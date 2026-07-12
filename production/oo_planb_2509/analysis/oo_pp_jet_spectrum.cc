#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

namespace {

constexpr std::array<double, 4> kRadii{0.1, 0.2, 0.4, 0.8};

struct Options {
  std::string output;
  std::vector<double> bins;
  double jet_abs_eta_max = 2.0;
  int progress_every = 50;
};

struct RunAccumulator {
  std::int64_t run_id = -1;
  std::uint64_t event_count = 0;
  double weight_sum = 0.0;
  double sigma_gen = std::numeric_limits<double>::quiet_NaN();
  std::vector<std::vector<std::uint64_t>> raw_counts;
  std::vector<std::vector<double>> weighted_density;
};

std::vector<double> parse_bins(const std::string &text) {
  std::vector<double> bins;
  std::istringstream stream(text);
  std::string token;
  while (std::getline(stream, token, ',')) {
    if (!token.empty()) {
      bins.push_back(std::stod(token));
    }
  }
  if (bins.size() < 2 || !std::is_sorted(bins.begin(), bins.end())) {
    throw std::runtime_error("--bins requires at least two strictly increasing values");
  }
  for (std::size_t index = 1; index < bins.size(); ++index) {
    if (!std::isfinite(bins[index - 1]) || bins[index] <= bins[index - 1]) {
      throw std::runtime_error("--bins requires finite, strictly increasing values");
    }
  }
  return bins;
}

Options parse_options(int argc, char **argv) {
  Options options;
  for (int index = 1; index < argc; ++index) {
    const std::string argument = argv[index];
    auto require_value = [&]() -> std::string {
      if (++index >= argc) {
        throw std::runtime_error("missing value after " + argument);
      }
      return argv[index];
    };
    if (argument == "--output") {
      options.output = require_value();
    } else if (argument == "--bins") {
      options.bins = parse_bins(require_value());
    } else if (argument == "--jet-abs-eta-max") {
      options.jet_abs_eta_max = std::stod(require_value());
    } else if (argument == "--progress-every") {
      options.progress_every = std::stoi(require_value());
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }
  if (options.output.empty() || options.bins.empty()) {
    throw std::runtime_error("--output and --bins are required");
  }
  if (!(options.jet_abs_eta_max > 0.0) || !std::isfinite(options.jet_abs_eta_max)) {
    throw std::runtime_error("--jet-abs-eta-max must be finite and positive");
  }
  return options;
}

class Analyzer {
 public:
  explicit Analyzer(Options options)
      : options_(std::move(options)), output_(options_.output) {
    if (!output_) {
      throw std::runtime_error("could not create " + options_.output);
    }
    fastjet::ClusterSequence::set_fastjet_banner_stream(nullptr);
    output_ << "run_id\tevent_count\tweight_sum\tsigma_gen\tradius\tbin_index"
               "\tpt_low\tpt_high\traw_jet_count\tweighted_density\n";
    output_ << std::setprecision(17);
  }

  void consume(std::istream &input) {
    std::string line;
    while (std::getline(input, line)) {
      if (!line.empty() && line.back() == '\r') {
        line.pop_back();
      }
      if (line.rfind("# OORUN ", 0) == 0) {
        start_run(std::stoll(line.substr(8)));
      } else if (line == "# OOENDRUN") {
        finish_run();
      } else if (line.rfind("# event", 0) == 0) {
        start_event();
      } else if (line.rfind("weight", 0) == 0) {
        parse_weight(line);
      } else if (line == "end") {
        finish_event();
      } else if (!line.empty() && line.front() != '#') {
        parse_particle(line);
      }
    }
    if (event_active_ || run_active_) {
      throw std::runtime_error("input ended before # OOENDRUN");
    }
  }

 private:
  void start_run(std::int64_t run_id) {
    if (run_active_ || event_active_) {
      throw std::runtime_error("nested # OORUN marker");
    }
    run_ = RunAccumulator{};
    run_.run_id = run_id;
    run_.raw_counts.assign(kRadii.size(),
                           std::vector<std::uint64_t>(options_.bins.size() - 1, 0));
    run_.weighted_density.assign(kRadii.size(),
                                 std::vector<double>(options_.bins.size() - 1, 0.0));
    run_active_ = true;
  }

  void start_event() {
    if (!run_active_) {
      throw std::runtime_error("# event outside # OORUN");
    }
    if (event_active_) {
      finish_event();
    }
    event_active_ = true;
    event_weight_ = std::numeric_limits<double>::quiet_NaN();
    event_sigma_gen_ = std::numeric_limits<double>::quiet_NaN();
    particles_.clear();
  }

  void parse_weight(const std::string &line) {
    if (!event_active_) {
      throw std::runtime_error("weight outside event");
    }
    std::istringstream stream(line);
    std::string key;
    stream >> key >> event_weight_;
    while (stream >> key) {
      if (key == "cross") {
        stream >> event_sigma_gen_;
        break;
      }
      std::string ignored;
      stream >> ignored;
    }
  }

  void parse_particle(const std::string &line) {
    if (!event_active_) {
      return;
    }
    double px = 0.0;
    double py = 0.0;
    double pz = 0.0;
    double mass = 0.0;
    int pdg_id = 0;
    int raw_label = 0;
    std::istringstream stream(line);
    if (!(stream >> px >> py >> pz >> mass >> pdg_id >> raw_label)) {
      return;
    }
    if (raw_label == -2) {
      return;
    }
    if (raw_label != 0) {
      throw std::runtime_error("pp input contains a nonzero final-state wake label");
    }
    const double pt = std::hypot(px, py);
    if (!(pt > 1e-12)) {
      return;
    }
    const double energy = std::sqrt(std::max(0.0, px * px + py * py + pz * pz + mass * mass));
    particles_.emplace_back(px, py, pz, energy);
  }

  int bin_index(double pt) const {
    for (std::size_t index = 0; index + 1 < options_.bins.size(); ++index) {
      if (pt > options_.bins[index] && pt <= options_.bins[index + 1]) {
        return static_cast<int>(index);
      }
    }
    return -1;
  }

  void finish_event() {
    if (!event_active_) {
      return;
    }
    if (!std::isfinite(event_weight_) || !(event_weight_ > 0.0) ||
        !std::isfinite(event_sigma_gen_) || !(event_sigma_gen_ > 0.0)) {
      throw std::runtime_error("event is missing a positive weight or sigmaGen");
    }
    for (std::size_t radius_index = 0; radius_index < kRadii.size(); ++radius_index) {
      const fastjet::JetDefinition definition(fastjet::antikt_algorithm,
                                              kRadii[radius_index], fastjet::E_scheme,
                                              fastjet::Best);
      const fastjet::ClusterSequence sequence(particles_, definition);
      const auto jets = fastjet::sorted_by_pt(sequence.inclusive_jets(0.0));
      for (const fastjet::PseudoJet &jet : jets) {
        if (!std::isfinite(jet.pt()) || !std::isfinite(jet.eta()) ||
            std::abs(jet.eta()) >= options_.jet_abs_eta_max) {
          continue;
        }
        const int index = bin_index(jet.pt());
        if (index < 0) {
          continue;
        }
        const std::size_t bin = static_cast<std::size_t>(index);
        ++run_.raw_counts[radius_index][bin];
        run_.weighted_density[radius_index][bin] +=
            event_weight_ / (options_.bins[bin + 1] - options_.bins[bin]);
      }
    }
    ++run_.event_count;
    run_.weight_sum += event_weight_;
    run_.sigma_gen = event_sigma_gen_;
    event_active_ = false;
    particles_.clear();
  }

  void finish_run() {
    if (!run_active_) {
      throw std::runtime_error("# OOENDRUN without # OORUN");
    }
    finish_event();
    if (run_.event_count == 0 || !(run_.weight_sum > 0.0) ||
        !std::isfinite(run_.sigma_gen)) {
      throw std::runtime_error("empty or invalid pp run");
    }
    for (std::size_t radius_index = 0; radius_index < kRadii.size(); ++radius_index) {
      for (std::size_t bin = 0; bin + 1 < options_.bins.size(); ++bin) {
        output_ << run_.run_id << '\t' << run_.event_count << '\t' << run_.weight_sum
                << '\t' << run_.sigma_gen << '\t' << kRadii[radius_index] << '\t'
                << bin << '\t' << options_.bins[bin] << '\t' << options_.bins[bin + 1]
                << '\t' << run_.raw_counts[radius_index][bin] << '\t'
                << run_.weighted_density[radius_index][bin] << '\n';
      }
    }
    ++completed_runs_;
    if (options_.progress_every > 0 && completed_runs_ % options_.progress_every == 0) {
      std::cerr << "processed_pp_runs=" << completed_runs_ << '\n';
    }
    run_active_ = false;
  }

  Options options_;
  std::ofstream output_;
  RunAccumulator run_;
  bool run_active_ = false;
  bool event_active_ = false;
  double event_weight_ = std::numeric_limits<double>::quiet_NaN();
  double event_sigma_gen_ = std::numeric_limits<double>::quiet_NaN();
  std::vector<fastjet::PseudoJet> particles_;
  std::uint64_t completed_runs_ = 0;
};

}  // namespace

int main(int argc, char **argv) {
  try {
    Analyzer analyzer(parse_options(argc, argv));
    analyzer.consume(std::cin);
    return 0;
  } catch (const std::exception &error) {
    std::cerr << "error: " << error.what() << '\n';
    return 1;
  }
}
