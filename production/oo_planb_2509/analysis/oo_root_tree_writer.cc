#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <exception>
#include <fstream>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

#include "Compression.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TNamed.h"
#include "TParameter.h"
#include "TTree.h"

#include "fastjet/ClusterSequence.hh"
#include "fastjet/JetDefinition.hh"
#include "fastjet/PseudoJet.hh"

namespace {

constexpr std::array<char, 8> kPairMagic{'O', 'O', 'P', 'A', 'I', 'R', '1', '\0'};
constexpr const char *kSchemaVersion = "oo-paired-root-v7";
constexpr double kHbarCGeVFm = 0.19732698;
constexpr double kFormationTimeCorrectedPtMin = 30.0;

#pragma pack(push, 1)
struct PairHeader {
  char magic[8];
  std::uint64_t pair_id;
  std::int32_t source_index;
  std::int64_t chunk_id;
  std::int32_t hydro_index;
  std::int64_t seed;
  std::int32_t event_number;
  double event_weight;
  double sigma_gen;
  double hard_x;
  double hard_y;
  std::uint32_t no_prehydro_count;
  std::uint32_t with_prehydro_count;
};

struct ParticleRecord {
  double px;
  double py;
  double pz;
  double mass;
  std::int32_t pdg_id;
  std::int32_t raw_label;
};
#pragma pack(pop)

static_assert(sizeof(PairHeader) == 84, "PairHeader protocol changed");
static_assert(sizeof(ParticleRecord) == 40, "ParticleRecord protocol changed");

struct Options {
  std::string output;
  double raw_jet_pt_min = 1.0;
  double jet_abs_eta_max = 5.0;
  double z_cut = 0.1;
  double beta = 0.0;
  double match_dr_fraction = 0.5;
  std::vector<std::pair<std::string, std::string>> sources;
};

struct EventMeta {
  std::uint64_t pair_id = 0;
  int source_index = -1;
  std::int64_t chunk_id = -1;
  int hydro_index = -1;
  std::int64_t seed = -1;
  int event_number = -1;
  double event_weight = std::numeric_limits<double>::quiet_NaN();
  double sigma_gen = std::numeric_limits<double>::quiet_NaN();
  double hard_x = std::numeric_limits<double>::quiet_NaN();
  double hard_y = std::numeric_limits<double>::quiet_NaN();
};

struct EventSummary {
  int n_hadrons = 0;
  int n_normal = 0;
  int n_positive_wake = 0;
  int n_negative_wake = 0;
  int n_negative_thermal = 0;
  int n_hadronized_holes = 0;
  int n_hard_markers = 0;
  double sum_pt = 0.0;
  double signed_sum_pt = 0.0;
};

struct HadronInfo : public fastjet::PseudoJet::UserInfoBase {
  HadronInfo(int raw_label_in, int pdg_id_in, int input_index_in)
      : raw_label(raw_label_in), pdg_id(pdg_id_in), input_index(input_index_in) {}
  int raw_label;
  int pdg_id;
  int input_index;
};

double particle_energy(const ParticleRecord &particle) {
  return std::sqrt(std::max(0.0, particle.px * particle.px + particle.py * particle.py +
                                    particle.pz * particle.pz + particle.mass * particle.mass));
}

double pseudo_eta(double px, double py, double pz) {
  const double pt = std::hypot(px, py);
  if (pt == 0.0) {
    if (pz == 0.0) {
      return std::numeric_limits<double>::quiet_NaN();
    }
    return std::copysign(std::numeric_limits<double>::infinity(), pz);
  }
  return std::asinh(pz / pt);
}

double signed_mass(double px, double py, double pz, double energy) {
  const double mass_squared = energy * energy - px * px - py * py - pz * pz;
  return std::copysign(std::sqrt(std::abs(mass_squared)), mass_squared);
}

double delta_phi(double first, double second) {
  return std::remainder(first - second, 2.0 * M_PI);
}

double delta_r(double first_y, double first_phi, double second_y, double second_phi) {
  return std::hypot(first_y - second_y, delta_phi(first_phi, second_phi));
}

template <typename T>
bool read_record(std::istream &input, T &record, bool allow_clean_eof = false) {
  input.read(reinterpret_cast<char *>(&record), sizeof(T));
  const auto bytes_read = input.gcount();
  if (bytes_read == 0 && allow_clean_eof && input.eof()) {
    return false;
  }
  if (bytes_read != static_cast<std::streamsize>(sizeof(T))) {
    throw std::runtime_error("truncated binary event stream");
  }
  return true;
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
    } else if (argument == "--raw-jet-pt-min") {
      options.raw_jet_pt_min = std::stod(require_value());
    } else if (argument == "--jet-abs-eta-max") {
      options.jet_abs_eta_max = std::stod(require_value());
    } else if (argument == "--z-cut") {
      options.z_cut = std::stod(require_value());
    } else if (argument == "--beta") {
      options.beta = std::stod(require_value());
    } else if (argument == "--match-dr-fraction") {
      options.match_dr_fraction = std::stod(require_value());
    } else if (argument == "--source") {
      const std::string source = require_value();
      const auto separator = source.find('\t');
      if (separator == std::string::npos) {
        throw std::runtime_error("--source must be NAME<TAB>PATH");
      }
      options.sources.emplace_back(source.substr(0, separator), source.substr(separator + 1));
    } else {
      throw std::runtime_error("unknown argument: " + argument);
    }
  }
  if (options.output.empty()) {
    throw std::runtime_error("--output is required");
  }
  if (options.raw_jet_pt_min < 0.0 || options.jet_abs_eta_max <= 0.0 || options.z_cut < 0.0 ||
      options.match_dr_fraction <= 0.0) {
    throw std::runtime_error("invalid jet configuration");
  }
  return options;
}

void branch_event_metadata(TTree *tree, EventMeta &metadata, bool &prehydro_enabled) {
  tree->Branch("pairId", &metadata.pair_id);
  tree->Branch("sourceIndex", &metadata.source_index);
  tree->Branch("chunkId", &metadata.chunk_id);
  tree->Branch("hydroIndex", &metadata.hydro_index);
  tree->Branch("seed", &metadata.seed);
  tree->Branch("eventNumber", &metadata.event_number);
  tree->Branch("prehydroEnabled", &prehydro_enabled);
  tree->Branch("eventWeight", &metadata.event_weight);
  tree->Branch("sigmaGen", &metadata.sigma_gen);
  tree->Branch("hardX", &metadata.hard_x);
  tree->Branch("hardY", &metadata.hard_y);
}

class HadronTree {
 public:
  HadronTree(TDirectory *directory, bool prehydro_enabled_in)
      : prehydro_enabled_(prehydro_enabled_in) {
    directory->cd();
    tree_ = new TTree("Hadrons", "One entry per accepted paired HYBRID event");
    tree_->SetAutoFlush(-10000000);
    branch_event_metadata(tree_, metadata_, prehydro_enabled_);
    tree_->Branch("nHadrons", &summary_.n_hadrons);
    tree_->Branch("nNormal", &summary_.n_normal);
    tree_->Branch("nPositiveWake", &summary_.n_positive_wake);
    tree_->Branch("nNegativeWake", &summary_.n_negative_wake);
    tree_->Branch("nNegativeThermal", &summary_.n_negative_thermal);
    tree_->Branch("nHadronizedHoles", &summary_.n_hadronized_holes);
    tree_->Branch("nHardMarkers", &summary_.n_hard_markers);
    tree_->Branch("sumPt", &summary_.sum_pt);
    tree_->Branch("signedSumPt", &summary_.signed_sum_pt);
    tree_->Branch("signedPx", &signed_px_);
    tree_->Branch("signedPy", &signed_py_);
    tree_->Branch("signedPz", &signed_pz_);
    tree_->Branch("signedE", &signed_energy_);
    tree_->Branch("hadronPt", &hadron_pt_);
    tree_->Branch("hadronEta", &hadron_eta_);
    tree_->Branch("hadronPhi", &hadron_phi_);
    tree_->Branch("hadronMass", &hadron_mass_);
    tree_->Branch("hadronE", &hadron_energy_);
    tree_->Branch("hadronPx", &hadron_px_);
    tree_->Branch("hadronPy", &hadron_py_);
    tree_->Branch("hadronPz", &hadron_pz_);
    tree_->Branch("hadronStatus", &hadron_status_);
    tree_->Branch("hadronRawLabel", &hadron_raw_label_);
    tree_->Branch("hadronID", &hadron_id_);
  }

  EventSummary fill(const EventMeta &metadata, const std::vector<ParticleRecord> &particles) {
    metadata_ = metadata;
    summary_ = EventSummary{};
    signed_px_ = 0.0;
    signed_py_ = 0.0;
    signed_pz_ = 0.0;
    signed_energy_ = 0.0;
    clear_vectors();
    hadron_pt_.reserve(particles.size());
    hadron_eta_.reserve(particles.size());
    hadron_phi_.reserve(particles.size());
    hadron_mass_.reserve(particles.size());
    hadron_energy_.reserve(particles.size());
    hadron_px_.reserve(particles.size());
    hadron_py_.reserve(particles.size());
    hadron_pz_.reserve(particles.size());
    hadron_status_.reserve(particles.size());
    hadron_raw_label_.reserve(particles.size());
    hadron_id_.reserve(particles.size());

    for (const ParticleRecord &particle : particles) {
      if (particle.raw_label == -2) {
        ++summary_.n_hard_markers;
        continue;
      }
      if (particle.raw_label < 0 || particle.raw_label > 3) {
        throw std::runtime_error("unsupported HYBRID particle label");
      }
      const bool negative = particle.raw_label == 2 || particle.raw_label == 3;
      const int status = negative ? -1 : (particle.raw_label == 1 ? 1 : 0);
      const double sign = negative ? -1.0 : 1.0;
      const double pt = std::hypot(particle.px, particle.py);
      const double energy = particle_energy(particle);

      ++summary_.n_hadrons;
      if (particle.raw_label == 0) {
        ++summary_.n_normal;
      } else if (particle.raw_label == 1) {
        ++summary_.n_positive_wake;
      } else {
        ++summary_.n_negative_wake;
        if (particle.raw_label == 2) {
          ++summary_.n_negative_thermal;
        } else {
          ++summary_.n_hadronized_holes;
        }
      }
      summary_.sum_pt += pt;
      summary_.signed_sum_pt += sign * pt;
      signed_px_ += sign * particle.px;
      signed_py_ += sign * particle.py;
      signed_pz_ += sign * particle.pz;
      signed_energy_ += sign * energy;

      hadron_pt_.push_back(static_cast<float>(pt));
      hadron_eta_.push_back(static_cast<float>(pseudo_eta(particle.px, particle.py, particle.pz)));
      hadron_phi_.push_back(static_cast<float>(std::atan2(particle.py, particle.px)));
      hadron_mass_.push_back(static_cast<float>(particle.mass));
      hadron_energy_.push_back(static_cast<float>(energy));
      hadron_px_.push_back(static_cast<float>(particle.px));
      hadron_py_.push_back(static_cast<float>(particle.py));
      hadron_pz_.push_back(static_cast<float>(particle.pz));
      hadron_status_.push_back(status);
      hadron_raw_label_.push_back(particle.raw_label);
      hadron_id_.push_back(particle.pdg_id);
    }
    tree_->Fill();
    return summary_;
  }

 private:
  void clear_vectors() {
    hadron_pt_.clear();
    hadron_eta_.clear();
    hadron_phi_.clear();
    hadron_mass_.clear();
    hadron_energy_.clear();
    hadron_px_.clear();
    hadron_py_.clear();
    hadron_pz_.clear();
    hadron_status_.clear();
    hadron_raw_label_.clear();
    hadron_id_.clear();
  }

  TTree *tree_ = nullptr;
  EventMeta metadata_;
  bool prehydro_enabled_ = false;
  EventSummary summary_;
  double signed_px_ = 0.0;
  double signed_py_ = 0.0;
  double signed_pz_ = 0.0;
  double signed_energy_ = 0.0;
  std::vector<float> hadron_pt_;
  std::vector<float> hadron_eta_;
  std::vector<float> hadron_phi_;
  std::vector<float> hadron_mass_;
  std::vector<float> hadron_energy_;
  std::vector<float> hadron_px_;
  std::vector<float> hadron_py_;
  std::vector<float> hadron_pz_;
  std::vector<int> hadron_status_;
  std::vector<int> hadron_raw_label_;
  std::vector<int> hadron_id_;
};

struct JetRecord {
  double eta = std::numeric_limits<double>::quiet_NaN();
  double y = std::numeric_limits<double>::quiet_NaN();
  double phi = std::numeric_limits<double>::quiet_NaN();
  double pt = std::numeric_limits<double>::quiet_NaN();
  double mass = std::numeric_limits<double>::quiet_NaN();
  double energy = std::numeric_limits<double>::quiet_NaN();
  double px = 0.0;
  double py = 0.0;
  double pz = 0.0;
  double raw_eta = std::numeric_limits<double>::quiet_NaN();
  double raw_y = std::numeric_limits<double>::quiet_NaN();
  double raw_phi = std::numeric_limits<double>::quiet_NaN();
  double raw_pt = 0.0;
  double raw_mass = 0.0;
  double zg = std::numeric_limits<double>::quiet_NaN();
  double rg = std::numeric_limits<double>::quiet_NaN();
  double sd_pt = std::numeric_limits<double>::quiet_NaN();
  double sd_mass = std::numeric_limits<double>::quiet_NaN();
  int soft_drop_valid = 0;
  int n_sd = 0;
  int mult = 0;
  int total_mult = 0;
  double pt_d = std::numeric_limits<double>::quiet_NaN();
  double effective_multiplicity = std::numeric_limits<double>::quiet_NaN();
  double leading_fraction = std::numeric_limits<double>::quiet_NaN();
  double normal_pt_d = std::numeric_limits<double>::quiet_NaN();
  double normal_effective_multiplicity = std::numeric_limits<double>::quiet_NaN();
  double leading_normal_fraction = std::numeric_limits<double>::quiet_NaN();
  double girth = std::numeric_limits<double>::quiet_NaN();
  double signed_girth = std::numeric_limits<double>::quiet_NaN();
  double max_kt = std::numeric_limits<double>::quiet_NaN();
  double max_kt_z = std::numeric_limits<double>::quiet_NaN();
  double max_kt_rg = std::numeric_limits<double>::quiet_NaN();
  double max_kt_primary = std::numeric_limits<double>::quiet_NaN();
  std::vector<float> formation_tau_f;
  std::vector<float> formation_tau_f_small_angle;
  std::vector<float> formation_z;
  std::vector<float> formation_theta;
  std::vector<float> formation_delta_r;
  std::vector<float> formation_kt;
  std::vector<float> formation_parent_energy;
  int formation_invalid_splits = 0;
  int formation_hardest_valid = 0;
  double formation_hardest_tau_f = std::numeric_limits<double>::quiet_NaN();
  double formation_hardest_tau_f_small_angle = std::numeric_limits<double>::quiet_NaN();
  double formation_hardest_z = std::numeric_limits<double>::quiet_NaN();
  double formation_hardest_theta = std::numeric_limits<double>::quiet_NaN();
  double formation_hardest_delta_r = std::numeric_limits<double>::quiet_NaN();
  double formation_hardest_kt = std::numeric_limits<double>::quiet_NaN();
  double formation_hardest_parent_energy = std::numeric_limits<double>::quiet_NaN();
  double normal_pt = 0.0;
  double positive_wake_pt = 0.0;
  double negative_wake_pt = 0.0;
  double wake_fraction = std::numeric_limits<double>::quiet_NaN();
  int n_normal = 0;
  int n_positive_wake = 0;
  int n_negative_wake = 0;
  int n_negative_thermal = 0;
  int n_hadronized_holes = 0;
  int hard_parton_id = 0;
  double hard_parton_pt = std::numeric_limits<double>::quiet_NaN();
  double hard_parton_dr = std::numeric_limits<double>::quiet_NaN();
  int pair_match_index = -1;
  double pair_match_dr = std::numeric_limits<double>::quiet_NaN();
  double pair_match_other_pt = std::numeric_limits<double>::quiet_NaN();
  int pair_match_other_hard_parton_id = 0;
};

struct SoftDropResult {
  bool valid = false;
  double zg = std::numeric_limits<double>::quiet_NaN();
  double rg = std::numeric_limits<double>::quiet_NaN();
  double pt = std::numeric_limits<double>::quiet_NaN();
  double mass = std::numeric_limits<double>::quiet_NaN();
  int n_sd = 0;
  double max_kt = std::numeric_limits<double>::quiet_NaN();
  double max_kt_z = std::numeric_limits<double>::quiet_NaN();
  double max_kt_rg = std::numeric_limits<double>::quiet_NaN();
  double max_kt_primary = std::numeric_limits<double>::quiet_NaN();
  std::vector<float> formation_tau_f;
  std::vector<float> formation_tau_f_small_angle;
  std::vector<float> formation_z;
  std::vector<float> formation_theta;
  std::vector<float> formation_delta_r;
  std::vector<float> formation_kt;
  std::vector<float> formation_parent_energy;
  int formation_invalid_splits = 0;
  int formation_hardest_index = -1;
  double formation_hardest_kt = -std::numeric_limits<double>::infinity();
};

struct FormationTimeRecord {
  double tau_f = std::numeric_limits<double>::quiet_NaN();
  double tau_f_small_angle = std::numeric_limits<double>::quiet_NaN();
  double z = std::numeric_limits<double>::quiet_NaN();
  double theta = std::numeric_limits<double>::quiet_NaN();
  double delta_r = std::numeric_limits<double>::quiet_NaN();
  double kt = std::numeric_limits<double>::quiet_NaN();
  double parent_energy = std::numeric_limits<double>::quiet_NaN();
};

bool calculate_formation_time(const fastjet::PseudoJet &parent,
                              const fastjet::PseudoJet &first,
                              const fastjet::PseudoJet &second,
                              FormationTimeRecord &record) {
  const double first_momentum =
      std::sqrt(first.px() * first.px() + first.py() * first.py() + first.pz() * first.pz());
  const double second_momentum = std::sqrt(second.px() * second.px() + second.py() * second.py() +
                                           second.pz() * second.pz());
  const double parent_energy = parent.e();
  if (!(parent_energy > 0.0 && first.e() > 0.0 && second.e() > 0.0 &&
        first_momentum > 0.0 && second_momentum > 0.0)) {
    return false;
  }
  const double z_first = first.e() / parent_energy;
  const double z_second = second.e() / parent_energy;
  if (!(z_first > 0.0 && z_second > 0.0)) {
    return false;
  }
  const double cosine = std::clamp(
      (first.px() * second.px() + first.py() * second.py() + first.pz() * second.pz()) /
          (first_momentum * second_momentum),
      -1.0, 1.0);
  const double one_minus_cosine = 1.0 - cosine;
  const double delta_r = first.delta_R(second);
  const double exact_denominator =
      parent_energy * z_first * z_second * one_minus_cosine;
  const double small_angle_denominator =
      0.5 * parent_energy * z_first * z_second * delta_r * delta_r;
  if (!(exact_denominator > 0.0 && small_angle_denominator > 0.0)) {
    return false;
  }
  record.tau_f = kHbarCGeVFm / exact_denominator;
  record.tau_f_small_angle = kHbarCGeVFm / small_angle_denominator;
  record.z = std::min(z_first, z_second);
  record.theta = std::acos(cosine);
  record.delta_r = delta_r;
  record.kt = std::min(first.pt(), second.pt()) * delta_r;
  record.parent_energy = parent_energy;
  return std::isfinite(record.tau_f) && std::isfinite(record.tau_f_small_angle) &&
         std::isfinite(record.z) && std::isfinite(record.theta) &&
         std::isfinite(record.delta_r) && std::isfinite(record.kt) &&
         std::isfinite(record.parent_energy);
}

SoftDropResult calculate_declustering(const std::vector<fastjet::PseudoJet> &constituents,
                                      double radius, double z_cut, double beta,
                                      bool retain_formation_time) {
  SoftDropResult result;
  if (constituents.size() < 2) {
    return result;
  }
  const fastjet::JetDefinition ca_definition(fastjet::cambridge_algorithm, 2.0 * radius + 1e-6,
                                             fastjet::E_scheme, fastjet::Best);
  const fastjet::ClusterSequence ca_sequence(constituents, ca_definition);
  auto ca_jets = fastjet::sorted_by_pt(ca_sequence.inclusive_jets(0.0));
  if (ca_jets.empty()) {
    return result;
  }
  fastjet::PseudoJet root = ca_jets.front();

  std::vector<fastjet::PseudoJet> stack{root};
  while (!stack.empty()) {
    const fastjet::PseudoJet node = stack.back();
    stack.pop_back();
    fastjet::PseudoJet first;
    fastjet::PseudoJet second;
    if (!node.has_parents(first, second)) {
      continue;
    }
    const double denominator = first.pt() + second.pt();
    const double z = denominator > 0.0 ? std::min(first.pt(), second.pt()) / denominator : 0.0;
    const double rg = first.delta_R(second);
    const double kt = std::min(first.pt(), second.pt()) * rg;
    if (!std::isfinite(result.max_kt) || kt > result.max_kt) {
      result.max_kt = kt;
      result.max_kt_z = z;
      result.max_kt_rg = rg;
    }
    if (retain_formation_time) {
      FormationTimeRecord formation;
      if (calculate_formation_time(node, first, second, formation)) {
        result.formation_tau_f.push_back(static_cast<float>(formation.tau_f));
        result.formation_tau_f_small_angle.push_back(
            static_cast<float>(formation.tau_f_small_angle));
        result.formation_z.push_back(static_cast<float>(formation.z));
        result.formation_theta.push_back(static_cast<float>(formation.theta));
        result.formation_delta_r.push_back(static_cast<float>(formation.delta_r));
        result.formation_kt.push_back(static_cast<float>(formation.kt));
        result.formation_parent_energy.push_back(static_cast<float>(formation.parent_energy));
        const int candidate = static_cast<int>(result.formation_kt.size()) - 1;
        if (formation.kt > result.formation_hardest_kt) {
          result.formation_hardest_index = candidate;
          result.formation_hardest_kt = formation.kt;
        }
      } else {
        ++result.formation_invalid_splits;
      }
    }
    stack.push_back(first);
    stack.push_back(second);
  }
  if (retain_formation_time &&
      result.formation_tau_f.size() +
              static_cast<std::size_t>(result.formation_invalid_splits) !=
          constituents.size() - 1) {
    throw std::runtime_error("C/A formation-time tree does not have Nconstituent-1 splits");
  }

  fastjet::PseudoJet current = root;
  while (true) {
    fastjet::PseudoJet first;
    fastjet::PseudoJet second;
    if (!current.has_parents(first, second)) {
      break;
    }
    if (second.pt() > first.pt()) {
      std::swap(first, second);
    }
    const double denominator = first.pt() + second.pt();
    const double z = denominator > 0.0 ? second.pt() / denominator : 0.0;
    const double rg = first.delta_R(second);
    const double kt = second.pt() * rg;
    if (!std::isfinite(result.max_kt_primary) || kt > result.max_kt_primary) {
      result.max_kt_primary = kt;
    }
    const double angular_factor = beta == 0.0 ? 1.0 : std::pow(rg / radius, beta);
    if (z > z_cut * angular_factor) {
      ++result.n_sd;
      if (!result.valid) {
        result.valid = true;
        result.zg = z;
        result.rg = rg;
        result.pt = current.pt();
        result.mass = current.m();
      }
    }
    current = first;
  }
  return result;
}

std::vector<JetRecord> make_jets(const std::vector<ParticleRecord> &particles, double radius,
                                 const Options &options) {
  constexpr double kGhostScale = 1e-100;
  std::vector<fastjet::PseudoJet> clustering_inputs;
  clustering_inputs.reserve(particles.size());
  for (std::size_t index = 0; index < particles.size(); ++index) {
    const ParticleRecord &particle = particles[index];
    if (particle.raw_label == -2) {
      continue;
    }
    const double pt = std::hypot(particle.px, particle.py);
    if (pt <= 1e-12) {
      continue;
    }
    const bool negative = particle.raw_label == 2 || particle.raw_label == 3;
    const double scale = negative ? kGhostScale : 1.0;
    fastjet::PseudoJet pseudojet(scale * particle.px, scale * particle.py, scale * particle.pz,
                                scale * particle_energy(particle));
    pseudojet.set_user_info(new HadronInfo(particle.raw_label, particle.pdg_id,
                                           static_cast<int>(index)));
    if (particle.raw_label >= 0 && particle.raw_label <= 3) {
      clustering_inputs.push_back(pseudojet);
    }
  }
  if (clustering_inputs.empty()) {
    return {};
  }

  const fastjet::JetDefinition anti_kt_definition(fastjet::antikt_algorithm, radius,
                                                  fastjet::E_scheme, fastjet::Best);
  const fastjet::ClusterSequence cluster_sequence(clustering_inputs, anti_kt_definition);
  const auto raw_jets = fastjet::sorted_by_pt(cluster_sequence.inclusive_jets(0.0));

  std::vector<JetRecord> records;
  records.reserve(raw_jets.size());
  for (std::size_t jet_index = 0; jet_index < raw_jets.size(); ++jet_index) {
    const fastjet::PseudoJet &raw_jet = raw_jets[jet_index];
    if (raw_jet.pt() < options.raw_jet_pt_min ||
        std::abs(raw_jet.eta()) > options.jet_abs_eta_max) {
      continue;
    }
    JetRecord record;
    record.raw_eta = raw_jet.eta();
    record.raw_y = raw_jet.rap();
    record.raw_phi = raw_jet.phi_std();
    record.raw_pt = raw_jet.pt();
    record.raw_mass = raw_jet.m();

    const auto all_constituents = raw_jet.constituents();
    std::vector<fastjet::PseudoJet> constituents;
    std::vector<const ParticleRecord *> assigned_negative;
    constituents.reserve(all_constituents.size());
    assigned_negative.reserve(all_constituents.size());
    for (const fastjet::PseudoJet &constituent : all_constituents) {
      const auto &info = constituent.user_info<HadronInfo>();
      if (info.raw_label == 0 || info.raw_label == 1) {
        constituents.push_back(constituent);
      } else {
        assigned_negative.push_back(&particles.at(static_cast<std::size_t>(info.input_index)));
      }
    }
    record.mult = static_cast<int>(constituents.size());
    double scalar_pt = 0.0;
    double sum_pt_squared = 0.0;
    double leading_pt = 0.0;
    double normal_sum_pt_squared = 0.0;
    double leading_normal_pt = 0.0;
    double girth_numerator = 0.0;
    for (const fastjet::PseudoJet &constituent : constituents) {
      const auto &info = constituent.user_info<HadronInfo>();
      const double pt = constituent.pt();
      scalar_pt += pt;
      sum_pt_squared += pt * pt;
      leading_pt = std::max(leading_pt, pt);
      girth_numerator += pt * constituent.delta_R(raw_jet);
      if (info.raw_label == 0) {
        ++record.n_normal;
        record.normal_pt += pt;
        normal_sum_pt_squared += pt * pt;
        leading_normal_pt = std::max(leading_normal_pt, pt);
      } else if (info.raw_label == 1) {
        ++record.n_positive_wake;
        record.positive_wake_pt += pt;
      }
    }
    if (scalar_pt > 0.0) {
      record.pt_d = std::sqrt(sum_pt_squared) / scalar_pt;
      record.effective_multiplicity = scalar_pt * scalar_pt / sum_pt_squared;
      record.leading_fraction = leading_pt / scalar_pt;
      record.girth = girth_numerator / scalar_pt;
      record.wake_fraction = record.positive_wake_pt / scalar_pt;
    }
    if (record.normal_pt > 0.0 && normal_sum_pt_squared > 0.0) {
      record.normal_pt_d = std::sqrt(normal_sum_pt_squared) / record.normal_pt;
      record.normal_effective_multiplicity =
          record.normal_pt * record.normal_pt / normal_sum_pt_squared;
      record.leading_normal_fraction = leading_normal_pt / record.normal_pt;
    }

    double negative_px = 0.0;
    double negative_py = 0.0;
    double negative_pz = 0.0;
    double negative_energy = 0.0;
    double signed_girth_numerator = girth_numerator;
    for (const ParticleRecord *negative_particle : assigned_negative) {
      const ParticleRecord &negative = *negative_particle;
      ++record.n_negative_wake;
      if (negative.raw_label == 2) {
        ++record.n_negative_thermal;
      } else {
        ++record.n_hadronized_holes;
      }
      const double negative_pt = std::hypot(negative.px, negative.py);
      const fastjet::PseudoJet negative_pseudojet(negative.px, negative.py, negative.pz,
                                                  particle_energy(negative));
      record.negative_wake_pt += negative_pt;
      negative_px += negative.px;
      negative_py += negative.py;
      negative_pz += negative.pz;
      negative_energy += negative_pseudojet.e();
      signed_girth_numerator -= negative_pt * negative_pseudojet.delta_R(raw_jet);
    }
    record.total_mult =
        record.n_normal + record.n_positive_wake - record.n_negative_wake;
    const double signed_scalar_pt = scalar_pt - record.negative_wake_pt;
    if (signed_scalar_pt != 0.0) {
      record.signed_girth = signed_girth_numerator / signed_scalar_pt;
    }
    if (scalar_pt > 0.0) {
      record.wake_fraction =
          (record.positive_wake_pt - record.negative_wake_pt) / scalar_pt;
    }

    record.px = raw_jet.px() - negative_px;
    record.py = raw_jet.py() - negative_py;
    record.pz = raw_jet.pz() - negative_pz;
    record.energy = raw_jet.e() - negative_energy;
    record.pt = std::hypot(record.px, record.py);
    record.eta = pseudo_eta(record.px, record.py, record.pz);
    record.phi = std::atan2(record.py, record.px);
    record.mass = signed_mass(record.px, record.py, record.pz, record.energy);
    if (record.energy > std::abs(record.pz)) {
      record.y = 0.5 * std::log((record.energy + record.pz) / (record.energy - record.pz));
    }

    const bool retain_formation_time =
        radius >= 0.4 - 1e-12 && record.pt > kFormationTimeCorrectedPtMin;
    const SoftDropResult declustering = calculate_declustering(
        constituents, radius, options.z_cut, options.beta, retain_formation_time);
    record.soft_drop_valid = declustering.valid ? 1 : 0;
    record.zg = declustering.zg;
    record.rg = declustering.rg;
    record.sd_pt = declustering.pt;
    record.sd_mass = declustering.mass;
    record.n_sd = declustering.n_sd;
    record.max_kt = declustering.max_kt;
    record.max_kt_z = declustering.max_kt_z;
    record.max_kt_rg = declustering.max_kt_rg;
    record.max_kt_primary = declustering.max_kt_primary;
    if (retain_formation_time) {
      record.formation_tau_f = declustering.formation_tau_f;
      record.formation_tau_f_small_angle = declustering.formation_tau_f_small_angle;
      record.formation_z = declustering.formation_z;
      record.formation_theta = declustering.formation_theta;
      record.formation_delta_r = declustering.formation_delta_r;
      record.formation_kt = declustering.formation_kt;
      record.formation_parent_energy = declustering.formation_parent_energy;
      record.formation_invalid_splits = declustering.formation_invalid_splits;
      if (declustering.formation_hardest_index >= 0) {
        const std::size_t hardest =
            static_cast<std::size_t>(declustering.formation_hardest_index);
        record.formation_hardest_valid = 1;
        record.formation_hardest_tau_f = declustering.formation_tau_f.at(hardest);
        record.formation_hardest_tau_f_small_angle =
            declustering.formation_tau_f_small_angle.at(hardest);
        record.formation_hardest_z = declustering.formation_z.at(hardest);
        record.formation_hardest_theta = declustering.formation_theta.at(hardest);
        record.formation_hardest_delta_r = declustering.formation_delta_r.at(hardest);
        record.formation_hardest_kt = declustering.formation_kt.at(hardest);
        record.formation_hardest_parent_energy =
            declustering.formation_parent_energy.at(hardest);
      }
    }
    records.push_back(record);
  }
  return records;
}

void match_hard_partons(std::vector<JetRecord> &jets,
                        const std::vector<ParticleRecord> &particles, double radius) {
  struct Candidate {
    double distance;
    int jet_index;
    int particle_index;
  };
  std::vector<std::size_t> marker_indices;
  for (std::size_t particle_index = 0; particle_index < particles.size(); ++particle_index) {
    if (particles[particle_index].raw_label == -2) {
      marker_indices.push_back(particle_index);
    }
  }
  std::vector<Candidate> candidates;
  for (std::size_t jet_index = 0; jet_index < jets.size(); ++jet_index) {
    for (const std::size_t particle_index : marker_indices) {
      const ParticleRecord &particle = particles[particle_index];
      const double marker_eta = pseudo_eta(particle.px, particle.py, particle.pz);
      const double marker_phi = std::atan2(particle.py, particle.px);
      const double distance =
          delta_r(jets[jet_index].raw_eta, jets[jet_index].raw_phi, marker_eta, marker_phi);
      if (std::isfinite(distance) && distance < radius) {
        candidates.push_back(Candidate{distance, static_cast<int>(jet_index),
                                       static_cast<int>(particle_index)});
      }
    }
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &first, const Candidate &second) {
              return first.distance < second.distance;
            });
  std::vector<bool> jet_used(jets.size(), false);
  std::vector<bool> particle_used(particles.size(), false);
  for (const Candidate &candidate : candidates) {
    if (jet_used[candidate.jet_index] || particle_used[candidate.particle_index]) {
      continue;
    }
    jet_used[candidate.jet_index] = true;
    particle_used[candidate.particle_index] = true;
    JetRecord &jet = jets[candidate.jet_index];
    const ParticleRecord &particle = particles[candidate.particle_index];
    jet.hard_parton_id = particle.pdg_id;
    jet.hard_parton_pt = std::hypot(particle.px, particle.py);
    jet.hard_parton_dr = candidate.distance;
  }
}

void match_jets(std::vector<JetRecord> &no_prehydro, std::vector<JetRecord> &with_prehydro,
                double radius, double max_fraction) {
  struct Candidate {
    double distance;
    int no_index;
    int with_index;
  };
  std::vector<Candidate> candidates;
  for (std::size_t no_index = 0; no_index < no_prehydro.size(); ++no_index) {
    for (std::size_t with_index = 0; with_index < with_prehydro.size(); ++with_index) {
      const double distance =
          delta_r(no_prehydro[no_index].raw_y, no_prehydro[no_index].raw_phi,
                  with_prehydro[with_index].raw_y, with_prehydro[with_index].raw_phi);
      if (distance < max_fraction * radius) {
        candidates.push_back(
            Candidate{distance, static_cast<int>(no_index), static_cast<int>(with_index)});
      }
    }
  }
  std::sort(candidates.begin(), candidates.end(),
            [](const Candidate &first, const Candidate &second) {
              return first.distance < second.distance;
            });
  std::vector<bool> no_used(no_prehydro.size(), false);
  std::vector<bool> with_used(with_prehydro.size(), false);
  for (const Candidate &candidate : candidates) {
    if (no_used[candidate.no_index] || with_used[candidate.with_index]) {
      continue;
    }
    no_used[candidate.no_index] = true;
    with_used[candidate.with_index] = true;
    JetRecord &no_jet = no_prehydro[candidate.no_index];
    JetRecord &with_jet = with_prehydro[candidate.with_index];
    no_jet.pair_match_index = candidate.with_index;
    no_jet.pair_match_dr = candidate.distance;
    no_jet.pair_match_other_pt = with_jet.pt;
    no_jet.pair_match_other_hard_parton_id = with_jet.hard_parton_id;
    with_jet.pair_match_index = candidate.no_index;
    with_jet.pair_match_dr = candidate.distance;
    with_jet.pair_match_other_pt = no_jet.pt;
    with_jet.pair_match_other_hard_parton_id = no_jet.hard_parton_id;
  }
}

class RadiusBranches {
 public:
  RadiusBranches(TTree *tree, std::string prefix, bool store_formation_time)
      : prefix_(std::move(prefix)), store_formation_time_(store_formation_time) {
    branch(tree, "Eta", eta_);
    branch(tree, "Y", y_);
    branch(tree, "Phi", phi_);
    branch(tree, "Pt", pt_);
    branch(tree, "mass", mass_);
    branch(tree, "E", energy_);
    branch(tree, "Px", px_);
    branch(tree, "Py", py_);
    branch(tree, "Pz", pz_);
    branch(tree, "RawEta", raw_eta_);
    branch(tree, "RawY", raw_y_);
    branch(tree, "RawPhi", raw_phi_);
    branch(tree, "RawPt", raw_pt_);
    branch(tree, "RawMass", raw_mass_);
    branch(tree, "Zg", zg_);
    branch(tree, "Rg", rg_);
    branch(tree, "SDPt", sd_pt_);
    branch(tree, "SDmass", sd_mass_);
    branch(tree, "SoftDropValid", soft_drop_valid_);
    branch(tree, "NSD", n_sd_);
    branch(tree, "Mult", mult_);
    branch(tree, "TotalMult", total_mult_);
    branch(tree, "PtD", pt_d_);
    branch(tree, "EffectiveMultiplicity", effective_multiplicity_);
    branch(tree, "LeadingFraction", leading_fraction_);
    branch(tree, "NormalPtD", normal_pt_d_);
    branch(tree, "NormalEffectiveMultiplicity", normal_effective_multiplicity_);
    branch(tree, "LeadingNormalFraction", leading_normal_fraction_);
    branch(tree, "G", girth_);
    branch(tree, "SignedG", signed_girth_);
    branch(tree, "MaxKt", max_kt_);
    branch(tree, "MaxKtZ", max_kt_z_);
    branch(tree, "MaxKtRg", max_kt_rg_);
    branch(tree, "MaxKtPrimary", max_kt_primary_);
    branch(tree, "NormalPt", normal_pt_);
    branch(tree, "PositiveWakePt", positive_wake_pt_);
    branch(tree, "NegativeWakePt", negative_wake_pt_);
    branch(tree, "WakeFraction", wake_fraction_);
    branch(tree, "NNormal", n_normal_);
    branch(tree, "NPositiveWake", n_positive_wake_);
    branch(tree, "NNegativeWake", n_negative_wake_);
    branch(tree, "NNegativeThermal", n_negative_thermal_);
    branch(tree, "NHadronizedHoles", n_hadronized_holes_);
    branch(tree, "HardPartonId", hard_parton_id_);
    branch(tree, "HardPartonPt", hard_parton_pt_);
    branch(tree, "HardPartonDR", hard_parton_dr_);
    branch(tree, "PairMatchIndex", pair_match_index_);
    branch(tree, "PairMatchDR", pair_match_dr_);
    branch(tree, "PairMatchOtherPt", pair_match_other_pt_);
    branch(tree, "PairMatchOtherHardPartonId", pair_match_other_hard_parton_id_);
    if (store_formation_time_) {
      branch(tree, "FormationTauF", formation_tau_f_);
      branch(tree, "FormationTauFSmallAngle", formation_tau_f_small_angle_);
      branch(tree, "FormationZ", formation_z_);
      branch(tree, "FormationTheta", formation_theta_);
      branch(tree, "FormationDeltaR", formation_delta_r_);
      branch(tree, "FormationKt", formation_kt_);
      branch(tree, "FormationParentE", formation_parent_energy_);
      branch(tree, "FormationOffset", formation_offset_);
      branch(tree, "FormationInvalidSplits", formation_invalid_splits_);
      branch(tree, "FormationHardestValid", formation_hardest_valid_);
      branch(tree, "FormationHardestTauF", formation_hardest_tau_f_);
      branch(tree, "FormationHardestTauFSmallAngle", formation_hardest_tau_f_small_angle_);
      branch(tree, "FormationHardestZ", formation_hardest_z_);
      branch(tree, "FormationHardestTheta", formation_hardest_theta_);
      branch(tree, "FormationHardestDeltaR", formation_hardest_delta_r_);
      branch(tree, "FormationHardestKt", formation_hardest_kt_);
      branch(tree, "FormationHardestParentE", formation_hardest_parent_energy_);
    }
  }

  void assign(const std::vector<JetRecord> &records) {
    clear();
    reserve(records.size());
    if (store_formation_time_) {
      formation_offset_.push_back(0);
    }
    for (const JetRecord &record : records) {
      push(eta_, record.eta);
      push(y_, record.y);
      push(phi_, record.phi);
      push(pt_, record.pt);
      push(mass_, record.mass);
      push(energy_, record.energy);
      push(px_, record.px);
      push(py_, record.py);
      push(pz_, record.pz);
      push(raw_eta_, record.raw_eta);
      push(raw_y_, record.raw_y);
      push(raw_phi_, record.raw_phi);
      push(raw_pt_, record.raw_pt);
      push(raw_mass_, record.raw_mass);
      push(zg_, record.zg);
      push(rg_, record.rg);
      push(sd_pt_, record.sd_pt);
      push(sd_mass_, record.sd_mass);
      soft_drop_valid_.push_back(record.soft_drop_valid);
      n_sd_.push_back(record.n_sd);
      mult_.push_back(record.mult);
      total_mult_.push_back(record.total_mult);
      push(pt_d_, record.pt_d);
      push(effective_multiplicity_, record.effective_multiplicity);
      push(leading_fraction_, record.leading_fraction);
      push(normal_pt_d_, record.normal_pt_d);
      push(normal_effective_multiplicity_, record.normal_effective_multiplicity);
      push(leading_normal_fraction_, record.leading_normal_fraction);
      push(girth_, record.girth);
      push(signed_girth_, record.signed_girth);
      push(max_kt_, record.max_kt);
      push(max_kt_z_, record.max_kt_z);
      push(max_kt_rg_, record.max_kt_rg);
      push(max_kt_primary_, record.max_kt_primary);
      push(normal_pt_, record.normal_pt);
      push(positive_wake_pt_, record.positive_wake_pt);
      push(negative_wake_pt_, record.negative_wake_pt);
      push(wake_fraction_, record.wake_fraction);
      n_normal_.push_back(record.n_normal);
      n_positive_wake_.push_back(record.n_positive_wake);
      n_negative_wake_.push_back(record.n_negative_wake);
      n_negative_thermal_.push_back(record.n_negative_thermal);
      n_hadronized_holes_.push_back(record.n_hadronized_holes);
      hard_parton_id_.push_back(record.hard_parton_id);
      push(hard_parton_pt_, record.hard_parton_pt);
      push(hard_parton_dr_, record.hard_parton_dr);
      pair_match_index_.push_back(record.pair_match_index);
      push(pair_match_dr_, record.pair_match_dr);
      push(pair_match_other_pt_, record.pair_match_other_pt);
      pair_match_other_hard_parton_id_.push_back(record.pair_match_other_hard_parton_id);
      if (store_formation_time_) {
        formation_tau_f_.insert(formation_tau_f_.end(), record.formation_tau_f.begin(),
                                record.formation_tau_f.end());
        formation_tau_f_small_angle_.insert(
            formation_tau_f_small_angle_.end(), record.formation_tau_f_small_angle.begin(),
            record.formation_tau_f_small_angle.end());
        formation_z_.insert(formation_z_.end(), record.formation_z.begin(),
                            record.formation_z.end());
        formation_theta_.insert(formation_theta_.end(), record.formation_theta.begin(),
                                record.formation_theta.end());
        formation_delta_r_.insert(formation_delta_r_.end(), record.formation_delta_r.begin(),
                                  record.formation_delta_r.end());
        formation_kt_.insert(formation_kt_.end(), record.formation_kt.begin(),
                             record.formation_kt.end());
        formation_parent_energy_.insert(formation_parent_energy_.end(),
                                        record.formation_parent_energy.begin(),
                                        record.formation_parent_energy.end());
        formation_offset_.push_back(static_cast<int>(formation_tau_f_.size()));
        formation_invalid_splits_.push_back(record.formation_invalid_splits);
        formation_hardest_valid_.push_back(record.formation_hardest_valid);
        push(formation_hardest_tau_f_, record.formation_hardest_tau_f);
        push(formation_hardest_tau_f_small_angle_,
             record.formation_hardest_tau_f_small_angle);
        push(formation_hardest_z_, record.formation_hardest_z);
        push(formation_hardest_theta_, record.formation_hardest_theta);
        push(formation_hardest_delta_r_, record.formation_hardest_delta_r);
        push(formation_hardest_kt_, record.formation_hardest_kt);
        push(formation_hardest_parent_energy_, record.formation_hardest_parent_energy);
      }
    }
  }

 private:
  template <typename T>
  void branch(TTree *tree, const std::string &suffix, std::vector<T> &value) {
    tree->Branch((prefix_ + suffix).c_str(), &value);
  }

  static void push(std::vector<float> &destination, double value) {
    destination.push_back(static_cast<float>(value));
  }

  void reserve(std::size_t size) {
    for (auto *values : float_vectors()) {
      values->reserve(size);
    }
    for (auto *values : int_vectors()) {
      values->reserve(size);
    }
    if (store_formation_time_) {
      formation_offset_.reserve(size + 1);
    }
  }

  void clear() {
    for (auto *values : float_vectors()) {
      values->clear();
    }
    for (auto *values : int_vectors()) {
      values->clear();
    }
    for (auto *values : formation_float_vectors()) {
      values->clear();
    }
    formation_offset_.clear();
  }

  std::vector<std::vector<float> *> float_vectors() {
    return {&eta_,          &y_,                 &phi_,
            &pt_,           &mass_,              &energy_,
            &px_,           &py_,                &pz_,
            &raw_eta_,      &raw_y_,             &raw_phi_,
            &raw_pt_,       &raw_mass_,           &zg_,
            &rg_,           &sd_pt_,              &sd_mass_,
            &pt_d_,         &effective_multiplicity_, &leading_fraction_,
            &normal_pt_d_,  &normal_effective_multiplicity_,
            &leading_normal_fraction_, &girth_,   &signed_girth_,
            &max_kt_,       &max_kt_z_,           &max_kt_rg_,
            &max_kt_primary_, &normal_pt_,        &positive_wake_pt_,
            &negative_wake_pt_, &wake_fraction_,  &pair_match_dr_,
            &pair_match_other_pt_, &hard_parton_pt_, &hard_parton_dr_,
            &formation_hardest_tau_f_, &formation_hardest_tau_f_small_angle_,
            &formation_hardest_z_, &formation_hardest_theta_,
            &formation_hardest_delta_r_, &formation_hardest_kt_,
            &formation_hardest_parent_energy_};
  }

  std::vector<std::vector<int> *> int_vectors() {
    return {&soft_drop_valid_, &n_sd_,              &mult_,
            &total_mult_,
            &n_normal_,        &n_positive_wake_,   &n_negative_wake_,
            &n_negative_thermal_, &n_hadronized_holes_, &hard_parton_id_,
            &pair_match_index_, &pair_match_other_hard_parton_id_,
            &formation_invalid_splits_, &formation_hardest_valid_};
  }

  std::vector<std::vector<float> *> formation_float_vectors() {
    return {&formation_tau_f_, &formation_tau_f_small_angle_, &formation_z_,
            &formation_theta_, &formation_delta_r_, &formation_kt_,
            &formation_parent_energy_};
  }

  std::string prefix_;
  bool store_formation_time_ = false;
  std::vector<float> eta_, y_, phi_, pt_, mass_, energy_, px_, py_, pz_;
  std::vector<float> raw_eta_, raw_y_, raw_phi_, raw_pt_, raw_mass_;
  std::vector<float> zg_, rg_, sd_pt_, sd_mass_;
  std::vector<int> soft_drop_valid_, n_sd_, mult_, total_mult_;
  std::vector<float> pt_d_, effective_multiplicity_, leading_fraction_;
  std::vector<float> normal_pt_d_, normal_effective_multiplicity_;
  std::vector<float> leading_normal_fraction_, girth_, signed_girth_;
  std::vector<float> max_kt_, max_kt_z_, max_kt_rg_, max_kt_primary_;
  std::vector<float> normal_pt_, positive_wake_pt_, negative_wake_pt_, wake_fraction_;
  std::vector<int> n_normal_, n_positive_wake_, n_negative_wake_, n_negative_thermal_;
  std::vector<int> n_hadronized_holes_, hard_parton_id_, pair_match_index_;
  std::vector<int> pair_match_other_hard_parton_id_;
  std::vector<float> hard_parton_pt_, hard_parton_dr_;
  std::vector<float> pair_match_dr_, pair_match_other_pt_;
  std::vector<float> formation_tau_f_, formation_tau_f_small_angle_;
  std::vector<float> formation_z_, formation_theta_, formation_delta_r_;
  std::vector<float> formation_kt_, formation_parent_energy_;
  std::vector<int> formation_offset_;
  std::vector<int> formation_invalid_splits_, formation_hardest_valid_;
  std::vector<float> formation_hardest_tau_f_, formation_hardest_tau_f_small_angle_;
  std::vector<float> formation_hardest_z_, formation_hardest_theta_;
  std::vector<float> formation_hardest_delta_r_, formation_hardest_kt_;
  std::vector<float> formation_hardest_parent_energy_;
};

class JetTree {
 public:
  JetTree(TDirectory *directory, bool prehydro_enabled_in)
      : prehydro_enabled_(prehydro_enabled_in) {
    directory->cd();
    tree_ = new TTree("Jets", "One entry per accepted paired HYBRID event");
    tree_->SetAutoFlush(-10000000);
    branch_event_metadata(tree_, metadata_, prehydro_enabled_);
    tree_->Branch("nJet1", &n_jet1_);
    tree_->Branch("nJet2", &n_jet2_);
    tree_->Branch("nJet4", &n_jet4_);
    tree_->Branch("nJet8", &n_jet8_);
    jet1_ = std::make_unique<RadiusBranches>(tree_, "jet1", false);
    jet2_ = std::make_unique<RadiusBranches>(tree_, "jet2", false);
    jet4_ = std::make_unique<RadiusBranches>(tree_, "jet4", true);
    jet8_ = std::make_unique<RadiusBranches>(tree_, "jet8", true);
  }

  void fill(const EventMeta &metadata, const std::vector<JetRecord> &jets1,
            const std::vector<JetRecord> &jets2,
            const std::vector<JetRecord> &jets4,
            const std::vector<JetRecord> &jets8) {
    metadata_ = metadata;
    n_jet1_ = static_cast<int>(jets1.size());
    n_jet2_ = static_cast<int>(jets2.size());
    n_jet4_ = static_cast<int>(jets4.size());
    n_jet8_ = static_cast<int>(jets8.size());
    jet1_->assign(jets1);
    jet2_->assign(jets2);
    jet4_->assign(jets4);
    jet8_->assign(jets8);
    tree_->Fill();
  }

 private:
  TTree *tree_ = nullptr;
  EventMeta metadata_;
  bool prehydro_enabled_ = false;
  int n_jet1_ = 0;
  int n_jet2_ = 0;
  int n_jet4_ = 0;
  int n_jet8_ = 0;
  std::unique_ptr<RadiusBranches> jet1_;
  std::unique_ptr<RadiusBranches> jet2_;
  std::unique_ptr<RadiusBranches> jet4_;
  std::unique_ptr<RadiusBranches> jet8_;
};

class PairTree {
 public:
  explicit PairTree(TDirectory *directory) {
    directory->cd();
    tree_ = new TTree("Pairs", "One entry per strict no/with-prehydro event pair");
    tree_->SetAutoFlush(-10000000);
    tree_->Branch("pairId", &metadata_.pair_id);
    tree_->Branch("sourceIndex", &metadata_.source_index);
    tree_->Branch("chunkId", &metadata_.chunk_id);
    tree_->Branch("hydroIndex", &metadata_.hydro_index);
    tree_->Branch("seed", &metadata_.seed);
    tree_->Branch("eventNumber", &metadata_.event_number);
    tree_->Branch("eventWeight", &metadata_.event_weight);
    tree_->Branch("sigmaGen", &metadata_.sigma_gen);
    tree_->Branch("hardX", &metadata_.hard_x);
    tree_->Branch("hardY", &metadata_.hard_y);
    tree_->Branch("noPrehydroNHadrons", &no_hadrons_);
    tree_->Branch("withPrehydroNHadrons", &with_hadrons_);
    tree_->Branch("deltaNHadrons", &delta_hadrons_);
    tree_->Branch("noPrehydroSignedSumPt", &no_signed_sum_pt_);
    tree_->Branch("withPrehydroSignedSumPt", &with_signed_sum_pt_);
    tree_->Branch("deltaSignedSumPt", &delta_signed_sum_pt_);
    tree_->Branch("noPrehydroNJet1", &no_jet1_);
    tree_->Branch("withPrehydroNJet1", &with_jet1_);
    tree_->Branch("noPrehydroNJet2", &no_jet2_);
    tree_->Branch("withPrehydroNJet2", &with_jet2_);
    tree_->Branch("noPrehydroNJet4", &no_jet4_);
    tree_->Branch("withPrehydroNJet4", &with_jet4_);
    tree_->Branch("noPrehydroNJet8", &no_jet8_);
    tree_->Branch("withPrehydroNJet8", &with_jet8_);
  }

  void fill(const EventMeta &metadata, const EventSummary &no_summary,
            const EventSummary &with_summary, std::size_t no_jet1, std::size_t with_jet1,
            std::size_t no_jet2, std::size_t with_jet2,
            std::size_t no_jet4, std::size_t with_jet4,
            std::size_t no_jet8, std::size_t with_jet8) {
    metadata_ = metadata;
    no_hadrons_ = no_summary.n_hadrons;
    with_hadrons_ = with_summary.n_hadrons;
    delta_hadrons_ = with_hadrons_ - no_hadrons_;
    no_signed_sum_pt_ = no_summary.signed_sum_pt;
    with_signed_sum_pt_ = with_summary.signed_sum_pt;
    delta_signed_sum_pt_ = with_signed_sum_pt_ - no_signed_sum_pt_;
    no_jet1_ = static_cast<int>(no_jet1);
    with_jet1_ = static_cast<int>(with_jet1);
    no_jet2_ = static_cast<int>(no_jet2);
    with_jet2_ = static_cast<int>(with_jet2);
    no_jet4_ = static_cast<int>(no_jet4);
    with_jet4_ = static_cast<int>(with_jet4);
    no_jet8_ = static_cast<int>(no_jet8);
    with_jet8_ = static_cast<int>(with_jet8);
    tree_->Fill();
  }

 private:
  TTree *tree_ = nullptr;
  EventMeta metadata_;
  int no_hadrons_ = 0;
  int with_hadrons_ = 0;
  int delta_hadrons_ = 0;
  double no_signed_sum_pt_ = 0.0;
  double with_signed_sum_pt_ = 0.0;
  double delta_signed_sum_pt_ = 0.0;
  int no_jet1_ = 0;
  int with_jet1_ = 0;
  int no_jet2_ = 0;
  int with_jet2_ = 0;
  int no_jet4_ = 0;
  int with_jet4_ = 0;
  int no_jet8_ = 0;
  int with_jet8_ = 0;
};

EventMeta metadata_from_header(const PairHeader &header) {
  EventMeta metadata;
  metadata.pair_id = header.pair_id;
  metadata.source_index = header.source_index;
  metadata.chunk_id = header.chunk_id;
  metadata.hydro_index = header.hydro_index;
  metadata.seed = header.seed;
  metadata.event_number = header.event_number;
  metadata.event_weight = header.event_weight;
  metadata.sigma_gen = header.sigma_gen;
  metadata.hard_x = header.hard_x;
  metadata.hard_y = header.hard_y;
  return metadata;
}

void write_metadata(TFile &output, const Options &options, std::uint64_t pair_count,
                    std::uint64_t no_hadron_count, std::uint64_t with_hadron_count,
                    std::uint64_t no_jet1_count, std::uint64_t with_jet1_count,
                    std::uint64_t no_jet2_count, std::uint64_t with_jet2_count,
                    std::uint64_t no_jet4_count, std::uint64_t with_jet4_count,
                    std::uint64_t no_jet8_count, std::uint64_t with_jet8_count) {
  TDirectory *directory = output.mkdir("metadata");
  directory->cd();
  TNamed schema("schemaVersion", kSchemaVersion);
  schema.Write();
  TNamed status_mapping(
      "hadronStatusDefinition",
      "0=normal(raw 0); +1=positive wake(raw 1); -1=negative wake/hole(raw 2 or 3); raw -2 excluded");
  status_mapping.Write();
  TNamed jet_definition(
      "jetDefinition",
      "anti-kt E-scheme; normal+positive-wake clustering; negative wake/hole ghost association; negative four-vector subtraction (4MomSub)");
  jet_definition.Write();
  TNamed substructure_definition(
      "substructureDefinition",
      "positive constituents only; Cambridge/Aachen reclustering with R_CA=2*R_antiKt+1e-6, E-scheme, Best strategy; Soft Drop first passing hardest-branch split with R0=R_antiKt; MaxKt over full C/A tree; Mult=Nnormal+Npositive; TotalMult=Nnormal+Npositive-Nnegative; normal-only effective multiplicity excludes wake hadrons");
  substructure_definition.Write();
  TNamed formation_time_definition(
      "formationTimeDefinition",
      "R=0.4,0.8 full C/A trees with R_CA=2*R_antiKt+1e-6, E-scheme, Best strategy; tau_f=2 hbarc Eparent/Qparent^2=hbarc/[Eparent z1 z2 (1-cos(theta12))], hbarc=0.19732698 GeV fm, zi=Ei/Eparent, exact three-dimensional theta12; small-angle audit replaces 1-cos(theta12) by DeltaR12^2/2; hardest split maximizes min(pT1,pT2)*DeltaR12 over the full tree");
  formation_time_definition.Write();
  TNamed formation_time_interpretation(
      "formationTimeInterpretation",
      "final-hadron Cambridge/Aachen formation-time estimator; not generator-level parton-shower history; negative wake and hadronized holes affect corrected jet selection through 4MomSub but are excluded from the nonlinear constituent tree");
  formation_time_interpretation.Write();
  TNamed hard_parton_definition(
      "hardPartonMatchDefinition",
      "one-to-one nearest-axis matching of raw-label -2 outgoing hard-parton markers to raw jet axes within DeltaR<R");
  hard_parton_definition.Write();
  TParameter<double>("rawJetPtMin", options.raw_jet_pt_min).Write();
  TParameter<double>("jetAbsEtaMax", options.jet_abs_eta_max).Write();
  TParameter<double>("softDropZCut", options.z_cut).Write();
  TParameter<double>("softDropBeta", options.beta).Write();
  TParameter<double>("pairMatchDRFraction", options.match_dr_fraction).Write();
  TParameter<double>("formationTimeCorrectedPtMin", kFormationTimeCorrectedPtMin).Write();

  int source_index = 0;
  std::string source_name;
  std::string source_path;
  TTree sources("Sources", "Input source index mapping");
  sources.Branch("sourceIndex", &source_index);
  sources.Branch("sourceName", &source_name);
  sources.Branch("sourcePath", &source_path);
  for (std::size_t index = 0; index < options.sources.size(); ++index) {
    source_index = static_cast<int>(index);
    source_name = options.sources[index].first;
    source_path = options.sources[index].second;
    sources.Fill();
  }
  sources.Write();

  TTree conversion("Conversion", "Conversion totals");
  conversion.Branch("pairCount", &pair_count);
  conversion.Branch("noPrehydroHadronCount", &no_hadron_count);
  conversion.Branch("withPrehydroHadronCount", &with_hadron_count);
  conversion.Branch("noPrehydroJet1Count", &no_jet1_count);
  conversion.Branch("withPrehydroJet1Count", &with_jet1_count);
  conversion.Branch("noPrehydroJet2Count", &no_jet2_count);
  conversion.Branch("withPrehydroJet2Count", &with_jet2_count);
  conversion.Branch("noPrehydroJet4Count", &no_jet4_count);
  conversion.Branch("withPrehydroJet4Count", &with_jet4_count);
  conversion.Branch("noPrehydroJet8Count", &no_jet8_count);
  conversion.Branch("withPrehydroJet8Count", &with_jet8_count);
  conversion.Fill();
  conversion.Write();
}

int run(const Options &options) {
  fastjet::ClusterSequence::set_fastjet_banner_stream(nullptr);
  TFile output(options.output.c_str(), "RECREATE", "OO paired event and jet trees");
  if (output.IsZombie()) {
    throw std::runtime_error("could not create ROOT output " + options.output);
  }
  output.SetCompressionAlgorithm(ROOT::RCompressionSetting::EAlgorithm::kZSTD);
  output.SetCompressionLevel(5);

  TDirectory *no_directory = output.mkdir("noPrehydro");
  TDirectory *with_directory = output.mkdir("withPrehydro");
  HadronTree no_hadron_tree(no_directory, false);
  JetTree no_jet_tree(no_directory, false);
  HadronTree with_hadron_tree(with_directory, true);
  JetTree with_jet_tree(with_directory, true);
  PairTree pair_tree(&output);

  std::uint64_t pair_count = 0;
  std::uint64_t no_hadron_count = 0;
  std::uint64_t with_hadron_count = 0;
  std::uint64_t no_jet1_count = 0;
  std::uint64_t with_jet1_count = 0;
  std::uint64_t no_jet2_count = 0;
  std::uint64_t with_jet2_count = 0;
  std::uint64_t no_jet4_count = 0;
  std::uint64_t with_jet4_count = 0;
  std::uint64_t no_jet8_count = 0;
  std::uint64_t with_jet8_count = 0;

  while (true) {
    PairHeader header{};
    if (!read_record(std::cin, header, true)) {
      break;
    }
    if (!std::equal(kPairMagic.begin(), kPairMagic.end(), header.magic)) {
      throw std::runtime_error("invalid event-stream magic");
    }
    constexpr std::uint32_t kMaximumParticlesPerVariant = 10000000;
    if (header.no_prehydro_count > kMaximumParticlesPerVariant ||
        header.with_prehydro_count > kMaximumParticlesPerVariant) {
      throw std::runtime_error("unreasonable particle count in event stream");
    }
    std::vector<ParticleRecord> no_particles(header.no_prehydro_count);
    std::vector<ParticleRecord> with_particles(header.with_prehydro_count);
    for (ParticleRecord &particle : no_particles) {
      read_record(std::cin, particle);
    }
    for (ParticleRecord &particle : with_particles) {
      read_record(std::cin, particle);
    }

    const EventMeta metadata = metadata_from_header(header);
    const EventSummary no_summary = no_hadron_tree.fill(metadata, no_particles);
    const EventSummary with_summary = with_hadron_tree.fill(metadata, with_particles);
    auto no_jets1 = make_jets(no_particles, 0.1, options);
    auto with_jets1 = make_jets(with_particles, 0.1, options);
    auto no_jets2 = make_jets(no_particles, 0.2, options);
    auto with_jets2 = make_jets(with_particles, 0.2, options);
    auto no_jets4 = make_jets(no_particles, 0.4, options);
    auto with_jets4 = make_jets(with_particles, 0.4, options);
    auto no_jets8 = make_jets(no_particles, 0.8, options);
    auto with_jets8 = make_jets(with_particles, 0.8, options);
    match_hard_partons(no_jets1, no_particles, 0.1);
    match_hard_partons(with_jets1, with_particles, 0.1);
    match_hard_partons(no_jets2, no_particles, 0.2);
    match_hard_partons(with_jets2, with_particles, 0.2);
    match_hard_partons(no_jets4, no_particles, 0.4);
    match_hard_partons(with_jets4, with_particles, 0.4);
    match_hard_partons(no_jets8, no_particles, 0.8);
    match_hard_partons(with_jets8, with_particles, 0.8);
    match_jets(no_jets1, with_jets1, 0.1, options.match_dr_fraction);
    match_jets(no_jets2, with_jets2, 0.2, options.match_dr_fraction);
    match_jets(no_jets4, with_jets4, 0.4, options.match_dr_fraction);
    match_jets(no_jets8, with_jets8, 0.8, options.match_dr_fraction);
    no_jet_tree.fill(metadata, no_jets1, no_jets2, no_jets4, no_jets8);
    with_jet_tree.fill(metadata, with_jets1, with_jets2, with_jets4, with_jets8);
    pair_tree.fill(metadata, no_summary, with_summary, no_jets1.size(), with_jets1.size(),
                   no_jets2.size(), with_jets2.size(),
                   no_jets4.size(), with_jets4.size(),
                   no_jets8.size(), with_jets8.size());

    ++pair_count;
    no_hadron_count += no_summary.n_hadrons;
    with_hadron_count += with_summary.n_hadrons;
    no_jet1_count += no_jets1.size();
    with_jet1_count += with_jets1.size();
    no_jet2_count += no_jets2.size();
    with_jet2_count += with_jets2.size();
    no_jet4_count += no_jets4.size();
    with_jet4_count += with_jets4.size();
    no_jet8_count += no_jets8.size();
    with_jet8_count += with_jets8.size();
  }

  write_metadata(output, options, pair_count, no_hadron_count, with_hadron_count,
                 no_jet1_count, with_jet1_count,
                 no_jet2_count, with_jet2_count, no_jet4_count, with_jet4_count,
                 no_jet8_count, with_jet8_count);
  output.Write();
  output.Close();
  std::cout << "pairs=" << pair_count << " no_hadrons=" << no_hadron_count
            << " with_hadrons=" << with_hadron_count << " no_jet1=" << no_jet1_count
            << " with_jet1=" << with_jet1_count << " no_jet2=" << no_jet2_count
            << " with_jet2=" << with_jet2_count << " no_jet4=" << no_jet4_count
            << " with_jet4=" << with_jet4_count << " no_jet8=" << no_jet8_count
            << " with_jet8=" << with_jet8_count << '\n';
  return 0;
}

}  // namespace

int main(int argc, char **argv) {
  try {
    return run(parse_options(argc, argv));
  } catch (const std::exception &error) {
    std::cerr << "oo_root_tree_writer: " << error.what() << '\n';
    return 1;
  }
}
