#include "HYBRID.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <sstream>
#include <stdexcept>
#include "vector_operators.h"
#include "ResponseLedger.h"

namespace {
constexpr int kShowerSeedOffset = 33;
constexpr int kHybridSeedOffset = 1346;
constexpr int kLundSeedOffset = 2337;
constexpr int kElasticSeedOffset = 3779;

int getSeedBase(const Config &cfg) {
    return cfg.getIntOr("seed_base", cfg.getIntOr("njob", 0));
}

const Config &validateConfig(const Config &cfg) {
    cfg.validateKnownKeys({
        "njob", "seed_base", "Nev", "cent", "kappa", "alpha", "tmethod", "mode",
        "do_quench", "do_wake", "do_source", "do_elastic", "do_lres",
        "do_Moliere_on_unresolved_partons", "do_Moliere_dynamic_unresolved_resolution",
        "do_Moliere_dynamic_daughter_unresolved_resolution",
        "do_Moliere_recursive_unresolved_resolution",
        "allow_modee_nonconserving_opening", "modee_max_opening_relative_residual",
        "moliere_unresolved_resolution_c", "lres_rpower", "rpower",
        "dump_hybrid_evolution_history", "hybrid_evolution_history_file",
        "doEventDisplay", "eventDisplayFile", "compat_moliere_legacy_hydro",
        "use_prehydro", "prehydro_file", "use_fixed_xy", "fixed_x", "fixed_y",
        "ebe_hydro", "hadro_type", "tables_path", "output_base",
        "use_trigger", "trigger_id", "trigger_pt", "trigger_eta",
        "pythia_cmnd", "max_tree_attempts", "max_event_attempts"});

    const int Nev = cfg.getIntOr("Nev", 1);
    const int max_tree_attempts = cfg.getIntOr("max_tree_attempts", 100000);
    const int max_event_attempts = cfg.getIntOr("max_event_attempts", 1000000);
    const int tmethod = cfg.getIntOr("tmethod", 0);
    const int mode = cfg.getIntOr("mode", 0);
    const int ebe_hydro = cfg.getIntOr("ebe_hydro", 0);
    const bool do_elastic = cfg.getBoolOr("do_elastic", false);
    const bool do_lres = cfg.getBoolOr("do_lres", false);
    const int hadro_type = cfg.getIntOr("hadro_type", do_elastic ? 1 : 0);
    if (Nev <= 0 || max_tree_attempts <= 0 || max_event_attempts <= 0) {
        throw std::invalid_argument("Nev and generation-attempt limits must be positive");
    }
    if (tmethod < 0 || tmethod > 1 || mode < 0 || mode > 2 ||
        ebe_hydro < 0 || ebe_hydro > 1 || hadro_type < 0 || hadro_type > 1) {
        throw std::invalid_argument("Invalid tmethod, mode, ebe_hydro, or hadro_type value");
    }
    const double lres_rpower = cfg.getDoubleOr(
        "lres_rpower", cfg.getDoubleOr("rpower", 2.0));
    if (cfg.getDoubleOr("kappa", 1.0) < 0. || cfg.getDoubleOr("alpha", 1.0) < 0. ||
        lres_rpower < 0. || cfg.getDoubleOr("moliere_unresolved_resolution_c", 1.0) < 0. ||
        cfg.getDoubleOr("modee_max_opening_relative_residual", 1.e-6) < 0.) {
        throw std::invalid_argument("Physics parameters and residual thresholds must be nonnegative");
    }

    const bool mode_b = cfg.getBoolOr("do_Moliere_on_unresolved_partons", false);
    const bool mode_c = cfg.getBoolOr("do_Moliere_dynamic_unresolved_resolution", false);
    const bool mode_d = cfg.getBoolOr("do_Moliere_dynamic_daughter_unresolved_resolution", false);
    const bool mode_e = cfg.getBoolOr("do_Moliere_recursive_unresolved_resolution", false);
    const bool allow_modee = cfg.getBoolOr("allow_modee_nonconserving_opening", false);
    const int unresolved_mode_count = static_cast<int>(mode_b) + static_cast<int>(mode_c) +
                                      static_cast<int>(mode_d) + static_cast<int>(mode_e);
    if (unresolved_mode_count > 1) {
        throw std::invalid_argument("Select at most one unresolved Moliere mode (B-E)");
    }
    if (unresolved_mode_count > 0 && (!do_elastic || !do_lres)) {
        throw std::invalid_argument("Unresolved Moliere modes require do_elastic=true and do_lres=true");
    }
    if (allow_modee && !mode_e) {
        throw std::invalid_argument("allow_modee_nonconserving_opening is only valid for Mode E");
    }
    if (mode_e && !allow_modee) {
        throw std::invalid_argument(
            "Mode E is validation-only: its live-parent opening is not four-momentum "
            "conserving. Set allow_modee_nonconserving_opening=true only for diagnostics.");
    }
    return cfg;
}
}

HYBRID::HYBRID(const Config &cfg) :
      do_quench_(validateConfig(cfg).getBoolOr("do_quench", true)),
      do_wake_(cfg.getBoolOr("do_wake", true)),
      do_source_(cfg.getBoolOr("do_source", false)),
      do_elastic_(cfg.getBoolOr("do_elastic", false)),
      do_lres_(cfg.getBoolOr("do_lres", false)),
      do_moliere_on_unresolved_partons_(cfg.getBoolOr("do_Moliere_on_unresolved_partons", false)),
      do_moliere_dynamic_unresolved_resolution_(cfg.getBoolOr("do_Moliere_dynamic_unresolved_resolution", false)),
      do_moliere_dynamic_daughter_unresolved_resolution_(cfg.getBoolOr("do_Moliere_dynamic_daughter_unresolved_resolution", false)),
      do_moliere_recursive_unresolved_resolution_(cfg.getBoolOr("do_Moliere_recursive_unresolved_resolution", false)),
      allow_modee_nonconserving_opening_(cfg.getBoolOr("allow_modee_nonconserving_opening", false)),
      dump_hybrid_evolution_history_(cfg.getBoolOr("dump_hybrid_evolution_history", false)),
      do_event_display_(cfg.getBoolOr("doEventDisplay", false)),
      use_fixed_xy_(cfg.getBoolOr("use_fixed_xy", false)),
      compat_moliere_legacy_hydro_(cfg.getBoolOr("compat_moliere_legacy_hydro", true)),
      use_prehydro_(cfg.getBoolOr("use_prehydro", false)),
      njob_(cfg.getIntOr("njob", 0)),
      Nev_(cfg.getIntOr("Nev", 1)),
      cent_(cfg.getStringOr("cent", "0-5")),
      kappa_(cfg.getDoubleOr("kappa", 1.0)),
      alpha_(cfg.getDoubleOr("alpha", 1.0)),
      tmethod_(cfg.getIntOr("tmethod", 0)),
      mode_(cfg.getIntOr("mode", 0)),
      ebe_hydro_(cfg.getIntOr("ebe_hydro", 0)),
      hadro_type_(cfg.getIntOr("hadro_type", cfg.getBoolOr("do_elastic", false) ? 1 : 0)),
      max_tree_attempts_(cfg.getIntOr("max_tree_attempts", 100000)),
      max_event_attempts_(cfg.getIntOr("max_event_attempts", 1000000)),
      lres_rpower_(cfg.getDoubleOr("lres_rpower", cfg.getDoubleOr("rpower", 2.0))),
      moliere_unresolved_resolution_c_(cfg.getDoubleOr("moliere_unresolved_resolution_c", 1.0)),
      modee_max_opening_relative_residual_(cfg.getDoubleOr("modee_max_opening_relative_residual", 1.e-6)),
      seed_base_(getSeedBase(cfg)),
      shower_seed_(seed_base_ + kShowerSeedOffset),
      hybrid_seed_(seed_base_ + kHybridSeedOffset),
      lund_seed_(seed_base_ + kLundSeedOffset),
      elastic_seed_(seed_base_ + kElasticSeedOffset),
      fixed_x_(cfg.getDoubleOr("fixed_x", 0.0)),
      fixed_y_(cfg.getDoubleOr("fixed_y", 0.0)),
      tables_path_(cfg.getStringOr("tables_path", "")),
      prehydro_file_(cfg.getStringOr("prehydro_file", "prehydro_table.tsv")),
      hybrid_evolution_history_file_(cfg.getStringOr("hybrid_evolution_history_file", "")),
      event_display_file_(cfg.getStringOr("eventDisplayFile", "eventDisplay.root")),
      pythia_cmnd_(cfg.getStringOr("pythia_cmnd", "setup_pythia.cmnd")),
      generated_event_attempts_(0),
      tree_failures_(0),
      hadronization_failures_(0),
      nr_(hybrid_seed_),
      tree_gen_(std::make_unique<TreeGenerator>()),
      hydro_profile_(std::make_unique<HydroProfile>()),
      wake_gen_(std::make_unique<WakeGenerator>()),
      lund_gen_(std::make_unique<LundGenerator>()),
      glauber_model_(std::make_unique<GlauberModel>()),
      energy_loss_(std::make_unique<EnergyLoss>(nr_, kappa_, alpha_, tmethod_, mode_,
                                                ebe_hydro_, do_elastic_, do_lres_,
                                                do_moliere_on_unresolved_partons_,
                                                do_moliere_dynamic_unresolved_resolution_,
                                                do_moliere_dynamic_daughter_unresolved_resolution_,
                                                do_moliere_recursive_unresolved_resolution_,
                                                modee_max_opening_relative_residual_,
                                                moliere_unresolved_resolution_c_,
                                                lres_rpower_,
                                                dump_hybrid_evolution_history_,
                                                hybrid_evolution_history_file_,
                                                do_event_display_,
                                                event_display_file_,
                                                compat_moliere_legacy_hydro_,
                                                elastic_seed_,
                                                tables_path_,
                                                *hydro_profile_)) {
    // Open output files
    const std::string out_base = cfg.getStringOr("output_base", "HYBRID");
    // Match the legacy executables: each run owns a fresh output file.
    // Appending across reruns in the same directory can duplicate event blocks.
    hjt_file_.open(out_base + "_Hadrons.out", std::ios_base::out | std::ios_base::trunc);
    pjt_file_.open(out_base + "_Partons.out", std::ios_base::out | std::ios_base::trunc);
    if (!hjt_file_ || !pjt_file_) {
        throw std::runtime_error("Failed to open output files for base: " + out_base);
    }

    std::cout << "Seed base= " << seed_base_
              << " shower= " << shower_seed_
              << " hybrid= " << hybrid_seed_
              << " lund= " << lund_seed_
              << " elastic= " << elastic_seed_ << std::endl;
    if (do_quench_) {
        std::cout << "MMLI heavy-flavor transport is disabled; use MMLHI for charm or bottom physics"
                  << std::endl;
    }
    if (do_elastic_) {
        std::cout << "Elastic scattering requested"
                  << " hadro_type= " << hadro_type_
                  << " compat_moliere_legacy_hydro= " << compat_moliere_legacy_hydro_;
        if (!tables_path_.empty()) {
            std::cout << " tables_path= " << tables_path_;
        }
        std::cout << std::endl;
    }
    if (do_lres_) {
        std::cout << "Finite LRES requested"
                  << " rpower= " << lres_rpower_
                  << " (do_elastic= " << do_elastic_ << ")"
                  << " (do_Moliere_on_unresolved_partons= "
                  << do_moliere_on_unresolved_partons_ << ")"
                  << " (do_Moliere_dynamic_unresolved_resolution= "
                  << do_moliere_dynamic_unresolved_resolution_ << ")"
                  << " (do_Moliere_dynamic_daughter_unresolved_resolution= "
                  << do_moliere_dynamic_daughter_unresolved_resolution_ << ")"
                  << " (do_Moliere_recursive_unresolved_resolution= "
                  << do_moliere_recursive_unresolved_resolution_ << ")"
                  << " (moliere_unresolved_resolution_c= "
                  << moliere_unresolved_resolution_c_ << ")"
                  << std::endl;
    }
    if (dump_hybrid_evolution_history_) {
        std::cout << "Hybrid evolution history dump enabled";
        if (!hybrid_evolution_history_file_.empty()) {
            std::cout << " file= " << hybrid_evolution_history_file_;
        }
        std::cout << std::endl;
    }
    if (do_event_display_) {
        std::cout << "Event-display ROOT segment dump enabled file= "
                  << event_display_file_ << std::endl;
    }
    if (use_fixed_xy_) {
        std::cout << "Using fixed production vertex"
                  << " X= " << fixed_x_
                  << " Y= " << fixed_y_
                  << std::endl;
    }
    if (use_prehydro_) {
        std::cout << "arXiv:2509.19430 pre-hydro enabled file= " << prehydro_file_ << std::endl;
    }

    if (cfg.getBoolOr("use_trigger", false)) {
        tree_gen_->setTrigger(
            cfg.getDoubleOr("trigger_pt", 0.0),
            cfg.getDoubleOr("trigger_eta", 2.4),
            cfg.getIntOr("trigger_id", 22));
    }

    if (do_quench_) {
        if (ebe_hydro_ == 0) read_nuclear();
        else Ncollsize_ = read_nuclear_ipsat();
        read_hydro();
    }
}

HYBRID::~HYBRID() {
    std::cout << "HYBRID run diagnostics: generated_event_attempts="
              << generated_event_attempts_ << " tree_failures=" << tree_failures_
              << " hadronization_failures=" << hadronization_failures_ << std::endl;
    hjt_file_.close();
    pjt_file_.close();
}

void HYBRID::run() {
    int count = 0;
    bool lund_initialized = false;
    bool tree_initialized = false;

    std::vector<Parton> partons;
    std::vector<Quench> quenched;
    std::vector<Quench> recoiled;
    std::vector<Wake> wake;
    std::vector<Hadron> vhadrons;
    std::vector<Hadron> qhadrons;
    std::vector<Hadron> hhadrons;

    while (count < Nev_) {
        // Generate PYTHIA tree
        partons.clear();
        if (!tree_initialized) {
            init_tree();
            tree_initialized = true;
        }
        double weight = 0.;
        double cross = 0.;
        double cross_err = 0.;
        int tree_attempts = 0;
        while (!do_tree(partons, weight, cross, cross_err)) {
            ++tree_failures_;
            if (++tree_attempts >= max_tree_attempts_) {
                throw std::runtime_error("Exceeded max_tree_attempts while generating or triggering an event");
            }
        }
        ++generated_event_attempts_;
        if (generated_event_attempts_ > max_event_attempts_) {
            throw std::runtime_error("Exceeded max_event_attempts after repeated event rejection");
        }
        // Create vector of quenched partons initially equal to vacuum partons
        quenched.clear();
        quenched.reserve(partons.size());
        for (const auto &parton : partons) {
            quenched.emplace_back(parton);
        }
        recoiled.clear();

        double x = 0., y = 0.;
        if (do_quench_) {
            // Generate x,y
            if (use_fixed_xy_) {
                x = fixed_x_;
                y = fixed_y_;
            } else if (ebe_hydro_ == 0) {
                gxy(x, y);
            } else {
                gxy_ipsat(x, y);
            }
            std::cout << " xcre= " << x << " ycre= " << y << std::endl;

            std::ofstream source_file;
            if (do_source_) {
                source_file.open("SOURCE.dat", std::ios::app);
                source_file << "# event " << count << "\n";
                source_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << "\n";
                for (const auto &p : partons) {
                    if (p.GetOrig() == "hs") {
                        auto pp = p.vGetP();
                        source_file << "p_x = " << pp[0] << " p_y= " << pp[1] << " p_z = " << pp[2] << " p_e = " << pp[3] << "\n";
                    }
                }
            }

            do_eloss(partons, quenched, recoiled, x, y);

            if (do_source_) {
                source_file << "end\n";
            }
        }

        if (do_wake_) {
            // Do back-reaction
            wake.clear();
            do_wake(quenched, partons, recoiled, wake);
            std::cout << "Wake size= " << wake.size() << std::endl;
        }

        // Hadronize
        if (!lund_initialized) {
            init_lund();
            lund_initialized = true;
        }

        vhadrons.clear();
        qhadrons.clear();
        hhadrons.clear();
        std::string hadronization_failure_reason;
        if (!do_lund(partons, quenched, recoiled, vhadrons, qhadrons, hhadrons,
                     hadronization_failure_reason)) {
            ++hadronization_failures_;
            std::cout << "EVENT_REJECT generated_event=" << generated_event_attempts_
                      << " output_slot=" << count
                      << " reason=" << hadronization_failure_reason << std::endl;
            continue;
        }
        std::cout << " Vac Hadron size= " << vhadrons.size() << " Med Hadron size= " << qhadrons.size() << std::endl;

        // 4-momentum conservation check
        {
            std::vector<double> sum_vac_partons(4, 0.);
            std::vector<double> sum_med_partons(4, 0.);
            std::vector<double> sum_vac_hadrons(4, 0.);
            std::vector<double> sum_med_hadrons(4, 0.);
            std::vector<double> sum_med_hadrons_wake(4, 0.);

            for (const auto &p : partons) {
                if (p.GetD1() == -1) {
                    auto pp = p.vGetP();
                    for (int j = 0; j < 4; j++) sum_vac_partons[j] += pp[j];
                }
            }
            for (const auto &q : quenched) {
                if (q.GetD1() == -1 && q.vGetP()[3] != 0.) {
                    auto qp = q.vGetP();
                    for (int j = 0; j < 4; j++) sum_med_partons[j] += qp[j];
                }
            }
            for (const auto &q : recoiled) {
                if (q.GetD1() != -1 || q.vGetP()[3] == 0.) continue;
                auto qp = q.vGetP();
                if (q.GetOrig() == "recoiler") {
                    for (int j = 0; j < 4; j++) sum_med_partons[j] += qp[j];
                } else if (q.GetOrig() == "hole") {
                    for (int j = 0; j < 4; j++) sum_med_partons[j] -= qp[j];
                }
            }
            for (const auto &h : vhadrons) {
                auto hp = h.vGetP();
                for (int j = 0; j < 4; j++) sum_vac_hadrons[j] += hp[j];
            }
            for (const auto &h : qhadrons) {
                auto hp = h.vGetP();
                for (int j = 0; j < 4; j++) sum_med_hadrons[j] += hp[j];
            }
            for (const auto &h : hhadrons) {
                auto hp = h.vGetP();
                for (int j = 0; j < 4; j++) sum_med_hadrons[j] -= hp[j];
            }
            sum_med_hadrons_wake = sum_med_hadrons;
            for (const auto &w : wake) {
                auto wp = w.vGetP();
                double stat = w.GetStatus();
                for (int j = 0; j < 4; j++) sum_med_hadrons_wake[j] += wp[j] * stat;
            }

            std::cout << " === 4-momentum conservation (event " << count << ") ===" << std::endl;
            std::cout << "  Vac partons  (px,py,pz,E): "
                      << sum_vac_partons[0] << " " << sum_vac_partons[1] << " "
                      << sum_vac_partons[2] << " " << sum_vac_partons[3] << std::endl;
            std::cout << "  Med partons  (px,py,pz,E): "
                      << sum_med_partons[0] << " " << sum_med_partons[1] << " "
                      << sum_med_partons[2] << " " << sum_med_partons[3] << std::endl;
            std::cout << "  Vac hadrons  (px,py,pz,E): "
                      << sum_vac_hadrons[0] << " " << sum_vac_hadrons[1] << " "
                      << sum_vac_hadrons[2] << " " << sum_vac_hadrons[3] << std::endl;
            std::cout << "  Med hadrons  (px,py,pz,E): "
                      << sum_med_hadrons[0] << " " << sum_med_hadrons[1] << " "
                      << sum_med_hadrons[2] << " " << sum_med_hadrons[3] << std::endl;
            std::cout << "  Med+wake     (px,py,pz,E): "
                      << sum_med_hadrons_wake[0] << " " << sum_med_hadrons_wake[1] << " "
                      << sum_med_hadrons_wake[2] << " " << sum_med_hadrons_wake[3] << std::endl;
        }

        // Output
        output_event(count, partons, quenched, recoiled, vhadrons, qhadrons, hhadrons, wake, weight, cross, x, y);

        ++count;
    }
}

// Implement the private methods by adapting from main.cc
// For brevity, I'll sketch them; in practice, copy the implementations

void HYBRID::read_nuclear() {
    // Use GlauberModel method
    glauber_model_->readNuclear(0, cent_);
}

int HYBRID::read_nuclear_ipsat() {
    return glauber_model_->readNuclearIPSAT(0, cent_);
}

void HYBRID::read_hydro() {
    hydro_profile_->loadHydro(ebe_hydro_, cent_);
    if (use_prehydro_) {
        hydro_profile_->loadPreHydroTable(prehydro_file_);
    }
}

void HYBRID::init_tree() {
    tree_gen_->init(shower_seed_, pythia_cmnd_);
}

bool HYBRID::do_tree(std::vector<Parton> &partons, double &weight, double &cross, double &cross_err) {
    return tree_gen_->nextEvent(partons, weight, cross, cross_err);
}

void HYBRID::gxy(double &x, double &y) {
    glauber_model_->sampleXY(x, y, nr_);
}

void HYBRID::gxy_ipsat(double &x, double &y) {
    glauber_model_->sampleXYIPSAT(x, y, Ncollsize_, nr_);
}

void HYBRID::do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched,
                      std::vector<Quench> &recoiled, double x, double y) {
    energy_loss_->do_eloss(partons, quenched, x, y, &recoiled);
}

void HYBRID::do_wake(const std::vector<Quench> &quenched, const std::vector<Parton> &partons,
                       const std::vector<Quench> &recoiled, std::vector<Wake> &wake) {
    const auto deposits = build_response_ledger(quenched, partons, recoiled);
    wake_gen_->generate(deposits, wake, nr_);
}

void HYBRID::init_lund() {
    lund_gen_->init(lund_seed_);
}

bool HYBRID::do_lund(const std::vector<Parton> &partons,
                     const std::vector<Quench> &quenched,
                     const std::vector<Quench> &recoiled,
                     std::vector<Hadron> &vhadrons,
                     std::vector<Hadron> &qhadrons,
                     std::vector<Hadron> &hhadrons,
                     std::string &failure_reason) {
    lund_gen_->hadronizeVacuum(partons, vhadrons);
    std::vector<Quench> quenchandrecoil = quenched;
    std::vector<Quench> holes;
    if (do_elastic_) {
        for (const auto &q : recoiled) {
            if (q.GetOrig() == "recoiler") {
                quenchandrecoil.push_back(q);
            } else if (q.GetOrig() == "hole") {
                holes.push_back(q);
            }
        }
    }
    int had_counter = 0;
    constexpr int had_counter_max = 5;
    bool had_is_ok = false;
    do {
        qhadrons.clear();
        had_is_ok = lund_gen_->hadronizeMedium(quenchandrecoil, qhadrons, hadro_type_);
        had_counter += 1;
    } while (!had_is_ok && had_counter < had_counter_max);
    if (had_counter > 1) {
        std::cout << "Had Counter = " << had_counter << " and had_is_ok= " << had_is_ok << std::endl;
    }
    if (!had_is_ok) {
        failure_reason = "medium_hadronization_failed";
        return false;
    }
    if (!holes.empty()) {
        hhadrons.clear();
        if (!lund_gen_->hadronizeMedium(holes, hhadrons, hadro_type_)) {
            qhadrons.clear();
            hhadrons.clear();
            failure_reason = "hole_hadronization_failed";
            return false;
        }
    }
    failure_reason.clear();
    return true;
}

void HYBRID::output_event(int count,
                           const std::vector<Parton> &partons,
                           const std::vector<Quench> &quenched,
                           const std::vector<Quench> &recoiled,
                           const std::vector<Hadron> &vhadrons,
                           const std::vector<Hadron> &qhadrons,
                           const std::vector<Hadron> &hhadrons,
                           const std::vector<Wake> &wake,
                           double weight,
                           double cross,
                           double x,
                           double y) {
    // Partonic output
    pjt_file_ << "# event " << count << std::endl;
    pjt_file_ << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << std::endl;
    for (const auto &p : partons) {
        if (p.GetOrig() == "hs") {
            auto pp = p.vGetP();
            pjt_file_ << pp[0] << " " << pp[1] << " " << pp[2] << " " << p.GetMass() << " " << p.GetId() << " " << -2 << std::endl;
        }
    }

    if (do_quench_) {
        for (const auto &q : quenched) {
            if (q.GetD1() == -1 && q.vGetP()[3] != 0.) {
                auto qp = q.vGetP();
                pjt_file_ << qp[0] << " " << qp[1] << " " << qp[2] << " " << q.GetMass() << " " << q.GetId() << " " << q.hadScattering() << std::endl;
            }
        }
        for (const auto &q : recoiled) {
            if (q.GetD1() == -1 && q.vGetP()[3] != 0.) {
                int recid = -1000;
                if (q.GetOrig() == "recoiler") recid = 3;
                else if (q.GetOrig() == "hole") recid = 4;
                else continue;
                auto qp = q.vGetP();
                pjt_file_ << qp[0] << " " << qp[1] << " " << qp[2] << " " << q.GetMass() << " " << q.GetId() << " " << recid << std::endl;
            }
        }
    } else {
        for (const auto &p : partons) {
            if (p.GetD1() == -1) {
                auto pp = p.vGetP();
                pjt_file_ << pp[0] << " " << pp[1] << " " << pp[2] << " " << p.GetMass() << " " << p.GetId() << " " << 0 << std::endl;
            }
        }
    }
    pjt_file_ << "end" << std::endl;

    // Hadronic output
    hjt_file_ << "# event " << count << std::endl;
    hjt_file_ << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << std::endl;
    for (const auto &p : partons) {
        if (p.GetOrig() == "hs") {
            auto pp = p.vGetP();
            hjt_file_ << pp[0] << " " << pp[1] << " " << pp[2] << " " << p.GetMass() << " " << p.GetId() << " " << -2 << std::endl;
        }
    }

    if (do_quench_) {
        for (const auto &h : qhadrons) {
            auto hp = h.vGetP();
            hjt_file_ << hp[0] << " " << hp[1] << " " << hp[2] << " " << h.GetMass() << " " << h.GetId() << " " << 0 << std::endl;
        }
        for (const auto &h : hhadrons) {
            auto hp = h.vGetP();
            hjt_file_ << hp[0] << " " << hp[1] << " " << hp[2] << " " << h.GetMass() << " " << h.GetId() << " " << 3 << std::endl;
        }

        if (do_wake_) {
            for (const auto &w : wake) {
                int ide_jt = (int(w.GetStatus()) == 1) ? 1 : 2;
                int wake_id;
                double wake_ch = w.GetCharge();
                if (w.GetMass() < 0.5) {
                    if (wake_ch == 0.) wake_id = 111;
                    else if (wake_ch == 1.) wake_id = 211;
                    else wake_id = -211;
                } else {
                    if (wake_ch == 1.) wake_id = 2212;
                    else wake_id = -2212;
                }
                auto wp = w.vGetP();
                hjt_file_ << wp[0] << " " << wp[1] << " " << wp[2] << " " << w.GetMass() << " " << wake_id << " " << ide_jt << std::endl;
            }
        }
    } else {
        for (const auto &h : vhadrons) {
            auto hp = h.vGetP();
            hjt_file_ << hp[0] << " " << hp[1] << " " << hp[2] << " " << h.GetMass() << " " << h.GetId() << " " << 0 << std::endl;
        }
    }

    hjt_file_ << "end" << std::endl;
}
