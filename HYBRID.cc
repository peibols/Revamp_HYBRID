#include "HYBRID.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <sstream>
#include "vector_operators.h"

namespace {
constexpr int kShowerSeedOffset = 33;
constexpr int kHybridSeedOffset = 1346;
constexpr int kLundSeedOffset = 2337;

int getSeedBase(const Config &cfg) {
    return cfg.getIntOr("seed_base", cfg.getIntOr("njob", 0));
}
}

HYBRID::HYBRID(const Config &cfg) :
      do_quench_(cfg.getBoolOr("do_quench", true)),
      do_wake_(cfg.getBoolOr("do_wake", true)),
      do_source_(cfg.getBoolOr("do_source", false)),
      do_elastic_(cfg.getBoolOr("do_elastic", false)),
      njob_(cfg.getIntOr("njob", 0)),
      Nev_(cfg.getIntOr("Nev", 1)),
      cent_(cfg.getStringOr("cent", "0-5")),
      kappa_(cfg.getDoubleOr("kappa", 1.0)),
      alpha_(cfg.getDoubleOr("alpha", 1.0)),
      tmethod_(cfg.getIntOr("tmethod", 0)),
      mode_(cfg.getIntOr("mode", 0)),
      ebe_hydro_(cfg.getIntOr("ebe_hydro", 0)),
      hadro_type_(cfg.getIntOr("hadro_type", cfg.getBoolOr("do_elastic", false) ? 1 : 0)),
      seed_base_(getSeedBase(cfg)),
      shower_seed_(seed_base_ + kShowerSeedOffset),
      hybrid_seed_(seed_base_ + kHybridSeedOffset),
      lund_seed_(seed_base_ + kLundSeedOffset),
      tables_path_(cfg.getStringOr("tables_path", "")),
      nr_(hybrid_seed_),
      tree_gen_(std::make_unique<TreeGenerator>()),
      hydro_profile_(std::make_unique<HydroProfile>()),
      wake_gen_(std::make_unique<WakeGenerator>()),
      lund_gen_(std::make_unique<LundGenerator>()),
      glauber_model_(std::make_unique<GlauberModel>()),
      energy_loss_(std::make_unique<EnergyLoss>(nr_, kappa_, alpha_, tmethod_, mode_,
                                                ebe_hydro_, do_elastic_, tables_path_,
                                                *hydro_profile_)) {

    // Open output files
    const std::string out_base = cfg.getStringOr("output_base", "HYBRID");
    hjt_file_.open(out_base + "_Hadrons.out", std::ios_base::app);
    pjt_file_.open(out_base + "_Partons.out", std::ios_base::app);

    std::cout << "Seed base= " << seed_base_
              << " shower= " << shower_seed_
              << " hybrid= " << hybrid_seed_
              << " lund= " << lund_seed_ << std::endl;
    if (do_elastic_) {
        std::cout << "Elastic scattering requested"
                  << " hadro_type= " << hadro_type_;
        if (!tables_path_.empty()) {
            std::cout << " tables_path= " << tables_path_;
        }
        std::cout << std::endl;
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
    hjt_file_.close();
    pjt_file_.close();
}

void HYBRID::run() {
    int count = 0;
    bool lund_initialized = false;

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
        if (count == 0) init_tree();
        double weight = 0.;
        double cross = 0.;
        double cross_err = 0.;
        while (true) {
          if (do_tree(partons, weight, cross, cross_err)) break;
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
            if (ebe_hydro_ == 0) {
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
            do_wake(quenched, partons, wake);
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
        if (!do_lund(partons, quenched, recoiled, vhadrons, qhadrons, hhadrons)) {
            std::cout << "Skipping event " << count << " after medium hadronization failure" << std::endl;
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
}

void HYBRID::init_tree() {
    tree_gen_->init(shower_seed_);
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

void HYBRID::do_wake(const std::vector<Quench> &quenched, const std::vector<Parton> &partons, std::vector<Wake> &wake) {
    wake_gen_->generate(quenched, partons, wake, nr_);
}

void HYBRID::init_lund() {
    lund_gen_->init(lund_seed_);
}

bool HYBRID::do_lund(const std::vector<Parton> &partons,
                     const std::vector<Quench> &quenched,
                     const std::vector<Quench> &recoiled,
                     std::vector<Hadron> &vhadrons,
                     std::vector<Hadron> &qhadrons,
                     std::vector<Hadron> &hhadrons) {
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
    if (had_is_ok && !holes.empty()) {
        hhadrons.clear();
        lund_gen_->hadronizeMedium(holes, hhadrons, hadro_type_);
    }
    return had_is_ok;
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
