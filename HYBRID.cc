#include "HYBRID.h"

#include <cmath>
#include <cstdlib>
#include <iostream>
#include <cassert>
#include <algorithm>
#include <sstream>
#include "vector_operators.h"

HYBRID::HYBRID(const Config &cfg) :
      do_quench_(cfg.getBoolOr("do_quench", true)),
      do_wake_(cfg.getBoolOr("do_wake", true)),
      do_source_(cfg.getBoolOr("do_source", false)),
      njob_(cfg.getIntOr("njob", 0)),
      Nev_(cfg.getIntOr("Nev", 1)),
      cent_(cfg.getStringOr("cent", "0-5")),
      kappa_(cfg.getDoubleOr("kappa", 1.0)),
      alpha_(cfg.getDoubleOr("alpha", 1.0)),
      tmethod_(cfg.getIntOr("tmethod", 0)),
      mode_(cfg.getIntOr("mode", 0)),
      ebe_hydro_(cfg.getIntOr("ebe_hydro", 0)),
      nr_(1346 + cfg.getIntOr("njob", 0)),
      tree_gen_(std::make_unique<TreeGenerator>()),
      hydro_profile_(std::make_unique<HydroProfile>()),
      wake_gen_(std::make_unique<WakeGenerator>()),
      lund_gen_(std::make_unique<LundGenerator>()),
      glauber_model_(std::make_unique<GlauberModel>()),
      energy_loss_(std::make_unique<EnergyLoss>(nr_, kappa_, alpha_, tmethod_, mode_, ebe_hydro_, *hydro_profile_)) {

    // Open output files
    const std::string out_base = cfg.getStringOr("output_base", "HYBRID");
    hjt_file_.open(out_base + "_Hadrons.out", std::ios_base::app);
    pjt_file_.open(out_base + "_Partons.out", std::ios_base::app);

    if (do_quench_) {
        if (ebe_hydro_ == 0) read_nuclear();
        else Ncollsize_ = read_nuclear_ipsat();
        read_hydro_ipsat();
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
    std::vector<Wake> wake;
    std::vector<Hadron> vhadrons;
    std::vector<Hadron> qhadrons;

    while (count < Nev_) {
        // Generate PYTHIA tree
        partons.clear();
        if (count == 0) init_tree();
        double weight = 0.;
        double cross = 0.;
        double cross_err = 0.;
        do_tree(partons, weight, cross, cross_err);

        // Create vector of quenched partons initially equal to vacuum partons
        quenched.clear();
        quenched.reserve(partons.size());
        for (const auto &parton : partons) {
            quenched.emplace_back(parton);
        }

        double x = 0., y = 0.;
        if (do_quench_) {
            // Generate x,y
            if (ebe_hydro_ == 0) {
                gxy(x, y);
            } else {
                int randNcoll = int(Ncollsize_ * double(rand()) / double(RAND_MAX));
                if (randNcoll == Ncollsize_) randNcoll = Ncollsize_ - 1;
                gxy_ipsat(x, y, randNcoll);
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

            do_eloss(partons, quenched, x, y);

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
        do_lund(partons, quenched, vhadrons, qhadrons);
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
            for (const auto &h : vhadrons) {
                auto hp = h.vGetP();
                for (int j = 0; j < 4; j++) sum_vac_hadrons[j] += hp[j];
            }
            for (const auto &h : qhadrons) {
                auto hp = h.vGetP();
                for (int j = 0; j < 4; j++) sum_med_hadrons[j] += hp[j];
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
        output_event(count, partons, quenched, vhadrons, qhadrons, wake, weight, cross, x, y);

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

void HYBRID::read_hydro_ipsat() {
    hydro_profile_->loadHydro(ebe_hydro_, cent_);
}

void HYBRID::init_tree() {
    tree_gen_->init(njob_);
}

void HYBRID::do_tree(std::vector<Parton> &partons, double &weight, double &cross, double &cross_err) {
    bool have_trigger;
    tree_gen_->nextEvent(partons, weight, cross, cross_err, have_trigger);
}

void HYBRID::gxy(double &x, double &y) {
    glauber_model_->sampleXY(x, y, nr_);
}

void HYBRID::gxy_ipsat(double &x, double &y, int randNcoll) {
    glauber_model_->sampleXYIPSAT(x, y, randNcoll, nr_);
}

void HYBRID::do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched, double x, double y) {
    energy_loss_->do_eloss(partons, quenched, x, y);
}

void HYBRID::do_wake(const std::vector<Quench> &quenched, const std::vector<Parton> &partons, std::vector<Wake> &wake) {
    wake_gen_->generate(quenched, partons, wake, nr_);
}

void HYBRID::init_lund() {
    lund_gen_->init();
}

void HYBRID::do_lund(const std::vector<Parton> &partons, const std::vector<Quench> &quenched, std::vector<Hadron> &vhadrons, std::vector<Hadron> &qhadrons) {
    lund_gen_->hadronizeVacuum(partons, vhadrons);
    lund_gen_->hadronizeMedium(quenched, qhadrons);
}

void HYBRID::output_event(int count,
                           const std::vector<Parton> &partons,
                           const std::vector<Quench> &quenched,
                           const std::vector<Hadron> &vhadrons,
                           const std::vector<Hadron> &qhadrons,
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
                pjt_file_ << qp[0] << " " << qp[1] << " " << qp[2] << " " << q.GetMass() << " " << q.GetId() << " " << 0 << std::endl;
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
