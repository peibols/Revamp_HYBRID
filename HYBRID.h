#pragma once

#include <string>
#include <vector>
#include <fstream>
#include <memory>

#include "Config.h"
#include "Parton.h"
#include "Quench.h"
#include "Hadron.h"
#include "Wake.h"
#include "Random.h"
#include "TreeGenerator.h"
#include "HydroProfile.h"
#include "WakeGenerator.h"
#include "LundGenerator.h"
#include "GlauberModel.h"
#include "EnergyLoss.h"

class HYBRID {
private:
    // Configuration flags
    bool do_quench_;
    bool do_wake_;
    bool do_source_;
    bool do_elastic_;

    // Parameters
    int njob_;
    int Nev_;
    std::string cent_;
    double kappa_;
    double alpha_;
    int tmethod_;
    int mode_;
    int ebe_hydro_;
    int hadro_type_;
    int seed_base_;
    int shower_seed_;
    int hybrid_seed_;
    int lund_seed_;
    std::string tables_path_;

    // Random number generator
    numrand nr_;

    // Output files
    std::ofstream hjt_file_;
    std::ofstream pjt_file_;

    // Nuclear and hydro data
    int Ncollsize_;

    // Generator instances (replace global state)
    std::unique_ptr<TreeGenerator> tree_gen_;
    std::unique_ptr<HydroProfile> hydro_profile_;
    std::unique_ptr<WakeGenerator> wake_gen_;
    std::unique_ptr<LundGenerator> lund_gen_;
    std::unique_ptr<GlauberModel> glauber_model_;
    std::unique_ptr<EnergyLoss> energy_loss_;

    // Private methods for each step
    void read_nuclear();
    int read_nuclear_ipsat();
    void read_hydro();

    void init_tree();
    bool do_tree(std::vector<Parton> &partons, double &weight, double &cross, double &cross_err);

    void gxy(double &x, double &y);
    void gxy_ipsat(double &x, double &y);

    void do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched,
                  std::vector<Quench> &recoiled, double x, double y);

    void do_wake(const std::vector<Quench> &quenched, const std::vector<Parton> &partons, std::vector<Wake> &wake);

    void init_lund();
    bool do_lund(const std::vector<Parton> &partons,
                 const std::vector<Quench> &quenched,
                 const std::vector<Quench> &recoiled,
                 std::vector<Hadron> &vhadrons,
                 std::vector<Hadron> &qhadrons);

    void output_event(int count,
                      const std::vector<Parton> &partons,
                      const std::vector<Quench> &quenched,
                      const std::vector<Quench> &recoiled,
                      const std::vector<Hadron> &vhadrons,
                      const std::vector<Hadron> &qhadrons,
                      const std::vector<Wake> &wake,
                      double weight,
                      double cross,
                      double x,
                      double y);

public:
    // Constructor
    HYBRID(int njob, int Nev, const std::string &cent, double kappa, double alpha, int tmethod, bool do_quench, int mode, int ebe_hydro, bool do_wake = true, bool do_source = false);

    // Construct from a configuration file
    explicit HYBRID(const Config &cfg);

    // Destructor
    ~HYBRID();

    // Main run method
    void run();

    // Setters for additional options if needed
    void set_do_wake(bool do_wake) { do_wake_ = do_wake; }
    void set_do_source(bool do_source) { do_source_ = do_source; }
};
