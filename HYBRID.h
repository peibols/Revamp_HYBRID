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
#include "HeavyQuarkEnergyLoss.h"

class HYBRID {
private:
    // Configuration flags
    bool do_quench_;
    bool do_wake_;
    bool do_source_;
    bool do_elastic_;
    bool do_lres_;
    bool do_moliere_on_unresolved_partons_;
    bool do_moliere_dynamic_unresolved_resolution_;
    bool do_moliere_dynamic_daughter_unresolved_resolution_;
    bool do_moliere_recursive_unresolved_resolution_;
    bool allow_modee_nonconserving_opening_;
    bool dump_hybrid_evolution_history_;
    bool do_event_display_;
    bool use_fixed_xy_;
    bool compat_moliere_legacy_hydro_;
    bool use_prehydro_;

    // Parameters
    int njob_;
    int Nev_;
    std::string cent_;
    double kappa_;
    double alpha_;
    int tmethod_;
    int mode_;
    heavy_quark::Parameters heavy_quark_parameters_;
    int ebe_hydro_;
    int hadro_type_;
    int max_tree_attempts_;
    int max_event_attempts_;
    double lres_rpower_;
    double moliere_unresolved_resolution_c_;
    double modee_max_opening_relative_residual_;
    int seed_base_;
    int shower_seed_;
    int hybrid_seed_;
    int lund_seed_;
    int elastic_seed_;
    double fixed_x_;
    double fixed_y_;
    std::string tables_path_;
    std::string prehydro_file_;
    std::string hybrid_evolution_history_file_;
    std::string event_display_file_;
    std::string pythia_cmnd_;
    long long generated_event_attempts_;
    long long tree_failures_;
    long long hadronization_failures_;

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

    void do_wake(const std::vector<Quench> &quenched, const std::vector<Parton> &partons,
                 const std::vector<Quench> &recoiled, std::vector<Wake> &wake);

    void init_lund();
    bool do_lund(const std::vector<Parton> &partons,
                 const std::vector<Quench> &quenched,
                 const std::vector<Quench> &recoiled,
                 std::vector<Hadron> &vhadrons,
                 std::vector<Hadron> &qhadrons,
                 std::vector<Hadron> &hhadrons,
                 std::string &failure_reason);

    void output_event(int count,
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
