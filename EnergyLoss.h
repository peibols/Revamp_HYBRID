#pragma once

#include <array>
#include <string>
#include <vector>
#include "Parton.h"
#include "Quench.h"
#include "Random.h"
#include "HydroProfile.h"

class EnergyLoss {
public:
    EnergyLoss(numrand &nr, double kappa, double alpha, int tmethod, int mode,
               int ebe_hydro, bool do_elastic, bool do_lres,
               bool do_moliere_on_unresolved_partons,
               bool do_moliere_dynamic_unresolved_resolution,
               bool do_moliere_dynamic_daughter_unresolved_resolution,
               bool do_moliere_recursive_unresolved_resolution,
               double moliere_unresolved_resolution_c,
               double lres_rpower,
               bool dump_hybrid_evolution_history,
               const std::string &hybrid_evolution_history_file,
               bool do_event_display,
               const std::string &event_display_file,
               bool compat_moliere_legacy_hydro,
               const std::string &tables_path,
               const HydroProfile &hydro_profile);
    ~EnergyLoss();

    // Perform energy loss on the given partons
    void do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched,
                  double x, double y, std::vector<Quench> *recoiled = nullptr);

private:
    numrand &nr_;
    double kappa_;
    double alpha_;
    int tmethod_;
    int mode_;
    int ebe_hydro_;
    bool do_elastic_;
    bool do_lres_;
    bool do_moliere_on_unresolved_partons_;
    bool do_moliere_dynamic_unresolved_resolution_;
    bool do_moliere_dynamic_daughter_unresolved_resolution_;
    bool do_moliere_recursive_unresolved_resolution_;
    double moliere_unresolved_resolution_c_;
    bool dump_hybrid_evolution_history_;
    std::string hybrid_evolution_history_file_;
    bool do_event_display_;
    std::string event_display_file_;
    int history_event_counter_;
    long long n_unresolved_segments_dynamic_;
    long long n_unresolved_candidate_scatters_;
    long long n_unresolved_coherent_scatters_;
    long long n_unresolved_resolving_scatters_;
    long long n_unresolved_pairs_elastically_decohered_;
    double sum_qperp_dperp_unresolved_candidates_;
    long long n_unresolved_segments_recursive_;
    long long n_recursive_frontier_candidates_;
    long long n_recursive_frontier_probe_batches_;
    long long n_recursive_frontier_probe_objects_;
    long long n_recursive_frontier_permutation_checks_;
    long long n_recursive_frontier_permutation_mismatches_;
    long long n_recursive_inner_resolutions_;
    long long n_recursive_outer_resolutions_;
    long long n_recursive_coherent_applications_;
    long long n_recursive_tree_updates_;
    long long n_recursive_opening_closure_checks_;
    double sum_recursive_opening_spatial_residual_;
    double max_recursive_opening_spatial_residual_;
    double sum_recursive_opening_energy_residual_;
    double max_recursive_opening_energy_residual_;
    long long n_recursive_live_dperp_tests_;
    long long n_recursive_vacuum_dperp_fallbacks_;
    bool compat_moliere_legacy_hydro_;
    double lres_rpower_;
    std::string tables_path_;
    const HydroProfile &hydro_profile_;

#ifdef HAVE_ROOT
    class TFile *event_display_root_file_;
    class TTree *event_display_tree_;
    class TTree *event_display_detail_tree_;
    int ed_event_id_;
    int ed_segment_id_;
    int ed_parton_index_;
    int ed_pdg_id_;
    int ed_parent_index_;
    int ed_d1_;
    int ed_d2_;
    int ed_is_colored_;
    int ed_is_unresolved_;
    int ed_had_scattering_;
    double ed_t_start_;
    double ed_x_start_;
    double ed_y_start_;
    double ed_z_start_;
    double ed_tau_start_;
    double ed_px_start_;
    double ed_py_start_;
    double ed_pz_start_;
    double ed_e_start_;
    double ed_t_end_;
    double ed_x_end_;
    double ed_y_end_;
    double ed_z_end_;
    double ed_tau_end_;
    double ed_px_end_;
    double ed_py_end_;
    double ed_pz_end_;
    double ed_e_end_;
    double ed_length_;
    double ed_tlength_;
    double ed_qperp_;
    std::string ed_segment_type_;
    int ed_record_id_;
    int ed_record_parton_index_;
    int ed_record_pdg_id_;
    int ed_record_parent_index_;
    int ed_record_d1_;
    int ed_record_d2_;
    int ed_record_related_index_;
    int ed_record_is_unresolved_;
    int ed_record_in_medium_;
    double ed_record_t_;
    double ed_record_x_;
    double ed_record_y_;
    double ed_record_z_;
    double ed_record_tau_;
    double ed_record_px_;
    double ed_record_py_;
    double ed_record_pz_;
    double ed_record_e_;
    double ed_record_t_end_;
    double ed_record_x_end_;
    double ed_record_y_end_;
    double ed_record_z_end_;
    double ed_record_tau_end_;
    double ed_record_px_end_;
    double ed_record_py_end_;
    double ed_record_pz_end_;
    double ed_record_e_end_;
    double ed_record_qperp_;
    double ed_record_temperature_;
    double ed_record_length_;
    double ed_record_tlength_;
    std::string ed_record_type_;
    std::string ed_record_label_;
#endif

    // Private member functions for energy loss calculations
    void do_eloss_impl(const std::vector<Parton> &partons, std::vector<Quench> &quenched, double x, double y);
    void do_lres_eloss_impl(const std::vector<Parton> &partons, std::vector<Quench> &quenched,
                            double x, double y, std::vector<Quench> *recoiled);
    void loss_rate(std::array<double,4> &p, std::array<double,4> &pos, double tof, int id,
                   double &length, double &tlength,
                   int event_id = -1, int *record_id = nullptr,
                   int parton_index = -1, int parent_index = -1,
                   int d1 = -1, int d2 = -1, bool is_unresolved = false);
    double resolution_time(double parent_e, double parent_px, double parent_py, double parent_pz,
                           double x, double y, double z,
                           double dvx, double dvy, double dvz,
                           double t0) const;
    void append_history_records(const std::vector<std::string> &records);
    void init_event_display();
    void close_event_display();
    void fill_event_display_segment(int event_id, int segment_id, int parton_index,
                                    const Parton &parton, int parent_index, int d1, int d2,
                                    bool is_unresolved, int had_scattering,
                                    const std::array<double,4> &pos_start,
                                    const std::array<double,4> &p_start,
                                    const std::array<double,4> &pos_end,
                                    const std::array<double,4> &p_end,
                                    double length, double tlength,
                                    const std::string &segment_type);
    void fill_event_display_record(int event_id, int record_id, int parton_index,
                                   int pdg_id, int parent_index, int d1, int d2,
                                   int related_index, bool is_unresolved, bool in_medium,
                                   const std::array<double,4> &pos_start,
                                   const std::array<double,4> &p_start,
                                   const std::array<double,4> &pos_end,
                                   const std::array<double,4> &p_end,
                                   double qperp, double temperature,
                                   double length, double tlength,
                                   const std::string &record_type,
                                   const std::string &label);
    void get_source_evol(double &tau_ev, double& x_f, double& y_f, double& vx_f, double& vy_f, double tau_ini, double x_ini, double y_ini, double Tc);
    double call_gT(double tau, double x, double y, int comp) const;
    void quenched_sons(const std::array<double,4> &p, const std::array<double,4> &qp, std::array<double,4> &d1, std::array<double,4> &d2);
    double normalise(std::array<double,4> &p);
    std::array<double,4> vec_prod(const std::array<double,4> &a, const std::array<double,4> &b);
    void trans_kick(const std::array<double,4> &w, double w2, const std::array<double,4> &v, std::array<double,4> &p, double temp, double vscalw, double lore, double step, double kappa);
};
