#include "EnergyLoss.h"
#include <array>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cstdint>
#include <cmath>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <utility>
#include "MoliereTables.h"
#include "MoliereElastic.h"
#include "vector_operators.h"
#ifdef HAVE_ROOT
#include "TFile.h"
#include "TTree.h"
#endif

namespace {
constexpr double kLresFinalFlightTime = 10000.;
constexpr double kLresInfiniteTime = 100000000000.;

std::uint64_t splitmix64(std::uint64_t x) {
    x += 0x9e3779b97f4a7c15ULL;
    x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
    x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
    return x ^ (x >> 31);
}

std::uint64_t hashCombine(std::uint64_t seed, std::uint64_t value) {
    return splitmix64(seed ^ (value + 0x9e3779b97f4a7c15ULL + (seed << 6) + (seed >> 2)));
}

int modeEBranchLocalSeed(int base_seed, int event_id, int segment_id,
                         int iteration, int active_ancestor, int probe_id) {
    std::uint64_t seed = 0x6d6f6465455f524eULL;  // "modeE_RN" tag
    seed = hashCombine(seed, static_cast<std::uint64_t>(base_seed));
    seed = hashCombine(seed, static_cast<std::uint64_t>(event_id));
    seed = hashCombine(seed, static_cast<std::uint64_t>(segment_id));
    seed = hashCombine(seed, static_cast<std::uint64_t>(iteration));
    seed = hashCombine(seed, static_cast<std::uint64_t>(active_ancestor + 1000003));
    seed = hashCombine(seed, static_cast<std::uint64_t>(probe_id + 2000003));

    // numrand stores an int seed. Keep the seed positive and nonzero while
    // retaining deterministic branch-local independence for Mode E probes.
    return static_cast<int>(seed % 2147483646ULL) + 1;
}

int modeECoherentSourceSeed(int base_seed, int event_id, int segment_id,
                            int iteration, int active_ancestor) {
    std::uint64_t seed = 0x6d6f6465455f434fULL;  // "modeE_CO" tag
    seed = hashCombine(seed, static_cast<std::uint64_t>(base_seed));
    seed = hashCombine(seed, static_cast<std::uint64_t>(event_id));
    seed = hashCombine(seed, static_cast<std::uint64_t>(segment_id));
    seed = hashCombine(seed, static_cast<std::uint64_t>(iteration));
    seed = hashCombine(seed, static_cast<std::uint64_t>(active_ancestor + 1000003));
    return static_cast<int>(seed % 2147483646ULL) + 1;
}

bool isColored(int id) {
    return std::abs(id) <= 6 || id == 21;
}

double safeFormationTime(const Parton &p) {
    const double q = p.GetQ();
    if (q == 0.) return 0.;
    return 0.2 * 2. * p.vGetP()[3] / (q * q);
}

std::array<double,4> velocity(const std::array<double,4> &p) {
    if (p[3] == 0.) return {0., 0., 0., 1.};
    return {p[0] / p[3], p[1] / p[3], p[2] / p[3], 1.};
}

std::array<double,4> causalVelocity(const std::array<double,4> &p) {
    auto v = velocity(p);
    const double vperp = std::sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
    if (vperp > 1. - 1.e-9 && vperp > 0.) {
        const double scale = (1. - 1.e-9) / vperp;
        v[0] *= scale;
        v[1] *= scale;
        v[2] *= scale;
    }
    v[3] = 1.;
    return v;
}

std::array<double,4> orientationFor(const std::array<double,4> &p) {
    if (p[3] == 0.) return {0., 0., 0., 1.};
    return {p[0] / p[3], p[1] / p[3], p[2] / p[3], 1.};
}

std::array<double,3> spatialUnit(const std::array<double,4> &p) {
    const double norm = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    if (norm <= 0.) return {0., 0., 1.};
    return {p[0] / norm, p[1] / norm, p[2] / norm};
}

double dot3(const std::array<double,3> &a, const std::array<double,3> &b) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

std::array<double,3> cross3(const std::array<double,3> &a,
                            const std::array<double,3> &b) {
    return {a[1] * b[2] - a[2] * b[1],
            a[2] * b[0] - a[0] * b[2],
            a[0] * b[1] - a[1] * b[0]};
}

double norm3(const std::array<double,3> &a) {
    return std::sqrt(dot3(a, a));
}

double spatialMomentumResidual(const std::array<double,4> &parent_p,
                               const std::array<double,4> &child1_p,
                               const std::array<double,4> &child2_p) {
    const double dx = parent_p[0] - child1_p[0] - child2_p[0];
    const double dy = parent_p[1] - child1_p[1] - child2_p[1];
    const double dz = parent_p[2] - child1_p[2] - child2_p[2];
    return std::sqrt(dx * dx + dy * dy + dz * dz);
}

double energyResidual(const std::array<double,4> &parent_p,
                      const std::array<double,4> &child1_p,
                      const std::array<double,4> &child2_p) {
    return std::abs(parent_p[3] - child1_p[3] - child2_p[3]);
}

std::array<double,3> normalized3(const std::array<double,3> &a,
                                 const std::array<double,3> &fallback) {
    const double n = norm3(a);
    if (n <= 0.) return fallback;
    return {a[0] / n, a[1] / n, a[2] / n};
}

std::array<double,4> rotateSpatialFromParentAxis(const std::array<double,4> &child_p,
                                                 const std::array<double,4> &vac_parent_p,
                                                 const std::array<double,4> &live_parent_p) {
    auto out = child_p;
    const auto from = spatialUnit(vac_parent_p);
    const auto to = spatialUnit(live_parent_p);
    const auto axis_cross = cross3(from, to);
    const double sin_angle = norm3(axis_cross);
    const double cos_angle = std::max(-1., std::min(1., dot3(from, to)));
    const std::array<double,3> r = {child_p[0], child_p[1], child_p[2]};

    std::array<double,3> rotated = r;
    if (sin_angle > 1.e-12) {
        const std::array<double,3> k = {axis_cross[0] / sin_angle,
                                        axis_cross[1] / sin_angle,
                                        axis_cross[2] / sin_angle};
        const auto k_cross_r = cross3(k, r);
        const double k_dot_r = dot3(k, r);
        rotated = {r[0] * cos_angle + k_cross_r[0] * sin_angle + k[0] * k_dot_r * (1. - cos_angle),
                   r[1] * cos_angle + k_cross_r[1] * sin_angle + k[1] * k_dot_r * (1. - cos_angle),
                   r[2] * cos_angle + k_cross_r[2] * sin_angle + k[2] * k_dot_r * (1. - cos_angle)};
    } else if (cos_angle < 0.) {
        std::array<double,3> trial_axis = std::abs(from[0]) < 0.9
                                             ? std::array<double,3>{1., 0., 0.}
                                             : std::array<double,3>{0., 1., 0.};
        const auto axis = normalized3(cross3(from, trial_axis), {0., 0., 1.});
        const double axis_dot_r = dot3(axis, r);
        rotated = {-r[0] + 2. * axis[0] * axis_dot_r,
                   -r[1] + 2. * axis[1] * axis_dot_r,
                   -r[2] + 2. * axis[2] * axis_dot_r};
    }

    out[0] = rotated[0];
    out[1] = rotated[1];
    out[2] = rotated[2];
    return out;
}

std::array<double,4> scaledMasslessMomentum(const std::array<double,4> &axis,
                                               double energy,
                                               const std::array<double,4> &fallback_axis) {
    auto direction = spatialUnit(axis);
    if (axis[0] == 0. && axis[1] == 0. && axis[2] == 0.) direction = spatialUnit(fallback_axis);
    return {direction[0] * energy, direction[1] * energy, direction[2] * energy, energy};
}

void mapChildMomentaFromLiveParent(const std::array<double,4> &vac_parent_p,
                                   const std::array<double,4> &live_parent_p,
                                   const std::array<double,4> &vac_child1_p,
                                   const std::array<double,4> &vac_child2_p,
                                   std::array<double,4> &p1,
                                   std::array<double,4> &p2) {
    const double vac_child_sum_e = vac_child1_p[3] + vac_child2_p[3];
    const double share1 = vac_child_sum_e > 0. ? vac_child1_p[3] / vac_child_sum_e : 0.5;
    const double share2 = 1. - share1;
    const double e1 = std::max(0., live_parent_p[3] * share1);
    const double e2 = std::max(0., live_parent_p[3] * share2);

    // A massless coherent parent cannot be split into two separated on-shell
    // daughters while preserving its full four-vector exactly.  For Mode E we
    // preserve the live parent deflection and total energy, then validate the
    // residual spatial mismatch through event-level studies.
    const auto rotated1 = rotateSpatialFromParentAxis(vac_child1_p, vac_parent_p, live_parent_p);
    const auto rotated2 = rotateSpatialFromParentAxis(vac_child2_p, vac_parent_p, live_parent_p);
    p1 = scaledMasslessMomentum(rotated1, e1, live_parent_p);
    p2 = scaledMasslessMomentum(rotated2, e2, live_parent_p);
}

std::array<double,4> separatedChildPosition(const std::array<double,4> &parent_pos,
                                            const std::array<double,4> &parent_p,
                                            const std::array<double,4> &child_p,
                                            double split_time,
                                            double target_time) {
    auto pos = parent_pos;
    const double dt = std::max(0., target_time - split_time);
    const auto vp = causalVelocity(parent_p);
    const auto vc = causalVelocity(child_p);
    pos[0] += (vc[0] - vp[0]) * dt;
    pos[1] += (vc[1] - vp[1]) * dt;
    pos[2] += (vc[2] - vp[2]) * dt;
    pos[3] = target_time;
    return pos;
}

double transverseSeparation(const std::array<double,4> &pos1,
                            const std::array<double,4> &pos2) {
    const double dx = pos1[0] - pos2[0];
    const double dy = pos1[1] - pos2[1];
    return std::sqrt(dx * dx + dy * dy);
}

double phiFromP(const std::array<double,4> &p) {
    return std::atan2(p[1], p[0]);
}

double deltaPhiXYFromP(const std::array<double,4> &a, const std::array<double,4> &b) {
    double dphi = phiFromP(a) - phiFromP(b);
    while (dphi > M_PI) dphi -= 2. * M_PI;
    while (dphi < -M_PI) dphi += 2. * M_PI;
    return std::abs(dphi);
}

std::vector<double> toVector4(const std::array<double,4> &p) {
    return {p[0], p[1], p[2], p[3]};
}

double properTimeFromPos(const std::array<double,4> &pos) {
    const double tau2 = pos[3] * pos[3] - pos[2] * pos[2];
    return tau2 > 0. ? std::sqrt(tau2) : 0.;
}

double compute_unresolved_pair_dperp(const std::array<double,4> &p1_vac,
                                     const std::array<double,4> &p2_vac,
                                     double split_time,
                                     double scattering_time) {
    const double dt = std::max(0., scattering_time - split_time);
    const auto v1 = velocity(p1_vac);
    const auto v2 = velocity(p2_vac);
    const double dx = (v1[0] - v2[0]) * dt;
    const double dy = (v1[1] - v2[1]) * dt;
    return std::sqrt(dx * dx + dy * dy);
}

bool passes_dynamic_moliere_resolution_test(double qperp, double dperp, double c_res) {
    return qperp * dperp > c_res;
}

moliere::ScatteringDecision apply_coherent_unresolved_kick() {
    return moliere::ScatteringDecision::Apply;
}

void apply_resolved_daughter_kick(const moliere::ScatteringCandidate &candidate,
                                  std::array<double,4> &p,
                                  std::vector<Quench> &new_particles,
                                  int &had_scattering,
                                  std::array<double,4> &orient) {
    for (int i = 0; i < 4; ++i) {
        p[i] += candidate.p_after[i] - candidate.p_before[i];
    }
    if (p[3] <= 0.) {
        p = {0., 0., 0., 0.};
    }
    had_scattering = 1;
    orient = orientationFor(p);

    new_particles.emplace_back(Parton(toVector4(candidate.recoiler_p), 100000000., 0., 0,
                                      -1, -1, candidate.recoiler_id, "recoiler", 0, 0, false));
    new_particles.back().vSetRi(candidate.pos);
    new_particles.emplace_back(Parton(toVector4(candidate.hole_p), 100000000., 0., 0,
                                      -1, -1, candidate.hole_id, "hole", 0, 0, true));
    new_particles.back().vSetRi(candidate.pos);
}

struct LresLifetime {
    double resolve = 0.;
    double creation = 0.;
    double finish = 0.;
    double resolve_abs = 0.;
    double effective_time = 0.;
    std::array<double,4> ri = {0., 0., 0., 0.};
    std::array<double,4> rf = {0., 0., 0., 0.};
};

struct LresState {
    std::array<double,4> p = {0., 0., 0., 0.};
    std::array<double,4> r = {0., 0., 0., 0.};
};
}

EnergyLoss::EnergyLoss(numrand &nr, double kappa, double alpha, int tmethod, int mode,
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
                       const HydroProfile &hydro_profile)
    : nr_(nr), kappa_(kappa), alpha_(alpha), tmethod_(tmethod), mode_(mode),
      ebe_hydro_(ebe_hydro), do_elastic_(do_elastic), do_lres_(do_lres),
      do_moliere_on_unresolved_partons_(do_moliere_on_unresolved_partons),
      do_moliere_dynamic_unresolved_resolution_(do_moliere_dynamic_unresolved_resolution),
      do_moliere_dynamic_daughter_unresolved_resolution_(do_moliere_dynamic_daughter_unresolved_resolution),
      do_moliere_recursive_unresolved_resolution_(do_moliere_recursive_unresolved_resolution),
      moliere_unresolved_resolution_c_(moliere_unresolved_resolution_c),
      dump_hybrid_evolution_history_(dump_hybrid_evolution_history),
      hybrid_evolution_history_file_(hybrid_evolution_history_file),
      do_event_display_(do_event_display),
      event_display_file_(event_display_file),
      history_event_counter_(0),
      n_unresolved_segments_dynamic_(0),
      n_unresolved_candidate_scatters_(0),
      n_unresolved_coherent_scatters_(0),
      n_unresolved_resolving_scatters_(0),
      n_unresolved_pairs_elastically_decohered_(0),
      sum_qperp_dperp_unresolved_candidates_(0.),
      n_unresolved_segments_recursive_(0),
      n_recursive_frontier_candidates_(0),
      n_recursive_frontier_probe_batches_(0),
      n_recursive_frontier_probe_objects_(0),
      n_recursive_frontier_permutation_checks_(0),
      n_recursive_frontier_permutation_mismatches_(0),
      n_recursive_inner_resolutions_(0),
      n_recursive_outer_resolutions_(0),
      n_recursive_coherent_applications_(0),
      n_recursive_failed_daughter_vetoes_(0),
      n_recursive_coherent_resample_requests_(0),
      n_recursive_coherent_resample_candidates_(0),
      n_recursive_coherent_candidate_vetoes_(0),
      n_recursive_coherent_candidate_accepts_(0),
      n_recursive_coherent_resample_exhausted_(0),
      n_recursive_color_neutral_parent_skips_(0),
      n_recursive_tree_updates_(0),
      n_recursive_opening_closure_checks_(0),
      sum_recursive_opening_spatial_residual_(0.),
      max_recursive_opening_spatial_residual_(0.),
      sum_recursive_opening_energy_residual_(0.),
      max_recursive_opening_energy_residual_(0.),
      n_recursive_live_dperp_tests_(0),
      n_recursive_vacuum_dperp_fallbacks_(0),
      compat_moliere_legacy_hydro_(compat_moliere_legacy_hydro),
      lres_rpower_(lres_rpower), tables_path_(tables_path),
      hydro_profile_(hydro_profile)
#ifdef HAVE_ROOT
      , event_display_root_file_(nullptr), event_display_tree_(nullptr),
      event_display_detail_tree_(nullptr),
      ed_event_id_(0), ed_segment_id_(0), ed_parton_index_(0), ed_pdg_id_(0),
      ed_parent_index_(0), ed_d1_(0), ed_d2_(0), ed_is_colored_(0),
      ed_is_unresolved_(0), ed_had_scattering_(0), ed_t_start_(0.),
      ed_x_start_(0.), ed_y_start_(0.), ed_z_start_(0.), ed_tau_start_(0.),
      ed_px_start_(0.), ed_py_start_(0.), ed_pz_start_(0.), ed_e_start_(0.),
      ed_t_end_(0.), ed_x_end_(0.), ed_y_end_(0.), ed_z_end_(0.),
      ed_tau_end_(0.), ed_px_end_(0.), ed_py_end_(0.), ed_pz_end_(0.),
      ed_e_end_(0.), ed_length_(0.), ed_tlength_(0.), ed_qperp_(0.),
      ed_segment_type_(""), ed_record_id_(0), ed_record_parton_index_(0),
      ed_record_pdg_id_(0), ed_record_parent_index_(0), ed_record_d1_(0),
      ed_record_d2_(0), ed_record_related_index_(0), ed_record_is_unresolved_(0),
      ed_record_in_medium_(0), ed_record_t_(0.), ed_record_x_(0.),
      ed_record_y_(0.), ed_record_z_(0.), ed_record_tau_(0.),
      ed_record_px_(0.), ed_record_py_(0.), ed_record_pz_(0.), ed_record_e_(0.),
      ed_record_t_end_(0.), ed_record_x_end_(0.), ed_record_y_end_(0.),
      ed_record_z_end_(0.), ed_record_tau_end_(0.), ed_record_px_end_(0.),
      ed_record_py_end_(0.), ed_record_pz_end_(0.), ed_record_e_end_(0.),
      ed_record_qperp_(0.), ed_record_temperature_(0.), ed_record_length_(0.),
      ed_record_tlength_(0.), ed_record_type_(""), ed_record_label_("")
#endif
{
    if (do_elastic_) {
        MoliereTables::ensureLoaded(tables_path_);
    }
    init_event_display();
}

EnergyLoss::~EnergyLoss() {
    close_event_display();
    if (n_unresolved_segments_dynamic_ > 0) {
        const double avg_qd =
            n_unresolved_candidate_scatters_ > 0
                ? sum_qperp_dperp_unresolved_candidates_ /
                      static_cast<double>(n_unresolved_candidate_scatters_)
                : 0.;
        std::cout << "Dynamic unresolved Moliere diagnostics:"
                  << " n_unresolved_segments_dynamic= " << n_unresolved_segments_dynamic_
                  << " n_unresolved_candidate_scatters= " << n_unresolved_candidate_scatters_
                  << " n_unresolved_coherent_scatters= " << n_unresolved_coherent_scatters_
                  << " n_unresolved_resolving_scatters= " << n_unresolved_resolving_scatters_
                  << " n_unresolved_pairs_elastically_decohered= "
                  << n_unresolved_pairs_elastically_decohered_
                  << " average_qperp_dperp_for_candidate_scatters= " << avg_qd
                  << std::endl;
    }
    if (n_unresolved_segments_recursive_ > 0) {
        const long long classified_recursive_candidates =
            n_recursive_coherent_applications_ +
            n_unresolved_resolving_scatters_ +
            n_recursive_failed_daughter_vetoes_ +
            n_recursive_coherent_candidate_vetoes_;
        const long long recursive_candidate_accounting_delta =
            n_unresolved_candidate_scatters_ - classified_recursive_candidates;
        const long long coherent_candidate_accounting_delta =
            n_recursive_coherent_resample_candidates_ -
            n_recursive_coherent_candidate_accepts_ -
            n_recursive_coherent_candidate_vetoes_;
        std::cout << "Recursive unresolved Moliere diagnostics:"
                  << " n_unresolved_segments_recursive= " << n_unresolved_segments_recursive_
                  << " n_recursive_frontier_candidates= " << n_recursive_frontier_candidates_
                  << " n_recursive_frontier_probe_batches= "
                  << n_recursive_frontier_probe_batches_
                  << " n_recursive_frontier_probe_objects= "
                  << n_recursive_frontier_probe_objects_
                  << " n_recursive_frontier_permutation_checks= "
                  << n_recursive_frontier_permutation_checks_
                  << " n_recursive_frontier_permutation_mismatches= "
                  << n_recursive_frontier_permutation_mismatches_
                  << " n_recursive_inner_resolutions= " << n_recursive_inner_resolutions_
                  << " n_recursive_outer_resolutions= " << n_recursive_outer_resolutions_
                  << " n_recursive_coherent_applications= " << n_recursive_coherent_applications_
                  << " n_recursive_failed_daughter_vetoes= "
                  << n_recursive_failed_daughter_vetoes_
                  << " n_recursive_coherent_resample_requests= "
                  << n_recursive_coherent_resample_requests_
                  << " n_recursive_coherent_resample_candidates= "
                  << n_recursive_coherent_resample_candidates_
                  << " n_recursive_coherent_candidate_vetoes= "
                  << n_recursive_coherent_candidate_vetoes_
                  << " n_recursive_coherent_candidate_accepts= "
                  << n_recursive_coherent_candidate_accepts_
                  << " n_recursive_coherent_resample_exhausted= "
                  << n_recursive_coherent_resample_exhausted_
                  << " n_recursive_color_neutral_parent_skips= "
                  << n_recursive_color_neutral_parent_skips_
                  << " recursive_candidate_accounting_delta= "
                  << recursive_candidate_accounting_delta
                  << " coherent_candidate_accounting_delta= "
                  << coherent_candidate_accounting_delta
                  << " n_recursive_tree_updates= " << n_recursive_tree_updates_
                  << " n_recursive_opening_closure_checks= "
                  << n_recursive_opening_closure_checks_
                  << " avg_recursive_opening_spatial_residual= "
                  << (n_recursive_opening_closure_checks_ > 0
                          ? sum_recursive_opening_spatial_residual_ /
                                static_cast<double>(n_recursive_opening_closure_checks_)
                          : 0.)
                  << " max_recursive_opening_spatial_residual= "
                  << max_recursive_opening_spatial_residual_
                  << " avg_recursive_opening_energy_residual= "
                  << (n_recursive_opening_closure_checks_ > 0
                          ? sum_recursive_opening_energy_residual_ /
                                static_cast<double>(n_recursive_opening_closure_checks_)
                          : 0.)
                  << " max_recursive_opening_energy_residual= "
                  << max_recursive_opening_energy_residual_
                  << " n_recursive_live_dperp_tests= " << n_recursive_live_dperp_tests_
                  << " n_recursive_vacuum_dperp_fallbacks= "
                  << n_recursive_vacuum_dperp_fallbacks_
                  << std::endl;
    }
}

void EnergyLoss::init_event_display() {
    if (!do_event_display_) return;
    if (event_display_file_.empty()) event_display_file_ = "eventDisplay.root";
#ifdef HAVE_ROOT
    event_display_root_file_ = TFile::Open(event_display_file_.c_str(), "RECREATE");
    if (event_display_root_file_ == nullptr || event_display_root_file_->IsZombie()) {
        std::cerr << "Event-display ROOT output disabled: could not open "
                  << event_display_file_ << std::endl;
        do_event_display_ = false;
        return;
    }
    event_display_tree_ = new TTree("PartonSegments",
                                    "Parton location, energy, and momentum by propagation segment");
    event_display_tree_->Branch("event_id", &ed_event_id_);
    event_display_tree_->Branch("segment_id", &ed_segment_id_);
    event_display_tree_->Branch("parton_index", &ed_parton_index_);
    event_display_tree_->Branch("pdg_id", &ed_pdg_id_);
    event_display_tree_->Branch("parent_index", &ed_parent_index_);
    event_display_tree_->Branch("d1", &ed_d1_);
    event_display_tree_->Branch("d2", &ed_d2_);
    event_display_tree_->Branch("is_colored", &ed_is_colored_);
    event_display_tree_->Branch("is_unresolved", &ed_is_unresolved_);
    event_display_tree_->Branch("had_scattering", &ed_had_scattering_);
    event_display_tree_->Branch("t_start", &ed_t_start_);
    event_display_tree_->Branch("x_start", &ed_x_start_);
    event_display_tree_->Branch("y_start", &ed_y_start_);
    event_display_tree_->Branch("z_start", &ed_z_start_);
    event_display_tree_->Branch("tau_start", &ed_tau_start_);
    event_display_tree_->Branch("px_start", &ed_px_start_);
    event_display_tree_->Branch("py_start", &ed_py_start_);
    event_display_tree_->Branch("pz_start", &ed_pz_start_);
    event_display_tree_->Branch("e_start", &ed_e_start_);
    event_display_tree_->Branch("t_end", &ed_t_end_);
    event_display_tree_->Branch("x_end", &ed_x_end_);
    event_display_tree_->Branch("y_end", &ed_y_end_);
    event_display_tree_->Branch("z_end", &ed_z_end_);
    event_display_tree_->Branch("tau_end", &ed_tau_end_);
    event_display_tree_->Branch("px_end", &ed_px_end_);
    event_display_tree_->Branch("py_end", &ed_py_end_);
    event_display_tree_->Branch("pz_end", &ed_pz_end_);
    event_display_tree_->Branch("e_end", &ed_e_end_);
    event_display_tree_->Branch("length", &ed_length_);
    event_display_tree_->Branch("tlength", &ed_tlength_);
    event_display_tree_->Branch("qperp", &ed_qperp_);
    event_display_tree_->Branch("segment_type", &ed_segment_type_);
    event_display_detail_tree_ = new TTree("DetailedRecords",
                                           "Time-ordered shower, propagation, Moliere, and LRES diagnostic records");
    event_display_detail_tree_->Branch("event_id", &ed_event_id_);
    event_display_detail_tree_->Branch("record_id", &ed_record_id_);
    event_display_detail_tree_->Branch("parton_index", &ed_record_parton_index_);
    event_display_detail_tree_->Branch("pdg_id", &ed_record_pdg_id_);
    event_display_detail_tree_->Branch("parent_index", &ed_record_parent_index_);
    event_display_detail_tree_->Branch("d1", &ed_record_d1_);
    event_display_detail_tree_->Branch("d2", &ed_record_d2_);
    event_display_detail_tree_->Branch("related_index", &ed_record_related_index_);
    event_display_detail_tree_->Branch("is_unresolved", &ed_record_is_unresolved_);
    event_display_detail_tree_->Branch("in_medium", &ed_record_in_medium_);
    event_display_detail_tree_->Branch("t", &ed_record_t_);
    event_display_detail_tree_->Branch("x", &ed_record_x_);
    event_display_detail_tree_->Branch("y", &ed_record_y_);
    event_display_detail_tree_->Branch("z", &ed_record_z_);
    event_display_detail_tree_->Branch("tau", &ed_record_tau_);
    event_display_detail_tree_->Branch("px", &ed_record_px_);
    event_display_detail_tree_->Branch("py", &ed_record_py_);
    event_display_detail_tree_->Branch("pz", &ed_record_pz_);
    event_display_detail_tree_->Branch("e", &ed_record_e_);
    event_display_detail_tree_->Branch("t_end", &ed_record_t_end_);
    event_display_detail_tree_->Branch("x_end", &ed_record_x_end_);
    event_display_detail_tree_->Branch("y_end", &ed_record_y_end_);
    event_display_detail_tree_->Branch("z_end", &ed_record_z_end_);
    event_display_detail_tree_->Branch("tau_end", &ed_record_tau_end_);
    event_display_detail_tree_->Branch("px_end", &ed_record_px_end_);
    event_display_detail_tree_->Branch("py_end", &ed_record_py_end_);
    event_display_detail_tree_->Branch("pz_end", &ed_record_pz_end_);
    event_display_detail_tree_->Branch("e_end", &ed_record_e_end_);
    event_display_detail_tree_->Branch("qperp", &ed_record_qperp_);
    event_display_detail_tree_->Branch("temperature", &ed_record_temperature_);
    event_display_detail_tree_->Branch("length", &ed_record_length_);
    event_display_detail_tree_->Branch("tlength", &ed_record_tlength_);
    event_display_detail_tree_->Branch("record_type", &ed_record_type_);
    event_display_detail_tree_->Branch("label", &ed_record_label_);
#else
    std::cerr << "Event-display ROOT output requested but this binary was built without ROOT. "
              << "Rebuild with root-config available." << std::endl;
    do_event_display_ = false;
#endif
}

void EnergyLoss::close_event_display() {
#ifdef HAVE_ROOT
    if (event_display_root_file_ != nullptr) {
        event_display_root_file_->cd();
        if (event_display_tree_ != nullptr) event_display_tree_->Write();
        if (event_display_detail_tree_ != nullptr) event_display_detail_tree_->Write();
        event_display_root_file_->Close();
        delete event_display_root_file_;
        event_display_root_file_ = nullptr;
        event_display_tree_ = nullptr;
        event_display_detail_tree_ = nullptr;
    }
#endif
}

void EnergyLoss::fill_event_display_segment(int event_id, int segment_id, int parton_index,
                                            const Parton &parton, int parent_index, int d1, int d2,
                                            bool is_unresolved, int had_scattering,
                                            const std::array<double,4> &pos_start,
                                            const std::array<double,4> &p_start,
                                            const std::array<double,4> &pos_end,
                                            const std::array<double,4> &p_end,
                                            double length, double tlength,
                                            const std::string &segment_type) {
    if (!do_event_display_) return;
#ifdef HAVE_ROOT
    if (event_display_tree_ == nullptr) return;
    ed_event_id_ = event_id;
    ed_segment_id_ = segment_id;
    ed_parton_index_ = parton_index;
    ed_pdg_id_ = parton.GetId();
    ed_parent_index_ = parent_index;
    ed_d1_ = d1;
    ed_d2_ = d2;
    ed_is_colored_ = isColored(parton.GetId()) ? 1 : 0;
    ed_is_unresolved_ = is_unresolved ? 1 : 0;
    ed_had_scattering_ = had_scattering;
    ed_t_start_ = pos_start[3];
    ed_x_start_ = pos_start[0];
    ed_y_start_ = pos_start[1];
    ed_z_start_ = pos_start[2];
    ed_tau_start_ = properTimeFromPos(pos_start);
    ed_px_start_ = p_start[0];
    ed_py_start_ = p_start[1];
    ed_pz_start_ = p_start[2];
    ed_e_start_ = p_start[3];
    ed_t_end_ = pos_end[3];
    ed_x_end_ = pos_end[0];
    ed_y_end_ = pos_end[1];
    ed_z_end_ = pos_end[2];
    ed_tau_end_ = properTimeFromPos(pos_end);
    ed_px_end_ = p_end[0];
    ed_py_end_ = p_end[1];
    ed_pz_end_ = p_end[2];
    ed_e_end_ = p_end[3];
    ed_length_ = length;
    ed_tlength_ = tlength;
    ed_qperp_ = std::sqrt((p_end[0] - p_start[0]) * (p_end[0] - p_start[0]) +
                          (p_end[1] - p_start[1]) * (p_end[1] - p_start[1]));
    ed_segment_type_ = segment_type;
    event_display_tree_->Fill();
#else
    (void)event_id;
    (void)segment_id;
    (void)parton_index;
    (void)parton;
    (void)parent_index;
    (void)d1;
    (void)d2;
    (void)is_unresolved;
    (void)had_scattering;
    (void)pos_start;
    (void)p_start;
    (void)pos_end;
    (void)p_end;
    (void)length;
    (void)tlength;
    (void)segment_type;
#endif
}

void EnergyLoss::fill_event_display_record(int event_id, int record_id, int parton_index,
                                           int pdg_id, int parent_index, int d1, int d2,
                                           int related_index, bool is_unresolved, bool in_medium,
                                           const std::array<double,4> &pos_start,
                                           const std::array<double,4> &p_start,
                                           const std::array<double,4> &pos_end,
                                           const std::array<double,4> &p_end,
                                           double qperp, double temperature,
                                           double length, double tlength,
                                           const std::string &record_type,
                                           const std::string &label) {
    if (!do_event_display_) return;
#ifdef HAVE_ROOT
    if (event_display_detail_tree_ == nullptr) return;
    ed_event_id_ = event_id;
    ed_record_id_ = record_id;
    ed_record_parton_index_ = parton_index;
    ed_record_pdg_id_ = pdg_id;
    ed_record_parent_index_ = parent_index;
    ed_record_d1_ = d1;
    ed_record_d2_ = d2;
    ed_record_related_index_ = related_index;
    ed_record_is_unresolved_ = is_unresolved ? 1 : 0;
    ed_record_in_medium_ = in_medium ? 1 : 0;
    ed_record_t_ = pos_start[3];
    ed_record_x_ = pos_start[0];
    ed_record_y_ = pos_start[1];
    ed_record_z_ = pos_start[2];
    ed_record_tau_ = properTimeFromPos(pos_start);
    ed_record_px_ = p_start[0];
    ed_record_py_ = p_start[1];
    ed_record_pz_ = p_start[2];
    ed_record_e_ = p_start[3];
    ed_record_t_end_ = pos_end[3];
    ed_record_x_end_ = pos_end[0];
    ed_record_y_end_ = pos_end[1];
    ed_record_z_end_ = pos_end[2];
    ed_record_tau_end_ = properTimeFromPos(pos_end);
    ed_record_px_end_ = p_end[0];
    ed_record_py_end_ = p_end[1];
    ed_record_pz_end_ = p_end[2];
    ed_record_e_end_ = p_end[3];
    ed_record_qperp_ = qperp;
    ed_record_temperature_ = temperature;
    ed_record_length_ = length;
    ed_record_tlength_ = tlength;
    ed_record_type_ = record_type;
    ed_record_label_ = label;
    event_display_detail_tree_->Fill();
#else
    (void)event_id; (void)record_id; (void)parton_index; (void)pdg_id;
    (void)parent_index; (void)d1; (void)d2; (void)related_index;
    (void)is_unresolved; (void)in_medium; (void)pos_start; (void)p_start;
    (void)pos_end; (void)p_end; (void)qperp; (void)temperature;
    (void)length; (void)tlength; (void)record_type; (void)label;
#endif
}

void EnergyLoss::do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched,
                          double x, double y, std::vector<Quench> *recoiled) {
    if (recoiled != nullptr) {
        recoiled->clear();
    }
    // Finite LRES must own the timeline. If elastic is also enabled, Moliere is
    // applied to the currently resolved effective object inside this path.
    if (do_lres_ && lres_rpower_ < 1.e9) {
        do_lres_eloss_impl(partons, quenched, x, y, recoiled);
        return;
    }
    if (do_elastic_) {
        const int event_id = do_event_display_ ? history_event_counter_++ : -1;
        int event_display_segment_id = 0;
        int event_display_record_id = 0;
        std::vector<std::array<double,4>> p_before;
        std::vector<std::array<double,4>> pos_before;
        if (do_event_display_) {
            for (size_t ip = 0; ip < partons.size(); ++ip) {
                fill_event_display_record(event_id, event_display_record_id++,
                                          static_cast<int>(ip), partons[ip].GetId(),
                                          partons[ip].GetMom(), partons[ip].GetD1(), partons[ip].GetD2(),
                                          -1, false, false,
                                          partons[ip].GetRi(), partons[ip].vGetP(),
                                          partons[ip].GetRi(), partons[ip].vGetP(),
                                          0., 0., 0., 0.,
                                          "pythia_splitting_node", "shower_record");
            }
            p_before.reserve(quenched.size());
            pos_before.reserve(quenched.size());
            for (const auto &q : quenched) {
                p_before.push_back(q.vGetP());
                pos_before.push_back(q.GetRi());
            }
        }
        auto callback_factory =
            [&](int parton_index, int pdg_id, int parent_index, int d1, int d2) {
                moliere::ScatteringCallback scattering_callback =
                    [&, parton_index, pdg_id, parent_index, d1, d2]
                    (const moliere::ScatteringCandidate &candidate) {
                        fill_event_display_record(event_id, event_display_record_id++,
                                                  parton_index, pdg_id, parent_index, d1, d2,
                                                  -1, false, true,
                                                  candidate.pos, candidate.p_before,
                                                  candidate.pos, candidate.p_after,
                                                  candidate.qperp, 0., 0., 0.,
                                                  "moliere_scattering",
                                                  "accepted_elastic_scattering");
                        fill_event_display_record(event_id, event_display_record_id++,
                                                  -1, candidate.recoiler_id, parton_index, -1, -1,
                                                  parton_index, false, true,
                                                  candidate.pos, candidate.recoiler_p,
                                                  candidate.pos, candidate.recoiler_p,
                                                  0., 0., 0., 0.,
                                                  "medium_response",
                                                  "moliere_recoiler");
                        fill_event_display_record(event_id, event_display_record_id++,
                                                  -1, candidate.hole_id, parton_index, -1, -1,
                                                  parton_index, false, true,
                                                  candidate.pos, candidate.hole_p,
                                                  candidate.pos, candidate.hole_p,
                                                  0., 0., 0., 0.,
                                                  "medium_response",
                                                  "moliere_hole");
                        return moliere::ScatteringDecision::Apply;
                    };
                moliere::PropagationStepCallback step_callback =
                    [&, parton_index, pdg_id, parent_index, d1, d2]
                    (const moliere::PropagationStep &step) {
                        fill_event_display_record(event_id, event_display_record_id++,
                                                  parton_index, pdg_id, parent_index, d1, d2,
                                                  -1, false, step.in_medium != 0,
                                                  step.pos_before, step.p_before,
                                                  step.pos_after, step.p_after,
                                                  std::sqrt((step.p_after[0] - step.p_before[0]) *
                                                            (step.p_after[0] - step.p_before[0]) +
                                                            (step.p_after[1] - step.p_before[1]) *
                                                            (step.p_after[1] - step.p_before[1])),
                                                  step.temperature, step.step, 0.,
                                                  "moliere_propagation_step",
                                                  "moliere_internal_step");
                    };
                return std::make_pair(scattering_callback, step_callback);
            };
        if (recoiled == nullptr) {
            std::vector<Quench> local_recoiled;
            moliere::do_eloss(partons, quenched, x, y, nr_, kappa_, alpha_, tmethod_, mode_,
                              ebe_hydro_, compat_moliere_legacy_hydro_, hydro_profile_, local_recoiled,
                              do_event_display_ ? callback_factory : moliere::PartonCallbackFactory());
        } else {
            moliere::do_eloss(partons, quenched, x, y, nr_, kappa_, alpha_, tmethod_, mode_,
                              ebe_hydro_, compat_moliere_legacy_hydro_, hydro_profile_, *recoiled,
                              do_event_display_ ? callback_factory : moliere::PartonCallbackFactory());
        }
        if (do_event_display_) {
            for (size_t i = 0; i < quenched.size() && i < p_before.size() && i < partons.size(); ++i) {
                if (quenched[i].GetOrig() == "rem") continue;
                if (!quenched[i].GetIsDone()) continue;
                fill_event_display_segment(event_id, event_display_segment_id++,
                                           static_cast<int>(i), partons[i], quenched[i].GetMom(),
                                           quenched[i].GetD1(), quenched[i].GetD2(),
                                           false, quenched[i].hadScattering(),
                                           pos_before[i], p_before[i],
                                           quenched[i].GetRf(), quenched[i].vGetP(),
                                           0.0, 0.0, "moliere_elastic_segment");
            }
        }
        return;
    }
    do_eloss_impl(partons, quenched, x, y);
}

void EnergyLoss::do_eloss_impl(const std::vector<Parton> &partons, std::vector<Quench> &quenched, double xcre, double ycre) {
    const int event_id = do_event_display_ ? history_event_counter_++ : -1;
    int event_display_segment_id = 0;
    int event_display_record_id = 0;

    if (do_event_display_) {
        for (size_t ip = 0; ip < partons.size(); ++ip) {
            fill_event_display_record(event_id, event_display_record_id++,
                                      static_cast<int>(ip), partons[ip].GetId(),
                                      partons[ip].GetMom(), partons[ip].GetD1(), partons[ip].GetD2(),
                                      -1, false, false,
                                      partons[ip].GetRi(), partons[ip].vGetP(),
                                      partons[ip].GetRi(), partons[ip].vGetP(),
                                      0., 0., 0., 0.,
                                      "pythia_splitting_node", "shower_record");
        }
    }

    // Tag final particles
    std::vector<int> FinId;
    for (size_t i = 0; i < quenched.size(); i++) {
        // Exclude remnants
        if (quenched[i].GetD1() == -1 && quenched[i].GetOrig() != "rem") {
            FinId.push_back(i);
        }
    }

    // Energy Loss Loop: select a final particle, find its oldest undone parent, climb down the family chain. Iterate until all final particles are done
    for (size_t i = 0; i < FinId.size(); i++) {
        int ind = FinId[i];  // Start loop with final particle
        std::vector<int> Fam;  // Family chain array
        Fam.push_back(ind);
        double inhe = 1.;    // Used to see whether mother was completely quenched
        // Family chain loop
        int found = 0;
        do {
            int mom = quenched[ind].GetMom();
            // If "ind" is not a parent parton
            if (mom != -1) {
                // If mother is done
                if (quenched[mom].GetIsDone() == true) {
                    inhe = quenched[mom].vGetP()[3];
                    // End of family chain
                    found = 1;
                }
                // If it is not done
                else {
                    Fam.push_back(mom);
                    ind = mom;
                }
            }
            // If "ind" is a parent parton
            else {
                quenched[ind].SetRi(xcre, ycre, 0., 0.);
                // End of family chain
                found = 1;
            }
        } while (found == 0);

        // Apply energy loss chronologically
        for (size_t w = Fam.size(); w > 0; w--) {
            int tp = Fam[w - 1];
            // If first done mother was totally quenched, set all descendance quenched and done, and exit loop
            if (inhe == 0.) {
                for (size_t j = w; j > 0; j--) {
                    tp = Fam[j - 1];
                    quenched[tp].SetP(0., 0., 0., 0.);
                    quenched[tp].SetIsDone(true);
                }
                break;
            }
            auto p = quenched[tp].vGetP();
            double q = quenched[tp].GetQ();
            auto pos = quenched[tp].GetRi();
            const auto p_before = p;
            const auto pos_before = pos;
            // Time of flight (from formation time argument)
            double tof = 0.2 * 2. * p[3] / (q * q); // in fm
            // If final particle, fly arbitrarily far
            if (w == 1) tof = 10000000000.;
            double length = 0.; // length in QGP
            double tlength = 0.; // temperature weighted length in QGP
            // If colored particle
            if (abs(quenched[tp].GetId()) <= 6 || quenched[tp].GetId() == 21) {
                loss_rate(p, pos, tof, quenched[tp].GetId(), length, tlength,
                          event_id, &event_display_record_id, tp, quenched[tp].GetMom(),
                          quenched[tp].GetD1(), quenched[tp].GetD2(), false);
            } else {
                // If not colored particle, don't do energy loss, but propagate position and time manually
                pos += p / p[3] * tof;
                fill_event_display_record(event_id, event_display_record_id++,
                                          tp, quenched[tp].GetId(), quenched[tp].GetMom(),
                                          quenched[tp].GetD1(), quenched[tp].GetD2(),
                                          -1, false, false,
                                          pos_before, p_before, pos, p,
                                          0., 0., 0., 0.,
                                          "free_stream_step", "legacy_non_colored_free_stream");
            }
            // Update mother momenta and set positions to final
            quenched[tp].vSetP(p);
            quenched[tp].vSetRf(pos);
            quenched[tp].SetIsDone(true);
            quenched[tp].AddLength(length, tlength);
            fill_event_display_segment(event_id, event_display_segment_id++,
                                       tp, partons[tp], quenched[tp].GetMom(),
                                       quenched[tp].GetD1(), quenched[tp].GetD2(),
                                       false, quenched[tp].hadScattering(),
                                       pos_before, p_before, pos, p,
                                       length, tlength,
                                       isColored(quenched[tp].GetId())
                                           ? "legacy_energy_loss_segment"
                                           : "legacy_free_stream_segment");
            // If it got fully quenched, quenched descendance and exit
            if (p[3] == 0.) {
                for (size_t j = w; j > 0; j--) {
                    tp = Fam[j - 1];
                    quenched[tp].SetP(0., 0., 0., 0.);
                    quenched[tp].SetIsDone(true);
                }
                break;
            }
            // If not final particle, propagate quenching to son in chain, and store results for other son
            if (w != 1) {
                // Find two daughters
                int d1 = Fam[w - 2];
                int d2;
                if (quenched[tp].GetD1() == d1) d2 = quenched[tp].GetD2();
                else d2 = quenched[tp].GetD1();
                // Find new momenta for sons: rotate and quench tri-momentum, quench energy
                auto m_p = partons[tp].vGetP();
                auto d1_p = partons[d1].vGetP();
                auto d2_p = partons[d2].vGetP();
                quenched_sons(m_p, p, d1_p, d2_p);
                // Propagate momenta and positions
                quenched[d1].vSetP(d1_p);
                quenched[d1].vSetInhP(d1_p);
                quenched[d1].vSetRi(pos);
                quenched[d1].AddLength(length, tlength);
                quenched[d2].vSetP(d2_p);
                quenched[d2].vSetInhP(d2_p);
                quenched[d2].vSetRi(pos);
                quenched[d2].AddLength(length, tlength);
            }
        }
        Fam.clear();
    }
    FinId.clear();
}

void EnergyLoss::do_lres_eloss_impl(const std::vector<Parton> &partons, std::vector<Quench> &quenched,
                                    double xcre, double ycre, std::vector<Quench> *recoiled) {
    const size_t n = quenched.size();
    if (n == 0) return;

    std::vector<int> final_ids;
    final_ids.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        if (quenched[i].GetD1() == -1 && quenched[i].GetOrig() != "rem") {
            final_ids.push_back(static_cast<int>(i));
        }
    }

    std::vector<LresLifetime> life(n);
    std::vector<bool> timeline_done(n, false);

    // First build the vacuum formation timeline used by the finite-resolution rule.
    for (int final_id : final_ids) {
        int ind = final_id;
        std::vector<int> family;
        family.push_back(ind);

        double x = xcre;
        double y = ycre;
        double z = 0.;
        double t = 0.;

        while (true) {
            const int mom = quenched[ind].GetMom();
            if (mom >= 0 && mom < static_cast<int>(n)) {
                if (timeline_done[mom]) {
                    x = life[mom].rf[0];
                    y = life[mom].rf[1];
                    z = life[mom].rf[2];
                    t = life[mom].finish;
                    break;
                }
                family.push_back(mom);
                ind = mom;
                continue;
            }
            break;
        }

        for (auto it = family.rbegin(); it != family.rend(); ++it) {
            const int idx = *it;
            const auto p = partons[idx].vGetP();
            const double e = p[3];
            const double tof = (idx == final_id) ? kLresFinalFlightTime : safeFormationTime(partons[idx]);
            const double previous = t;

            life[idx].creation = previous;
            life[idx].ri = {x, y, z, previous};

            if (e != 0.) {
                x += p[0] / e * tof;
                y += p[1] / e * tof;
                z += p[2] / e * tof;
            }
            t = previous + tof;

            life[idx].finish = t;
            life[idx].rf = {x, y, z, t};
            timeline_done[idx] = true;
        }
    }

    std::vector<int> brother(n, -1);
    std::vector<int> effective_mom(n, -1);
    for (size_t i = 0; i < n; ++i) {
        effective_mom[i] = quenched[i].GetMom();
    }

    for (size_t mom = 0; mom < n; ++mom) {
        int first = -1;
        int second = -1;
        for (size_t child = 0; child < n; ++child) {
            if (quenched[child].GetMom() != static_cast<int>(mom)) continue;
            if (first == -1) first = static_cast<int>(child);
            else {
                second = static_cast<int>(child);
                break;
            }
        }
        if (first != -1 && second != -1) {
            brother[first] = second;
            brother[second] = first;
        }
    }

    for (size_t i = 0; i < n; ++i) {
        const int mom = quenched[i].GetMom();
        const int sib = brother[i];
        if (sib < 0 || mom < 0 || mom >= static_cast<int>(n)) continue;

        const auto parent_p = partons[mom].vGetP();
        const auto p_i = partons[i].vGetP();
        const auto p_s = partons[sib].vGetP();
        const auto v_i = velocity(p_i);
        const auto v_s = velocity(p_s);

        const double store = resolution_time(parent_p[3], parent_p[0], parent_p[1], parent_p[2],
                                             life[i].ri[0], life[i].ri[1], life[i].ri[2],
                                             v_s[0] - v_i[0], v_s[1] - v_i[1], v_s[2] - v_i[2],
                                             life[i].creation);
        life[i].resolve = store;
        life[i].resolve_abs = life[i].creation + store;
        life[i].effective_time = life[i].resolve_abs;
    }

    // Build the finite-LRES effective shower tree.  A daughter pair should not
    // become visible to the medium before an unresolved ancestor is visible, and
    // the two daughters of one splitting must become visible together.  The
    // iteration below pulls ancestor/sibling effective times earlier until this
    // ordered medium-resolution tree is self-consistent.
    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < n; ++i) {
            const int mom = effective_mom[i];
            if (mom >= 0 && mom < static_cast<int>(n) && effective_mom[mom] >= 0) {
                if (life[i].effective_time < life[mom].effective_time) {
                    life[mom].effective_time = life[i].effective_time;
                    changed = true;
                }
            }

            const int sib = brother[i];
            if (sib >= 0 && sib < static_cast<int>(n) && mom >= 0 && mom < static_cast<int>(n)) {
                const double shared = std::min(life[i].effective_time, life[sib].effective_time);
                if (life[i].effective_time != shared || life[sib].effective_time != shared) {
                    life[i].effective_time = shared;
                    life[sib].effective_time = shared;
                    changed = true;
                }
            }
        }
    }

    // If a node and its parent have the same effective time, that parent would
    // have zero lifetime as a separate medium object.  Collapse such nodes so
    // propagation is done only for finite-duration effective charges.
    changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < n; ++i) {
            const int mom = effective_mom[i];
            if (mom < 0 || mom >= static_cast<int>(n)) continue;
            const int grand = effective_mom[mom];
            if (grand < 0 || grand >= static_cast<int>(n)) continue;
            if (life[i].effective_time == life[mom].effective_time) {
                effective_mom[i] = grand;
                changed = true;
            }
        }
    }

    // live_time is the duration for which each effective object propagates
    // before the medium resolves its effective daughter.  Final leaves are
    // handled later with a long final-flight time.
    std::vector<int> effective_daughter(n, -1);
    for (size_t i = 0; i < n; ++i) {
        const int mom = effective_mom[i];
        if (mom >= 0 && mom < static_cast<int>(n)) {
            effective_daughter[mom] = static_cast<int>(i);
        }
    }

    std::vector<double> live_time(n, 0.);
    for (size_t i = 0; i < n; ++i) {
        const int daughter = effective_daughter[i];
        if (daughter >= 0 && daughter < static_cast<int>(n)) {
            live_time[i] = life[daughter].effective_time - life[i].effective_time;
            if (live_time[i] < 0.) live_time[i] = 0.;
        }
    }

    std::vector<LresState> qstate(n);
    std::vector<bool> done(n, false);
    std::vector<bool> modeE_seeded(n, false);
    std::vector<std::array<double,4>> qorient(n);
    std::vector<int> qhad(n, 0);
    for (size_t i = 0; i < n; ++i) {
        qorient[i] = quenched[i].orient();
        qhad[i] = quenched[i].hadScattering();
    }
    std::vector<Quench> lres_moliere_particles;
    std::vector<Quench> local_recoiled;
    std::vector<std::string> history_records;
    const int event_id = history_event_counter_++;
    int event_display_segment_id = 0;
    int event_display_record_id = 0;
    const std::string history_mode =
        do_moliere_recursive_unresolved_resolution_
            ? "recursive_unresolved_resolution"
            : (do_moliere_dynamic_daughter_unresolved_resolution_
                   ? "dynamic_daughter_unresolved_resolution"
                   : (do_moliere_dynamic_unresolved_resolution_
                   ? "dynamic_unresolved_resolution"
                   : (do_moliere_on_unresolved_partons_ ? "individual_unresolved" : "coherent_unresolved")));
    auto emit = [&](const std::string &type, int parton_id, int parent_id, int d1, int d2, double t,
                    const std::array<double,4> &pos, const std::array<double,4> &p, double qperp,
                    const std::string &label, const std::string &note) {
        if (do_event_display_) {
            const int pdg_id =
                (parton_id >= 0 && parton_id < static_cast<int>(partons.size()))
                    ? partons[parton_id].GetId()
                    : 0;
            fill_event_display_record(event_id, event_display_record_id++,
                                      parton_id, pdg_id, parent_id, d1, d2,
                                      -1, type.find("unresolved") != std::string::npos,
                                      false, pos, p, pos, p, qperp, 0., 0., 0.,
                                      type, label + ":" + note);
        }
        if (dump_hybrid_evolution_history_) {
            std::ostringstream os;
            os << std::setprecision(10)
               << event_id << '\t'
               << history_mode << '\t'
               << type << '\t'
               << parton_id << '\t'
               << parent_id << '\t'
               << d1 << '\t'
               << d2 << '\t'
               << t << '\t'
               << pos[0] << '\t' << pos[1] << '\t' << pos[2] << '\t'
               << p[0] << '\t' << p[1] << '\t' << p[2] << '\t' << p[3] << '\t'
               << qperp << '\t'
               << label << '\t'
               << note;
            history_records.push_back(os.str());
        }
    };
    auto make_moliere_scattering_callback =
        [&](int parton_id, int pdg_id, int parent_id, int d1, int d2,
            bool is_unresolved, const std::string &label) {
            return moliere::ScatteringCallback(
                [&, parton_id, pdg_id, parent_id, d1, d2, is_unresolved, label]
                (const moliere::ScatteringCandidate &candidate) {
                    fill_event_display_record(event_id, event_display_record_id++,
                                              parton_id, pdg_id, parent_id, d1, d2,
                                              -1, is_unresolved, true,
                                              candidate.pos, candidate.p_before,
                                              candidate.pos, candidate.p_after,
                                              candidate.qperp, 0., 0., 0.,
                                              "moliere_scattering", label);
                    fill_event_display_record(event_id, event_display_record_id++,
                                              -1, candidate.recoiler_id, parton_id, -1, -1,
                                              parton_id, is_unresolved, true,
                                              candidate.pos, candidate.recoiler_p,
                                              candidate.pos, candidate.recoiler_p,
                                              0., 0., 0., 0.,
                                              "medium_response", "moliere_recoiler");
                    fill_event_display_record(event_id, event_display_record_id++,
                                              -1, candidate.hole_id, parton_id, -1, -1,
                                              parton_id, is_unresolved, true,
                                              candidate.pos, candidate.hole_p,
                                              candidate.pos, candidate.hole_p,
                                              0., 0., 0., 0.,
                                              "medium_response", "moliere_hole");
                    return moliere::ScatteringDecision::Apply;
                });
        };
    auto make_moliere_step_callback =
        [&](int parton_id, int pdg_id, int parent_id, int d1, int d2,
            bool is_unresolved, const std::string &label) {
            return moliere::PropagationStepCallback(
                [&, parton_id, pdg_id, parent_id, d1, d2, is_unresolved, label]
                (const moliere::PropagationStep &step) {
                    fill_event_display_record(event_id, event_display_record_id++,
                                              parton_id, pdg_id, parent_id, d1, d2,
                                              -1, is_unresolved, step.in_medium != 0,
                                              step.pos_before, step.p_before,
                                              step.pos_after, step.p_after,
                                              std::sqrt((step.p_after[0] - step.p_before[0]) *
                                                        (step.p_after[0] - step.p_before[0]) +
                                                        (step.p_after[1] - step.p_before[1]) *
                                                        (step.p_after[1] - step.p_before[1])),
                                              step.temperature, step.step, 0.,
                                              "moliere_propagation_step", label);
                });
        };

    for (size_t i = 0; i < n; ++i) {
        const int mom = quenched[i].GetMom();
        emit("split", static_cast<int>(i), mom, quenched[i].GetD1(), quenched[i].GetD2(),
             life[i].creation, life[i].ri, partons[i].vGetP(), 0.0, "shower_branch", "formation_time");
        if (mom >= 0 && mom < static_cast<int>(n)) {
            const int sib = brother[i];
            if (sib >= 0 && sib < static_cast<int>(n)) {
                const double dr = deltaPhiXYFromP(partons[i].vGetP(), partons[sib].vGetP());
                emit("opening_angle", static_cast<int>(i), mom, static_cast<int>(i), sib,
                     life[i].creation, life[i].ri, partons[i].vGetP(), dr, "deltaPhi_xy", "at_creation");
            }
        }
    }

    for (size_t i = 0; i < n; ++i) {
        const int d1 = quenched[i].GetD1();
        const int d2 = quenched[i].GetD2();
        const bool has_unresolved_pair =
            d1 >= 0 && d2 >= 0 &&
            d1 < static_cast<int>(n) && d2 < static_cast<int>(n) &&
            effective_mom[d1] == static_cast<int>(i) && effective_mom[d2] == static_cast<int>(i) &&
            life[d1].effective_time > life[i].effective_time &&
            life[d2].effective_time > life[i].effective_time;
        if (!has_unresolved_pair) continue;
        const double t0 = life[i].effective_time;
        const double t1 = std::min(life[d1].effective_time, life[d2].effective_time);
        if (t1 <= t0) continue;
        emit("unresolved_region_start", static_cast<int>(i), quenched[i].GetMom(), d1, d2,
             t0, life[i].ri, partons[i].vGetP(), 0.0, "deltaR<Lres", "unresolved_by_Lres");
        emit("resolution", static_cast<int>(i), quenched[i].GetMom(), d1, d2,
             t1, life[i].rf, partons[i].vGetP(), 0.0, "deltaR>=Lres", "resolved_by_medium");
    }

    for (int final_id : final_ids) {
        int ind = final_id;
        std::vector<int> family;
        family.push_back(ind);
        bool mother_zero = false;

        while (true) {
            if (modeE_seeded[ind] && !done[ind]) {
                break;
            }
            const int mom = effective_mom[ind];
            if (mom >= 0 && mom < static_cast<int>(n)) {
                if (done[mom]) {
                    const double parent_e = partons[mom].vGetP()[3];
                    const double frac = (parent_e != 0.) ? qstate[mom].p[3] / parent_e : 0.;
                    qstate[ind].p = partons[ind].vGetP() * frac;
                    if (do_elastic_) {
                        qorient[ind] = orientationFor(qstate[ind].p);
                        if (qhad[mom] == 1 || qhad[mom] == 2) qhad[ind] = 2;
                    }
                    const auto p = partons[ind].vGetP();
                    const double e = p[3];
                    const double dt = life[ind].effective_time - life[ind].creation;
                    qstate[ind].r = life[ind].ri;
                    if (e != 0.) {
                        qstate[ind].r[0] += p[0] / e * dt;
                        qstate[ind].r[1] += p[1] / e * dt;
                        qstate[ind].r[2] += p[2] / e * dt;
                    }
                    qstate[ind].r[3] = life[ind].effective_time;
                    if (qstate[mom].p[3] == 0.) mother_zero = true;
                    break;
                }
                family.push_back(mom);
                ind = mom;
                continue;
            }

            qstate[ind].p = partons[ind].vGetP();
            qstate[ind].r = {xcre, ycre, 0., life[ind].effective_time};
            if (do_elastic_) {
                qorient[ind] = orientationFor(qstate[ind].p);
            }
            break;
        }

        for (auto it = family.rbegin(); it != family.rend(); ++it) {
            const int idx = *it;
            if (mother_zero) {
                for (auto zero_it = it; zero_it != family.rend(); ++zero_it) {
                    qstate[*zero_it].p = {0., 0., 0., 0.};
                    done[*zero_it] = true;
                }
                break;
            }

            auto p = qstate[idx].p;
            const auto vac_p = partons[idx].vGetP();
            const double vac_e = vac_p[3];
            std::array<double,4> pos = qstate[idx].r;
            if (!modeE_seeded[idx]) {
                const double dt = life[idx].effective_time - life[idx].creation;
                pos = life[idx].ri;
                // Start the effective segment at the vacuum shower location
                // corresponding to its effective resolution time.  The momentum p
                // carries the accumulated quenching inherited from its ancestor.
                if (vac_e != 0.) {
                    pos[0] += vac_p[0] / vac_e * dt;
                    pos[1] += vac_p[1] / vac_e * dt;
                    pos[2] += vac_p[2] / vac_e * dt;
                }
                pos[3] = life[idx].effective_time;
            }
            if (pos[3] * pos[3] < pos[2] * pos[2]) {
                pos[3] = std::abs(pos[2]) + 1.e-9;
            }
            double tof = live_time[idx];
            if (modeE_seeded[idx] && idx != final_id) {
                const double segment_end = life[idx].effective_time + live_time[idx];
                tof = std::max(0., segment_end - pos[3]);
            }
            if (idx == final_id) tof = kLresFinalFlightTime;
            modeE_seeded[idx] = false;

            double length = 0.;
            double tlength = 0.;
            const auto p_segment_start = p;
            const auto pos_segment_start = pos;
            const int d1_segment = quenched[idx].GetD1();
            const int d2_segment = quenched[idx].GetD2();
            const bool segment_has_unresolved_pair =
                d1_segment >= 0 && d2_segment >= 0 &&
                d1_segment < static_cast<int>(n) && d2_segment < static_cast<int>(n) &&
                effective_mom[d1_segment] == idx && effective_mom[d2_segment] == idx &&
                life[d1_segment].effective_time > life[idx].effective_time &&
                life[d2_segment].effective_time > life[idx].effective_time;
            const bool modee_unresolved_segment =
                do_elastic_ && do_moliere_recursive_unresolved_resolution_ &&
                segment_has_unresolved_pair;
            std::string segment_type = "lres_zero_momentum_segment";
            // Mode E must still inspect colored descendants of an unresolved
            // color-neutral parent (for example gamma -> q qbar).  Such a parent
            // is never used as a coherent Moliere source below.
            if (isColored(partons[idx].GetId()) || modee_unresolved_segment) {
                if (do_elastic_) {
                    const auto p_before = p;
                    const auto pos_before = pos;
                    const int d1 = d1_segment;
                    const int d2 = d2_segment;
                    const bool has_unresolved_pair = segment_has_unresolved_pair;

                    // Unresolved Moliere mode precedence:
                    //   E: recursive daughter-candidate dipole resolution through nested LRES trees;
                    //   D: daughter-candidate dynamic dipole resolution for one active sibling pair;
                    //   C: parent-candidate dynamic dipole resolution;
                    //   B: propagate unresolved daughters independently;
                    //   A: default coherent unresolved parent.
                    if (do_moliere_recursive_unresolved_resolution_ && has_unresolved_pair) {
                        // Mode E recursively tests the already-formed colored frontier below each
                        // active coherent object. A daughter candidate is walked upward through its
                        // enclosing sibling dipoles. If it resolves one level, it acts on that branch.
                        // If it resolves none, the daughter proposal is vetoed completely and the
                        // active colored coherent object is sampled with its own identity and rate.
                        // Parent proposals are accepted only when they remain unresolved; resolving
                        // parent proposals are vetoed and the parent stream continues. This keeps one
                        // recoil/hole source per accepted scattering and forbids coherent kicks on a
                        // color-neutral parent. Shower-formation times split the search into ordered
                        // windows, while the finite-LRES timeline itself remains unchanged.
                        segment_type = "modeE_recursive_unresolved";
                        ++n_unresolved_segments_dynamic_;
                        ++n_unresolved_segments_recursive_;
                        const double total_end = pos[3] + tof;
                        const int modee_segment_id = event_display_segment_id;
                        int modee_iteration = 0;

                        std::vector<int> active_groups{idx};
                        qstate[idx].p = p;
                        qstate[idx].r = pos;

                        auto child_pair = [&](int node, int &c1, int &c2) {
                            c1 = (node >= 0 && node < static_cast<int>(n)) ? quenched[node].GetD1() : -1;
                            c2 = (node >= 0 && node < static_cast<int>(n)) ? quenched[node].GetD2() : -1;
                            return c1 >= 0 && c2 >= 0 &&
                                   c1 < static_cast<int>(n) && c2 < static_cast<int>(n);
                        };
                        auto is_active = [&](int node) {
                            return std::find(active_groups.begin(), active_groups.end(), node) != active_groups.end();
                        };
                        auto descends_from = [&](int node, int ancestor) {
                            int cur = node;
                            while (cur >= 0 && cur < static_cast<int>(n)) {
                                if (cur == ancestor) return true;
                                cur = quenched[cur].GetMom();
                            }
                            return false;
                        };
                        auto split_child_momenta = [&](int parent, const std::array<double,4> &parent_p,
                                                       int c1, int c2,
                                                       std::array<double,4> &p1,
                                                       std::array<double,4> &p2) {
                            mapChildMomentaFromLiveParent(partons[parent].vGetP(), parent_p,
                                                           partons[c1].vGetP(), partons[c2].vGetP(),
                                                           p1, p2);
                        };
                        auto split_child_positions = [&](int c1, int c2,
                                                         const std::array<double,4> &parent_p,
                                                         const std::array<double,4> &parent_pos,
                                                         const std::array<double,4> &p1,
                                                         const std::array<double,4> &p2,
                                                         std::array<double,4> &pos1,
                                                         std::array<double,4> &pos2) {
                            const double target_time = parent_pos[3];
                            pos1 = separatedChildPosition(parent_pos, parent_p, p1,
                                                          life[c1].creation, target_time);
                            pos2 = separatedChildPosition(parent_pos, parent_p, p2,
                                                          life[c2].creation, target_time);
                        };
                        auto split_active_parent = [&](int parent, const std::string &note) {
                            if (!is_active(parent)) return false;
                            int c1 = -1;
                            int c2 = -1;
                            if (!child_pair(parent, c1, c2)) return false;
                            std::array<double,4> p1;
                            std::array<double,4> p2;
                            split_child_momenta(parent, qstate[parent].p, c1, c2, p1, p2);
                            const double spatial_residual =
                                spatialMomentumResidual(qstate[parent].p, p1, p2);
                            const double energy_residual =
                                energyResidual(qstate[parent].p, p1, p2);
                            const double relative_spatial_residual =
                                qstate[parent].p[3] > 0. ? spatial_residual / qstate[parent].p[3] : 0.;
                            ++n_recursive_opening_closure_checks_;
                            sum_recursive_opening_spatial_residual_ += spatial_residual;
                            max_recursive_opening_spatial_residual_ =
                                std::max(max_recursive_opening_spatial_residual_, spatial_residual);
                            sum_recursive_opening_energy_residual_ += energy_residual;
                            max_recursive_opening_energy_residual_ =
                                std::max(max_recursive_opening_energy_residual_, energy_residual);
                            std::ostringstream closure_note;
                            closure_note << note
                                         << ":spatial_abs=" << spatial_residual
                                         << ":energy_abs=" << energy_residual;
                            emit("recursive_opening_closure", parent, quenched[parent].GetMom(), c1, c2,
                                 qstate[parent].r[3], qstate[parent].r, qstate[parent].p,
                                 relative_spatial_residual, "relative_spatial_residual",
                                 closure_note.str());
                            std::array<double,4> pos1;
                            std::array<double,4> pos2;
                            split_child_positions(c1, c2, qstate[parent].p, qstate[parent].r,
                                                  p1, p2, pos1, pos2);
                            qstate[c1].p = p1;
                            qstate[c2].p = p2;
                            qstate[c1].r = pos1;
                            qstate[c2].r = pos2;
                            qorient[c1] = orientationFor(qstate[c1].p);
                            qorient[c2] = orientationFor(qstate[c2].p);
                            if (qhad[parent] == 1 || qhad[parent] == 2) {
                                qhad[c1] = 2;
                                qhad[c2] = 2;
                            }
                            active_groups.erase(std::remove(active_groups.begin(), active_groups.end(), parent),
                                                active_groups.end());
                            active_groups.push_back(c1);
                            active_groups.push_back(c2);
                            ++n_recursive_tree_updates_;
                            emit("recursive_tree_update", parent, quenched[parent].GetMom(), c1, c2,
                                 qstate[parent].r[3], qstate[parent].r, qstate[parent].p,
                                 0., "split_active_group", note);
                            return true;
                        };
                        std::function<bool(int)> ensure_active = [&](int target) {
                            if (is_active(target)) return true;
                            std::vector<int> path;
                            int cur = target;
                            while (cur >= 0 && cur < static_cast<int>(n) && !is_active(cur)) {
                                path.push_back(cur);
                                cur = quenched[cur].GetMom();
                            }
                            if (cur < 0 || cur >= static_cast<int>(n) || !is_active(cur)) return false;
                            for (auto it_path = path.rbegin(); it_path != path.rend(); ++it_path) {
                                const int child = *it_path;
                                const int parent = quenched[child].GetMom();
                                if (!is_active(child)) {
                                    if (!split_active_parent(parent, "modeE_dynamic_decoherence")) return false;
                                }
                            }
                            return is_active(target);
                        };
                        auto find_active_ancestor = [&](int node) {
                            int cur = node;
                            while (cur >= 0 && cur < static_cast<int>(n)) {
                                if (is_active(cur)) return cur;
                                cur = quenched[cur].GetMom();
                            }
                            return -1;
                        };

                        struct RecursiveProjectedState {
                            bool ok = false;
                            std::array<double,4> p = {0., 0., 0., 0.};
                            std::array<double,4> pos = {0., 0., 0., 0.};
                            int had = 0;
                            std::array<double,4> orient = {0., 0., 0., 1.};
                        };
                        auto project_from_ancestor_state =
                            [&](int node, int ancestor,
                                const std::array<double,4> &ancestor_p,
                                const std::array<double,4> &ancestor_pos,
                                int ancestor_had, double target_time) {
                            RecursiveProjectedState state;
                            if (ancestor < 0 || !descends_from(node, ancestor)) return state;
                            std::vector<int> path;
                            int cur = node;
                            while (cur != ancestor && cur >= 0 && cur < static_cast<int>(n)) {
                                path.push_back(cur);
                                cur = quenched[cur].GetMom();
                            }
                            if (cur != ancestor) return state;
                            state.ok = true;
                            state.p = ancestor_p;
                            state.pos = ancestor_pos;
                            if (target_time > state.pos[3] + 1.e-9 && state.p[3] > 0.) {
                                state.pos += causalVelocity(state.p) * (target_time - state.pos[3]);
                            }
                            state.pos[3] = target_time;
                            state.had = qhad[node];
                            for (auto it_path = path.rbegin(); it_path != path.rend(); ++it_path) {
                                const int child = *it_path;
                                const int parent = quenched[child].GetMom();
                                int c1 = -1;
                                int c2 = -1;
                                if (!child_pair(parent, c1, c2)) {
                                    state.ok = false;
                                    return state;
                                }
                                std::array<double,4> p1;
                                std::array<double,4> p2;
                                std::array<double,4> pos1;
                                std::array<double,4> pos2;
                                split_child_momenta(parent, state.p, c1, c2, p1, p2);
                                split_child_positions(c1, c2, state.p, state.pos,
                                                      p1, p2, pos1, pos2);
                                state.p = (child == c1) ? p1 : p2;
                                state.pos = (child == c1) ? pos1 : pos2;
                            }
                            if (ancestor_had == 1 || ancestor_had == 2) state.had = 2;
                            state.orient = orientationFor(state.p);
                            return state;
                        };
                        auto project_from_active = [&](int node, double target_time) {
                            const int ancestor = find_active_ancestor(node);
                            if (ancestor < 0) return RecursiveProjectedState{};
                            return project_from_ancestor_state(
                                node, ancestor, qstate[ancestor].p, qstate[ancestor].r,
                                qhad[ancestor], target_time);
                        };
                        struct ModeEPairDperp {
                            double dperp = 0.;
                            bool used_live_positions = false;
                        };
                        auto modeE_pair_dperp = [&](int child, int sibling, double target_time) {
                            ModeEPairDperp result;
                            const auto child_state = project_from_active(child, target_time);
                            const auto sibling_state = project_from_active(sibling, target_time);
                            if (child_state.ok && sibling_state.ok) {
                                result.dperp = transverseSeparation(child_state.pos, sibling_state.pos);
                                result.used_live_positions = true;
                                return result;
                            }
                            result.dperp = compute_unresolved_pair_dperp(
                                partons[child].vGetP(), partons[sibling].vGetP(),
                                life[child].creation, target_time);
                            return result;
                        };
                        auto modeE_pair_dperp_from_source =
                            [&](int child, int sibling, int ancestor,
                                const std::array<double,4> &source_p,
                                const std::array<double,4> &source_pos,
                                int source_had, double target_time) {
                                ModeEPairDperp result;
                                const auto child_state = project_from_ancestor_state(
                                    child, ancestor, source_p, source_pos,
                                    source_had, target_time);
                                const auto sibling_state = project_from_ancestor_state(
                                    sibling, ancestor, source_p, source_pos,
                                    source_had, target_time);
                                if (child_state.ok && sibling_state.ok) {
                                    result.dperp = transverseSeparation(
                                        child_state.pos, sibling_state.pos);
                                    result.used_live_positions = true;
                                    return result;
                                }
                                result.dperp = compute_unresolved_pair_dperp(
                                    partons[child].vGetP(), partons[sibling].vGetP(),
                                    life[child].creation, target_time);
                                return result;
                            };
                        std::function<void(int, double, std::vector<int>&)> collect_probe_frontier =
                            [&](int node, double t, std::vector<int> &frontier) {
                                int c1 = -1;
                                int c2 = -1;
                                if (child_pair(node, c1, c2) &&
                                    life[c1].creation <= t + 1.e-9 &&
                                    life[c2].creation <= t + 1.e-9) {
                                    collect_probe_frontier(c1, t, frontier);
                                    collect_probe_frontier(c2, t, frontier);
                                    return;
                                }
                                if (isColored(partons[node].GetId()) &&
                                    std::find(frontier.begin(), frontier.end(), node) == frontier.end()) {
                                    frontier.push_back(node);
                                }
                            };
                        std::function<void(int, double, double&)> tighten_next_formation_boundary =
                            [&](int node, double t, double &boundary) {
                                int c1 = -1;
                                int c2 = -1;
                                if (!child_pair(node, c1, c2)) return;
                                for (int child : {c1, c2}) {
                                    const double creation = life[child].creation;
                                    if (creation > t + 1.e-9) {
                                        boundary = std::min(boundary, creation);
                                    } else {
                                        tighten_next_formation_boundary(child, t, boundary);
                                    }
                                }
                            };
                        auto propagate_active_groups = [&](double target_time) {
                            for (int group : active_groups) {
                                const double dt_group = std::max(0., target_time - qstate[group].r[3]);
                                if (dt_group <= 1.e-9 || qstate[group].p[3] <= 0.) continue;
                                int gd1 = -1;
                                int gd2 = -1;
                                child_pair(group, gd1, gd2);
                                if (isColored(partons[group].GetId())) {
                                    loss_rate(qstate[group].p, qstate[group].r, dt_group,
                                              partons[group].GetId(), length, tlength,
                                              event_id, &event_display_record_id, group,
                                              quenched[group].GetMom(), gd1, gd2, true);
                                } else {
                                    qstate[group].r += qstate[group].p / qstate[group].p[3] * dt_group;
                                }
                            }
                        };
                        auto recombine_active_groups = [&]() {
                            std::array<double,4> sum_p = {0., 0., 0., 0.};
                            std::array<double,4> sum_pos = {0., 0., 0., 0.};
                            for (int group : active_groups) sum_p += qstate[group].p;
                            if (sum_p[3] > 0.) {
                                for (int group : active_groups) {
                                    sum_pos[0] += qstate[group].p[3] * qstate[group].r[0];
                                    sum_pos[1] += qstate[group].p[3] * qstate[group].r[1];
                                    sum_pos[2] += qstate[group].p[3] * qstate[group].r[2];
                                }
                                sum_pos[0] /= sum_p[3];
                                sum_pos[1] /= sum_p[3];
                                sum_pos[2] /= sum_p[3];
                            }
                            double max_t = pos[3];
                            for (int group : active_groups) max_t = std::max(max_t, qstate[group].r[3]);
                            sum_pos[3] = max_t;
                            p = sum_p;
                            pos = sum_pos;
                            qstate[idx].p = p;
                            qstate[idx].r = pos;
                            qorient[idx] = orientationFor(p);
                            for (int group : active_groups) {
                                if (qhad[group] == 1 || qhad[group] == 2) qhad[idx] = 2;
                            }
                        };

                        struct RecursiveProbe {
                            bool found = false;
                            int probe = -1;
                            int active_ancestor = -1;
                            moliere::ScatteringCandidate candidate;
                            std::array<double,4> p = {0., 0., 0., 0.};
                            std::array<double,4> pos = {0., 0., 0., 0.};
                            int had = 0;
                            std::array<double,4> orient = {0., 0., 0., 1.};
                            numrand rng;
                        };
                        auto probe_object = [&](int probe_id, const RecursiveProjectedState &start,
                                                double probe_end, numrand rng_start,
                                                int elastic_seed) {
                            RecursiveProbe probe;
                            probe.probe = probe_id;
                            probe.active_ancestor = find_active_ancestor(probe_id);
                            probe.p = start.p;
                            probe.pos = start.pos;
                            probe.had = start.had;
                            probe.orient = start.orient;
                            probe.rng = rng_start;
                            const double remaining = std::max(0., probe_end - probe.pos[3]);
                            auto callback = [&](const moliere::ScatteringCandidate &candidate) {
                                probe.found = true;
                                probe.candidate = candidate;
                                return moliere::ScatteringDecision::StopBeforeApply;
                            };
                            // Moliere gen_particles uses the global Distributions.hpp generator,
                            // so branch-local probes must scope that generator in addition to numrand.
                            const auto elastic_rng_state = moliere::elastic_generator_state();
                            moliere::seed_elastic_generator(static_cast<unsigned int>(elastic_seed));
                            if (remaining > 1.e-9 && probe.p[3] > 0.) {
                                moliere::propagate_segment_with_scattering_callback(
                                    probe.p, probe.pos, remaining, partons[probe_id].GetId(),
                                    probe.rng, kappa_, alpha_, tmethod_, mode_, ebe_hydro_,
                                    compat_moliere_legacy_hydro_, hydro_profile_,
                                    lres_moliere_particles, probe.had, probe.orient, callback);
                            }
                            moliere::set_elastic_generator_state(elastic_rng_state);
                            return probe;
                        };
                        std::function<void(int, double, std::vector<std::pair<int,int>>&)> collect_formed_pairs =
                            [&](int node, double t, std::vector<std::pair<int,int>> &pairs) {
                                int c1 = -1;
                                int c2 = -1;
                                if (!child_pair(node, c1, c2)) return;
                                if (life[c1].creation > t + 1.e-9 ||
                                    life[c2].creation > t + 1.e-9) {
                                    return;
                                }
                                pairs.emplace_back(c1, c2);
                                collect_formed_pairs(c1, t, pairs);
                                collect_formed_pairs(c2, t, pairs);
                            };

                        auto resample_coherent_source = [&](int active_ancestor, double probe_end,
                                                            int iteration_id) {
                            RecursiveProbe probe;
                            probe.probe = active_ancestor;
                            probe.active_ancestor = active_ancestor;
                            probe.p = qstate[active_ancestor].p;
                            probe.pos = qstate[active_ancestor].r;
                            probe.had = qhad[active_ancestor];
                            probe.orient = qorient[active_ancestor];

                            ++n_recursive_coherent_resample_requests_;
                            if (!isColored(partons[active_ancestor].GetId())) {
                                ++n_recursive_color_neutral_parent_skips_;
                                emit("recursive_coherent_resample_skip", active_ancestor,
                                     quenched[active_ancestor].GetMom(),
                                     quenched[active_ancestor].GetD1(),
                                     quenched[active_ancestor].GetD2(), probe.pos[3],
                                     probe.pos, probe.p, 0., "color_neutral_parent",
                                     "modeE_no_coherent_moliere_source");
                                return probe;
                            }

                            const int coherent_seed = modeECoherentSourceSeed(
                                nr_.GetIr(), event_id, modee_segment_id, iteration_id,
                                active_ancestor);
                            probe.rng = numrand(coherent_seed);
                            const double remaining = std::max(0., probe_end - probe.pos[3]);
                            auto callback = [&](const moliere::ScatteringCandidate &candidate) {
                                ++n_recursive_coherent_resample_candidates_;
                                ++n_unresolved_candidate_scatters_;
                                std::vector<std::pair<int,int>> pairs;
                                collect_formed_pairs(active_ancestor, candidate.pos[3], pairs);

                                bool resolves_any = false;
                                double max_qd = 0.;
                                for (const auto &pair : pairs) {
                                    // The coherent probe has already accumulated all
                                    // continuous updates up to this candidate. Project
                                    // its daughters from that candidate-time state, not
                                    // from the stale state at the failed daughter probe.
                                    const auto dperp_result = modeE_pair_dperp_from_source(
                                        pair.first, pair.second, active_ancestor,
                                        candidate.p_before, candidate.pos, probe.had,
                                        candidate.pos[3]);
                                    if (dperp_result.used_live_positions) ++n_recursive_live_dperp_tests_;
                                    else ++n_recursive_vacuum_dperp_fallbacks_;
                                    const double qd = candidate.qperp * dperp_result.dperp;
                                    max_qd = std::max(max_qd, qd);
                                    resolves_any = resolves_any ||
                                        passes_dynamic_moliere_resolution_test(
                                            candidate.qperp, dperp_result.dperp,
                                            moliere_unresolved_resolution_c_);
                                }
                                sum_qperp_dperp_unresolved_candidates_ += max_qd;

                                std::ostringstream test_note;
                                test_note << (resolves_any
                                                  ? "modeE_parent_candidate_veto_resolving"
                                                  : "modeE_parent_candidate_accept_coherent")
                                          << ":candidate_qperp=" << candidate.qperp;
                                emit("recursive_coherent_resample_test", active_ancestor,
                                     quenched[active_ancestor].GetMom(),
                                     quenched[active_ancestor].GetD1(),
                                     quenched[active_ancestor].GetD2(), candidate.pos[3],
                                     candidate.pos, candidate.p_after, max_qd,
                                     "max_qperp_dperp",
                                     test_note.str());
                                if (resolves_any) {
                                    ++n_recursive_coherent_candidate_vetoes_;
                                    return moliere::ScatteringDecision::VetoAndContinue;
                                }

                                ++n_recursive_coherent_candidate_accepts_;
                                probe.found = true;
                                probe.candidate = candidate;
                                return moliere::ScatteringDecision::StopBeforeApply;
                            };

                            const auto elastic_rng_state = moliere::elastic_generator_state();
                            moliere::seed_elastic_generator(static_cast<unsigned int>(coherent_seed));
                            if (remaining > 1.e-9 && probe.p[3] > 0.) {
                                moliere::propagate_segment_with_scattering_callback(
                                    probe.p, probe.pos, remaining,
                                    partons[active_ancestor].GetId(), probe.rng,
                                    kappa_, alpha_, tmethod_, mode_, ebe_hydro_,
                                    compat_moliere_legacy_hydro_, hydro_profile_,
                                    lres_moliere_particles, probe.had, probe.orient,
                                    callback);
                            }
                            moliere::set_elastic_generator_state(elastic_rng_state);
                            if (!probe.found) {
                                ++n_recursive_coherent_resample_exhausted_;
                                emit("recursive_coherent_resample_exhausted", active_ancestor,
                                     quenched[active_ancestor].GetMom(),
                                     quenched[active_ancestor].GetD1(),
                                     quenched[active_ancestor].GetD2(), probe.pos[3],
                                     probe.pos, probe.p, 0., "no_accepted_candidate",
                                     "modeE_reached_lres_interval_end");
                            }
                            return probe;
                        };

                        while (p[3] > 0. && pos[3] < total_end - 1.e-9) {
                            const int modee_iteration_id = modee_iteration++;
                            double window_end = total_end;
                            for (int group : active_groups) {
                                tighten_next_formation_boundary(group, pos[3], window_end);
                            }
                            std::vector<int> probe_frontier;
                            for (int group : active_groups) {
                                collect_probe_frontier(group, pos[3], probe_frontier);
                            }
                            std::sort(probe_frontier.begin(), probe_frontier.end());
                            probe_frontier.erase(std::unique(probe_frontier.begin(), probe_frontier.end()),
                                                 probe_frontier.end());
                            ++n_recursive_frontier_probe_batches_;
                            n_recursive_frontier_probe_objects_ +=
                                static_cast<long long>(probe_frontier.size());

                            const double probe_time = pos[3];
                            auto earlier_probe = [](const RecursiveProbe &probe,
                                                    const RecursiveProbe &chosen) {
                                if (!probe.found) return false;
                                if (!chosen.found) return true;
                                if (probe.candidate.pos[3] < chosen.candidate.pos[3] - 1.e-12) {
                                    return true;
                                }
                                if (std::abs(probe.candidate.pos[3] - chosen.candidate.pos[3]) <= 1.e-12 &&
                                    probe.probe < chosen.probe) {
                                    return true;
                                }
                                return false;
                            };
                            auto same_probe_choice = [](const RecursiveProbe &a,
                                                        const RecursiveProbe &b) {
                                if (a.found != b.found) return false;
                                if (!a.found) return true;
                                return a.probe == b.probe &&
                                       a.active_ancestor == b.active_ancestor &&
                                       std::abs(a.candidate.pos[3] - b.candidate.pos[3]) <= 1.e-12 &&
                                       std::abs(a.candidate.qperp - b.candidate.qperp) <= 1.e-12;
                            };
                            auto select_frontier_candidate = [&](const std::vector<int> &frontier_order,
                                                                bool record_candidates) {
                                RecursiveProbe selected;
                                for (int probe_id : frontier_order) {
                                    const auto start = project_from_active(probe_id, probe_time);
                                    if (!start.ok || start.p[3] <= 0.) continue;
                                    const int active_ancestor_for_probe = find_active_ancestor(probe_id);
                                    const int probe_seed = modeEBranchLocalSeed(
                                        nr_.GetIr(), event_id, modee_segment_id, modee_iteration_id,
                                        active_ancestor_for_probe, probe_id);
                                    RecursiveProbe probe = probe_object(
                                        probe_id, start, window_end, numrand(probe_seed), probe_seed);
                                    if (record_candidates && probe.found) {
                                        emit("recursive_probe_candidate", probe.probe, probe.active_ancestor,
                                             quenched[probe.probe].GetD1(), quenched[probe.probe].GetD2(),
                                             probe.candidate.pos[3], probe.candidate.pos,
                                             probe.candidate.p_after, probe.candidate.qperp,
                                             "q_perp",
                                             "modeE_branch_local_seed=" + std::to_string(probe_seed));
                                    }
                                    if (earlier_probe(probe, selected)) {
                                        selected = probe;
                                    }
                                }
                                return selected;
                            };

                            RecursiveProbe chosen = select_frontier_candidate(probe_frontier, true);
                            if ((dump_hybrid_evolution_history_ || do_event_display_) &&
                                probe_frontier.size() > 1) {
                                ++n_recursive_frontier_permutation_checks_;
                                std::vector<int> reversed_frontier = probe_frontier;
                                std::reverse(reversed_frontier.begin(), reversed_frontier.end());
                                const RecursiveProbe reversed_chosen =
                                    select_frontier_candidate(reversed_frontier, false);
                                if (!same_probe_choice(chosen, reversed_chosen)) {
                                    ++n_recursive_frontier_permutation_mismatches_;
                                    const RecursiveProbe &reported = chosen.found ? chosen : reversed_chosen;
                                    emit("recursive_probe_permutation_mismatch",
                                         reported.probe, reported.active_ancestor,
                                         reported.probe >= 0 ? quenched[reported.probe].GetD1() : -1,
                                         reported.probe >= 0 ? quenched[reported.probe].GetD2() : -1,
                                         reported.found ? reported.candidate.pos[3] : probe_time,
                                         reported.found ? reported.candidate.pos : pos,
                                         reported.found ? reported.candidate.p_after : p,
                                         reported.found ? reported.candidate.qperp : 0.,
                                         "modeE_frontier_permutation",
                                         "sorted_probe=" + std::to_string(chosen.probe) +
                                         ":reversed_probe=" + std::to_string(reversed_chosen.probe));
                                }
                            }

                            if (!chosen.found) {
                                propagate_active_groups(window_end);
                                recombine_active_groups();
                                if (window_end >= total_end - 1.e-9) break;
                                continue;
                            }

                            const double candidate_time = std::min(chosen.candidate.pos[3], total_end);
                            propagate_active_groups(candidate_time);

                            int final_apply = chosen.probe;
                            int current_object = chosen.probe;
                            int child_for_test = chosen.probe;
                            bool resolved_by_candidate = false;
                            int tested_dipoles = 0;
                            double last_qd = 0.;
                            double candidate_qd_for_average = 0.;
                            const int active_ancestor = chosen.active_ancestor;
                            ++n_unresolved_candidate_scatters_;
                            ++n_recursive_frontier_candidates_;
                            while (true) {
                                const int parent = quenched[child_for_test].GetMom();
                                if (parent < 0 || parent >= static_cast<int>(n) ||
                                    !descends_from(parent, active_ancestor)) {
                                    final_apply = current_object;
                                    break;
                                }
                                const int sibling = brother[child_for_test];
                                if (sibling < 0 || sibling >= static_cast<int>(n)) {
                                    final_apply = current_object;
                                    break;
                                }
                                const auto dperp_result = modeE_pair_dperp(
                                    child_for_test, sibling, chosen.candidate.pos[3]);
                                ++tested_dipoles;
                                const double dperp = dperp_result.dperp;
                                if (dperp_result.used_live_positions) ++n_recursive_live_dperp_tests_;
                                else ++n_recursive_vacuum_dperp_fallbacks_;
                                const double qd = chosen.candidate.qperp * dperp;
                                last_qd = qd;
                                candidate_qd_for_average = qd;
                                const bool resolves = passes_dynamic_moliere_resolution_test(
                                    chosen.candidate.qperp, dperp, moliere_unresolved_resolution_c_);
                                const std::string dperp_label = dperp_result.used_live_positions
                                                                    ? "qperp_dperp_live"
                                                                    : "qperp_dperp_vac_fallback";
                                emit("recursive_resolution_test", chosen.probe, parent,
                                     child_for_test, sibling, chosen.candidate.pos[3],
                                     chosen.candidate.pos, chosen.candidate.p_after, qd,
                                     dperp_label,
                                     resolves ? "modeE_resolves_tested_dipole"
                                              : "modeE_coherent_tested_dipole");
                                if (resolves) {
                                    final_apply = current_object;
                                    resolved_by_candidate = true;
                                    if (parent == active_ancestor) ++n_recursive_outer_resolutions_;
                                    else ++n_recursive_inner_resolutions_;
                                    break;
                                }
                                current_object = parent;
                                child_for_test = parent;
                            }
                            sum_qperp_dperp_unresolved_candidates_ += candidate_qd_for_average;

                            if (!resolved_by_candidate && tested_dipoles == 0 &&
                                chosen.probe == active_ancestor) {
                                // Before the first daughter formation, the active object
                                // itself is the frontier. Its proposal already comes from
                                // the correct coherent source and needs no veto/resampling.
                                auto applied_candidate = chosen.candidate;
                                applied_candidate.pos = qstate[active_ancestor].r;
                                applied_candidate.pos[3] = candidate_time;
                                apply_resolved_daughter_kick(
                                    applied_candidate, qstate[active_ancestor].p,
                                    lres_moliere_particles, qhad[active_ancestor],
                                    qorient[active_ancestor]);
                                qstate[active_ancestor].r = applied_candidate.pos;
                                ++n_unresolved_coherent_scatters_;
                                ++n_recursive_coherent_applications_;
                                emit("moliere_kick", active_ancestor,
                                     quenched[active_ancestor].GetMom(),
                                     quenched[active_ancestor].GetD1(),
                                     quenched[active_ancestor].GetD2(), candidate_time,
                                     applied_candidate.pos, qstate[active_ancestor].p,
                                     applied_candidate.qperp, "q_perp",
                                     "modeE_active_coherent_source_kick");
                                recombine_active_groups();
                                if (pos[3] <= candidate_time) {
                                    pos[3] = std::min(window_end, candidate_time + 1.e-6);
                                    for (int group : active_groups) {
                                        if (qstate[group].r[3] <= candidate_time) {
                                            qstate[group].r[3] = pos[3];
                                        }
                                    }
                                }
                                continue;
                            }

                            if (!resolved_by_candidate) {
                                ++n_recursive_failed_daughter_vetoes_;
                                emit("recursive_failed_daughter_veto", chosen.probe,
                                     active_ancestor, quenched[chosen.probe].GetD1(),
                                     quenched[chosen.probe].GetD2(),
                                     chosen.candidate.pos[3], chosen.candidate.pos,
                                     chosen.candidate.p_before, candidate_qd_for_average,
                                     "qperp_dperp",
                                     "modeE_discard_daughter_qperp_recoil_hole");

                                const bool colored_parent =
                                    isColored(partons[active_ancestor].GetId());
                                // Dani's failed-probe rule: once the daughter proposal is
                                // rejected, follow the coherent source until its first
                                // genuinely unresolving proposal or the LRES interval end.
                                // Formation times crossed along this search are included by
                                // collect_formed_pairs at each parent-candidate timestamp.
                                RecursiveProbe coherent = resample_coherent_source(
                                    active_ancestor, total_end, modee_iteration_id);
                                if (coherent.found) {
                                    const double coherent_time =
                                        std::min(coherent.candidate.pos[3], total_end);
                                    propagate_active_groups(coherent_time);
                                    auto applied_candidate = coherent.candidate;
                                    applied_candidate.pos = qstate[active_ancestor].r;
                                    applied_candidate.pos[3] = coherent_time;
                                    apply_resolved_daughter_kick(
                                        applied_candidate, qstate[active_ancestor].p,
                                        lres_moliere_particles, qhad[active_ancestor],
                                        qorient[active_ancestor]);
                                    qstate[active_ancestor].r = applied_candidate.pos;
                                    ++n_unresolved_coherent_scatters_;
                                    ++n_recursive_coherent_applications_;
                                    emit("moliere_kick", active_ancestor,
                                         quenched[active_ancestor].GetMom(),
                                         quenched[active_ancestor].GetD1(),
                                         quenched[active_ancestor].GetD2(), coherent_time,
                                         applied_candidate.pos, qstate[active_ancestor].p,
                                         applied_candidate.qperp, "q_perp",
                                         "modeE_resampled_coherent_parent_kick");
                                    recombine_active_groups();
                                    if (pos[3] <= coherent_time) {
                                        pos[3] = std::min(total_end, coherent_time + 1.e-6);
                                        for (int group : active_groups) {
                                            if (qstate[group].r[3] <= coherent_time) {
                                                qstate[group].r[3] = pos[3];
                                            }
                                        }
                                    }
                                    continue;
                                }

                                if (colored_parent) {
                                    // No acceptable coherent-source scattering occurred
                                    // before this finite-LRES interval ended.
                                    propagate_active_groups(total_end);
                                } else {
                                    const double next_time =
                                        std::min(window_end, candidate_time + 1.e-6);
                                    for (int group : active_groups) {
                                        if (qstate[group].r[3] <= candidate_time) {
                                            qstate[group].r[3] = next_time;
                                        }
                                    }
                                }
                                recombine_active_groups();
                                if (window_end >= total_end - 1.e-9 &&
                                    pos[3] >= window_end - 1.e-9) {
                                    break;
                                }
                                continue;
                            }

                            if (resolved_by_candidate) {
                                ++n_unresolved_resolving_scatters_;
                                ++n_unresolved_pairs_elastically_decohered_;
                            }

                            if (final_apply != active_ancestor) {
                                ensure_active(final_apply);
                            }
                            if (!is_active(final_apply)) {
                                final_apply = active_ancestor;
                            }

                            auto applied_candidate = chosen.candidate;
                            // A resolving daughter proposal acts at the first enclosing
                            // dipole that it resolves. Failed proposals have already
                            // taken the separate veto/coherent-resampling path above.
                            applied_candidate.pos = qstate[final_apply].r;
                            applied_candidate.pos[3] = chosen.candidate.pos[3];
                            apply_resolved_daughter_kick(applied_candidate, qstate[final_apply].p,
                                                         lres_moliere_particles, qhad[final_apply],
                                                         qorient[final_apply]);
                            qstate[final_apply].r = applied_candidate.pos;
                            emit("moliere_kick", final_apply, quenched[final_apply].GetMom(),
                                 quenched[final_apply].GetD1(), quenched[final_apply].GetD2(),
                                 applied_candidate.pos[3], applied_candidate.pos,
                                 qstate[final_apply].p, applied_candidate.qperp,
                                 "q_perp",
                                 "modeE_recursive_resolving_kick");
                            emit("dynamic_resolution", final_apply, active_ancestor,
                                 quenched[final_apply].GetD1(), quenched[final_apply].GetD2(),
                                 applied_candidate.pos[3], applied_candidate.pos,
                                 qstate[final_apply].p, last_qd, "qperp_dperp",
                                 "modeE_recursive_tree_decoherence");
                            recombine_active_groups();
                            if (pos[3] <= candidate_time) {
                                pos[3] = std::min(total_end, candidate_time + 1.e-6);
                                for (int group : active_groups) {
                                    if (qstate[group].r[3] <= candidate_time) {
                                        qstate[group].r[3] = pos[3];
                                    }
                                }
                            }
                        }

                        propagate_active_groups(total_end);
                        if (active_groups.size() == 1 && active_groups.front() == idx) {
                            int c1_boundary = -1;
                            int c2_boundary = -1;
                            if (child_pair(idx, c1_boundary, c2_boundary) &&
                                effective_mom[c1_boundary] == idx &&
                                effective_mom[c2_boundary] == idx) {
                                // If no elastic candidate resolved the coherent object,
                                // the normal finite-LRES boundary still opens the pair.
                                // Seed the daughters from the live kicked parent so the
                                // coherent deflection is inherited by all later descendants.
                                split_active_parent(idx, "modeE_lres_boundary_resolution");
                            }
                        }
                        recombine_active_groups();
                        for (int group : active_groups) {
                            if (group != idx) {
                                modeE_seeded[group] = true;
                            }
                        }
                    } else if (do_moliere_dynamic_daughter_unresolved_resolution_ && has_unresolved_pair) {
                        segment_type = "modeD_dynamic_daughter_unresolved";
                        ++n_unresolved_segments_dynamic_;
                        const double total_end = pos[3] + tof;
                        bool elastically_decohered = false;

                        // Mode D: during a geometrically unresolved LRES interval, use
                        // daughter-level Moliere candidates to test whether a hard kick
                        // resolves the dipole.  Until a candidate passes q_perp*d_perp
                        // > c_res, continuous energy loss and unresolved kicks act on
                        // the coherent parent.
                        while (!elastically_decohered && p[3] > 0. && pos[3] < total_end - 1.e-9) {
                            const double remaining = std::max(0., total_end - pos[3]);
                            const double vac_parent_e = partons[idx].vGetP()[3];
                            const double frac = (vac_parent_e != 0.) ? p[3] / vac_parent_e : 0.;
                            auto p1 = partons[d1].vGetP() * frac;
                            auto p2 = partons[d2].vGetP() * frac;
                            const double daughter_e_sum = p1[3] + p2[3];
                            const double share1 = (daughter_e_sum > 0.) ? p1[3] / daughter_e_sum : 0.5;
                            const double share2 = 1. - share1;
                            const auto parent_mismatch = p - (p1 + p2);
                            p1 += parent_mismatch * share1;
                            p2 += parent_mismatch * share2;

                            auto pos1 = pos;
                            auto pos2 = pos;
                            int had1 = qhad[d1];
                            int had2 = qhad[d2];
                            if (qhad[idx] == 1 || qhad[idx] == 2) {
                                had1 = 2;
                                had2 = 2;
                            }
                            auto orient1 = orientationFor(p1);
                            auto orient2 = orientationFor(p2);

                            struct DaughterProbe {
                                bool found = false;
                                int daughter = -1;
                                moliere::ScatteringCandidate candidate;
                                std::array<double,4> p = {0., 0., 0., 0.};
                                std::array<double,4> pos = {0., 0., 0., 0.};
                                int had = 0;
                                std::array<double,4> orient = {0., 0., 0., 1.};
                                numrand rng;
                            };

                            auto probe_daughter = [&](int daughter, std::array<double,4> p_start,
                                                      std::array<double,4> pos_start, int had_start,
                                                      std::array<double,4> orient_start,
                                                      numrand rng_start) {
                                DaughterProbe probe;
                                probe.daughter = daughter;
                                probe.p = p_start;
                                probe.pos = pos_start;
                                probe.had = had_start;
                                probe.orient = orient_start;
                                probe.rng = rng_start;
                                auto callback = [&](const moliere::ScatteringCandidate &candidate) {
                                    probe.found = true;
                                    probe.candidate = candidate;
                                    return moliere::ScatteringDecision::StopBeforeApply;
                                };
                                moliere::propagate_segment_with_scattering_callback(
                                    probe.p, probe.pos, remaining, partons[daughter].GetId(),
                                    probe.rng, kappa_, alpha_, tmethod_, mode_, ebe_hydro_,
                                    compat_moliere_legacy_hydro_, hydro_profile_,
                                    lres_moliere_particles, probe.had, probe.orient, callback);
                                return probe;
                            };

                            DaughterProbe probe1 = probe_daughter(d1, p1, pos1, had1, orient1, nr_);
                            DaughterProbe probe2 = probe_daughter(d2, p2, pos2, had2, orient2, probe1.rng);
                            DaughterProbe chosen;
                            nr_ = probe2.rng;
                            if (probe1.found && probe2.found) {
                                chosen = probe1.candidate.pos[3] <= probe2.candidate.pos[3] ? probe1 : probe2;
                            } else if (probe1.found) {
                                chosen = probe1;
                            } else if (probe2.found) {
                                chosen = probe2;
                            } else {
                                loss_rate(p, pos, remaining, partons[idx].GetId(), length, tlength,
                                          event_id, &event_display_record_id, idx, quenched[idx].GetMom(),
                                          d1, d2, true);
                                break;
                            }

                            const double candidate_time = chosen.candidate.pos[3];
                            const double coherent_tof =
                                std::max(0., std::min(candidate_time, total_end) - pos[3]);
                            // Before the first candidate scattering, the pair is
                            // still unresolved by LRES, so continuous HYBRID
                            // energy loss acts on the coherent parent.
                            if (coherent_tof > 1.e-9 && p[3] > 0.) {
                                loss_rate(p, pos, coherent_tof, partons[idx].GetId(), length, tlength,
                                          event_id, &event_display_record_id, idx, quenched[idx].GetMom(),
                                          d1, d2, true);
                            }
                            if (pos[3] < candidate_time) {
                                pos = chosen.candidate.pos;
                            }

                            const double dperp = compute_unresolved_pair_dperp(
                                partons[d1].vGetP(), partons[d2].vGetP(),
                                life[d1].creation, chosen.candidate.pos[3]);
                            const double qd = chosen.candidate.qperp * dperp;
                            ++n_unresolved_candidate_scatters_;
                            sum_qperp_dperp_unresolved_candidates_ += qd;
                            const bool resolves = passes_dynamic_moliere_resolution_test(
                                chosen.candidate.qperp, dperp, moliere_unresolved_resolution_c_);
                            emit("dynamic_resolution_test", chosen.daughter, idx, d1, d2,
                                 chosen.candidate.pos[3], chosen.candidate.pos,
                                 chosen.candidate.p_after, qd, "qperp_dperp",
                                 resolves ? "daughter_candidate_resolves_dipole"
                                          : "daughter_candidate_coherent_dipole");

                            if (!resolves) {
                                ++n_unresolved_coherent_scatters_;
                                const double previous_segment_time = pos[3];
                                // The candidate was generated from a daughter
                                // probe, but its wavelength is too long to
                                // resolve the dipole.  Convert it into one
                                // coherent parent kick and one recoil/hole
                                // source, then continue looking for later
                                // candidate scatterings.
                                apply_resolved_daughter_kick(chosen.candidate, p,
                                                             lres_moliere_particles,
                                                             qhad[idx], qorient[idx]);
                                pos = chosen.candidate.pos;
                                if (pos[3] <= previous_segment_time) {
                                    pos[3] = std::min(total_end, previous_segment_time + 1.e-6);
                                }
                                emit("moliere_kick", idx, quenched[idx].GetMom(), d1, d2,
                                     chosen.candidate.pos[3], chosen.candidate.pos,
                                     p, chosen.candidate.qperp,
                                     "q_perp", "dynamic_daughter_coherent_parent_kick");
                                continue;
                            }

                            ++n_unresolved_resolving_scatters_;
                            ++n_unresolved_pairs_elastically_decohered_;
                            elastically_decohered = true;

                            // A resolving candidate breaks elastic coherence.
                            // Rebuild daughters from the updated parent after
                            // pre-candidate coherent energy loss, then apply
                            // only the sampled kick delta to the struck daughter.
                            const double vac_parent_e_after_loss = partons[idx].vGetP()[3];
                            const double frac_after_loss =
                                (vac_parent_e_after_loss != 0.) ? p[3] / vac_parent_e_after_loss : 0.;
                            p1 = partons[d1].vGetP() * frac_after_loss;
                            p2 = partons[d2].vGetP() * frac_after_loss;
                            const double daughter_e_sum_after_loss = p1[3] + p2[3];
                            const double share1_after_loss =
                                (daughter_e_sum_after_loss > 0.) ? p1[3] / daughter_e_sum_after_loss : 0.5;
                            const double share2_after_loss = 1. - share1_after_loss;
                            const auto parent_mismatch_after_loss = p - (p1 + p2);
                            p1 += parent_mismatch_after_loss * share1_after_loss;
                            p2 += parent_mismatch_after_loss * share2_after_loss;
                            pos1 = pos;
                            pos2 = pos;

                            if (chosen.daughter == d1) {
                                apply_resolved_daughter_kick(chosen.candidate, p1,
                                                             lres_moliere_particles, had1, orient1);
                            } else {
                                apply_resolved_daughter_kick(chosen.candidate, p2,
                                                             lres_moliere_particles, had2, orient2);
                            }
                            emit("moliere_kick", chosen.daughter, idx, d1, d2,
                                 chosen.candidate.pos[3], chosen.candidate.pos,
                                 (chosen.daughter == d1) ? p1 : p2,
                                 chosen.candidate.qperp,
                                 "q_perp", "dynamic_daughter_resolving_kick");
                            emit("dynamic_resolution", chosen.daughter, idx, d1, d2,
                                 chosen.candidate.pos[3], chosen.candidate.pos,
                                 p1 + p2, qd, "qperp_dperp",
                                 "daughter_candidate_elastic_decoherence");

                            const double remaining_after = std::max(0., total_end - chosen.candidate.pos[3]);
                            if (remaining_after > 0. && p1[3] > 0.) {
                                auto cb = make_moliere_scattering_callback(
                                    d1, partons[d1].GetId(), idx, d1, d2, true,
                                    "modeD_post_decoherence_daughter_scattering");
                                auto step_cb = make_moliere_step_callback(
                                    d1, partons[d1].GetId(), idx, d1, d2, true,
                                    "modeD_post_decoherence_daughter_step");
                                moliere::propagate_segment_with_scattering_callback(
                                    p1, pos1, remaining_after, partons[d1].GetId(),
                                    nr_, kappa_, alpha_, tmethod_, mode_, ebe_hydro_,
                                    compat_moliere_legacy_hydro_, hydro_profile_,
                                    lres_moliere_particles, had1, orient1, cb, step_cb);
                            }
                            if (remaining_after > 0. && p2[3] > 0.) {
                                auto cb = make_moliere_scattering_callback(
                                    d2, partons[d2].GetId(), idx, d1, d2, true,
                                    "modeD_post_decoherence_daughter_scattering");
                                auto step_cb = make_moliere_step_callback(
                                    d2, partons[d2].GetId(), idx, d1, d2, true,
                                    "modeD_post_decoherence_daughter_step");
                                moliere::propagate_segment_with_scattering_callback(
                                    p2, pos2, remaining_after, partons[d2].GetId(),
                                    nr_, kappa_, alpha_, tmethod_, mode_, ebe_hydro_,
                                    compat_moliere_legacy_hydro_, hydro_profile_,
                                    lres_moliere_particles, had2, orient2, cb, step_cb);
                            }

                            qhad[d1] = had1;
                            qhad[d2] = had2;
                            qorient[d1] = orient1;
                            qorient[d2] = orient2;
                            // The surrounding finite-LRES code still expects
                            // one effective object until the original LRES
                            // segment boundary, so recombine the independently
                            // propagated daughters at the end of this segment.
                            p = p1 + p2;
                            if (p[3] > 0.) {
                                pos[0] = (p1[3] * pos1[0] + p2[3] * pos2[0]) / p[3];
                                pos[1] = (p1[3] * pos1[1] + p2[3] * pos2[1]) / p[3];
                                pos[2] = (p1[3] * pos1[2] + p2[3] * pos2[2]) / p[3];
                            } else {
                                pos[0] = 0.;
                                pos[1] = 0.;
                                pos[2] = 0.;
                            }
                            pos[3] = std::max(pos1[3], pos2[3]);
                            qhad[idx] = (had1 == 1 || had1 == 2 || had2 == 1 || had2 == 2) ? 2 : qhad[idx];
                            qorient[idx] = orientationFor(p);
                        }
                    } else if (do_moliere_dynamic_unresolved_resolution_ && has_unresolved_pair) {
                        segment_type = "modeC_dynamic_parent_unresolved";
                        ++n_unresolved_segments_dynamic_;
                        const double total_end = pos[3] + tof;
                        bool elastically_decohered = false;
                        moliere::ScatteringCandidate resolving_candidate;
                        double resolving_qd = 0.;

                        auto dynamic_callback = [&](const moliere::ScatteringCandidate &candidate) {
                            // Mode C samples a candidate on the coherent parent.
                            // The q_perp*d_perp test decides whether this parent
                            // scattering should remain coherent or should
                            // elastically decohere the unresolved daughter pair.
                            const double dperp = compute_unresolved_pair_dperp(
                                partons[d1].vGetP(), partons[d2].vGetP(),
                                life[d1].creation, candidate.pos[3]);
                            const double qd = candidate.qperp * dperp;
                            ++n_unresolved_candidate_scatters_;
                            sum_qperp_dperp_unresolved_candidates_ += qd;

                            const bool resolves = passes_dynamic_moliere_resolution_test(
                                candidate.qperp, dperp, moliere_unresolved_resolution_c_);
                            emit("dynamic_resolution_test", idx, quenched[idx].GetMom(), d1, d2,
                                 candidate.pos[3], candidate.pos, candidate.p_after,
                                 qd, "qperp_dperp",
                                 resolves ? "resolves_dipole" : "coherent_dipole");

                            if (!resolves) {
                                ++n_unresolved_coherent_scatters_;
                                return apply_coherent_unresolved_kick();
                            }

                            ++n_unresolved_resolving_scatters_;
                            elastically_decohered = true;
                            resolving_candidate = candidate;
                            resolving_qd = qd;
                            return moliere::ScatteringDecision::StopBeforeApply;
                        };

                        moliere::propagate_segment_with_scattering_callback(
                            p, pos, tof, partons[idx].GetId(), nr_, kappa_, alpha_,
                            tmethod_, mode_, ebe_hydro_, compat_moliere_legacy_hydro_,
                            hydro_profile_, lres_moliere_particles, qhad[idx],
                            qorient[idx], dynamic_callback,
                            make_moliere_step_callback(idx, partons[idx].GetId(), quenched[idx].GetMom(),
                                                       d1, d2, true, "modeC_dynamic_parent_step"));

                        if (elastically_decohered) {
                            ++n_unresolved_pairs_elastically_decohered_;

                            const double vac_parent_e = partons[idx].vGetP()[3];
                            const double frac = (vac_parent_e != 0.) ? p[3] / vac_parent_e : 0.;
                            auto p1 = partons[d1].vGetP() * frac;
                            auto p2 = partons[d2].vGetP() * frac;
                            const double daughter_e_sum = p1[3] + p2[3];
                            const double share1 = (daughter_e_sum > 0.) ? p1[3] / daughter_e_sum : 0.5;
                            const double share2 = 1. - share1;
                            const auto parent_mismatch = p - (p1 + p2);
                            p1 += parent_mismatch * share1;
                            p2 += parent_mismatch * share2;
                            auto pos1 = pos;
                            auto pos2 = pos;

                            int had1 = qhad[d1];
                            int had2 = qhad[d2];
                            if (qhad[idx] == 1 || qhad[idx] == 2) {
                                had1 = 2;
                                had2 = 2;
                            }
                            auto orient1 = orientationFor(p1);
                            auto orient2 = orientationFor(p2);

                            // Because the candidate was generated from the
                            // coherent parent, Mode C does not know which
                            // daughter was struck microscopically.  Assign the
                            // resolving kick to a daughter with probability
                            // proportional to the positive daughter energy.
                            const double delta_e =
                                resolving_candidate.p_after[3] - resolving_candidate.p_before[3];
                            double prob_d1 = 0.5;
                            const double positive_e_sum = std::max(0., p1[3]) + std::max(0., p2[3]);
                            if (positive_e_sum > 0.) {
                                prob_d1 = std::max(0., p1[3]) / positive_e_sum;
                            }
                            int struck = (nr_.rando() < prob_d1) ? d1 : d2;
                            if (struck == d1 && p1[3] + delta_e <= 0. && p2[3] + delta_e > 0.) struck = d2;
                            if (struck == d2 && p2[3] + delta_e <= 0. && p1[3] + delta_e > 0.) struck = d1;

                            if (struck == d1) {
                                apply_resolved_daughter_kick(resolving_candidate, p1,
                                                             lres_moliere_particles, had1, orient1);
                                emit("moliere_kick", d1, idx, d1, d2,
                                     resolving_candidate.pos[3], resolving_candidate.pos,
                                     p1, resolving_candidate.qperp,
                                     "q_perp", "dynamic_resolving_daughter_kick");
                            } else {
                                apply_resolved_daughter_kick(resolving_candidate, p2,
                                                             lres_moliere_particles, had2, orient2);
                                emit("moliere_kick", d2, idx, d1, d2,
                                     resolving_candidate.pos[3], resolving_candidate.pos,
                                     p2, resolving_candidate.qperp,
                                     "q_perp", "dynamic_resolving_daughter_kick");
                            }
                            emit("dynamic_resolution", struck, idx, d1, d2,
                                 resolving_candidate.pos[3], resolving_candidate.pos,
                                 p1 + p2, resolving_qd,
                                 "qperp_dperp", "elastic_decoherence");

                            const double remaining_after = std::max(0., total_end - pos[3]);
                            if (remaining_after > 0. && p1[3] > 0.) {
                                auto cb = make_moliere_scattering_callback(
                                    d1, partons[d1].GetId(), idx, d1, d2, true,
                                    "modeC_post_decoherence_daughter_scattering");
                                auto step_cb = make_moliere_step_callback(
                                    d1, partons[d1].GetId(), idx, d1, d2, true,
                                    "modeC_post_decoherence_daughter_step");
                                moliere::propagate_segment_with_scattering_callback(
                                    p1, pos1, remaining_after, partons[d1].GetId(),
                                    nr_, kappa_, alpha_, tmethod_, mode_, ebe_hydro_,
                                    compat_moliere_legacy_hydro_, hydro_profile_,
                                    lres_moliere_particles, had1, orient1, cb, step_cb);
                            }
                            if (remaining_after > 0. && p2[3] > 0.) {
                                auto cb = make_moliere_scattering_callback(
                                    d2, partons[d2].GetId(), idx, d1, d2, true,
                                    "modeC_post_decoherence_daughter_scattering");
                                auto step_cb = make_moliere_step_callback(
                                    d2, partons[d2].GetId(), idx, d1, d2, true,
                                    "modeC_post_decoherence_daughter_step");
                                moliere::propagate_segment_with_scattering_callback(
                                    p2, pos2, remaining_after, partons[d2].GetId(),
                                    nr_, kappa_, alpha_, tmethod_, mode_, ebe_hydro_,
                                    compat_moliere_legacy_hydro_, hydro_profile_,
                                    lres_moliere_particles, had2, orient2, cb, step_cb);
                            }

                            qhad[d1] = had1;
                            qhad[d2] = had2;
                            qorient[d1] = orient1;
                            qorient[d2] = orient2;

                            p = p1 + p2;
                            if (p[3] > 0.) {
                                pos[0] = (p1[3] * pos1[0] + p2[3] * pos2[0]) / p[3];
                                pos[1] = (p1[3] * pos1[1] + p2[3] * pos2[1]) / p[3];
                                pos[2] = (p1[3] * pos1[2] + p2[3] * pos2[2]) / p[3];
                            } else {
                                pos[0] = 0.;
                                pos[1] = 0.;
                                pos[2] = 0.;
                            }
                            pos[3] = std::max(pos1[3], pos2[3]);
                            qhad[idx] = (had1 == 1 || had1 == 2 || had2 == 1 || had2 == 2) ? 2 : qhad[idx];
                            qorient[idx] = orientationFor(p);
                        }
                    } else if (do_moliere_on_unresolved_partons_ && has_unresolved_pair) {
                        segment_type = "modeB_individual_unresolved_daughters";
                        // Mode B keeps the LRES timeline unchanged but lets
                        // Moliere act on the unresolved daughters separately.
                        // The parent energy sets the daughter energy scale;
                        // after propagation the daughters are recombined into
                        // the effective parent required by the LRES segment.
                        const double vac_parent_e = partons[idx].vGetP()[3];
                        const double frac = (vac_parent_e != 0.) ? p[3] / vac_parent_e : 0.;

                        auto p1 = partons[d1].vGetP() * frac;
                        auto p2 = partons[d2].vGetP() * frac;
                        auto pos1 = pos;
                        auto pos2 = pos;
                        const auto p1_before = p1;
                        const auto p2_before = p2;
                        const auto pos1_before = pos1;
                        const auto pos2_before = pos2;

                        int had1 = qhad[d1];
                        int had2 = qhad[d2];
                        auto orient1 = qorient[d1];
                        auto orient2 = qorient[d2];

                        auto cb1 = make_moliere_scattering_callback(
                            d1, partons[d1].GetId(), idx, d1, d2, true,
                            "modeB_unresolved_daughter_scattering");
                        auto step_cb1 = make_moliere_step_callback(
                            d1, partons[d1].GetId(), idx, d1, d2, true,
                            "modeB_unresolved_daughter_step");
                        moliere::propagate_segment_with_scattering_callback(
                            p1, pos1, tof, partons[d1].GetId(), nr_, kappa_, alpha_,
                            tmethod_, mode_, ebe_hydro_, compat_moliere_legacy_hydro_,
                            hydro_profile_, lres_moliere_particles, had1, orient1, cb1, step_cb1);
                        auto cb2 = make_moliere_scattering_callback(
                            d2, partons[d2].GetId(), idx, d1, d2, true,
                            "modeB_unresolved_daughter_scattering");
                        auto step_cb2 = make_moliere_step_callback(
                            d2, partons[d2].GetId(), idx, d1, d2, true,
                            "modeB_unresolved_daughter_step");
                        moliere::propagate_segment_with_scattering_callback(
                            p2, pos2, tof, partons[d2].GetId(), nr_, kappa_, alpha_,
                            tmethod_, mode_, ebe_hydro_, compat_moliere_legacy_hydro_,
                            hydro_profile_, lres_moliere_particles, had2, orient2, cb2, step_cb2);

                        qhad[d1] = had1;
                        qhad[d2] = had2;
                        qorient[d1] = orient1;
                        qorient[d2] = orient2;

                        const double qperp1 = std::sqrt((p1[0] - p1_before[0]) * (p1[0] - p1_before[0]) +
                                                        (p1[1] - p1_before[1]) * (p1[1] - p1_before[1]));
                        const double qperp2 = std::sqrt((p2[0] - p2_before[0]) * (p2[0] - p2_before[0]) +
                                                        (p2[1] - p2_before[1]) * (p2[1] - p2_before[1]));
                        if (qperp1 > 1.e-12) {
                            emit("moliere_kick", d1, idx, d1, d2,
                                 pos1_before[3], pos1_before, p1, qperp1,
                                 "q_perp", "unresolved_daughter_kick");
                        }
                        if (qperp2 > 1.e-12) {
                            emit("moliere_kick", d2, idx, d1, d2,
                                 pos2_before[3], pos2_before, p2, qperp2,
                                 "q_perp", "unresolved_daughter_kick");
                        }

                        p = p1 + p2;
                        if (p[3] > 0.) {
                            pos[0] = (p1[3] * pos1[0] + p2[3] * pos2[0]) / p[3];
                            pos[1] = (p1[3] * pos1[1] + p2[3] * pos2[1]) / p[3];
                            pos[2] = (p1[3] * pos1[2] + p2[3] * pos2[2]) / p[3];
                        } else {
                            pos[0] = 0.;
                            pos[1] = 0.;
                            pos[2] = 0.;
                        }
                        pos[3] = std::max(pos1[3], pos2[3]);
                        qhad[idx] = (had1 == 1 || had1 == 2 || had2 == 1 || had2 == 2) ? 2 : qhad[idx];
                        qorient[idx] = orientationFor(p);
                    } else {
                        segment_type = has_unresolved_pair ? "modeA_coherent_unresolved_parent"
                                                           : "resolved_moliere_segment";
                        // Mode A is the backward-compatible default.  If the
                        // segment contains an unresolved pair, Moliere sees only
                        // the coherent effective parent.  If no unresolved pair
                        // is present, this is ordinary resolved-parton Moliere
                        // propagation.
                        auto cb = make_moliere_scattering_callback(
                            idx, partons[idx].GetId(), quenched[idx].GetMom(), d1, d2,
                            has_unresolved_pair, has_unresolved_pair
                                                     ? "modeA_coherent_parent_scattering"
                                                     : "resolved_parton_scattering");
                        auto step_cb = make_moliere_step_callback(
                            idx, partons[idx].GetId(), quenched[idx].GetMom(), d1, d2,
                            has_unresolved_pair, has_unresolved_pair
                                                     ? "modeA_coherent_parent_step"
                                                     : "resolved_parton_step");
                        moliere::propagate_segment_with_scattering_callback(
                            p, pos, tof, partons[idx].GetId(), nr_, kappa_, alpha_,
                            tmethod_, mode_, ebe_hydro_, compat_moliere_legacy_hydro_,
                            hydro_profile_, lres_moliere_particles, qhad[idx], qorient[idx],
                            cb, step_cb);
                    }
                    const double qperp = std::sqrt((p[0] - p_before[0]) * (p[0] - p_before[0]) +
                                                   (p[1] - p_before[1]) * (p[1] - p_before[1]));
                    if (qperp > 1.e-12) {
                        const std::string kick_note =
                            (do_moliere_recursive_unresolved_resolution_ && has_unresolved_pair)
                                ? "recursive_unresolved_net_kick"
                                : ((do_moliere_dynamic_daughter_unresolved_resolution_ && has_unresolved_pair)
                                       ? "dynamic_daughter_unresolved_net_kick"
                                       : ((do_moliere_dynamic_unresolved_resolution_ && has_unresolved_pair)
                                              ? "dynamic_unresolved_net_kick"
                                              : (has_unresolved_pair ? "coherent_unresolved_kick" : "resolved_parton_kick")));
                        emit("moliere_kick", idx, quenched[idx].GetMom(), d1, d2,
                             pos_before[3], pos_before, p, qperp, "q_perp", kick_note);
                    }
                } else {
                    segment_type = "lres_energy_loss_segment";
                    loss_rate(p, pos, tof, partons[idx].GetId(), length, tlength,
                              event_id, &event_display_record_id, idx, quenched[idx].GetMom(),
                              d1_segment, d2_segment, segment_has_unresolved_pair);
                }
            } else if (p[3] != 0.) {
                segment_type = "lres_free_stream_segment";
                pos += p / p[3] * tof;
                fill_event_display_record(event_id, event_display_record_id++,
                                          idx, partons[idx].GetId(), quenched[idx].GetMom(),
                                          d1_segment, d2_segment,
                                          -1, segment_has_unresolved_pair, false,
                                          pos_segment_start, p_segment_start, pos, p,
                                          0., 0., 0., 0.,
                                          "free_stream_step", "lres_free_stream");
            }

            fill_event_display_segment(event_id, event_display_segment_id++,
                                       idx, partons[idx], quenched[idx].GetMom(),
                                       d1_segment, d2_segment,
                                       segment_has_unresolved_pair, qhad[idx],
                                       pos_segment_start, p_segment_start,
                                       pos, p, length, tlength, segment_type);

            qstate[idx].p = p;
            qstate[idx].r = pos;
            done[idx] = true;
            quenched[idx].AddLength(length, tlength);

            if (p[3] <= 0.) {
                for (auto zero_it = it; zero_it != family.rend(); ++zero_it) {
                    qstate[*zero_it].p = {0., 0., 0., 0.};
                    done[*zero_it] = true;
                }
                break;
            }

            const double frac = (partons[idx].vGetP()[3] != 0.) ? qstate[idx].p[3] / partons[idx].vGetP()[3] : 0.;
            auto next_it = it;
            ++next_it;
            if (next_it != family.rend()) {
                const int daughter = *next_it;
                if (!modeE_seeded[daughter]) {
                    qstate[daughter].p = partons[daughter].vGetP() * frac;
                    if (do_elastic_) {
                        qorient[daughter] = orientationFor(qstate[daughter].p);
                        if (qhad[idx] == 1 || qhad[idx] == 2) qhad[daughter] = 2;
                    }
                }
            }
        }
    }

    if (do_elastic_) {
        std::vector<Quench> &recoiled_out = recoiled != nullptr ? *recoiled : local_recoiled;
        moliere::process_recoilers(lres_moliere_particles, nr_, kappa_, alpha_, tmethod_, mode_,
                                   ebe_hydro_, compat_moliere_legacy_hydro_, hydro_profile_, recoiled_out);
        for (const auto &rp : recoiled_out) {
            const std::string label = (rp.GetOrig() == "recoiler" || rp.GetOrig() == "hole") ? "response_parton" : "other";
            emit("medium_response", -1, -1, -1, -1, rp.GetRi()[3], rp.GetRi(), rp.vGetP(), 0.0, label, rp.GetOrig());
        }
    }

    for (size_t i = 0; i < n; ++i) {
        if (!done[i]) continue;
        quenched[i].vSetP(qstate[i].p);
        quenched[i].vSetRi(life[i].ri);
        quenched[i].vSetRf(qstate[i].r);
        quenched[i].vSetInhP(qstate[i].p);
        if (do_elastic_) {
            quenched[i].setOrient(qorient[i]);
            quenched[i].setHadScattering(qhad[i]);
        }
        quenched[i].SetIsDone(true);
    }

    append_history_records(history_records);
}

void EnergyLoss::append_history_records(const std::vector<std::string> &records) {
    if (!dump_hybrid_evolution_history_ || hybrid_evolution_history_file_.empty() || records.empty()) return;
    std::ofstream out(hybrid_evolution_history_file_, std::ios::app);
    if (!out.is_open()) return;
    if (out.tellp() == 0) {
        out << "event_id\tmode\trecord_type\tparton_id\tparent_id\td1\td2\t"
            << "time\tx\ty\tz\tpx\tpy\tpz\tE\tqperp\tlabel\tnote\n";
    }
    for (const auto &r : records) out << r << '\n';
}

void EnergyLoss::loss_rate(std::array<double,4> &p, std::array<double,4> &pos, double tof, int id,
                           double &length, double &tlength,
                           int event_id, int *record_id,
                           int parton_index, int parent_index,
                           int d1, int d2, bool is_unresolved) {
    double Tc;
    if (tmethod_ == 0) Tc = 0.170;
    else Tc = 0.145;
    constexpr double charm_mass = 1.25;
    constexpr double b_mass = 4.2;

    double tot = pos[3] + tof;    // Final time

    double tau0h = 0.6;             //Ave hydro
    if (ebe_hydro_ == 1) tau0h = 0.4;   //ebe hydro

    double ei = p[3];      // Initial energy

    double f_dist = 0.;    // Traversed distance in Fluid Frame

    double virt_f_dist = 0.;

    double CF;
    if (id == 21) {
        if (mode_ == 0) CF = pow(9. / 4., 1. / 3.);  // If gluon, color charge dependence is ratio of casimirs to power 1/3
        else CF = 9. / 4.;
    } else CF = 1.;

    int marker = 0;    // If one, exit loop
    double step = 0.1;  // Time step in LAB frame

    auto w = p / p[3];  // 4-velocity

    do {
        const auto pos_step_start = pos;
        const auto p_step_start = p;
        double temp_for_record = 0.;
        double step_length_before = length;
        double step_tlength_before = tlength;
        bool in_medium_for_record = false;
#ifdef DO_SOURCE
        // Keep 4momentum before applying quenching this step
        auto p_prev = p;
#endif
        auto p_pre_floor = p;

        if (pos[3] == tot) marker = 1;
        if (pos[3] > tot) std::cout << " Warning: Went beyond tot= " << tot << " t= " << pos[3] << std::endl;

        // Proper time
        double tau = sqrt(pos[3] * pos[3] - pos[2] * pos[2]);
        if (tau != tau) {
            std::cout << " TAU Not a number z= " << pos[2] << " t= " << pos[3] << " wz= " << w[2] << " en = " << p[3] << " pz= " << p[2] << "\n";
            std::cout << " Id= " << id << std::endl;
            exit(1);
        }

        // Rapidity
        double eta = 1. / 2. * log((pos[3] + pos[2]) / (pos[3] - pos[2]));
        if (eta != eta && tau > 0.) {
            std::cout << " Eta is NaN= " << eta << " t= " << pos[3] << " z= " << pos[2] << std::endl;
            exit(1);
        }

        int will_hot = 0;  // Advance variable (to reach hot zones)
        double vx = 0.;
        double vy = 0.;
        if (tau >= tau0h || hydro_profile_.hasPreHydroAt(tau)) {  // Hydro or optional pre-hydro profile
            double temp = 0.;
            hydro_profile_.getValues(tau, pos[0], pos[1], temp, vx, vy);
            temp_for_record = temp;
            in_medium_for_record = temp >= Tc;

            double vz = pos[2] / pos[3];
            double frap = atanh(vz);
            vx /= cosh(frap);
            vy /= cosh(frap);

            std::array<double,4> v = {vx, vy, vz, 1.};

            double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
            double w2 = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
            double vscalw = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
            if (v2 >= 1.) v2 = 0.999999999;
            double lore = 1. / sqrt(1. - v2);

            double f_lore = w2 + lore * lore * (v2 - 2. * vscalw + vscalw * vscalw);
            if (f_lore < 0.) {
                f_lore = 0.;
            }
            double f_step = step * sqrt(f_lore);
            f_dist += f_step;

            // temp is already available from getValues
            
            if (temp > Tc) {
                // In-medium distance tracked here (for potential future use)
            }

            // Safe way to exit the plasma: check whether temperature will be above Tc in the next 1000 steps
            if (temp < Tc) {
                // Check whether temperature increases in its way
                for (unsigned int j = 1; j < 1000; j++) {
                    double step_j = step * double(j);
                    double tpos0 = pos[0] + w[0]*step_j;
                    double tpos1 = pos[1] + w[1]*step_j;
                    double tpos2 = pos[2] + w[2]*step_j;
                    double tpos3 = pos[3] + w[3]*step_j;
                    if (tpos3 > tot) break;
                    tau = sqrt(tpos3*tpos3 - tpos2*tpos2);
                    eta = 1. / 2. * log((tpos3 + tpos2) / (tpos3 - tpos2));
                    double ctemp = call_gT(tau, tpos0, tpos1, 0);
                    if (ctemp > Tc) { //It will get to hot
                        will_hot = int(j);
                        break;
                    }
                }
                if (will_hot == 0) { // It will not get to hot
                    pos += w * (tot - pos[3]);
                    marker = 1;
                }
            }

            // Now broad&quench
            if (p[3] > 0. && temp >= Tc && f_step != 0.) {
                if (mode_ == 0) {
                    length += f_step;
                    tlength += temp / 0.2 * f_step;
                } else {
                    length += f_step;
                    tlength += (temp/0.2)*(temp/0.2) * f_step;
                }

                // Broadening
                if (kappa_ != 0.) {
                    trans_kick(w, w2, v, p, temp, vscalw, lore, step, kappa_);
                }

                p_pre_floor = p;
                bool doquench = true;
                if (std::abs(id) == 4 && p[3] <= charm_mass) doquench = false;
                if (std::abs(id) == 5 && p[3] <= b_mass) doquench = false;

                // Strong coupling
                if (alpha_ != 0. && mode_ == 0 && doquench) {
                    double Efs = ei * lore * (1. - vscalw);
                    double tstop = 0.2 * pow(Efs, 1. / 3.) / (2. * pow(temp, 4. / 3.) * alpha_) / CF;
                    double beta = tstop / f_dist;
                    if (beta > 1.) {
                        double intpiece = Efs * step * 4. / (3.141592) * (1. / (beta * tstop * sqrt(beta * beta - 1.)));
                        double quench = (p[3] - intpiece) / p[3];
                        p *= quench;
                    } else {
                        p[3] = 0.;
                    }
                }

                // Radiative
                if (alpha_ != 0. && mode_ == 1) {
                    double intpiece = CF * (step / 0.2) * alpha_ * temp * temp * temp * (f_dist / 0.2);
                    double quench = (p[3] - intpiece) / p[3];
                    p *= quench;
                }

                // Collisional
                if (alpha_ != 0. && mode_ == 2) {
                    double intpiece = CF * (step / 0.2) * alpha_ * temp * temp;
                    double quench = (p[3] - intpiece) / p[3];
                    p *= quench;
                }
                
            }
        }

        if (std::abs(id) == 4 && p[3] < charm_mass) {
            p[3] = charm_mass;
            double pmod = std::sqrt(p_pre_floor[0] * p_pre_floor[0] + p_pre_floor[1] * p_pre_floor[1] +
                                    p_pre_floor[2] * p_pre_floor[2]);
            if (pmod == 0.) pmod = 1.;
            p[0] = p_pre_floor[0] / pmod * p[3];
            p[1] = p_pre_floor[1] / pmod * p[3];
            p[2] = p_pre_floor[2] / pmod * p[3];
        }
        if (std::abs(id) == 5 && p[3] < b_mass) {
            p[3] = b_mass;
            double pmod = std::sqrt(p_pre_floor[0] * p_pre_floor[0] + p_pre_floor[1] * p_pre_floor[1] +
                                    p_pre_floor[2] * p_pre_floor[2]);
            if (pmod == 0.) pmod = 1.;
            p[0] = p_pre_floor[0] / pmod * p[3];
            p[1] = p_pre_floor[1] / pmod * p[3];
            p[2] = p_pre_floor[2] / pmod * p[3];
        }

        //This is to check travelled distance, not used
        if (tof < 1000000) {
            double vz = pos[2] / std::max(pos[3], 0.000001);
            double v2 = vz * vz;
            double w2 = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
            double vscalw = vz * w[2];
            if (v2 >= 1.) v2 = 0.999999999;
            double lore = 1. / sqrt(1. - v2);
            double f_lore = w2 + lore * lore * (v2 - 2. * vscalw + vscalw * vscalw);
            if (f_lore < 0.) {
                std::cout << " craazy f_lore= " << f_lore << std::endl;
                f_lore = 1.;
            }
            virt_f_dist += step * sqrt(f_lore);
        }

        // If parton gets totally quenched, exit
        if (p[3] <= 0.) {
            marker = 1;
            for (unsigned int i = 0; i < 4; i++) p[i] = 0.;
        } else {
            // Manually protect very soft particles from getting kicks that yield velocities greater than 1
            for (unsigned int i = 0; i < 3; i++) {
                if (p[i] > p[3]) {
                    std::cout << " Got crazy kick in i= " << i << "p[i]= " << p[i] << " and p[3]= " << p[3] << std::endl;
                    p[i] = 0.99999 * p[3];
                }
            }
            // Update kinematical quantities, with the possibility of advancing to hot regions
            w = p / p[3];
            double tstep = std::max(double(will_hot), 1.) * step;
            if (pos[3] + tstep > tot) {
                tstep = tot - pos[3];
            }
            if (marker != 1) pos += w * tstep;
        }

        if (do_event_display_ && record_id != nullptr && event_id >= 0) {
            fill_event_display_record(event_id, (*record_id)++, parton_index, id,
                                      parent_index, d1, d2, -1,
                                      is_unresolved, in_medium_for_record,
                                      pos_step_start, p_step_start, pos, p,
                                      std::sqrt((p[0] - p_step_start[0]) * (p[0] - p_step_start[0]) +
                                                (p[1] - p_step_start[1]) * (p[1] - p_step_start[1])),
                                      temp_for_record,
                                      length - step_length_before,
                                      tlength - step_tlength_before,
                                      "energy_loss_step",
                                      "hybrid_integration_step");
        }

#ifdef DO_SOURCE
        // Fill source file, for new wake purposes
        if (p[3] != p_prev[3]) {
            // Get tau ev, x_f, y_f and vx_f and vy_f for source file
            double tau_ev, x_f, y_f, vx_f, vy_f;
            get_source_evol(tau_ev, x_f, y_f, vx_f, vy_f, tau, pos[0], pos[1], Tc);

            std::ofstream source_file("SOURCE.dat", std::ios_base::app);
            source_file << tau << " " << pos[0] << " " << pos[1] << " " << eta << " " << tau_ev << " "
                        << -p[3] + p_prev[3] << " " << -p[0] + p_prev[0] << " " << -p[1] + p_prev[1] << " " << -p[2] + p_prev[2] << " "
                        << vx << " " << vy << " " << vx_f << " " << vy_f << " " << x_f << " " << y_f << std::endl;
        }
#endif

    } while (marker == 0);

    //Just to check the distance travelled, not used
    if (virt_f_dist == 0.) virt_f_dist = -1;
}

double EnergyLoss::call_gT(double tau, double x, double y, int comp) const {
    if (tau < 0.) return 0.0;

    switch (comp) {
        case 0:
            return hydro_profile_.temperature(tau, x, y);
        case 1:
            return hydro_profile_.velocityX(tau, x, y);
        case 2:
            return hydro_profile_.velocityY(tau, x, y);
        default:
            std::cerr << "Wrong comp= " << comp << " in call_gT" << std::endl;
            return 0.0;
    }
}

double EnergyLoss::resolution_time(double parent_e, double parent_px, double parent_py, double parent_pz,
                                   double x, double y, double z,
                                   double dvx, double dvy, double dvz,
                                   double t0) const {
    if (lres_rpower_ < 0.00001) return kLresInfiniteTime;
    if (parent_e == 0.) return 0.;

    const double Tc = (tmethod_ == 0) ? 0.170 : 0.145;
    const double scale = 1. / 3.14 / lres_rpower_;
    const double wx = parent_px / parent_e;
    const double wy = parent_py / parent_e;
    const double wz = parent_pz / parent_e;

    double ti = t0;
    double xprime = 0.;
    double yprime = 0.;
    double zprime = 0.;
    double step_res = 0.1;
    bool refine_step = false;

    for (unsigned int i = 0; i <= 100000000; ++i) {
        const double proper_sq = ti * ti - z * z;
        const double tau = proper_sq > 0. ? std::sqrt(proper_sq) : 0.;

        if (tau >= 0.6 || hydro_profile_.hasPreHydroAt(tau)) {
            const double temp = call_gT(tau, x, y, 0);
            const double sep = std::sqrt(xprime * xprime + yprime * yprime + zprime * zprime);
            if (sep >= scale / (temp / 0.2) || temp <= Tc) {
                if (step_res == 0.1) {
                    if (i == 0) return ti - t0;
                    refine_step = true;
                } else {
                    return ti - t0;
                }
            }
        }

        if (!refine_step) {
            x += wx * step_res;
            y += wy * step_res;
            z += wz * step_res;
            xprime += dvx * step_res;
            yprime += dvy * step_res;
            zprime += dvz * step_res;
            ti += step_res;
        } else {
            x -= wx * step_res;
            y -= wy * step_res;
            z -= wz * step_res;
            xprime -= dvx * step_res;
            yprime -= dvy * step_res;
            zprime -= dvz * step_res;
            ti -= step_res;
            step_res = 0.01;
            refine_step = false;
        }
    }

    return kLresInfiniteTime;
}

void EnergyLoss::get_source_evol(double &tau_ev, double& x_f, double& y_f, double& vx_f, double& vy_f, double tau_ini, double x_ini, double y_ini, double Tc) {
    tau_ev = 0.;
    double dtau = 0.1;
    x_f = x_ini;
    y_f = y_ini;
    double tau_now = tau_ini;
    double vx_f_local = 0., vy_f_local = 0.;

    while (true) {
        vx_f_local = call_gT(tau_now, x_f, y_f, 1);
        vy_f_local = call_gT(tau_now, x_f, y_f, 2);

        double localT;
        localT = call_gT(tau_now, x_f, y_f, 0);
        if (localT < Tc) break;

        x_f += vx_f_local * dtau;
        y_f += vy_f_local * dtau;
        tau_ev += dtau;
        tau_now += dtau;
    }

    vx_f = vx_f_local;
    vy_f = vy_f_local;
}

void EnergyLoss::trans_kick(const std::array<double,4> &w, double w2, const std::array<double,4> &v, std::array<double,4> &p, double temp, double vscalw, double lore, double step, double kappa) {
    if (vscalw == 1.) return;

    auto e1 = vec_prod(w, v);
    double Ne1 = normalise(e1);
    if (Ne1 == 0.) {
        double b = 0.5;
        double c = 0.2;
        double a = (-b * w[1] - c * w[2]) / w[0];
        e1 = {a, b, c, 0.};
        Ne1 = normalise(e1);
    }

    double Nw = sqrt(w2);
    auto l = vec_prod(w, e1) / Nw;

    double uscalW = lore * (1. - vscalw);
    double uscall = lore * (-v[0] * l[0] - v[1] * l[1] - v[2] * l[2]);
    double W2 = 1. - w2;

    auto Wp = w - v * (lore * W2 / uscalW);

    double Nalpha = -uscall * uscalW / (uscalW*uscalW - W2);
    if (Nalpha != Nalpha || std::isinf(Nalpha)) return;
    double NN = 1. + W2 * uscall*uscall / (-(uscalW*uscalW) + W2);
    // In some rare situations, this norm squared can be negative. Only do kick otherwise
    if (sqrt(NN) != sqrt(NN)) std::cout << " negative NN " << std::endl;
    else {
        auto e2 = (l + Wp * Nalpha) / sqrt(NN);

        double Ef = p[3] * lore * (1. - vscalw);
        double lore_1mv = lore * (1. - vscalw);
        double wf2 = 1. - W2 / (lore_1mv * lore_1mv);
        double DelQ2 = kappa * temp*temp*temp * lore * (1. - vscalw) * step * 5.;

        double qfac = 0.;
        double qbeta = 0.;
        const bool valid_rest_frame_energy =
            std::isfinite(Ef) && Ef > 0. && std::isfinite(wf2) && wf2 > 0.;
        // A zero local-rest-frame energy is a zero kick, not a valid 0/0 limit.
        if (valid_rest_frame_energy) {
            const double cutoff = Ef * sqrt(wf2);
            qfac = sqrt(-1. * log(nr_.rando())) * sqrt(DelQ2); // Box-Muller method
            if (!std::isfinite(qfac)) {
                qfac = 0.;
            } else if (qfac > cutoff) {
                qfac = std::nextafter(cutoff, 0.);
            }

            const double beta_argument = 1. - qfac * qfac / (Ef * Ef * wf2);
            if (std::isfinite(beta_argument)) {
                qbeta = sqrt(std::clamp(beta_argument, 0., 1.)) - 1.;
            } else {
                qfac = 0.;
            }
        }

        double qphi = 2. * 3.141592654 * nr_.rando();

        auto e = e1 * cos(qphi) + e2 * sin(qphi);

        auto Wt = (w - v * (uscalW * lore)) / lore / (1. - vscalw);

        // Update 4momentum
        if (lore == 1.) {
            e2 = vec_prod(e1, w);
            e = e1 * cos(qphi) + e2 * sin(qphi);
        }
        p += Wt * qbeta * Ef + e * qfac;
        if (p[3] != p[3]) {
            std::cout << " p in bro= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
            std::cout << " qbeta= " << qbeta << " qfac= " << qfac << std::endl;
            std::cout << " Wt= " << Wt[0] << " " << Wt[1] << " " << Wt[2] << " " << Wt[3] << std::endl;
            std::cout << " e1= " << e1[0] << " " << e1[1] << " " << e1[2] << " " << e1[3] << std::endl;
            std::cout << " e2= " << e2[0] << " " << e2[1] << " " << e2[2] << " " << e2[3] << std::endl;
            std::cout << " Nalpha= " << Nalpha << std::endl;
            std::cout << " l= " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
            std::cout << " Wp= " << Wp[0] << " " << Wp[1] << " " << Wp[2] << " " << Wp[3] << std::endl;
            std::cout << " uscalW= " << uscalW << std::endl;
            std::cout << " Nw= " << Nw << " uscall= " << uscall << std::endl;
            std::cout << " W2= " << W2 << std::endl;
            exit(1);
        }
    }
}

void EnergyLoss::quenched_sons(const std::array<double,4> &p, const std::array<double,4> &qp, std::array<double,4> &d1, std::array<double,4> &d2) {
    // Mutable local copies for normalisation
    std::array<double,4> p_n = p;
    std::array<double,4> qp_n = qp;
    // Normalise 3-momentum
    double qmod = normalise(qp_n);
    double modmom = normalise(p_n);
    double modthis = normalise(d1);
    double modoson = normalise(d2);
    // Define transverse axis for rotation and normalise
    auto axis = vec_prod(p_n, qp_n);
    normalise(axis);
    // Find angle in plane
    double angle = 0.;
    if (p_n[0]*qp_n[0] + p_n[1]*qp_n[1] + p_n[2]*qp_n[2] >= 1.) angle = 0.;
    else angle = acos(p_n[0]*qp_n[0] + p_n[1]*qp_n[1] + p_n[2]*qp_n[2]);
    // Perform Rodrigues rotation
    double thisscal = axis[0]*d1[0] + axis[1]*d1[1] + axis[2]*d1[2];
    auto use = d1 * cos(angle) + vec_prod(axis, d1) * sin(angle) + axis * (thisscal * (1. - cos(angle)));
    double oscal = axis[0]*d2[0] + axis[1]*d2[1] + axis[2]*d2[2];
    auto ouse = d2 * cos(angle) + vec_prod(axis, d2) * sin(angle) + axis * (oscal * (1. - cos(angle)));
    // Update momenta
    double lamp = qmod / modmom;
    double lambda = qp[3] / p[3];
    for (unsigned int i = 0; i < 3; i++) {
        d1[i] = use[i] * modthis * lamp;
        d2[i] = ouse[i] * modoson * lamp;
    }
    d1[3] *= lambda;
    d2[3] *= lambda;
}

double EnergyLoss::normalise(std::array<double,4> &p) {
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    if (norm == 0.) return norm;
    p[0] /= norm;
    p[1] /= norm;
    p[2] /= norm;
    return norm;
}

std::array<double,4> EnergyLoss::vec_prod(const std::array<double,4> &a, const std::array<double,4> &b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0], 0.};
}
