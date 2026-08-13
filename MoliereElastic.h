#pragma once

#include <array>
#include <functional>
#include <random>
#include <vector>

#include "HydroProfile.h"
#include "Parton.h"
#include "Quench.h"
#include "Random.h"
#include "TransportState.h"

namespace moliere {

struct ScatteringCandidate {
    std::array<double,4> pos = {0., 0., 0., 0.};
    std::array<double,4> p_before = {0., 0., 0., 0.};
    std::array<double,4> p_after = {0., 0., 0., 0.};
    std::array<double,4> recoiler_p = {0., 0., 0., 0.};
    std::array<double,4> hole_p = {0., 0., 0., 0.};
    double qperp = 0.;
    int recoiler_id = 0;
    int hole_id = 0;
};

struct PropagationStep {
    std::array<double,4> pos_before = {0., 0., 0., 0.};
    std::array<double,4> pos_after = {0., 0., 0., 0.};
    std::array<double,4> p_before = {0., 0., 0., 0.};
    std::array<double,4> p_after = {0., 0., 0., 0.};
    double temperature = 0.;
    double tau = 0.;
    double step = 0.;
    int in_medium = 0;
};

enum class ScatteringDecision {
    // Commit the sampled hard scattering and its recoil/hole source.
    Apply,
    // Return the sampled candidate to the caller without committing it.
    StopBeforeApply,
    // Reject this hard candidate, produce no recoil/hole, and keep advancing
    // the same source so the callback can inspect later candidates.
    VetoAndContinue
};

using ScatteringCallback = std::function<ScatteringDecision(const ScatteringCandidate&)>;
using PropagationStepCallback = std::function<void(const PropagationStep&)>;
using PartonCallbackFactory = std::function<std::pair<ScatteringCallback, PropagationStepCallback>(
    int parton_index, int pdg_id, int parent_index, int d1, int d2)>;

void propagate_segment(std::array<double,4> &p,
                       std::array<double,4> &pos,
                       double tof,
                       int id,
                       numrand &nr,
                       double kappa,
                       double alpha,
                       int tmethod,
                       int model,
                       int ebe_hydro,
                       bool compat_moliere_legacy_hydro,
                       const HydroProfile &hydro_profile,
                       std::default_random_engine &elastic_rng,
                       std::vector<Quench> &new_particles,
                       int &had_scattering,
                       std::array<double,4> &orient,
                       const PropagationStepCallback &step_callback = nullptr,
                       TransportState *transport_state = nullptr);

void propagate_segment_with_scattering_callback(std::array<double,4> &p,
                                                std::array<double,4> &pos,
                                                double tof,
                                                int id,
                                                numrand &nr,
                                                double kappa,
                                                double alpha,
                                                int tmethod,
                                                int model,
                                                int ebe_hydro,
                                                bool compat_moliere_legacy_hydro,
                                                const HydroProfile &hydro_profile,
                                                std::default_random_engine &elastic_rng,
                                                std::vector<Quench> &new_particles,
                                                int &had_scattering,
                                                std::array<double,4> &orient,
                                                const ScatteringCallback &callback,
                                                const PropagationStepCallback &step_callback = nullptr,
                                                TransportState *transport_state = nullptr);

void process_recoilers(std::vector<Quench> &new_particles,
                       numrand &nr,
                       double kappa,
                       double alpha,
                       int tmethod,
                       int model,
                       int ebe_hydro,
                       bool compat_moliere_legacy_hydro,
                       const HydroProfile &hydro_profile,
                       std::default_random_engine &elastic_rng,
                       std::vector<Quench> &recoiled);

void do_eloss(const std::vector<Parton> &partons,
              std::vector<Quench> &quenched,
              double xcre,
              double ycre,
              numrand &nr,
              double kappa,
              double alpha,
              int tmethod,
              int model,
              int ebe_hydro,
              bool compat_moliere_legacy_hydro,
              const HydroProfile &hydro_profile,
              std::default_random_engine &elastic_rng,
              std::vector<Quench> &recoiled,
              const PartonCallbackFactory &callback_factory = nullptr);

}
