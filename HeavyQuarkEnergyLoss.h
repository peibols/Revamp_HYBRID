#pragma once

#include <array>
#include <string>

#include "Random.h"

namespace heavy_quark {

// Runtime modes are deliberately separate from the light-parton HYBRID
// `mode` option. Disabled performs no arithmetic and consumes no RNG state.
enum class Mode {
    Disabled = 0,
    Drag = 1,
    DragAndDiffusion = 2,
    DiffusionOnly = 3
};

struct Parameters {
    Mode mode = Mode::Disabled;
    double t_hooft_lambda = 0.;
    double charm_mass = 1.25;
    double bottom_mass = 4.2;
};

struct StepInput {
    int pdg_id = 0;
    double temperature = 0.;
    // Local-fluid-frame distance traversed during this integration step.
    double fluid_path_length_fm = 0.;
    // Candidate light-parton HYBRID loss used by the drag crossover.
    bool baseline_available = false;
    double baseline_energy_loss_fluid = 0.;
};

enum class StepDecision {
    NotApplicable,
    UseBaseline,
    AppliedDrag,
    AppliedDragAndDiffusion,
    AppliedDiffusionOnly
};

struct StepResult {
    StepDecision decision = StepDecision::NotApplicable;
    double effective_mass = 0.;
    double energy_before = 0.;
    double energy_after = 0.;
};

struct Diagnostics {
    long long n_heavy_steps = 0;
    long long n_baseline_steps = 0;
    long long n_drag_steps = 0;
    long long n_diffusion_steps = 0;
    long long n_diffusion_only_steps = 0;
    long long n_invalid_steps = 0;
    double sum_energy_change = 0.;
};

Parameters make_parameters(int mode, double t_hooft_lambda,
                           double charm_mass = 1.25,
                           double bottom_mass = 4.2);
void validate_for_energy_loss_model(const Parameters &parameters, int energy_loss_model);
std::string mode_name(Mode mode);
bool is_heavy_quark(int pdg_id);
double mass_for_pdg(int pdg_id, const Parameters &parameters);
bool applies_heavy_update(StepDecision decision);

StepResult apply_step(std::array<double,4> &p_lab,
                      const std::array<double,4> &fluid_velocity,
                      const StepInput &input,
                      const Parameters &parameters,
                      numrand &nr,
                      Diagnostics *diagnostics = nullptr);

}  // namespace heavy_quark
