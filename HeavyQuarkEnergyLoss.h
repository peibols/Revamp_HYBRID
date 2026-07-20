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
    bool add_generic_broadening_with_diffusion = true;
    bool enable_hard_moliere = true;
};

struct StepInput {
    int pdg_id = 0;
    double temperature = 0.;
    // Local-fluid-frame elapsed time during this integration step:
    // dt* = gamma_flow * (1 - v_flow dot v_parton) * dt.  The historical
    // field name says "path length"; it equals a spatial path only in the
    // ultrarelativistic limit.
    double fluid_path_length_fm = 0.;
    // Candidate light-parton HYBRID loss used by the drag crossover.
    bool baseline_available = false;
    double baseline_energy_loss_fluid = 0.;
};

enum class StepDecision {
    NotApplicable,
    UseBaseline,
    AppliedBaselineOnShell,
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
    long long n_hard_scattered_heavy = 0;
    long long n_hard_heavy_mass_shell_failures = 0;
    long long n_hard_heavy_momentum_closure_failures = 0;
    double max_hard_heavy_mass_shell_residual = 0.;
    double max_hard_heavy_momentum_closure_residual = 0.;
    double sum_energy_change = 0.;
};

struct HardScatteringCheck {
    double expected_mass2 = 0.;
    double outgoing_mass2 = 0.;
    double mass_shell_residual = 0.;
    double momentum_closure_residual = 0.;
    bool finite = false;
    bool outgoing_on_mass_shell = false;
    bool four_momentum_conserved = false;
};

Parameters make_parameters(int mode, double t_hooft_lambda,
                           double charm_mass = 1.25,
                           double bottom_mass = 4.2,
                           bool add_generic_broadening_with_diffusion = true,
                           bool enable_hard_moliere = true);
void validate_for_energy_loss_model(const Parameters &parameters, int energy_loss_model);
std::string mode_name(Mode mode);
bool is_heavy_quark(int pdg_id);
double mass_for_pdg(int pdg_id, const Parameters &parameters);
bool applies_heavy_update(StepDecision decision);
bool apply_generic_broadening(int pdg_id, const Parameters &parameters);
bool apply_hard_moliere(int pdg_id, const Parameters &parameters);

// Inspect the immediate 2->2 hard-scattering record before any later soft
// broadening or energy loss. The expected outgoing mass is the larger of the
// configured heavy-quark floor and the incoming invariant mass.
HardScatteringCheck check_hard_scattering_kinematics(
    const std::array<double,4> &projectile_before,
    const std::array<double,4> &projectile_after,
    const std::array<double,4> &recoiler_after,
    const std::array<double,4> &thermal_hole_before,
    double mass_floor);
void record_hard_scattering_check(const HardScatteringCheck &check,
                                  Diagnostics *diagnostics);

StepResult apply_step(std::array<double,4> &p_lab,
                      const std::array<double,4> &fluid_velocity,
                      const StepInput &input,
                      const Parameters &parameters,
                      numrand &nr,
                      Diagnostics *diagnostics = nullptr);

}  // namespace heavy_quark
