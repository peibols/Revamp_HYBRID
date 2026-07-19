#include "HeavyQuarkEnergyLoss.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace heavy_quark {
namespace {

constexpr double kPi = 3.14159265358979323846;
constexpr double kHbarCGeVFm = 0.2;

double spatial_norm2(const std::array<double,4> &p) {
    return p[0] * p[0] + p[1] * p[1] + p[2] * p[2];
}

bool finite_four_vector(const std::array<double,4> &p) {
    return std::all_of(p.begin(), p.end(), [](double value) { return std::isfinite(value); });
}

std::array<double,4> boost_to_fluid(const std::array<double,4> &p,
                                    const std::array<double,4> &v) {
    const double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if (v2 <= 0.) return p;
    const double bounded_v2 = std::min(v2, 1. - 1.e-12);
    const double gamma = 1. / std::sqrt(1. - bounded_v2);
    const double vdotp = v[0] * p[0] + v[1] * p[1] + v[2] * p[2];
    const double coefficient = (gamma - 1.) * vdotp / bounded_v2 - gamma * p[3];
    return {p[0] + coefficient * v[0],
            p[1] + coefficient * v[1],
            p[2] + coefficient * v[2],
            gamma * (p[3] - vdotp)};
}

std::array<double,4> boost_from_fluid(const std::array<double,4> &p,
                                      const std::array<double,4> &v) {
    const double v2 = v[0] * v[0] + v[1] * v[1] + v[2] * v[2];
    if (v2 <= 0.) return p;
    const double bounded_v2 = std::min(v2, 1. - 1.e-12);
    const double gamma = 1. / std::sqrt(1. - bounded_v2);
    const double vdotp = v[0] * p[0] + v[1] * p[1] + v[2] * p[2];
    const double coefficient = (gamma - 1.) * vdotp / bounded_v2 + gamma * p[3];
    return {p[0] + coefficient * v[0],
            p[1] + coefficient * v[1],
            p[2] + coefficient * v[2],
            gamma * (p[3] + vdotp)};
}

double effective_mass(const std::array<double,4> &p, double mass_floor) {
    const double invariant_mass2 = p[3] * p[3] - spatial_norm2(p);
    if (invariant_mass2 <= mass_floor * mass_floor) return mass_floor;
    return std::sqrt(invariant_mass2);
}

void put_on_shell(std::array<double,4> &p, double mass) {
    p[3] = std::sqrt(std::max(0., spatial_norm2(p) + mass * mass));
}

double open_uniform(numrand &nr) {
    double value = nr.rando();
    while (value <= 0.) value = nr.rando();
    return value;
}

double normal_draw(numrand &nr, double variance) {
    if (variance <= 0.) return 0.;
    const double radius = std::sqrt(-2. * std::log(open_uniform(nr)) * variance);
    return radius * std::cos(2. * kPi * nr.rando());
}

void apply_diffusion(std::array<double,4> &p_fluid, double mass,
                     double temperature, double path_length_fm,
                     double t_hooft_lambda, numrand &nr) {
    if (t_hooft_lambda <= 0. || temperature <= 0. || path_length_fm <= 0.) return;
    put_on_shell(p_fluid, mass);
    const double gamma = std::max(1., p_fluid[3] / mass);
    // Legacy convention: kappa_HQ = pi*sqrt(lambda), with
    // <Delta p_i^2> = kappa_HQ*gamma*T^3*Delta x for each component.
    const double kappa_hq = kPi * std::sqrt(t_hooft_lambda);
    const double path_length_inverse_gev = path_length_fm / kHbarCGeVFm;
    const double variance = kappa_hq * gamma * std::pow(temperature, 3.) *
                            path_length_inverse_gev;
    p_fluid[0] += normal_draw(nr, variance);
    p_fluid[1] += normal_draw(nr, variance);
    p_fluid[2] += normal_draw(nr, variance);
    put_on_shell(p_fluid, mass);
}

void record_result(const StepResult &result, Diagnostics *diagnostics) {
    if (diagnostics == nullptr || result.decision == StepDecision::NotApplicable) return;
    ++diagnostics->n_heavy_steps;
    switch (result.decision) {
        case StepDecision::UseBaseline:
            ++diagnostics->n_baseline_steps;
            break;
        case StepDecision::AppliedDrag:
            ++diagnostics->n_drag_steps;
            break;
        case StepDecision::AppliedDragAndDiffusion:
            ++diagnostics->n_drag_steps;
            ++diagnostics->n_diffusion_steps;
            break;
        case StepDecision::AppliedDiffusionOnly:
            ++diagnostics->n_diffusion_steps;
            ++diagnostics->n_diffusion_only_steps;
            break;
        case StepDecision::NotApplicable:
            break;
    }
    if (std::isfinite(result.energy_before) && std::isfinite(result.energy_after)) {
        diagnostics->sum_energy_change += result.energy_before - result.energy_after;
    }
}

}  // namespace

Parameters make_parameters(int mode, double t_hooft_lambda,
                           double charm_mass, double bottom_mass) {
    if (mode < 0 || mode > 3) {
        throw std::invalid_argument("heavy_quark_eloss_mode must be in [0,3]");
    }
    if (!std::isfinite(t_hooft_lambda) || t_hooft_lambda < 0.) {
        throw std::invalid_argument("heavy_quark_lambda must be finite and non-negative");
    }
    if (!std::isfinite(charm_mass) || charm_mass <= 0. ||
        !std::isfinite(bottom_mass) || bottom_mass <= 0.) {
        throw std::invalid_argument("heavy-quark masses must be finite and positive");
    }
    Parameters parameters;
    parameters.mode = static_cast<Mode>(mode);
    parameters.t_hooft_lambda = t_hooft_lambda;
    parameters.charm_mass = charm_mass;
    parameters.bottom_mass = bottom_mass;
    return parameters;
}

void validate_for_energy_loss_model(const Parameters &parameters, int energy_loss_model) {
    if ((parameters.mode == Mode::Drag || parameters.mode == Mode::DragAndDiffusion) &&
        energy_loss_model != 0) {
        throw std::invalid_argument(
            "heavy-quark drag modes 1 and 2 require the HYBRID strong-coupling mode=0");
    }
}

std::string mode_name(Mode mode) {
    switch (mode) {
        case Mode::Disabled: return "disabled";
        case Mode::Drag: return "drag_crossover";
        case Mode::DragAndDiffusion: return "drag_crossover_plus_diffusion";
        case Mode::DiffusionOnly: return "diffusion_only";
    }
    return "invalid";
}

bool is_heavy_quark(int pdg_id) {
    const int abs_id = std::abs(pdg_id);
    return abs_id == 4 || abs_id == 5;
}

double mass_for_pdg(int pdg_id, const Parameters &parameters) {
    return std::abs(pdg_id) == 5 ? parameters.bottom_mass : parameters.charm_mass;
}

bool applies_heavy_update(StepDecision decision) {
    return decision == StepDecision::AppliedDrag ||
           decision == StepDecision::AppliedDragAndDiffusion ||
           decision == StepDecision::AppliedDiffusionOnly;
}

StepResult apply_step(std::array<double,4> &p_lab,
                      const std::array<double,4> &fluid_velocity,
                      const StepInput &input,
                      const Parameters &parameters,
                      numrand &nr,
                      Diagnostics *diagnostics) {
    StepResult result;
    result.energy_before = p_lab[3];
    result.energy_after = p_lab[3];
    if (parameters.mode == Mode::Disabled || !is_heavy_quark(input.pdg_id)) return result;

    const double mass_floor = mass_for_pdg(input.pdg_id, parameters);
    const double fluid_v2 = fluid_velocity[0] * fluid_velocity[0] +
                            fluid_velocity[1] * fluid_velocity[1] +
                            fluid_velocity[2] * fluid_velocity[2];
    if (!finite_four_vector(p_lab) || p_lab[3] <= 0. ||
        !finite_four_vector(fluid_velocity) || fluid_v2 >= 1. ||
        !std::isfinite(input.temperature) || input.temperature < 0. ||
        !std::isfinite(input.fluid_path_length_fm) || input.fluid_path_length_fm < 0. ||
        (input.baseline_available &&
         (std::isnan(input.baseline_energy_loss_fluid) ||
          input.baseline_energy_loss_fluid < 0.))) {
        result.decision = StepDecision::UseBaseline;
        if (diagnostics != nullptr) ++diagnostics->n_invalid_steps;
        record_result(result, diagnostics);
        return result;
    }
    result.effective_mass = effective_mass(p_lab, mass_floor);

    if (parameters.t_hooft_lambda == 0. || input.temperature == 0. ||
        input.fluid_path_length_fm == 0.) {
        if (parameters.mode == Mode::DiffusionOnly) {
            result.decision = StepDecision::AppliedDiffusionOnly;
        } else if (parameters.mode == Mode::DragAndDiffusion) {
            result.decision = StepDecision::AppliedDragAndDiffusion;
        } else {
            result.decision = StepDecision::AppliedDrag;
        }
        record_result(result, diagnostics);
        return result;
    }

    if (parameters.mode == Mode::DiffusionOnly) {
        auto p_fluid = boost_to_fluid(p_lab, fluid_velocity);
        put_on_shell(p_fluid, result.effective_mass);
        apply_diffusion(p_fluid, result.effective_mass, input.temperature,
                        input.fluid_path_length_fm, parameters.t_hooft_lambda, nr);
        p_lab = boost_from_fluid(p_fluid, fluid_velocity);
        result.decision = StepDecision::AppliedDiffusionOnly;
        result.energy_after = p_lab[3];
        record_result(result, diagnostics);
        return result;
    }

    if (result.effective_mass <= input.temperature || input.fluid_path_length_fm <= 0.) {
        result.decision = StepDecision::UseBaseline;
        record_result(result, diagnostics);
        return result;
    }

    auto p_fluid_before = boost_to_fluid(p_lab, fluid_velocity);
    put_on_shell(p_fluid_before, result.effective_mass);
    auto p_fluid_drag = p_fluid_before;
    const double eta_drag = 0.5 * kPi * std::sqrt(parameters.t_hooft_lambda) *
                            input.temperature * input.temperature / result.effective_mass;
    const double path_length_inverse_gev = input.fluid_path_length_fm / kHbarCGeVFm;
    // This is the legacy first-order drag update for the spatial momentum.
    // Recomputing E from the mass shell removes its O(step^2) inconsistency.
    const double drag_factor = std::max(0., 1. - eta_drag * path_length_inverse_gev);
    p_fluid_drag[0] *= drag_factor;
    p_fluid_drag[1] *= drag_factor;
    p_fluid_drag[2] *= drag_factor;
    put_on_shell(p_fluid_drag, result.effective_mass);
    const double drag_energy_loss = std::max(0., p_fluid_before[3] - p_fluid_drag[3]);

    const bool baseline_would_cross_mass =
        input.baseline_available &&
        p_fluid_before[3] - input.baseline_energy_loss_fluid < result.effective_mass;
    const bool use_drag = !input.baseline_available || baseline_would_cross_mass ||
                          input.baseline_energy_loss_fluid > drag_energy_loss;
    if (!use_drag) {
        result.decision = StepDecision::UseBaseline;
        record_result(result, diagnostics);
        return result;
    }

    if (parameters.mode == Mode::DragAndDiffusion) {
        apply_diffusion(p_fluid_drag, result.effective_mass, input.temperature,
                        input.fluid_path_length_fm, parameters.t_hooft_lambda, nr);
        result.decision = StepDecision::AppliedDragAndDiffusion;
    } else {
        result.decision = StepDecision::AppliedDrag;
    }
    p_lab = boost_from_fluid(p_fluid_drag, fluid_velocity);
    result.energy_after = p_lab[3];
    record_result(result, diagnostics);
    return result;
}

}  // namespace heavy_quark
