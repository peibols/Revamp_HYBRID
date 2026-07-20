#include "HeavyQuarkEnergyLoss.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <stdexcept>

namespace {

bool close(double a, double b, double tolerance = 1.e-11) {
    return std::abs(a - b) <= tolerance * std::max({1., std::abs(a), std::abs(b)});
}

bool same(const std::array<double,4> &a, const std::array<double,4> &b,
          double tolerance = 1.e-11) {
    for (int i = 0; i < 4; ++i) {
        if (!close(a[i], b[i], tolerance)) return false;
    }
    return true;
}

double invariant_mass2(const std::array<double,4> &p) {
    return p[3] * p[3] - p[0] * p[0] - p[1] * p[1] - p[2] * p[2];
}

std::array<double,4> on_shell(double px, double py, double pz, double mass) {
    return {px, py, pz, std::sqrt(px * px + py * py + pz * pz + mass * mass)};
}

}  // namespace

int main() {
    const std::array<double,4> static_fluid = {0., 0., 0., 1.};
    heavy_quark::StepInput charm_step;
    charm_step.pdg_id = 4;
    charm_step.temperature = 0.3;
    charm_step.fluid_path_length_fm = 0.1;

    {
        auto p = on_shell(10., 1., 0.5, 1.25);
        const auto before = p;
        numrand tested_rng(1234);
        numrand reference_rng(1234);
        heavy_quark::Parameters disabled;
        const auto result = heavy_quark::apply_step(
            p, static_fluid, charm_step, disabled, tested_rng);
        assert(result.decision == heavy_quark::StepDecision::NotApplicable);
        assert(same(p, before));
        assert(tested_rng.rando() == reference_rng.rando());
    }

    {
        auto p = on_shell(10., 1., 0.5, 0.);
        const auto before = p;
        auto parameters = heavy_quark::make_parameters(2, 5.5);
        auto light_step = charm_step;
        light_step.pdg_id = 1;
        numrand tested_rng(55);
        numrand reference_rng(55);
        const auto result = heavy_quark::apply_step(
            p, static_fluid, light_step, parameters, tested_rng);
        assert(result.decision == heavy_quark::StepDecision::NotApplicable);
        assert(same(p, before));
        assert(tested_rng.rando() == reference_rng.rando());
    }

    {
        auto p = on_shell(10., 1., 0.5, 1.25);
        const auto before = p;
        auto parameters = heavy_quark::make_parameters(1, 5.5);
        heavy_quark::Diagnostics diagnostics;
        numrand nr(10);
        const auto result = heavy_quark::apply_step(
            p, static_fluid, charm_step, parameters, nr, &diagnostics);
        assert(result.decision == heavy_quark::StepDecision::AppliedDrag);
        assert(p[3] < before[3]);
        assert(close(invariant_mass2(p), result.effective_mass * result.effective_mass));
        assert(diagnostics.n_drag_steps == 1);
        assert(diagnostics.n_diffusion_steps == 0);

        const double eta_drag = 0.5 * 3.14159265358979323846 * std::sqrt(5.5) *
                                charm_step.temperature * charm_step.temperature / 1.25;
        const double expected_px = before[0] *
                                   (1. - eta_drag * charm_step.fluid_path_length_fm / 0.2);
        assert(close(p[0], expected_px));
    }

    {
        auto p = on_shell(10., 1., 0.5, 1.25);
        const auto before = p;
        auto parameters = heavy_quark::make_parameters(1, 5.5);
        auto baseline_step = charm_step;
        baseline_step.baseline_available = true;
        baseline_step.baseline_energy_loss_fluid = 0.;
        numrand tested_rng(17);
        numrand reference_rng(17);
        const auto result = heavy_quark::apply_step(
            p, static_fluid, baseline_step, parameters, tested_rng);
        assert(result.decision == heavy_quark::StepDecision::AppliedBaselineOnShell);
        assert(same(p, before));
        assert(tested_rng.rando() == reference_rng.rando());
    }

    {
        auto p = on_shell(10., 1., 0.5, 1.25);
        const auto before = p;
        auto parameters = heavy_quark::make_parameters(1, 5.5);
        auto baseline_step = charm_step;
        baseline_step.baseline_available = true;
        baseline_step.baseline_energy_loss_fluid = 0.05;
        numrand nr(18);
        const auto result = heavy_quark::apply_step(
            p, static_fluid, baseline_step, parameters, nr);
        assert(result.decision == heavy_quark::StepDecision::AppliedBaselineOnShell);
        assert(p[3] < before[3]);
        assert(close(before[3] - p[3], 0.05));
        assert(close(invariant_mass2(p), result.effective_mass * result.effective_mass));
    }

    {
        auto p1 = on_shell(10., 1., 0.5, 1.25);
        auto p2 = p1;
        auto parameters = heavy_quark::make_parameters(2, 5.5);
        numrand nr1(99);
        numrand nr2(99);
        const auto result1 = heavy_quark::apply_step(
            p1, static_fluid, charm_step, parameters, nr1);
        const auto result2 = heavy_quark::apply_step(
            p2, static_fluid, charm_step, parameters, nr2);
        assert(result1.decision == heavy_quark::StepDecision::AppliedDragAndDiffusion);
        assert(result2.decision == result1.decision);
        assert(same(p1, p2));
        assert(close(invariant_mass2(p1), result1.effective_mass * result1.effective_mass));
    }

    {
        auto p = on_shell(10., 1., 0.5, 1.25);
        const auto before = p;
        auto parameters = heavy_quark::make_parameters(2, 5.5);
        auto baseline_step = charm_step;
        baseline_step.baseline_available = true;
        baseline_step.baseline_energy_loss_fluid = 0.;
        numrand tested_rng(101);
        numrand reference_rng(101);
        const auto result = heavy_quark::apply_step(
            p, static_fluid, baseline_step, parameters, tested_rng);
        assert(result.decision == heavy_quark::StepDecision::AppliedBaselineOnShell);
        assert(same(p, before));
        assert(tested_rng.rando() == reference_rng.rando());
    }

    {
        auto p = on_shell(7., 0.5, 0.2, 1.25);
        const auto before = p;
        auto parameters = heavy_quark::make_parameters(3, 5.5);
        numrand nr(808);
        const auto result = heavy_quark::apply_step(
            p, static_fluid, charm_step, parameters, nr);
        assert(result.decision == heavy_quark::StepDecision::AppliedDiffusionOnly);
        assert(!same(p, before));
        assert(close(invariant_mass2(p), result.effective_mass * result.effective_mass));
    }

    {
        constexpr int samples = 50000;
        const auto before = on_shell(7., 0.5, 0.2, 1.25);
        const auto parameters = heavy_quark::make_parameters(3, 5.5);
        const double gamma = before[3] / 1.25;
        const double expected_variance = 3.14159265358979323846 * std::sqrt(5.5) *
                                         gamma * std::pow(charm_step.temperature, 3.) *
                                         charm_step.fluid_path_length_fm / 0.2;
        std::array<double,3> sum = {0., 0., 0.};
        std::array<double,3> sum2 = {0., 0., 0.};
        numrand nr(809);
        for (int sample = 0; sample < samples; ++sample) {
            auto p = before;
            const auto result = heavy_quark::apply_step(
                p, static_fluid, charm_step, parameters, nr);
            assert(result.decision == heavy_quark::StepDecision::AppliedDiffusionOnly);
            for (int component = 0; component < 3; ++component) {
                const double kick = p[component] - before[component];
                sum[component] += kick;
                sum2[component] += kick * kick;
            }
        }
        for (int component = 0; component < 3; ++component) {
            const double mean = sum[component] / samples;
            const double variance = sum2[component] / samples - mean * mean;
            assert(std::abs(mean) < 0.02 * std::sqrt(expected_variance));
            assert(std::abs(variance / expected_variance - 1.) < 0.03);
        }
    }

    {
        auto p = on_shell(7., 0.5, 0.2, 1.25);
        const auto before = p;
        auto parameters = heavy_quark::make_parameters(3, 0.);
        numrand tested_rng(808);
        numrand reference_rng(808);
        const auto result = heavy_quark::apply_step(
            p, static_fluid, charm_step, parameters, tested_rng);
        assert(result.decision == heavy_quark::StepDecision::AppliedDiffusionOnly);
        assert(same(p, before));
        assert(tested_rng.rando() == reference_rng.rando());
    }

    {
        auto p = on_shell(14., 2., 1., 4.2);
        const std::array<double,4> moving_fluid = {0.25, -0.1, 0.05, 1.};
        auto parameters = heavy_quark::make_parameters(1, 5.5);
        auto bottom_step = charm_step;
        bottom_step.pdg_id = 5;
        numrand nr(42);
        const auto result = heavy_quark::apply_step(
            p, moving_fluid, bottom_step, parameters, nr);
        assert(result.decision == heavy_quark::StepDecision::AppliedDrag);
        assert(result.effective_mass >= 4.2);
        assert(close(invariant_mass2(p), result.effective_mass * result.effective_mass, 1.e-10));
    }

    {
        bool bad_mode = false;
        try {
            (void)heavy_quark::make_parameters(4, 5.5);
        } catch (const std::invalid_argument &) {
            bad_mode = true;
        }
        assert(bad_mode);

        bool bad_model = false;
        try {
            const auto parameters = heavy_quark::make_parameters(1, 5.5);
            heavy_quark::validate_for_energy_loss_model(parameters, 1);
        } catch (const std::invalid_argument &) {
            bad_model = true;
        }
        assert(bad_model);
    }

    {
        heavy_quark::Parameters disabled;
        assert(heavy_quark::apply_generic_broadening(4, disabled));

        const auto drag = heavy_quark::make_parameters(1, 5.5, 1.25, 4.2, false);
        assert(heavy_quark::apply_generic_broadening(4, drag));

        const auto matched_diffusion =
            heavy_quark::make_parameters(2, 5.5, 1.25, 4.2, false);
        assert(!heavy_quark::apply_generic_broadening(4, matched_diffusion));
        assert(!heavy_quark::apply_generic_broadening(5, matched_diffusion));
        assert(heavy_quark::apply_generic_broadening(1, matched_diffusion));

        const auto legacy_additive =
            heavy_quark::make_parameters(2, 5.5, 1.25, 4.2, true);
        assert(heavy_quark::apply_generic_broadening(4, legacy_additive));

        const auto hard_disabled =
            heavy_quark::make_parameters(2, 5.5, 1.25, 4.2, false, false);
        assert(!heavy_quark::apply_hard_moliere(4, hard_disabled));
        assert(!heavy_quark::apply_hard_moliere(5, hard_disabled));
        assert(heavy_quark::apply_hard_moliere(1, hard_disabled));
        assert(heavy_quark::apply_hard_moliere(21, hard_disabled));
        assert(heavy_quark::apply_hard_moliere(4, legacy_additive));
    }

    {
        auto p = on_shell(4., 1., 0., 1.25);
        const auto before = p;
        auto parameters = heavy_quark::make_parameters(1, 5.5);
        auto invalid_step = charm_step;
        invalid_step.temperature = std::nan("");
        heavy_quark::Diagnostics diagnostics;
        numrand tested_rng(909);
        numrand reference_rng(909);
        const auto result = heavy_quark::apply_step(
            p, static_fluid, invalid_step, parameters, tested_rng, &diagnostics);
        assert(result.decision == heavy_quark::StepDecision::UseBaseline);
        assert(same(p, before));
        assert(tested_rng.rando() == reference_rng.rando());
        assert(diagnostics.n_invalid_steps == 1);
        assert(std::isfinite(diagnostics.sum_energy_change));
    }

    {
        const auto projectile = on_shell(10., 1., 0.5, 1.25);
        const auto thermal = on_shell(0.2, -0.1, 0.3, 0.);
        const auto exact = heavy_quark::check_hard_scattering_kinematics(
            projectile, projectile, thermal, thermal, 1.25);
        assert(exact.finite);
        assert(exact.outgoing_on_mass_shell);
        assert(exact.four_momentum_conserved);

        auto massless_projectile = projectile;
        massless_projectile[3] = std::sqrt(
            massless_projectile[0] * massless_projectile[0] +
            massless_projectile[1] * massless_projectile[1] +
            massless_projectile[2] * massless_projectile[2]);
        const auto legacy_table_like = heavy_quark::check_hard_scattering_kinematics(
            projectile, massless_projectile, thermal, thermal, 1.25);
        assert(legacy_table_like.finite);
        assert(!legacy_table_like.outgoing_on_mass_shell);
        assert(!legacy_table_like.four_momentum_conserved);

        heavy_quark::Diagnostics diagnostics;
        heavy_quark::record_hard_scattering_check(legacy_table_like, &diagnostics);
        assert(diagnostics.n_hard_scattered_heavy == 1);
        assert(diagnostics.n_hard_heavy_mass_shell_failures == 1);
        assert(diagnostics.n_hard_heavy_momentum_closure_failures == 1);
        assert(diagnostics.max_hard_heavy_mass_shell_residual > 0.);
        assert(diagnostics.max_hard_heavy_momentum_closure_residual > 0.);
    }

    std::cout << "heavy-quark kernel tests passed" << std::endl;
    return 0;
}
