#include "MoliereResolution.h"
#include "TransportState.h"

#include <array>
#include <cassert>
#include <cmath>
#include <iostream>

namespace {
bool close(double a, double b, double tolerance = 1.e-12) {
    return std::abs(a - b) <= tolerance;
}
}

int main() {
    assert(close(transport_step_duration(1., 1.), 0.));
    assert(close(transport_step_duration(1., 1.03), 0.03));
    assert(close(transport_step_duration(1., 1.10), 0.10));
    assert(close(transport_step_duration(1., 1.13), 0.10));

    TransportState state;
    state.initialize(100.);
    state.fluid_distance = 0.7;
    state.initialize(5.);
    assert(close(state.initial_energy, 100.));
    assert(close(state.fluid_distance, 0.7));

    // Static-medium loss increments must not depend on scheduler cuts.
    const double initial_energy = 100.;
    const double stopping_distance = 2.;
    const double whole_strong = strong_coupling_step_loss(
        initial_energy, stopping_distance, 0., 0.2);
    const double split_strong =
        strong_coupling_step_loss(initial_energy, stopping_distance, 0., 0.07) +
        strong_coupling_step_loss(initial_energy, stopping_distance, 0.07, 0.2);
    assert(close(whole_strong, split_strong));

    const double whole_weight = fluid_distance_weighted_step(0., 0.2, 0.2);
    const double split_weight =
        fluid_distance_weighted_step(0., 0.07, 0.07) +
        fluid_distance_weighted_step(0.07, 0.2, 0.13);
    assert(close(whole_weight, split_weight));
    assert(close(strong_coupling_step_loss(
                     initial_energy, stopping_distance, 0.2, 0.2), 0.));

    const std::array<double,4> origin = {0., 0., 0., 1.};
    const std::array<double,4> z_separated = {0., 0., 2., 1.};
    const std::array<double,4> source_x = {10., 0., 0., 10.};
    const std::array<double,4> source_z = {0., 0., 10., 10.};
    assert(close(moliere_resolution::transverse_separation(
                     origin, z_separated, source_x), 2.));
    assert(close(moliere_resolution::transverse_separation(
                     origin, z_separated, source_z), 0.));

    const double q = moliere_resolution::kHbarCGeVFm;
    assert(close(moliere_resolution::qperp_dperp(q, 1.), 1.));
    assert(!moliere_resolution::passes(0., 1., 0.));
    assert(!moliere_resolution::passes(1., 0., 0.));
    assert(moliere_resolution::passes(q, 1.01, 1.));

    // Rigidly rotate both the dipole and source from z/x to x/z.
    const std::array<double,4> x_separated = {2., 0., 0., 1.};
    assert(close(moliere_resolution::transverse_separation(
                     origin, x_separated, source_z), 2.));

    std::cout << "transport/resolution unit tests passed\n";
    return 0;
}
