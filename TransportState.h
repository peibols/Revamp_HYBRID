#pragma once

#include <algorithm>
#include <cmath>

// State that must survive when a physical trajectory is split into scheduler
// windows.  The transport kernels initialize it once for each live object.
struct TransportState {
    bool initialized = false;
    double initial_energy = 0.;
    double fluid_distance = 0.;
    double lab_distance = 0.;
    double virtual_fluid_distance = 0.;

    void initialize(double energy) {
        if (initialized) return;
        initialized = true;
        initial_energy = std::max(0., energy);
    }
};

inline double transport_step_duration(double current_time, double end_time,
                                      double nominal_step = 0.1) {
    if (end_time <= current_time || nominal_step <= 0.) return 0.;
    return std::min(nominal_step, end_time - current_time);
}

// Integral of the strong-coupling stopping-law rate from zero to distance x.
// Using differences of this primitive makes a static-medium trajectory
// independent of scheduler window boundaries.
inline double strong_coupling_cumulative_loss(double initial_fluid_energy,
                                              double stopping_distance,
                                              double fluid_distance) {
    if (initial_fluid_energy <= 0. || stopping_distance <= 0. || fluid_distance <= 0.) {
        return 0.;
    }
    constexpr double pi = 3.14159265358979323846;
    const double u = std::clamp(fluid_distance / stopping_distance, 0., 1.);
    return initial_fluid_energy * (2. / pi) *
           (std::asin(u) - u * std::sqrt(std::max(0., 1. - u * u)));
}

inline double strong_coupling_step_loss(double initial_fluid_energy,
                                        double stopping_distance,
                                        double fluid_distance_before,
                                        double fluid_distance_after) {
    return std::max(0., strong_coupling_cumulative_loss(
                             initial_fluid_energy, stopping_distance, fluid_distance_after) -
                         strong_coupling_cumulative_loss(
                             initial_fluid_energy, stopping_distance, fluid_distance_before));
}

inline double fluid_distance_weighted_step(double fluid_distance_before,
                                           double fluid_distance_after,
                                           double lab_duration) {
    if (lab_duration <= 0.) return 0.;
    return 0.5 * (fluid_distance_before + fluid_distance_after) * lab_duration;
}
