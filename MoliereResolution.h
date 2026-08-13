#pragma once

#include <array>
#include <algorithm>
#include <cmath>

namespace moliere_resolution {

// hbar*c converts q_perp [GeV] times d_perp [fm] to a dimensionless number.
constexpr double kHbarCGeVFm = 0.1973269804;

inline std::array<double,3> source_direction(const std::array<double,4> &p) {
    const double norm = std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    if (norm <= 0.) return {0., 0., 1.};
    return {p[0] / norm, p[1] / norm, p[2] / norm};
}

inline double transverse_separation(const std::array<double,4> &pos1,
                                    const std::array<double,4> &pos2,
                                    const std::array<double,4> &source_p) {
    const std::array<double,3> dr = {
        pos1[0] - pos2[0], pos1[1] - pos2[1], pos1[2] - pos2[2]};
    const auto n = source_direction(source_p);
    const double parallel = dr[0] * n[0] + dr[1] * n[1] + dr[2] * n[2];
    const double d2 = dr[0] * dr[0] + dr[1] * dr[1] + dr[2] * dr[2]
                    - parallel * parallel;
    return std::sqrt(std::max(0., d2));
}

inline double qperp_dperp(double qperp_gev, double dperp_fm) {
    if (qperp_gev <= 0. || dperp_fm <= 0.) return 0.;
    return qperp_gev * dperp_fm / kHbarCGeVFm;
}

inline bool passes(double qperp_gev, double dperp_fm, double c_res) {
    return qperp_dperp(qperp_gev, dperp_fm) > c_res;
}

}  // namespace moliere_resolution
