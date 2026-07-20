#pragma once

#include <array>
#include <vector>

#include "HydroProfile.h"
#include "Parton.h"
#include "Quench.h"
#include "Random.h"

namespace moliere {

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
                       std::vector<Quench> &new_particles,
                       int &had_scattering,
                       std::array<double,4> &orient);

void process_recoilers(std::vector<Quench> &new_particles,
                       numrand &nr,
                       double kappa,
                       double alpha,
                       int tmethod,
                       int model,
                       int ebe_hydro,
                       bool compat_moliere_legacy_hydro,
                       const HydroProfile &hydro_profile,
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
              std::vector<Quench> &recoiled);

}
