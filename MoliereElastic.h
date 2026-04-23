#pragma once

#include <vector>

#include "HydroProfile.h"
#include "Parton.h"
#include "Quench.h"
#include "Random.h"

namespace moliere {

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
              const HydroProfile &hydro_profile,
              std::vector<Quench> &recoiled);

}
