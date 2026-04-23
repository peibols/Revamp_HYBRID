#include <algorithm>
#include <array>
#include <cmath>
#include <fstream>
#include <iostream>
#include <vector>

#include "Parton.h"
#include "Quench.h"
#include "Random.h"
#include "Wake.h"

using std::vector;
using namespace std;

namespace {
constexpr int kProposalBudgetBase = 200000;
constexpr int kProposalBudgetStep = 20000;
constexpr int kMaxProposalRestarts = 100;
using Four = std::array<double,4>;

double transcut = 0.;
double basesig = 0.65;
int Nrun = 800000;
double maxptsq = 3.5 * 3.5;
double maxrap = 2.5;
double tole = 0.25;
double masspi = 0.1396;
double masspro = 0.938;
double masstra[2] = {masspi, masspro};
double normcoop[2] = {30., 30.};
double maxcooper[2] = {0., 0.};
int toomuch = 0;

Four to_array(const vector<double> &p) {
  return {p[0], p[1], p[2], p[3]};
}

vector<double> to_vector(const Four &p) {
  return {p[0], p[1], p[2], p[3]};
}

Four wake_array(const Wake &w) {
  return to_array(w.vGetP());
}

Four add(const Four &a, const Four &b) {
  return {a[0] + b[0], a[1] + b[1], a[2] + b[2], a[3] + b[3]};
}

Four sub(const Four &a, const Four &b) {
  return {a[0] - b[0], a[1] - b[1], a[2] - b[2], a[3] - b[3]};
}

Four scale(const Four &a, double c) {
  return {a[0] * c, a[1] * c, a[2] * c, a[3] * c};
}

void add_inplace(Four &a, const Four &b) {
  a[0] += b[0];
  a[1] += b[1];
  a[2] += b[2];
  a[3] += b[3];
}

double rapid(double pt, double pz) {
  if (pt == 0.) return 0.;
  return atanh(pz / sqrt(pt * pt + pz * pz));
}

int set_charge(int spe, numrand &nr) {
  if (spe == 0) {
    double r = nr.rando();
    if (r < 0.5) return 0;
    else if (r < 0.75) return 1;
    else return -1;
  } else {
    double r = nr.rando();
    if (r < 0.5) return 1;
    else return -1;
  }
}

double thermal(int spe, double ptrand) {
  double temp;
  if (spe == 0) {
    temp = 0.211501 * pow(ptrand, 0.275362);
    if (temp > 0.4) temp = 0.4;
    if (temp < 0.19) temp = 0.19;
  } else {
    temp = 0.33 * pow(ptrand, 0.3);
    if (temp > 0.4) temp = 0.4;
    if (temp < 0.15) temp = 0.15;
  }
  return temp;
}

Four vec_abs(const Four &p) {
  return {std::fabs(p[0]), std::fabs(p[1]), std::fabs(p[2]), std::fabs(p[3])};
}

void one_body(vector<Wake> &wake,
              const Four &delta,
              const Four &momback,
              double ptlost,
              double mtlost,
              double raplost,
              numrand &nr,
              int spe,
              int mode) {
  double mc = 0.;
  double cooper = 0.;
  double randian = 1.;
  double pxrand, pyrand, raprand, ptrand;
  double phirand, mtrand, phidif, rapdif;

  do {
    ptrand = std::max(sqrt(maxptsq * nr.rando()), 0.000001);
    phirand = 2. * 3.141592654 * nr.rando();
    raprand = maxrap * (-1. + 2. * nr.rando());
    pxrand = ptrand * cos(phirand);
    pyrand = ptrand * sin(phirand);
    mtrand = sqrt(ptrand * ptrand + masstra[spe] * masstra[spe]);

    double inang = (delta[0] * pxrand + delta[1] * pyrand) / (ptlost * ptrand);
    if (inang > 1.) inang = 1.;
    if (inang < -1.) inang = -1.;
    phidif = acos(inang);
    rapdif = raprand;

    double Temp = thermal(spe, ptrand);
    double cosh_rad = cosh(rapdif);
    double T5 = Temp * Temp * Temp * Temp * Temp;
    cooper = exp(-mtrand / Temp * cosh_rad) * mtrand / T5 * cosh_rad *
             (ptrand * 3. * ptlost / mtlost * cos(phidif) + mtrand * cosh_rad) /
             normcoop[spe];

    if (cooper != cooper) cooper = 0.;

    if (fabs(cooper) > maxcooper[spe]) maxcooper[spe] = fabs(cooper);
    if (fabs(cooper) > 1.) {
      normcoop[spe] *= (fabs(cooper) + 0.0001);
    }

    if (mode == -1) mc = fabs(cooper);
    else mc = cooper * wake[mode].GetStatus();

    randian = nr.rando();
  } while (mc < randian);

  double pt_cosh = ptrand * cosh(raprand + raplost);
  Four p = {pxrand, pyrand,
            ptrand * sinh(raprand + raplost),
            sqrt(pt_cosh * pt_cosh + masstra[spe] * masstra[spe])};

  int charge = set_charge(spe, nr);
  double status = (cooper > 0.) ? 1. : -1.;
  wake.emplace_back(to_vector(p), masstra[spe], charge, spe, status);
}

bool try_add_residual_wake(vector<Wake> &wake,
                           const Four &delta,
                           Four &momback,
                           numrand &nr) {
  double difx = -momback[0] + delta[0];
  double dify = -momback[1] + delta[1];
  double difz = -momback[2] + delta[2];
  double dife = -momback[3] + delta[3];
  double remass2 = dife * dife - difx * difx - dify * dify - difz * difz;

  if (!(fabs(difx) < 3. * tole && fabs(dify) < 3. * tole && fabs(difz) < 8. * tole && remass2 > 0.)) {
    return false;
  }

  int spe = 0;
  if (fabs(remass2 - masstra[0] * masstra[0]) >= fabs(remass2 - masstra[1] * masstra[1])) {
    spe = 1;
  }

  int charge = set_charge(spe, nr);
  double stat = (dife < 0.) ? -1. : 1.;
  double tdife = sqrt(difx * difx + dify * dify + difz * difz + masstra[spe] * masstra[spe]);
  Four p = {stat * difx, stat * dify, stat * difz, tdife};
  wake.emplace_back(to_vector(p), masstra[spe], charge, spe, stat);
  add_inplace(momback, scale(wake_array(wake.back()), stat));
  return true;
}
}  // namespace

void do_wake(vector<Quench> quenched,
             vector<Parton> partons,
             vector<Wake> &wake,
             numrand &nr,
             vector<vector<double>> all_wakes) {
  (void)all_wakes;

  double total_injected = 0.;
  int wake_source_idx = 0;
  for (size_t i = 0; i < quenched.size(); ++i) {
    if (partons[i].GetD1() != -1) continue;
    Four delta = sub(to_array(partons[i].vGetP()), to_array(quenched[i].vGetP()));

    total_injected += delta[3];

    double ptlost = sqrt(delta[0] * delta[0] + delta[1] * delta[1]);
    if (ptlost <= 0.) {
      ++wake_source_idx;
      continue;
    }
    double raplost = rapid(ptlost, delta[2]);
    if (raplost != raplost) {
      std::cout << " RapLost NaN: Dx= " << delta[0] << " Dy= " << delta[1] << " Dz= " << delta[2] << " De= " << delta[3] << std::endl;
      ++wake_source_idx;
      continue;
    }
    double mtlost = delta[3] / cosh(raplost);

    if (fabs(delta[3]) >= 0. && delta[3] < 10000000. && fabs(mtlost) > transcut) {
      std::vector<Wake> pwake;
      Four momback = {};
      Four pmomback = {};
      Four dif = {};
      double msigma, pass, newpass;
      int spe, mode;
      int runi, encallao;
      int numenc = 0;
      bool restart_metropolis = true;
      int restart_count = 0;

      while (restart_metropolis) {
        restart_metropolis = false;

        pwake.clear();
        momback = {};
        runi = 0; encallao = 0;
        int proposal_count = 0;
        const int proposal_budget = kProposalBudgetBase + restart_count * kProposalBudgetStep;

        do {
          if (nr.rando() <= 0.05) spe = 1;
          else spe = 0;
          one_body(pwake, delta, momback, ptlost, mtlost, raplost, nr, spe, -1);
          add_inplace(momback, scale(wake_array(pwake.back()), pwake.back().GetStatus()));
        } while (fabs(momback[3]) < fabs(delta[3]));

        dif = vec_abs(sub(delta, momback));

        double dif_sq = dif[0]*dif[0] + dif[1]*dif[1] + dif[2]*dif[2] + dif[3]*dif[3];
        msigma = sqrt(dif_sq) / sqrt(log(2.));
        double msigma2 = msigma * msigma;

        do {
          mode = int(double(pwake.size()) * nr.rando());
          spe = pwake[mode].GetId();

          one_body(pwake, delta, momback, ptlost, mtlost, raplost, nr, spe, mode);

          pass = exp(-dif_sq / msigma2);

          pmomback = add(momback, scale(sub(wake_array(pwake.back()), wake_array(pwake[mode])), pwake[mode].GetStatus()));

          dif = vec_abs(sub(delta, pmomback));
          dif_sq = dif[0]*dif[0] + dif[1]*dif[1] + dif[2]*dif[2] + dif[3]*dif[3];

          newpass = exp(-dif_sq / msigma2);

          if (newpass > pass) {
            std::swap(pwake[mode], pwake.back());
            pwake.pop_back();
            runi += 1;
            encallao = 0;
            momback = pmomback;
          } else {
            pwake.pop_back();
            dif = vec_abs(sub(delta, momback));
            dif_sq = dif[0]*dif[0] + dif[1]*dif[1] + dif[2]*dif[2] + dif[3]*dif[3];
            encallao += 1;
          }

          if (encallao > 50000) {
            std::cout << " ENCALLAO !!!! \n \n";
            numenc += 1;
            if (numenc > 5) {
              toomuch += 1;
              break;
            }
            restart_metropolis = true;
            break;
          }

          ++proposal_count;
          if (proposal_count >= proposal_budget) {
            if (try_add_residual_wake(pwake, delta, momback, nr)) {
              dif = vec_abs(sub(delta, momback));
              break;
            }

            ++restart_count;
            if (restart_count > kMaxProposalRestarts) {
              toomuch += 1;
              break;
            }
            restart_metropolis = true;
            break;
          }

        } while (runi < Nrun && (dif[0] > tole || dif[1] > tole || dif[2] > tole || dif[3] > tole));

        if (!restart_metropolis) {
          if (runi >= Nrun) toomuch += 1;

          if (dif[0] != 0.) {
            try_add_residual_wake(pwake, delta, momback, nr);
            dif = vec_abs(sub(delta, momback));
          }

          for (size_t k = 0; k < pwake.size(); ++k) {
            pwake[k].SetMom(wake_source_idx);
            int pdg = -1000;
            double charge = pwake[k].GetCharge();
            int particle_spe = pwake[k].GetId();
            if (particle_spe == 0) {
              if (charge == 0) pdg = 111;
              else if (charge == 1.) pdg = 211;
              else if (charge == -1.) pdg = -211;
            } else if (particle_spe == 1) {
              if (charge == 1.) pdg = 2212;
              else if (charge == -1.) pdg = -2212;
            }
            pwake[k].SetId(pdg);
            wake.push_back(pwake[k]);
          }
          pwake.clear();
        }
      }
    }
    ++wake_source_idx;
  }

  std::ofstream outfile;
  outfile.open("Injected.dat", std::ios_base::app);
  outfile << total_injected << endl;
}
