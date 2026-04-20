#include "WakeGenerator.h"

#include <array>
#include <cmath>
#include <iostream>
#include "vector_operators.h"

using std::vector;

namespace {
// The old wake Metropolis loop used wall-clock time to trigger restarts,
// which made same-seed runs diverge under different machine load.
constexpr int kProposalBudgetBase = 200000;
constexpr int kProposalBudgetStep = 20000;
constexpr int kMaxProposalRestarts = 100;
}

WakeGenerator::WakeGenerator()
    : transcut_(0.),
      basesig_(0.65),
      Nrun_(800000),
      maxptsq_(3.5 * 3.5),
      maxrap_(2.5),
      tole_(0.25),
      masspi_(0.1396),
      masspro_(0.938),
      toomuch_(0) {
    masstra_[0] = masspi_;
    masstra_[1] = masspro_;
    normcoop_[0] = 30.;
    normcoop_[1] = 30.;
    maxcooper_[0] = 0.;
    maxcooper_[1] = 0.;
}

WakeGenerator::~WakeGenerator() = default;

double WakeGenerator::rapid(double pt, double pz) const {
    if (pt == 0.) return 0.;
    return atanh(pz / sqrt(pt * pt + pz * pz));
}

int WakeGenerator::set_charge(int spe, numrand &nr) const {
    // spe: 0 = pion, 1 = proton
    if (spe == 0) {
        double r = nr.rando();
        if (r < 0.5) return 0;      // pi0
        else if (r < 0.75) return 1;  // pi+
        else return -1;                // pi-
    } else {
        double r = nr.rando();
        if (r < 0.5) return 1;         // proton
        else return -1;                // antiproton
    }
}

double WakeGenerator::thermal(int spe, double ptrand) const {
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

void WakeGenerator::one_body(vector<Wake> &wake, 
                             const std::array<double,4>& delta, 
                             const std::array<double,4>& momback, 
                             double ptlost, double mtlost, double raplost, 
                             numrand &nr, int spe, int mode) {
    double mc = 0.;
    double cooper = 0.;
    double randian = 1.;
    double pxrand, pyrand, raprand, ptrand;
    double phirand, mtrand, phidif, rapdif;
    
    do {
        ptrand = std::max(sqrt(maxptsq_ * nr.rando()), 0.000001);
        phirand = 2. * 3.141592654 * nr.rando();
        raprand = maxrap_ * (-1. + 2. * nr.rando());
        pxrand = ptrand * cos(phirand);
        pyrand = ptrand * sin(phirand);
        mtrand = sqrt(ptrand*ptrand + masstra_[spe]*masstra_[spe]);
        
        double inang = (delta[0] * pxrand + delta[1] * pyrand) / (ptlost * ptrand);
        if (inang > 1.) inang = 1.;
        if (inang < -1.) inang = -1.;
        phidif = acos(inang);
        rapdif = raprand;
        
        double Temp = thermal(spe, ptrand);
        // mtrand already computed above, no need to recompute
        
        double cosh_rad = cosh(rapdif);
        double T5 = Temp*Temp*Temp*Temp*Temp;
        cooper = exp(-mtrand / Temp * cosh_rad) * mtrand / T5 * cosh_rad *
                 (ptrand * 3. * ptlost / mtlost * cos(phidif) + mtrand * cosh_rad) / normcoop_[spe];
        
        if (cooper != cooper) {
            cooper = 0.;
        }
        
        if (fabs(cooper) > maxcooper_[spe]) maxcooper_[spe] = fabs(cooper);
        if (fabs(cooper) > 1.) {
            normcoop_[spe] *= (fabs(cooper) + 0.0001);
        }
        
        if (mode == -1) mc = fabs(cooper);
        else mc = cooper * wake[mode].GetStatus();
        
        randian = nr.rando();
    } while (mc < randian);
    
    // Create wake particle
    double pt_cosh = ptrand * cosh(raprand + raplost);
    std::array<double,4> p = {pxrand, pyrand,
                              ptrand * sinh(raprand + raplost),
                              sqrt(pt_cosh*pt_cosh + masstra_[spe]*masstra_[spe])};
    
    int charge = set_charge(spe, nr);
    double status;
    if (cooper > 0.) status = 1.;
    else status = -1.;
    wake.emplace_back(p, masstra_[spe], charge, spe, status);
}

bool WakeGenerator::tryAddResidualWake(std::vector<Wake> &wake,
                                       const std::array<double,4>& delta,
                                       std::array<double,4>& momback,
                                       numrand &nr) const {
    double difx = -momback[0] + delta[0];
    double dify = -momback[1] + delta[1];
    double difz = -momback[2] + delta[2];
    double dife = -momback[3] + delta[3];
    double remass2 = dife*dife - difx*difx - dify*dify - difz*difz;

    if (!(fabs(difx) < 3. * tole_ && fabs(dify) < 3. * tole_ && fabs(difz) < 8. * tole_ && remass2 > 0.)) {
        return false;
    }

    int spe = 0;
    if (fabs(remass2 - masstra_[0] * masstra_[0]) >= fabs(remass2 - masstra_[1] * masstra_[1])) {
        spe = 1;
    }

    int charge = set_charge(spe, nr);
    double stat = (dife < 0.) ? -1. : 1.;
    double tdife = sqrt(difx*difx + dify*dify + difz*difz + masstra_[spe]*masstra_[spe]);
    std::array<double,4> p = {stat * difx, stat * dify, stat * difz, tdife};
    wake.push_back(Wake(p, masstra_[spe], charge, spe, stat));
    momback += wake.back().vGetP() * stat;
    return true;
}

void WakeGenerator::generate(const std::vector<Quench> &quenched, 
                             const std::vector<Parton> &partons, 
                             std::vector<Wake> &wake, 
                             numrand &nr) {
    int wake_source_idx = 0;  // sequential index of processed final partons
    for (size_t i = 0; i < quenched.size(); ++i) {
        if (partons[i].GetD1() != -1) continue;
        std::array<double,4> delta = partons[i].vGetP() - quenched[i].vGetP();

        double ptlost = sqrt(delta[0]*delta[0] + delta[1]*delta[1]);
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

        if (fabs(delta[3]) >= 0. && delta[3] < 10000000. && fabs(mtlost) > transcut_) {
            std::vector<Wake> pwake;
            std::array<double,4> momback = {};
            std::array<double,4> pmomback = {};
            std::array<double,4> dif = {};
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
                    momback += pwake.back().vGetP() * pwake.back().GetStatus();
                } while (fabs(momback[3]) < fabs(delta[3]));

                dif = vec_abs(delta - momback);

                double dif_sq = dif[0]*dif[0] + dif[1]*dif[1] + dif[2]*dif[2] + dif[3]*dif[3];
                msigma = sqrt(dif_sq) / sqrt(log(2.));
                double msigma2 = msigma * msigma;

                do {
                    mode = int(double(pwake.size()) * nr.rando());
                    spe = pwake[mode].GetId();

                    one_body(pwake, delta, momback, ptlost, mtlost, raplost, nr, spe, mode);

                    pass = exp(-dif_sq / msigma2);

                    pmomback = momback + (pwake.back().vGetP() - pwake[mode].vGetP()) * pwake[mode].GetStatus();

                    dif = vec_abs(delta - pmomback);
                    dif_sq = dif[0]*dif[0] + dif[1]*dif[1] + dif[2]*dif[2] + dif[3]*dif[3];

                    newpass = exp(-dif_sq / msigma2);

                    if (newpass > pass) {
                        // O(1) erase: swap target with back, then pop
                        std::swap(pwake[mode], pwake.back());
                        pwake.pop_back();
                        runi += 1;
                        encallao = 0;
                        momback = pmomback;
                    } else {
                        pwake.pop_back();
                        // restore dif to pre-proposal state
                        dif = vec_abs(delta - momback);
                        dif_sq = dif[0]*dif[0] + dif[1]*dif[1] + dif[2]*dif[2] + dif[3]*dif[3];
                        encallao += 1;
                    }

                    if (encallao > 50000) {
                        std::cout << " ENCALLAO !!!! \n \n";
                        numenc += 1;
                        if (numenc > 5) {
                            toomuch_ += 1;
                            break;
                        }
                        restart_metropolis = true;
                        break;
                    }

                    ++proposal_count;
                    if (proposal_count >= proposal_budget) {
                        if (tryAddResidualWake(pwake, delta, momback, nr)) {
                            dif = vec_abs(delta - momback);
                            break;
                        }

                        ++restart_count;
                        if (restart_count > kMaxProposalRestarts) {
                            toomuch_ += 1;
                            break;
                        }
                        restart_metropolis = true;
                        break;
                    }

                } while (runi < Nrun_ && (dif[0] > tole_ || dif[1] > tole_ || dif[2] > tole_ || dif[3] > tole_));

                if (!restart_metropolis) {
                    if (runi >= Nrun_) toomuch_ += 1;

                    if (dif[0] != 0.) {
                        tryAddResidualWake(pwake, delta, momback, nr);
                        dif = vec_abs(delta - momback);
                    }

                    for (size_t k = 0; k < pwake.size(); k++) {
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
}

std::array<double,4> WakeGenerator::vec_abs(const std::array<double,4>& p) const {
    return {std::fabs(p[0]), std::fabs(p[1]), std::fabs(p[2]), std::fabs(p[3])};
}
