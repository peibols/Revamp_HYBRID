#include "WakeGenerator.h"

#include <cmath>
#include <iostream>
#include <ctime>
#include "vector_operators.h"

using std::vector;

WakeGenerator::WakeGenerator()
    : transcut_(0.),
      basesig_(0.65),
      Nrun_(800000),
      maxptsq_(pow(3.5, 2.)),
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
                             vector<double> delta, 
                             vector<double> momback, 
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
        mtrand = sqrt(pow(ptrand, 2.) + pow(masstra_[spe], 2.));
        
        double inang = (delta[0] * pxrand + delta[1] * pyrand) / (ptlost * ptrand);
        if (inang > 1.) inang = 1.;
        if (inang < -1.) inang = -1.;
        phidif = acos(inang);
        rapdif = raprand;
        
        double Temp = thermal(spe, ptrand);
        mtrand = sqrt(pow(ptrand, 2.) + pow(masstra_[spe], 2.));
        
        cooper = exp(-mtrand / Temp * cosh(rapdif)) * mtrand / pow(Temp, 5.) * cosh(rapdif) *
                 (ptrand * 3. * ptlost / mtlost * cos(phidif) + mtrand * cosh(rapdif)) / normcoop_[spe];
        
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
    vector<double> p(4);
    p[0] = pxrand;
    p[1] = pyrand;
    p[2] = ptrand * sinh(raprand + raplost);
    p[3] = sqrt(pow(ptrand * cosh(raprand + raplost), 2.) + pow(masstra_[spe], 2.));
    
    int charge = set_charge(spe, nr);
    double status;
    if (cooper > 0.) status = 1.;
    else status = -1.;
    wake.emplace_back(p, masstra_[spe], charge, spe, status);
}

void WakeGenerator::generate(const std::vector<Quench> &quenched, 
                             const std::vector<Parton> &partons, 
                             std::vector<Wake> &wake, 
                             numrand &nr) {
    // Compute all_wakes
    std::vector<std::vector<double>> all_wakes;
    for (size_t i = 0; i < quenched.size(); ++i) {
        if (partons[i].GetD1() != -1) continue;
        std::vector<double> delta = partons[i].vGetP() - quenched[i].vGetP();
        all_wakes.push_back(delta);
    }

    for (size_t i = 0; i < all_wakes.size(); ++i) {
        std::vector<double> delta = all_wakes[i];

        double ptlost = sqrt(pow(delta[0], 2.) + pow(delta[1], 2.));
        if (ptlost <= 0.) {
            //std::cout << " ptlost= " << ptlost << " elost= " << delta[3] << std::endl;
            continue;
        }
        double raplost = rapid(ptlost, delta[2]);
        if (raplost != raplost) {
            std::cout << " RapLost NaN: Dx= " << delta[0] << " Dy= " << delta[1] << " Dz= " << delta[2] << " De= " << delta[3] << std::endl;
            continue;
        }
        double mtlost = delta[3] / cosh(raplost);

        if (fabs(delta[3]) >= 0. && delta[3] < 10000000. && fabs(mtlost) > transcut_) {
            std::vector<Wake> pwake;
            std::vector<double> momback(4, 0.);
            std::vector<double> pmomback(4, 0.);
            std::vector<double> dif(4, 0.);
            std::vector<double> old_dif(4, 0.);
            double msigma, pass, newpass;
            int spe, mode;
            int runi, encallao;
            int numenc = 0;
            double clocklim = 0.1;
            int tooclock = 0;
            bool restart_metropolis = true;

            while (restart_metropolis) {
                restart_metropolis = false;  // Will be set to true if we need to restart
                clock_t startClock = clock();

                pwake.clear();
                for (unsigned int j = 0; j < 4; j++) momback[j] = 0.;
                runi = 0, encallao = 0;

                do {
                    if (nr.rando() <= 0.05) spe = 1;
                    else spe = 0;
                    one_body(pwake, delta, momback, ptlost, mtlost, raplost, nr, spe, -1);
                    momback += pwake[pwake.size() - 1].vGetP() * pwake[pwake.size() - 1].GetStatus();
                } while (fabs(momback[3]) < fabs(delta[3]));

                dif = vec_abs(delta - momback);

                msigma = sqrt(pow(dif[0], 2.) + pow(dif[1], 2.) + pow(dif[2], 2.) + pow(dif[3], 2.)) / sqrt(log(2.));

                do {
                    mode = int(double(pwake.size()) * nr.rando());

                    spe = pwake[mode].GetId();

                    one_body(pwake, delta, momback, ptlost, mtlost, raplost, nr, spe, mode);

                    pass = exp((-pow(dif[0], 2.) - pow(dif[1], 2.) - pow(dif[2], 2.) - pow(dif[3], 2.)) / pow(msigma, 2.));

                    pmomback = momback + (pwake[pwake.size() - 1].vGetP() - pwake[mode].vGetP()) * pwake[mode].GetStatus();

                    dif = vec_abs(delta - pmomback);

                    newpass = exp((-pow(dif[0], 2.) - pow(dif[1], 2.) - pow(dif[2], 2.) - pow(dif[3], 2.)) / pow(msigma, 2.));

                    if (newpass > pass) {
                        pwake.erase(pwake.begin() + mode);
                        runi += 1;
                        encallao = 0;
                        momback = pmomback;
                    } else {
                        pwake.erase(pwake.begin() + pwake.size() - 1);
                        dif = vec_abs(delta - momback);
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

                    clock_t endClock = clock();
                    if (double((endClock - startClock)) / CLOCKS_PER_SEC > clocklim) {
                        clocklim += 0.02;
                        tooclock += 1;

                        double difx = -momback[0] + delta[0];
                        double dify = -momback[1] + delta[1];
                        double difz = -momback[2] + delta[2];
                        double dife = -momback[3] + delta[3];
                        double remass2 = dife * dife - difx * difx - dify * dify - difz * difz;

                        if (fabs(difx) < 3. * tole_ && fabs(dify) < 3. * tole_ && fabs(difz) < 8. * tole_ && remass2 > 0.) {
                            if (fabs(remass2 - masstra_[0] * masstra_[0]) < fabs(remass2 - masstra_[1] * masstra_[1])) spe = 0;
                            else spe = 1;

                            int charge = set_charge(spe, nr);

                            double stat;
                            if (dife < 0.) stat = -1.;
                            else stat = 1.;

                            std::vector<double> p;
                            double tdife = sqrt(difx * difx + dify * dify + difz * difz + masstra_[spe] * masstra_[spe]);
                            p.push_back(stat * difx), p.push_back(stat * dify), p.push_back(stat * difz), p.push_back(tdife);
                            pwake.push_back(Wake(p, masstra_[spe], charge, spe, stat));
                            momback += pwake[pwake.size() - 1].vGetP() * stat;
                            dif = vec_abs(delta - momback);
                            break;
                        }

                        if (tooclock > 100) {
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
                        double difx = -momback[0] + delta[0];
                        double dify = -momback[1] + delta[1];
                        double difz = -momback[2] + delta[2];
                        double dife = -momback[3] + delta[3];
                        double remass2 = dife * dife - difx * difx - dify * dify - difz * difz;

                        if (fabs(remass2 - masstra_[0] * masstra_[0]) < fabs(remass2 - masstra_[1] * masstra_[1])) spe = 0;
                        else spe = 1;

                        int charge = set_charge(spe, nr);

                        double stat;
                        if (dife < 0.) stat = -1.;
                        else stat = 1.;

                        std::vector<double> p;
                        double tdife = sqrt(difx * difx + dify * dify + difz * difz + masstra_[spe] * masstra_[spe]);
                        p.push_back(stat * difx), p.push_back(stat * dify), p.push_back(stat * difz), p.push_back(tdife);
                        pwake.push_back(Wake(p, masstra_[spe], charge, spe, stat));
                        momback += pwake[pwake.size() - 1].vGetP() * stat;
                        dif = vec_abs(delta - momback);
                    }

                    for (size_t k = 0; k < pwake.size(); k++) {
                        pwake[k].SetMom(i);
                        int pdg = -1000;
                        double charge = pwake[k].GetCharge();
                        int particle_spe = pwake[k].GetId();  // Get stored species
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
    }
}

std::vector<double> WakeGenerator::vec_abs(std::vector<double> p) const {
    for (auto &val : p) val = std::fabs(val);
    return p;
}
