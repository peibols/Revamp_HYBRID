#include "EnergyLoss.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iomanip>
#include "vector_operators.h"

// Forward declarations for global functions that remain
void quenched_sons(std::vector<double> p, std::vector<double> qp, std::vector<double> &d1, std::vector<double> &d2);

EnergyLoss::EnergyLoss(numrand &nr, double kappa, double alpha, int tmethod, int mode, int ebe_hydro, const HydroProfile &hydro_profile)
    : nr_(nr), kappa_(kappa), alpha_(alpha), tmethod_(tmethod), mode_(mode), ebe_hydro_(ebe_hydro), hydro_profile_(hydro_profile) {
}

EnergyLoss::~EnergyLoss() {
    // No special cleanup needed
}

void EnergyLoss::do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched, double x, double y) {
    do_eloss_impl(partons, quenched, x, y);
}

void EnergyLoss::do_eloss_impl(const std::vector<Parton> &partons, std::vector<Quench> &quenched, double xcre, double ycre) {
    // Tag final particles
    std::vector<int> FinId;
    for (size_t i = 0; i < quenched.size(); i++) {
        // Exclude remnants
        if (quenched[i].GetD1() == -1 && quenched[i].GetOrig() != "rem") {
            FinId.push_back(i);
        }
    }

    // Energy Loss Loop: select a final particle, find its oldest undone parent, climb down the family chain. Iterate until all final particles are done
    for (size_t i = 0; i < FinId.size(); i++) {
        int ind = FinId[i];  // Start loop with final particle
        std::vector<int> Fam;  // Family chain array
        Fam.push_back(ind);
        double inhe = 1.;    // Used to see whether mother was completely quenched
        // Family chain loop
        int found = 0;
        do {
            int mom = quenched[ind].GetMom();
            // If "ind" is not a parent parton
            if (mom != -1) {
                // If mother is done
                if (quenched[mom].GetIsDone() == true) {
                    inhe = quenched[mom].vGetP()[3];
                    // End of family chain
                    found = 1;
                }
                // If it is not done
                else {
                    Fam.push_back(mom);
                    ind = mom;
                }
            }
            // If "ind" is a parent parton
            else {
                quenched[ind].SetRi(xcre, ycre, 0., 0.);
                // End of family chain
                found = 1;
            }
        } while (found == 0);

        // Apply energy loss chronologically
        for (size_t w = Fam.size(); w > 0; w--) {
            int tp = Fam[w - 1];
            // If first done mother was totally quenched, set all descendance quenched and done, and exit loop
            if (inhe == 0.) {
                for (size_t j = w; j > 0; j--) {
                    tp = Fam[j - 1];
                    quenched[tp].SetP(0., 0., 0., 0.);
                    quenched[tp].SetIsDone(true);
                }
                break;
            }
            std::vector<double> p = quenched[tp].vGetP();
            double q = quenched[tp].GetQ();
            std::vector<double> pos = quenched[tp].GetRi();
            // Time of flight (from formation time argument)
            double tof = 0.2 * 2. * p[3] / pow(q, 2.); // in fm
            // If final particle, fly arbitrarily far
            if (w == 1) tof = 10000000000.;
            double length = 0.; // length in QGP
            double tlength = 0.; // temperature weighted length in QGP
            // If colored particle
            if (abs(quenched[tp].GetId()) <= 6 || quenched[tp].GetId() == 21) {
                loss_rate(p, pos, tof, quenched[tp].GetId(), length, tlength);
            } else {
                // If not colored particle, don't do energy loss, but propagate position and time manually
                pos += p / p[3] * tof;
            }
            // Update mother momenta and set positions to final
            quenched[tp].vSetP(p);
            quenched[tp].vSetRf(pos);
            quenched[tp].SetIsDone(true);
            quenched[tp].AddLength(length, tlength);
            // If it got fully quenched, quenched descendance and exit
            if (p[3] == 0.) {
                for (size_t j = w; j > 0; j--) {
                    tp = Fam[j - 1];
                    quenched[tp].SetP(0., 0., 0., 0.);
                    quenched[tp].SetIsDone(true);
                }
                break;
            }
            // If not final particle, propagate quenching to son in chain, and store results for other son
            if (w != 1) {
                // Find two daughters
                int d1 = Fam[w - 2];
                int d2;
                if (quenched[tp].GetD1() == d1) d2 = quenched[tp].GetD2();
                else d2 = quenched[tp].GetD1();
                // Find new momenta for sons: rotate and quench tri-momentum, quench energy
                std::vector<double> m_p = partons[tp].vGetP();
                std::vector<double> d1_p = partons[d1].vGetP();
                std::vector<double> d2_p = partons[d2].vGetP();
                quenched_sons(m_p, p, d1_p, d2_p);
                // Propagate momenta and positions
                quenched[d1].vSetP(d1_p);
                quenched[d1].vSetInhP(d1_p);
                quenched[d1].vSetRi(pos);
                quenched[d1].AddLength(length, tlength);
                quenched[d2].vSetP(d2_p);
                quenched[d2].vSetInhP(d2_p);
                quenched[d2].vSetRi(pos);
                quenched[d2].AddLength(length, tlength);
            }
        }
        Fam.clear();
    }
    FinId.clear();
}

void EnergyLoss::loss_rate(std::vector<double> &p, std::vector<double> &pos, double tof, int id, double &length, double &tlength) {
    double Tc;
    if (tmethod_ == 0) Tc = 0.170;
    else Tc = 0.145;

    double tot = pos[3] + tof;    // Final time

    double tau0h = 0.6;             //Ave hydro
    if (ebe_hydro_ == 1) tau0h = 0.4;   //ebe hydro

    double ei = p[3];      // Initial energy

    double f_dist = 0.;    // Traversed distance in Fluid Frame

    double virt_f_dist = 0.;

    double CF;
    if (id == 21) {
        if (mode_ == 0) CF = pow(9. / 4., 1. / 3.);  // If gluon, color charge dependence is ratio of casimirs to power 1/3
        else CF = 9. / 4.;
    } else CF = 1.;

    int marker = 0;    // If one, exit loop
    double step = 0.1;  // Time step in LAB frame

    std::vector<double> w = p / p[3];  // 4-velocity

    do {
        // Keep 4momentum before applying quenching this step
        std::vector<double> p_prev = p;

        if (pos[3] == tot) marker = 1;
        if (pos[3] > tot) std::cout << " Warning: Went beyond tot= " << tot << " t= " << pos[3] << std::endl;

        // Proper time
        double tau = sqrt(pos[3] * pos[3] - pos[2] * pos[2]);
        if (tau != tau) {
            std::cout << " TAU Not a number z= " << pos[2] << " t= " << pos[3] << " wz= " << w[2] << " en = " << p[3] << " pz= " << p[2] << "\n";
            std::cout << " Id= " << id << std::endl;
            exit(1);
        }

        // Rapidity
        double eta = 1. / 2. * log((pos[3] + pos[2]) / (pos[3] - pos[2]));
        if (eta != eta && tau > 0.) {
            std::cout << " Eta is NaN= " << eta << " t= " << pos[3] << " z= " << pos[2] << std::endl;
            exit(1);
        }

        int will_hot = 0;  // Advance variable (to reach hot zones)
        double vx = 0.;
        double vy = 0.;
        if (tau >= tau0h) {  // Hydro profile starting time
            std::vector<double> v;
            vx = call_gT(tau, pos[0], pos[1], 1);
            vy = call_gT(tau, pos[0], pos[1], 2);

            double vz = pos[2] / pos[3];
            double frap = atanh(vz);
            vx /= cosh(frap);
            vy /= cosh(frap);

            v.push_back(vx);
            v.push_back(vy);
            v.push_back(vz);
            v.push_back(1.);

            double v2 = pow(v[0], 2.) + pow(v[1], 2.) + pow(v[2], 2.);
            double w2 = pow(w[0], 2.) + pow(w[1], 2.) + pow(w[2], 2.);
            double vscalw = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
            if (v2 >= 1.) v2 = 0.999999999;
            double lore = 1. / sqrt(1. - v2);

            double f_lore = w2 + lore * lore * (v2 - 2. * vscalw + vscalw * vscalw);
            if (f_lore < 0.) {
                f_lore = 0.;
            }
            double f_step = step * sqrt(f_lore);
            f_dist += f_step;

            double temp = call_gT(tau, pos[0], pos[1], 0);
            
            if (temp > Tc) {
                // In-medium distance tracked here (for potential future use)
            }

            // Safe way to exit the plasma: check whether temperature will be above Tc in the next 1000 steps
            if (temp < Tc) {
                // Check whether temperature increases in its way
                for (unsigned int j = 1; j < 1000; j++) {
                    std::vector<double> tpos = pos + w * step * double(j);
                    if (tpos[3] > tot) break;
                    tau = sqrt(tpos[3] * tpos[3] - tpos[2] * tpos[2]);
                    eta = 1. / 2. * log((tpos[3] + tpos[2]) / (tpos[3] - tpos[2]));
                    double ctemp = call_gT(tau, tpos[0], tpos[1], 0);
                    if (ctemp > Tc) { //It will get to hot
                        will_hot = int(j);
                        break;
                    }
                }
                if (will_hot == 0) { // It will not get to hot
                    pos += w * (tot - pos[3]);
                    marker = 1;
                }
            }

            // Now broad&quench
            if (p[3] > 0. && temp >= Tc && f_step != 0.) {
                if (mode_ == 0) {
                    length += f_step;
                    tlength += temp / 0.2 * f_step;
                } else {
                    length += f_step;
                    tlength += pow(temp / 0.2, 2.) * f_step;
                }

                // Broadening
                if (kappa_ != 0.) {
                    trans_kick(w, w2, v, p, temp, vscalw, lore, step, kappa_);
                }

                // Strong coupling
                if (alpha_ != 0. && mode_ == 0) {
                    double Efs = ei * lore * (1. - vscalw);
                    double tstop = 0.2 * pow(Efs, 1. / 3.) / (2. * pow(temp, 4. / 3.) * alpha_) / CF;
                    double beta = tstop / f_dist;
                    if (beta > 1.) {
                        double intpiece = Efs * step * 4. / (3.141592) * (1. / (beta * tstop * sqrt(beta * beta - 1.)));
                        double quench = (p[3] - intpiece) / p[3];
                        p *= quench;
                    } else {
                        p[3] = 0.;
                    }
                }

                // Radiative
                if (alpha_ != 0. && mode_ == 1) {
                    double intpiece = CF * (step / 0.2) * alpha_ * temp * temp * temp * (f_dist / 0.2);
                    double quench = (p[3] - intpiece) / p[3];
                    p *= quench;
                }

                // Collisional
                if (alpha_ != 0. && mode_ == 2) {
                    double intpiece = CF * (step / 0.2) * alpha_ * temp * temp;
                    double quench = (p[3] - intpiece) / p[3];
                    p *= quench;
                }
                
            }
        }

        //This is to check travelled distance, not used
        if (tof < 1000000) {
            double vz = pos[2] / std::max(pos[3], 0.000001);
            double v2 = vz * vz;
            double w2 = pow(w[0], 2.) + pow(w[1], 2.) + pow(w[2], 2.);
            double vscalw = vz * w[2];
            if (v2 >= 1.) v2 = 0.999999999;
            double lore = 1. / sqrt(1. - v2);
            double f_lore = w2 + lore * lore * (v2 - 2. * vscalw + vscalw * vscalw);
            if (f_lore < 0.) {
                std::cout << " craazy f_lore= " << f_lore << std::endl;
                f_lore = 1.;
            }
            virt_f_dist += step * sqrt(f_lore);
        }

        // If parton gets totally quenched, exit
        if (p[3] <= 0.) {
            marker = 1;
            for (unsigned int i = 0; i < 4; i++) p[i] = 0.;
        } else {
            // Manually protect very soft particles from getting kicks that yield velocities greater than 1
            for (unsigned int i = 0; i < 3; i++) {
                if (p[i] > p[3]) {
                    std::cout << " Got crazy kick in i= " << i << "p[i]= " << p[i] << " and p[3]= " << p[3] << std::endl;
                    p[i] = 0.99999 * p[3];
                }
            }
            // Update kinematical quantities, with the possibility of advancing to hot regions
            w = p / p[3];
            double tstep = std::max(double(will_hot), 1.) * step;
            if (pos[3] + tstep > tot) {
                tstep = tot - pos[3];
            }
            if (marker != 1) pos += w * tstep;
        }

#ifdef DO_SOURCE
        // Fill source file, for new wake purposes
        if (p[3] != p_prev[3]) {
            // Get tau ev, x_f, y_f and vx_f and vy_f for source file
            double tau_ev, x_f, y_f, vx_f, vy_f;
            get_source_evol(tau_ev, x_f, y_f, vx_f, vy_f, tau, pos[0], pos[1], Tc);

            std::ofstream source_file("SOURCE.dat", std::ios_base::app);
            source_file << tau << " " << pos[0] << " " << pos[1] << " " << eta << " " << tau_ev << " "
                        << -p[3] + p_prev[3] << " " << -p[0] + p_prev[0] << " " << -p[1] + p_prev[1] << " " << -p[2] + p_prev[2] << " "
                        << vx << " " << vy << " " << vx_f << " " << vy_f << " " << x_f << " " << y_f << std::endl;
        }
#endif

    } while (marker == 0);

    //Just to check the distance travelled, not used
    if (virt_f_dist == 0.) virt_f_dist = -1;
}

double EnergyLoss::call_gT(double tau, double x, double y, int comp) const {
    if (tau < 0.) return 0.0;

    switch (comp) {
        case 0:
            return hydro_profile_.temperature(tau, x, y);
        case 1:
            return hydro_profile_.velocityX(tau, x, y);
        case 2:
            return hydro_profile_.velocityY(tau, x, y);
        default:
            std::cerr << "Wrong comp= " << comp << " in call_gT" << std::endl;
            return 0.0;
    }
}

void EnergyLoss::get_source_evol(double &tau_ev, double& x_f, double& y_f, double& vx_f, double& vy_f, double tau_ini, double x_ini, double y_ini, double Tc) {
    tau_ev = 0.;
    double dtau = 0.1;
    x_f = x_ini;
    y_f = y_ini;
    double tau_now = tau_ini;
    double vx_f_local = 0., vy_f_local = 0.;

    while (true) {
        vx_f_local = call_gT(tau_now, x_f, y_f, 1);
        vy_f_local = call_gT(tau_now, x_f, y_f, 2);

        double localT;
        localT = call_gT(tau_now, x_f, y_f, 0);
        if (localT < Tc) break;

        x_f += vx_f_local * dtau;
        y_f += vy_f_local * dtau;
        tau_ev += dtau;
        tau_now += dtau;
    }

    vx_f = vx_f_local;
    vy_f = vy_f_local;
}

void EnergyLoss::trans_kick(const std::vector<double> &w, double w2, const std::vector<double> &v, std::vector<double> &p, double temp, double vscalw, double lore, double step, double kappa) {
    if (vscalw == 1.) return;

    std::vector<double> e1 = vec_prod(w, v);
    double Ne1 = normalise(e1);
    if (Ne1 == 0.) {
        double b = 0.5;
        double c = 0.2;
        double a = (-b * w[1] - c * w[2]) / w[0];
        e1 = {a, b, c, 0.};
        Ne1 = normalise(e1);
    }

    double Nw = sqrt(w2);
    std::vector<double> l = vec_prod(w, e1) / Nw;

    double uscalW = lore * (1. - vscalw);
    double uscall = lore * (-v[0] * l[0] - v[1] * l[1] - v[2] * l[2]);
    double W2 = 1. - w2;

    std::vector<double> Wp = w - v * lore * W2 / uscalW;

    double Nalpha = -uscall * uscalW / (pow(uscalW, 2.) - W2);
    if (Nalpha != Nalpha || std::isinf(Nalpha)) return;
    double NN = 1. + W2 * pow(uscall, 2.) / (-pow(uscalW, 2.) + W2);
    // In some rare situations, this norm squared can be negative. Only do kick otherwise
    if (sqrt(NN) != sqrt(NN)) std::cout << " negative NN " << std::endl;
    else {
        std::vector<double> e2 = (l + Wp * Nalpha) / sqrt(NN);

        double Ef = p[3] * lore * (1. - vscalw);
        double wf2 = 1. - W2 / pow(lore * (1. - vscalw), 2.);
        double DelQ2 = kappa * pow(temp, 3.) * lore * (1. - vscalw) * step * 5.;

        double qfac = 0.;
        // Only do kick if energy is greater than temperature
        if (Ef * sqrt(wf2) > 0.) {
            qfac = sqrt(-1. * log(nr_.rando())) * sqrt(DelQ2); // Box-Muller method
            if (qfac > Ef * sqrt(wf2)) {
                qfac = Ef * sqrt(wf2) - 0.00000001;
            }
        }

        double qbeta;
        if (wf2 > 0.) qbeta = sqrt(1. - qfac * qfac / Ef / Ef / wf2) - 1.;
        else qbeta = 0.;
        if (qbeta != qbeta) {
            std::cout << " qbeta= " << std::setprecision(6) << qbeta << " qfac= " << qfac << " Ef= " << Ef << " wf2= " << wf2 << " cutoff= " << Ef * sqrt(wf2) << std::endl;
            std::cout << " W2= " << std::setprecision(6) << W2 << " lore= " << lore << " vscalw= " << vscalw << std::endl;
        }

        double qphi = 2. * 3.141592654 * nr_.rando();

        std::vector<double> e = e1 * cos(qphi) + e2 * sin(qphi);

        std::vector<double> Wt = (w - v * uscalW * lore) / lore / (1. - vscalw);

        // Update 4momentum
        if (lore == 1.) {
            e2 = vec_prod(e1, w);
            e = e1 * cos(qphi) + e2 * sin(qphi);
        }
        p += Wt * qbeta * Ef + e * qfac;
        if (p[3] != p[3]) {
            std::cout << " p in bro= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << std::endl;
            std::cout << " qbeta= " << qbeta << " qfac= " << qfac << std::endl;
            std::cout << " Wt= " << Wt[0] << " " << Wt[1] << " " << Wt[2] << " " << Wt[3] << std::endl;
            std::cout << " e1= " << e1[0] << " " << e1[1] << " " << e1[2] << " " << e1[3] << std::endl;
            std::cout << " e2= " << e2[0] << " " << e2[1] << " " << e2[2] << " " << e2[3] << std::endl;
            std::cout << " Nalpha= " << Nalpha << std::endl;
            std::cout << " l= " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << std::endl;
            std::cout << " Wp= " << Wp[0] << " " << Wp[1] << " " << Wp[2] << " " << Wp[3] << std::endl;
            std::cout << " uscalW= " << uscalW << std::endl;
            std::cout << " Nw= " << Nw << " uscall= " << uscall << std::endl;
            std::cout << " W2= " << W2 << std::endl;
            exit(1);
        }
    }
}

void EnergyLoss::quenched_sons(const std::vector<double> &p, const std::vector<double> &qp, std::vector<double> &d1, std::vector<double> &d2) {
    // Normalise 3-momentum
    double qmod = normalise(const_cast<std::vector<double>&>(qp));
    double modmom = normalise(const_cast<std::vector<double>&>(p));
    double modthis = normalise(d1);
    double modoson = normalise(d2);
    // Define transverse axis for rotation and normalise
    std::vector<double> axis = vec_prod(p, qp);
    normalise(axis);
    // Find angle in plane
    double angle = 0.;
    if (p[0] * qp[0] + p[1] * qp[1] + p[2] * qp[2] >= 1.) angle = 0.;
    else angle = acos(p[0] * qp[0] + p[1] * qp[1] + p[2] * qp[2]);
    // Perform Rodrigues rotation
    double thisscal = axis[0] * d1[0] + axis[1] * d1[1] + axis[2] * d1[2];
    std::vector<double> use = d1 * cos(angle) + vec_prod(axis, d1) * sin(angle) + axis * thisscal * (1. - cos(angle));
    double oscal = axis[0] * d2[0] + axis[1] * d2[1] + axis[2] * d2[2];
    std::vector<double> ouse = d2 * cos(angle) + vec_prod(axis, d2) * sin(angle) + axis * oscal * (1. - cos(angle));
    // Update momenta
    double lamp = qmod / modmom;
    double lambda = qp[3] / p[3];
    for (unsigned int i = 0; i < 3; i++) {
        d1[i] = use[i] * modthis * lamp;
        d2[i] = ouse[i] * modoson * lamp;
    }
    d1[3] *= lambda;
    d2[3] *= lambda;
}

double EnergyLoss::normalise(std::vector<double> &p) {
    double norm = sqrt(pow(p[0], 2.) + pow(p[1], 2.) + pow(p[2], 2.));
    if (norm == 0.) return norm;
    p[0] /= norm;
    p[1] /= norm;
    p[2] /= norm;
    return norm;
}

std::vector<double> EnergyLoss::vec_prod(const std::vector<double> &a, const std::vector<double> &b) {
    std::vector<double> rot(4, 0.);
    rot[0] = a[1] * b[2] - a[2] * b[1];
    rot[1] = a[2] * b[0] - a[0] * b[2];
    rot[2] = a[0] * b[1] - a[1] * b[0];
    rot[3] = 0.;
    return rot;
}