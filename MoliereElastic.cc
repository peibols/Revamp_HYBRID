#include "MoliereElastic.h"

#include <cstdlib>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

#include "FourVector.h"
#include "gsl/gsl_integration.h"
#include "vector_operators.h"

using std::vector;
using namespace std;

#include "Distributions.hpp"

namespace moliere {
namespace {

struct GslIntegrationWorkspaces {
    GslIntegrationWorkspaces()
        : wdk(gsl_integration_workspace_alloc(100000)),
          wkcm(gsl_integration_workspace_alloc(100000)),
          wx(gsl_integration_workspace_alloc(100000)) {}

    ~GslIntegrationWorkspaces() {
        gsl_integration_workspace_free(wdk);
        gsl_integration_workspace_free(wkcm);
        gsl_integration_workspace_free(wx);
    }

    GslIntegrationWorkspaces(const GslIntegrationWorkspaces&) = delete;
    GslIntegrationWorkspaces& operator=(const GslIntegrationWorkspaces&) = delete;

    gsl_integration_workspace *wdk;
    gsl_integration_workspace *wkcm;
    gsl_integration_workspace *wx;
};

GslIntegrationWorkspaces &integration_workspaces() {
    static thread_local GslIntegrationWorkspaces workspaces;
    return workspaces;
}

vector<double> to_vec(const std::array<double,4> &a) {
    return {a[0], a[1], a[2], a[3]};
}

std::array<double,4> to_arr(const vector<double> &v) {
    return {v[0], v[1], v[2], v[3]};
}

double transverse_momentum_transfer(const vector<double> &p_before, const vector<double> &p_after) {
    const double qx = p_after[0] - p_before[0];
    const double qy = p_after[1] - p_before[1];
    const double qz = p_after[2] - p_before[2];
    const double pmag2 = p_before[0] * p_before[0] +
                         p_before[1] * p_before[1] +
                         p_before[2] * p_before[2];
    const double qmag2 = qx * qx + qy * qy + qz * qz;
    if (pmag2 <= 0.) return std::sqrt(qmag2);

    const double qdotp = qx * p_before[0] + qy * p_before[1] + qz * p_before[2];
    const double qparallel2 = qdotp * qdotp / pmag2;
    return std::sqrt(std::max(0., qmag2 - qparallel2));
}

double gT(const HydroProfile &hydro_profile, double tau, double x, double y) {
    return hydro_profile.temperature(tau, x, y);
}

double gVx(const HydroProfile &hydro_profile, double tau, double x, double y) {
    return hydro_profile.velocityX(tau, x, y);
}

double gVy(const HydroProfile &hydro_profile, double tau, double x, double y) {
    return hydro_profile.velocityY(tau, x, y);
}

double gT_legacy_elastic(const HydroProfile &hydro_profile, double tau, double x, double y, double eta) {
    return hydro_profile.temperatureElasticLegacy(tau, x, y, eta);
}

double gVx_legacy_elastic(const HydroProfile &hydro_profile, double tau, double x, double y, double eta) {
    return hydro_profile.velocityXElasticLegacy(tau, x, y, eta);
}

double gVy_legacy_elastic(const HydroProfile &hydro_profile, double tau, double x, double y, double eta) {
    return hydro_profile.velocityYElasticLegacy(tau, x, y, eta);
}

double normalise(vector<double> &p) {
    double norm = std::sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    if (norm == 0.) return norm;
    p[0] /= norm;
    p[1] /= norm;
    p[2] /= norm;
    return norm;
}

vector<double> vec_prod(vector<double> a, vector<double> b) {
    return {a[1]*b[2] - b[1]*a[2],
            -a[0]*b[2] + b[0]*a[2],
            a[0]*b[1] - b[0]*a[1],
            0.};
}

void quenched_sons(vector<double> p, vector<double> qp, vector<double> &d1, vector<double> &d2) {
    double qmod = normalise(qp);
    double modmom = normalise(p);
    double modthis = normalise(d1);
    double modoson = normalise(d2);
    vector<double> axis = vec_prod(p, qp);
    normalise(axis);
    double angle = 0.;
    if (p[0]*qp[0] + p[1]*qp[1] + p[2]*qp[2] >= 1.) angle = 0.;
    else angle = std::acos(p[0]*qp[0] + p[1]*qp[1] + p[2]*qp[2]);
    double thisscal = axis[0]*d1[0] + axis[1]*d1[1] + axis[2]*d1[2];
    vector<double> use = d1*std::cos(angle) + vec_prod(axis,d1)*std::sin(angle) + axis*thisscal*(1.-std::cos(angle));
    double oscal = axis[0]*d2[0] + axis[1]*d2[1] + axis[2]*d2[2];
    vector<double> ouse = d2*std::cos(angle) + vec_prod(axis,d2)*std::sin(angle) + axis*oscal*(1.-std::cos(angle));
    double lamp = qmod/modmom;
    double lambda = qp[3]/p[3];
    for (unsigned int i=0; i<3; ++i) {
        d1[i] = use[i]*modthis*lamp;
        d2[i] = ouse[i]*modoson*lamp;
    }
    d1[3] *= lambda;
    d2[3] *= lambda;
}

void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp,
                double vscalw, double lore, double step, double kappa, numrand &nr) {
    if (vscalw == 1.) return;

    vector<double> e1 = vec_prod(w, v);
    double Ne1 = normalise(e1);
    if (Ne1 == 0.) {
        double b = 0.5;
        double c = 0.2;
        double a = (-b*w[1] - c*w[2]) / w[0];
        e1 = {a, b, c, 0.};
        Ne1 = normalise(e1);
    }

    double Nw = std::sqrt(w2);
    vector<double> l = vec_prod(w, e1) / Nw;

    double uscalW = lore*(1.-vscalw);
    double uscall = lore*(-v[0]*l[0] - v[1]*l[1] - v[2]*l[2]);
    double W2 = 1.-w2;

    vector<double> Wp = w - v*lore*W2/uscalW;

    double Nalpha = -uscall*uscalW/(std::pow(uscalW,2.) - W2);
    if (Nalpha != Nalpha || std::isinf(Nalpha)) return;
    double NN = 1. + W2*std::pow(uscall,2.)/(-std::pow(uscalW,2.) + W2);
    if (std::sqrt(NN) != std::sqrt(NN)) {
        std::cout << " negative NN " << std::endl;
        return;
    }

    vector<double> e2 = (l + Wp*Nalpha)/std::sqrt(NN);

    double Ef = p[3]*lore*(1.-vscalw);
    double wf2 = 1. - W2/std::pow(lore*(1.-vscalw),2.);
    double DelQ2 = kappa*std::pow(temp,3.)*lore*(1.-vscalw)*step*5.;

    double qfac = 0.;
    if (Ef*std::sqrt(wf2) > 0.) {
        qfac = std::sqrt(-1.*std::log(nr.rando()))*std::sqrt(DelQ2);
        if (qfac > Ef*std::sqrt(wf2)) {
            qfac = Ef*std::sqrt(wf2) - 0.00000001;
        }
    }

    double qbeta;
    if (wf2 > 0.) qbeta = std::sqrt(1. - qfac*qfac/Ef/Ef/wf2) - 1.;
    else qbeta = 0.;
    if (qbeta != qbeta) {
        std::cout << " qbeta= " << std::setprecision(6) << qbeta
                  << " qfac= " << qfac << " Ef= " << Ef << " wf2= " << wf2
                  << " cutoff= " << Ef*std::sqrt(wf2) << std::endl;
        std::cout << " W2= " << std::setprecision(6) << W2
                  << " lore= " << lore << " vscalw= " << vscalw << std::endl;
    }

    double qphi = 2.*3.141592654*nr.rando();
    vector<double> e = e1*std::cos(qphi) + e2*std::sin(qphi);
    vector<double> Wt = (w - v*uscalW*lore)/lore/(1.-vscalw);

    if (lore == 1.) {
        e2 = vec_prod(e1, w);
        normalise(e2);
        e = e1*std::cos(qphi) + e2*std::sin(qphi);
    }
    p += Wt*qbeta*Ef + e*qfac;
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
        std::exit(1);
    }
}

FourVector Boost(double b[3], FourVector p) {
    double betamod = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    double gamma = 1.0/std::sqrt(1.0-betamod*betamod);
    double pscalb = p.x()*b[0] + p.y()*b[1] + p.z()*b[2];
    double spat = (pscalb*gamma/(1.0+gamma) - p.t())*gamma;
    double ptt = gamma*(p.t() - pscalb);
    double ptx = p.x() + b[0]*spat;
    double pty = p.y() + b[1]*spat;
    double ptz = p.z() + b[2]*spat;
    return FourVector(ptx, pty, ptz, ptt);
}

FourVector BoostBack(double b[3], FourVector p) {
    double betamod = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
    double gamma = 1.0/std::sqrt(1.0-betamod*betamod);
    double pscalb = p.x()*b[0] + p.y()*b[1] + p.z()*b[2];
    double spat = (-pscalb*gamma/(1.0+gamma) - p.t())*gamma;
    double bt = gamma*(p.t() + pscalb);
    double bx = p.x() - b[0]*spat;
    double by = p.y() - b[1]*spat;
    double bz = p.z() - b[2]*spat;
    return FourVector(bx, by, bz, bt);
}

void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa,
               double alpha, int tmethod, int model, int ebe_hydro,
               bool compat_moliere_legacy_hydro, const HydroProfile &hydro_profile, vector<Quench> &new_particles,
               int &had_scattering, vector<double> &orient,
               const ScatteringCallback *scattering_callback = nullptr,
               const PropagationStepCallback *step_callback = nullptr) {
    auto &workspaces = integration_workspaces();
    gsl_integration_workspace *wdk = workspaces.wdk;
    gsl_integration_workspace *wkcm = workspaces.wkcm;
    gsl_integration_workspace *wx = workspaces.wx;

    double Tc = (tmethod == 0) ? 0.170 : 0.145;
    double charm_mass = 1.25;
    double b_mass = 4.2;

    double tot = pos[3] + tof;
    double ei = p[3];
    double f_dist = 0.;
    double l_dist = 0.;

    double CF;
    if (id == 21) {
        if (model == 0) CF = std::pow(9./4.,1./3.);
        else CF = 9./4.;
    } else CF = 1.;

    int marker = 0;
    double step = 0.1;

    vector<double> p_prev;
    vector<double> p_ini = p;
    vector<double> w = p/p[3];
    vector<double> o_in = orient;

    do {
        const vector<double> pos_step_start = pos;
        const vector<double> p_step_start = p;
        double temp_for_record = 0.;
        double tau_for_record = 0.;
        int in_medium_for_record = 0;
        p_prev = p;

        if (pos[3] == tot) marker = 1;
        if (pos[3] > tot) std::cout << " Warning: Went beyond tot= " << tot << " t= " << pos[3] << std::endl;

        double tau = std::sqrt(pos[3]*pos[3]-pos[2]*pos[2]);
        tau_for_record = tau;
        if (tau != tau) {
            std::cout << " TAU Not a number z= " << pos[2] << " t= " << pos[3]
                      << " wz= " << w[2] << " en = " << p[3] << " pz= " << p[2] << "\n";
            std::cout << " Id= " << id << std::endl;
            std::exit(1);
        }

        double eta = 0.;
        if (tau > 0.) eta = 0.5*std::log((pos[3]+pos[2])/(pos[3]-pos[2]));
        if (eta != eta && tau > 0.) eta = 0.;

        int will_hot = 0;
        double vx = 0., vy = 0., vz = 0.;
        double tau0h = (ebe_hydro == 1) ? 0.4 : 0.6;
        if (tau >= tau0h || hydro_profile.hasPreHydroAt(tau)) {
            vector<double> v;
            vx = compat_moliere_legacy_hydro
                ? gVx_legacy_elastic(hydro_profile, tau, pos[0], pos[1], eta)
                : gVx(hydro_profile, tau, pos[0], pos[1]);
            vy = compat_moliere_legacy_hydro
                ? gVy_legacy_elastic(hydro_profile, tau, pos[0], pos[1], eta)
                : gVy(hydro_profile, tau, pos[0], pos[1]);
            vz = pos[2]/pos[3];
            double frap = std::atanh(vz);
            vx /= std::cosh(frap);
            vy /= std::cosh(frap);
            v.push_back(vx); v.push_back(vy); v.push_back(vz); v.push_back(1.);

            double v2 = std::pow(v[0],2.) + std::pow(v[1],2.) + std::pow(v[2],2.);
            double w2 = std::pow(w[0],2.) + std::pow(w[1],2.) + std::pow(w[2],2.);
            double vscalw = v[0]*w[0] + v[1]*w[1] + v[2]*w[2];
            if (v2 >= 1.) v2 = 0.999999999;
            double lore = 1./std::sqrt(1.-v2);

            double temp = compat_moliere_legacy_hydro
                ? gT_legacy_elastic(hydro_profile, tau, pos[0], pos[1], eta)
                : gT(hydro_profile, tau, pos[0], pos[1]);
            temp_for_record = temp;
            in_medium_for_record = temp >= Tc ? 1 : 0;

            l_dist += step;
            double f_lore = w2 + lore*lore*(v2 - 2.*vscalw + vscalw*vscalw);
            if (f_lore < 0.) f_lore = 0.;
            double f_step = step*std::sqrt(f_lore);
            f_dist += f_step;

            if (temp < Tc) {
                for (unsigned int j=1; j<1000; ++j) {
                    vector<double> tpos = pos + w*step*double(j);
                    if (tpos[3] > tot) break;
                    tau = std::sqrt(tpos[3]*tpos[3]-tpos[2]*tpos[2]);
                    eta = 0.5*std::log((tpos[3]+tpos[2])/(tpos[3]-tpos[2]));
                    const double ttemp = compat_moliere_legacy_hydro
                        ? gT_legacy_elastic(hydro_profile, tau, tpos[0], tpos[1], eta)
                        : gT(hydro_profile, tau, tpos[0], tpos[1]);
                    if (ttemp > Tc) {
                        will_hot = int(j);
                        break;
                    }
                }
                if (will_hot == 0) {
                    pos += w*(tot-pos[3]);
                    marker = 1;
                }
            }

            if (p[3] > 0. && temp >= Tc) {
                FourVector pp;
                pp.Set(p[0], p[1], p[2], std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));
                double beta[3] = {v[0], v[1], v[2]};
                pp = Boost(beta, pp);

                double pin = std::sqrt(pp.x()*pp.x()+pp.y()*pp.y()+pp.z()*pp.z())/temp;
                if ((id == 21 || std::abs(id) <= 4) && pin < 1500.) {
                    double deltime = f_step*temp/0.2;
                    vector<double> elscat = gen_particles(pp.x()/temp, pp.y()/temp, pp.z()/temp, id, deltime, wdk, wkcm, wx);
                    if (elscat[0] == 1) {
                        const vector<double> p_before_scattering = p;
                        FourVector pf(elscat[2]*temp, elscat[3]*temp, elscat[4]*temp,
                                      temp*std::sqrt(elscat[2]*elscat[2]+elscat[3]*elscat[3]+elscat[4]*elscat[4]));
                        pf = BoostBack(beta, pf);
                        int pf_id = int(elscat[5]);

                        FourVector kf(elscat[6]*temp, elscat[7]*temp, elscat[8]*temp,
                                      temp*std::sqrt(elscat[6]*elscat[6]+elscat[7]*elscat[7]+elscat[8]*elscat[8]));
                        kf = BoostBack(beta, kf);
                        int kf_id = int(elscat[9]);

                        FourVector ff(elscat[10]*temp, elscat[11]*temp, elscat[12]*temp,
                                      temp*std::sqrt(elscat[10]*elscat[10]+elscat[11]*elscat[11]+elscat[12]*elscat[12]));
                        ff = BoostBack(beta, ff);
                        int ff_id = int(elscat[13]);

                        pp = BoostBack(beta, pp);
                        vector<double> vpos{pos[0], pos[1], pos[2], pos[3]};
                        double edif = 100000000000.;
                        int iclose = -1000;
                        if (std::abs(id) != 4) {
                            if (std::fabs(pp.t()-pf.t()) < edif) { edif = std::fabs(pp.t()-pf.t()); iclose = 0; }
                            if (std::fabs(pp.t()-ff.t()) < edif) { edif = std::fabs(pp.t()-ff.t()); iclose = 1; }
                        } else {
                            if (pf_id == id && ff_id == id) {
                                std::cout << " Both charms!?" << std::endl;
                                std::exit(1);
                            }
                            if (pf_id == id) iclose = 0;
                            else if (ff_id == id) iclose = 1;
                            else {
                                std::cout << " No match for charm!?" << std::endl;
                                std::exit(1);
                            }
                        }

                        vector<double> p_after_scattering;
                        vector<double> recoiler_p;
                        int recoiler_id = 0;
                        if (iclose == 0) {
                            p_after_scattering = {pf.x(), pf.y(), pf.z(), pf.t()};
                            recoiler_p = {ff.x(), ff.y(), ff.z(), ff.t()};
                            recoiler_id = ff_id;
                        } else if (iclose == 1) {
                            p_after_scattering = {ff.x(), ff.y(), ff.z(), ff.t()};
                            recoiler_p = {pf.x(), pf.y(), pf.z(), pf.t()};
                            recoiler_id = pf_id;
                        } else {
                            std::cout << "No match for normal particle" << std::endl;
                            std::exit(1);
                        }

                        vector<double> vkf{kf.x(), kf.y(), kf.z(), kf.t()};
                        ScatteringDecision scattering_decision = ScatteringDecision::Apply;
                        if (scattering_callback != nullptr) {
                            ScatteringCandidate candidate;
                            candidate.pos = to_arr(vpos);
                            candidate.p_before = to_arr(p_before_scattering);
                            candidate.p_after = to_arr(p_after_scattering);
                            candidate.recoiler_p = to_arr(recoiler_p);
                            candidate.hole_p = to_arr(vkf);
                            candidate.qperp = transverse_momentum_transfer(
                                p_before_scattering, p_after_scattering);
                            candidate.recoiler_id = recoiler_id;
                            candidate.hole_id = kf_id;
                            scattering_decision = (*scattering_callback)(candidate);
                            if (scattering_decision == ScatteringDecision::StopBeforeApply) {
                                marker = 1;
                                break;
                            }
                        }

                        if (scattering_decision == ScatteringDecision::Apply) {
                            had_scattering = 1;
                            p = p_after_scattering;
                            new_particles.emplace_back(Parton(recoiler_p, 100000000., 0., 0, -1, -1,
                                                              recoiler_id, "recoiler", 0, 0, false));
                            new_particles.back().vSetRi(to_arr(vpos));
                            new_particles.emplace_back(Parton(vkf, 100000000., 0., 0, -1, -1, kf_id, "hole", 0, 0, true));
                            new_particles.back().vSetRi(to_arr(vpos));
                            p_prev = p;
                        }
                    } else if (elscat[0] == -1) {
                        std::cout << " Elscat problem!" << std::endl;
                    }
                }

                if (marker == 1) continue;

                if (kappa != 0. && step != 0.) trans_kick(w, w2, v, p, temp, vscalw, lore, step, kappa, nr);

                orient[0] = p[0]/p[3];
                orient[1] = p[1]/p[3];
                orient[2] = p[2]/p[3];

                bool doquench = true;
                if (std::abs(id) == 4 && p[3] <= charm_mass) doquench = false;
                if (std::abs(id) == 5 && p[3] <= b_mass) doquench = false;

                p_prev = p;

                if (alpha != 0. && model == 0 && doquench) {
                    double Efs = ei*lore*(1.-vscalw);
                    double tstop = 0.2*std::pow(Efs,1./3.)/(2.*std::pow(temp,4./3.)*alpha)/CF;
                    double beta_s = tstop/f_dist;
                    if (beta_s > 1.) {
                        double intpiece = Efs*step*4./(3.141592)*(1./(beta_s*tstop*std::sqrt(beta_s*beta_s-1.)));
                        double quench = (p[3]-intpiece)/p[3];
                        p *= quench;
                    } else {
                        p[3] = 0.;
                    }
                }
                if (alpha != 0. && model == 1) {
                    double intpiece = CF*(step/0.2)*alpha*temp*temp*temp*(f_dist/0.2);
                    double quench = (p[3]-intpiece)/p[3];
                    p *= quench;
                }
                if (alpha != 0. && model == 2) {
                    double intpiece = CF*(step/0.2)*alpha*temp*temp;
                    double quench = (p[3]-intpiece)/p[3];
                    p *= quench;
                }
            }
        }

        if (std::abs(id) == 4 && p[3] < charm_mass) {
            p[3] = charm_mass;
            double pmod = std::sqrt(p_prev[0]*p_prev[0]+p_prev[1]*p_prev[1]+p_prev[2]*p_prev[2]);
            if (pmod == 0.) pmod = 1.;
            p[0] = p_prev[0]/pmod*p[3];
            p[1] = p_prev[1]/pmod*p[3];
            p[2] = p_prev[2]/pmod*p[3];
        }
        if (std::abs(id) == 5 && p[3] < b_mass) {
            p[3] = b_mass;
            double pmod = std::sqrt(p_prev[0]*p_prev[0]+p_prev[1]*p_prev[1]+p_prev[2]*p_prev[2]);
            if (pmod == 0.) pmod = 1.;
            p[0] = p_prev[0]/pmod*p[3];
            p[1] = p_prev[1]/pmod*p[3];
            p[2] = p_prev[2]/pmod*p[3];
        }

        if (p[3] <= 0.) {
            marker = 1;
            for (unsigned int i=0; i<4; ++i) p[i] = 0.;
        } else {
            for (unsigned int i=0; i<3; ++i) {
                if (p[i] > p[3]) {
                    std::cout << " Got crazy kick in i= " << i << "p[i]= " << p[i] << " and p[3]= " << p[3] << std::endl;
                    p[i] = 0.99999*p[3];
                }
            }
            w = p/p[3];
            double tstep = std::max(double(will_hot),1.)*step;
            if (pos[3] + tstep > tot) tstep = tot-pos[3];
            if (marker != 1) pos += w*tstep;
        }

        if (step_callback != nullptr) {
            PropagationStep step_record;
            step_record.pos_before = to_arr(pos_step_start);
            step_record.pos_after = to_arr(pos);
            step_record.p_before = to_arr(p_step_start);
            step_record.p_after = to_arr(p);
            step_record.temperature = temp_for_record;
            step_record.tau = tau_for_record;
            step_record.step = pos[3] - pos_step_start[3];
            step_record.in_medium = in_medium_for_record;
            (*step_callback)(step_record);
        }
    } while (marker == 0);

    double scalpprev = orient[0]*o_in[0] + orient[1]*o_in[1] + orient[2]*o_in[2];
    if (scalpprev > 1.) scalpprev = 1.;
    (void)std::acos(scalpprev);
}

}  // namespace

std::default_random_engine elastic_generator_state() {
    return ::generator;
}

void set_elastic_generator_state(const std::default_random_engine &state) {
    ::generator = state;
}

void seed_elastic_generator(unsigned int seed) {
    ::generator.seed(seed);
}

void propagate_segment(std::array<double,4> &p, std::array<double,4> &pos, double tof, int id,
                       numrand &nr, double kappa, double alpha, int tmethod, int model, int ebe_hydro,
                       bool compat_moliere_legacy_hydro, const HydroProfile &hydro_profile,
                       std::vector<Quench> &new_particles, int &had_scattering,
                       std::array<double,4> &orient,
                       const PropagationStepCallback &step_callback) {
    vector<double> p_vec = to_vec(p);
    vector<double> pos_vec = to_vec(pos);
    vector<double> orient_vec = to_vec(orient);
    loss_rate(p_vec, pos_vec, tof, id, nr, kappa, alpha, tmethod, model,
              ebe_hydro, compat_moliere_legacy_hydro, hydro_profile,
              new_particles, had_scattering, orient_vec, nullptr,
              step_callback ? &step_callback : nullptr);
    p = to_arr(p_vec);
    pos = to_arr(pos_vec);
    orient = to_arr(orient_vec);
}

void propagate_segment_with_scattering_callback(std::array<double,4> &p, std::array<double,4> &pos,
                                                double tof, int id,
                                                numrand &nr, double kappa, double alpha,
                                                int tmethod, int model, int ebe_hydro,
                                                bool compat_moliere_legacy_hydro,
                                                const HydroProfile &hydro_profile,
                                                std::vector<Quench> &new_particles,
                                                int &had_scattering,
                                                std::array<double,4> &orient,
                                                const ScatteringCallback &callback,
                                                const PropagationStepCallback &step_callback) {
    vector<double> p_vec = to_vec(p);
    vector<double> pos_vec = to_vec(pos);
    vector<double> orient_vec = to_vec(orient);
    loss_rate(p_vec, pos_vec, tof, id, nr, kappa, alpha, tmethod, model,
              ebe_hydro, compat_moliere_legacy_hydro, hydro_profile,
              new_particles, had_scattering, orient_vec, &callback,
              step_callback ? &step_callback : nullptr);
    p = to_arr(p_vec);
    pos = to_arr(pos_vec);
    orient = to_arr(orient_vec);
}

void process_recoilers(std::vector<Quench> &new_particles, numrand &nr, double kappa, double alpha,
                       int tmethod, int model, int ebe_hydro, bool compat_moliere_legacy_hydro,
                       const HydroProfile &hydro_profile, std::vector<Quench> &recoiled) {
    while (true) {
        vector<Quench> current_particles = new_particles;
        new_particles.clear();
        for (unsigned int ip = 0; ip < current_particles.size(); ++ip) {
            if (current_particles[ip].GetOrig() == "hole") {
                current_particles[ip].SetIsDone(true);
                recoiled.push_back(current_particles[ip]);
                continue;
            }
            vector<double> p = to_vec(current_particles[ip].vGetP());
            vector<double> pos = to_vec(current_particles[ip].GetRi());
            double tof = current_particles[ip].GetQ();
            vector<double> orient = to_vec(current_particles[ip].orient());
            std::array<double,4> orig_en = current_particles[ip].vGetP();
            int had_scattering = 0;
            loss_rate(p, pos, tof, current_particles[ip].GetId(), nr, kappa, alpha, tmethod, model,
                      ebe_hydro, compat_moliere_legacy_hydro, hydro_profile, new_particles, had_scattering, orient);
            current_particles[ip].setOrigEn(orig_en);
            current_particles[ip].vSetP(p);
            current_particles[ip].vSetRf(to_arr(pos));
            current_particles[ip].setOrient(to_arr(orient));
            current_particles[ip].setHadScattering(had_scattering);
            current_particles[ip].SetIsDone(true);
            if (p[3] == 0.) current_particles[ip].SetP(0.,0.,0.,0.);
            recoiled.push_back(current_particles[ip]);
        }
        if (new_particles.empty()) break;
        std::cout << " RESCATTERING! " << std::endl;
    }
}

void do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched, double xcre, double ycre,
              numrand &nr, double kappa, double alpha, int tmethod, int model, int ebe_hydro,
              bool compat_moliere_legacy_hydro, const HydroProfile &hydro_profile, std::vector<Quench> &recoiled,
              const PartonCallbackFactory &callback_factory) {
    vector<int> FinId;
    for (unsigned int i = 0; i < quenched.size(); ++i) {
        if (quenched[i].GetD1() == -1 && quenched[i].GetOrig() != "rem") {
            FinId.push_back(i);
        }
    }

    vector<Quench> new_particles;
    for (unsigned int i = 0; i < FinId.size(); ++i) {
        int ind = FinId[i];
        vector<int> Fam;
        Fam.push_back(ind);
        double inhe = 1.;
        int found = 0;
        do {
            int mom = quenched[ind].GetMom();
            if (mom != -1) {
                if (quenched[mom].GetIsDone() == true) {
                    inhe = quenched[mom].vGetP()[3];
                    found = 1;
                } else {
                    Fam.push_back(mom);
                    ind = mom;
                }
            } else {
                quenched[ind].SetRi(xcre, ycre, 0., 0.);
                found = 1;
            }
        } while (found == 0);

        for (unsigned int w = Fam.size(); w > 0; --w) {
            int tp = Fam[w-1];
            if (inhe == 0.) {
                for (unsigned int j = w; j > 0; --j) {
                    tp = Fam[j-1];
                    quenched[tp].SetP(0.,0.,0.,0.);
                    quenched[tp].SetIsDone(true);
                }
                break;
            }
            vector<double> p = to_vec(quenched[tp].vGetP());
            vector<double> pos = to_vec(quenched[tp].GetRi());
            double tof = quenched[tp].GetQ();
            double q = quenched[tp].GetQ();
            if (q > 0.) tof = 0.2 * 2. * p[3] / (q * q);
            if (w == 1) tof = 10000000000.;

            int had_scattering = quenched[tp].hadScattering();
            vector<double> orient = to_vec(quenched[tp].orient());
            if (std::abs(quenched[tp].GetId()) <= 6 || quenched[tp].GetId() == 21) {
                ScatteringCallback scattering_callback;
                PropagationStepCallback step_callback;
                if (callback_factory) {
                    auto callbacks = callback_factory(tp, quenched[tp].GetId(), quenched[tp].GetMom(),
                                                      quenched[tp].GetD1(), quenched[tp].GetD2());
                    scattering_callback = callbacks.first;
                    step_callback = callbacks.second;
                }
                loss_rate(p, pos, tof, quenched[tp].GetId(), nr, kappa, alpha, tmethod, model,
                          ebe_hydro, compat_moliere_legacy_hydro, hydro_profile, new_particles, had_scattering, orient,
                          scattering_callback ? &scattering_callback : nullptr,
                          step_callback ? &step_callback : nullptr);
            } else {
                pos += p/p[3]*tof;
            }

            quenched[tp].vSetP(p);
            quenched[tp].vSetRf(to_arr(pos));
            quenched[tp].setOrient(to_arr(orient));
            quenched[tp].setHadScattering(had_scattering);
            quenched[tp].SetIsDone(true);

            if (p[3] == 0.) {
                for (unsigned int j = w; j > 0; --j) {
                    tp = Fam[j-1];
                    quenched[tp].SetP(0.,0.,0.,0.);
                    quenched[tp].SetIsDone(true);
                }
                break;
            }
            if (w != 1) {
                int d1 = Fam[w-2];
                int d2;
                if (quenched[tp].GetD1() == d1) d2 = quenched[tp].GetD2();
                else d2 = quenched[tp].GetD1();
                vector<double> m_p = to_vec(partons[tp].vGetP());
                vector<double> d1_p = to_vec(partons[d1].vGetP());
                vector<double> d2_p = to_vec(partons[d2].vGetP());
                quenched_sons(m_p, p, d1_p, d2_p);

                vector<double> orient_d1{d1_p[0]/d1_p[3], d1_p[1]/d1_p[3], d1_p[2]/d1_p[3], d1_p[3]/d1_p[3]};
                quenched[d1].setOrient(to_arr(orient_d1));
                quenched[d1].vSetP(d1_p);
                quenched[d1].vSetInhP(to_arr(d1_p));
                quenched[d1].vSetRi(to_arr(pos));

                vector<double> orient_d2{d2_p[0]/d2_p[3], d2_p[1]/d2_p[3], d2_p[2]/d2_p[3], d2_p[3]/d2_p[3]};
                quenched[d2].setOrient(to_arr(orient_d2));
                quenched[d2].vSetP(d2_p);
                quenched[d2].vSetInhP(to_arr(d2_p));
                quenched[d2].vSetRi(to_arr(pos));
                if (had_scattering == 1 || had_scattering == 2) {
                    quenched[d1].setHadScattering(2);
                    quenched[d2].setHadScattering(2);
                }
            }
        }
        Fam.clear();
    }
    FinId.clear();

    process_recoilers(new_particles, nr, kappa, alpha, tmethod, model, ebe_hydro,
                      compat_moliere_legacy_hydro, hydro_profile, recoiled);
}

}  // namespace moliere
