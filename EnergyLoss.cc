#include "EnergyLoss.h"
#include <array>
#include <vector>
#include <iostream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <iomanip>
#include "MoliereTables.h"
#include "MoliereElastic.h"
#include "vector_operators.h"

namespace {
constexpr double kLresFinalFlightTime = 10000.;
constexpr double kLresInfiniteTime = 100000000000.;

bool isColored(int id) {
    return std::abs(id) <= 6 || id == 21;
}

double safeFormationTime(const Parton &p) {
    const double q = p.GetQ();
    if (q == 0.) return 0.;
    return 0.2 * 2. * p.vGetP()[3] / (q * q);
}

std::array<double,4> velocity(const std::array<double,4> &p) {
    if (p[3] == 0.) return {0., 0., 0., 1.};
    return {p[0] / p[3], p[1] / p[3], p[2] / p[3], 1.};
}

std::array<double,4> orientationFor(const std::array<double,4> &p) {
    if (p[3] == 0.) return {0., 0., 0., 1.};
    return {p[0] / p[3], p[1] / p[3], p[2] / p[3], 1.};
}

struct LresLifetime {
    double resolve = 0.;
    double creation = 0.;
    double finish = 0.;
    double resolve_abs = 0.;
    double effective_time = 0.;
    std::array<double,4> ri = {0., 0., 0., 0.};
    std::array<double,4> rf = {0., 0., 0., 0.};
};

struct LresState {
    std::array<double,4> p = {0., 0., 0., 0.};
    std::array<double,4> r = {0., 0., 0., 0.};
};
}

EnergyLoss::EnergyLoss(numrand &nr, double kappa, double alpha, int tmethod, int mode,
                       int ebe_hydro, bool do_elastic, bool do_lres, double lres_rpower,
                       bool compat_moliere_legacy_hydro,
                       const std::string &tables_path,
                       const HydroProfile &hydro_profile)
    : nr_(nr), kappa_(kappa), alpha_(alpha), tmethod_(tmethod), mode_(mode),
      ebe_hydro_(ebe_hydro), do_elastic_(do_elastic), do_lres_(do_lres),
      compat_moliere_legacy_hydro_(compat_moliere_legacy_hydro),
      lres_rpower_(lres_rpower), tables_path_(tables_path),
      hydro_profile_(hydro_profile) {
    if (do_elastic_) {
        MoliereTables::ensureLoaded(tables_path_);
    }
}

EnergyLoss::~EnergyLoss() {
    // No special cleanup needed
}

void EnergyLoss::do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched,
                          double x, double y, std::vector<Quench> *recoiled) {
    if (recoiled != nullptr) {
        recoiled->clear();
    }
    // Finite LRES must own the timeline. If elastic is also enabled, Moliere is
    // applied to the currently resolved effective object inside this path.
    if (do_lres_ && lres_rpower_ < 1.e9) {
        do_lres_eloss_impl(partons, quenched, x, y, recoiled);
        return;
    }
    if (do_elastic_) {
        if (recoiled == nullptr) {
            std::vector<Quench> local_recoiled;
            moliere::do_eloss(partons, quenched, x, y, nr_, kappa_, alpha_, tmethod_, mode_,
                              ebe_hydro_, compat_moliere_legacy_hydro_, hydro_profile_, local_recoiled);
        } else {
            moliere::do_eloss(partons, quenched, x, y, nr_, kappa_, alpha_, tmethod_, mode_,
                              ebe_hydro_, compat_moliere_legacy_hydro_, hydro_profile_, *recoiled);
        }
        return;
    }
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
            auto p = quenched[tp].vGetP();
            double q = quenched[tp].GetQ();
            auto pos = quenched[tp].GetRi();
            // Time of flight (from formation time argument)
            double tof = 0.2 * 2. * p[3] / (q * q); // in fm
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
                auto m_p = partons[tp].vGetP();
                auto d1_p = partons[d1].vGetP();
                auto d2_p = partons[d2].vGetP();
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

void EnergyLoss::do_lres_eloss_impl(const std::vector<Parton> &partons, std::vector<Quench> &quenched,
                                    double xcre, double ycre, std::vector<Quench> *recoiled) {
    const size_t n = quenched.size();
    if (n == 0) return;

    std::vector<int> final_ids;
    final_ids.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        if (quenched[i].GetD1() == -1 && quenched[i].GetOrig() != "rem") {
            final_ids.push_back(static_cast<int>(i));
        }
    }

    std::vector<LresLifetime> life(n);
    std::vector<bool> timeline_done(n, false);

    // First build the vacuum formation timeline used by the finite-resolution rule.
    for (int final_id : final_ids) {
        int ind = final_id;
        std::vector<int> family;
        family.push_back(ind);

        double x = xcre;
        double y = ycre;
        double z = 0.;
        double t = 0.;

        while (true) {
            const int mom = quenched[ind].GetMom();
            if (mom >= 0 && mom < static_cast<int>(n)) {
                if (timeline_done[mom]) {
                    x = life[mom].rf[0];
                    y = life[mom].rf[1];
                    z = life[mom].rf[2];
                    t = life[mom].finish;
                    break;
                }
                family.push_back(mom);
                ind = mom;
                continue;
            }
            break;
        }

        for (auto it = family.rbegin(); it != family.rend(); ++it) {
            const int idx = *it;
            const auto p = partons[idx].vGetP();
            const double e = p[3];
            const double tof = (idx == final_id) ? kLresFinalFlightTime : safeFormationTime(partons[idx]);
            const double previous = t;

            life[idx].creation = previous;
            life[idx].ri = {x, y, z, previous};

            if (e != 0.) {
                x += p[0] / e * tof;
                y += p[1] / e * tof;
                z += p[2] / e * tof;
            }
            t = previous + tof;

            life[idx].finish = t;
            life[idx].rf = {x, y, z, t};
            timeline_done[idx] = true;
        }
    }

    std::vector<int> brother(n, -1);
    std::vector<int> effective_mom(n, -1);
    for (size_t i = 0; i < n; ++i) {
        effective_mom[i] = quenched[i].GetMom();
    }

    for (size_t mom = 0; mom < n; ++mom) {
        int first = -1;
        int second = -1;
        for (size_t child = 0; child < n; ++child) {
            if (quenched[child].GetMom() != static_cast<int>(mom)) continue;
            if (first == -1) first = static_cast<int>(child);
            else {
                second = static_cast<int>(child);
                break;
            }
        }
        if (first != -1 && second != -1) {
            brother[first] = second;
            brother[second] = first;
        }
    }

    for (size_t i = 0; i < n; ++i) {
        const int mom = quenched[i].GetMom();
        const int sib = brother[i];
        if (sib < 0 || mom < 0 || mom >= static_cast<int>(n)) continue;

        const auto parent_p = partons[mom].vGetP();
        const auto p_i = partons[i].vGetP();
        const auto p_s = partons[sib].vGetP();
        const auto v_i = velocity(p_i);
        const auto v_s = velocity(p_s);

        const double store = resolution_time(parent_p[3], parent_p[0], parent_p[1], parent_p[2],
                                             life[i].ri[0], life[i].ri[1], life[i].ri[2],
                                             v_s[0] - v_i[0], v_s[1] - v_i[1], v_s[2] - v_i[2],
                                             life[i].creation);
        life[i].resolve = store;
        life[i].resolve_abs = life[i].creation + store;
        life[i].effective_time = life[i].resolve_abs;
    }

    bool changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < n; ++i) {
            const int mom = effective_mom[i];
            if (mom >= 0 && mom < static_cast<int>(n) && effective_mom[mom] >= 0) {
                if (life[i].effective_time < life[mom].effective_time) {
                    life[mom].effective_time = life[i].effective_time;
                    changed = true;
                }
            }

            const int sib = brother[i];
            if (sib >= 0 && sib < static_cast<int>(n) && mom >= 0 && mom < static_cast<int>(n)) {
                const double shared = std::min(life[i].effective_time, life[sib].effective_time);
                if (life[i].effective_time != shared || life[sib].effective_time != shared) {
                    life[i].effective_time = shared;
                    life[sib].effective_time = shared;
                    changed = true;
                }
            }
        }
    }

    changed = true;
    while (changed) {
        changed = false;
        for (size_t i = 0; i < n; ++i) {
            const int mom = effective_mom[i];
            if (mom < 0 || mom >= static_cast<int>(n)) continue;
            const int grand = effective_mom[mom];
            if (grand < 0 || grand >= static_cast<int>(n)) continue;
            if (life[i].effective_time == life[mom].effective_time) {
                effective_mom[i] = grand;
                changed = true;
            }
        }
    }

    std::vector<int> effective_daughter(n, -1);
    for (size_t i = 0; i < n; ++i) {
        const int mom = effective_mom[i];
        if (mom >= 0 && mom < static_cast<int>(n)) {
            effective_daughter[mom] = static_cast<int>(i);
        }
    }

    std::vector<double> live_time(n, 0.);
    for (size_t i = 0; i < n; ++i) {
        const int daughter = effective_daughter[i];
        if (daughter >= 0 && daughter < static_cast<int>(n)) {
            live_time[i] = life[daughter].effective_time - life[i].effective_time;
            if (live_time[i] < 0.) live_time[i] = 0.;
        }
    }

    std::vector<LresState> qstate(n);
    std::vector<bool> done(n, false);
    std::vector<std::array<double,4>> qorient(n);
    std::vector<int> qhad(n, 0);
    for (size_t i = 0; i < n; ++i) {
        qorient[i] = quenched[i].orient();
        qhad[i] = quenched[i].hadScattering();
    }
    std::vector<Quench> lres_moliere_particles;
    std::vector<Quench> local_recoiled;

    for (int final_id : final_ids) {
        int ind = final_id;
        std::vector<int> family;
        family.push_back(ind);
        bool mother_zero = false;

        while (true) {
            const int mom = effective_mom[ind];
            if (mom >= 0 && mom < static_cast<int>(n)) {
                if (done[mom]) {
                    const double parent_e = partons[mom].vGetP()[3];
                    const double frac = (parent_e != 0.) ? qstate[mom].p[3] / parent_e : 0.;
                    qstate[ind].p = partons[ind].vGetP() * frac;
                    if (do_elastic_) {
                        qorient[ind] = orientationFor(qstate[ind].p);
                        if (qhad[mom] == 1 || qhad[mom] == 2) qhad[ind] = 2;
                    }
                    const auto p = partons[ind].vGetP();
                    const double e = p[3];
                    const double dt = life[ind].effective_time - life[ind].creation;
                    qstate[ind].r = life[ind].ri;
                    if (e != 0.) {
                        qstate[ind].r[0] += p[0] / e * dt;
                        qstate[ind].r[1] += p[1] / e * dt;
                        qstate[ind].r[2] += p[2] / e * dt;
                    }
                    qstate[ind].r[3] = life[ind].effective_time;
                    if (qstate[mom].p[3] == 0.) mother_zero = true;
                    break;
                }
                family.push_back(mom);
                ind = mom;
                continue;
            }

            qstate[ind].p = partons[ind].vGetP();
            qstate[ind].r = {xcre, ycre, 0., life[ind].effective_time};
            if (do_elastic_) {
                qorient[ind] = orientationFor(qstate[ind].p);
            }
            break;
        }

        for (auto it = family.rbegin(); it != family.rend(); ++it) {
            const int idx = *it;
            if (mother_zero) {
                for (auto zero_it = it; zero_it != family.rend(); ++zero_it) {
                    qstate[*zero_it].p = {0., 0., 0., 0.};
                    done[*zero_it] = true;
                }
                break;
            }

            auto p = qstate[idx].p;
            const auto vac_p = partons[idx].vGetP();
            const double vac_e = vac_p[3];
            const double dt = life[idx].effective_time - life[idx].creation;
            std::array<double,4> pos = life[idx].ri;
            if (vac_e != 0.) {
                pos[0] += vac_p[0] / vac_e * dt;
                pos[1] += vac_p[1] / vac_e * dt;
                pos[2] += vac_p[2] / vac_e * dt;
            }
            pos[3] = life[idx].effective_time;
            if (pos[3] * pos[3] < pos[2] * pos[2]) {
                pos[3] = std::abs(pos[2]) + 1.e-9;
            }
            double tof = live_time[idx];
            if (idx == final_id) tof = kLresFinalFlightTime;

            double length = 0.;
            double tlength = 0.;
            if (isColored(partons[idx].GetId())) {
                if (do_elastic_) {
                    // LRES has already selected idx as the currently resolved
                    // color object. Moliere scatters that object until its
                    // finite-resolution lifetime ends; daughters enter later.
                    moliere::propagate_segment(p, pos, tof, partons[idx].GetId(), nr_, kappa_, alpha_,
                                               tmethod_, mode_, ebe_hydro_, compat_moliere_legacy_hydro_,
                                               hydro_profile_, lres_moliere_particles, qhad[idx], qorient[idx]);
                } else {
                    loss_rate(p, pos, tof, partons[idx].GetId(), length, tlength);
                }
            } else if (p[3] != 0.) {
                pos += p / p[3] * tof;
            }

            qstate[idx].p = p;
            qstate[idx].r = pos;
            done[idx] = true;
            quenched[idx].AddLength(length, tlength);

            if (p[3] <= 0.) {
                for (auto zero_it = it; zero_it != family.rend(); ++zero_it) {
                    qstate[*zero_it].p = {0., 0., 0., 0.};
                    done[*zero_it] = true;
                }
                break;
            }

            const double frac = (partons[idx].vGetP()[3] != 0.) ? qstate[idx].p[3] / partons[idx].vGetP()[3] : 0.;
            auto next_it = it;
            ++next_it;
            if (next_it != family.rend()) {
                const int daughter = *next_it;
                qstate[daughter].p = partons[daughter].vGetP() * frac;
                if (do_elastic_) {
                    qorient[daughter] = orientationFor(qstate[daughter].p);
                    if (qhad[idx] == 1 || qhad[idx] == 2) qhad[daughter] = 2;
                }
            }
        }
    }

    if (do_elastic_) {
        std::vector<Quench> &recoiled_out = recoiled != nullptr ? *recoiled : local_recoiled;
        moliere::process_recoilers(lres_moliere_particles, nr_, kappa_, alpha_, tmethod_, mode_,
                                   ebe_hydro_, compat_moliere_legacy_hydro_, hydro_profile_, recoiled_out);
    }

    for (size_t i = 0; i < n; ++i) {
        if (!done[i]) continue;
        quenched[i].vSetP(qstate[i].p);
        quenched[i].vSetRi(life[i].ri);
        quenched[i].vSetRf(qstate[i].r);
        quenched[i].vSetInhP(qstate[i].p);
        if (do_elastic_) {
            quenched[i].setOrient(qorient[i]);
            quenched[i].setHadScattering(qhad[i]);
        }
        quenched[i].SetIsDone(true);
    }
}

void EnergyLoss::loss_rate(std::array<double,4> &p, std::array<double,4> &pos, double tof, int id, double &length, double &tlength) {
    double Tc;
    if (tmethod_ == 0) Tc = 0.170;
    else Tc = 0.145;
    constexpr double charm_mass = 1.25;
    constexpr double b_mass = 4.2;

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

    auto w = p / p[3];  // 4-velocity

    do {
#ifdef DO_SOURCE
        // Keep 4momentum before applying quenching this step
        auto p_prev = p;
#endif
        auto p_pre_floor = p;

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
            double temp = 0.;
            hydro_profile_.getValues(tau, pos[0], pos[1], temp, vx, vy);

            double vz = pos[2] / pos[3];
            double frap = atanh(vz);
            vx /= cosh(frap);
            vy /= cosh(frap);

            std::array<double,4> v = {vx, vy, vz, 1.};

            double v2 = v[0]*v[0] + v[1]*v[1] + v[2]*v[2];
            double w2 = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
            double vscalw = v[0] * w[0] + v[1] * w[1] + v[2] * w[2];
            if (v2 >= 1.) v2 = 0.999999999;
            double lore = 1. / sqrt(1. - v2);

            double f_lore = w2 + lore * lore * (v2 - 2. * vscalw + vscalw * vscalw);
            if (f_lore < 0.) {
                f_lore = 0.;
            }
            double f_step = step * sqrt(f_lore);
            f_dist += f_step;

            // temp is already available from getValues
            
            if (temp > Tc) {
                // In-medium distance tracked here (for potential future use)
            }

            // Safe way to exit the plasma: check whether temperature will be above Tc in the next 1000 steps
            if (temp < Tc) {
                // Check whether temperature increases in its way
                for (unsigned int j = 1; j < 1000; j++) {
                    double step_j = step * double(j);
                    double tpos0 = pos[0] + w[0]*step_j;
                    double tpos1 = pos[1] + w[1]*step_j;
                    double tpos2 = pos[2] + w[2]*step_j;
                    double tpos3 = pos[3] + w[3]*step_j;
                    if (tpos3 > tot) break;
                    tau = sqrt(tpos3*tpos3 - tpos2*tpos2);
                    eta = 1. / 2. * log((tpos3 + tpos2) / (tpos3 - tpos2));
                    double ctemp = call_gT(tau, tpos0, tpos1, 0);
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
                    tlength += (temp/0.2)*(temp/0.2) * f_step;
                }

                // Broadening
                if (kappa_ != 0.) {
                    trans_kick(w, w2, v, p, temp, vscalw, lore, step, kappa_);
                }

                p_pre_floor = p;
                bool doquench = true;
                if (std::abs(id) == 4 && p[3] <= charm_mass) doquench = false;
                if (std::abs(id) == 5 && p[3] <= b_mass) doquench = false;

                // Strong coupling
                if (alpha_ != 0. && mode_ == 0 && doquench) {
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

        if (std::abs(id) == 4 && p[3] < charm_mass) {
            p[3] = charm_mass;
            double pmod = std::sqrt(p_pre_floor[0] * p_pre_floor[0] + p_pre_floor[1] * p_pre_floor[1] +
                                    p_pre_floor[2] * p_pre_floor[2]);
            if (pmod == 0.) pmod = 1.;
            p[0] = p_pre_floor[0] / pmod * p[3];
            p[1] = p_pre_floor[1] / pmod * p[3];
            p[2] = p_pre_floor[2] / pmod * p[3];
        }
        if (std::abs(id) == 5 && p[3] < b_mass) {
            p[3] = b_mass;
            double pmod = std::sqrt(p_pre_floor[0] * p_pre_floor[0] + p_pre_floor[1] * p_pre_floor[1] +
                                    p_pre_floor[2] * p_pre_floor[2]);
            if (pmod == 0.) pmod = 1.;
            p[0] = p_pre_floor[0] / pmod * p[3];
            p[1] = p_pre_floor[1] / pmod * p[3];
            p[2] = p_pre_floor[2] / pmod * p[3];
        }

        //This is to check travelled distance, not used
        if (tof < 1000000) {
            double vz = pos[2] / std::max(pos[3], 0.000001);
            double v2 = vz * vz;
            double w2 = w[0]*w[0] + w[1]*w[1] + w[2]*w[2];
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

double EnergyLoss::resolution_time(double parent_e, double parent_px, double parent_py, double parent_pz,
                                   double x, double y, double z,
                                   double dvx, double dvy, double dvz,
                                   double t0) const {
    if (lres_rpower_ < 0.00001) return kLresInfiniteTime;
    if (parent_e == 0.) return 0.;

    const double Tc = (tmethod_ == 0) ? 0.170 : 0.145;
    const double scale = 1. / 3.14 / lres_rpower_;
    const double wx = parent_px / parent_e;
    const double wy = parent_py / parent_e;
    const double wz = parent_pz / parent_e;

    double ti = t0;
    double xprime = 0.;
    double yprime = 0.;
    double zprime = 0.;
    double step_res = 0.1;
    bool refine_step = false;

    for (unsigned int i = 0; i <= 100000000; ++i) {
        const double proper_sq = ti * ti - z * z;
        const double tau = proper_sq > 0. ? std::sqrt(proper_sq) : 0.;

        if (tau >= 0.6) {
            const double temp = call_gT(tau, x, y, 0);
            const double sep = std::sqrt(xprime * xprime + yprime * yprime + zprime * zprime);
            if (sep >= scale / (temp / 0.2) || temp <= Tc) {
                if (step_res == 0.1) {
                    if (i == 0) return ti - t0;
                    refine_step = true;
                } else {
                    return ti - t0;
                }
            }
        }

        if (!refine_step) {
            x += wx * step_res;
            y += wy * step_res;
            z += wz * step_res;
            xprime += dvx * step_res;
            yprime += dvy * step_res;
            zprime += dvz * step_res;
            ti += step_res;
        } else {
            x -= wx * step_res;
            y -= wy * step_res;
            z -= wz * step_res;
            xprime -= dvx * step_res;
            yprime -= dvy * step_res;
            zprime -= dvz * step_res;
            ti -= step_res;
            step_res = 0.01;
            refine_step = false;
        }
    }

    return kLresInfiniteTime;
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

void EnergyLoss::trans_kick(const std::array<double,4> &w, double w2, const std::array<double,4> &v, std::array<double,4> &p, double temp, double vscalw, double lore, double step, double kappa) {
    if (vscalw == 1.) return;

    auto e1 = vec_prod(w, v);
    double Ne1 = normalise(e1);
    if (Ne1 == 0.) {
        double b = 0.5;
        double c = 0.2;
        double a = (-b * w[1] - c * w[2]) / w[0];
        e1 = {a, b, c, 0.};
        Ne1 = normalise(e1);
    }

    double Nw = sqrt(w2);
    auto l = vec_prod(w, e1) / Nw;

    double uscalW = lore * (1. - vscalw);
    double uscall = lore * (-v[0] * l[0] - v[1] * l[1] - v[2] * l[2]);
    double W2 = 1. - w2;

    auto Wp = w - v * (lore * W2 / uscalW);

    double Nalpha = -uscall * uscalW / (uscalW*uscalW - W2);
    if (Nalpha != Nalpha || std::isinf(Nalpha)) return;
    double NN = 1. + W2 * uscall*uscall / (-(uscalW*uscalW) + W2);
    // In some rare situations, this norm squared can be negative. Only do kick otherwise
    if (sqrt(NN) != sqrt(NN)) std::cout << " negative NN " << std::endl;
    else {
        auto e2 = (l + Wp * Nalpha) / sqrt(NN);

        double Ef = p[3] * lore * (1. - vscalw);
        double lore_1mv = lore * (1. - vscalw);
        double wf2 = 1. - W2 / (lore_1mv * lore_1mv);
        double DelQ2 = kappa * temp*temp*temp * lore * (1. - vscalw) * step * 5.;

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

        auto e = e1 * cos(qphi) + e2 * sin(qphi);

        auto Wt = (w - v * (uscalW * lore)) / lore / (1. - vscalw);

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

void EnergyLoss::quenched_sons(const std::array<double,4> &p, const std::array<double,4> &qp, std::array<double,4> &d1, std::array<double,4> &d2) {
    // Mutable local copies for normalisation
    std::array<double,4> p_n = p;
    std::array<double,4> qp_n = qp;
    // Normalise 3-momentum
    double qmod = normalise(qp_n);
    double modmom = normalise(p_n);
    double modthis = normalise(d1);
    double modoson = normalise(d2);
    // Define transverse axis for rotation and normalise
    auto axis = vec_prod(p_n, qp_n);
    normalise(axis);
    // Find angle in plane
    double angle = 0.;
    if (p_n[0]*qp_n[0] + p_n[1]*qp_n[1] + p_n[2]*qp_n[2] >= 1.) angle = 0.;
    else angle = acos(p_n[0]*qp_n[0] + p_n[1]*qp_n[1] + p_n[2]*qp_n[2]);
    // Perform Rodrigues rotation
    double thisscal = axis[0]*d1[0] + axis[1]*d1[1] + axis[2]*d1[2];
    auto use = d1 * cos(angle) + vec_prod(axis, d1) * sin(angle) + axis * (thisscal * (1. - cos(angle)));
    double oscal = axis[0]*d2[0] + axis[1]*d2[1] + axis[2]*d2[2];
    auto ouse = d2 * cos(angle) + vec_prod(axis, d2) * sin(angle) + axis * (oscal * (1. - cos(angle)));
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

double EnergyLoss::normalise(std::array<double,4> &p) {
    double norm = sqrt(p[0]*p[0] + p[1]*p[1] + p[2]*p[2]);
    if (norm == 0.) return norm;
    p[0] /= norm;
    p[1] /= norm;
    p[2] /= norm;
    return norm;
}

std::array<double,4> EnergyLoss::vec_prod(const std::array<double,4> &a, const std::array<double,4> &b) {
    return {a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0], 0.};
}
