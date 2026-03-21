#pragma once

#include <vector>
#include "Parton.h"
#include "Quench.h"
#include "Random.h"
#include "HydroProfile.h"

class EnergyLoss {
public:
    EnergyLoss(numrand &nr, double kappa, double alpha, int tmethod, int mode, int ebe_hydro, const HydroProfile &hydro_profile);
    ~EnergyLoss();

    // Perform energy loss on the given partons
    void do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched, double x, double y);

private:
    numrand &nr_;
    double kappa_;
    double alpha_;
    int tmethod_;
    int mode_;
    int ebe_hydro_;
    const HydroProfile &hydro_profile_;

    // Private member functions for energy loss calculations
    void do_eloss_impl(const std::vector<Parton> &partons, std::vector<Quench> &quenched, double x, double y);
    void loss_rate(std::vector<double> &p, std::vector<double> &pos, double tof, int id, double &length, double &tlength);
    void get_source_evol(double &tau_ev, double& x_f, double& y_f, double& vx_f, double& vy_f, double tau_ini, double x_ini, double y_ini, double Tc);
    double call_gT(double tau, double x, double y, int comp) const;
    void quenched_sons(const std::vector<double> &p, const std::vector<double> &qp, std::vector<double> &d1, std::vector<double> &d2);
    double normalise(std::vector<double> &p);
    std::vector<double> vec_prod(const std::vector<double> &a, const std::vector<double> &b);
    void trans_kick(const std::vector<double> &w, double w2, const std::vector<double> &v, std::vector<double> &p, double temp, double vscalw, double lore, double step, double kappa);
};