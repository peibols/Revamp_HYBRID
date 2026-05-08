#pragma once

#include <array>
#include <string>
#include <vector>
#include "Parton.h"
#include "Quench.h"
#include "Random.h"
#include "HydroProfile.h"

class EnergyLoss {
public:
    EnergyLoss(numrand &nr, double kappa, double alpha, int tmethod, int mode,
               int ebe_hydro, bool do_elastic, const std::string &tables_path,
               const HydroProfile &hydro_profile);
    ~EnergyLoss();

    // Perform energy loss on the given partons
    void do_eloss(const std::vector<Parton> &partons, std::vector<Quench> &quenched,
                  double x, double y, std::vector<Quench> *recoiled = nullptr);

private:
    numrand &nr_;
    double kappa_;
    double alpha_;
    int tmethod_;
    int mode_;
    int ebe_hydro_;
    bool do_elastic_;
    std::string tables_path_;
    const HydroProfile &hydro_profile_;

    // Private member functions for energy loss calculations
    void do_eloss_impl(const std::vector<Parton> &partons, std::vector<Quench> &quenched, double x, double y);
    void loss_rate(std::array<double,4> &p, std::array<double,4> &pos, double tof, int id, double &length, double &tlength);
    void get_source_evol(double &tau_ev, double& x_f, double& y_f, double& vx_f, double& vy_f, double tau_ini, double x_ini, double y_ini, double Tc);
    double call_gT(double tau, double x, double y, int comp) const;
    void quenched_sons(const std::array<double,4> &p, const std::array<double,4> &qp, std::array<double,4> &d1, std::array<double,4> &d2);
    double normalise(std::array<double,4> &p);
    std::array<double,4> vec_prod(const std::array<double,4> &a, const std::array<double,4> &b);
    void trans_kick(const std::array<double,4> &w, double w2, const std::array<double,4> &v, std::array<double,4> &p, double temp, double vscalw, double lore, double step, double kappa);
};
