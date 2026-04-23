#pragma once

#include <array>
#include <vector>

#include "Parton.h"

using std::vector;

class Quench : public Parton
{
  private:
    std::array<double,4> _ri;
    std::array<double,4> _rf;

    std::array<double,4> _orient;
    std::array<double,4> _inh_p;
    std::array<double,4> _orig_en;

    int _had_scattering;

  public:
    Quench();
    Quench(Parton partons);
    Quench(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, vector<double> inh_p);
    virtual ~Quench();

    virtual void display() const;

    void setOrient(const std::array<double,4>& orient);
    const std::array<double,4>& orient() const;

    void vSetInhP(const std::array<double,4>& p);
    const std::array<double,4>& GetInhP() const;

    void setOrigEn(const std::array<double,4>& orig_en);
    const std::array<double,4>& origEn() const;

    double delta_R(const Quench& other_parton) const;

    void SetRi(double xi, double yi, double zi, double ti);
    void vSetRi(const std::array<double,4>& ri);
    const std::array<double,4>& GetRi() const;

    void SetRf(double xf, double yf, double zf, double tf);
    void vSetRf(const std::array<double,4>& rf);
    const std::array<double,4>& GetRf() const;

    void setHadScattering(int had_scattering);
    int hadScattering() const;
};
