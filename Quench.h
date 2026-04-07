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

    std::array<double,4> _inh_p;

  public:
    Quench();
    Quench(Parton partons);
    Quench(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, vector<double> inh_p);
    virtual ~Quench();

    virtual void display() const;

    void vSetInhP(const std::array<double,4>& p);
    const std::array<double,4>& GetInhP() const;

    void SetRi(double xi, double yi, double zi, double ti);
    void vSetRi(const std::array<double,4>& ri);
    const std::array<double,4>& GetRi() const;

    void SetRf(double xf, double yf, double zf, double tf);
    void vSetRf(const std::array<double,4>& rf);
    const std::array<double,4>& GetRf() const;
};
