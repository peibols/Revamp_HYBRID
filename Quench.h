#pragma once

#include <vector>

#include "Parton.h"

using std::vector;

class Quench : public Parton
{
  private:
    vector<double> _ri;
    vector<double> _rf;

    vector<double> _inh_p;

  public:
    Quench();
    Quench(Parton partons);
    Quench(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, vector<double> inh_p);
    ~Quench();

    virtual void display() const;

    void vSetInhP(vector<double> p);
    vector<double> GetInhP() const;

    void SetRi(double xi, double yi, double zi, double ti);
    void vSetRi(vector<double> ri);
    vector<double> GetRi() const;

    void SetRf(double xf, double yf, double zf, double tf);
    void vSetRf(vector<double> rf);
    vector<double> GetRf() const;
};
