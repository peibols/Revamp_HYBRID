#pragma once

#include <vector>

#include "Parton.h"

using std::vector;

class Quench : public Parton
{
  private:
    vector<double> _ri;
    vector<double> _rf;

    vector<double> _orient;
    
    vector<double> _inh_p;

    int _had_scattering;

    vector<double> _orig_en;

  public:
    Quench();
    Quench(Parton partons);
    Quench(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, vector<double> inh_p);
    ~Quench();

    virtual void display() const;

    void set_orient(vector<double> orient) { _orient = orient; }
    vector<double> orient() const { return _orient; }
    
    void vSetInhP(vector<double> p);
    vector<double> GetInhP() const;

    void SetOrigEn(vector<double> orig_en) { _orig_en=orig_en; }
    vector<double> orig_en() { return _orig_en; }

    double delta_R(Quench other_parton);

    void SetRi(double xi, double yi, double zi, double ti);
    void vSetRi(vector<double> ri);
    vector<double> GetRi() const;

    void SetRf(double xf, double yf, double zf, double tf);
    void vSetRf(vector<double> rf);
    vector<double> GetRf() const;

    void set_had_scattering(int had_scattering) { _had_scattering = had_scattering; }
    int had_scattering() { return _had_scattering; }
};
