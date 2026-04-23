#include "Quench.h"
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>

using std::vector;

Quench::Quench()
{
}

Quench::Quench(Parton partons) : Parton(partons)
{
  _ri.fill(0.);
  _rf.fill(0.);
  _inh_p.fill(0.);
  _orig_en.fill(0.);
  _orient = {0., 0., 0., 1.};
  if (_p[3] != 0.) {
    _orient = {_p[0] / _p[3], _p[1] / _p[3], _p[2] / _p[3], 1.};
  }
  _isdone=false;
  _had_scattering = 0;
}

//Do I actually use this constructor?
Quench::Quench(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, vector<double> inh_p) : Parton(partons)
{
  _ri = {xi, yi, zi, ti};
  _rf = {xf, yf, zf, tf};
  _inh_p[0]=inh_p[0]; _inh_p[1]=inh_p[1]; _inh_p[2]=inh_p[2]; _inh_p[3]=inh_p[3];
  _orig_en.fill(0.);
  _orient = {0., 0., 0., 1.};
  if (_p[3] != 0.) {
    _orient = {_p[0] / _p[3], _p[1] / _p[3], _p[2] / _p[3], 1.};
  }
  _isdone=false;
  _had_scattering = 0;
}

Quench::~Quench()
{
  //std::cout << "Quench destructor called" << std::endl;
}

void Quench::display() const
{
  Parton::display();
  std::cout << " inh_px= " << _inh_p[0] << " inh_py= " << _inh_p[1] << " inh_pz= " << _inh_p[2] << " inh_en= " << _inh_p[3] << std::endl;
  std::cout << " xi= " << _ri[0] << " yi= " << _ri[1] << " zi= " << _ri[2] << " ti= " << _ri[3] << std::endl;
  std::cout << " xf= " << _rf[0] << " yf= " << _rf[1] << " zf= " << _rf[2] << " tf= " << _rf[3] << std::endl;
}

void Quench::setOrient(const std::array<double,4>& orient)
{
  _orient = orient;
}

const std::array<double,4>& Quench::orient() const
{
  return _orient;
}

void Quench::vSetInhP(const std::array<double,4>& p)
{
  _inh_p=p;
}
const std::array<double,4>& Quench::GetInhP() const
{
  return _inh_p; 
}

void Quench::setOrigEn(const std::array<double,4>& orig_en)
{
  _orig_en = orig_en;
}

const std::array<double,4>& Quench::origEn() const
{
  return _orig_en;
}

double Quench::delta_R(const Quench& other_parton) const
{
  double delta_eta = GetEta() - other_parton.GetEta();
  double this_pt = GetPt();
  double other_pt = other_parton.GetPt();
  double delta_phi = 0.;
  if (this_pt > 0. && other_pt > 0.) {
    double cos_dphi = (_p[0] * other_parton.vGetP()[0] + _p[1] * other_parton.vGetP()[1]) / this_pt / other_pt;
    cos_dphi = std::max(-1.0, std::min(1.0, cos_dphi));
    delta_phi = std::acos(cos_dphi);
  }
  return std::sqrt(delta_phi * delta_phi + delta_eta * delta_eta);
}

void Quench::vSetRi(const std::array<double,4>& ri)
{
  _ri=ri;
}
void Quench::SetRi(double xi, double yi, double zi, double ti)
{
  _ri[0]=xi;
  _ri[1]=yi;
  _ri[2]=zi;
  _ri[3]=ti;
}
const std::array<double,4>& Quench::GetRi() const
{
  return _ri;
}

void Quench::vSetRf(const std::array<double,4>& rf)
{
  _rf=rf;
}
void Quench::SetRf(double xf, double yf, double zf, double tf) 
{
  _rf[0]=xf;
  _rf[1]=yf;
  _rf[2]=zf;
  _rf[3]=tf;
}
const std::array<double,4>& Quench::GetRf() const
{
  return _rf;
}

void Quench::setHadScattering(int had_scattering)
{
  _had_scattering = had_scattering;
}

int Quench::hadScattering() const
{
  return _had_scattering;
}
