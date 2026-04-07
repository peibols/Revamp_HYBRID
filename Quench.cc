#include "Quench.h"
#include <iostream>
#include <vector>

using std::vector;

Quench::Quench()
{
}

Quench::Quench(Parton partons) : Parton(partons)
{
  _ri.fill(0.); _rf.fill(0.); _inh_p.fill(0.);
  _isdone=false;
}

//Do I actually use this constructor?
Quench::Quench(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, vector<double> inh_p) : Parton(partons)
{
  _ri = {xi, yi, zi, ti};
  _rf = {xf, yf, zf, tf};
  _inh_p[0]=inh_p[0]; _inh_p[1]=inh_p[1]; _inh_p[2]=inh_p[2]; _inh_p[3]=inh_p[3];
  _isdone=false;
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

void Quench::vSetInhP(const std::array<double,4>& p)
{
  _inh_p=p;
}
const std::array<double,4>& Quench::GetInhP() const
{
  return _inh_p; 
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
