#include "Hadron.h"
#include <iostream>
#include <vector>

using std::vector;

Hadron::Hadron()
{
}

Hadron::Hadron(Parton partons) : Parton(partons)
{
  _ri.fill(0.); _rf.fill(0.);
  
  _charge=0.;
  _width=0.;
  
  _isdone=false;
}

//Preferred constructor for now
Hadron::Hadron(Parton partons, double charge, double width) : Parton(partons)
{
  _ri.fill(0.); _rf.fill(0.);

  _charge=charge;
  _width=width;
}

Hadron::Hadron(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, double charge, double width) : Parton(partons)
{
  _ri = {xi, yi, zi, ti};
  _rf = {xf, yf, zf, tf};
  
  _isdone=false;

  _charge=charge;
  _width=width;
}

Hadron::~Hadron()
{
  //std::cout << "Hadron destructor called" << std::endl;
}

void Hadron::display() const
{
  Parton::display();
  //std::cout << " xi= " << _ri[0] << " yi= " << _ri[1] << " zi= " << _ri[2] << " ti= " << _ri[3] << std::endl;
  //std::cout << " xf= " << _rf[0] << " yf= " << _rf[1] << " zf= " << _rf[2] << " tf= " << _rf[3] << std::endl;
  std::cout << " charge= " << _charge << " width= " << _width << std::endl;
}

void Hadron::SetCharge(double charge)
{
  _charge=charge;
}
double Hadron::GetCharge() const
{
  return _charge;
}

void Hadron::SetWidth(double width)
{
  _width=width;
}
double Hadron::GetWidth() const
{
  return _width;
}

void Hadron::vSetRi(const std::array<double,4>& ri)
{
  _ri=ri;
}
void Hadron::SetRi(double xi, double yi, double zi, double ti)
{
  _ri[0]=xi;
  _ri[1]=yi;
  _ri[2]=zi;
  _ri[3]=ti;
}
const std::array<double,4>& Hadron::GetRi() const
{
  return _ri;
}

void Hadron::vSetRf(const std::array<double,4>& rf)
{
  _rf=rf;
}
void Hadron::SetRf(double xf, double yf, double zf, double tf) 
{
  _rf[0]=xf;
  _rf[1]=yf;
  _rf[2]=zf;
  _rf[3]=tf;
}
const std::array<double,4>& Hadron::GetRf() const
{
  return _rf;
}
