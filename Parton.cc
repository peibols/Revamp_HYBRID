#include "Parton.h"
#include <iostream>
#include <vector>
#include <cmath>

using std::vector;

Parton::Parton()
{

}

Parton::Parton(vector<double> p, double q, double mass, int mom, int d1, int d2, int id, std::string orig, int col, int acol, bool isdone)
{
  _p=p;

  _q=q;
  _mass=mass;
  
  _mom=mom;
  _d1=d1;
  _d2=d2;

  _id=id;
  _orig=orig;

  _col=col;
  _acol=acol;

  _isdone=isdone;

  _length=0.;
  _tlength=0.;

  for (unsigned in=0; in<4; in++) _ri.push_back(-1.);
}

Parton::~Parton()
{
  //std::cout << "Parton destructor called" << std::endl;
}

void Parton::display() const
{
  std::cout << " Px = " << _p[0] << " Py= " << _p[1] << " Pz = " << _p[2] << " En= " << _p[3] << std::endl;
  std::cout << " mom= " << _mom << " d1= " << _d1 << " d2= " << _d2 << " id= " << _id << " orig= " << _orig << " IsDone?= " << _isdone << std::endl;
}

void Parton::vSetP(vector<double> p)
{
  _p=p;
}
void Parton::SetP(double px, double py, double pz, double en)
{
  _p[0]=px;
  _p[1]=py;
  _p[2]=pz;
  _p[3]=en;
}
vector<double> Parton::vGetP() const
{
  return _p;
}

double Parton::GetPt() const
{
  double pt=sqrt(pow(_p[0],2.)+pow(_p[1],2.));
  return pt;
}

void Parton::vSetRi(vector<double> ri)
{
  _ri=ri;
}
vector<double> Parton::GetRi() const
{
  return _ri;
}

double Parton::GetEta() const
{
  double pt=sqrt(pow(_p[0],2.)+pow(_p[1],2.));
  double eta=1./2.*log((sqrt(pow(pt,2.)+pow(_p[2],2.))+_p[2])/(sqrt(pow(pt,2.)+pow(_p[2],2.))-_p[2]));
  return eta;
}

void Parton::SetQ(double q)
{
  _q=q;
}
double Parton::GetQ() const
{
  return _q;
}

void Parton::SetMass(double mass)
{
  _mass=mass;
}
double Parton::GetMass() const
{
  return _mass;
}

void Parton::SetMom(int mom)
{
  _mom=mom;
}
int Parton::GetMom() const
{
  return _mom;
}

void Parton::SetD1(int d1)
{
  _d1=d1;
}
int Parton::GetD1() const
{
  return _d1;
}

void Parton::SetD2(int d2)
{
  _d2=d2;
}
int Parton::GetD2() const
{
  return _d2;
}

void Parton::SetId(int id)
{
  _id=id;
}
int Parton::GetId() const
{
  return _id;
}

void Parton::SetOrig(std::string orig)
{
  _orig=orig;
}
std::string Parton::GetOrig() const
{
  return _orig;
}

void Parton::SetCol(int col)
{
  _col=col;
}
int Parton::GetCol() const
{
  return _col;
}

void Parton::SetAcol(int acol)
{
  _acol=acol;
}
int Parton::GetAcol() const
{
  return _acol;
}

void Parton::SetIsDone(bool isdone)
{
  _isdone=isdone;
}
bool Parton::GetIsDone() const
{
  return _isdone;
}
