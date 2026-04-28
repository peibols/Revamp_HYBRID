#include "Wake.h"
#include <iostream>

#include <vector>
using std::vector;

Wake::Wake()
{

}

Wake::Wake(vector<double> p, double mass, int charge, int id, double status)
{
  _p=p;
  
  _mass=mass;
  _charge=charge;

  _id=id;
  _status=status;

  _mom=-1;
}

Wake::~Wake()
{
  //std::cout << "Wake destructor called" << std::endl;
}

void Wake::display() const
{
  std::cout << " Px = " << _p[0] << " Py= " << _p[1] << " Pz = " << _p[2] << " En= " << _p[3] << std::endl;
  std::cout << " Mass= " << _mass << " Charge= " << _charge << " ";
  std::cout << " id= " << _id << " status= " << _status << std::endl;
}

void Wake::vSetP(vector<double> p)
{
  _p=p;
}
vector<double> Wake::vGetP() const
{
  return _p;
}

void Wake::SetMass(double mass)
{
  _mass=mass;
}
double Wake::GetMass() const
{
  return _mass;
}

void Wake::SetCharge(int charge)
{
  _charge=charge;
}
int Wake::GetCharge() const
{
  return _charge;
}

void Wake::SetId(int id)
{
  _id=id;
}
int Wake::GetId() const
{
  return _id;
}

void Wake::SetMom(int mom)
{
  _mom=mom;
}
int Wake::GetMom() const
{
  return _mom;
}

void Wake::SetStatus(double status)
{
  _status=status;
}
double Wake::GetStatus() const
{
  return _status;
}
