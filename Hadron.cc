#include "Hadron.h"
#include <iostream>
#include <vector>

using std::vector;

Hadron::Hadron()
{
}

Hadron::Hadron(Parton partons) : Parton(partons)
{
	for (int i = 0; i < 4; i++) _ri.push_back(0.), _rf.push_back(0.);
	
	_charge=0.;
	_width=0.;
	
	_isdone=false;
}

//Preferred constructor for now
Hadron::Hadron(Parton partons, double charge, double width) : Parton(partons)
{
	for (int i = 0; i < 4; i++) _ri.push_back(0.), _rf.push_back(0.);

        _charge=charge;
        _width=width;
}

Hadron::Hadron(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, double charge, double width) : Parton(partons)
{
	_ri.push_back(xi);
	_ri.push_back(yi);
	_ri.push_back(zi);
	_ri.push_back(ti);

	_rf.push_back(xf);
        _rf.push_back(yf);
        _rf.push_back(zf);
        _rf.push_back(tf);
	
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

void Hadron::vSetRi(vector<double> ri)
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
vector<double> Hadron::GetRi() const
{
	return _ri;
}

void Hadron::vSetRf(vector<double> rf)
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
vector<double> Hadron::GetRf() const
{
        return _rf;
}
