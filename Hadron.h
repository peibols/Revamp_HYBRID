#pragma once

#include <vector>

#include "Parton.h"

using std::vector;

class Hadron : public Parton
{
	private:
		//Keep position information in case one wants to do hadron rescattering
		vector<double> _ri;
		vector<double> _rf;

		double _charge;
		double _width;

	public:
		Hadron();
		Hadron(Parton partons);
		Hadron(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, double charge, double width);
		Hadron(Parton partons, double charge, double width);
		~Hadron();

		virtual void display() const;

		void SetRi(double xi, double yi, double zi, double ti);
		void vSetRi(vector<double> ri);
		vector<double> GetRi() const;

		void SetRf(double xf, double yf, double zf, double tf);
               	void vSetRf(vector<double> rf);
		vector<double> GetRf() const;

		void SetCharge(double charge);
                double GetCharge() const;

		void SetWidth(double width);
                double GetWidth() const;
};
