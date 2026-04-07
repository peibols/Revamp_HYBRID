#pragma once

#include <array>
#include <vector>

#include "Parton.h"

using std::vector;

class Hadron : public Parton
{
  private:
    //Keep position information in case one wants to do hadron rescattering
    std::array<double,4> _ri;
    std::array<double,4> _rf;

    double _charge;
    double _width;

  public:
    Hadron();
    Hadron(Parton partons);
    Hadron(Parton partons, double xi, double yi, double zi, double ti, double xf, double yf, double zf, double tf, double charge, double width);
    Hadron(Parton partons, double charge, double width);
    virtual ~Hadron();

    virtual void display() const;

    void SetRi(double xi, double yi, double zi, double ti);
    void vSetRi(const std::array<double,4>& ri);
    const std::array<double,4>& GetRi() const;

    void SetRf(double xf, double yf, double zf, double tf);
    void vSetRf(const std::array<double,4>& rf);
    const std::array<double,4>& GetRf() const;

    void SetCharge(double charge);
    double GetCharge() const;

    void SetWidth(double width);
    double GetWidth() const;
};
