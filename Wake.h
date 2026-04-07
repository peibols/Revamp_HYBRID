#pragma once

#include <array>
#include <vector>
#include <string>

using std::vector;

class Wake
{
  private:
    std::array<double,4> _p;
    
    double _mass;
    int _charge;  

    int _id;
    double _status;

    int _mom;

  public:
    Wake();
    Wake(std::array<double,4> p, double mass, int charge, int id, double status);
    ~Wake();

    void display() const;

    void vSetP(const std::array<double,4>& p);
    const std::array<double,4>& vGetP() const;

    void SetMass(double mass);
    double GetMass() const;

    void SetCharge(int charge);
    int GetCharge() const;

    void SetId(int id);
    int GetId() const;

    void SetMom(int mom);
    int GetMom() const;
    
    void SetStatus(double status);
    double GetStatus() const;
};
