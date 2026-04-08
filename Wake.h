#pragma once

#include <vector>
#include <string>

using std::vector;

class Wake
{
  private:
    vector<double> _p;
    
    double _mass;
    int _charge;  

    int _id;
    double _status;

    int _mom;

  public:
    Wake();
    Wake(vector<double> p, double mass, int charge, int id, double status);
    ~Wake();

    void display() const;

    void vSetP(vector<double> p);
    vector<double> vGetP() const;

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
