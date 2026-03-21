#pragma once

#include <vector>
#include <string>

using std::vector;

class Parton
{
  protected:
    vector<double> _p;

    double _q;
    double _mass;
  
    vector<double> _ri;

    int _mom;
    int _d1;
    int _d2;

    int _id;
    std::string _orig;

    int _col;
    int _acol;

    double _length;
    double _tlength;

    bool _isdone;

  public:
    Parton();
    Parton(vector<double> p, double q, double mass, int mom, int d1, int d2, int id, std::string orig, int col, int acol, bool isdone);
    virtual ~Parton();

    virtual void display() const;

    void AddLength(double length, double tlength) { _length+=length; _tlength+=tlength; }
    double length() { return _length; }
    double tlength() { return _tlength; }

    void vSetP(vector<double> p);
    void SetP(double px, double py, double pz, double en);
    vector<double> vGetP() const;

    double GetPt() const;

    double GetEta() const;

    void SetQ(double q);
    double GetQ() const;

    void SetMass(double mass);
    double GetMass() const;

    void SetMom(int mom);
    int GetMom() const;

    void SetD1(int d1);
    int GetD1() const;

    void SetD2(int d2);
    int GetD2() const;

    void SetId(int id);
    int GetId() const;

    void vSetRi(vector<double> ri);
    vector<double> GetRi() const;

    void SetOrig(std::string orig);
    std::string GetOrig() const;

    void SetCol(int col);
    int GetCol() const;

    void SetAcol(int acol);
    int GetAcol() const;

    void SetIsDone(bool isdone);
    bool GetIsDone() const;
};
