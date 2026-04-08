#include "Parton.h"
#include "Quench.h"
#include "Random.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <functional>

#include "vector_operators.h"

using std::vector;
using namespace std;

void do_eloss(vector<Parton> partons, vector<Quench> &quenched, double xcre, double ycre, numrand &nr, double kappa, double alpha, int tmethod, int mode, int ebe_hydro);
void quenched_sons(vector<double> p, vector<double> qp, vector<double> &d1, vector<double> &d2);
void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod, int mode, double &length, double &tlength, int ebe_hydro);

void do_eloss(vector<Parton> partons, vector<Quench> &quenched, double xcre, double ycre, numrand &nr, double kappa, double alpha, int tmethod, int mode, int ebe_hydro)
{
  //Tag final particles
  vector<int> FinId;
  for (unsigned int i = 0; i < quenched.size(); i++) {
    //Exclude remnants
    if (quenched[i].GetD1()==-1 && quenched[i].GetOrig()!="rem")
    {
      FinId.push_back(i);
    }
  }

  //Energy Loss Loop: select a final particle, find its oldest undone parent, climb down the family chain. Iterate until all final particles are done
  for (unsigned int i = 0; i < FinId.size(); i++)
  {
    int ind = FinId[i];  //Start loop with final particle
    vector<int> Fam;  //Family chain array
    Fam.push_back(ind);
    double inhe=1.;    //Used to see whether mother was completely quenched
    //Family chain loop
    int found=0;
    do {
      int mom = quenched[ind].GetMom();
      //If "ind" is not a parent parton
      if (mom!=-1)
      {
        //If mother is done
        if (quenched[mom].GetIsDone()==true)
        {
          inhe=quenched[mom].vGetP()[3];
          //End of family chain
          found=1;
        }
        //If it is not done
        else
        {
          Fam.push_back(mom);
          ind=mom;  
        }
      }
      //If "ind" is a parent parton
      else
      {
        quenched[ind].SetRi(xcre,ycre,0.,0.);
        //End of family chain
        found=1;  
      }
    } while (found==0);
      
    //Apply energy loss chronologically
    for (unsigned int w = Fam.size(); w>0; w--)
    {
      int tp = Fam[w-1];
      //If first done mother was totally quenched, set all descendance quenched and done, and exit loop
      if (inhe==0.)
      {
        for (unsigned int j = w; j>0; j--)
        {
          tp = Fam[j-1];
          quenched[tp].SetP(0.,0.,0.,0.);
          quenched[tp].SetIsDone(true);
        }
        break;
      }
      vector<double> p = quenched[tp].vGetP();
      double q = quenched[tp].GetQ();
      vector<double> pos=quenched[tp].GetRi();
      //Time of flight (from formation time argument)
      double tof = 0.2*2.*p[3]/pow(q,2.); // in fm
      //If final particle, fly arbitrarily far
      if (w==1) tof = 10000000000.;
      double length=0.; //length in QGP
      double tlength=0.; //temperature weighted length in QGP
      //If colored particle
      if (abs(quenched[tp].GetId())<=6 || quenched[tp].GetId()==21)
      {
        //cout << "Quenching" << endl;
	//cout << " En bef= " << p[3] << endl;
        loss_rate(p, pos, tof, quenched[tp].GetId(), nr, kappa, alpha, tmethod, mode, length, tlength, ebe_hydro);
	//cout << " En aft= " << p[3] << endl;
      }
      else
      {
        //If not colored particle, don't do energy loss, but propagate position and time manually
        pos+=p/p[3]*tof;
      }
      //Update mother momenta and set positions to final
      quenched[tp].vSetP(p);
      quenched[tp].vSetRf(pos);
      quenched[tp].SetIsDone(true);
      quenched[tp].AddLength(length, tlength);
      //If mother got fully quenched, quenched descendance and exit
      if (p[3]==0.)
      {
        for (unsigned int j = w; j>0; j--)
        {
          tp = Fam[j-1];
          quenched[tp].SetP(0.,0.,0.,0.);
          quenched[tp].SetIsDone(true);
        }
        break;
      }
      //If not final particle, propagate quenching to son in chain, and store results for other son
      if (w!=1) {
        //Find two daughters
        int d1 = Fam[w-2];
        int d2;
        if (quenched[tp].GetD1()==d1) d2 = quenched[tp].GetD2();
        else d2 = quenched[tp].GetD1();
        //Find new momenta for sons: rotate and quench tri-momentum, quench energy
        vector<double> m_p = partons[tp].vGetP();
        vector<double> d1_p = partons[d1].vGetP();
        vector<double> d2_p = partons[d2].vGetP();
        quenched_sons(m_p, p, d1_p, d2_p);
        //Propagate momenta and positions
        quenched[d1].vSetP(d1_p);
        quenched[d1].vSetInhP(d1_p);
        quenched[d1].vSetRi(pos);
        quenched[d1].AddLength(length, tlength);
        quenched[d2].vSetP(d2_p);
        quenched[d2].vSetInhP(d2_p);
        quenched[d2].vSetRi(pos);
        quenched[d2].AddLength(length, tlength);
      }
    }
    Fam.clear();
  }
  FinId.clear();

  return;
//End of function
}
