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

void do_eloss(vector<Parton> partons, vector<Quench> &quenched, double xcre, double ycre, numrand &nr, double kappa, double alpha, int tmethod, int model, vector<Quench> &recoiled);
void quenched_sons(vector<double> p, vector<double> qp, vector<double> &d1, vector<double> &d2);
void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod, int mode, vector<Quench> &new_particles, int &had_scattering, vector<double> &orient, double wake_dep[4], double wake_rec[4]);

void do_eloss(vector<Parton> partons, vector<Quench> &quenched, double xcre, double ycre, numrand &nr, double kappa, double alpha, int tmethod, int model, vector<Quench> &recoiled, vector<vector<double>> &all_wakes)
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

  vector<Quench> new_particles;
  //Energy Loss Loop: select a final particle, find its oldest undone parent, climb down the family chain. Iterate until all final particles are done
  for (unsigned int i = 0; i < FinId.size(); i++)
  {
    int ind = FinId[i];			//Start loop with final particle
    vector<int> Fam;			//Family chain array
    Fam.push_back(ind);
    double inhe=1.;			//Used to see whether mother was completely quenched
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
      vector<double> pos=quenched[tp].GetRi();
      //Time of flight (from formation time argument)
      double tof = quenched[tp].GetTf();
      //If final particle, fly arbitrarily far
      if (w==1) tof = 100000000.;
      //If colored particle
      int had_scattering=quenched[tp].had_scattering();
      vector<double> orient=quenched[tp].orient();
      double ein=p[3];
      double wake_dep[4]={0.};
      double wake_rec[4]={0.};
      if (abs(quenched[tp].GetId())<=6 || quenched[tp].GetId()==21)
      {
        loss_rate(p, pos, tof, quenched[tp].GetId(), nr, kappa, alpha, tmethod, model, new_particles, had_scattering, orient, wake_dep, wake_rec);
        //cout << " Eloss= " << ein-p[3] << endl;  
        vector<double> wakey;
	for (int oo=0; oo<4; oo++) {
          wakey.push_back(wake_dep[oo]);
	}
	all_wakes.push_back(wakey);
	for (int oo=0; oo<4; oo++) {
          wakey[oo]=wake_rec[oo];
	}
	all_wakes.push_back(wakey);
      }
      else
      {
        //If not colored particle, don't do energy loss, but propagate position and time manually
	pos+=p/p[3]*tof;
      }
      //cout << " length= " << length << " tof= " << tof << endl;
      //Update mother momenta and set positions to final
      quenched[tp].vSetP(p);
      quenched[tp].vSetRf(pos);
      quenched[tp].set_orient(orient);
      quenched[tp].set_had_scattering(had_scattering);
      quenched[tp].SetIsDone(true);
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
	vector<double> orient_d1;
	orient_d1.push_back(d1_p[0]/d1_p[3]), orient_d1.push_back(d1_p[1]/d1_p[3]), orient_d1.push_back(d1_p[2]/d1_p[3]), orient_d1.push_back(d1_p[3]/d1_p[3]);
	quenched[d1].set_orient(orient_d1);
	quenched[d1].vSetP(d1_p);
	quenched[d1].vSetInhP(d1_p);
	quenched[d1].vSetRi(pos);
	vector<double> orient_d2;
	orient_d2.push_back(d2_p[0]/d2_p[3]), orient_d2.push_back(d2_p[1]/d2_p[3]), orient_d2.push_back(d2_p[2]/d2_p[3]), orient_d2.push_back(d2_p[3]/d2_p[3]);
	quenched[d2].set_orient(orient_d2);
	quenched[d2].vSetP(d2_p);
        quenched[d2].vSetInhP(d2_p);
        quenched[d2].vSetRi(pos);
	if (had_scattering == 1 || had_scattering == 2) {
          quenched[d1].set_had_scattering(2);
          quenched[d2].set_had_scattering(2);
	}
      }
    }
    Fam.clear();
  }
  FinId.clear();

  //Ignore recoils!
  //new_particles.clear();

  //Quench new particles
  //cout << " new_particles.size= " << new_particles.size() << endl;
  while (true) {
    vector<Quench> current_particles = new_particles;
    new_particles.clear();
    for (unsigned int ip=0; ip<current_particles.size(); ip++) {
     //cout << " Dealing with recoils " << endl;
     if (current_particles[ip].GetOrig()=="hole") {
        current_particles[ip].SetIsDone(true);
        recoiled.push_back(current_particles[ip]);
	continue;
      }
      vector<double> p = current_particles[ip].vGetP();
      vector<double> pos = current_particles[ip].GetRi();
      //Time of flight (from formation time argument)
      double tof = current_particles[ip].GetTf();
      vector<double> orient=current_particles[ip].orient();
      vector<double> orig_en = p;
      //double ebef=p[3];
      int had_scattering=0;
      double wake_dep[4]={0.};
      double wake_rec[4]={0.};
      loss_rate(p, pos, tof, current_particles[ip].GetId(), nr, kappa, alpha, tmethod, model, new_particles, had_scattering, orient, wake_dep, wake_rec);
      vector<double> wakey;
      for (int oo=0; oo<4; oo++) {
        wakey.push_back(wake_dep[oo]);
      }
      all_wakes.push_back(wakey);
      for (int oo=0; oo<4; oo++) {
        wakey[oo]=wake_rec[oo];
      }
      all_wakes.push_back(wakey);
      //cout << "energy aft= " << p[3] << endl;
      //double qfac=ebef;
      current_particles[ip].SetOrigEn(orig_en);
      //Update mother momenta and set positions to final
      current_particles[ip].vSetP(p);
      current_particles[ip].vSetRf(pos);
      current_particles[ip].set_orient(orient);
      current_particles[ip].vSetRf(pos);
      current_particles[ip].set_had_scattering(had_scattering);
      current_particles[ip].SetIsDone(true);
      if (p[3]==0.) current_particles[ip].SetP(0.,0.,0.,0.);
      recoiled.push_back(current_particles[ip]);
    }
    if (new_particles.size()==0) break;
    else {
      cout << " RESCATTERING! " << endl;
    }
  }

}
