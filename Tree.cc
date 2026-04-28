#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include <iostream>
#include <vector>
#include <fstream>

#include "Parton.h"

#include "vector_operators.h"

using std::vector;
using namespace std;
using namespace Pythia8;

void init_tree(int, string);
bool do_tree(vector<Parton>&, double&, double&, double&, int);

Pythia pythia;

void init_tree(int njob, string pythiafile)
{
  //Read cmnd file
  ostringstream pythiaset;
  pythiaset << pythiafile;
  pythia.readFile(pythiaset.str());

  int seed=33+njob;
  ostringstream seedstring;
  seedstring << "Random:seed = " << seed;
  pythia.readString(seedstring.str().c_str());

  pythia.init();

  return;
}

bool do_tree(vector<Parton> &partons, double &weight, double &cross, double &cross_err, int biased_tree)
{
  const Info& info = pythia.info; //Added const to fix error
  // Info& info = pythia.info;
  bool hasHard = pythia.settings.hasHardProc();

  if (!pythia.next()) {
    return false;
  }
  // /*
  if (!hasHard){
    double pTHat  = info.pTHat();
    if (pTHat>20 && info.isNonDiffractive()){//Shouldn't be hardcoded forever
      weight=info.weight();
      cross=info.sigmaGen();
      cross_err=info.sigmaErr();
      return false;
    }
  }
  // */

  if (biased_tree == 1) { //Biased towards b-quarks
  bool has_b=false;
    for (int i = 0; i < pythia.event.size(); i++)
    {
      if (pythia.event[i].isFinal())
      {
        //Chech for b-parton
        if (abs(pythia.event[i].id())==5)
        {
          has_b=true;
          break;
        }
      }
    }
    if (!has_b) {
      weight=info.weight();
      cross=info.sigmaGen();
      cross_err=info.sigmaErr();
      return false;
    }
  } else if (biased_tree == 2) { //Biased towards c-quarks
  bool has_c=false;
    for (int i = 0; i < pythia.event.size(); i++)
    {
      if (pythia.event[i].isFinal())
      {
        //Chech for c-parton
        if (abs(pythia.event[i].id())==4)
        {
          has_c=true;
          break;
        }
      }
    }
    if (!has_c) {
      weight=info.weight();
      cross=info.sigmaGen();
      cross_err=info.sigmaErr();
      return false;
    }
  } else if (biased_tree == 3) { //Biased towards b,c-quarks
  bool has_bc=false;
    for (int i = 0; i < pythia.event.size(); i++)
    {
      if (pythia.event[i].isFinal())
      {
        //Chech for c-parton
        if (abs(pythia.event[i].id())==4 || abs(pythia.event[i].id())==5)
        {
          has_bc=true;
          break;
        }
      }
    }
    if (!has_bc) {
      weight=info.weight();
      cross=info.sigmaGen();
      cross_err=info.sigmaErr();
      return false;
    }
  }
  //Find Final Particles
  for (int i = 0; i < pythia.event.size(); i++)
  {
    if (pythia.event[i].isFinal())
    {
      vector<double> p;
      for (unsigned int j=1; j<4; j++) p.push_back(pythia.event[i].p()[j]);
      p.push_back(pythia.event[i].p()[0]);

      //Simply store remnants (virt 0, mother 0, daugthers -1)
      if (pythia.event[i].status() == 63)
      {
        partons.push_back ( Parton ( p, 0., pythia.event[i].m(), 0, -1, -1, pythia.event[i].id(), "rem", pythia.event[i].col(), pythia.event[i].acol(), true, false ) );
        continue;
      }

      //Find first non-trivial mother, skipping carbon copies
      int use = i;
      int m1 = 0;
      int m2 = 0;
      do
      {
        m1 = pythia.event[use].mother1();
        m2 = pythia.event[use].mother2();
        if (m1==m2) use = m1;
      } while (m1==m2);
		
      //Compute virtuality
      double virt=sqrt(fabs(pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2()));
      //cout << " virt= " << virt << endl;	
      // if (pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2()<0) {
      //   cout << " virt= " << virt << " signed virt= " << pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2() << " e= " << pythia.event[i].e() << " px= " << pythia.event[i].px() << " py= " << pythia.event[i].py() << " pz= " << pythia.event[i].pz() << " m2= " << pythia.event[i].m2() << " id= " << pythia.event[i].id() << endl;
      // }
      //cout << " virt2= " << pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2() << endl;
	
      //Add it to partons array
      partons.push_back ( Parton ( p, virt, pythia.event[i].m(), m1, -1, -1, pythia.event[i].id(), "ps", pythia.event[i].col(), pythia.event[i].acol(), false, false ) );
	
      //If mother is ISR initiator, say mother is -1
      if (pythia.event[m1].status()==-41)
      {
        partons[partons.size()-1].SetMom(-1);
        partons[partons.size()-1].SetOrig("isr");
        partons[partons.size()-1].SetIsDone(true);
      }
		
      //If mother is hard scattering initiator, say mother is -1
      if (pythia.event[m1].status()==-21)
      {
        //cout << " ahora " << endl;
        partons[partons.size()-1].SetMom(-1);
        partons[partons.size()-1].SetOrig("hs");
        partons[partons.size()-1].SetIsDone(true);
      }
    }
  }
	
  //Reconstruct tree, using momentum of daughters to get the one for mothers
  int changes=0;
  do
  {
    changes=1;
    unsigned int ps = partons.size();
    for (unsigned int i = 0; i < ps; i++)
    {
      if (partons[i].GetIsDone()==false)
      {
        for (unsigned int j = 0; j < ps; j++)
        {
	  if (partons[i].GetMom()==partons[j].GetMom() && i!=j && partons[j].GetIsDone()==false)
	  {
	    int mom=partons[i].GetMom();
		
            //Set mother momentum and virtuality
            vector<double> p=partons[i].vGetP()+partons[j].vGetP();
	    double virt=sqrt(fabs(pow(p[3],2.)-pow(p[0],2.)-pow(p[1],2.)-pow(p[2],2.)-pow(pythia.event[mom].m(),2.)));
      // if (pow(p[3],2.)-pow(p[0],2.)-pow(p[1],2.)-pow(p[2],2.)-pow(pythia.event[mom].m(),2.)<0) {
      //   cout << " virt= " << virt << " signed virt= " << pow(p[3],2.)-pow(p[0],2.)-pow(p[1],2.)-pow(p[2],2.)-pow(pythia.event[mom].m(),2.) << " e= " << p[3] << " px= " << p[0] << " py= " << p[1] << " pz= " << p[2] << " m2= " << pythia.event[mom].m() << " id= " << pythia.event[mom].id() << endl;
      // }
		
            //Find first non-trivial mother of the mother
            int use = mom;
            int m1 = 0;
            int m2 = 0;
            do {
              m1 = pythia.event[use].mother1();
              m2 = pythia.event[use].mother2();
              if (m1==m2) use = m1;
            } while (m1==m2);
		
	    //Fill it in partons array
	    partons.push_back ( Parton ( p, virt, pythia.event[mom].m(), m1, i, j, pythia.event[mom].id(), "ps", pythia.event[mom].col(), pythia.event[mom].acol(), false, false ) );

	    //Update mother of daughters to point to position in partons array instead of pythia list, and declare as done
	    partons[i].SetMom(partons.size()-1);
	    partons[j].SetMom(partons.size()-1);
	    partons[i].SetIsDone(true);
	    partons[j].SetIsDone(true);
		
	    //If mother is ISR initiator, say mother of mother is -1
	    if (pythia.event[m1].status()==-41)
	    {
	      partons[partons.size()-1].SetMom(-1);
	      partons[partons.size()-1].SetOrig("isr");
              partons[partons.size()-1].SetIsDone(true);
            }
		
	    //If mother is hard scattering initiator, say mother of mother is -1
	    if (pythia.event[m1].status()==-21)
            {
	      partons[partons.size()-1].SetMom(-1);
              partons[partons.size()-1].SetOrig("hs");
	      partons[partons.size()-1].SetIsDone(true);
	    }

	    changes=0;
            break;
	  }
        }
      }
    }
  } while (changes==0);
  //cout << " Partons size= " << partons.size() << endl;

  //Extract weight and cross section of the event
  weight=info.weight();
  cross=info.sigmaGen();
  cross_err=info.sigmaErr();

  return true;
}
