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

bool safe_jet(vector<Parton>, double);
void init_tree(int);
void do_tree(vector<Parton>&, double&, double&, double&, bool&, double, double, int);

Pythia pythia;

void init_tree(int seed)
{
  //Read cmnd file
  ostringstream pythiaset;
  //pythiaset << "setup_pythia_" << pdf << "_" << system << ".cmnd";
  pythiaset << "setup_pythia.cmnd";
  pythia.readFile(pythiaset.str());

  //Set Random Seed
  pythia.readString("Random:setSeed = on");
  ostringstream seedstring;
  seedstring << "Random:seed = " << seed;
  pythia.readString(seedstring.str().c_str());

  pythia.init();
}

void do_tree(vector<Parton> &partons, double &weight, double &cross, double &cross_err, bool &have_trigger, double trigger_pt, double trigger_eta, int trigger_id)
{
  //Event& event      = pythia.event;
  //ParticleData& pdt = pythia.particleData;
  
  if (!pythia.next()) {
    return;
  }

  //Look for trigger
  have_trigger = 0;
  for (int i = 0; i < pythia.event.size(); i++)
  {
    if (!pythia.event[i].isFinal()) continue;
    if (abs(pythia.event[i].id()) != trigger_id) continue;
    if (pythia.event[i].pT() < trigger_pt) continue;
    if (fabs(pythia.event[i].eta()) > trigger_eta) continue;
    have_trigger = 1;
    break;
  }
  if (!have_trigger) return;

  //Find Final Particles
  for (int i = 0; i < pythia.event.size(); i++)
  {
    if (pythia.event[i].isFinal())
    {
      vector<double> p;
      for (unsigned int j=1; j<4; j++) p.push_back(pythia.event[i].p()[j]);
      p.push_back(pythia.event[i].p()[0]);

      //Simply store remnants (tf 0, mother 0, daugthers -1)
      if (pythia.event[i].status() == 63)
      {
        partons.push_back ( Parton ( p, 0., pythia.event[i].m(), 0, -1, -1, pythia.event[i].id(), "rem", pythia.event[i].col(), pythia.event[i].acol(), true ) );
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
      //double virt=sqrt(abs(pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2()));
      double tf = 100000000.;
      //cout << " virt= " << virt << endl;	
      //cout << " virt2= " << pow(pythia.event[i].e(),2.)-pow(pythia.event[i].px(),2.)-pow(pythia.event[i].py(),2.)-pow(pythia.event[i].pz(),2.)-pythia.event[i].m2() << endl;
	
      //Add it to partons array
      partons.push_back ( Parton ( p, tf, pythia.event[i].m(), m1, -1, -1, pythia.event[i].id(), "ps", pythia.event[i].col(), pythia.event[i].acol(), false ) );
	
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
	    double virt=sqrt(abs(pow(p[3],2.)-pow(p[0],2.)-pow(p[1],2.)-pow(p[2],2.)-pow(pythia.event[mom].m(),2.)));
		
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
	    double tf = 2.*p[3]/virt/virt*0.2;
	    partons.push_back ( Parton ( p, tf, pythia.event[mom].m(), m1, i, j, pythia.event[mom].id(), "ps", pythia.event[mom].col(), pythia.event[mom].acol(), false ) );

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

  //Set on-shell, assuming massless
  for (unsigned int ip=0; ip<partons.size(); ip++) {
    vector<double> p = partons[ip].vGetP();
    partons[ip].SetP(p[0],p[1],p[2],sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));
  }

  //Extract weight and cross section of the event
  weight=pythia.info.weight();
  cross=pythia.info.sigmaGen();
  cross_err=pythia.info.sigmaErr();
	
}
