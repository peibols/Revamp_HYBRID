#include "Parton.h"
#include "Quench.h"
#include "Hadron.h"
#include "Random.h"
#include "Wake.h"
#include <cmath>
#include <cstdlib> 
#include <fstream>
#include <iostream>
#include <assert.h>
#include <vector>
#include <algorithm>
#include <sstream>

#include "vector_operators.h"

#define DO_WAKE

using std::vector;
using namespace std;

void read_nuclear(int, std::string);
void read_hydro(int, std::string);

void init_tree(int njob);
void do_tree(vector<Parton> &partons, double &weight, double &cross, double &cross_err, bool &have_trigger, double trigger_pt, double trigger_eta, int trigger_id);

void gxy(double &, double &, numrand &);

void do_eloss(vector<Parton>, vector<Quench> &, double, double, numrand &, double kappa, double alpha, int tmethod, int mode);

void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr);

void init_lund();
void do_lund(vector<Parton> partons, vector<Quench> quenched, vector<Hadron> &vhadrons, vector<Hadron> &qhadrons);

int main(int argc, char** argv)
{

  assert(argc==12);

  int njob=atoi(argv[1]);
  int Nev=atoi(argv[2]);
  std::string cent=argv[3];
  double kappa=atof(argv[4]);
  double alpha=atof(argv[5]);
  int tmethod=atoi(argv[6]);
  std::stringstream ss(argv[7]);	
  bool do_quench;
  if (!(ss >> std::boolalpha >> do_quench)) { cout << " Wrong bool variable for do_quench "; exit(0); }
  int mode=atoi(argv[8]);
  double trigger_pt = atof(argv[9]);
  double trigger_eta = atof(argv[10]);
  int trigger_id = atoi(argv[11]);
	
  cout << " njob= " << njob << " N= " << Nev << " cent= " << cent << " kappa= " << kappa << " alpha= " << alpha << " tmethod= " << tmethod << endl;

  char JToutHad[100];
  sprintf(JToutHad,"./HYBRID_Hadrons.out");
  ofstream hjt_file;
  hjt_file.open (JToutHad);
	
  char JToutPart[100];
  sprintf(JToutPart,"./HYBRID_Partons.out");
  ofstream pjt_file;
  pjt_file.open (JToutPart);

  //Initialize Random Seed
  numrand nr(1346+njob);

  if (do_quench) {
    //Read initial energy density
    read_nuclear(0, cent);

    //Read hydro file, event averaged
    read_hydro(0, cent);
  }	

  //Hydro Event Loop
  int count=0;
  do {

    //Generate PYTHIA tree
    //Declare partons vector
    vector<Parton> partons;
    if (count==0) init_tree(njob);
    double weight, cross, cross_err;
    
    bool have_trigger = 0;
    while (true) {

      do_tree(partons,weight,cross,cross_err,have_trigger,trigger_pt,trigger_eta,trigger_id);
      if (have_trigger) break;
      //cout << " partons.size= " << partons.size() << endl;

    }
	
    //Create vector of quenched partons initially equal to vacuum partons
    vector<Quench> quenched;	
    for (unsigned int i = 0; i<partons.size(); i++)
    {
      quenched.push_back ( Quench ( partons[i] ) );
    }

    double x,y;
    ofstream source_file;	
    if (do_quench) {	
      //Generate x,y
      gxy(x, y, nr);
      cout << " xcre= " << x << " ycre= " << y << endl;

      do_eloss(partons, quenched, x, y, nr, kappa, alpha, tmethod, mode);
		
    }
    else {
      x=0., y=0.;	
    }		

  #ifdef DO_WAKE
    //Do back-reaction, return a vector of wake hadrons
    vector<Wake> wake;
    do_wake(quenched,partons,wake,nr);
    cout << "Wake size= " << wake.size() << endl;
  #endif

    //Hadronize in pythia, return a vector of pythia hadrons
    vector<Hadron> vhadrons, qhadrons;
    if (count==0) init_lund();
    do_lund(partons,quenched,vhadrons,qhadrons);
    cout << " Vac Hadron size= " << vhadrons.size() << " Med Hadron size= " << qhadrons.size() << endl;
		
    //Print partonic
    pjt_file << "# event " << count << endl; 
    pjt_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << endl; 
    for (unsigned int i=0; i<partons.size(); i++) {
      //Print original pair
      if (partons[i].GetOrig()=="hs") pjt_file << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].GetMass() << " " << partons[i].GetId() << " " << -2 << endl;
    }			
      
    if (do_quench==1) {
      for (unsigned int i=0; i<quenched.size(); i++) {
        if (quenched[i].GetD1()==-1 && quenched[i].vGetP()[3]!=0.) pjt_file << quenched[i].vGetP()[0] << " " << quenched[i].vGetP()[1] << " " << quenched[i].vGetP()[2] << " " << quenched[i].GetMass() << " " << quenched[i].GetId() << " " << 0 << endl;
      }
    }
    else {
      for (unsigned int i=0; i<partons.size(); i++) {
        if (partons[i].GetD1()==-1) pjt_file << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].GetMass() << " " << partons[i].GetId() << " " << 0 << endl;
      }
    }
    pjt_file << "end" << endl;
		
    //Print hadronic
    hjt_file << "# event " << count << endl; 
    hjt_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << endl; 
    for (unsigned int i=0; i<partons.size(); i++) {
      if (partons[i].GetOrig()=="hs") hjt_file << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].GetMass() << " " << partons[i].GetId() << " " << -2 << endl;
    }
			
    if (do_quench==1) {
      for (unsigned int i=0; i<qhadrons.size(); i++) {
	hjt_file << qhadrons[i].vGetP()[0] << " " << qhadrons[i].vGetP()[1] << " " << qhadrons[i].vGetP()[2] << " " << qhadrons[i].GetMass() << " " << qhadrons[i].GetId() << " " << 0 << endl;
      }

    #ifdef DO_WAKE
      for (unsigned int i=0; i<wake.size(); i++)
      {
        int ide_jt;
	if (int(wake[i].GetStatus())==1) ide_jt=1;
	else ide_jt=2;
	int wake_id;
	double wake_ch=wake[i].GetCharge();
	if (wake[i].GetMass()<0.5) {
	  if (wake_ch==0.) wake_id=111;
	  else if (wake_ch==1.) wake_id=211;
	  else wake_id=-211;
	}
	else {
	  if (wake_ch==1.) wake_id=2212;
	  else wake_id=-2212;
	}
	hjt_file << wake[i].vGetP()[0] << " " << wake[i].vGetP()[1] << " " << wake[i].vGetP()[2] << " ";
	hjt_file << wake[i].GetMass() << " " << wake_id << " " << ide_jt << endl; 
      }
    #endif

    }
      
    else {
      for (unsigned int i=0; i<vhadrons.size(); i++) {
	hjt_file << vhadrons[i].vGetP()[0] << " " << vhadrons[i].vGetP()[1] << " " << vhadrons[i].vGetP()[2] << " " << vhadrons[i].GetMass() << " " << vhadrons[i].GetId() << " " << 0 << endl;
      }
    }
    hjt_file << "end" << endl;

    vhadrons.clear();
    qhadrons.clear();
    partons.clear();
    quenched.clear();	
  #ifdef DO_WAKE
    wake.clear();
  #endif

    count+=1;
	
  } while (count<Nev);

  pjt_file.close();
  hjt_file.close();
	
  return 0;
}
