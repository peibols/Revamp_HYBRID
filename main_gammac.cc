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

#include "global.h"
#include "vector_operators.h"

#define DO_WAKE
#define DO_LUND

char Sfile[100];        //Source File name: one file per event

using std::vector;
using namespace std;

namespace {
constexpr int kShowerSeedOffset = 33;
constexpr int kHybridSeedOffset = 1346;
constexpr int kLundSeedOffset = 2337;
}

void read_nuclear(int, std::string);
void read_hydro(int, std::string);

void init_tree(int seed);
void do_tree(vector<Parton> &partons, double &weight, double &cross, double &cross_err);

void gxy(double &, double &, numrand &);

void do_eloss(vector<Parton>, vector<Quench> &, double, double, numrand &, double kappa, double alpha, int tmethod, int model, vector<Quench> &);

void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr);

void init_lund(int seed);
void do_lund(vector<Parton> partons, vector<Quench> quenched, vector<Quench> recoiled, vector<Hadron> &vhadrons, vector<Hadron> &qhadrons, int hadro_type);

int main(int argc, char** argv)
{
  
  if (argc != 10 && argc != 11) {
    cout << "Usage: " << argv[0] << " njob Nev cent kappa alpha tmethod do_quench model hadro_type [seed_base]" << endl;
    return 1;
  }
  
  int njob=atoi(argv[1]);
  int Nev=atoi(argv[2]);
  std::string cent=argv[3];
  double kappa=atof(argv[4]);
  double alpha=atof(argv[5]);
  int tmethod=atoi(argv[6]);
  std::stringstream ss(argv[7]);	
  bool do_quench;
  if (!(ss >> std::boolalpha >> do_quench)) { cout << " Wrong bool variable for do_quench "; exit(0); }
  int model=atoi(argv[8]);
  int hadro_type=atoi(argv[9]);
  bool use_explicit_seed_base = (argc == 11);
  int seed_base = use_explicit_seed_base ? atoi(argv[10]) : njob;
  int shower_seed = seed_base + kShowerSeedOffset;
  int hybrid_seed = seed_base + kHybridSeedOffset;
  int lund_seed = use_explicit_seed_base ? (seed_base + kLundSeedOffset) : 0;
	
  cout << " njob= " << njob << " N= " << Nev << " cent= " << cent << " kappa= " << kappa << " alpha= " << alpha << " tmethod= " << tmethod << endl;
  cout << " seed_base= " << seed_base << " shower= " << shower_seed << " hybrid= " << hybrid_seed << " lund= " << lund_seed;
  if (!use_explicit_seed_base) cout << " (legacy default)";
  cout << endl;

#ifdef DO_LUND
  char JToutHad[100];
  sprintf(JToutHad,"./HYBRID_Hadrons_%i.out",njob);
  ofstream hjt_file;
  hjt_file.open (JToutHad);
#endif	

  char JToutPart[100];
  sprintf(JToutPart,"./HYBRID_Partons_%i.out",njob);
  ofstream pjt_file;
  pjt_file.open (JToutPart);

  //Initialize Random Seed
  numrand nr(hybrid_seed);
  //cout << " rando= " << nr.rando() << endl;

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
    if (count==0) init_tree(shower_seed);
    double weight, cross, cross_err;
    bool gotc=0;
    int ntry=0;
    do {
      partons.clear();
      do_tree(partons,weight,cross,cross_err);
      for (unsigned int i=0; i<partons.size(); i++) {
        if (partons[i].GetOrig()=="hs" && abs(partons[i].GetId())==4) gotc=1;
      }
      ntry++;
    } while (gotc==0);
    cout << "ntry= " << ntry << endl;
    //cout << " partons.size= " << partons.size() << endl;
	
    //Create vector of quenched partons initially equal to vacuum partons
    vector<Quench> quenched;	
    for (unsigned int i = 0; i<partons.size(); i++)
    {
      quenched.push_back ( Quench ( partons[i] ) );
    }

    // Shower Analysis
    //shower_analysis(partons);

    double x,y;
    vector<Quench> recoiled;
    if (do_quench) {	
      
      //Generate x,y
      gxy(x, y, nr);
      cout << " xcre= " << x << " ycre= " << y << endl;

      do_eloss(partons, quenched, x, y, nr, kappa, alpha, tmethod, model, recoiled);

    }
    else {
      x=0., y=0.;	
    }		

  #ifdef DO_WAKE
    //Do back-reaction, return a vector of wake hadrons
    vector<Wake> wake;
    do_wake(quenched,partons,wake,nr);
    cout << "wake size from jet partons= " << wake.size() << endl;
    do_wake(recoiled,partons,wake,nr);
    cout << "wake size after recoil partons= " << wake.size() << endl;
  std::ofstream outfile;
  outfile.open("Injected.dat", std::ios_base::app);
  outfile << "NEXT" << endl;
  #endif

  #ifdef DO_LUND
    //Hadronize in pythia, return a vector of pythia hadrons
    vector<Hadron> vhadrons, qhadrons;
    if (count==0) init_lund(lund_seed);
    do_lund(partons,quenched,recoiled,vhadrons,qhadrons,hadro_type);
    cout << " Vac Hadron size= " << vhadrons.size() << " Med Hadron size= " << qhadrons.size() << endl;
  #endif	

    //Print partonic
    pjt_file << "# event " << count << endl; 
    pjt_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << endl; 
    for (unsigned int i=0; i<partons.size(); i++) {
      if (partons[i].GetOrig()=="hs") pjt_file << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].GetMass() << " " << partons[i].GetId() << " " << -2 << endl;
    }
    if (do_quench==1) {
      for (unsigned int i=0; i<quenched.size(); i++) {
        if (quenched[i].GetD1()==-1 && quenched[i].vGetP()[3]!=0.) pjt_file << quenched[i].vGetP()[0] << " " << quenched[i].vGetP()[1] << " " << quenched[i].vGetP()[2] << " " << quenched[i].GetMass() << " " << quenched[i].GetId() << " " << quenched[i].had_scattering() << endl;
        //if (quenched[i].GetD1()==-1) pjt_file << quenched[i].vGetP()[0] << " " << quenched[i].vGetP()[1] << " " << quenched[i].vGetP()[2] << " " << quenched[i].GetMass() << " " << quenched[i].GetId() << " " << quenched[i].had_scattering() << endl;
      }
      for (unsigned int i=0; i<recoiled.size(); i++) {
        int recid=-1000;
	if (recoiled[i].GetOrig()=="recoiler") recid=3;
	else if (recoiled[i].GetOrig()=="hole") recid=4;
	else {
          cout << " Unrecognised orig= " << recoiled[i].GetOrig() << endl;
	  exit(1);
	}
        if (recoiled[i].GetD1()==-1 && recoiled[i].vGetP()[3]!=0.) pjt_file << recoiled[i].vGetP()[0] << " " << recoiled[i].vGetP()[1] << " " << recoiled[i].vGetP()[2] << " " << recoiled[i].GetMass() << " " << recoiled[i].GetId() << " " << recid << endl;
      }
    }
    else {
      for (unsigned int i=0; i<partons.size(); i++) {
        if (partons[i].GetD1()==-1) pjt_file << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].GetMass() << " " << partons[i].GetId() << " " << 0 << endl;
      }
    }
    pjt_file << "end" << endl;

  #ifdef DO_LUND    
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
      for (unsigned int i=0; i<wake.size(); i++) {
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
  #endif

    //Compute jet observables, choose part or had
    //std::string part_or_had = "part";
    //jet_obs(partons,quenched,vhadrons,qhadrons,part_or_had,cent,tmethod,alpha,0.3,system,weight,cross,pdf);

  #ifdef DO_LUND
    vhadrons.clear();
    qhadrons.clear();
  #endif
    partons.clear();
    quenched.clear();
    recoiled.clear();
  #ifdef DO_WAKE
    wake.clear();
  #endif
  
    count+=1;
  } while (count<Nev);
	
  pjt_file.close();
#ifdef DO_LUND
  hjt_file.close();
#endif

  return 0;
}
