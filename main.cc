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

#include "read_tables.hpp"

#define DO_WAKE
#define DO_LUND
#define DO_ELASTIC

char Sfile[100];        //Source File name: one file per event

using std::vector;
using namespace std;

void read_tables(std::string tables_path);

void read_nuclear(int, std::string);
void read_hydro(int, std::string);

void init_tree(int njob);
void do_tree(vector<Parton> &partons, double &weight, double &cross, double &cross_err);

void gxy(double &, double &, numrand &);

void do_eloss(vector<Parton>, vector<Quench> &, double, double, numrand &, double kappa, double alpha, int tmethod, int model, vector<Quench> &, vector<vector<double>> &all_wakes);

void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr, vector<vector<double>> all_wakes);

void init_lund();
void do_lund_vac(vector<Parton> partons, vector<Hadron> &hadrons, int hadro_type);
bool do_lund_med(vector<Quench> partons, vector<Hadron> &hadrons, int hadro_type);

int main(int argc, char** argv)
{
  
  assert(argc==10);
  
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

  cout << " njob= " << njob << " N= " << Nev << " cent= " << cent << " kappa= " << kappa << " alpha= " << alpha << " tmethod= " << tmethod << endl;

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
  numrand nr(1346+njob);
  //cout << " rando= " << nr.rando() << endl;

  if (do_quench) {	
    //Read initial energy density
    read_nuclear(0, cent);
		
    //Read hydro file, event averaged
    read_hydro(0, cent);

  }

#ifdef DO_ELASTIC
  //Get moliere tables
  std::string tables_path="/mnt/lustre/scratch/nlsas/home/usc/ie/dpa/moliere_tables/a10_tables/";
  read_tables(tables_path);
  use_tables = 1;
#endif

  //Hydro Event Loop
  int count=0;
  do {
    
    //Generate PYTHIA tree
    //Declare partons vector
    vector<Parton> partons;
    if (count==0) init_tree(njob);
    double weight, cross, cross_err;
    do_tree(partons,weight,cross,cross_err);
    //cout << " partons.size= " << partons.size() << endl;

    //if (count<80) { count++; continue; }

    //Create vector of quenched partons initially equal to vacuum partons
    vector<Quench> quenched;	
    for (unsigned int i = 0; i<partons.size(); i++)
    {
      quenched.push_back ( Quench ( partons[i] ) );
      //quenched.push_back ( partons[i] );
    }

    // Shower Analysis
    //shower_analysis(partons);

    double x,y;
    vector<Quench> recoiled;
    vector<vector<double>> all_wakes;
    ofstream source_file;
    if (do_quench) {	
      
      //Generate x,y
      gxy(x, y, nr);
      cout << " xcre= " << x << " ycre= " << y << endl;

      /*
      //Print creation point, and four momentum of hs partons in Source file
      sprintf(Sfile,"./SOURCE_Ev%i.dat",count);
      source_file.open(Sfile);
      source_file << "X_cre= " << x << " Y_cre= " << y << endl;
      for (unsigned int i = 0; i<partons.size(); i++)
      {
        if (partons[i].GetOrig()=="hs") source_file << "p_x = " << partons[i].vGetP()[0] << " p_y= " << partons[i].vGetP()[1] << " p_z = " << partons[i].vGetP()[2] << " p_e = " << partons[i].vGetP()[3] << endl;
      }
      */

      do_eloss(partons, quenched, x, y, nr, kappa, alpha, tmethod, model, recoiled, all_wakes);

    }
    else {
      x=0., y=0.;	
    }		

  #ifdef DO_WAKE
    //Do back-reaction, return a vector of wake hadrons
    vector<Wake> wake;
    do_wake(quenched,partons,wake,nr,all_wakes);
    cout << "wake size from jet partons= " << wake.size() << endl;
    //do_wake(recoiled,partons,wake,nr,all_wakes);
    //cout << "wake size after recoil partons= " << wake.size() << endl;
    std::ofstream outfile;
    outfile.open("Injected.dat", std::ios_base::app);
    outfile << "NEXT" << endl;
  #endif
 
  //EM check
  cout << "Event " << count << endl;
  vector<double> vac_p(4,0.);
  for (int iv=0; iv<partons.size(); iv++) {
    if (partons[iv].GetD1()==-1) vac_p += partons[iv].vGetP();
  }
  cout << " Vac p = " << vac_p[0] << " " << vac_p[1] << " " << vac_p[2] << " " << vac_p[3] << endl;
  vector<double> med_p(4,0.);
  for (int iv=0; iv<quenched.size(); iv++) {
    if (quenched[iv].GetD1()==-1) med_p += quenched[iv].vGetP();
  }
  cout << " Med p = " << med_p[0] << " " << med_p[1] << " " << med_p[2] << " " << med_p[3] << endl;
  for (int iv=0; iv<recoiled.size(); iv++) {
    if (recoiled[iv].GetD1()!=-1) { cout << "WTF RECOIL NOT D1=-1" << endl; exit(1); }
    if (recoiled[iv].GetOrig()=="recoiler") med_p += recoiled[iv].vGetP();
    else if (recoiled[iv].GetOrig()=="hole") med_p -= recoiled[iv].vGetP();
    else { cout << "WTF RECOIL ORIG= " << recoiled[iv].GetOrig() << endl; exit(1); }
  }
  cout << " W recoiled Med p = " << med_p[0] << " " << med_p[1] << " " << med_p[2] << " " << med_p[3] << endl;

#ifdef DO_WAKE

  vector<double> wmed_p = med_p;
  for (int iv=0; iv<all_wakes.size(); iv++) {
    wmed_p += all_wakes[iv];
  }
  cout << " With Bef Wakes Med p = " << wmed_p[0] << " " << wmed_p[1] << " " << wmed_p[2] << " " << wmed_p[3] << endl;
  double pt_med = sqrt(wmed_p[0]*wmed_p[0]+wmed_p[1]*wmed_p[1]);
  if (pt_med>10.) {
    cout << endl << "NON CONSERVATION PARTONIC" << endl << endl;
    for (int iv=0; iv<partons.size(); iv++) { cout << " i= " << iv << " "; partons[iv].display(); cout << endl;}
    for (int iv=0; iv<quenched.size(); iv++) { cout << " i = " << iv << " "; quenched[iv].display(); cout << endl;}
    for (int iv=0; iv<all_wakes.size(); iv++) {
      cout << " i= " << iv << " " << all_wakes[iv][0] << " " << all_wakes[iv][1] << " " << all_wakes[iv][2] << " " << all_wakes[iv][3] << endl;
    }
    exit(1);
  }
  
  for (int iv=0; iv<wake.size(); iv++) {
    med_p += wake[iv].vGetP()*double(wake[iv].GetStatus());
  }
  cout << " With Wake Med p = " << med_p[0] << " " << med_p[1] << " " << med_p[2] << " " << med_p[3] << endl;
  pt_med = sqrt(med_p[0]*med_p[0]+med_p[1]*med_p[1]);
  //if (pt_med>10.) {
    //cout << endl << "NON CONSERVATION PARTONIC" << endl << endl;
    //exit(1);
  //}
  //
#endif

  #ifdef DO_LUND
    //Hadronize in pythia, return a vector of pythia hadrons
    if (count==0) init_lund();
    //if (count==80) init_lund();
    
    vector<Hadron> vhadrons, qhadrons, hhadrons;
    
    int vac_hadro_type=0;
    do_lund_vac(partons,vhadrons,vac_hadro_type);

    vector<Quench> quenchandrecoil=quenched;
    vector<Quench> holes;
    for (unsigned int i=0; i<recoiled.size(); i++) {
      if (recoiled[i].GetOrig()=="recoiler") quenchandrecoil.push_back(recoiled[i]);
      else if (recoiled[i].GetOrig()=="hole") holes.push_back(recoiled[i]);
      else {
        cout << " WTF orig= " << recoiled[i].GetOrig() << endl;
	exit(1);
      }
    }
 
    //for (int iv=0; iv<quenched.size(); iv++) { cout << " i = " << iv << " "; quenched[iv].display(); cout << endl;}
    //for (int iv=0; iv<holes.size(); iv++) { cout << " i = " << iv << " "; holes[iv].display(); cout << endl;}

    //Check if med hadron went ok
    int had_counter=0;
    bool had_is_ok=0;
    int had_counter_max=5;
    do {
      qhadrons.clear();
      had_is_ok = do_lund_med(quenchandrecoil,qhadrons,hadro_type);
      had_counter++;
    } while (!had_is_ok && had_counter<had_counter_max);
    if (had_counter>0) {
      cout << "Had Counter = " << had_counter << " and had_is_ok= " << had_is_ok << endl;
    }
    if (!had_is_ok) {
      continue; //Skip event and do not print it
    }

    if (!holes.empty()) do_lund_med(holes,hhadrons,hadro_type);
    cout << " Vac Hadron size= " << vhadrons.size() << " Med Hadron size= " << qhadrons.size() << endl;
 
    
    if (qhadrons.size()==0) {
      for (int iv=0; iv<partons.size(); iv++) {
        cout << " Parton = " << iv << " ";
	partons[iv].display();
	cout << endl;
      }
      cout << endl;
      for (int iv=0; iv<quenched.size(); iv++) {
        cout << " Quench = " << iv << " ";
	quenched[iv].display();
	cout << endl;
      }
      cout << endl;
      for (int iv=0; iv<quenchandrecoil.size(); iv++) {
        cout << " Quench&Recoil = " << iv << " ";
	quenchandrecoil[iv].display();
	cout << endl;
      }
      cout << endl;
      for (int iv=0; iv<holes.size(); iv++) {
        cout << " Holes = " << iv << " ";
	holes[iv].display();
	cout << endl;
      }
      cout << endl;
      for (int iv=0; iv<hhadrons.size(); iv++) {
        cout << " Hole Hadrons = " << iv << " ";
	hhadrons[iv].display();
	cout << endl;
      }
      cout << endl;
    }
    

    //EM hadron level check
    vector<double> hvac_p(4,0.);
    for (int iv=0; iv<vhadrons.size(); iv++) {
      hvac_p += vhadrons[iv].vGetP();
    }
    cout << " Had Vac p = " << hvac_p[0] << " " << hvac_p[1] << " " << hvac_p[2] << " " << hvac_p[3] << endl;
    vector<double> hmed_p(4,0.);
    for (int iv=0; iv<qhadrons.size(); iv++) {
      hmed_p += qhadrons[iv].vGetP();
    }
    cout << " Had Med p = " << hmed_p[0] << " " << hmed_p[1] << " " << hmed_p[2] << " " << hmed_p[3] << endl;

#ifdef DO_WAKE
    for (int iv=0; iv<wake.size(); iv++) {
      hmed_p += wake[iv].vGetP()*double(wake[iv].GetStatus());
    }
    cout << " With Wake Had Med p = " << hmed_p[0] << " " << hmed_p[1] << " " << hmed_p[2] << " " << hmed_p[3] << endl;
#endif

    for (int iv=0; iv<hhadrons.size(); iv++) {
      hmed_p -= hhadrons[iv].vGetP();
    }
    cout << " With Holes Had Med p = " << hmed_p[0] << " " << hmed_p[1] << " " << hmed_p[2] << " " << hmed_p[3] << endl;
    double hmed_pt=sqrt(hmed_p[0]*hmed_p[0]+hmed_p[1]*hmed_p[1]);
    if (hmed_pt>10.) {
      cout << endl << "WATCH OUT" << endl;
    }

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
      for (unsigned int i=0; i<hhadrons.size(); i++) {
        hjt_file << hhadrons[i].vGetP()[0] << " " << hhadrons[i].vGetP()[1] << " " << hhadrons[i].vGetP()[2] << " " << hhadrons[i].GetMass() << " " << hhadrons[i].GetId() << " " << 3 << endl;
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
 
    if (do_quench) {
      source_file.close();
    }

    count+=1;
  } while (count<Nev);
	
  pjt_file.close();
#ifdef DO_LUND
  hjt_file.close();
#endif

  return 0;
}
