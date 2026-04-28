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

 #define DO_WAKE //ONLY COMMENT OUT WHEN WE WANT FAST RUNNING FOR TESTING

using std::vector;
using namespace std;

void read_nuclear(int, std::string);
void read_hydro(int, std::string);

void init_tree(int njob, string pythiafile);
bool do_tree(vector<Parton> &partons, double &weight, double &cross, double &cross_err, int biased_tree);

void gxy(double &, double &, numrand &);

void do_eloss(vector<Parton>, vector<Quench> &, double, double, numrand &, double kappa, double alpha, double lambda, int tmethod, int mode, int new_eloss);

void do_wake(vector<Quench> quenched, vector<Parton> partons, vector<Wake> &wake, numrand &nr);

void init_lund();
void do_lund_vac(vector<Parton> partons, vector<Hadron> &vhadrons, int hadro_type);
bool do_lund_med(vector<Quench> quenched, vector<Hadron> &qhadrons, int hadro_type);

int main(int argc, char** argv)
{

  assert(argc==15);

  int njob=atoi(argv[1]);
  int Nev=atoi(argv[2]);
  std::string cent=argv[3];
  double kappa=atof(argv[4]);
  double alpha=atof(argv[5]);
  double lambda=atof(argv[6]);
  int tmethod=atoi(argv[7]);
  std::stringstream ss(argv[8]);	
  string pythiafile = argv[9];
  string outfileprefix = argv[10];
  bool do_quench;
  if (!(ss >> std::boolalpha >> do_quench)) { cout << " Wrong bool variable for do_quench "; exit(0); }
  int mode=atoi(argv[11]);

  int biased_tree=atoi(argv[12]); //0 for no bias, 1 for b bias, 2 for c bias
  int new_eloss=atoi(argv[13]); //0 for no heavy quark energy loss, 1 for centaur formula, 2 for centaur and heavy broad, 3 for only heavy broad
  int outputmode=atoi(argv[14]); //0 for usual, 1 for streamline D's (for doing v2), 2 for normal+charms, 3 for only partons for coalesence, 4 for streamlined D's and coalesence partons, 5 for jets for coalescence

  cout << " njob= " << njob << " N= " << Nev << " cent= " << cent << " kappa= " << kappa << " alpha= " << alpha << " tmethod= " << tmethod << endl;

  char JToutHad[100];
  sprintf(JToutHad,"%s_Hadrons.out",outfileprefix.c_str());
  ofstream hjt_file;
  hjt_file.open (JToutHad);
	
  ofstream pjt_file;
  if (outputmode==0 || outputmode==2){
    char JToutPart[100];
    sprintf(JToutPart,"%s_Partons.out",outfileprefix.c_str()); //Use same prefix as hadrons
    pjt_file.open (JToutPart);
  } else if (outputmode==3 || outputmode==4 || outputmode==5){
    char JToutPart[100];
    sprintf(JToutPart,"%s_PreCoalesencePartons.out",outfileprefix.c_str()); //Use same prefix as hadrons
    pjt_file.open (JToutPart);
  }

  ofstream chjt_file;
  if (outputmode == 5){
    char JToutCHad[100];
    sprintf(JToutCHad,"%s_NoCoalHadrons.out",outfileprefix.c_str()); //Use same prefix as hadrons
    chjt_file.open (JToutCHad);
  }

  char JToutSkip[100];
  sprintf(JToutSkip,"%s_Skipped.out",outfileprefix.c_str()); //Use same prefix as hadrons
  ofstream sjt_file;
  sjt_file.open (JToutSkip);

  
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
    cout << " Event " << count << " of " << Nev << endl;
    //Generate PYTHIA tree
    //Declare partons vector
    vector<Parton> partons;
    if (count==0) init_tree(njob, pythiafile);
    double weight, cross, cross_err;
    
    bool has_b;
    while (true)
    {
      partons.clear();
      has_b=do_tree(partons,weight,cross,cross_err,biased_tree);
      if (has_b) {break;}
      else {sjt_file << "weight " << weight << " cross " << cross  << endl;}
    }

    double highest_c_pt = 0.;
    for (unsigned int i=0; i<partons.size(); i++) {
      //Print original pair
      // if (partons[i].GetOrig()=="hs") pjt_file << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].GetMass() << " " << partons[i].GetId() << " " << -2 << endl;
      if (abs(partons[i].GetId())==4){
        if (highest_c_pt < sqrt(partons[i].vGetP()[0]*partons[i].vGetP()[0] + partons[i].vGetP()[1]*partons[i].vGetP()[1])) {
          highest_c_pt = sqrt(partons[i].vGetP()[0]*partons[i].vGetP()[0] + partons[i].vGetP()[1]*partons[i].vGetP()[1]);
        }
      }
    }	
    
    //cout << " partons.size= " << partons.size() << endl;
	
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
      do_eloss(partons, quenched, x, y, nr, kappa, alpha, lambda, tmethod, mode, new_eloss);
    }
    else {
      x=0., y=0.;	
    }

    #ifdef DO_WAKE
      //Do back-reaction, return a vector of wake hadrons
      vector<Wake> wake;
      do_wake(quenched,partons,wake,nr);
      // cout << "Wake size= " << wake.size() << endl;
    #endif

    //Hadronize in pythia, return a vector of pythia hadrons
    vector<Hadron> vhadrons, qhadrons;
    if (count==0) init_lund();
    int hadro_type = 0; //0 for vacuum, 1 for medium
    if (outputmode==5) hadro_type = 1; 
    do_lund_vac(partons,vhadrons, hadro_type);
    do_lund_med(quenched,qhadrons, hadro_type);
		vector<Parton> charmlesspartons;
    vector<Quench> charmlessquenched;
    vector<Hadron> charlessqhadrons;
    if (outputmode==5) {
      for (unsigned int i=0; i<partons.size(); i++) {
        if (abs(partons[i].GetId())!=4 && abs(partons[i].GetId())!=5) charmlesspartons.push_back(partons[i]);
      }
      for (unsigned int i=0; i<quenched.size(); i++) {
        if (abs(quenched[i].GetId())!=4 && abs(quenched[i].GetId())!=5) charmlessquenched.push_back(quenched[i]);
      }

      do_lund_med(charmlessquenched, charlessqhadrons, 1); 
    }

    if (outputmode==0 || outputmode==2){
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
    }

    //Print heavy parton ready for coalescence
    if ((outputmode==3 || outputmode==4 || outputmode==5) && do_quench==1){
      for (unsigned int i=0; i<quenched.size(); i++) {
        if (quenched[i].GetD1()==-1 && (quenched[i].vGetP()[0]!=0. || quenched[i].vGetP()[1]!=0. ||quenched[i].vGetP()[2]!=0.) && ((abs(quenched[i].GetId())==4 && biased_tree==2) || (abs(quenched[i].GetId())==5 && biased_tree==1)))
        {
            pjt_file << quenched[i].vGetP()[0] << " " << quenched[i].vGetP()[1] << " " << quenched[i].vGetP()[2] << " " << quenched[i].GetMass() << " " << quenched[i].GetId() << " " << 0 << " " << highest_c_pt << " " << quenched[i].GetFluidVx() << " " << quenched[i].GetFluidVy() << " " << quenched[i].GetFluidVz() << " " << weight << " " << cross << " " << quenched[i].GetFreezoutCrosser() << " " << count << endl;
        }
      }
    }
		
    //Print hadronic
    if (outputmode==0 || outputmode==2 || outputmode==5){
      hjt_file << "# event " << count << endl; 
      hjt_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << " cpt " << highest_c_pt << endl; 
      for (unsigned int i=0; i<partons.size(); i++) {
        if (partons[i].GetOrig()=="hs") hjt_file << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].GetMass() << " " << partons[i].GetId() << " " << -2 << endl;
      }
    } else if (outputmode==3 || outputmode==4){
    }
    
    if (outputmode==5){
      double largest_phase_dist = 0.;
      for (unsigned int i=0; i<qhadrons.size(); i++) {
        if (qhadrons[i].GetOrig() == "promptD") {
          largest_phase_dist = max(largest_phase_dist, qhadrons[i].GetPhaseSpaceDistance());
        }
      }
      chjt_file << "# event " << count << endl;
      chjt_file << "weight " << weight << " cross " << cross << " X " << x << " Y " << y << " cpt " << highest_c_pt << " largest_phase_dist " << largest_phase_dist << endl;
      for (unsigned int i=0; i<partons.size(); i++) {
        if (partons[i].GetOrig()=="hs") chjt_file << partons[i].vGetP()[0] << " " << partons[i].vGetP()[1] << " " << partons[i].vGetP()[2] << " " << partons[i].GetMass() << " " << partons[i].GetId() << " " << -2 << endl;
      }
    }


    if (do_quench==1) {
      for (unsigned int i=0; i<qhadrons.size(); i++) {
        int deca_label = 0;
	      if (qhadrons[i].GetOrig()=="non-promptD") deca_label = 10;
	      if (qhadrons[i].GetOrig()=="promptD") deca_label = 11;
        if (outputmode==1 || outputmode==4){
          if (deca_label != 0 && !(deca_label==10 && biased_tree==2) && !(deca_label==11 && biased_tree==1)){
            double pt= sqrt(qhadrons[i].vGetP()[0]*qhadrons[i].vGetP()[0] + qhadrons[i].vGetP()[1]*qhadrons[i].vGetP()[1]);
            double v2= (qhadrons[i].vGetP()[0]*qhadrons[i].vGetP()[0] - qhadrons[i].vGetP()[1]*qhadrons[i].vGetP()[1])/(pt * pt);
            double yrap= 1./2.*log((sqrt(qhadrons[i].vGetP()[0]*qhadrons[i].vGetP()[0] + qhadrons[i].vGetP()[1]*qhadrons[i].vGetP()[1] + qhadrons[i].vGetP()[2]*qhadrons[i].vGetP()[2] + qhadrons[i].GetMass()*qhadrons[i].GetMass()) + qhadrons[i].vGetP()[2])/(sqrt(qhadrons[i].vGetP()[0]*qhadrons[i].vGetP()[0] + qhadrons[i].vGetP()[1]*qhadrons[i].vGetP()[1] + qhadrons[i].vGetP()[2]*qhadrons[i].vGetP()[2] + qhadrons[i].GetMass()*qhadrons[i].GetMass()) - qhadrons[i].vGetP()[2]));
            double eta = 1./2.*log((sqrt(qhadrons[i].vGetP()[0]*qhadrons[i].vGetP()[0] + qhadrons[i].vGetP()[1]*qhadrons[i].vGetP()[1] + qhadrons[i].vGetP()[2]*qhadrons[i].vGetP()[2]                                              ) + qhadrons[i].vGetP()[2])/(sqrt(qhadrons[i].vGetP()[0]*qhadrons[i].vGetP()[0] + qhadrons[i].vGetP()[1]*qhadrons[i].vGetP()[1] + qhadrons[i].vGetP()[2]*qhadrons[i].vGetP()[2]                                              ) - qhadrons[i].vGetP()[2]));
            double phase_space_distance = qhadrons[i].GetPhaseSpaceDistance();
            hjt_file << weight << " " << pt << " " << v2 << " " << yrap << " " << qhadrons[i].GetId() << " " << deca_label << " " << cross << " " << highest_c_pt << " " << qhadrons[i].GetFreezoutCrosser() << " " << eta << " " << qhadrons[i].GetMass() << " " << phase_space_distance << endl;
          }
        } else {
	        hjt_file << qhadrons[i].vGetP()[0] << " " << qhadrons[i].vGetP()[1] << " " << qhadrons[i].vGetP()[2] << " " << qhadrons[i].GetMass() << " " << qhadrons[i].GetId() << " " << deca_label << endl;
        }
      }
      if (outputmode==5){
        for (unsigned int i=0; i<charlessqhadrons.size(); i++) {
          int deca_label = 0;
          if (charlessqhadrons[i].GetOrig()=="non-promptD") deca_label = 10;
          if (charlessqhadrons[i].GetOrig()=="promptD") deca_label = 11;
          if (deca_label!=0){
            cout << "WARNING, weird deca_label" << endl;
          }
          chjt_file << charlessqhadrons[i].vGetP()[0] << " " << charlessqhadrons[i].vGetP()[1] << " " << charlessqhadrons[i].vGetP()[2] << " " << charlessqhadrons[i].GetMass() << " " << charlessqhadrons[i].GetId() << " " << deca_label << endl;
        }
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
        int deca_label = 0;
	      if (vhadrons[i].GetOrig()=="non-promptD") deca_label = 10;
	      if (vhadrons[i].GetOrig()=="promptD") deca_label = 11;
        if (outputmode==1 || outputmode==4){
          if (deca_label != 0){
            double pt= sqrt(vhadrons[i].vGetP()[0]*vhadrons[i].vGetP()[0] + vhadrons[i].vGetP()[1]*vhadrons[i].vGetP()[1]);
            double v2= (vhadrons[i].vGetP()[0]*vhadrons[i].vGetP()[0] - vhadrons[i].vGetP()[1]*vhadrons[i].vGetP()[1])/(pt * pt);
            double yrap= 1./2.*log((sqrt(vhadrons[i].vGetP()[0]*vhadrons[i].vGetP()[0] + vhadrons[i].vGetP()[1]*vhadrons[i].vGetP()[1] + vhadrons[i].vGetP()[2]*vhadrons[i].vGetP()[2] + vhadrons[i].GetMass()*vhadrons[i].GetMass()) + vhadrons[i].vGetP()[2])/(sqrt(vhadrons[i].vGetP()[0]*vhadrons[i].vGetP()[0] + vhadrons[i].vGetP()[1]*vhadrons[i].vGetP()[1] + vhadrons[i].vGetP()[2]*vhadrons[i].vGetP()[2] + vhadrons[i].GetMass()*vhadrons[i].GetMass()) - vhadrons[i].vGetP()[2]));
            double eta = 1./2.*log((sqrt(vhadrons[i].vGetP()[0]*vhadrons[i].vGetP()[0] + vhadrons[i].vGetP()[1]*vhadrons[i].vGetP()[1] + vhadrons[i].vGetP()[2]*vhadrons[i].vGetP()[2]                                              ) + vhadrons[i].vGetP()[2])/(sqrt(vhadrons[i].vGetP()[0]*vhadrons[i].vGetP()[0] + vhadrons[i].vGetP()[1]*vhadrons[i].vGetP()[1] + vhadrons[i].vGetP()[2]*vhadrons[i].vGetP()[2]                                              ) - vhadrons[i].vGetP()[2]));
            double phase_space_distance = vhadrons[i].GetPhaseSpaceDistance();
            hjt_file << weight << " " << pt << " " << v2 << " " << yrap << " " << vhadrons[i].GetId() << " " << deca_label << " " << cross << " " << highest_c_pt << " " << vhadrons[i].GetFreezoutCrosser() << " " << eta << " " << vhadrons[i].GetMass() << " " << phase_space_distance << endl;
          }
        } else {
	        hjt_file << vhadrons[i].vGetP()[0] << " " << vhadrons[i].vGetP()[1] << " " << vhadrons[i].vGetP()[2] << " " << vhadrons[i].GetMass() << " " << vhadrons[i].GetId() << " " << deca_label << endl;
        }
      }
    }
    if (outputmode==0 || outputmode==2 || outputmode==5){hjt_file << "end" << endl;}

    vhadrons.clear();
    qhadrons.clear();
    partons.clear();
    quenched.clear();	
  #ifdef DO_WAKE
    wake.clear();
  #endif

    count+=1;
	
  } while (count<Nev);

  if (outputmode==0 || outputmode==2 || outputmode==3 || outputmode==4 || outputmode==5){pjt_file.close();}
  hjt_file.close();
  sjt_file.close();
	
  return 0;
}
