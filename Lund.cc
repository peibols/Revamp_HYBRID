#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include <iostream>
#include <vector>
#include <fstream>

#include "Parton.h"
#include "Quench.h"
#include "Hadron.h"

using std::vector;
using namespace std;
using namespace Pythia8;

void do_lund_vac(vector<Parton> partons, vector<Hadron> &vhadrons, int hadro_type);
bool do_lund_med(vector<Quench> quenched, vector<Hadron> &qhadrons, int hadro_type);
bool isQuarkOrDiquark(const Particle& p);
int findFragmentationVertex(const Event& event, int index);
int climbToInitialQuark(const Event& event, int index);
double phaseSpaceDistance(const Event& event, int index);
void init_lund();

Pythia hpythia;

void init_lund()
{
  hpythia.readString("Random:setSeed = on");
  hpythia.readString("Random:seed = 0");

  // 1) Check whether “decays.cmnd” is in the working directory
  std::ifstream fin("decays.cmnd");
  if (fin) {
    // If the file opened successfully, load it
    hpythia.readFile("decays.cmnd");
    fin.close();
    std::cout << "[init_lund] Loaded decay‐switches from decays.cmnd\n";
  } else {
    // Otherwise, skip it—everything stays at Pythia’s default
    std::cout << "[init_lund] decays.cmnd not found; using default decay settings\n";
  }
  
  hpythia.readString("Print:quiet = off");
  hpythia.readString("111:mayDecay = off"); //No pi0 decay
  hpythia.readString("23:mayDecay = off");  //No Z0 decay
  hpythia.readString("ProcessLevel:all = off");

  hpythia.init();
}

void do_lund_vac(vector<Parton> partons, vector<Hadron> &vhadrons, int hadro_type)
{
  Event& event      = hpythia.event;
  ParticleData& pdt = hpythia.particleData;

  //Hadronize vacuum
  event.reset();
  int colsum=0;
  for (unsigned int i = 0; i < partons.size(); i++) {
    //If final, introduce in PYTHIA
    if (partons[i].GetD1()==-1) {
      int ide=partons[i].GetId();
      vector<double> p=partons[i].vGetP();
      double px=p[0];
      double py=p[1];
      double pz=p[2];
      double mm=pdt.m0(ide);
      int col=partons[i].GetCol();
      int acol=partons[i].GetAcol();
      double ee=sqrtpos(px*px+py*py+pz*pz+mm*mm);
      event.append(int(ide),23,int(col),int(acol),px,py,pz,ee,mm);
      colsum+=col-acol;
    }
  }
  if (colsum!=0.) cout << " Sumadecolor= " << colsum << endl;
  if (!hpythia.next()) {
    cout << " Event generation aborted prematurely, owing to error" << endl;
  }

    for (int i=0; i<hpythia.event.size(); ++i) {
    if (hpythia.event[i].isFinal()) {
      //Fill vhadrons vector
      vector<double> p;
      for (unsigned int j=1; j<4; j++) p.push_back(hpythia.event[i].p()[j]);
      p.push_back(hpythia.event[i].p()[0]);
      //cout <<"Label= " << i << " Id = " << hpythia.event[i].id() << " t prod= " << hpythia.event[i].tProd() << endl;
      string had_type="lund";
      double phasespace_distance = 0.;
      if ((abs(hpythia.event[i].id())>=400 && abs(hpythia.event[i].id())<=600) || 
      (abs(hpythia.event[i].id())>=4000 && abs(hpythia.event[i].id())<=6000)) {
        //cout <<"Label= " << i << " Id = " << hpythia.event[i].id() << " t prod= " << hpythia.event[i].tProd() << endl;
        int m1 = hpythia.event[i].mother1();
        had_type = "promptD";
        phasespace_distance = phaseSpaceDistance(hpythia.event, i);
        while (true) {
          if (hpythia.event[m1].status()<=-71 && hpythia.event[m1].status()>=-79) break;
          if (hpythia.event[m1].status()==-23) break;
          int idm1 = abs(hpythia.event[m1].id());
          if (idm1==0) break;
          if (
            (idm1 >= 500 && idm1 <= 600) ||
            (idm1 >= 5000 && idm1 <= 6000)
            )
          {
            had_type = "non-promptD";
            break;
          }
          m1 = hpythia.event[m1].mother1();
        }
      }
      vhadrons.push_back( Hadron ( Parton ( p, 0., hpythia.event[i].m(), 0, -1, -1, hpythia.event[i].id(), had_type, 0, 0, true, false), hpythia.event[i].charge(), -1, phasespace_distance) );
    }
  }
}

bool do_lund_med(vector<Quench> quenched, vector<Hadron> &qhadrons, int hadro_type)
{
  Event& event      = hpythia.event;
  ParticleData& pdt = hpythia.particleData;

  // Check if event had a final charm or bottom cross the freezout surface
  bool crosser_event = false;
  for (unsigned int i = 0; i < quenched.size(); i++) {
    if (quenched[i].GetD1()==-1) {
      int ide=quenched[i].GetId();
      if (abs(ide)==4 || abs(ide)==5) {
        if (quenched[i].GetFreezoutCrosser()) crosser_event = true;
        break;
      }
    }
  } //SHOULD WE RATHER CHECK THAT ALL CROSSED THE SURFACE?

  //Hadronize medium
  event.reset();

  if (hadro_type==0) {
    int colsum=0;
    for (unsigned int i = 0; i < quenched.size(); i++) {
      //If final, introduce in PYTHIA
      if (quenched[i].GetD1()==-1) {
        int ide=quenched[i].GetId();
        vector<double> p=quenched[i].vGetP();
        double px=p[0];
        double py=p[1];
        double pz=p[2];
        double mm=pdt.m0(ide);
        int col=quenched[i].GetCol();
        int acol=quenched[i].GetAcol();
        double ee=sqrtpos(px*px+py*py+pz*pz+mm*mm);
        event.append(int(ide),23,int(col),int(acol),px,py,pz,ee,mm);
        colsum+=col-acol;
      }
    }
    if (colsum!=0.) cout << " Sumadecolor= " << colsum << endl;
    if (!hpythia.next()) {
      cout << " Event generation aborted prematurely, owing to error" << endl;
    }

  }
  else if (hadro_type == 1)
  {
    double rempx=0.2;
    double rempy=0.2;
    double p_fake=2500.;
    double rempz=p_fake;
    double reme=std::sqrt(std::pow(rempx,2.)+std::pow(rempy,2.)+std::pow(rempz,2.));

    vector<Quench> pIn, uncolored;
    for (unsigned int i=0; i<quenched.size(); i++) {
      if (quenched[i].GetD1() == -1 && quenched[i].vGetP()[3]!=0. && quenched[i].GetOrig()!="rem") {
        if (quenched[i].GetId()==21 || abs(quenched[i].GetId())<=6) pIn.push_back(quenched[i]);
	      else uncolored.push_back(quenched[i]);
      }
    } 
    
    int col[pIn.size()+2], acol[pIn.size()+2], isdone[pIn.size()+2];
    memset( col, 0, (pIn.size()+2)*sizeof(int) ), memset( acol, 0, (pIn.size()+2)*sizeof(int) ), memset( isdone, 0, (pIn.size()+2)*sizeof(int) );
  
    cout << " pin size= " << pIn.size() << endl;

    // Find number of quarks
    int nquarks=0;
    int isquark[pIn.size()+2];
    memset( isquark, 0, (pIn.size()+2)*sizeof(int) );
    for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {
      if (abs(pIn[ipart].GetId())<=6) {
          isquark[nquarks]=ipart;
          nquarks+=1;  
      }
    }
  
    // Find number of strings
    int nstrings=max(int(double(nquarks)/2.+0.6),1);
    
    // If there are no quarks, need to attach two of them
    int istring=0;
    int one_end[nstrings], two_end[nstrings];
    if (nquarks==0) { // Only attach remnants if event is not empty
      // First quark
	    vector<double> pfirstq={rempx,rempy,rempz,reme};
      pIn.push_back ( Parton ( pfirstq, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true, false ) );
      isquark[nquarks]=pIn.size()-1;
      nquarks+=1;
      isdone[pIn.size()-1]=1;
      one_end[0]=pIn.size()-1;
      // Second quark
      vector<double> psecq={rempx,rempy,-rempz,reme};
      pIn.push_back ( Parton ( psecq, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true, false ) );
      isquark[nquarks]=pIn.size()-1;
      nquarks+=1;
      isdone[pIn.size()-1]=1;
      two_end[istring]=pIn.size()-1;
    }

    // Assign ends of strings (order matters in this algo)
    for(int iquark=0; iquark<nquarks; iquark++) {
      if (isdone[isquark[iquark]]==0) {
        isdone[isquark[iquark]]=1;
        one_end[istring]=isquark[iquark];
        double min_delR=1000000.;
        int partner=-2;
        for(int jquark=0; jquark<nquarks; jquark++) {  
          if (iquark==jquark) continue;
          int d_jquark=isquark[jquark];
          if (isdone[d_jquark]==0) {
            double delR = pIn[isquark[iquark]].delta_R(pIn[d_jquark]);
            if (delR<min_delR) min_delR=delR, partner=jquark;
          }
        }
        if (partner!=-2) {
          isdone[isquark[partner]]=1;
          two_end[istring]=isquark[partner];
          istring+=1;
        }
        else {
          vector<double> pfirstq={rempx,rempy,rempz,reme};
          pIn.push_back ( Parton ( pfirstq, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true, false ) );
          isquark[nquarks]=pIn.size()-1;
          nquarks+=1;
          isdone[pIn.size()-1]=1;
          two_end[istring]=pIn.size()-1;
          cout << "Attached quark remnant flying down +Pz beam" << endl;
        }
      }
    }

    // Assign gluons to a certain string
    int my_string[pIn.size()];
    memset( my_string, 0, pIn.size()*sizeof(int) );
    for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
    {
      if (pIn[ipart].GetId()==21) {
        double min_delR=100000.;
        for (int ns=0; ns<nstrings; ns++)
        {         
          int fq=one_end[ns];
          int sq=two_end[ns];
          double f_delR = pIn[ipart].delta_R(pIn[fq]);
          double s_delR = pIn[ipart].delta_R(pIn[sq]);
          double delR=(f_delR+s_delR)/2.;
          if (delR<min_delR) my_string[ipart]=ns, min_delR=delR;
        }
      }
    }

    // Build up chain using gluons assigned to each string, in a closest pair order
    int lab_col=102;
    for (int ns=0; ns<nstrings; ns++)
    {
      // First end
      int tquark=one_end[ns];
      if (pIn[tquark].GetId()>0) col[tquark]=lab_col;
      else acol[tquark]=lab_col;
      lab_col+=1;
      int link=tquark;

      // Feed into PYTHIA
      int ide=pIn[tquark].GetId();
      double px=pIn[tquark].vGetP()[0];
      double py=pIn[tquark].vGetP()[1];
      double pz=pIn[tquark].vGetP()[2];
      double ee=pIn[tquark].vGetP()[3];
      double mm=pdt.m0(int(ide));
      ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
      if (col[tquark]==0 && acol[tquark]==0 && (ide==21 || abs(ide)<=6)) {
        cout<<"Stopping because of colorless parton trying to be introduced in PYTHIA string";
        exit(0);
      }
      event.append(int(ide),23,col[tquark],acol[tquark],px,py,pz,ee,mm);
        
      int changes=1;
      do {
        changes=0;
        double min_delR=100000.;
        int next_link=0;
        for(unsigned int ipart=0; ipart <  pIn.size(); ++ipart)
        {
          if (pIn[ipart].GetId()==21 && isdone[ipart]==0 && my_string[ipart]==ns)
          {
            changes=1;
            double delR = pIn[link].delta_R(pIn[ipart]);
            if (delR<min_delR) min_delR=delR, next_link=ipart;
          }
        }
        if (changes==1)
        {
          isdone[next_link]=1;
          if (col[link]==lab_col-1) col[next_link]=lab_col, acol[next_link]=lab_col-1;
          else col[next_link]=lab_col-1, acol[next_link]=lab_col;
          lab_col+=1;
          link=next_link;
		
	        // Feed into PYTHIA
          ide=pIn[next_link].GetId();
          px=pIn[next_link].vGetP()[0];
          py=pIn[next_link].vGetP()[1];
          pz=pIn[next_link].vGetP()[2];
          mm=pdt.m0(int(ide));
          ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
          if (col[next_link]==0 && acol[next_link]==0 && (ide==21 || abs(ide)<=6)) {
            cout<<"Stopping because of colorless parton trying to be introduced in PYTHIA string";
            exit(0);
          }
          event.append(int(ide),23,col[next_link],acol[next_link],px,py,pz,ee,mm);
	      }
      } while (changes==1);

      // Attach second end
      if (col[link]==lab_col-1) col[two_end[ns]]=0, acol[two_end[ns]]=lab_col-1;
      else col[two_end[ns]]=lab_col-1, acol[two_end[ns]]=0;
	        
      if (col[two_end[ns]]!=0) { 
        if (pIn[two_end[ns]].GetId()<0) pIn[two_end[ns]].SetId(-pIn[two_end[ns]].GetId());
      } else {
        if (pIn[two_end[ns]].GetId()>0) pIn[two_end[ns]].SetId(-pIn[two_end[ns]].GetId());
      }

      // Feed into PYTHIA
      ide=pIn[two_end[ns]].GetId();
      px=pIn[two_end[ns]].vGetP()[0];
      py=pIn[two_end[ns]].vGetP()[1];
      pz=pIn[two_end[ns]].vGetP()[2];
      mm=pdt.m0(int(ide));
      ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
      if (col[two_end[ns]]==0 && acol[two_end[ns]]==0 && (ide==21 || abs(ide)<=6)) {
        cout<<"Stopping because of colorless parton trying to be introduced in PYTHIA string";
        exit(0);
      }
      event.append(int(ide),23,col[two_end[ns]],acol[two_end[ns]],px,py,pz,ee,mm);
    }
   
    //Introduce uncolored
    for (unsigned int i=0; i<uncolored.size(); i++) {
      int ide=uncolored[i].GetId();
      double px=uncolored[i].vGetP()[0];
      double py=uncolored[i].vGetP()[1];
      double pz=uncolored[i].vGetP()[2];
      double mm=pdt.m0(int(ide));
      double ee=std::sqrt(px*px+py*py+pz*pz+mm*mm);
      event.append(int(ide),23,0,0,px,py,pz,ee,mm);
    }
    
    if (!hpythia.next()) {
      cout << " GOT ISSUES " << endl;
      return 0;
    }

  }

  for (int i=0; i<hpythia.event.size(); ++i) {
    if (hpythia.event[i].isFinal()) {
      //Fill qhadrons vector
      vector<double> p;
      for (unsigned int j=1; j<4; j++) p.push_back(hpythia.event[i].p()[j]);
      p.push_back(hpythia.event[i].p()[0]);
      string had_type="lund";
      double phasespace_distance = 0.;
      if ((abs(hpythia.event[i].id())>=400 && abs(hpythia.event[i].id())<=600) ||
          (abs(hpythia.event[i].id())>=4000 && abs(hpythia.event[i].id())<=6000)) {
        //cout <<"Label= " << i << " Id = " << hpythia.event[i].id() << " t prod= " << hpythia.event[i].tProd() << endl;
        int m1 = hpythia.event[i].mother1();
        had_type = "promptD";
        phasespace_distance = phaseSpaceDistance(hpythia.event, i);
	      while (true) {
          if (hpythia.event[m1].status()<=-71 && hpythia.event[m1].status()>=-79) break;
	        if (hpythia.event[m1].status()==-23) break;
          int idm1 = abs(hpythia.event[m1].id());
          if (idm1==0) break;
	        if (
            (idm1 >= 500 && idm1 <= 600) ||
            (idm1 >= 5000 && idm1 <= 6000)
	        )
	        {
            had_type = "non-promptD";
            break;
          }
	        m1 = hpythia.event[m1].mother1();
        }
      }
      // if (!(abs(hpythia.event[i].id())>=21 && abs(hpythia.event[i].id())<=40)) { // Gauge boson
      //   phasespace_distance = phaseSpaceDistance(hpythia.event, i);
      // }
      qhadrons.push_back( Hadron ( Parton ( p, 0., hpythia.event[i].m(), 0, -1, -1, hpythia.event[i].id(), had_type, 0, 0, true, crosser_event), hpythia.event[i].charge(), -1, phasespace_distance) );
    }
  }

  return 1;

}


// Helper: Check if a particle is a quark or diquark using PDG id ranges.
bool isQuarkOrDiquark(const Particle& p) {
  int absId = std::abs(p.id());
  // Accept light quarks (1-6) or diquarks (1103-5554).
  return ((absId >= 1 && absId <= 6) || (absId >= 1103 && absId <= 5554));
}

// Recursively climb until a fragmentation vertex is found,
// i.e. where the particle has two valid and distinct mothers.
int findFragmentationVertex(const Event& event, int index) {
  if (index <= 0) return index;  // invalid index
  int m1 = event[index].mother1();
  int m2 = event[index].mother2();
  if (m1 > 0 && m2 > 0 && m1 != m2) {
      return index;
  } else if (m1 > 0) {
      return findFragmentationVertex(event, m1);
  }
  cout << "Warning: No valid fragmentation vertex found for index " << index << endl;
  cout << " Ids: " << event[index].id() << endl;
  cout << " Statuses: " << event[index].status() << endl;
  cout << " Mothers: " << event[index].mother1() << " " << event[index].mother2() << endl;
  std::cerr << "Error: No valid fragmentation vertex found for index " << index << std::endl;
  return index;
}

// Climb up a single branch choosing the branch that leads to a quark/diquark endpoint.
// This function climbs until an "initial" particle is reached (both mothers equal to 0).
int climbToInitialQuark(const Event& event, int index) {
  while (true) {
      int m1 = event[index].mother1();
      int m2 = event[index].mother2();
      // Check if this particle is initial.
      if (m1 == 0 && m2 == 0)
          break;
      
      bool m1IsQuark = (m1 > 0 && isQuarkOrDiquark(event[m1]));
      bool m2IsQuark = (m2 > 0 && isQuarkOrDiquark(event[m2]));
      
      if (m1IsQuark && !m2IsQuark) {
          index = m1;
      } else if (m2IsQuark && !m1IsQuark) {
          index = m2;
      } else if (m1IsQuark && m2IsQuark) {
          if (m1==m2){
            index = m1;
          } else {
            cout << "Warning: m1 and m2 are both quarks but not equal" << endl;
            cout << " Ids: " << event[m1].id() << " " << event[m2].id() << endl;
            exit(1);
          }
      } else {
          cout << "Warning: m1 and m2 are not quarks" << endl;
          cout << " Ids: " << event[m1].id() << " " << event[m2].id() << endl;
          cout << " Statuses: " << event[m1].status() << " " << event[m2].status() << endl;
          cout << " Mothers: " << event[m1].mother1() << " " << event[m1].mother2() << " " << event[m2].mother1() << " " << event[m2].mother2() << endl;
          
          event.list();

          exit(1);
      }
  }
  return index;
}

// Given the event and index of a particle (e.g. a hadron),
// this function finds the two string endpoints (initial quarks/diquarks)
// and calculates the phase space distance between them.
double phaseSpaceDistance(const Event& event, int index) {
  // 1. Find the fragmentation vertex.
  int fragIndex = findFragmentationVertex(event, index);
  if (fragIndex <= 0) {
      cout << "Warning: Invalid fragmentation vertex for index " << index << endl;
      cout << " Ids: " << event[index].id() << endl;
      cout << " Statuses: " << event[index].status() << endl;
      cout << " Mothers: " << event[index].mother1() << " " << event[index].mother2() << endl;
      std::cerr << "Error: Invalid fragmentation vertex for index " << index << std::endl;
      return 0.0;
  }
  
  int m1 = event[fragIndex].mother1();
  int m2 = event[fragIndex].mother2();
  if (m1 <= 0 || m2 <= 0) {
      cout << "Warning: Fragmentation vertex does not have two valid mothers." << endl;
      cout << " Ids: " << event[m1].id() << " " << event[m2].id() << endl;
      cout << " Statuses: " << event[m1].status() << " " << event[m2].status() << endl;
      cout << " Mothers: " << event[m1].mother1() << " " << event[m1].mother2() << " " << event[m2].mother1() << " " << event[m2].mother2() << endl;
      cout << " Orig ID: " << event[index].id() << endl;
      cout << " Fragmentation vertex ID: " << event[fragIndex].id() << endl;
      cout << " Fragmentation vertex status: " << event[fragIndex].status() << endl;
      cout << " Fragmentation vertex mothers: " << event[fragIndex].mother1() << " " << event[fragIndex].mother2() << endl;
      std::cerr << "Error: Fragmentation vertex does not have two valid mothers." << std::endl;
      return 0.0;
  }
  
  // 2. Climb each branch to get the quark/diquark endpoints.
  // int endpoint1 = climbToInitialQuark(event, m1);
  // int endpoint2 = climbToInitialQuark(event, m2);
  int endpoint1 = m1;
  int endpoint2 = m2;
  
  // 3. Compute the Minkowski phase space distance:
  //    Δs² = (E1 - E2)² - (px1 - px2)² - (py1 - py2)² - (pz1 - pz2)².
  // double dE  = event[endpoint1].e()  - event[endpoint2].e();
  double dpx = event[endpoint1].px() - event[endpoint2].px();
  double dpy = event[endpoint1].py() - event[endpoint2].py();
  double dpz = event[endpoint1].pz() - event[endpoint2].pz();
  
  // return dE*dE - dpx*dpx - dpy*dpy - dpz*dpz;
  return sqrt( dpx*dpx + dpy*dpy + dpz*dpz );
}