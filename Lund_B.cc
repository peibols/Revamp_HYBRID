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
void init_lund(int seed);

Pythia hpythia;

void init_lund(int seed)
{
  hpythia.readString("Random:setSeed = on");
  ostringstream seedstring;
  seedstring << "Random:seed = " << seed;
  hpythia.readString(seedstring.str().c_str());


  /*
  //Ezra decays
  hpythia.readString("111:mayDecay = on");
  hpythia.readString("310:mayDecay = off");
  hpythia.readString("3122:mayDecay = off");
  hpythia.readString("3112:mayDecay = off");
  hpythia.readString("3222:mayDecay = off");
  hpythia.readString("3312:mayDecay = off");
  hpythia.readString("3322:mayDecay = off");
  hpythia.readString("3334:mayDecay = off");
  */

  /*
  //Decays
  hpythia.readString("411:mayDecay = off");  //D+
  hpythia.readString("421:mayDecay = off");  //D0
  hpythia.readString("431:mayDecay = off");  //D+S
  */

  hpythia.readString("511:mayDecay = off");  //B0
  hpythia.readString("521:mayDecay = off");  //B+
  hpythia.readString("531:mayDecay = off");  //Bs0
  hpythia.readString("541:mayDecay = off");  //Bc+

  hpythia.readString("ProcessLevel:all = off");

  hpythia.init();
  //cout << " LUND INITIALISED \n";
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
      vhadrons.push_back( Hadron ( Parton ( p, 0., hpythia.event[i].m(), 0, -1, -1, hpythia.event[i].id(), "lund", 0, 0, true), hpythia.event[i].charge(), -1 ) );
    }
  }

//End program
}

bool do_lund_med(vector<Quench> quenched, vector<Hadron> &qhadrons, int hadro_type)
{

  Event& event      = hpythia.event;
  ParticleData& pdt = hpythia.particleData;
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
        //if (px*px+py*py+pz*pz==0.) {
          //cout << " Putting 0 momentum particle!" << endl;
          //cout << " energy= " << ee << endl;
        //}
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

    //double Lambda_QCD=0.2;
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
      //  cout << pIn[ipart].vGetP()[0] << " "
	//     << pIn[ipart].vGetP()[1] << " "
	//     << pIn[ipart].vGetP()[2] << " "
	//     << pIn[ipart].vGetP()[3] << endl;

        if (abs(pIn[ipart].GetId())<=6) {
            isquark[nquarks]=ipart;
            nquarks+=1;  
        }
    }

    //cout << " nquarks = " << nquarks << endl;
  
    // Find number of strings
    int nstrings=max(int(double(nquarks)/2.+0.6),1);
    
    // If there are no quarks, need to attach two of them
    int istring=0;
    int one_end[nstrings], two_end[nstrings];
    if (nquarks==0) { // Only attach remnants if event is not empty
        // First quark
	vector<double> pfirstq={rempx,rempy,rempz,reme};
        pIn.push_back ( Parton ( pfirstq, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true ) );
        isquark[nquarks]=pIn.size()-1;
        nquarks+=1;
        isdone[pIn.size()-1]=1;
        one_end[0]=pIn.size()-1;
        // Second quark
	vector<double> psecq={rempx,rempy,-rempz,reme};
        pIn.push_back ( Parton ( psecq, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true ) );
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
		    //cout << " delR = " << delR << endl;
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
                pIn.push_back ( Parton ( pfirstq, 0., 0., 0, -1, -1, 1, "rem", 0, 0, true ) );
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
      //event.reset();

      //cout << " NEXT STRING " << endl;
    
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
      }
      else {
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
      qhadrons.push_back( Hadron ( Parton ( p, 0., hpythia.event[i].m(), 0, -1, -1, hpythia.event[i].id(), "lund", 0, 0, true), hpythia.event[i].charge(), -1 ) );
    }
  }

  return 1;

//End program
}
