#include "Pythia8/Pythia.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

using namespace std;
using namespace Pythia8;

int main(int argc, char** argv) {

  assert(argc==3);

  int core=atoi(argv[1]);
  int nEvent=atoi(argv[2]);

  Rndm randa;
  randa.init(0);

  //Input Files
  char InpF[100];
  sprintf(InpF,"%s.out",argv[1]);
  ifstream infile;
  infile.open (InpF);

  ostringstream WakeF;
  WakeF << "wake_" << core << ".out";
  std::ifstream infileWake(WakeF.str().c_str());

  ofstream outfile("HYBRID_Hadrons.out");
  ofstream parfile("HYBRID_Partons.out");
 
  Pythia pythia;
  Event& event      = pythia.event;
  ParticleData& pdt = pythia.particleData;

  pythia.readString("Print:quiet = on");
  pythia.readString("111:mayDecay = off");
  pythia.readString("ProcessLevel:all = off");

  pythia.init();

  for (int iEvent = 0; iEvent < nEvent; ++iEvent) {

      event.reset();
      double colsum=0.;
      double mm, ee;
      double px, py, pz, en;
      double qpx, qpy, qpz, qe, col, acol, ide;
      std::string sQx,sQy,sQz,sQe;
      std::string spx;
      double weight, cross, cross_err, pthat;
      double xcre, ycre;

      outfile << "# event " << iEvent << endl;
      parfile << "# event " << iEvent << endl;

      double p_px[1000], p_py[1000], p_pz[1000], p_e[1000], p_m[1000], p_ide[1000];
      int k=0;
      while (infile >> spx >> py >> pz >> en >> sQx >> sQy >> sQz >> sQe >> col >> acol >> ide) {
        
        qpx = strtod(sQx.c_str(), NULL);
        qpy = strtod(sQy.c_str(), NULL);
        qpz = strtod(sQz.c_str(), NULL);
        qe = strtod(sQe.c_str(), NULL);
        if (qpx!=qpx || qpy!=qpy || qpz!=qpz || qe!=qe ) {
          cout << " BUG STILL THEEEEERE\n";
        }
        if (qpx!=qpx) qpx=0.;
        if (qpy!=qpy) qpy=0.;
        if (qpz!=qpz) qpz=0.;
        if (qe!=qe) qe=0.;
        if (spx=="NEXT") {
          xcre=en;
          ycre=qpx;
	  cross=acol;
          weight=ide;
	  cross_err=col;
	  pthat=qpz;
          break;
        }
        else {
          px = strtod(spx.c_str(), NULL);
          if (qe==0.) {
            qpx=px/100000.;
            qpy=py/100000.;
            qpz=pz/100000.;
          }
	  else {
            p_px[k]=qpx;
            p_py[k]=qpy;
            p_pz[k]=qpz;
            if (ide>-1000) mm=pdt.m0(int(ide));
	    else mm=0.;
            ee=sqrtpos(qpx*qpx+qpy*qpy+qpz*qpz+mm*mm);
	    p_m[k]=mm;
	    p_ide[k]=ide;
            p_e[k]=ee;
	    k++;
	    if (ide<=-1000) continue;
	  }

          mm=pdt.m0(int(ide));
          ee=sqrtpos(qpx*qpx+qpy*qpy+qpz*qpz+mm*mm);
          event.append(int(ide),23,int(col),int(acol),qpx,qpy,qpz,ee,mm);
          colsum+=col-acol;

        }

      }

      if (colsum!=0.) cout << " Sumadecolor= " << colsum << " at count= " << iEvent << " ";
      if (!pythia.next()) {
        cout << " Event generation aborted prematurely, owing to error at count= " << iEvent << "\n";
        //break;
      }

      outfile << "weight " << weight << " cross " << cross << " X " << xcre << " Y " << ycre <<
              " cross_err " << cross_err << " pthat " << pthat << endl;
      parfile << "weight " << weight << " cross " << cross << " X " << xcre << " Y " << ycre <<
              " cross_err " << cross_err << " pthat " << pthat << endl;

      for (int ik=0; ik<k; ik++) {
        int tide=p_ide[ik];
	int tag;
	if (tide==-1000) tide=1, tag=-2;
	else if (tide==-2000) tide=2, tag=-2;
	else if (tide==-3000) tide=3, tag=-2;
	else tag=0;
        parfile << p_px[ik] << " " << p_py[ik] << " " << p_pz[ik] << " " << p_m[ik] << " " << tide << " " << tag << endl;
        if (tag==-2) outfile << p_px[ik] << " " << p_py[ik] << " " << p_pz[ik] << " " << p_m[ik] << " " << tide << " " << tag << endl;
      }

      double Px, Py, Pz, E, mass;
      int iide, charge;
      for (int i=0; i<pythia.event.size(); ++i) {
        if (pythia.event[i].isFinal()) {
          Px=pythia.event[i].px();
          Py=pythia.event[i].py();
          Pz=pythia.event[i].pz();
          E=pythia.event[i].e();
          mass=sqrt(pow(E,2.)-pow(Px,2.)-pow(Py,2.)-pow(Pz,2.));
          if (mass!=mass) mass=0.;
	  iide=pythia.event[i].id();
          charge=pythia.event[i].charge();
          if (Px!=Px) continue;
          outfile << Px << " " << Py << " " << Pz << " " << mass << " " << iide << " " << 0 << endl;
        }
     }

     int stat;
     while (infileWake >> Px >> Py >> Pz >> E >> stat)
     {
       if (Px==0.000123 && Py==0. && Pz==0. && E==0.) {
          break;
        }
        else {
          double masstwo=pow(E,2.)-pow(Px,2.)-pow(Py,2.)-pow(Pz,2.);
        
	  double randchar=randa.flat();
          if (masstwo<0.5*0.5) {
            if (randchar<=1./3.) charge=-1, iide=-211;
            else if (randchar<=2./3.) charge=0, iide=111;
            else charge=1, iide=211;
          }
          else {
            if (randchar>1./2.) charge=1, iide=2212;
            else charge=-1., iide=-2212;
          }
          int label;
          if (stat==1) label=1;
          else if (stat==-1) label=2;
	  mass=pdt.m0(int(iide));
          outfile << Px << " " << Py << " " << Pz << " " << mass << " " << iide << " " << label << endl;
       }

     }

     outfile << "end" << endl;
     parfile << "end" << endl;

  }

  outfile.close();
  parfile.close();

}

     
