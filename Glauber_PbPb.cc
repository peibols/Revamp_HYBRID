#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <assert.h>

#include "Random.h"

using std::vector;
using namespace std;

double glaub[200][200];
int g_maxx=200;
int g_maxy=200;
double g_deltax=0.098650;
double g_deltay=0.098650;

double TA[4000], step;
double norm=1.;

double bmin, bmax;

//LHC
//centrality: 		0 5 10 20 30 40 50 60 70
//impact parameter: 	0 3.5 4.94 6.98 8.55 9.88 11.04 12.09 13.05

//RHC
//centrality:		0 5 10 20 30 40 50 60 70
//impact parameter:     0 3.5 4.7 6.7 8.2 9.4 10.6 11.6 12.5

void read_nuclear(int nhyd, std::string cent);
void gxy(double &x, double &y, numrand &nr);
double gTAA(double x, double y, double b);

//nhyd obsolete, used when using different hydro profiles
void read_nuclear(int nhyd, std::string cent)
{

  char glauFile[200];
  sprintf(glauFile,"./TAb2LL.dat");
  ifstream initial (glauFile);

  assert(!initial.fail());

  cout << " Reading Initial Energy Density..." << endl;

  if (cent=="0-5") bmin=0., bmax=3.5;
  else if (cent=="5-10") bmin=3.5, bmax=4.94;
  else if (cent=="10-20") bmin=4.94, bmax=6.98;
  else if (cent=="20-30") bmin=6.98, bmax=8.55;
  else if (cent=="30-40") bmin=8.55, bmax=9.88;
  else if (cent=="40-50") bmin=9.88, bmax=11.04;
  else if (cent=="50-60") bmin=11.04, bmax=12.09;
  else if (cent=="60-70") bmin=12.09, bmax=13.05;
  else { cout << " Unrecognized cent= " << cent.c_str() << endl; exit(0); }
  cout << " Bmin= " << bmin << " Bmax= " << bmax << endl;
	
  //Reading nuclear overlap function as a function of impact parameter
  double b2;
  for (unsigned a=0; a<4000; a++)
  {       
    initial >> b2 >> TA[a];
    if (a == 1) step = b2;
  }	

}

//------------------gxy-----------------------------//
void gxy(double &x, double &y, numrand &nr) {

  double rho,phi;
  double P;
  double b;

  naiguels:
  b=sqrt((bmax*bmax-bmin*bmin)*nr.rando()+bmin*bmin);
  norm=1.;
  norm=gTAA(0.,0.,bmin);
  rho=sqrt(150.*nr.rando());
  phi=2.*3.141592654*nr.rando();
  x=rho*cos(phi);
  y=rho*sin(phi);
  P=nr.rando();
  if(P>gTAA(x,y,b)) goto naiguels;

  return;
}

//-----------------gTAA------------------------------//
double gTAA(double x, double y, double b)
{
  int il, irr;
  double rho2, use;

  rho2=pow(x+b/2.,2.)+y*y;
  il=int(rho2/step);
  rho2=pow(x-b/2.,2.)+y*y;
  irr=int(rho2/step);
  use=0.;
  if(il<4000 && irr<4000) {
    use=TA[il]*TA[irr]/norm;
  }

  return use;
}
