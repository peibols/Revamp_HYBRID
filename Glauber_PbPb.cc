#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <assert.h>

#include "Random.h"

#include "global.h"

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

//0-5% LHC and RHIC
//double bmin=0.;
//double bmax=3.5;

//50-60% RHIC
//double bmin=10.;
//double bmax=11.6;

//double bmin=3.5;
//double bmax=4.7;

//LHC
//centrality: 		0 5 10 20 30 40 50 60 70
//impact parameter: 	0 3.5 4.94 6.98 8.55 9.88 11.04 12.09 13.05

//RHC
//centrality:		0 5 10 20 30 40 50 60 70
//impact parameter:     0 3.5 4.7 6.7 8.2 9.4 10.6 11.6 12.5
//50-60% LHC
//double bmin=11.04;
//double bmax=12.09;

void read_nuclear(int nhyd, std::string cent);
void gxy(double &x, double &y, numrand &nr);
//double gGlaub(double x, double y);
double gTAA(double x, double y, double b);

void read_nuclear(int nhyd, std::string cent)
{
	char glauFile[200];
        //sprintf(glauFile,"/gs/project/cqn-654-ad/mayank/events_for_jets/PbPb_%s_2p76/job-%i/results/u_field_1.dat",cent.c_str(),nhyd);
        //sprintf(glauFile,"/home/peibols/projects/rrg-jeon-ac/group_writable/hydro_from_Mayank/PbPb_%s_2p76/job-%i/results/u_field_1.dat",cent.c_str(),nhyd);
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
	
	//Reading T(t)
        double b2;
        for (unsigned a=0; a<4000; a++)
        {       
                initial >> b2 >> TA[a];
                if (a == 1) step = b2;
        }	
/*
	string s, sub;
	//First line: some crap
	getline(initial,s);
        //Rest of lines: crap, x, y, edensity(MeV/fm^3), 7 x crap
	double crap, x, y, edens;
	double maxedens=0.;
	int ix=0, iy=0;
	do {
		getline(initial,s);
		if (s.empty()) continue;
		istringstream iss(s);
		iss >> crap >> x >> y >> edens;
		if (edens>maxedens) maxedens=edens;
		//cout << " x= " << x << " y= " << y << " edens= " << edens << endl;
		glaub[ix][iy]=edens;
		iy+=1;
		if (iy==g_maxy) ix+=1, iy=0;	
	} while(ix<g_maxx);
	cout << " Finish Reading Initial Energy Density." << endl;
	cout << " Max Energy Density = " << maxedens << endl;
	//Set Maximum to 1
	for (int i=0; i<g_maxx; i++)
	{
		for (int j=0; j<g_maxy; j++)
		{
			//double xt = -9.855135e+00+g_deltax*double(i);
			//double yt = -9.855135e+00+g_deltay*double(j);
			//cout << xt << " " << yt << " " << gGlaub(xt,yt) << endl;
			glaub[i][j]/=maxedens;
		}
	}
*/

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
        if (il >= 4000) il = 3999;
        if (irr >= 4000) irr = 3999;
        use=TA[il]*TA[irr]/TA[0];
        return use;
}

/*
void gxy(double &x, double &y, numrand &nr) {

        double rho,phi;
	double P;

        naiguels:
        rho=sqrt(150.*nr.rando());
        phi=2.*3.141592654*nr.rando();
        x=rho*cos(phi);
        y=rho*sin(phi);
	P=nr.rando();
        if(P>gGlaub(x,y)) goto naiguels;
}

double gGlaub(double x, double y)
{
	double gdens=0.;
	
	int ix, dx, iy, dy;
	
	if (x>=0.) {
		ix = int(x/g_deltax)+(g_maxx)/2;
		dx = (x - double(ix-(g_maxx)/2)*g_deltax)/g_deltax;
	}
	else {
		ix = int(x/g_deltax)+(g_maxx)/2-1;
                dx = (x - double(ix-(g_maxx)/2)*g_deltax)/g_deltax;
	}
	
	if (y>=0.) {
                iy = int(y/g_deltay)+(g_maxy)/2;
                dy = (y - double(iy-(g_maxy)/2)*g_deltay)/g_deltay;
        }
        else {
                iy = int(y/g_deltay)+(g_maxy)/2-1;
                dy = (y - double(iy-(g_maxy)/2)*g_deltay)/g_deltay;
        }

	if (ix<0 || ix>=g_maxx-1 || iy<0 || iy>=g_maxy-1) return gdens; 
	gdens=glaub[ix][iy]*(1.-dx)*(1.-dy);
	gdens+=glaub[ix][iy+1]*(1.-dx)*dy;
	gdens+=glaub[ix+1][iy]*dx*(1.-dy);
	gdens+=glaub[ix+1][iy+1]*dx*dy;

	if (gdens>1.) cout << " gGlaub not properly normalised: gdens = " << gdens << endl;
	return gdens;
}
*/
