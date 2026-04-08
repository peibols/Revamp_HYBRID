#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <assert.h>
#include <stdlib.h>

#include "Random.h"

using std::vector;
using namespace std;

double glaub[401][401];
double g_deltax=0.05;
double g_deltay=0.05;
int g_maxx=401;
int g_maxy=401;

void read_nuclear(int, std::string);
void gxy(double &x, double &y, numrand &nr);
double gGlaub(double x, double y);

void read_nuclear(int nhyd, std::string cent)
{
	char glauFile[100];
	//if (cent!="0-100") sprintf(glauFile,"/gs/project/cqn-654-ad/cshen/MUSIC_pA_3+1d/pPb/initial_conditions_MCGlbpPb5020C%s_Chanwook/rho_binary_event_%i_block.dat",cent.c_str(),nhyd);	
	if (cent!="0-100") sprintf(glauFile,"/home/peibols/projects/rrg-jeon-ac/group_writable/hydro_from_Chun/initial_conditions_MCGlbpPb5020C%s_Chanwook/rho_binary_event_%i_block.dat",cent.c_str(),nhyd);
	else {
		cout << " Not ready for pPb MIN BIAS !!!";
		exit(0);
	}
	//else sprintf(glauFile,"/gs/project/cqn-654-ad/peibols/hybrid_ebe/chun_hydro/initial_conditions_MCGlbpPb5020C0-100_resized/rho_binary_event_%i_block.dat",nhyd);
	ifstream initial (glauFile);

	cout << "glauber= " << glauFile << endl;
	assert(!initial.fail());

	//cout << " Reading Initial Energy Density..." << endl;
	//401x401 Matrix
	string s, sub;
	double crap, x, y, bdens;
	double maxbdens=0.;
	int ix=0, iy=0;
	do {
		getline(initial,s);
		istringstream iss(s);
		do {
			iss >> bdens;
			glaub[ix][iy]=bdens;
			iy+=1;
			if (bdens>maxbdens) maxbdens=bdens;
		} while(iy<g_maxy);
		iy=0;
		ix+=1;
	} while(ix<g_maxx);
	//cout << " Finish Reading Initial Energy Density." << endl;
	//cout << " Max Energy Density = " << maxbdens << endl;
	//Set Maximum to 1
	double ncoll=0.;
	for (int i=0; i<g_maxx; i++)
	{
		for (int j=0; j<g_maxy; j++)
		{
			ncoll+=glaub[i][j]*g_deltax*g_deltay;
			glaub[i][j]/=maxbdens;
		}
	}
	//cout << " Ncoll= " << ncoll << endl;
}

void gxy(double &x, double &y, numrand &nr) {

        double rho,phi;
        int ix, iy;
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
		ix = int(x/g_deltax)+(g_maxx-1)/2;
		dx = (x - double(ix-(g_maxx-1)/2)*g_deltax)/g_deltax;
	}
	else {
		ix = int(x/g_deltax)+(g_maxx-1)/2-1;
                dx = (x - double(ix-(g_maxx-1)/2)*g_deltax)/g_deltax;
	}
	
	if (y>=0.) {
                iy = int(y/g_deltay)+(g_maxy-1)/2;
                dy = (y - double(iy-(g_maxy-1)/2)*g_deltay)/g_deltay;
        }
        else {
                iy = int(y/g_deltay)+(g_maxy-1)/2-1;
                dy = (y - double(iy-(g_maxy-1)/2)*g_deltay)/g_deltay;
        }

	if (ix<0 || ix>=g_maxx-1 || iy<0 || iy>=g_maxy-1) return gdens; 
	gdens=glaub[ix][iy]*(1.-dx)*(1.-dy);
	gdens+=glaub[ix][iy+1]*(1.-dx)*dy;
	gdens+=glaub[ix+1][iy]*dx*(1.-dy);
	gdens+=glaub[ix+1][iy+1]*dx*dy;

	if (gdens>1.) cout << " gGlaub not properly normalised: gdens = " << gdens << endl;
	return gdens;
}
