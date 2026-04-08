#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>
#include <stdlib.h>

#include "global.h"

double hydrot[297][64][101][101];
double hydrox[297][64][101][101];
double hydroy[297][64][101][101];
double hydroz[297][64][101][101];

int maxx=101;
int maxy=101;
int maxeta=64;

double deltax=0.2;
double deltay=0.2;
double deltaeta=0.2;
double deltat=0.06;

double tau0=0.6;	//Initial time
double tau1;		//Final time
double eta1=6.4;	//Extreme abs(eta)

double dt, dx, dy, deta;
int it, ix, iy, ieta;

void getGrid(double tau, double x, double y, double eta, double* dt, double* dx, double* dy, double* deta, int* it, int* ix, int* iy, int* ieta);
double gT(double tau, double x, double y, double eta);
double gVx(double tau, double x, double y, double eta);
double gVy(double tau, double x, double y, double eta);
double gVz(double tau, double x, double y, double eta);
void read_hydro(int, std::string);

using std::vector;
using namespace std;

/*
template<typename T>
std::istream & binary_read(std::istream& stream, T& value){
    return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}
*/

void read_hydro(int nhyd, std::string cent)
{
	char hydFile[200];
	//if (cent!="0-100") sprintf(hydFile,"/gs/project/cqn-654-ad/cshen/MUSIC_pA_3+1d/pPb/pPb_C%s-with_evo_for_Chanwook_collected/event_%i/evolution_xyeta.dat",cent.c_str(),nhyd);
	if (cent!="0-100") sprintf(hydFile,"/home/peibols/projects/rrg-jeon-ac/group_writable/hydro_from_Chun/pPb_C%s-with_evo_for_Chanwook_collected/event_%i/evolution_xyeta.dat",cent.c_str(),nhyd);
	else {
		cout << " No hydro yet for pPb MIN BIAS !!";
		exit(0);
	}
	//else sprintf(hydFile,"/gs/project/cqn-654-ad/peibols/hybrid_ebe/chun_hydro/pPb_C0-100_with_evo_for_Chanwook_collected/event_%i/evolution_xyeta.dat",nhyd);
	std::FILE *hydro;
	hydro = std::fopen(hydFile, "rb");

	cout << " hydro= " << hydFile << endl;
	assert(hydro != NULL);

	cout << " Reading Hydro Event " << nhyd << endl;
	clock_t startClock = clock();
	int t=0;
	double T, QGPfrac, vx, vy, vz;
	int size = sizeof(double);
	int do_exit=0;
	do {
		for (int l=0; l<maxeta; l++)
		{
			for (int k=0; k<maxy; k++)
			{
				for (int j=0; j<maxx; j++) {
					int status = 0;
              				status = std::fread(&T, size, 1, hydro);
              				status += std::fread(&QGPfrac, size, 1, hydro);
              				status += std::fread(&vx, size, 1, hydro);
              				status += std::fread(&vy, size, 1, hydro);
              				status += std::fread(&vz, size, 1, hydro);
					if (status!=5) { do_exit=1; break; }
					hydrot[t][l][k][j]=T;
					hydrox[t][l][k][j]=vx;
					hydroy[t][l][k][j]=vy;
					hydroz[t][l][k][j]=vz;
					//if (column[0]==column[4]) cout << " MADAFACKA col0= " << column[0] << "\n";
					//cout << " c0= " << column[0] << " c1= " << column[1] << " c2= " << column[2] << " ";
					//cout << " c3= " << column[3] << " c4= " << column[4] << endl;
					//cout << endl;
				}
				if (do_exit==1) break;
			}
			if (do_exit==1) break;
		}
		//cout << " t= " << t << endl;
		t+=1;
	} while (do_exit==0);
	tau1=tau0+double(t-1)*deltat;
	clock_t endClock = clock();
	cout << " Finish Hydro Read in " << double((endClock - startClock)) / CLOCKS_PER_SEC << " secs. \n";

	/*	
	double step=0.01;
	cout << " T at center= " << hydrot[0][32][50][50] << endl;
	for (unsigned int i=0; i<1000; i++)
	{
		//double x=0.8*double(i)*step;
		double temp=gT(tau0+double(i)*step,double(i)*step,0.,0.);
		double vx=gVx(tau0+double(i)*step,double(i)*step,0.,0.);
		double vy=gVy(tau0+double(i)*step,double(i)*step,0.,0.);
		double vz=gVz(tau0+double(i)*step,double(i)*step,0.,0.);
		double v2=pow(vx,2.)+pow(vy,2.)+pow(vz,2.);
		cout << tau0+double(i)*step << " temp= " << temp << " vx= " << vx << " vy= " << vy << " vz= " << vz << " v2= " << v2 << endl;
		if (temp==0.) break;
	}
	*/
	
}

void getGrid(double tau, double x, double y, double eta, double* udt, double* udx, double* udy, double* udeta, int* uit, int* uix, int* uiy, int* uieta)
{
        *uit=int((tau-tau0)/deltat);
        *udt=(tau-tau0-double(*uit)*deltat)/deltat;

        if (y>=0.) {
                *uiy = int(y/deltay)+(maxy-1)/2;
                *udy = (y - double(*uiy-(maxy-1)/2)*deltay)/deltay;
        }
        else {
                *uiy = int(y/deltay)+(maxy-1)/2-1;
                *udy = (y - double(*uiy-(maxy-1)/2)*deltay)/deltay;
        }

        if (x>=0.) {
                *uix = int(x/deltax)+(maxx-1)/2;
                *udx = (x - double(*uix-(maxx-1)/2)*deltax)/deltax;
        }
        else {
                *uix = int(x/deltax)+(maxx-1)/2-1;
                *udx = (x - double(*uix-(maxx-1)/2)*deltax)/deltax;
        }

	if (eta>=0.) {
                *uieta = int(eta/deltaeta)+maxeta/2;
                *udeta = (eta - double(*uieta-maxeta/2)*deltaeta)/deltaeta;
        }
        else {
        	*uieta = int(eta/deltaeta)+maxeta/2-1;
                *udeta = (eta - double(*uieta-maxeta/2)*deltaeta)/deltaeta;
	}
}

double gT(double tau, double x, double y, double eta)
{
        double gete=0.;

        if (tau>=tau1 || abs(eta)>=eta1 || tau<tau0) {
		//cout << " No Hydro here: tau= " << tau << ", eta= " << eta << endl;
                return gete;
        }

	getGrid(tau,x,y,eta,&dt,&dx,&dy,&deta,&it,&ix,&iy,&ieta);

	if (ix<0 || ix>=maxx-1 || iy<0 || iy>=maxy-1 || ieta<0 || ieta>=maxeta-1) return gete;
        gete=hydrot[it][ieta][iy][ix]*(1.-dt)*(1.-dx)*(1.-dy)*(1.-deta);
        gete+=hydrot[it+1][ieta][iy][ix]*(1.-dx)*(1.-dy)*dt*(1.-deta);
        gete+=hydrot[it][ieta][iy][ix+1]*(1.-dt)*(1.-dy)*dx*(1.-deta);
        gete+=hydrot[it][ieta][iy+1][ix]*(1.-dx)*(1.-dt)*dy*(1.-deta);
        gete+=hydrot[it+1][ieta][iy][ix+1]*dx*(1.-dy)*dt*(1.-deta);
        gete+=hydrot[it][ieta][iy+1][ix+1]*(1.-dt)*dy*dx*(1.-deta);
        gete+=hydrot[it+1][ieta][iy+1][ix]*(1.-dx)*dt*dy*(1.-deta);
        gete+=hydrot[it+1][ieta][iy+1][ix+1]*dx*dt*dy*(1.-deta);
	gete+=hydrot[it][ieta+1][iy][ix]*(1.-dt)*(1.-dx)*(1.-dy)*deta;
	gete+=hydrot[it+1][ieta+1][iy][ix]*(1.-dx)*(1.-dy)*dt*deta;
        gete+=hydrot[it][ieta+1][iy][ix+1]*(1.-dt)*(1.-dy)*dx*deta;
        gete+=hydrot[it][ieta+1][iy+1][ix]*(1.-dx)*(1.-dt)*dy*deta;
        gete+=hydrot[it+1][ieta+1][iy][ix+1]*dx*(1.-dy)*dt*deta;
        gete+=hydrot[it][ieta+1][iy+1][ix+1]*(1.-dt)*dy*dx*deta;
        gete+=hydrot[it+1][ieta+1][iy+1][ix]*(1.-dx)*dt*dy*deta;
        gete+=hydrot[it+1][ieta+1][iy+1][ix+1]*dx*dt*dy*deta;

	return gete;
}

double gVx(double tau, double x, double y, double eta)
{
        double gete=0.;

        if (tau>=tau1 || abs(eta)>=eta1 || tau<tau0) {
                //cout << " No Hydro here: tau= " << tau << ", eta= " << eta << endl;
                return gete;
        }

        getGrid(tau,x,y,eta,&dt,&dx,&dy,&deta,&it,&ix,&iy,&ieta);

        if (ix<0 || ix>=maxx-1 || iy<0 || iy>=maxy-1 || ieta<0 || ieta>=maxeta-1) return gete;
        gete=hydrox[it][ieta][iy][ix]*(1.-dt)*(1.-dx)*(1.-dy)*(1.-deta);
        gete+=hydrox[it+1][ieta][iy][ix]*(1.-dx)*(1.-dy)*dt*(1.-deta);
        gete+=hydrox[it][ieta][iy][ix+1]*(1.-dt)*(1.-dy)*dx*(1.-deta);
        gete+=hydrox[it][ieta][iy+1][ix]*(1.-dx)*(1.-dt)*dy*(1.-deta);
        gete+=hydrox[it+1][ieta][iy][ix+1]*dx*(1.-dy)*dt*(1.-deta);
        gete+=hydrox[it][ieta][iy+1][ix+1]*(1.-dt)*dy*dx*(1.-deta);
        gete+=hydrox[it+1][ieta][iy+1][ix]*(1.-dx)*dt*dy*(1.-deta);
        gete+=hydrox[it+1][ieta][iy+1][ix+1]*dx*dt*dy*(1.-deta);
        gete+=hydrox[it][ieta+1][iy][ix]*(1.-dt)*(1.-dx)*(1.-dy)*deta;
        gete+=hydrox[it+1][ieta+1][iy][ix]*(1.-dx)*(1.-dy)*dt*deta;
        gete+=hydrox[it][ieta+1][iy][ix+1]*(1.-dt)*(1.-dy)*dx*deta;
        gete+=hydrox[it][ieta+1][iy+1][ix]*(1.-dx)*(1.-dt)*dy*deta;
        gete+=hydrox[it+1][ieta+1][iy][ix+1]*dx*(1.-dy)*dt*deta;
        gete+=hydrox[it][ieta+1][iy+1][ix+1]*(1.-dt)*dy*dx*deta;
        gete+=hydrox[it+1][ieta+1][iy+1][ix]*(1.-dx)*dt*dy*deta;
        gete+=hydrox[it+1][ieta+1][iy+1][ix+1]*dx*dt*dy*deta;

        return gete;
}

double gVy(double tau, double x, double y, double eta)
{
        double gete=0.;

        if (tau>=tau1 || abs(eta)>=eta1 || tau<tau0) {
                //cout << " No Hydro here: tau= " << tau << ", eta= " << eta << endl;
                return gete;
        }

	//cout << " tau= " << tau << " eta= " << eta << " x= " << x << " y= " << y << endl;
        getGrid(tau,x,y,eta,&dt,&dx,&dy,&deta,&it,&ix,&iy,&ieta);

        if (ix<0 || ix>=maxx-1 || iy<0 || iy>=maxy-1 || ieta<0 || ieta>=maxeta-1) return gete;
        gete=hydroy[it][ieta][iy][ix]*(1.-dt)*(1.-dx)*(1.-dy)*(1.-deta);
        gete+=hydroy[it+1][ieta][iy][ix]*(1.-dx)*(1.-dy)*dt*(1.-deta);
        gete+=hydroy[it][ieta][iy][ix+1]*(1.-dt)*(1.-dy)*dx*(1.-deta);
        gete+=hydroy[it][ieta][iy+1][ix]*(1.-dx)*(1.-dt)*dy*(1.-deta);
        gete+=hydroy[it+1][ieta][iy][ix+1]*dx*(1.-dy)*dt*(1.-deta);
        gete+=hydroy[it][ieta][iy+1][ix+1]*(1.-dt)*dy*dx*(1.-deta);
        gete+=hydroy[it+1][ieta][iy+1][ix]*(1.-dx)*dt*dy*(1.-deta);
        gete+=hydroy[it+1][ieta][iy+1][ix+1]*dx*dt*dy*(1.-deta);
        gete+=hydroy[it][ieta+1][iy][ix]*(1.-dt)*(1.-dx)*(1.-dy)*deta;
        gete+=hydroy[it+1][ieta+1][iy][ix]*(1.-dx)*(1.-dy)*dt*deta;
        gete+=hydroy[it][ieta+1][iy][ix+1]*(1.-dt)*(1.-dy)*dx*deta;
        gete+=hydroy[it][ieta+1][iy+1][ix]*(1.-dx)*(1.-dt)*dy*deta;
        gete+=hydroy[it+1][ieta+1][iy][ix+1]*dx*(1.-dy)*dt*deta;
        gete+=hydroy[it][ieta+1][iy+1][ix+1]*(1.-dt)*dy*dx*deta;
        gete+=hydroy[it+1][ieta+1][iy+1][ix]*(1.-dx)*dt*dy*deta;
        gete+=hydroy[it+1][ieta+1][iy+1][ix+1]*dx*dt*dy*deta;

        return gete;
}

double gVz(double tau, double x, double y, double eta)
{
        double gete=0.;

        if (tau>=tau1 || abs(eta)>=eta1 || tau<tau0) {
                //cout << " No Hydro here: tau= " << tau << ", eta= " << eta << endl;
                return gete;
        }

        getGrid(tau,x,y,eta,&dt,&dx,&dy,&deta,&it,&ix,&iy,&ieta);

        if (ix<0 || ix>=maxx-1 || iy<0 || iy>=maxy-1 || ieta<0 || ieta>=maxeta-1) return gete;
        gete=hydroz[it][ieta][iy][ix]*(1.-dt)*(1.-dx)*(1.-dy)*(1.-deta);
        gete+=hydroz[it+1][ieta][iy][ix]*(1.-dx)*(1.-dy)*dt*(1.-deta);
        gete+=hydroz[it][ieta][iy][ix+1]*(1.-dt)*(1.-dy)*dx*(1.-deta);
        gete+=hydroz[it][ieta][iy+1][ix]*(1.-dx)*(1.-dt)*dy*(1.-deta);
        gete+=hydroz[it+1][ieta][iy][ix+1]*dx*(1.-dy)*dt*(1.-deta);
        gete+=hydroz[it][ieta][iy+1][ix+1]*(1.-dt)*dy*dx*(1.-deta);
        gete+=hydroz[it+1][ieta][iy+1][ix]*(1.-dx)*dt*dy*(1.-deta);
        gete+=hydroz[it+1][ieta][iy+1][ix+1]*dx*dt*dy*(1.-deta);
        gete+=hydroz[it][ieta+1][iy][ix]*(1.-dt)*(1.-dx)*(1.-dy)*deta;
        gete+=hydroz[it+1][ieta+1][iy][ix]*(1.-dx)*(1.-dy)*dt*deta;
        gete+=hydroz[it][ieta+1][iy][ix+1]*(1.-dt)*(1.-dy)*dx*deta;
        gete+=hydroz[it][ieta+1][iy+1][ix]*(1.-dx)*(1.-dt)*dy*deta;
        gete+=hydroz[it+1][ieta+1][iy][ix+1]*dx*(1.-dy)*dt*deta;
        gete+=hydroz[it][ieta+1][iy+1][ix+1]*(1.-dt)*dy*dx*deta;
        gete+=hydroz[it+1][ieta+1][iy+1][ix]*(1.-dx)*dt*dy*deta;
        gete+=hydroz[it+1][ieta+1][iy+1][ix+1]*dx*dt*dy*deta;

        return gete;
}
