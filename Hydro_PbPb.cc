//g++ -Wall -mcmodel=medium hydro_read.cc -o hydro_read
#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <assert.h>
#include <stdlib.h>

#include "global.h"

//double hydrot[190][64][125][125];
//double hydrox[190][64][125][125];
//double hydroy[190][64][125][125];
//double hydroz[190][64][125][125];

double ****hydrot;
double ****hydrox;
double ****hydroy;
double ****hydroz;

int maxx=100;
int maxy=100;
int maxeta=64;

double deltax=0.3;
double deltay=0.3;
double deltaeta=0.203125;
double deltat=0.1;

double tau0=0.6;	//Initial time
double tau1=18.5;	//Final time
double eta1=6.5;	//Extreme abs(eta)

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

	//Allocate hydro array
	int nt=200, neta=64, nx=101, ny=101;
	hydrot = (double ****)malloc(nt * sizeof(double ***));
	hydrox = (double ****)malloc(nt * sizeof(double ***));
	hydroy = (double ****)malloc(nt * sizeof(double ***));
	hydroz = (double ****)malloc(nt * sizeof(double ***));
	if (hydrot==NULL)
        {
          fprintf(stderr, "out of memory\n");
	  exit(0);
        }
	for (int i=0; i<nt; i++)
	{         
          hydrot[i]=(double ***)malloc(neta * sizeof(double **));
          hydrox[i]=(double ***)malloc(neta * sizeof(double **));
	  hydroy[i]=(double ***)malloc(neta * sizeof(double **));
	  hydroz[i]=(double ***)malloc(neta * sizeof(double **));
	  if (hydrot[i]==NULL)
          {
            fprintf(stderr, "out of memory\n");
            exit(0);
          }
	}
	for (int i=0; i<nt; i++)
        {
	  for (int j=0; j<neta; j++) 
	  {            
            hydrot[i][j]=(double **)malloc(nx * sizeof(double *));
            hydrox[i][j]=(double **)malloc(nx * sizeof(double *));
	    hydroy[i][j]=(double **)malloc(nx * sizeof(double *));
	    hydroz[i][j]=(double **)malloc(nx * sizeof(double *));
	    if (hydrot[i][j]==NULL)    
            {   
              fprintf(stderr, "out of memory\n");         
              exit(0);
            }
          }
        }	  
	for (int i=0; i<nt; i++)
        {
	  for (int j=0; j<neta; j++) 
	  {            
	    for (int k=0; k<nx; k++)        
            {
              hydrot[i][j][k]=(double *)malloc(ny * sizeof(double));
              hydrox[i][j][k]=(double *)malloc(ny * sizeof(double));
	      hydroy[i][j][k]=(double *)malloc(ny * sizeof(double));
	      hydroz[i][j][k]=(double *)malloc(ny * sizeof(double));
              if (hydrot[i][j][k]==NULL)    
              {   
                fprintf(stderr, "out of memory\n");         
                exit(0);
              }
            }
          }
        }	  

	//hydrot = new double[190][64][125][125];
	//hydrox = new double[190][64][125][125];
	//hydroy = new double[190][64][125][125];
	//hydroz = new double[190][64][125][125];

	char hydFile[200];
	//sprintf(hydFile,"/gs/project/cqn-654-ad/mayank/events_for_jets/PbPb_%s_2p76/job-%i/results/evolution_xyeta.dat",cent.c_str(),nhyd);
	//sprintf(hydFile,"/home/peibols/projects/rrg-jeon-ac/group_writable/hydro_from_Mayank/PbPb_%s_2p76/job-%i/results/evolution_xyeta.dat",cent.c_str(),nhyd);
	sprintf(hydFile,"./hydroinfoPlaintxtHuichaoFormat.dat");
	std::ifstream hydro(hydFile);
        assert(!hydro.fail());
	//std::FILE *hydro;
        //hydro = std::fopen(hydFile, "rb");
	
	cout << " Reading Hydro..." << endl;
	clock_t startClock = clock();
        double enedat, tdat, vxdat, vydat;
        double tou, hor, ver;
        maxx=100;
        maxy=100;
        tau0=0.6;
        deltat=0.1;
        deltax=0.3;
        deltay=0.3;
        while (hydro >> hor >> ver >> tou >> enedat >> tdat >> vxdat >> vydat)//horizontal,vetical, time, energy density, temperature in fluid rest, vel in fluid restframe of plasma
        {        
                it = int((tou+deltat/2.-tau0)/deltat);//discretizing the plasma
                ix = int((hor+deltax*maxx/2.+deltax/2.)/deltax);//discretizing the plasma
                iy = int((ver+deltay*maxy/2.+deltay/2.)/deltay);//discretizing the plasma
       
                for (unsigned il=0; il<64; il++)
		{
                  hydrot[it][il][iy][ix]=tdat;
                  hydrox[it][il][iy][ix]=vxdat;
                  hydroy[it][il][iy][ix]=vydat;
        	}
	}

/*
	int t=0;
        double T, QGPfrac, vx, vy, vz;
        int size = sizeof(double);
        int do_exit=0;
	do {
		//for (int l=0; l<maxeta; l++)
		{
			int l=0; //Boost invariant smooth hydro
			for (int j=0; j<maxx; j++)
			{
				for (int k=0; k<maxy; k++)
				{
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
				}
				if (do_exit==1) break;
			}
			if (do_exit==1) break;
		}
		//cout << " t= " << t << endl;
		t+=1;
	} while (do_exit==0);
	tau1=tau0+double(t-1)*deltat;
	cout << " Tau1= " << tau1 << endl;
*/
	clock_t endClock = clock();
	cout << " Finish Hydro Read in " << double((endClock - startClock)) / CLOCKS_PER_SEC << " secs. \n";


	/*	
	double step=0.01;
	cout << " T at center= " << hydrot[0][32][50][50] << endl;
	for (unsigned int i=0; i<1000; i++)
	{
		//double x=0.8*double(i)*step;
		double temp=gT(0.6+double(i)*step,0.,0.,0.);
		cout << 0.6+double(i)*step << " " << temp << endl;
		if (temp==0.) break;
	}
	*/
}

void getGrid(double tau, double x, double y, double eta, double* udt, double* udx, double* udy, double* udeta, int* uit, int* uix, int* uiy, int* uieta)
{
        *uit=int((tau-tau0)/deltat);
        *udt=(tau-tau0-double(*uit)*deltat)/deltat;

        if (y>=0.) {
                *uiy = int(y/deltay)+(maxy)/2;
                *udy = (y - double(*uiy-(maxy)/2)*deltay)/deltay;
        }
        else {
                *uiy = int(y/deltay)+(maxy)/2-1;
                *udy = (y - double(*uiy-(maxy)/2)*deltay)/deltay;
        }

        if (x>=0.) {
                *uix = int(x/deltax)+(maxx)/2;
                *udx = (x - double(*uix-(maxx)/2)*deltax)/deltax;
        }
        else {
                *uix = int(x/deltax)+(maxx)/2-1;
                *udx = (x - double(*uix-(maxx)/2)*deltax)/deltax;
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

	if (ix<0 || ix>=maxx || iy<0 || iy>=maxy || ieta<0 || ieta>=maxeta-1) return gete;
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

	return gete*0.2;
}

double gVx(double tau, double x, double y, double eta)
{
        double gete=0.;

        if (tau>=tau1 || abs(eta)>=eta1 || tau<tau0) {
                //cout << " No Hydro here: tau= " << tau << ", eta= " << eta << endl;
                return gete;
        }

        getGrid(tau,x,y,eta,&dt,&dx,&dy,&deta,&it,&ix,&iy,&ieta);

        if (ix<0 || ix>=maxx || iy<0 || iy>=maxy || ieta<0 || ieta>=maxeta-1) return gete;
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

        getGrid(tau,x,y,eta,&dt,&dx,&dy,&deta,&it,&ix,&iy,&ieta);

        if (ix<0 || ix>=maxx || iy<0 || iy>=maxy || ieta<0 || ieta>=maxeta-1) return gete;
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

        if (ix<0 || ix>=maxx || iy<0 || iy>=maxy || ieta<0 || ieta>=maxeta-1) return gete;
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
