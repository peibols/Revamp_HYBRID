#include <fstream>
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <assert.h>
#include <stdlib.h>

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

void getGrid(double tau, double x, double y, double eta, double* dt, double* dx, double* dy, double* deta, int* it, int* ix, int* iy, int* ieta);
double gT(double tau, double x, double y, double eta);
double gVx(double tau, double x, double y, double eta);
double gVy(double tau, double x, double y, double eta);
double gVz(double tau, double x, double y, double eta);
void read_hydro(int, std::string);

using std::vector;
using namespace std;

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

  char hydFile[200];
  sprintf(hydFile,"./hydroinfoPlaintxtHuichaoFormat.dat");
  std::ifstream hydro(hydFile);
  assert(!hydro.fail());

  //For realistic wake
  double min_freeze_tau = 2.;
  double freeze_tau_cut = 9.;
  double tc_low = 0.142;
  double tc_high = 0.148;
  double f_vx[37][37][2][2]={{{{0.}}}};
  double f_vy[37][37][2][2]={{{{0.}}}};
  double f_tau[37][37][2][2]={{{{0.}}}};
  double fdeltax=0.6;
  double fdeltay=0.6;
  int fnx=37;
  int fny=37;
  double fmax_valx=11.;
  double fmax_valy=11.;

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
  double max_valx=15.;
  double max_valy=15.;
  int it, ix, iy;
  while (hydro >> hor >> ver >> tou >> enedat >> tdat >> vxdat >> vydat)//horizontal,vetical, time, energy density, temperature in fluid rest, vel
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

    //Store info for realistic wake
    if (tou < min_freeze_tau) continue;
    if (tdat*0.2 < tc_low || tdat*0.2 > tc_high) continue;
    int fix = int((hor+fdeltax*(fnx-1)/2.+fdeltax/2.)/fdeltax);
    int fiy = int((ver+fdeltay*(fny-1)/2.+fdeltay/2.)/fdeltay);
    if (fix<0 || fiy<0 || fix>=fnx-1 || fiy>=fny-1) continue;
    int ib = 0;
    if (tou < freeze_tau_cut) ib = 1;
    f_vx[fix][fiy][0][ib] += vxdat;
    f_vx[fix][fiy][1][ib] += 1;
    f_vy[fix][fiy][0][ib] += vydat;
    f_vy[fix][fiy][1][ib] += 1;
    f_tau[fix][fiy][0][ib] += tou;
    f_tau[fix][fiy][1][ib] += 1;

  }
 
  clock_t endClock = clock();
  cout << " Finish Hydro Read in " << double((endClock - startClock)) / CLOCKS_PER_SEC << " secs. \n";

  //Output flow and dtauf at freezeout
  ofstream flowfile("freezeout_flow.dat");
  for (int i=0; i<fnx; i++) {
    double x_cen = -fmax_valx + double(i)*fdeltax;
    for (int j=0; j<fny; j++) {
      double y_cen = -fmax_valy + double(j)*fdeltay;
      flowfile << x_cen << " " << y_cen << " ";
      for (int ib=0; ib<2; ib++) {
        flowfile << f_vx[i][j][0][ib]/max(f_vx[i][j][1][ib],1.) << " "
	         << f_vy[i][j][0][ib]/max(f_vy[i][j][1][ib],1.) << " "
		 << f_tau[i][j][0][ib]/max(f_tau[i][j][1][ib],1.) << " ";
      }
      flowfile << endl;
    }
    flowfile << endl;
  }
  flowfile.close();
  
  ofstream gradfile("freezeout_dtauf.dat");
  for (int i=0; i<fnx; i++) {
    double x_cen = -fmax_valx + double(i)*fdeltax;
    for (int j=0; j<fny; j++) {
      double y_cen = -fmax_valy + double(j)*fdeltay;
      gradfile << x_cen << " " << y_cen << " ";
      for (int ib=0; ib<2; ib++) {
        if (i==fnx-1 || j==fny-1) { gradfile << 0. << " " << 0. << " "; continue; }
	double tau_now = f_tau[i][j][0][ib]/max(f_tau[i][j][1][ib],1.);
	double tau_dx = f_tau[i+1][j][0][ib]/max(f_tau[i+1][j][1][ib],1.);
	double tau_dy = f_tau[i][j+1][0][ib]/max(f_tau[i][j+1][1][ib],1.);
        if (tau_dx!=0. && tau_now!=0.) gradfile << std::setprecision(7) << (tau_dx-tau_now)/fdeltax << " ";
        else gradfile << 0. << " ";
        if (tau_dy!=0. && tau_now!=0.) gradfile << std::setprecision(7) << (tau_dy-tau_now)/fdeltay << " ";
        else gradfile << 0. << " ";
      }
      gradfile << endl;
    }
    gradfile << endl;
  }
  gradfile.close();

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

  return;

}

double gT(double tau, double x, double y, double eta)
{
  double gete=0.;

  if (tau>=tau1 || abs(eta)>=eta1 || tau<tau0) {
    return gete;
  }

  double dt, dx, dy, deta;
  int it, ix, iy, ieta;
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
    return gete;
  }

  double dt, dx, dy, deta;
  int it, ix, iy, ieta;
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
    return gete;
  }

  double dt, dx, dy, deta;
  int it, ix, iy, ieta;
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
    return gete;
  }


  double dt, dx, dy, deta;
  int it, ix, iy, ieta;
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
