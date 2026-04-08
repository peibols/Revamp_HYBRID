//g++ -Wall -mcmodel=medium hydro_read.cc -o hydro_read
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
#include <cmath>
#include <assert.h>
#include <stdlib.h>

double hydrot_ebe[100][128][128];
double hydrox_ebe[100][128][128];
double hydroy_ebe[100][128][128];

int ixmax;
int ietamax;
int itaumax;

double hydroDx;
double hydroXmax;
double hydroDeta;
double hydro_eta_max;
double hydroTau0;
double hydroDtau;
double hydroTauMax;

double etaovers=0.12;
double hbarc=0.197327;

double tauofwtilde(double Tau13Teq, double tau, double w);
double Eofw(double w);
double getTeff(double Thyd, double tau);
void getGrid(double tau, double x, double y, double* dt, double* dx, double* dy, int* it, int* ix, int* iy);
double gT_ebe(double tau, double x, double y, double T[100][128][128]);
double call_gT(double tau, double x, double y, int i);
void read_hydro_ipsat(int, std::string);

using std::vector;
using namespace std;

/*
template<typename T>
std::istream & binary_read(std::istream& stream, T& value){
    return stream.read(reinterpret_cast<char*>(&value), sizeof(T));
}
*/

void read_hydro_ipsat(int nhyd, std::string cent)
{

        std::ostringstream filename;
	//filename << "hydros_" << cent.c_str() << "/job-" << nhyd << "/evolution_all_xyeta.dat";
	filename << "./evolution_all_xyeta.dat";
        std::FILE *hydro;
	hydro = std::fopen(filename.str().c_str(), "rb");
        if (hydro == NULL) {
          cout << "Hydro open fail = " << filename.str() << endl;
          exit(1);
        }

	cout << " Reading Hydro..." << endl;
	clock_t startClock = clock();

	// Get header
        float header[16];
	int status = std::fread(&header, sizeof(float), 16, hydro);
        if (status == 0) {
            cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                 << "Can not read the evolution file header" << endl;
            exit(1);
        }

	hydroTau0 = double(header[0]);
        hydroDtau = double(header[1]);
        ixmax = static_cast<int>(header[2]);
        hydroDx = double(header[3]);
        hydroXmax = double(std::abs(header[4]));
        ietamax = static_cast<int>(header[8]);
        hydroDeta = double(header[9]);
        hydro_eta_max = double(std::abs(header[10]));
        int turn_on_rhob = static_cast<int>(header[11]);
        int turn_on_shear = static_cast<int>(header[12]);
        int turn_on_bulk = static_cast<int>(header[13]);
        int turn_on_diff = static_cast<int>(header[14]);
        const int nVar_per_cell = static_cast<int>(header[15]);

        cout << "tau0= " << hydroTau0
	     << " dtau= " << hydroDtau
	     << " dx= " << hydroDx
	     << " xmax= " << hydroXmax
	     << " ixmax= " << ixmax
	     << " deta= " << hydroDeta
	     << " etamax= " << hydro_eta_max
	     << " ietamax= " << ietamax
	     << endl;

	float cell_info[nVar_per_cell];

	int itau_max = 0;
        double maxtemp=0.;

        int ik = 0;
        while (true) {
            status = 0;
            status = std::fread(&cell_info, sizeof(float), nVar_per_cell, hydro);
            if (status == 0) break;
            if (status != nVar_per_cell) {
                cerr << "[Hydroinfo_MUSIC::readHydroData]: ERROR: "
                     << "the evolution file format is not correct" << endl;
                exit(1);
            }

            if (itau_max < static_cast<int>(cell_info[0]))
                itau_max = static_cast<int>(cell_info[0]);
            int itau = static_cast<int>(cell_info[0]);
            int ix   = static_cast<int>(cell_info[1]);
            int iy   = static_cast<int>(cell_info[2]);
            int ieta = static_cast<int>(cell_info[3]);
            float temperature = cell_info[6];
            float ed = cell_info[4];
            float pressure = cell_info[5];
            float ux = cell_info[8];
            float uy = cell_info[9];
            float uz = cell_info[10];
	    //gamma = sqrt(1+u^2)
	    float u2 = ux*ux+uy*uy+uz*uz;
	    float gamma = sqrt(1+u2);
	    float vx = ux/gamma;
	    float vy = uy/gamma;

	    hydrot_ebe[itau][iy][ix]=double(temperature);
	    hydrox_ebe[itau][iy][ix]=double(vx);
	    hydroy_ebe[itau][iy][ix]=double(vy);

	    if (temperature>maxtemp) maxtemp=temperature;

	    ik++;
            if (ik%50000 == 0)
                cout << "o" << flush;
        }	

        cout << "itau_max= " << itau_max << endl;
        cout << " Max temp= " << maxtemp << endl;

	itaumax = itau_max;
	hydroTauMax = hydroTau0 + hydroDtau*itaumax;

	clock_t endClock = clock();
	cout << " Finish Hydro Read in " << double((endClock - startClock)) / CLOCKS_PER_SEC << " secs. \n";
	
	/*	
	double step=0.1;
	cout << " T at center= " << hydrot[0][64][64] << endl;
	for (unsigned int i=0; i<1000; i++)
	{
		//double x=0.8*double(i)*step;
		double temp=call_gT(hydroTau0+double(i)*step,0.,0.,0);
		cout << "Temp test= " << hydroTau0+double(i)*step << " " << temp << endl;
		if (temp==0.) break;
	}
	*/
	
	
}

void getGrid_ebe(double tau, double x, double y, double* udt, double* udx, double* udy, int* uit, int* uix, int* uiy)
{
        *uit=int((tau-hydroTau0)/hydroDtau);
        *udt=(tau-hydroTau0-double(*uit)*hydroDtau)/hydroDtau;

	*uix = int((hydroXmax+x)/hydroDx);
	*udx = (x-(double(*uix)*hydroDx-hydroXmax))/hydroDx;

	*uiy = int((hydroXmax+y)/hydroDx);
	*udy = (y-(double(*uiy)*hydroDx-hydroXmax))/hydroDx;
}

double call_gT(double tau, double x, double y, int i)
{
  double var;

  double etau=tau;
  if (tau<hydroTau0) {
    etau=hydroTau0;
  }

  if (i==0) {
    var = gT_ebe(etau,x,y,hydrot_ebe);
  }
  else if (i==1) {
    var = gT_ebe(etau,x,y,hydrox_ebe);
  }
  else if (i==2) {
    var = gT_ebe(etau,x,y,hydroy_ebe);
  }
  else {
    cout << " Wrong prof= " << i << endl;
    exit(1);
  }

  if (tau<hydroTau0) {
    if (i==1 || i==2) {
      var = var * (tau/hydroTau0);
    }
    if (i==0) {
      //var = 0.6;
      var = getTeff(var,tau);
      //cout << tau << " " << var << endl;
    }
  }

  return var;
}

double gT_ebe(double tau, double x, double y, double T[100][128][128])
{
        double gete=0.;

        if (tau>=hydroTauMax || tau<hydroTau0) {
		//cout << " No Hydro here: tau= " << tau << ", eta= " << eta << endl;
                return gete;
        }

        double dt, dx, dy;
        int it, ix, iy;
	getGrid_ebe(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);

        if (ix<0 || ix>=ixmax-1 || iy<0 || iy>=ixmax-1 || it<0 || it>=itaumax-1) return gete;
        gete=T[it][iy][ix]*(1.-dt)*(1.-dx)*(1.-dy);
        gete+=T[it+1][iy][ix]*(1.-dx)*(1.-dy)*dt;
        gete+=T[it][iy][ix+1]*(1.-dt)*(1.-dy)*dx;
        gete+=T[it][iy+1][ix]*(1.-dx)*(1.-dt)*dy;
        gete+=T[it+1][iy][ix+1]*dx*(1.-dy)*dt;
        gete+=T[it][iy+1][ix+1]*(1.-dt)*dy*dx;
        gete+=T[it+1][iy+1][ix]*(1.-dx)*dt*dy;
        gete+=T[it+1][iy+1][ix+1]*dx*dt*dy;

        return gete;
}

double getTeff(double Thyd, double tau)
{

  double tau13Teq=pow(hydroTau0/hbarc,1./3.)*(Thyd+2./3.*etaovers/(hydroTau0/hbarc));

  double tole=0.0000001;
  double hi=6.29;
  double lo=0.0010001;

  double taumin=pow(4.*M_PI*etaovers*lo/tau13Teq,3./2.)/pow(Eofw(lo),3./8.);
  double taumax=pow(4.*M_PI*etaovers*hi/tau13Teq,3./2.)/pow(Eofw(hi),3./8.);

  if (tau<taumin) {
    //cout << "tau= " << tau << " taumin= " << taumin << endl;
    return 0;
  }
  if (tau>taumax) {
    cout << "tau= " << tau << " taumax= " << taumax << endl;
    return 0;
  }
  
  double fhi=tauofwtilde(tau13Teq,tau,hi);
  double flo=tauofwtilde(tau13Teq,tau,lo);
  double wtilde=-1000;
  double mid_prev=0.;
  while (true)
  {
    double mid=(hi+lo)/2.;
    double fmid=tauofwtilde(tau13Teq,tau,mid);
    //if (fabs(fmid)<tole) {
    if (fabs((mid-mid_prev)/mid)<tole) {
      wtilde=mid;
      break;
    }
    mid_prev = mid;
    if (fhi*fmid>0.) {
      fhi=fmid;
      hi=mid;
    }
    else {
      flo=fmid;
      lo=mid;
    }
    if (flo*fhi>0.) {
      cout << "Same sign!" << endl;
      exit(1);
    }
  }

  double theTemp = pow(tau13Teq,3./2.)*pow(Eofw(wtilde),3./8.)/pow(4.*M_PI*etaovers*wtilde,1./2.);
  return theTemp;

}

double tauofwtilde(double tau13Teq, double tau, double w)
{
  return tau/hbarc-pow(4.*M_PI*etaovers*w/tau13Teq,3./2.)/pow(Eofw(w),3./8.);
}

double Eofw(double w)
{
  double ah=1.25160583540026846;
  double bh=0.0167194722583442507;
  double ch=-0.375124344952100277;

  double al=1.04038563732240692;
  double bl=1.87248466419276916;
  double cl=-0.00567650680970884004;

  if (w>0.037) return ah*pow(tanh(w),1./5.)+bh*w+ch;
  else return al*pow(tanh(w),1./2.5)+bl*w*w+cl; 
}
