#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "Random.h"

#include "vector_operators.h"

//#define DO_SOURCE

using std::vector;
using namespace std;

void get_source_evol(double &tau_ev, double& x_f, double& y_f, double& vx_f, double& vy_f, double tau_ini, double x_ini, double y_ini, double Tc, int ebe_hydro);
void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod, int mode, double &length, double &tlength, int ebe_hydro);
double normalise(vector<double> &p);
vector<double> vec_prod(vector<double> a, vector<double> b);
double gQ(double Del, numrand &nr, double cutoff);
void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp, double vscalw, double lore, double step, double kappa, numrand &nr);
double gVx(double tau, double x, double y, double eta);
double gVy(double tau, double x, double y, double eta);
double gVz(double tau, double x, double y, double eta);
double gT(double tau, double x, double y, double eta);
double call_gT(double tau, double x, double y, int comp);

void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod, int mode, double &length, double &tlength, int ebe_hydro)
{
  double Tc;
  if (tmethod==0) Tc=0.170;
  else Tc=0.145;

  double tot=pos[3]+tof;    //Final time

  double tau0h=0.6;
  if (ebe_hydro==1) tau0h=0.4;

  double ei=p[3];      //Initial energy

  double f_dist=0.;    //Traversed distance in Fluid Frame
  double l_dist=0.;    //Traversed distance in Lab Frame
  
  double inmed_f_dist=0.;
  double virt_f_dist=0.;

  //std::ofstream lengthfile;
  //lengthfile.open("lengthdist.txt", std::ios::app);

  double CF;
  if (id==21) {
    if (mode==0) CF=pow(9./4.,1./3.);  //If gluon, color charge dependence is ratio of casimirs to power 1/3
    else CF=9./4.;
  }
  else CF=1.;

  int marker=0;    //If one, exit loop
  double step=0.1;  //Time step in LAB frame

  vector<double> w = p/p[3];  //4-velocity
  //cout << " x= " << pos[0] << " y= " << pos[1] << " z= " << pos[2] << " t= " << pos[3] << endl;
  //cout << " px= " << p[0] << " py= " << p[1] << " pz= " << p[2] << " en= " << p[3] << endl;
  //cout << " wx= " << w[0] << " wy= " << w[1] << " wz= " << w[2] << " we= " << w[3] << endl;
  do {

    //Keep 4momentum before applying quenching this step
    vector<double> p_prev = p;

    if (pos[3]==tot) marker=1;
    if (pos[3]>tot) cout << " Warning: Went beyond tot= " << tot << " t= " << pos[3] << endl;  
    
    //Proper time
    double tau=sqrt(pos[3]*pos[3]-pos[2]*pos[2]);
    if (tau!=tau) {
      tau=0.;
      cout << " TAU Not a number z= " << pos[2] << " t= " << pos[3] << " wz= " << w[2] << " en = " << p[3] << " pz= " << p[2] << "\n";
      cout << " Id= " << id << endl;
      exit(1);
    }

    //Rapidity
    double eta = 1./2.*log((pos[3]+pos[2])/(pos[3]-pos[2]));
    //if (abs(eta)>=10.) cout << " eta= " << eta << endl;
    if (eta!=eta && tau>0.) {
      cout << " Eta is NaN= " << eta << " t= " << pos[3] << " z= " << pos[2] << endl;
      eta=0.;
    }

    int will_hot=0;  //Advance variable (to reach hot zones)
    double vx=0.;
    double vy=0.;
    if (tau>=tau0h)  //Smooth profile starting time
    {
      vector<double> v;
      if (ebe_hydro==0) {
        vx=gVx(tau,pos[0],pos[1],eta);
        vy=gVy(tau,pos[0],pos[1],eta);
      }
      else {
        vx=call_gT(tau,pos[0],pos[1],1);
        vy=call_gT(tau,pos[0],pos[1],2);
      }
      //double vz=gVz(tau,pos[0],pos[1],eta);
      //**BOOST INVARIANT**
      double vz=pos[2]/pos[3];
      double frap=atanh(vz);
      vx/=cosh(frap);
      vy/=cosh(frap);
      //*******************
      v.push_back(vx), v.push_back(vy), v.push_back(vz), v.push_back(1.);
      
      double v2=pow(v[0],2.)+pow(v[1],2.)+pow(v[2],2.);
      double w2=pow(w[0],2.)+pow(w[1],2.)+pow(w[2],2.);
      double vscalw=v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
      if (v2>=1.) v2=0.999999999, cout << " V2 >= 1 \n";
      double lore=1./sqrt(1.-v2);

      l_dist+=step;
      double f_lore=w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw);
      if (f_lore < 0.) { 
        //cout << " WTFFFFFFFFFF f_lore " << "\n"; 
        //cout << w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw) << " lore= " << lore << " vscalw= " << vscalw << endl;
	//cout << " v= " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << endl;
	//cout << " w= " << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << endl;
	//cout << " p= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
	//cout << " w2= " << w2 << endl;	
	f_lore=0.;
      }
      double f_step = step*sqrt(f_lore);
      f_dist+=step*sqrt(f_lore);

      if (f_dist > l_dist ) {
        //cout << " f_dist= " << f_dist << " l_dist=" << l_dist << endl;
      }

      double temp;
      if (ebe_hydro==0) temp=gT(tau,pos[0],pos[1],eta);
      else temp=call_gT(tau,pos[0],pos[1],0);

      if (temp>Tc) {
        inmed_f_dist+=f_step;
      }

      //Safe way to exit the plasma: check whether temperature will be above Tc in the next 1000 steps
      if (temp<Tc)
      {
        //Check whether temperature increases in its way
        for (unsigned int j=1; j<1000; j++)
        {
          vector<double> tpos=pos+w*step*double(j);
          if (tpos[3]>tot) break;
          tau=sqrt(tpos[3]*tpos[3]-tpos[2]*tpos[2]);
          eta = 1./2.*log((tpos[3]+tpos[2])/(tpos[3]-tpos[2]));
          double ctemp;
	  if (ebe_hydro==0) ctemp=gT(tau,pos[0],pos[1],eta);
	  else ctemp=call_gT(tau,tpos[0],tpos[1],0);
	  if (ctemp>Tc)
          {
            //cout << " Will Hot! from time= " << pos[3] << " at j= " << j << endl;
            will_hot=int(j);
            break;
          }
        }
        if (will_hot==0)
        {
          pos+=w*(tot-pos[3]);
          marker=1;
        }
      }

      //Now broad&quench
      if (p[3]>0. && temp>=Tc && f_step!=0.)
      {

        if (mode==0) {
          length += f_step;
          //cout << " stepping length= " << length << endl;
          tlength += temp / 0.2 * f_step;
        }
        else {
          length += f_step;
          //cout << " stepping length= " << length << endl;
          tlength += pow(temp / 0.2,2.) * f_step;
        }
        //Broadening
        if (kappa!=0.) {
          trans_kick(w,w2,v,p,temp,vscalw,lore,step,kappa,nr);
	}
	//Strong coupling
        if (alpha!=0. && mode==0)
        {
          double Efs=ei*lore*(1.-vscalw);
          double tstop=0.2*pow(Efs,1./3.)/(2.*pow(temp,4./3.)*alpha)/CF;
          double beta=tstop/f_dist;
          if (beta>1.)
          {
            double intpiece=Efs*step*4./(3.141592)*(1./(beta*tstop*sqrt(beta*beta-1.)));
            //Update 4momentum
            double quench=(p[3]-intpiece)/p[3];
            p*=quench;
          }
          else
          {
            p[3]=0.;
          }
        }
        //Radiative
        if (alpha!=0. && mode==1)
        {
          double intpiece=CF*(step/0.2)*alpha*temp*temp*temp*(f_dist/0.2);
          double quench=(p[3]-intpiece)/p[3];
          p*=quench;
        }
        //Collisional
        if (alpha!=0. && mode==2)
        {
          double intpiece=CF*(step/0.2)*alpha*temp*temp;
          double quench=(p[3]-intpiece)/p[3];
          p*=quench;
        }
      } 
    }

    if (tof<1000000) {
      double vz=pos[2]/max(pos[3],0.000001);
      double frap=atanh(vz);
      double v2=vz*vz;
      double w2=pow(w[0],2.)+pow(w[1],2.)+pow(w[2],2.);
      double vscalw=vz*w[2];
      if (v2>=1.) v2=0.999999999, cout << " V2 >= 1 \n";
      double lore=1./sqrt(1.-v2);
      double f_lore=w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw);
      if (f_lore<0.) {
        cout << " craazy f_lore= " << f_lore << endl;
	f_lore=1.;
      }
      virt_f_dist+=step*sqrt(f_lore);
    }

    //If parton gets totally quenched, exit 
    if (p[3]<=0.)
    {
      marker=1;
      for (unsigned int i=0; i<4; i++) p[i]=0.;
    }
    else {
      //Manually protect very soft particles from getting kicks that yield velocities greater than 1
      for (unsigned int i=0; i<3; i++) {
        if (p[i]>p[3]) {
	  cout << " Got crazy kick in i= " << i << "p[i]= " << p[i] << " and p[3]= " << p[3] << endl;
	  p[i]=0.99999*p[3];
        }
      }
      //Update kinematical quantities, with the possibility of advancing to hot regions
      w=p/p[3];
      double tstep=max(double(will_hot),1.)*step;
      if (pos[3]+tstep>tot)
      {
        tstep=tot-pos[3];
      }
      if (marker!=1) pos+=w*tstep;
    }


  #ifdef DO_SOURCE
    //Fill source file
    if (p[3]!=p_prev[3]) {
      
      //Get tau ev, x_f, y_f and vx_f and vy_f for source file
      double tau_ev, x_f, y_f, vx_f, vy_f;
      get_source_evol(tau_ev,x_f,y_f,vx_f,vy_f,tau,pos[0],pos[1],Tc,ebe_hydro);
      
      ofstream source_file ("SOURCE.dat",std::ios_base::app);
      source_file << tau << " " << pos[0] << " " << pos[1] << " " << eta << " " << tau_ev << " "
	          << -p[3]+p_prev[3] << " " << -p[0]+p_prev[0] << " " << -p[1]+p_prev[1] << " " << -p[2]+p_prev[2] << " "
		  << vx << " " << vy << " " << vx_f << " " << vy_f << " " << x_f << " " << y_f << endl;
    }
  #endif

  } while (marker==0);

  if (virt_f_dist==0.) virt_f_dist=-1;
  //lengthfile << inmed_f_dist << " " << virt_f_dist << endl;
  //lengthfile.close();

  return;
}

void get_source_evol(double &tau_ev, double& x_f, double& y_f, double& vx_f, double& vy_f, double tau_ini, double x_ini, double y_ini, double Tc, int ebe_hydro)
{

  tau_ev = 0.;
  double dtau = 0.1;
  x_f = x_ini;
  y_f = y_ini;
  double tau_now = tau_ini;
  
  while (true) {

    double vx_f, vy_f;
    if (ebe_hydro==0) {
      vx_f = gVx(tau_now,x_f,y_f,0.);
      vy_f = gVy(tau_now,x_f,y_f,0.);
    }
    else {
      vx_f = call_gT(tau_now,x_f,y_f,1);
      vy_f = call_gT(tau_now,x_f,y_f,2);
    }
    
    double localT;
    if (ebe_hydro==0) localT=gT(tau_now,x_f,y_f,0.);
    else localT=call_gT(tau_now,x_f,y_f,0);
    if (localT < Tc) break;

    x_f += vx_f * dtau;
    y_f += vy_f * dtau;
    tau_ev += dtau;
    tau_now += dtau;
  
  }

  return;
}

void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp, double vscalw, double lore, double step, double kappa, numrand &nr)
{
  if (vscalw==1.) return;

  vector<double> e1 = vec_prod(w,v);
  double Ne1=normalise(e1);
  if (Ne1==0.) {
    double b=0.5;
    double c=0.2;
    double a=(-b*w[1]-c*w[2])/w[0];
    e1={a,b,c,0.};
    Ne1=normalise(e1);
  }

  double Nw=sqrt(w2);
  vector<double> l = vec_prod(w,e1)/Nw;

  double uscalW=lore*(1.-vscalw);
  double uscall=lore*(-v[0]*l[0]-v[1]*l[1]-v[2]*l[2]);
  double W2=1.-w2;

  vector<double> Wp = w-v*lore*W2/uscalW;
  //cout << " Wp scal w= " << w[3]*Wp[3]-w[0]*Wp[0]-w[1]*Wp[1]-w[2]*Wp[2] << endl;

  double Nalpha=-uscall*uscalW/(pow(uscalW,2.)-W2);
  if (Nalpha!=Nalpha || isinf(Nalpha)) return;
  double NN=1.+W2*pow(uscall,2.)/(-pow(uscalW,2.)+W2);
  //In some rare situations, this norm squared can be negative. Only do kick otherwise
  if (sqrt(NN)!=sqrt(NN)) cout << " negative NN " << endl;
  else {
    vector<double> e2 = (l+Wp*Nalpha)/sqrt(NN);
    //cout << " norm e2= " << e2[3]*e2[3]-e2[0]*e2[0]-e2[1]*e2[1]-e2[2]*e2[2] << endl;

    double Ef=p[3]*lore*(1.-vscalw);
    double wf2=1.-W2/pow(lore*(1.-vscalw),2.);
    double DelQ2=kappa*pow(temp,3.)*lore*(1.-vscalw)*step*5.;

    double qfac=0.;
    //Only do kick if energy is greater than temperature
    //if (Ef>temp && Ef*sqrt(wf2)>0.)
    if (Ef*sqrt(wf2)>0.)
    {
      //cout << " DelQ2= " << DelQ2 << " cutoff= " << Ef*sqrt(wf2) << endl;
      //qfac=gQ(DelQ2,nr,Ef*sqrt(wf2));
      qfac = sqrt(-1.*log(nr.rando()))*sqrt(DelQ2); //Box-Muller method
      if (qfac>Ef*sqrt(wf2)) {
        qfac=Ef*sqrt(wf2)-0.00000001;
      }
    }
  
    double qbeta;
    if (wf2>0.) qbeta=sqrt(1.-qfac*qfac/Ef/Ef/wf2)-1.;
    else qbeta=0.;
    if (qbeta!=qbeta) {
      cout << " qbeta= " << setprecision(6) << qbeta << " qfac= " << qfac << " Ef= " << Ef << " wf2= " << wf2 << " cutoff= " << Ef*sqrt(wf2) << endl;
      cout << " W2= " << setprecision(6) << W2 << " lore= " << lore << " vscalw= " << vscalw << endl;
    }

    double qphi=2.*3.141592654*nr.rando(); 

    vector<double> e = e1*cos(qphi)+e2*sin(qphi); 

    vector<double> Wt = (w-v*uscalW*lore)/lore/(1.-vscalw);

    //Update 4momentum
    if (lore==1.) {
      e2=vec_prod(e1,w);
      double Ne2=normalise(e2);
      e=e1*cos(qphi)+e2*sin(qphi);
    }
    p+=Wt*qbeta*Ef+e*qfac;
    if (p[3]!=p[3]) {
      cout << " p in bro= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
      cout << " qbeta= " << qbeta << " qfac= " << qfac << endl;
      cout << " Wt= " << Wt[0] << " " << Wt[1] << " " << Wt[2] << " " << Wt[3] << endl;
      cout << " e1= " << e1[0] << " " << e1[1] << " " << e1[2] << " " << e1[3] << endl;
      cout << " e2= " << e2[0] << " " << e2[1] << " " << e2[2] << " " << e2[3] << endl;
      cout << " Nalpha= " << Nalpha << endl;
      cout << " l= " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << endl;
      cout << " Wp= " << Wp[0] << " " << Wp[1] << " " << Wp[2] << " " << Wp[3] << endl;
      cout << " uscalW= " << uscalW << endl;
      cout << " Nw= " << Nw << " uscall= " << uscall << endl;
      cout << " W2= " << W2 << endl;
      exit(1);
    }

  }
}

//obselete, using BoxMuller now
double gQ(double Del, numrand &nr, double cutoff)
{
  double gaussq;
  double qfac;
  double gaussmax=sqrt(2./Del/exp(1.));
  if (cutoff<sqrt(Del/2.)) gaussmax=2.*cutoff/Del*exp(-cutoff*cutoff/Del);
  qhatelsen:
  qfac=min(4.*sqrt(Del),cutoff)*nr.rando();
  gaussq=2.*qfac/Del*exp(-qfac*qfac/Del)/gaussmax;
  double nrand=nr.rando();
  if (gaussq>1.) std::cout << "Wrong normalisation in gQ! gaussq= " << gaussq << std::endl << std::endl;
  if (nrand>gaussq) goto qhatelsen;

  return qfac;
}
