#include <iostream>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include <iomanip>
#include "Quench.h"
#include "Random.h"
#include "FourVector.h"

#include "global.h"
#include "vector_operators.h"
#include "Distributions.hpp"

#include "gsl/gsl_integration.h"

#define DO_ELASTIC

using std::vector;
using namespace std;

void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod, int model, vector<Quench> &new_particles, int &had_scattering, vector<double> &orient, double wake_dep[4], double wake_rec[4]);
double normalise(vector<double> &p);
vector<double> vec_prod(vector<double> a, vector<double> b);
double gQ(double Del, numrand &nr, double cutoff);
void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp, double vscalw, double lore, double step, double kappa, numrand &nr);
double gVx(double tau, double x, double y, double eta);
double gVy(double tau, double x, double y, double eta);
double gVz(double tau, double x, double y, double eta);
double gT(double tau, double x, double y, double eta);
FourVector Boost( double b[3], FourVector p);
FourVector BoostBack( double b[3], FourVector p);
double Delta(FourVector v1, FourVector v2);

void loss_rate(vector<double> &p, vector<double> &pos, double tof, int id, numrand &nr, double kappa, double alpha, int tmethod, int model, vector<Quench> &new_particles, int &had_scattering, vector<double> &orient, double wake_dep[4], double wake_rec[4])
{

#ifdef DO_ELASTIC
  gsl_integration_workspace *wdk = gsl_integration_workspace_alloc(100000);
  gsl_integration_workspace *wkcm = gsl_integration_workspace_alloc(100000);
  gsl_integration_workspace *wx = gsl_integration_workspace_alloc(100000);
#endif

  //Source File, append
  ofstream source_file;   
  source_file.open(Sfile, ios::app);

  double Tc;
  if (tmethod==0) Tc=0.170;
  else Tc=0.145;

  double charm_mass=1.25;
  double b_mass=4.2;

  double tot=pos[3]+tof;		//Final time

  double ei=p[3];			//Initial energy

  double f_dist=0.;		//Traversed distance in Fluid Frame
  double l_dist=0.;		//Traversed distance in Lab Frame

  double CF;
  if (id==21) {
    if (model==0) CF=pow(9./4.,1./3.);     //If gluon, color charge dependence is ratio of casimirs to power 1/3
    else CF=9./4.;
  }
  else CF=1.;

  int marker=0;		//If one, exit loop
  double step=0.1;	//Time step in LAB frame

  vector<double> p_ini = p;

  vector<double> w = p/p[3];	//4-velocity
  //cout << " x= " << pos[0] << " y= " << pos[1] << " z= " << pos[2] << " t= " << pos[3] << endl;
  //cout << " px= " << p[0] << " py= " << p[1] << " pz= " << p[2] << " en= " << p[3] << endl;
  //cout << " wx= " << w[0] << " wy= " << w[1] << " wz= " << w[2] << " we= " << w[3] << endl;
  vector<double> o_in = orient;
    
  double tauS=sqrt(pos[3]*pos[3]-pos[2]*pos[2]);
  if (tauS!=tauS) {
    tauS=0.;
    cout << " Start TAU Not a number z= " << pos[2] << " t= " << pos[3] << " wz= " << w[2] << " en = " << p[3] << " pz= " << p[2] << "\n";
    exit(1);
  }

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

    int will_hot=0;	//Advance variable (to reach hot zones)
    double vx, vy, vz;
    if (tau>=0.6)	//Smooth profile starting time
    {
      vector<double> v;
      vx=gVx(tau,pos[0],pos[1],eta);
      vy=gVy(tau,pos[0],pos[1],eta);
      //double vz=gVz(tau,pos[0],pos[1],eta);
      //**BOOST INVARIANT**
      vz=pos[2]/pos[3];
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

      double temp=gT(tau,pos[0],pos[1],eta);
      
      l_dist+=step;
      double f_lore=w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw);
      if (f_lore < 0.) { 
        //cout << " WTFFFFFFFFFF temp= " << temp << "\n"; 
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
	  if (gT(tau,tpos[0],tpos[1],eta)>Tc)
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
      if (p[3]>0. && temp>=Tc)
      {
      #ifdef DO_ELASTIC
        //std::cout << " use_tables = " << use_tables << endl;
        FourVector pp;
	pp.Set(p[0],p[1],p[2],std::sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));
        double beta[3]={v[0],v[1],v[2]};
	pp = Boost(beta,pp);

        //pin check
	double pin = sqrt(pp.x()*pp.x()+pp.y()*pp.y()+pp.z()*pp.z()) / temp;

        if (pin>1500.) {
          std::cout << " Would have skipped because too large energy, " << 
	  		" p_FLUID= " << sqrt(pp.x()*pp.x()+pp.y()*pp.y()+pp.z()*pp.z()) <<
	  		" p_LAB= " << sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]) << " temp = " << temp << endl;
	}

      	//if ((id==21 || abs(id)<=4) && pin<=1000. && pin>2.) {
      	if ((id==21 || abs(id)<=4) && pin<1500) {
          //Elastic scattering
          double deltime=f_step*temp/0.2;
	  vector<double> elscat = gen_particles(pp.x()/temp, pp.y()/temp, pp.z()/temp, id, deltime, wdk, wkcm, wx);
          
	  if (elscat[0]==1) {
            if (abs(id)==4) cout << "c scattering" << endl;
	    cout << "Elastic scattering!!" << endl;
	    //cout << "Elab= " << p[3] << " Efluid= " << pp.t() << endl;
	    //cout << " tof= " << tof << endl;
	    had_scattering = 1;

	    FourVector psumfluid = pp;
	    
	    FourVector pf(elscat[2]*temp,elscat[3]*temp,elscat[4]*temp,temp*std::sqrt(elscat[2]*elscat[2]+elscat[3]*elscat[3]+elscat[4]*elscat[4]));
            psumfluid -= pf;
	    pf = BoostBack(beta,pf);
	    int pf_id = elscat[5];
	    
	    FourVector kf(elscat[6]*temp,elscat[7]*temp,elscat[8]*temp,temp*std::sqrt(elscat[6]*elscat[6]+elscat[7]*elscat[7]+elscat[8]*elscat[8]));
            psumfluid += kf;
	    kf = BoostBack(beta,kf);
	    int kf_id = elscat[9];
	    
	    FourVector ff(elscat[10]*temp,elscat[11]*temp,elscat[12]*temp,temp*std::sqrt(elscat[10]*elscat[10]+elscat[11]*elscat[11]+elscat[12]*elscat[12]));
            psumfluid -= ff;
	    ff = BoostBack(beta,ff);
	    int ff_id = elscat[13];
	    
	    pp = BoostBack(beta,pp);
	    FourVector psum=pp;
	    psum += kf-pf-ff;

            if (fabs(psum.t())>0.0001 || fabs(psumfluid.t())>0.0001) {
	      cout << " Energy sum LAB= " << psum.t() << " Energy sum fluid= " << psumfluid.t() << endl;
	      cout << " Px sum LAB= " << psum.x() << " Px sum fluid= " << psumfluid.x() << endl;
	      cout << " Py sum LAB= " << psum.y() << " Py sum fluid= " << psumfluid.y() << endl;
	      cout << " Pz sum LAB= " << psum.z() << " Pz sum fluid= " << psumfluid.z() << endl;
	    }
	    /*
	    cout << " pin_x = " << pp.x() << " kf_x= " << kf.x() << " pf_x= " << pf.x() << " ff_x= " << ff.x() << endl;
	    cout << " pin_y = " << pp.y() << " kf_y= " << kf.y() << " pf_y= " << pf.y() << " ff_y= " << ff.y() << endl;
	    cout << " pin_z = " << pp.z() << " kf_z= " << kf.z() << " pf_z= " << pf.z() << " ff_z= " << ff.z() << endl;
	    cout << " pin_t = " << pp.t() << " kf_t= " << kf.t() << " pf_t= " << pf.t() << " ff_t= " << ff.t() << endl;
	    cout << " pin_id = " << id << " kf_id= " << kf_id << " pf_id= " << pf_id << " ff_id= " << ff_id << endl;
            */
	    vector<double> vpos{pos[0],pos[1],pos[2],pos[3]};
	    //Look for closest particle
	    double edif=100000000000.;
	    int iclose=-1000;
	    if (abs(id)!= 4) {
	      if (fabs(pp.t()-pf.t())<edif) edif=fabs(pp.t()-pf.t()), iclose=0;
	      if (fabs(pp.t()-ff.t())<edif) edif=fabs(pp.t()-ff.t()), iclose=1;
	    }
            else {
	      if (pf_id==id && ff_id==id) {
                cout << " Both charms!?" << endl;
		exit(1);
	      }
              if (pf_id==id) iclose=0;
	      else if (ff_id==id) iclose=1;
	      else {
                cout << " No match for charm!?" << endl;
	        cout << " Energy sum LAB= " << psum.t() << " Energu sum fluid= " << psumfluid.t() << endl;
	        cout << " pin_x = " << pp.x() << " kf_x= " << kf.x() << " pf_x= " << pf.x() << " ff_x= " << ff.x() << endl;
	        cout << " pin_y = " << pp.y() << " kf_y= " << kf.y() << " pf_y= " << pf.y() << " ff_y= " << ff.y() << endl;
	        cout << " pin_z = " << pp.z() << " kf_z= " << kf.z() << " pf_z= " << pf.z() << " ff_z= " << ff.z() << endl;
	        cout << " pin_t = " << pp.t() << " kf_t= " << kf.t() << " pf_t= " << pf.t() << " ff_t= " << ff.t() << endl;
	        cout << " pin_id = " << id << " kf_id= " << kf_id << " pf_id= " << pf_id << " ff_id= " << ff_id << endl;
		exit(1);
	      }
	    }
	    if (iclose==0) {
	      p[0]=pf.x();
	      p[1]=pf.y();
	      p[2]=pf.z();
	      p[3]=pf.t();
              //if (id!=pf_id) cout << " DIF IDS WITH ENERGY CRITERION!" << endl;
              vector<double> vff{ff.x(),ff.y(),ff.z(),ff.t()};
	      new_particles.push_back ( Parton ( vff, 100000000., 0, 0, -1, -1, ff_id, "recoiler", 0, 0, false ) );
	      new_particles[new_particles.size()-1].vSetRi(vpos);
	    }
	    else if (iclose==1) {
	      p[0]=ff.x();
	      p[1]=ff.y();
	      p[2]=ff.z();
	      p[3]=ff.t();
              //if (id!=ff_id) cout << " DIF IDS WITH ENERGY CRITERION!" << endl;
              vector<double> vpf{pf.x(),pf.y(),pf.z(),pf.t()};
	      new_particles.push_back ( Parton ( vpf, 100000000., 0, 0, -1, -1, pf_id, "recoiler", 0, 0, false ) );
	      new_particles[new_particles.size()-1].vSetRi(vpos);
	    }
	    else {
              cout << "No match for normal particle" << endl;
	      exit(1);
	    }
           
	    //Angle in lab frame
	    double elang=acos((pf.x()*ff.x()+pf.y()*ff.y()+pf.z()*ff.z())/pf.t()/ff.t());
	    //cout << "Products Angle= " << elang << " Pf pt= " << pf.pt() << " Ff pt= " << ff.pt() << endl;

            if (iclose==0) {
              double oriangle=acos((pf.x()*pp.x()+pf.y()*pp.y()+pf.z()*pp.z())/pf.t()/pp.t());
	      //if (oriangle>0.2) cout << "Ori Angle= " << oriangle << " Pp pt= " << pp.pt() << " Pf Fin pt= " << pf.pt() << endl;
	    }
            if (iclose==1) {
              double oriangle=acos((ff.x()*pp.x()+ff.y()*pp.y()+ff.z()*pp.z())/ff.t()/pp.t());
	      //if (oriangle>0.2) cout << "Ori Angle= " << oriangle << " Pp pt= " << pp.pt() << " Ff Fin pt= " << ff.pt() << endl;
	    }

            vector<double> vkf{kf.x(),kf.y(),kf.z(),kf.t()};
	    new_particles.push_back ( Parton ( vkf, 100000000., 0, 0, -1, -1, kf_id, "hole", 0, 0, true ) );
	    new_particles[new_particles.size()-1].vSetRi(vpos);
	    //cout << " orient size for elastic= " << new_particles[new_particles.size()-1].orient().size() << endl;
	    p_prev = p;
	  }
	  else if (elscat[0]==-1) {
	    // put a check to redo here in case you get -1, with smaller timestep
            cout << " Elscat problem!" << endl;
	  }
	}
      #endif
	//Broadening
	//p_prev[3]=sqrt(p_prev[0]*p_prev[0]+p_prev[1]*p_prev[1]+p_prev[2]*p_prev[2]);
	if (w2==0.) {
          cout << " w2= " << w2 << endl;
	  cout << " p= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
	}
	double p_prev_mod=sqrt(p_prev[0]*p_prev[0]+p_prev[1]*p_prev[1]+p_prev[2]*p_prev[2]);
	double p_mod=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
	double scalpprev=(p[0]*p_prev[0]+p[1]*p_prev[1]+p[2]*p_prev[2])/p_mod/p_prev_mod;
	if (scalpprev>1.) scalpprev=1.;
	double angbro=acos(scalpprev);
	//cout << endl;
	//cout << " p= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
	//cout << " p_prev= " << p_prev[0] << " " << p_prev[1] << " " << p_prev[2] << " " << p_prev[3] << endl;
        //cout << "BEF BRO P prev Pt= " << sqrt(p_prev[0]*p_prev[0]+p_prev[1]*p_prev[1]) << " P now Pt= " << sqrt(p[0]*p[0]+p[1]*p[1]) << " angle = " << angbro << endl;
	double virtbef=p[3]*p[3]-p[0]*p[0]-p[1]*p[1]-p[2]*p[2];
	//cout << " virtbef= " << virtbef << endl;
	
	if (kappa!=0. && step!=0.) trans_kick(w,w2,v,p,temp,vscalw,lore,step,kappa,nr);
	
	double virtaft=p[3]*p[3]-p[0]*p[0]-p[1]*p[1]-p[2]*p[2];
	//cout << " virtaft= " << virtaft << endl;
	//p[3]=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
	for (int oo=0; oo<4; oo++) {
          wake_rec[oo]+=p_prev[oo]-p[oo];
	}
	scalpprev=(p[0]*p_prev[0]+p[1]*p_prev[1]+p[2]*p_prev[2])/p_mod/p_prev_mod;
	if (scalpprev>1.) scalpprev=1.;
	angbro=acos(scalpprev);
        //cout << " AFT BRO P prev Pt= " << sqrt(p_prev[0]*p_prev[0]+p_prev[1]*p_prev[1]) << " P now Pt= " << sqrt(p[0]*p[0]+p[1]*p[1]) << " angle = " << angbro << endl;
	if (angbro!=angbro) {
          for (unsigned a=0; a<4; a++) {
            cout << " p prev " << a << " = " << p_prev[a] << endl;
            cout << " p " << a << " = " << p[a] << endl;
	  }
	  //double virt_p = p[3]*p[3]-p[0]*p[0]-p[1]*p[1]-p[2]*p[2];
	  //double virt_pprev = p_prev[3]*p_prev[3]-p_prev[0]*p_prev[0]-p_prev[1]*p_prev[1]-p_prev[2]*p_prev[2];
	  //cout << std::setprecision(5) << " virt= " << virt_p << " virt_prev= " << virt_pprev << endl;
	}
	orient[0]=p[0]/p[3];
	orient[1]=p[1]/p[3];
	orient[2]=p[2]/p[3];
	//Strong coupling
	bool doquench=1;
	if (abs(id)==4 && p[3]<=charm_mass) doquench=0;
	if (abs(id)==5 && p[3]<=b_mass) doquench=0;
	
	p_prev=p;

        if (alpha!=0. && model==0 && doquench)
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
	if (alpha!=0. && model==1)
	{
	  double intpiece=CF*(step/0.2)*alpha*temp*temp*temp*(f_dist/0.2);
	  double quench=(p[3]-intpiece)/p[3];
	  p*=quench;
	}
	//Collisional
	if (alpha!=0. && model==2)
	{
	  double intpiece=CF*(step/0.2)*alpha*temp*temp;
	  double quench=(p[3]-intpiece)/p[3];
	  p*=quench;
	}
      } 
    }

    //Make sure charm or b is not totally quenched
    if (abs(id)==4 && p[3]<charm_mass) {
      p[3]=charm_mass;
      double pmod=sqrt(p_prev[0]*p_prev[0]+p_prev[1]*p_prev[1]+p_prev[2]*p_prev[2]);
      if (pmod==0.) {
        cout << "pmod 0, why exactly?" << endl;
	pmod=1.;
      }
      p[0]=p_prev[0]/pmod*p[3];
      p[1]=p_prev[1]/pmod*p[3];
      p[2]=p_prev[2]/pmod*p[3];
    }
    if (abs(id)==5 && p[3]<b_mass) {
      p[3]=b_mass;
      double pmod=sqrt(p_prev[0]*p_prev[0]+p_prev[1]*p_prev[1]+p_prev[2]*p_prev[2]);
      if (pmod==0.) {
        cout << "pmod 0, why exactly?" << endl;
	pmod=1.;
      }
      p[0]=p_prev[0]/pmod*p[3];
      p[1]=p_prev[1]/pmod*p[3];
      p[2]=p_prev[2]/pmod*p[3];
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
    
    //Fill wake
    for (int oo=0; oo<4; oo++) {
      wake_dep[oo]+=p_prev[oo]-p[oo];
    }

    //Fill J_mu
    //if (p[3]!=p_prev[3]) {
      //source_file << tau << " " << pos[0] << " " << pos[1] << " " << eta << " ";
      //source_file << -p[3]+p_prev[3] << " " << -p[0]+p_prev[0] << " " << -p[1]+p_prev[1] << " " << -p[2]+p_prev[2] << " ";
      //source_file << p_prev[3] << " " << p_prev[0] << " " << p_prev[1] << " " << p_prev[2] << " ";
      //source_file << vx << " " << vy << " " << vz  << endl;
    //}

  } while (marker==0);

#ifdef DO_ELASTIC
  gsl_integration_workspace_free(wdk);
  gsl_integration_workspace_free(wx);
  gsl_integration_workspace_free(wkcm);
#endif

  //compare wake_dep and pini-p;
  /*
  for (int oo=0; oo<4; oo++) {
    if (fabs(wake_rec[oo]+wake_dep[oo]-(p_ini[oo]-p[oo]))>0.0001) { 
      cout << " wake_dep " << wake_dep[oo] << " wake_rec= " << wake_rec[oo] << " p_ini= " << p_ini[oo] << " p = " << p[oo] << endl;
      cout << " dif " << oo << " " << wake_rec[oo]+wake_dep[oo]-(p_ini[oo]-p[oo]) << endl; 
      cout << " id= " << id << endl;
    }
  }*/

  double scalpprev=(orient[0]*o_in[0]+orient[1]*o_in[1]+orient[2]*o_in[2]);
  if (scalpprev>1.) scalpprev=1.;
  double angbro=acos(scalpprev);
  //cout << " angbro final= " << angbro << endl;
  //cout << endl;

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

/*
void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp, double vscalw, double lore, double step, double kappa, numrand &nr)
{
  if (vscalw==1.) return; //Since it implies Ef=0
  
  vector<double> e1 = vec_prod(w,v);
  double Ne1=normalise(e1);
  if (Ne1==0.) {
    double b=0.5;
    double c=0.2;
    double a=(-b*w[1]-c*w[2])/w[0];
    e1={a,b,c,0.};
    Ne1=normalise(e1);
  }

  vector<double> l = vec_prod(w,e1)/sqrt(w2);

  double uscalW=lore*(1.-vscalw);
  double uscall=lore*(-v[0]*l[0]-v[1]*l[1]-v[2]*l[2]);
  double W2=1.-w2;

  vector<double> Wp = w-v*lore*W2/uscalW;

  double Nalpha=-uscall*uscalW/(pow(uscalW,2.)-W2);
  if (Nalpha!=Nalpha || __builtin_isinf(Nalpha)) return; //isinf replaced with __builtin_isinf to fix error

  vector<double> e2;
  if (lore==1.) {
    e2 = vec_prod(e1,w);
  } else {
    e2 = l+Wp*Nalpha;
  }
  double Ne2=normalise(e2);

  double Ef=p[3]*lore*(1.-vscalw);
  double wf2=1.-W2/pow(lore*(1.-vscalw),2.);
  double DelQ2=kappa*pow(temp,3.)*lore*(1.-vscalw)*step*5.;

  double qfac=0.;
  double phi=2.*3.141592654*nr.rando();
  if (Ef*sqrt(wf2)>0.)
  {
    qfac = sqrt(-1.*log(nr.rando()))*sqrt(DelQ2); //Box-Muller method
    if (qfac>Ef*sqrt(wf2)) {
      qfac=Ef*sqrt(wf2)-0.00000001;
    }
  }

  double qbeta;
  if (wf2>0.) qbeta=sqrt(1.-qfac*qfac/Ef/Ef/wf2)-1.;
  else qbeta=0.;
  if (qbeta!=qbeta) {
    cout << " qbeta= " << std::setprecision(6) << qbeta << " qfac= " << qfac << " Ef= " << Ef << " wf2= " << wf2 << " cutoff= " << Ef*sqrt(wf2) << endl;
    cout << " W2= " << std::setprecision(6) << W2 << " lore= " << lore << " vscalw= " << vscalw << endl;
  }

  vector<double> e = e1*cos(phi)+e2*sin(phi);
  vector<double> Wt = (w-v*uscalW*lore)/lore/(1.-vscalw);

  p+=Wt*qbeta*Ef+e*qfac;
  if (p[3]!=p[3] || p[2]!=p[2] || p[1]!=p[1] || p[0]!=p[0]) {
    cout << " p in bro= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
    cout << " qbeta= " << qbeta << " qfac= " << qfac << endl;
    cout << " Wt= " << Wt[0] << " " << Wt[1] << " " << Wt[2] << " " << Wt[3] << endl;
    cout << " e1= " << e1[0] << " " << e1[1] << " " << e1[2] << " " << e1[3] << endl;
    cout << " e2= " << e2[0] << " " << e2[1] << " " << e2[2] << " " << e2[3] << endl;
    cout << " Nalpha= " << Nalpha << endl;
    cout << " l= " << l[0] << " " << l[1] << " " << l[2] << " " << l[3] << endl;
    cout << " Wp= " << Wp[0] << " " << Wp[1] << " " << Wp[2] << " " << Wp[3] << endl;
    cout << " uscalW= " << uscalW << endl;
    cout << " uscall= " << uscall << endl;
    cout << " W2= " << W2 << endl;
    exit(1);
  }

}
*/

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

FourVector Boost( double b[3], FourVector p) {

  double betamod = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  double gamma = 1.0 / std::sqrt(1.0-betamod*betamod);
  double pscalb = p.x()*b[0] + p.y()*b[1] + p.z()*b[2];
  double spat = ( pscalb * gamma / (1.0 + gamma) - p.t() ) * gamma;

  double ptt = gamma * ( p.t() - pscalb );
  double ptx = p.x() + b[0] * spat;
  double pty = p.y() + b[1] * spat;
  double ptz = p.z() + b[2] * spat;

  return FourVector(ptx, pty, ptz, ptt);
}

FourVector BoostBack( double b[3], FourVector p) { //FIXME can be simplified BoostBack(v, p) == Boost(-v, p)

  double betamod = std::sqrt(b[0]*b[0] + b[1]*b[1] + b[2]*b[2]);
  double gamma = 1.0 / std::sqrt(1.0-betamod*betamod);
  double pscalb = p.x()*b[0] + p.y()*b[1] + p.z()*b[2];
  double spat = ( -pscalb * gamma / (1.0 + gamma) - p.t() ) * gamma;

  double bt = gamma * ( p.t() + pscalb );
  double bx = p.x() - b[0] * spat;
  double by = p.y() - b[1] * spat;
  double bz = p.z() - b[2] * spat;

  return FourVector(bx, by, bz, bt);
}

double Delta(FourVector v1, FourVector v2) {
  double deltaphi = fabs(v1.phi() - v2.phi());
  if (deltaphi>M_PI) {
    deltaphi = 2.*M_PI-std::max(v1.phi(),v2.phi())+std::min(v1.phi(),v2.phi());
  }
  //std::cout << " New deltaphi= " << deltaphi << std::endl;
  return std::sqrt(std::pow(v1.rapidity() - v2.rapidity(), 2.) +
                       std::pow(deltaphi, 2.));
}

