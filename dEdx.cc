#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <iostream>
#include <cmath>
#include "Random.h"

#include "vector_operators.h"

// #include <chrono> //FOR TESTING

using std::vector;
using namespace std;

void loss_rate(vector<double> &p, vector<double> &pos, double mass, double tof, int id, numrand &nr, double kappa, double alpha, double lambda, int tmethod, int mode, double &length, double &tlength, int new_eloss, double &fluidvx, double &fluidvy, double &fluidvz, bool &freezoutcrosser);
double normalise(vector<double> &p);
vector<double> vec_prod(vector<double> a, vector<double> b);
double gQ(double Del, numrand &nr, double cutoff);
void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp, double vscalw, double lore, double step, double kappa, numrand &nr);
void hqdiff(vector<double> &p, double temp, double step, double kappa, numrand &nr);
double gVx(double tau, double x, double y, double eta);
double gVy(double tau, double x, double y, double eta);
double gVz(double tau, double x, double y, double eta);
double gT(double tau, double x, double y, double eta);

void loss_rate(vector<double> &p, vector<double> &pos, double mass, double tof, int id, numrand &nr, double kappa, double alpha, double lambda, int tmethod, int mode, double &length, double &tlength, int new_eloss, double &fluidvx, double &fluidvy, double &fluidvz, bool &freezoutcrosser)
{
  if (p[3]*p[3]-p[0]*p[0]-p[1]*p[1]-p[2]*p[2]<0.) {
    // cout << " Initial state virtuality = " << p[3]*p[3]-p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-mass*mass << " mass= " << mass << " id=" << id << " p= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
    if (p[3]*p[3]-p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-mass*mass<-1.) { cout << "Large negative virtuality! " << p[3]*p[3]-p[0]*p[0]-p[1]*p[1]-p[2]*p[2]-mass*mass << " id=" << id << " p= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl; }
    p[3]=sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+mass*mass); //Force particle on shell if negative virtuality
    // double factor = sqrt((p[3]*p[3]-mass*mass)/(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]));
    // p[0]*=factor;
    // p[1]*=factor;
    // p[2]*=factor;
  }
  double Tc;
  if (tmethod==0) Tc=0.170;
  else Tc=0.145;

  double tot=pos[3]+tof;    //Final time
  //cout << " tot= " << tot << endl;

  double ei=p[3];      //Initial energy in lab frame
  // double etai;
  // if (pos[3]!=0. && pos[2]!=0.){
  //   etai =1./2.*log((pos[3]+pos[2])/(pos[3]-pos[2])); //Initial rapidity
  // } else {
  //   etai=0.;
  // }
  // double pti=sqrt(p[0]*p[0]+p[1]*p[1]);  //Initial transverse momentum

  double f_dist=0.;    //Traversed distance in Fluid Frame
  double l_dist=0.;    //Traversed distance in Lab Frame

  double CF;
  if (id==21) {
    if (mode==0) CF=pow(9./4.,1./3.);  //If gluon, color charge dependence is ratio of casimirs to power 1/3
    else CF=9./4.;
  }
  else CF=1.;

  int marker=0;    //If one, exit loop
  double step=0.02;  //Time step in LAB frame

  vector<double> w = p/p[3];  //4-velocity
  //cout << " x= " << pos[0] << " y= " << pos[1] << " z= " << pos[2] << " t= " << pos[3] << endl;
  //cout << " px= " << p[0] << " py= " << p[1] << " py= " << p[1] << " pz= " << p[2] << " en= " << p[3] << endl;
  //cout << " wx= " << w[0] << " wy= " << w[1] << " wz= " << w[2] << " we= " << w[3] << endl;
  
  // int nlight=0;
  // int nheavy=0;
  // int ncold=0;

  // auto start = std::chrono::high_resolution_clock::now(); //FOR TESTING

  bool washot=false;

  do {
    //Keep 4momentum before applying quenching this step
    vector<double> p_prev = p;

    if (pos[3]==tot) marker=1;
    if (pos[3]>tot) { cout << " Warning: Went beyond tot= " << tot << " t= " << pos[3] << endl; marker=1; } 
    //Proper time
    double tau=sqrt(pos[3]*pos[3]-pos[2]*pos[2]);
    if (tau!=tau) {
      tau=0.;
      cout << " TAU Not a number z= " << pos[2] << " t= " << pos[3] << " tot= " << tot << " marker= " << marker << " wz= " << w[2] << " en = " << p[3] << " pz= " << p[2] << " eta= " << 1./2.*log((pos[3]+pos[2])/(pos[3]-pos[2])) << endl;
      // exit(1);
    }

    //Rapidity
    double eta = 1./2.*log((pos[3]+pos[2])/(pos[3]-pos[2]));
    if (eta!=eta && tau>0.) {
      cout << " Eta is NaN= " << eta << " t= " << pos[3] << " z= " << pos[2] << endl;
      eta=0.;
    }
    if ((fabs(eta)>10. || eta != eta) && fabs(pos[2])>0. && fabs(pos[3])>0.) {
      cout << " Eta is too large to matter= " << eta << " t= " << pos[3] << " z= " << pos[2] << endl;
      marker=1;
    }

    int will_hot=0;  //Advance variable (to reach hot zones)
    if (tau>=0.6 && !(tau!=tau) && abs(eta) <= 10. && !(eta!=eta))  //Smooth profile starting time NOT NEEDED IN BRICK
    {
      vector<double> v;
      double vx=gVx(tau,pos[0],pos[1],eta); //0 FOR BRICK
      double vy=gVy(tau,pos[0],pos[1],eta); //0 FOR BRICK
      //**BOOST INVARIANT**
      double vz=pos[2]/pos[3]; //0 FOR BRICK
      double frap=atanh(vz);
      vx/=cosh(frap);
      vy/=cosh(frap);
      v.push_back(vx), v.push_back(vy), v.push_back(vz), v.push_back(1.); 

      
      
      double v2=pow(v[0],2.)+pow(v[1],2.)+pow(v[2],2.);
      double w2=pow(w[0],2.)+pow(w[1],2.)+pow(w[2],2.);
      if (w2>=1.){
        // cout << " W2 >= 1" << endl;
        w[0]/=w2;
        w[1]/=w2;
        w[2]/=w2;
        w2=1.;
      }
      if (v2>=1.){
        cout << " V2 >= 1" << endl;
        v[0]/=(v2/0.9999);
        v[1]/=(v2/0.9999);
        v[2]/=(v2/0.9999);
        v2=0.9999;
      }
      double vscalw=v[0]*w[0]+v[1]*w[1]+v[2]*w[2];
      // if (vscalw>=1.) vscalw=0.9999;
      if (vscalw>=1.) vscalw=0.99, cout << " VscalW >= 1 \n" << " v= " << v[0] << " " << v[1] << " " << v[2] << " w= " << w[0] << " " << w[1] << " " << w[2] << " virtuality= " << p[3]*p[3]-(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]+mass*mass) << " id=" << id << " en= " << p[3] << "rap=" << eta << endl;
      double lore=1./sqrt(1.-v2);

      //cout << " v2= " << v2 << " lore= " << lore << " vscalw= " << vscalw << endl;

      if (f_dist > l_dist ) {
        //cout << " f_dist= " << f_dist << " l_dist=" << l_dist << endl;
      }

      double temp=gT(tau,pos[0],pos[1],eta); //0 FOR BRICK
      // if (w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw) < 0.) cout << " WTFFFFFFFFFF temp= " << temp << " w2=" << w2 << " lore=" << lore << " v2=" << v2 << " vscalw=" << vscalw << " mass=" << mass << " total=" << w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw) << endl; 
      
      l_dist+=step;
      double f_step = step*sqrt(w2+lore*lore*(v2-2.*vscalw+vscalw*vscalw));
      if (f_step!=f_step) f_step=0.;
      f_dist+=f_step;

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
            // ncold += will_hot;
            break;
          }
        }
        if (will_hot==0)
        {
          pos+=w*(tot-pos[3]);
          marker=1;
          if (washot) freezoutcrosser=true;//(sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]) < 15);
          if (fluidvx==0. && fluidvy==0. && fluidvz==0.){
            fluidvx=v[0];
            fluidvy=v[1];
            fluidvz=v[2];
          }
        }
      } else {
        fluidvx=v[0];
        fluidvy=v[1];
        fluidvz=v[2];
      }

      //Now broad&quench
      if (p[3]>0. && temp>=Tc)
      {
        washot=true;
        if (mode==0) {
          length += f_step;
          tlength += temp / 0.2 * f_step;
        }
        else {
          length += f_step;
          tlength += pow(temp / 0.2,2.) * f_step;
        }
        //Broadening
        if (kappa!=0.) trans_kick(w,w2,v,p,temp,vscalw,lore,step,kappa,nr);
        // if (mass>temp && (new_eloss==2)) { //Heavy quark broadening
        //   trans_kick(w,w2,v,p,temp,vscalw,lore,step,3.141592/2.*sqrt(lambda),nr);
        // }
        //Strong coupling
        if (alpha!=0. && mode==0)
        {
          // Lorentz transform stuff
          double Efs=ei*lore*(1.-vscalw); //Initial energy in fluid frame
          vector<double> wfluid;
          vector<double> pfluid;
          double coef = (lore-1.0)*(p[0]*v[0]+p[1]*v[1]+p[2]*v[2])/v2;
          pfluid.push_back(p[0] + v[0]*coef - lore*p[3]*v[0]);
          pfluid.push_back(p[1] + v[1]*coef - lore*p[3]*v[1]);
          pfluid.push_back(p[2] + v[2]*coef - lore*p[3]*v[2]);
          pfluid.push_back(lore*( p[3] - (v[0]*p[0]+v[1]*p[1]+v[2]*p[2]) ));
          wfluid = pfluid/pfluid[3];
          double fluidstep=step*lore*(1.-vscalw);

          // cout << " p= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
          // cout << " pfluid= " << pfluid[0] << " " << pfluid[1] << " " << pfluid[2] << " " << pfluid[3] << endl;
          // cout << " fluidstep= " << fluidstep << " v= " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << " w= " << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << endl;

          // Light energyloss quantities
          double tstop=0.2*pow(Efs,1./3.)/(2.*pow(temp,4./3.)*alpha)/CF;
          double beta;
          if (f_dist==0) {
            beta=9999999999.;
            cout << " f_dist=0! " << " temp= " << temp << " mass= " << mass << " Efs= " << Efs << " beta= " << beta << " tstop= " << tstop << " f_dist= " << f_dist << endl;
          } else {
            beta=tstop/f_dist;
          }

          if (mass>temp && (new_eloss==1 || new_eloss==2) && (abs(id)==4 || abs(id)==5 || abs(id)==6)) { //Heavy quark energy loss            
            double etaD=3.141592/2*sqrt(lambda)*temp*temp/mass;
            double altintpiece=etaD*(pfluid[0]*pfluid[0]+pfluid[1]*pfluid[1]+pfluid[2]*pfluid[2])/pfluid[3]*(fluidstep*5);
            if (std::isinf(altintpiece) && std::signbit(altintpiece)) {
              cout << " altintpiece is inf! " << " etaD= " << etaD << " p= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << " vscalw= " << vscalw << " mass= " << mass << endl;
              altintpiece=0.;
            }

            bool treatasheavy=false;
            double intpiece;
            if (beta <= 1.) 
            {
              treatasheavy=true;
              // intpiece=999999.*fabs(altintpiece);// Just making sure it is larger so its not used. Could this cause overflow?
            } else {
              intpiece=Efs*4./(3.141592)*(1./(beta*tstop*sqrt(beta*beta-1.)))*fluidstep;
              if (intpiece>altintpiece) treatasheavy=true;
            }

            
            if (pfluid[3]-intpiece<mass) treatasheavy=true;

            if (treatasheavy)
            {
              pfluid[3]-=altintpiece;
              pfluid[0]-=etaD*pfluid[0]*(fluidstep*5);
              pfluid[1]-=etaD*pfluid[1]*(fluidstep*5);
              pfluid[2]-=etaD*pfluid[2]*(fluidstep*5);

              if (new_eloss==2) hqdiff(pfluid,temp,fluidstep,3.141592*sqrt(lambda),nr);

              if (p[3]!=p[3] || p[2]!=p[2] || p[1]!=p[1] || p[0]!=p[0]) {
                cout << " NaN Momenta! altintpiece" << endl;
                cout << " p in bro= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
                cout << " intpiece= " << intpiece << " altintpiece= " << altintpiece << endl;
                cout << " w= " << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << endl;
                cout << " v= " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << endl;
                cout << " lore= " << lore << " vscalw= " << vscalw << endl;
                cout << " temp= " << temp << " mass= " << mass << " Efs= " << Efs << " beta= " << beta << " tstop= " << tstop << " f_dist= " << f_dist << endl;
              }
            } else {
              double p2=pfluid[0]*pfluid[0]+pfluid[1]*pfluid[1]+pfluid[2]*pfluid[2];

              pfluid[0]-=pfluid[3]/p2*intpiece*pfluid[0];
              pfluid[1]-=pfluid[3]/p2*intpiece*pfluid[1];
              pfluid[2]-=pfluid[3]/p2*intpiece*pfluid[2];
              pfluid[3]-=intpiece;

              if (p[3]!=p[3] || p[2]!=p[2] || p[1]!=p[1] || p[0]!=p[0]) {
                cout << " NaN Momenta! intpiece" << endl;
                cout << " p in bro= " << p[0] << " " << p[1] << " " << p[2] << " " << p[3] << endl;
                cout << " intpiece= " << intpiece << " altintpiece= " << altintpiece << endl;
                cout << " w= " << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << endl;
                cout << " v= " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << endl;
                cout << " lore= " << lore << " vscalw= " << vscalw << endl;
                cout << " temp= " << temp << " mass= " << mass << " Efs= " << Efs << " beta= " << beta << " tstop= " << tstop << " f_dist= " << f_dist << endl;
              }
            }

            if (pfluid[3]<mass){
              cout << " pfluid[3] < mass! " << pfluid[3] << " mass= " << mass << " id=" << id << " eta= " << eta << " temp= " << temp << " Efs= " << Efs << " beta= " << beta << " tstop= " << tstop << " f_dist= " << f_dist << " intpiece= " << intpiece << " altintpiece= " << altintpiece << " treatasheavy= " << treatasheavy << endl;
              pfluid[3]=sqrt(pfluid[0]*pfluid[0]+pfluid[1]*pfluid[1]+pfluid[2]*pfluid[2]+mass*mass); // Force particle on shell if negative virtuality
              // exit(1);
            }

          } else if (new_eloss == 0 || mass <= temp || !(abs(id)==4 || abs(id)==5 || abs(id)==6)) //massless case
          {
            // double Efs=ei*lore*(1.-vscalw);
            // double tstop=0.2*pow(Efs,1./3.)/(2.*pow(temp,4./3.)*alpha)/CF;
            // double beta=tstop/f_dist;
            if (beta>1.)
            {
            double intpiece=Efs*fluidstep*4./(3.141592)*(1./(beta*tstop*sqrt(beta*beta-1.)));
            // //Update 4momentum
            double quench=(pfluid[3]-intpiece)/pfluid[3];
            pfluid*=quench;// SHOULD THIS ALSO INCLUDE THE SMALL MASS LIKE ABOVE??


              if (pfluid[3]!=pfluid[3] || pfluid[2]!=pfluid[2] || pfluid[1]!=pfluid[1] || pfluid[0]!=pfluid[0]) {
                cout << " NaN Momenta! intpiece" << endl;
                cout << " p in bro= " << pfluid[0] << " " << pfluid[1] << " " << pfluid[2] << " " << pfluid[3] << endl;
                cout << " intpiece= " << intpiece << endl;
                cout << " w= " << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << endl;
                cout << " v= " << v[0] << " " << v[1] << " " << v[2] << " " << v[3] << endl;
                cout << " lore= " << lore << " vscalw= " << vscalw << endl;
                cout << " temp= " << temp << " mass= " << mass << " Efs= " << Efs << " beta= " << beta << " tstop= " << tstop << " f_dist= " << f_dist << endl;
              }
            } else
            {
              pfluid[3]=0.;
            }
              if (pfluid[3]<0.) pfluid[3]=0.;
          } else {
            cout << " Unknown energy loss mode! " << new_eloss << endl;
            exit(1);
          }
          
          
          // Lorentz transform back
          coef = (lore-1.0)*(pfluid[0]*v[0]+pfluid[1]*v[1]+pfluid[2]*v[2])/v2;
          p[0]=pfluid[0] + v[0]*coef + lore*pfluid[3]*v[0];
          p[1]=pfluid[1] + v[1]*coef + lore*pfluid[3]*v[1];
          p[2]=pfluid[2] + v[2]*coef + lore*pfluid[3]*v[2];
          p[3]=lore*( pfluid[3] + (v[0]*pfluid[0]+v[1]*pfluid[1]+v[2]*pfluid[2]) );
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

    //If parton gets totally quenched, exit 
    if (p[3]<=0.)
    {
      marker=1;
      for (unsigned int i=0; i<3; i++) p[i]=0.;
      p[3]=0.;
    }
    else {
      //Manually protect very soft particles from getting kicks that yield velocities greater than 1
      // for (unsigned int i=0; i<3; i++) { if (p[i]>p[3]) p[i]=0.99999*p[3]; } THIS GETS DONE BELOW IN A SLIGHTLY MORE ROBUST WAY
      //Update kinematical quantities, with the possibility of advancing to hot regions
      w=p/p[3];
      if (w!=w) {
        cout << "w=" << w[0] << " " << w[1] << " " << w[2] << " " << w[3] << " E=" << p[3] << " p= " << p[0] << " " << p[1] << " " << p[2] << " " << endl;
        // exit(1);
      }
      double w2=w[0]*w[0]+w[1]*w[1]+w[2]*w[2];
      if (w2>=1.) {
        p[0]*=0.9999/sqrt(w2);
        p[1]*=0.9999/sqrt(w2);
        p[2]*=0.9999/sqrt(w2);

        w[0]=p[0]/p[3];
        w[1]=p[1]/p[3];
        w[2]=p[2]/p[3];
        // exit(1);
      }
      double tstep=max(double(will_hot),1.)*step;
      if (pos[3]+tstep>tot)
      {
        tstep=tot-pos[3];
      }
      if (marker!=1) pos+=w*tstep;
    }
    
  } while (marker==0);

  // if (p[3]>ei) { // This should happen sometimes, but if it happens too often, or if it happens with a large energy gain, it is a problem
  //   cout << "Gained energy! ei= " << ei << " ef= " << p[3] << " mass= " << mass << " id=" << id << " rap= " <<  1./2.*log((pos[3]+pos[2])/(pos[3]-pos[2])) << " rapi= " << etai << " pt= " << sqrt(p[0]*p[0]+p[1]*p[1]) << " pti= " << pti << endl;
  // }
  return;
}

void trans_kick(vector<double> w, double w2, vector<double> v, vector<double> &p, double temp, double vscalw, double lore, double step, double kappa, numrand &nr)
{
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
  // double Ne2=normalise(e2);

  double Ef=p[3]*lore*(1.-vscalw);
  double wf2=1.-W2/pow(lore*(1.-vscalw),2.);
  double DelQ2=kappa*pow(temp,3.)*lore*(1.-vscalw)*step*5.;

  double qfac=0.;
  double phi=2.*3.141592*nr.rando();
  if (Ef*sqrt(wf2)>0.)
  {
    qfac = sqrt(-2.*log(nr.rando()))*sqrt(DelQ2); //Box-Muller method
    if (qfac>Ef*sqrt(wf2)) qfac=Ef*sqrt(wf2);
  }

  double qbeta;
  if (wf2>0.) qbeta=sqrt(1.-qfac*qfac/Ef/Ef/wf2)-1.;
  else qbeta=0.;
  if (qbeta!=qbeta) {
    cout << " qbeta= " << std::setprecision(6) << qbeta << " qfac= " << qfac << " Ef= " << Ef << " wf2= " << wf2 << " cutoff= " << Ef*sqrt(wf2) << endl;
    cout << " W2= " << std::setprecision(6) << W2 << " lore= " << lore << " vscalw= " << vscalw << endl;
  }

  vector<double> Wt = (w-v*uscalW*lore)/lore/(1.-vscalw);
  vector<double> e = e1*qfac*cos(phi)+e2*qfac*sin(phi);
  
  p+=Wt*qbeta*Ef+ ;
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
  }
}

void hqdiff(vector<double> &p, double temp, double step, double kappa, numrand &nr)
{
    // broadening width
    double quarkv2 = (p[0]*p[0] + p[1]*p[1] + p[2]*p[2])/(p[3]*p[3]);
    if (quarkv2 >= 1.) {
        cout << "Quark velocity squared is greater than 1! " << quarkv2 << endl;
        quarkv2 = 0.9999;
    }
    double quarkgamma = 1.0 / sqrt(1.0 - quarkv2);
    if (quarkgamma != quarkgamma) {
        cout << "Quark gamma is NaN! " << quarkgamma << endl;
        quarkgamma = 1.0;
    }
    double DelQ2 = kappa*quarkgamma*pow(temp,3.)*(step*5.);
    
    // magnitudes of broadening using box-muller
    double qfac1 = sqrt(-2.*log(nr.rando()))*sqrt(DelQ2)*sin(2.*3.141592*nr.rando());
    double qfac2 = sqrt(-2.*log(nr.rando()))*sqrt(DelQ2)*sin(2.*3.141592*nr.rando());
    double qfac3 = sqrt(-2.*log(nr.rando()))*sqrt(DelQ2)*sin(2.*3.141592*nr.rando());

    // update momentum
    p[3] += sqrt(p[3]*p[3] + 2*(p[0]*qfac1 + p[1]*qfac2 + p[2]*qfac3) + qfac1*qfac1 + qfac2*qfac2 + qfac3*qfac3) - p[3];
    p[0] += qfac1;
    p[1] += qfac2;
    p[2] += qfac3;
}
