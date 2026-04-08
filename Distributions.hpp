// I hate c++,I really do. https://xkcd.com/353/
#include "helperfunctions.hpp"
#include <random>
#include <vector>
#include <complex>
#include "read_tables.hpp"

//requires gnuplot and gsl

//This header file contains useful functions for the Hybrid model large angle scattering
//Please input all parameters with units (i.e. pin) normalized by the temperature at the point,and all results (i.e. k,p,etc.) will come out normalized by temp.
//i.e. useful functions take arguments of t*T,x,pin_x/T,x,pin_y/T,x,pin_z/T,and pinPDG
//Please specify a pmin and a theta min below

//formulas are renormalized from the paper such that the integrals are with respect to k/T and p/T,and all variables are assumed dimensionless with temp normalization. (When T=1,the formulas agree)

double G=-1;//sign in the denominator of the FD/BE distributions. Bosons dont mind being all in the same state,so theres the minus sign in the denom. Knowledge!
double Q=+1;//yo same
//initialize needed values

double Nc=3;
double Nf=3;
double CF=(Nc*Nc-1)/(2*Nc);
double CA=Nc;
double dA=(Nc*Nc-1);
double dF=Nc;
double gs=2.25; //the strong coupling constant for the amplitudes 
double gsmd=1.5; //different gs for debye mass 
double vg=2*(Nc*Nc-1);
double vq=2*Nc;

//parameters
//double amD2=10*gsmd*gsmd/3*(Nc+Nf/2);
//double acutoff = 4.444444;
//double acutoff = 4.;
double acutoff = 10.;
double amD2=acutoff*gs*gs/3*(Nc+Nf/2);

double very_small_number = pow(10.,-50.);

int max_n[7]={0};

//^^^^^^^^^^^^^^^^^^^^^^^^


std::default_random_engine generator;
std::uniform_real_distribution<double> distcon(0.0,1.0);
std::uniform_int_distribution<> distint(0,1);
std::uniform_int_distribution<> distintNf(1,Nf);

double kmin_(double kcm, double x, double pin)
{
    if (x==0||kcm==0) return inf;
    return amD2/4/kcm/x;
}

double p_(double kcm,double x,double pin)
{
    if (x==0||kcm==0) return inf;
    return kmin_(kcm,x,pin)+pin-kcm;
}
double theta_(double kcm,double x,double pin)
{
    return acos(1-amD2/x/2/pin/p_(kcm,x,pin));
}

double getThermalFactor(double kcm, double x, double pin, double particleType) {
    // particleType: Q = +1 for fermion, G = -1 for boson
    double p = p_(kcm, x, pin);

    if (particleType == Q) { // Fermion
        return (1 - 1/(exp(p) + 1));
    } else if (particleType == G) { // Boson
        return (1 + 1/(exp(p) - 1));
    } else {
        std::cout << "Invalid particle type. Use Q for fermion and G for boson." << std::endl;
        exit(EXIT_FAILURE);
    }
}

double k_(double dk, double kcm,double x,double pin)
{
    return dk+kmin_(kcm,x,pin);
}
double f_(double kcm,double x,double pin)
{
    return pin+kmin_(kcm,x,pin);
}
double s_(double kcm,double x,double pin)
{
    return pin-kcm;
}

double q_(double kcm,double x,double pin)
{
    return kcm+kmin_(kcm,x,pin);
}
double kcmpm(double pm,double x,double pin)
{
    if (abs(pm)!=1)
    {
        std::cout<<"pm must=+-1"<<std::endl;
        exit(EXIT_FAILURE);
    }
    return (pin+pm*sqrt(pin*pin-amD2*(1+x)))/2./(1+x);
}

double dkpm(double pm,double kcm,double x,double pin)
{
    double f=f_(kcm,x,pin);
    double s=s_(kcm,x,pin);
    if (abs(pm) !=1)
    {
        std::cout<<"pm must=+-1"<<std::endl;
        exit(EXIT_FAILURE);
    }
    double out =pow(sqrt(s/f*(1+1./x))+pm,2)*f*x;
    if (out!=out)
    {
        std::cout<<"dkpm failed "<<pm<<"  "<<kcm<<"   "<< x<<"  "<<pin<<std::endl;
        exit(EXIT_FAILURE);
    }
    return out ;
}
double Cu_(double dk,double kcm,double x,double pin)
{
    double p=p_(kcm,x,pin);
    double s=s_(kcm,x,pin);
    double q=q_(kcm,x,pin);
    double out= ((p+pin)*dk/q/q+s/q);
    if (out!=out)
    {
        std::cout<<"Cu failed "<<dk<<"   "<<kcm<<"   "<< x<<"  "<<pin<<std::endl;
        exit(EXIT_FAILURE);
    }
    return out ;


    //return 4*kcm*(pin-kcm)*x/(amD2+4*kcm*kcm*x)+dk*(amD2-4*kcm*x*(kcm-2*pin))*4*kcm*x/pow(amD2+4*kcm*kcm*x,2);//(pin+p) * (2 * k+pin-p)/2/pow(q,2)-.5;
}
double C_(double dk,double kcm,double x,double pin)
{
    return Cu_(dk,kcm,x,pin)+1;
}
double D_(double dk,double kcm,double x,double pin)
{
    double f=f_(kcm,x,pin);
    double s=s_(kcm,x,pin);
    double q=q_(kcm,x,pin);
    if (dk==inf&&q!=inf) return inf;
    double out= 2*sqrt(f*s*dk*(dk+q))/pow(q,2);
    if (out!=out)
    {
        std::cout<<"D failed "<<f<<"   "<< s<<"  "<<dk<<std::endl;
        exit(EXIT_FAILURE);
    }
    return out ;
    //return 8*kcm*x*sqrt(dk*(pin-kcm)*(amD2+4*kcm*pin*x)*(amD2+4*kcm*(dk+kcm)*x))/pow(amD2+4*kcm*kcm*x,2);//sqrt((4 * pin * p+tt) * (4 * k * (k+(pin-p))+tt))/2/q^2;
}

double phiu_(double dk,double kcm,double x,double pin)
{
    if (dk==inf) return 0;
    // std::cout<<dk<<"   "<<kcm<<"   "<<x<<"   "<<pin<<"   "<<(Cu_(dk,kcm,x,pin)-x)<<" || "<<D_(dk,kcm,x,pin)<<std::endl;
    double tCu = Cu_(dk,kcm,x,pin);
    double tD = D_(dk,kcm,x,pin);
    double frac=(tCu-x)/tD;
    //std::cout << " frac= " << frac << std::endl;
    if (frac<-1) return M_PI;
    if (frac>1) return 0;
    if (tD == 0.) {
      if (tCu-x >= 0.) return 0;
      else return M_PI;
    }
    //double frac=(Cu_(dk,kcm,x,pin)-x)/D_(dk,kcm,x,pin);
    //if (frac<-1) return M_PI;
    //if (frac>1) return 0;

    return acos(frac);
}

double P_(double kcm,double x,double pin)
{
    double s=s_(kcm,x,pin);
    double p=p_(kcm,x,pin);
    double P=s+(p+pin)*x;
    if (P!=P||P==inf)
    {
        std::cout<<"P failed "<<kcm<<"   "<< x<<"  "<<pin<<std::endl;
        exit(EXIT_FAILURE);
    }
    return P;
}

double M_(double kcm,double x,double pin)
{
    if (x==inf) return inf;
    double f=f_(kcm,x,pin);
    double s=s_(kcm,x,pin);
    double M=2*sqrt(f*s*x*(x+1));
    if (M!=M)
    {
        std::cout<<"M failed "<<kcm<<"   "<< x<<"  "<<pin<<std::endl;
        exit(EXIT_FAILURE);
    }
    // std::cout<<M<<std::endl;
    return M;
}

double psi_(double dk,double kcm,double x,double pin)
{
    if (dk==inf) return 0;
    double P=P_(kcm,x,pin);
    double M=M_(kcm,x,pin);
    if (dk>P+M) return 0;
    if (dk<P-M) return M_PI;
    double psi=acos((dk-P)/M);
    if (P==inf&&dk==inf&&M==inf) return 0;
    if (psi!=psi)
    {
        std::cout<<"psi failed "<<dk<<"   "<< P<<"  "<<M<<std::endl;
        exit(EXIT_FAILURE);
    }
    return psi;
}

double Psi_(double dk,double kcm,double x,double pin)
{
    if (dk==inf) return 0;
    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double phiu=phiu_(dk,kcm,x,pin);
    if (phiu==0&&(C!=0||D!=0)) return 0;
    // std::cout<<C<<"  "<<D<<"   "<<phiu<<std::endl;
    if (phiu==M_PI) return M_PI;
    if (C/D<1)
    {
        std::cout << "Psi <1" << std::endl;
        exit(EXIT_FAILURE);
    }
    return h(phiu,C/D);
}

double Psiso_(double dk,double kcm,double x,double pin)
{
    if (dk==inf) return 0;
    double psi=psi_(dk,kcm,x,pin);
    //if (psi==M_PI) return 0;
    if (fabs(psi-M_PI)<small_number) return 0;
    double Ps=P_(kcm,x,pin)-s_(kcm,x,pin);
    double M=M_(kcm,x,pin);
    if (Ps>=M)
    {
        return 2*atan((Ps-M)/sqrt(Ps*Ps-M*M)*tan(psi/2))/sqrt(Ps*Ps-M*M);
    }
    if (M>Ps)
    {
        double mine=(M-Ps)/sqrt(M*M-Ps*Ps)*tan(psi/2);
        return log(abs((1+mine)/(1-mine)))/sqrt(M*M-Ps*Ps);
    }
    std::cout << "Psiso not complete" << std::endl;
    exit(EXIT_FAILURE);
}

double Psiu_(double dk,double kcm,double x,double pin)
{
    if (dk==inf) return 0;
    double Cu=Cu_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double phiu=phiu_(dk,kcm,x,pin);
    if (phiu==0&&(Cu!=0||D!=0)) return 0;
    if (phiu==M_PI) return M_PI;
    if (Cu/D<1-pow(10.,-7.))
    {
        std::cout << "Psiu <1" << std::endl;
	cout << " Cu= " << Cu << " D= " << D << " dk= " << dk << " kcm= " << kcm << " x= " << x << " pin= " << pin << " phiu= " << phiu << std::endl;
        exit(EXIT_FAILURE);
    }
    return h(phiu,Cu/D);
}

// double dkpsi_(double kcm, double x, double pin,double psi)
// {
//     double mine=P_(kcm,x,pin)+M_(kcm,x,pin)*cos(psi);
//     if (mine!=mine)
//     {
//         std::cout<<"dkpsi failed "<<kcm<<"   "<< x<<"  "<<psi<<std::endl;
//         exit(EXIT_FAILURE);
//     }

//     return mine;
// }



double front0_(double dk,double kcm,double x,double pin,double g,double d, double particleA)
{
    double kmin=kmin_(kcm,x,pin);
    double thermalFactor = getThermalFactor(kcm, x, pin, particleA);
    return thermalFactor * exp(-dk-kmin)*kmin/x/(1+g*exp(-dk-kmin))/ (1+d*exp(-dk-kcm));
}
double frontdk_(double dk,double kcm,double x,double pin,double g,double d, double particleA, double n)
{
    // if (x>1) return 0;
    double kmin=kmin_(kcm,x,pin);
    // std::cout<<kcm<<"   "<<x<<"  "<<pin<<std::endl;
    // if (abs(exp(-n*dk-n*kmin)*powdiff(-d*exp(-kcm+kmin),-g,n)*kmin/x)<pow(10,-50)) std::cout<<"frontdk "<<exp(-n*dk-n*kmin)*powdiff(-d*exp(-kcm+kmin),-g,n)*kmin/x<<std::endl;
    double thermalFactor = getThermalFactor(kcm, x, pin, particleA);
    return thermalFactor * exp(-kmin)*powdiff(-d*exp(-kcm),-g*exp(-kmin),n)*kmin/x;
}

double frontkcm_(double kcm, double x, double pin, double g, double d, double n)
{
    double kmin=kmin_(kcm,x,pin);
    return exp(-kmin)*powdiff(-d*exp(-kcm),-g*exp(-kmin),n);
}



double Vu(double phi,double dk,double kcm,double x,double pin,double middle=1)
{
    //this is just the base level reduced domain of integration
    double Cu=Cu_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double kp=dkpm(1,kcm,x,pin);
    double km=dkpm(-1,kcm,x,pin);
    double s=s_(kcm,x,pin);
    double q=q_(kcm,x,pin);
    return (dk>kp)||(((cos(phi)<(Cu-x)/D)&&dk<kp&&dk>km)&&middle)||(dk<km&&s>q*x);
}


// double dphidpsi(double dk,double kcm, double x, double pin)
// {
//     double q=q_(kcm,x,pin);
//     double s=s_(kcm,x,pin);
//     double f=f_(kcm,x,pin);
//     return -.5*((x*q-s)/dk+(x*q+f)/(dk+q));
// }

// double dPhidpsi(double dk,double kcm, double x, double pin)
// {
//     double q=q_(kcm,x,pin);
//     double s=s_(kcm,x,pin);
//     double f=f_(kcm,x,pin);
//     return -.5*(2+(x*q-s)/dk-(x*q+f)/(dk+q));
// }

// double dPhiu0dpsi(double dk,double kcm, double x, double pin)
// {
//     double q=q_(kcm,x,pin);
//     double s=s_(kcm,x,pin);
//     double f=f_(kcm,x,pin);
//     return -.5*(2-(x*q-s)/dk+(x*q+f)/(dk+q));
// }



// double Phi(double dk,double kcm, double x, double pin)
// {
//     double C=C_(dk,kcm,x,pin);
//     double D=D_(dk,kcm,x,pin);
//     double phiu=phiu_(dk,kcm,x,pin);
//     return h(phiu,C/D);
// }

// double Phiu(double dk,double kcm, double x, double pin)
// {
//     double Cu=Cu_(dk,kcm,x,pin);
//     double D=D_(dk,kcm,x,pin);
//     double phiu=phiu_(dk,kcm,x,pin);
//     return h(phiu,Cu/D);
// }

double phigamint(double kcm, double x, double pin, double dk, double n, double i,gsl_integration_workspace *w, double u = 0 )
{
    double psi=psi_(dk,kcm,x,pin);
    if (fabs(psi-M_PI)<small_number) return 0;
    //if (psi==M_PI) return 0;
    if(!w)
    {
        std::cout << "set the workspace in phigamint  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    
    // if (sn!=0||cs!=0)
    // {
    //     std::cout << "phigamint hand integrate pls" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
    //std::cout << "running phigamint " <<psi<< std::endl;
    // if (sn>0)
    // {
    //     // std::cout << "finished phigamint" << std::endl;
    //     //return 2*phigamint(kcm,x,pin,psi,n,i,u,cs+1,sn-2)-phigamint(kcm,x,pin,psi,n,i,u,cs+2,sn-2);
    //     return phigamint(kcm,x,pin,psi,n,i,u,cs,sn-2)-phigamint(kcm,x,pin,psi,n,i,u,cs+2,sn-2);
    // }
    // if (sn!=0)
    // {
    //     std::cout << "phigamint hand integrate pls" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
    //int from 0 to psi of (dphi dpsi, dPhi dpsi,dPhi_u dpsi)*gamma(j,n(K))*cos(x)^cs*sin(x)^sn if j=1,2,3
    double q=q_(kcm,x,pin);
    double P=P_(kcm,x,pin);
    double M=M_(kcm,x,pin);
    double s=s_(kcm,x,pin);
    double f=f_(kcm,x,pin);
    double p=p_(kcm,x,pin);

    double corr=1;
    if (i==1)
    {
        double p1=(.5)*(expKint(P,M,psi,-1,0,n,w)*(x*q-s)+expKint(P,M,psi,-1,q,n,w)*(x*q+f))*corr;
        if (p1 == inf || p1 != p1)
        {
            std::cout << "phigamint failure "<<   i<< std::endl;
            exit(EXIT_FAILURE);
        }
        // std::cout << "finished phigamint1" << std::endl;
        return p1;
    }
    
    if (i>1) 
    {
        // std::cout << "started phigamint "<< i<< std::endl;
        double p1=(i-1)*phigamint(kcm,x,pin,dk,n,i-1,w)+(.5)*(expKint(P,M,psi,i-1,0,n,w,0)*(x*q-s)+expKint(P,M,psi,i-1,0,n,w,q)*(x*q+f))*corr;
        // double p1=(1+n*(P+M))*(-.5)*(expKint(P,M,psi,-1,0,n,cs,sn)*(x*q-s)+expKint(P,M,psi,-1,q,n,cs,sn)*(x*q+f));
        // double p2=(-n*M)*(-.5)*(expKint(P,M,psi,-1,0,n,1+cs,sn)*(x*q-s)+expKint(P,M,psi,-1,q,n,1+cs,sn)*(x*q+f));
        // std::cout << "finished phigamint " <<i <<std::endl;
        if (p1 == inf || p1 != p1)
        {
            std::cout << "phigamint failure "<<   i<< std::endl;
            exit(EXIT_FAILURE);
        }
        return p1;
    }
    // if (i==3) 
    // {
    //     double p1=(n*n*(P+M)*(P+M)+2*n*(P+M)+2)*(-.5)*(expKint(P,M,psi,-1,0,n,cs,sn)*(x*q-s)+expKint(P,M,psi,-1,q,n,cs,sn)*(x*q+f));
    //     double p2=(-2*M*(1+n*(P+M)))*(-.5)*(expKint(P,M,psi,-1,0,n,1+cs,sn)*(x*q-s)+expKint(P,M,psi,-1,q,n,1+cs,sn)*(x*q+f));
    //     double p3=(M*M)*(-.5)*(expKint(P,M,psi,-1,0,n,2+cs,sn)*(x*q-s)+expKint(P,M,psi,-1,q,n,2+cs,sn)*(x*q+f));
    //     // std::cout << "finished phigamint" << std::endl;
    //     return p1+p2+p3;
    // }
    double a=q/2.+pow(-1,u)*(p+pin)/2.;
    if (i==0)
    {
         //std::cout << "finished phigamint3 "<< expintK(P,M,psi,1,0,0,a,n)<<"   "<<pow(-1,u)*expintK(P,M,psi,1,-1,0,a,n)*(x*q-s)<<"   "<<pow(-1,u)*expintK(P,M,psi,1,-1,q,a,n)*(x*q+f)<< std::endl;
        return (.5)*(2*expintK(P,M,psi,1,0,0,a,n,w)+pow(-1,u)*expintK(P,M,psi,1,-1,0,a,n,w)*(x*q-s)-pow(-1,u)*expintK(P,M,psi,1,-1,q,a,n,w)*(x*q+f))*corr;
    }
    if (i<0)
    {
        // std::cout << "entering <0" << std::endl;
        // double p1=(-.5)/i*(                       2*expintK(P,M,psi,i+1,0,0,a,n,cs,sn)+pow(-1,u)*expintK(P,M,psi,i+1,-1,0,a,n,cs,sn)*(x*q-s)-pow(-1,u)*expintK(P,M,psi,i+1,-1,q,a,n,cs,sn)*(x*q+f));
        
        double p1=1./i*phigamint(kcm,x,pin,dk,n,i+1,w,u);

        
        // std::cout << "p1 done" << std::endl;
        // std::cout << expKint(P,M,psi,i,a,n,cs,sn,inf) << std::endl;
        double p2=-pow(n,i)*(.5)/i*(2*expKint(P,M,psi,i,a,n,w,inf)+pow(-1,u)*(x*q-s)*expKint(P,M,psi,i,a,n,w,0)-pow(-1,u)*(x*q+f)*expKint(P,M,psi,i,a,n,w,q));
        // std::cout << "finished phigamint4" << std::endl;
        return p1+p2*corr;

    }
    std::cout<<"phigamint not complete"<<std::endl;
    exit(EXIT_FAILURE);
}



double m1dk(double dk,double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *w )
{
    
    if (dk!=dk||kcm!=kcm||x!=x) std::cout<<"m1dk bad input "<<dk<<"   "<<kcm<<"   "<<x<<std::endl;
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    if (fabs(x)<0.0001) return 0;
    if(!w)
    {
        std::cout << "set the workspace in m1dk "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    //std::cout << "starting m1dk" << std::endl;
    double q=q_(kcm,x,pin);
    double f=f_(kcm,x,pin);
    double s=s_(kcm,x,pin);
    double p = p_(kcm, x, pin);
    double km=dkpm(-1,kcm,x,pin);
    double phiu=phiu_(dk,kcm,x,pin);
    double psi=psi_(dk,kcm,x,pin);

    double P=P_(kcm,x,pin);
    double M=M_(kcm,x,pin);

    
    
    double frontdk=frontdk_(dk,kcm,x,pin,g,d, particleA, n);
    double m1dk_13=-((pow(pin+p,2)+2*f*s)/pow(q,4)*Gamma(3,n*dk)/pow(n,3)+(2*f/pow(q,3)*(pin+p+s)*Gamma(2,n*dk)/pow(n,2)+pow(f/q,2)*Gamma(1,n*dk)/n))*(1-phiu/M_PI);
    double m1dk_2=(pow(pin+p,2)+2*f*s)/pow(q,4)/pow(n,3)*phigamint(kcm,x,pin,dk,n,3,w)+2*f/pow(q,3)/pow(n,2)*(p+pin+s)*phigamint(kcm,x,pin,dk,n,2,w)+f*f/pow(q,2)/n*phigamint(kcm,x,pin,dk,n,1,w)+pow(M,2)/2/q*((3*f/q+3*P*(p+pin)/q/q+1+x)*expKint(P,M,psi,0,0,n,w,inf,0,2)+3*(p+pin)/q/q*M*expKint(P,M,psi,0,0,n,w,inf,1,2));
    double out =frontdk*(m1dk_13-m1dk_2);

    if (out == inf || out != out)
    {
        std::cout << "m1dk failure " <<expKint(P,M,psi,0,0,n,w,inf,0,2)<<"  "<<phigamint(kcm,x,pin,dk,n,1,w)<<"  "<<s<<"    "<<p<<"   "<<km<<"   "<< std::endl;
	std::cout << dk <<"  "<<kcm<<"  "<<x<<"    "<<pin<< std::endl;
        exit(EXIT_FAILURE);
    }
    return out; 
}
double m2dk(double dk,double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *w ) 
{
    if (dk!=dk||kcm!=kcm||x!=x) std::cout<<"m2dk bad input"<<std::endl;
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    if (fabs(x)<0.0001) return 0;
    if(!w)
    {
        std::cout << "set the workspace in m2dk  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    //std::cout << "starting m2dk" << std::endl;
    double q=q_(kcm,x,pin);
    double f=f_(kcm,x,pin);
    // double s=s_(kcm,x,pin);
    double p = p_(kcm, x, pin);
    //double km=dkpm(-1,kcm,x,pin);
    //double s=s_(kcm,x,pin);
    double phiu=phiu_(dk,kcm,x,pin);
    double psi=psi_(dk,kcm,x,pin);

    double P=P_(kcm,x,pin);
    double M=M_(kcm,x,pin);
    
    double frontdk=frontdk_(dk,kcm,x,pin,g,d, particleA, n);
    double m2dk_13=-((pin+p)/pow(q*n,2)*Gamma(2,n*dk)+f/q/n*Gamma(1,n*dk))*(1-phiu/M_PI);
    double m2dk_2=((pin+p)/pow(q*n,2)*phigamint(kcm,x,pin,dk,n,2,w)+f/q/n*phigamint(kcm,x,pin,dk,n,1,w))+pow(M,2)/q*expKint(P,M,psi,0,0,n,w,inf,0,2);
    double out =frontdk*(m2dk_13-m2dk_2);

    if (out == inf || out != out)
    {
        std::cout << "m2dk failure " <<phiu<<"  "<<dk<<"  "<<kcm<<"    "<<x<<"   "<<pin<<"   "<< std::endl;
        exit(EXIT_FAILURE);
    }
    return out; 
}
double m3dk(double dk,double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *w ) 
{
      

    if (dk!=dk||kcm!=kcm||x!=x) std::cout<<"m3dk bad input"<<std::endl;
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    if (fabs(x)<0.0001) return 0;
    if(!w)
    {
        std::cout << "set the workspace in m3dk  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    //std::cout << "starting m3dk" << std::endl;

    double km=dkpm(-1,kcm,x,pin);
    //double q=q_(kcm,x,pin);
    //double s=s_(kcm,x,pin);

    double phiu=phiu_(dk,kcm,x,pin);

    //double psi=psi_(dk,kcm,x,pin);

    double frontdk=frontdk_(dk,kcm,x,pin,g,d, particleA, n);
    double m3dk=-1./n*(1-phiu/M_PI)*exp(-n*dk);
    if (dk>km) m3dk-=1./n*phigamint(kcm,x,pin,dk,n,1,w);
    // std::cout<<"m3dk "<<phiu<<"    "<<psi<<"     "<<std::endl;
    double out =frontdk*m3dk;

    if (out == inf || out != out)
    {
        std::cout << "m3dk failure " <<phiu<<"  "<<dk<<"  "<<kcm<<"    "<<x<<"   "<<pin<<"   "<< std::endl;
        exit(EXIT_FAILURE);
    }
    return out; 
}
double m4dk(double dk,double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *w ) 
{
   
    if (dk!=dk||kcm!=kcm||x!=x) std::cout<<"m4dk bad input"<<std::endl;
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    if (fabs(x)<0.0001) return 0;
    if(!w)
    {
        std::cout << "set the workspace in m4dk  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    //std::cout << "starting m4dk" << std::endl;
    double q=q_(kcm,x,pin);
    double f=f_(kcm,x,pin);
    //double s=s_(kcm,x,pin);
    //double km=dkpm(-1,kcm,x,pin);
    //double psi=psi_(dk,kcm,x,pin);
    double Psi=Psi_(dk,kcm,x,pin);


    //std::cout<<Psi<<"   "<<kcm<<"   "<<x<<"   "<<dk<<"    "<<pin<<std::endl;
    //exit(EXIT_FAILURE);
    double frontdk=frontdk_(dk,kcm,x,pin,g,d, particleA, n);
    double m4dk_13=-q*expGamma(0,n*(dk+f))*(1-Psi/M_PI)*exp(-n*dk);
    double m4dk_2=q*phigamint(kcm,x,pin,dk,n,0,w);

    double out =frontdk*(m4dk_13-m4dk_2);

    if (out == inf || out != out)
    {
        std::cout << "m4dk failure " << frontdk << "  " << m4dk_2 << std::endl;
        exit(EXIT_FAILURE);
    }
    return out; 
}
double m5dk(double dk,double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *w ) 
{
    
    if (dk!=dk||kcm!=kcm||x!=x) std::cout<<"m5dk bad input"<<std::endl;
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    if (fabs(x)<0.0001) return 0;
    if(!w)
    {
        std::cout << "set the workspace in m5dk  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    //std::cout << "starting m5dk" << std::endl;
    double q=q_(kcm,x,pin);
    double f=f_(kcm,x,pin);
    double p=p_(kcm,x,pin);
    double s=s_(kcm,x,pin);
    //double km=dkpm(-1,kcm,x,pin);
    double psi=psi_(dk,kcm,x,pin);
    double P=P_(kcm,x,pin);
    double M=M_(kcm,x,pin);
    
    double frontdk=frontdk_(dk,kcm,x,pin,g,d, particleA, n);
    double m5dk_13=-q*((pin+p)*n*expGamma(-1,n*(dk+f))-2*f*s*n*n*expGamma(-2,n*(dk+f)))*(1-Psi_(dk,kcm,x,pin)/M_PI)*exp(-n*dk);
    double m5dk_2=q*((pin+p)*n*phigamint(kcm,x,pin,dk,n,-1,w)-2*f*s*n*n*phigamint(kcm,x,pin,dk,n,-2,w))-q*M*M/(1+x)*expKint(P,M,psi,-2,f,n,w,inf,0,2);
    double out =frontdk*(m5dk_13-m5dk_2);

    if (out == inf || out != out)
    {
        std::cout << "m5dk failure " << frontdk<< "  " << out << std::endl;
        exit(EXIT_FAILURE);
    }
    return out; 
}
double m6dk(double dk,double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *w ) 
{
    
    if (dk!=dk||kcm!=kcm||x!=x) std::cout<<"m6dk bad input"<<std::endl;
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    if (fabs(x)<0.0001) return 0;

    if(!w)
    {
        std::cout << "set the workspace in m6dk  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    //std::cout << "starting m6dk" << std::endl;
    double q=q_(kcm,x,pin);
    double s=s_(kcm,x,pin);
    //double km=dkpm(-1,kcm,x,pin);

    //double psi=psi_(dk,kcm,x,pin);
    double Psiu=Psiu_(dk,kcm,x,pin);

    double frontdk=frontdk_(dk,kcm,x,pin,g,d, particleA, n);
    double m6dk_13=(q*expGamma(0,n*(dk-s))*sgn(-(dk-s)))*(1-Psiu/M_PI)*exp(-n*dk);
    //if (abs(dk-s)<.5) std::cout << 1- Psiu/M_PI<<"  "<< std::endl;
    double m6dk_2=q*phigamint(kcm,x,pin,dk,n,0,w,1);
    double out =frontdk*(m6dk_13-m6dk_2);

    if (out == inf || out != out)
    {
        std::cout << "m6dk failure " << frontdk << "  " << out << std::endl;
        exit(EXIT_FAILURE);
    }
    return out; 
}

struct special_params { double P;double M;double s;double n;};

double special_f (double x, void * p)
{
    struct special_params * params = (struct special_params *)p;
    double P = (params->P);
    double M = (params->M);
    double n = (params->n);
    double s = (params->s);
    double K = P+M*cos(x);

    return  (exp(-n*K))/(K-s)/M_PI;
}

double special(double dk, double kcm, double x, double pin, double n,gsl_integration_workspace *w )
{
    
    double psi=psi_(dk,kcm,x,pin);
    //if (psi==M_PI) return 0;
    if (fabs(psi-M_PI)<small_number) return 0;
    
    double s=s_(kcm,x,pin);
    double P=P_(kcm,x,pin);
    double M=M_(kcm,x,pin);
    long double sum=0;
    long double inc = 1000;

    double psis=M_PI;
    // if (M==0) return exp(-n*P)/(P-s)*(psi-M_PI)/M_PI;
    if (P-s>M) return expKint(P,M,psi,0,0,n,w,-s);

    gsl_function F;
    struct special_params params = {P, M, s, n};

    F.function = &special_f;
    F.params = &params;
    double result1=0;
    double result2=0;
    double abserr=0;
    double epsabs=abs_res_dk;
    double epsrel=rel_res_dk;
    double limit=100000;
    double pole=acos(-(P-s)/M);
    gsl_set_error_handler_off();
    if (!w)
    {
        std::cout << "set the workspace in special  " << std::endl;
        exit(EXIT_FAILURE);
       
    }
    double myeps=.00000001;

    gsl_integration_qags(&F, M_PI, std::max(pole + myeps, psi), epsabs, epsrel, limit, w, &result1, &abserr);
    gsl_integration_qags(&F, pole - myeps, std::min(psi, pole - myeps), epsabs, epsrel, limit, w, &result2, &abserr);

    if (result1+result2!=result1+result2||abs(result1+result2)==inf) 
    {
        std::cout << "special failure1  "<< result1<<"   "<<result2<<std::endl;
        exit(EXIT_FAILURE);
    }
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result1+result2;


   

    // std::cout<<F21(.5,.5,1.5-10,(M+P-s)/2/M)<<std::endl;
    // if (psi==0)
    // {
    //     std::cout<<P<<"  "<<s<<"  "<<M<<std::endl;
    //     for (int r = 1; abs(inc)*100000000000000000/ep >= abs(sum); r++) // 100000000000
    //     {
    //            std::cout<<r<<"  "<<inc<<"  "<<sum<<"  "<<F21(.5,.5,.5+r,(M+P-s)/2/M)<<"   "<<F21(.5,.5,.5+r,(M+P-s)/2/M)<<"  "<<pow(M+P-s,r)<<"  "<<pow((P-s-M)/2/M,r)<<std::endl;
    //         //  inc=M_PI*xofact(-n*(M+P-s),r,-n*s)*pow((M+P-s),-1)*F21(1-r,.5,1.5-r,(P-s-M)/(M+P-s))*choose(r-1.5,r-1);
    //         //  inc=M_PI*xofact((M-P+s)*n,r,-n*s)*F21(.5,.5,1.5-r,(M+P-s)/2/M)/2/M*pow((M-P+s)/2/M,.5)*choose(r-.5,.5)/2;
    //         //  if (psi> acos((P-s)/M))
    //         // else inc= ;
    //         inc=-xofact(-n*(M+P-s),r,-n*s)*2/choose(r-.5,r-1)/sqrt((M*M-(P-s)*(P-s)))*(sqrt((M+s-P)/2/M)*F21(.5,.5,.5+r,(M+P-s)/2/M)-pow((P-s-M)/2/M,r)*F21(.5,.5,.5+r,(M+P-s)/2/M));
            
    //         sum += inc;

    //     }
    //     if (inc!=inc||inc==inf) 
    //     {
    //         std::cout << "special failure1  "<< inc<<std::endl;
    //         exit(EXIT_FAILURE);
    //     }
    //     // std::cout<<" finishing special "<<std::endl;
    //     //std::cout<<P<<"  "<<M<<"   "<<s<<" "<<psi<<"  "<<sum<<"  "<<atan2(sqrt(1-a*a),-a)<<std::endl;
    //     return sum;

    // }

    // if (psi==0)
    // {
    //     //  std::cout<<P<<"  "<<M<<"   "<<s<<" "<<psi<<"  "<<sum<<"  "<<std::endl;
    //     double sq=sqrt(M*M-(P-s)*(P-s));
    //     for (int k = 0; abs(inc)*10000000/ep >= abs(sum); k++) // 100000000000
    //     {
    //         double sum1=0;
    //         for (int l = 0; l<=k; l++) sum1+=xofact(n*M/2.,l,-n*P/2)*sin(l*atan2(-sq,P-s)+k*atan2(sq,P-s));
    //         inc=sum1*xofact(n*M/2.,k,-n*P/2);
    //         //  if (psi> acos((P-s)/M))
    //         // else inc= ;
    //         //    std::cout<<k<<"  "<<inc<<"  "<<sum1<<"  "<<sum<<"   "<<M*n<<"  "<<psi<<std::endl;
    //         if (inc!=inc||inc==inf) 
    //         {
    //             std::cout << "special failure0  "<< inc<<std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    //         sum += inc;
    //     }
    //     // std::cout<<P<<"  "<<M<<"   "<<s<<" "<<psi<<"  "<<sum<<"  "<<std::endl;
    //      return 2*sum/sq*M_PI;

    // }
    
    // std::cout<<" starting special "<<std::endl;
    double a=(P-s)/M;

    for (int r = 1; abs(inc)*100000000/ep >= abs(sum); r++) // 100000000000
    {
         inc=2*sin(r*atan2(sqrt(1-a*a),-a))*BesselI_scaled(r,n*M,psi)*exp(-n*P+n*M);
         //  if (psi> acos((P-s)/M))
         // else inc= ;
        //    std::cout<<r<<"  "<<inc<<"  "<<sum<<"  "<<BesselI_scaled(r,n*M,psi)<<"   "<<M*n<<"  "<<psi<<std::endl;
         sum += inc;

    }
    if (inc!=inc||abs(inc)==inf) 
    {
        std::cout << "special failure1  "<< inc<<std::endl;
        exit(EXIT_FAILURE);
    }
    // std::cout<<" finishing special "<<std::endl;
     //std::cout<<P<<"  "<<M<<"   "<<s<<" "<<psi<<"  "<<sum<<"  "<<atan2(sqrt(1-a*a),-a)<<std::endl;
    return sum/sqrt(M*M-(P-s)*(P-s));

    std::cout << "special failure got to end  " << inc << std::endl;
    exit(EXIT_FAILURE);

    ////

    //////
    // long double mysum=0;
    // double mypsi=psis;
    // double K=0;
    // double mystep=(psi-M_PI)/psi_res;
    // for (int l = 0; l < psi_res; l++)
    // {
    //     K=P+M*cos(mypsi);
    //     inc=(exp(-n*K)-exp(-n*s))/(K-s)*mystep;
    //     mypsi+=mystep;
    //     mysum+=inc;
    
    //     // std::cout<<inc<<std::endl;
    //     if (inc == inf || inc != inc)
    //         {
    //             // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //             // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //             // std::cout << "expKintincfailr "<<std::endl;
    //             std::cout << "specialint failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // //std::cout << mystep << std::endl;
    // //std::cout << inc << std::endl;
    // return mysum;
    /////////

    //
    // long double mysum=0;
    // double psia=psi;
    // double mypsi=M_PI;
    // double mystep=(M_PI-psia)/psi_res;
    // double K=P+M*cos(mypsi);
    // for (int l = 0; l < psi_res; l++)
    // {
    //     K=P+M*cos(mypsi);
    //     inc=-(exp(-n*K)-exp(-n*s))/(K-s)*mystep;
    //     mysum+=inc;
    //     mypsi-=mystep;
    //     //std::cout<<inc<<std::endl;
    //     if (inc == inf || inc != inc)
    //         {
    //             // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //             // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //             // std::cout << "expKintincfailr "<<std::endl;
    //             std::cout << "special failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // //std::cout << mystep << std::endl;
    // //std::cout << inc << std::endl;
    // return mysum;
    //




    //double N=0;
    /////
    // for (int r = 1; (abs(inc)*1000000000./ep >=abs(sum));r++)
    // {
    //     inc=0;
    //     if (P-s<M)
    //     {
    //         psis = acos(-(P - s) / M);
    //         double psil=std::max(psi,psis);
    //         // double psim = psi;
    //         // double psis = psi;
    //         // if (P - s - M * cos(psi) > 0)
    //         // {
    //         //     psis = acos(-(P - s) / M);
    //         //     psim = std::max(psi, psis);
    //         // }
    //         //  if((P-s -M)/(M + P-s)>1) std::cout << (P-s -M)/(M + P-s)<<std::endl;
            
    //         if (psis!=psil) inc+=-(1./r)*sqrt(1/(M*M - (P-s)*(P-s)))*F1(r, .5, .5, 1 + r, (P-s + M *cos(psil))/(-M + P-s), (P-s + M *cos(psil))/(M + P-s))*xofact(-n*(P-s + M*cos(psil)),r,-n*s);
    //         inc-=-(1./r)*sqrt(1/(M*M - (P-s)*(P-s)))*F1(r, .5, .5, 1 + r, 1, (P-s -M)/(M + P-s))*xofact(-n*(P-s -M),r,-n*s);

    //         // inc = 0;
    //         // for (int j = 0; j <= r - 1; j++)
    //         // {
    //         //     // inc+= -pow(P - s, r - 1.) *cos(psi)* choose(r - 1, j)*pow(M*cos(psi)/(P-s),j)*F21(.5, (1 + j)*.5, (3 + j)*.5, cos(psi)*cos(psi)) / (1. + j);
    //         //     inc += -pow(P - s - M, -1.) * choose(r - 1, j) * pow(M * (1 + cos(psi)) / (P - s - M), j) / (.5 + j) * cos(psi / 2) * F21(.5, .5 + j, 1.5 + j, cos(psi / 2) * cos(psi / 2)) * exp(-n * s) * xofact(-n * (P - s - M), r);
    //         // }
    //         //  std::cout << inc << "  " << exp(-n * s) * xofact(-n * (P - s - M), r) << "   " << r << std::endl;
    //     }

    //     if (inc!=inc||inc==inf) 
    //     {
    //         std::cout << "special failure1  "<< inc<<std::endl;
    //         exit(EXIT_FAILURE);
    //     }



    //     if (psi<psis) 
    //     {
    //         ////////// 
    //         // double mine=std::min(1.,M*(1-cos(psis))/(M+P-s));
    //         // double sum1 = 0;
    //         // double inc1 = 1000;
    //         // for (int j= 0; (abs(inc1)*100000./ep >= abs(sum1));j++)
    //         // {
    //         //     inc1=0;
    //         //     // std::cout << " special F1" << 2*M/(P-s+M)*sin(psi/2)*sin(psi/2) << "   "<<mine<<"  "<<std::endl;
    //         //     if (psi!=0) inc1=-sin(psi/2)/(j+.5)*xofact(n*M*(1-cos(psi)),j)*F1(j+.5,.5,1-r,1.5+j,sin(psi/2)*sin(psi/2),2*M/(P-s+M)*sin(psi/2)*sin(psi/2));
    //         //     inc1-=-sin(psis/2)/(j+.5)*xofact(n*M*(1-cos(psis)),j)*F1(j+.5,.5,1-r,1.5+j,sin(psis/2)*sin(psis/2),mine);

    //         //     // inc1=-sin(psi/2)/(P-s+M)/(j+.5)*xofact(n*M*(1-cos(psi)),j)*xofact(n*(P-s+M),r,-n*(P-s+M))*F1(j+.5,.5,1-r,1.5+j,sin(psi/2)*sin(psi/2),2*M/(P-s+M)*sin(psi/2)*sin(psi/2));
    //         //     // inc1-=-sin(psis/2)/(P-s+M)/(j+.5)*xofact(n*M*(1-cos(psis)),j)*xofact(n*(P-s+M),r,-n*(P-s+M))*F1(j+.5,.5,1-r,1.5+j,sin(psis/2)*sin(psis/2),2*M/(P-s+M)*sin(psis/2)*sin(psis/2));
    //         //     // inc1-=-1/(P-s+M)/(j+.5)*xofact(2*n*M,j)*xofact(n*(P-s+M),r,-n*(P-s+M))*F1(j+.5,.5,1-r,1.5+j,1,2*M/(P-s+M));
    //         //     sum1+=inc1;
    //         //     //  std::cout << " special " << inc1 << "   "<<sum1<<"  "<<n*M*(1-cos(psis))<<"  "<<xofact(n*M*(1-cos(psis)),j)<<std::endl;
    //         //     if(j==5&&inc1==0&&sum1==0) break;
           
            
    //         // }
    //         // inc+=sum1*xofact(n*(P-s+M),r,-n*(P+M))/(P-s+M);
    //         /////////// 

    //         /////////
    //         if (r==1)
    //         {
    //             long double mysum=0;
    //             double mypsi=psis;
    //             double K=0;
    //             double mystep=(psi-psis)/psi_res;
    //             for (int l = 0; l < psi_res; l++)
    //             {
    //                 mypsi+=mystep;
    //                 K=P+M*cos(mypsi);
    //                 inc=(exp(-n*K)-exp(-n*s))/(K-s)*mystep;
    //                 mysum+=inc;
                    
    //                 // std::cout<<inc<<std::endl;
    //                 if (inc == inf || inc != inc)
    //                     {
    //                         // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //                         // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //                         // std::cout << "expKintincfailr "<<std::endl;
    //                         std::cout << "special2 failure" << std::endl;
    //                         exit(EXIT_FAILURE);
    //                     }
    //             }
    //             //std::cout << mystep << std::endl;
    //             //std::cout << inc << std::endl;
    //             inc+=mysum;
    //         }
    //         //////////////

    //         //inc+=xofact(-n*(P-s+M),r,-n*s)*(2*F1(.5,.5,-r,1.5,sin(x/2)*sin(x/2),M*(1-cos(x))/(M+P-s))*sin(x/2)-M_PI*F21(.5,-r,1,2*M/(P-s+M)));
    //         // std::cout<<"special1 "<<r<<"  "<<psi<<"  "<<psis<<std::endl;
    //         // std::cout<<"special11 "<<F1(.5,.5,-r,1.5,sin(psis/2)*sin(psis/2),1)<<std::endl;
    //         //  std::cout<<"special22 "<<F1(.5,.5,-r,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(M+P-s))<<std::endl;
    //         // std::cout<<"special2 "<<std::endl;


    //         // double mine=std::min(1.,M*(1-cos(psis))/(M+P-s));
    //         // if (psi!=0) inc+=xofact(-n*(P-s+M),r,-n*s)*(2*F1(.5,.5,1-r,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(M+P-s))*sin(psi/2))/(P-s+M);
    //         // inc-=xofact(-n*(P-s+M),r,-n*s)*(2*F1(.5,.5,1-r,1.5,sin(psis/2)*sin(psis/2),mine)*sin(psis/2))/(P-s+M);
            
    //         //std::cout<<"special3 "<<std::endl;
            


            
    //     }
        
    //     //std::cout << inc  << std::endl;
    //     if (inc!=inc||inc==inf) 
    //     {
    //         std::cout << "special failure2  "<< inc<<"   "<<xofact(-n*(P-s+M),r,-n*s)<<std::endl;
    //         exit(EXIT_FAILURE);
    //     }
        
    //     sum += inc;
    //     //N+=1;
    //     if(r==5&&inc==0&&sum==0) break;
    //     // std::cout << " special2 " << inc<<"  "<<sum << std::endl;
    // }
    // // std::cout << "special "<<std::endl;
    // return sum;
    /////
}

double m7dk(double dk,double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *w ) 
{
    
    if (dk!=dk||kcm!=kcm||x!=x) std::cout<<"m7dk bad input"<<std::endl;
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (fabs(x)<0.0001) return 0;

    if(!w)
    {
        std::cout << "set the workspace in m7dk  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    // std::cout<<dk<<"   "<<kcm<<"   "<<x<<"   "<<pin<<"    "<<g<<std::endl;
    //std::cout << "starting m7dk" << std::endl;
    double q=q_(kcm,x,pin);
    double f=f_(kcm,x,pin);
    double s=s_(kcm,x,pin);
    double Z=f+s;
    double M=M_(kcm,x,pin);
    double P=P_(kcm,x,pin);
    //double km=dkpm(-1,kcm,x,pin);
    double psi=psi_(dk,kcm,x,pin);
    // double dkpsi=dkpsi_(kcm,x,pin,psi);
    double Psiu=Psiu_(dk,kcm,x,pin);
    //if (psi!=0&&psi!=M_PI) std::cout<<Psi<<"   "<<psi<<"   "<<(s>q*x)<<"   "<<dk<<"    "<<dkpsi<<std::endl;
    
    double frontdk=frontdk_(dk,kcm,x,pin,g,d, particleA, n);
    double m7dk_13=q*(Z*n*expGamma(-1,n*(dk-s))+2*f*s*n*n*expGamma(-2,n*(dk-s)))*sgn(-(dk-s))*(1-Psiu/M_PI)*exp(-n*dk);
    
    if (m7dk_13!=m7dk_13) {
      std::cout << "m7dk13 q= " << q  << " dk= " << dk << " s= " << s << " Psiu= " << Psiu << " n= " << n << " f= " << f << " kcm= " << kcm << " x= " << x << " pin= " << pin << std::endl;
    }
    
    double mine=-1;
    double m7dk_2=q/M_PI*(M*sin(psi)*exp(-n*dk)*.5/x/(dk-s)*expEn(2,n*(dk-s)))*mine;//+n/2*(Z*Z-q*q*(1+x))*Psiso_(dk,kcm,x,pin)*exp(n*(-s))
    if (dk==s) m7dk_2 = pow(10.,50.);
     if (m7dk_2!=m7dk_2) 
     {
         std::cout << "m7dk2 failure 0 " <<Psiso_(dk,kcm,x,pin) << std::endl;
         exit(EXIT_FAILURE);
     }
    m7dk_2-=q*(f*s*n*n-Z*n)*phigamint(kcm,x,pin,dk,n,0,w,1)*mine;
     if (m7dk_2 != m7dk_2)
     {
         std::cout << "m7dk2 failure 1" << std::endl;
         exit(EXIT_FAILURE);
     }
    m7dk_2+=q*(1-Z*n/2)*((x*q+f)*expKint(P,M,psi,0,0,n,w,q)-(x*q-s)*expKint(P,M,psi,0,0,n,w,0))/2*mine;
     if (m7dk_2!=m7dk_2) 
     {
         std::cout << "m7dk2 failure 2"<< std::endl;
         exit(EXIT_FAILURE);
     }
    m7dk_2-=q*(expKint(P,M,psi,0,0,n,w)/x-n/2/x*expintK(P,M,psi,1,1,-s,-s,n,w)+Z*n/2*expintK(P,M,psi,1,0,0,-s,n,w))*mine;
     if (m7dk_2 != m7dk_2)
     {
         std::cout << "m7dk2 failure 3" << std::endl;
         exit(EXIT_FAILURE);
     }
    m7dk_2+=q*q*n/2*phigamint(kcm,x,pin,dk,n,1,w)*mine;// (x*q+f)*expKint(P,M,psi,0,0,n,q)+(x*q-s)*expKint(P,M,psi,0,0,n,0))*.5*;
     if (m7dk_2 != m7dk_2)
     {
         std::cout << "m7dk2 failure 4" << std::endl;
         exit(EXIT_FAILURE);
     }
    m7dk_2+=q*n/2*(Z*Z-q*q*(1+x))*special(dk,kcm,x,pin,n,w)*mine;
     if (m7dk_2!=m7dk_2) 
     {
         std::cout << "m7dk2 failure 5"<< " special= " << special(dk,kcm,x,pin,n,w) << std::endl;
	 std::cout << " (Z*Z-q*q*(1+x))= " << (Z*Z-q*q*(1+x)) << " Z= " << Z << " q= " << q << " x= " << x << " mine= " << mine << std::endl;
	 std::cout << " kcm= " << kcm << " pin= " << pin << " f= " << f << " s= " << s << std::endl;
         exit(EXIT_FAILURE);
     }

    if (m7dk_2!=m7dk_2) {
      std::cout << "m7dk2 q= " << q  << " dk= " << dk << " s= " << s << " Psiu= " << Psiu << " n= " << n << " f= " << f << " kcm= " << kcm << " x= " << x << " pin= " << pin << std::endl;
    }
    
    double out =frontdk*(m7dk_13-m7dk_2);

    if (out == inf || out != out)
    {
        std::cout << "m7dk failure " << frontdk << "  " << (m7dk_13-m7dk_2) << std::endl;
        exit(EXIT_FAILURE);
    }
    return out; 
}

double m1phi(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA)//bounds 0-2pi
{
    //  std::cout << "starting m1phi" << std::endl;
    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    double kp=dkpm(1,kcm,x,pin);
    double km=dkpm(-1,kcm,x,pin);
    double n=2;
    // std::cout << "C " << C << std::endl;
    // std::cout << "x " << x << std::endl;
    // std::cout << "D "<<D << std::endl;
    // std::cout << "CDx "<<(C-1)/D-x/D << std::endl;
    double phim=phiu_(dk,kcm,x,pin);
    double m1phi_2=(phi>phim)*(ABo2pi(n,C,D,phi)-ABo2pi(n,C,D,phim))+(phi>2*M_PI-phim)*(ABo2pi(n,C,D,2*M_PI-phim)-ABo2pi(n,C,D,phi));
    double m1phi_13=ABo2pi(n,C,D,phi);
    // std::cout << "finishing m1phi "<<m1phi_2<<"  "<< m1phi_13 <<"  " << phim<<" "<<D<<std::endl;
    // std::cout << "finishing m1phi "<<front0*(Vu(0,dk,kcm,x,pin)*m1phi_13+(kp>dk)*(dk>km)*m1phi_2) << std::endl;
    return front0*(Vu(0,dk,kcm,x,pin,0)*m1phi_13+(kp>dk)*(dk>km)*m1phi_2); 
}
double m2phi(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA) //bounds 0-2pi
{
    //  std::cout << "starting m2phi" << std::endl;
    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    double kp=dkpm(1,kcm,x,pin);
    double km=dkpm(-1,kcm,x,pin);
    double n=1;
    double phim=phiu_(dk,kcm,x,pin);
    double m2phi_2=(phi>phim)*(ABo2pi(n,C,D,phi)-ABo2pi(n,C,D,phim))+(phi>2*M_PI-phim)*(ABo2pi(n,C,D,2*M_PI-phim)-ABo2pi(n,C,D,phi));
    double m2phi_13=ABo2pi(n,C,D,phi);
    return front0*(Vu(0,dk,kcm,x,pin,0)*m2phi_13+(kp>dk)*(dk>km)*m2phi_2);
}
double m3phi(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA) //bounds 0-2pi
{
    //  std::cout << "starting m3phi" << std::endl;
    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    double kp=dkpm(1,kcm,x,pin);
    double km=dkpm(-1,kcm,x,pin);
    double n=0;
    double phim=phiu_(dk,kcm,x,pin);
    double m3phi_2=(phi>phim)*(ABo2pi(n,C,D,phi)-ABo2pi(n,C,D,phim))+(phi>2*M_PI-phim)*(ABo2pi(n,C,D,2*M_PI-phim)-ABo2pi(n,C,D,phi));
    double m3phi_13=ABo2pi(n,C,D,phi);
    return front0*(Vu(0,dk,kcm,x,pin,0)*m3phi_13+(kp>dk)*(dk>km)*m3phi_2);
}
double m4phi(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA) //bounds 0-2pi
{
    //  std::cout << "starting m4phi" << std::endl;
    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    double kp=dkpm(1,kcm,x,pin);
    double km=dkpm(-1,kcm,x,pin);
    double n=-1;
    double phim=phiu_(dk,kcm,x,pin);
    double m4phi_2=(phi>phim)*(ABo2pi(n,C,D,phi)-ABo2pi(n,C,D,phim))+(phi>2*M_PI-phim)*(ABo2pi(n,C,D,2*M_PI-phim)-ABo2pi(n,C,D,phi));
    double m4phi_13=ABo2pi(n,C,D,phi);
    // std::cout<<front0*(Vu(0,dk,kcm,x,pin)*m4phi_13+(kp>dk)*(dk>km)*m4phi_2)<<std::endl;//"     "<<phim<<"   "<<C<<"   "<<D<<"   "<<ABo2pi(n,C,D,2*M_PI-phim)<<"   "<<ABo2pi(n,C,D,phim)<<std::endl;
    return front0*(Vu(0,dk,kcm,x,pin,0)*m4phi_13+(kp>dk)*(dk>km)*m4phi_2);
}
double m5phi(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA) //bounds 0-2pi
{
    //  std::cout << "starting m5phi" << std::endl;
    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    double kp=dkpm(1,kcm,x,pin);
    double km=dkpm(-1,kcm,x,pin);
    double n=-2;
    double phim=phiu_(dk,kcm,x,pin);
    double m5phi_2=(phi>phim)*(ABo2pi(n,C,D,phi)-ABo2pi(n,C,D,phim))+(phi>2*M_PI-phim)*(ABo2pi(n,C,D,2*M_PI-phim)-ABo2pi(n,C,D,phi));
    double m5phi_13=ABo2pi(n,C,D,phi);
    return front0*(Vu(0,dk,kcm,x,pin,0)*m5phi_13+(kp>dk)*(dk>km)*m5phi_2);
}
double m6phi(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA)//bounds 0-2pi
{
    //  std::cout << "starting m6phi" << std::endl;
    double Cu=Cu_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    double kp=dkpm(1,kcm,x,pin);
    double km=dkpm(-1,kcm,x,pin);
    double n=-1;
    double phim=phiu_(dk,kcm,x,pin);
    // std::cout << "m6phi" << Cu<<"  "<<D<< std::endl;
    double m6phi_2=(phi>phim)*(ABo2pi(n,Cu,D,phi)-ABo2pi(n,Cu,D,phim))+(phi>2*M_PI-phim)*(ABo2pi(n,Cu,D,2*M_PI-phim)-ABo2pi(n,Cu,D,phi));
    double m6phi_13=ABo2pi(n,Cu,D,phi);
    // std::cout << "finished m6phi" << front0*(Vu(0,dk,kcm,x,pin)*m6phi_13+(kp>dk)*(dk>km)*m6phi_2)<< std::endl;
    
    return front0*(Vu(0,dk,kcm,x,pin,0)*m6phi_13+(kp>dk)*(dk>km)*m6phi_2);
}
double m7phi(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA) //bounds 0-2pi
{
    //  std::cout << "starting m7phi" << std::endl;
    double Cu=Cu_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);

    double kp=dkpm(1,kcm,x,pin);
    double km=dkpm(-1,kcm,x,pin);

    double n=-2;
    double phim=phiu_(dk,kcm,x,pin);
    double m7phi_2=(phi>phim)*(ABo2pi(n,Cu,D,phi)-ABo2pi(n,Cu,D,phim))+(phi>2*M_PI-phim)*(ABo2pi(n,Cu,D,2*M_PI-phim)-ABo2pi(n,Cu,D,phi));
    double m7phi_13=ABo2pi(n,Cu,D,phi);
    // std::cout << "finished m7phi" << front0*(Vu(0,dk,kcm,x,pin)*m7phi_13+(kp>dk)*(dk>km)*m7phi_2)<< std::endl;
    return front0*(Vu(0,dk,kcm,x,pin,0)*m7phi_13+(kp>dk)*(dk>km)*m7phi_2);
    
}

double m1base(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA)
{

    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    //  std::cout << "finished m1base" << front0*pow(C-D * cos(phi),2)/ (2 * M_PI) * Vu(phi,dk,kcm,x,pin)<< std::endl;
    // cout<<B<<"  "<<A<<std::endl;
    return front0*pow(C-D * cos(phi),2)/ (2 * M_PI) * Vu(phi,dk,kcm,x,pin);
}
double m2base(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA)
{

    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);

    return front0 * (C-D * cos(phi)) / (2 * M_PI)*Vu(phi,dk,kcm,x,pin);
}
double m3base(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA)
{
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    return front0/(2 * M_PI)*Vu(phi,dk,kcm,x,pin);
}
double m4base(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA)
{

    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    
    return front0 / (2 * M_PI) / (C-D * cos(phi))*Vu(phi,dk,kcm,x,pin);
}
double m5base(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA)
{
    double C=C_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    
    return front0 / (2 * M_PI) / pow(C-D * cos(phi),2)*Vu(phi,dk,kcm,x,pin);
}
double m6base(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA)
{

    double Cu=Cu_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    
    return front0 / (2 * M_PI) / (Cu-D * cos(phi))*Vu(phi,dk,kcm,x,pin);
}
double m7base(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA)
{

    double Cu=Cu_(dk,kcm,x,pin);
    double D=D_(dk,kcm,x,pin);
    double front0=front0_(dk,kcm,x,pin,g,d, particleA);
    // std::cout<<front0 / (2 * M_PI) / pow(Cu-D * cos(phi),2)*Vu(phi,dk,kcm,x,pin)<<std::endl;
    return front0 / (2 * M_PI) / pow(Cu-D * cos(phi),2)*Vu(phi,dk,kcm,x,pin);
}



struct makcm_params {double x;double pin;double g;double d; double particleA; double n; double a;gsl_integration_workspace *w;};

double makcm_f (double y, void * p)
{
    struct makcm_params * params = (struct makcm_params *)p;
    double x = (params->x);
    double pin = (params->pin);
    double g = (params->g);
    double d = (params->d);
    double particleA = (params->particleA);
    double n = (params->n);
    double a= (params->a);
    gsl_integration_workspace* w = (params->w);

    if (x>1) return 0;
    if (y==0) return 0;
    if (x==0) return 0;
    if (x!=x) std::cout << "help x has gone bad "<<std::endl;
    if(!w)
    {
        std::cout << "set the workspace in makcm "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 

    if (a==1) return m1dk(inf,y,x,pin,g,d, particleA, n,w)-m1dk(0,y,x,pin,g,d, particleA, n,w);
    if (a==2) return m2dk(inf,y,x,pin,g,d, particleA, n,w)-m2dk(0,y,x,pin,g,d, particleA, n,w);
    if (a==3) return m3dk(inf,y,x,pin,g,d, particleA, n,w)-m3dk(0,y,x,pin,g,d, particleA, n,w);
    if (a==4) return m4dk(inf,y,x,pin,g,d, particleA, n,w)-m4dk(0,y,x,pin,g,d, particleA, n,w);
    if (a==5) return m5dk(inf,y,x,pin,g,d, particleA, n,w)-m5dk(0,y,x,pin,g,d, particleA, n,w);
    if (a==6) return m6dk(inf,y,x,pin,g,d, particleA, n,w)-m6dk(0,y,x,pin,g,d, particleA, n,w);
    if (a==7) return m7dk(inf,y,x,pin,g,d, particleA, n,w)-m7dk(0,y,x,pin,g,d, particleA, n,w);
}

/*
double do_kcm_interp(int iX, double x, double pin, double kcm, int g, int d, int n)
{

  if (n>17) {
    std::cout << "iX= " << iX << " n= " << n << " is too large for current tabulation! EXITING" << std::endl;
    return 0.;
    exit(1);
  }

  int iA = (n-1)*4;
  if (g == -1 && d == -1) iA += 0;
  if (g == -1 && d == 1) iA += 1;
  if (g == 1 && d == -1) iA += 2;
  if (g == 1 && d == 1) iA += 3;

  if (x > 1) {
    std::cout << " x > 1 in mIkcm, x = " << x << std::endl;
    exit(1);
    //x=1.; //CHECK
  }

  if (kcm > pin) {
    std::cout << " kcm > pin in mIkcm, kcm = " << kcm << " pin = " << pin << std::endl;
    exit(1);
  }
  
  int ux = int(x/step_x);
  double dux = (x - double(ux)*step_x)/step_x;
  int fux = 1;
  if (ux == SIZE_X-1) fux = 0, dux = 0.;

  int up = int((pin - pin_min)/pin_step);
  double dup = (pin - pin_min - double(up)*pin_step)/pin_step;
  int fup = 1;
  if (up == SIZE_PIN-1) fup = 0, dup = 0.;

  double step_kcm = pin_vals[up]/double(nbins_kcm); //step of kcm depends on pin... is this the best I can do?
  int uk = int(kcm/step_kcm);
  double duk = (kcm - double(uk)*step_kcm)/step_kcm;
  int fuk = 1;
  if (uk == SIZE_KCM-1) fuk = 0, duk = 0.;

  double result = 0.;

  result += mkcm_table[iX-1][iA][up][ux][uk]*(1.-dup)*(1.-dux)*(1.-duk);
  result += mkcm_table[iX-1][iA][up+fup][ux][uk]*dup*(1.-dux)*(1.-duk);
  result += mkcm_table[iX-1][iA][up][ux+fux][uk]*(1.-dup)*dux*(1.-duk);
  result += mkcm_table[iX-1][iA][up][ux][uk+fuk]*(1.-dup)*(1.-dux)*duk;
  result += mkcm_table[iX-1][iA][up+fup][ux+fux][uk]*dup*dux*(1.-duk);
  result += mkcm_table[iX-1][iA][up+fup][ux][uk+fuk]*dup*(1.-dux)*duk;
  result += mkcm_table[iX-1][iA][up][ux+fux][uk+fuk]*(1.-dup)*dux*duk;
  result += mkcm_table[iX-1][iA][up+fup][ux+fux][uk+fuk]*dup*dux*duk;

  return result;

}
*/

//////////had to numerically integrate
double m1kcm(double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk, gsl_integration_workspace *wkcm)
{
    if (x!=x)
    {
        std::cout << "x gone bad m1kcm" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;

    //if (use_kcm_tables) return do_kcm_interp(1, x, pin, kcm, int(g), int(d), int(n));

    gsl_function F;
    struct makcm_params params = {x,pin,g,d, particleA, n,1,wdk};

    F.function = &makcm_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m1kcm  "<<std::endl;
        exit(EXIT_FAILURE);
    }
    gsl_integration_qags(&F, 0, kcm, epsabs, epsrel, limit, wkcm, &result, &abserr);

    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;

    //std::cout << "starting m1kcm "<< std::endl;
    // exit(EXIT_FAILURE);
    // long double mysum=0;
    // double mystep=kcm/kcmres;
    // double mykcm=mystep;
    // double inc=0;
    // for (int l = 1; l < kcmres; l++)
    // {
    //     inc=(m1dk(inf,mykcm,x,pin,g,d,n)-m1dk(0,mykcm,x,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     mykcm+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m1kcm failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // //std::cout << "finishing m1kcm "<< mysum<< std::endl;
    // return mysum;
}
double m2kcm(double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk, gsl_integration_workspace *wkcm )
{
    // if (kcmres<kcm)
    // {
    //     std::cout << "m2res too short failure" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    //std::cout << "starting m2kcm "<< std::endl;
    // exit(EXIT_FAILURE);

    //if (use_kcm_tables) return do_kcm_interp(2, x, pin, kcm, int(g), int(d), int(n));
    
    gsl_function F;
    struct makcm_params params = {x,pin,g,d, particleA, n,2,wdk};

    F.function = &makcm_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m2kcm  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, kcm, epsabs, epsrel, limit, wkcm, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=kcm/kcmres;
    // double mykcm=mystep;
    // double inc=0;
    // for (int l = 1; l < kcmres; l++)
    // {
    //     inc=(m2dk(inf,mykcm,x,pin,g,d,n)-m2dk(0,mykcm,x,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     mykcm+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m2kcm failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // //std::cout << "finishing m2kcm "<< mysum<< std::endl;
    // return mysum;
}
double m3kcm(double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk, gsl_integration_workspace *wkcm )
{
    // if (kcmres<kcm)
    // {
    //     std::cout << "m3res too short failure" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    //std::cout << "starting m3kcm "<< std::endl;
    // exit(EXIT_FAILURE);
    
    //if (use_kcm_tables) return do_kcm_interp(3, x, pin, kcm, int(g), int(d), int(n));

    gsl_function F;
    struct makcm_params params = {x,pin,g,d, particleA, n,3,wdk};

    F.function = &makcm_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m3kcm  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, kcm, epsabs, epsrel, limit, wkcm, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=kcm/kcmres;
    // double mykcm=mystep;
    // double inc=0;
    // for (int l = 1; l < kcmres; l++)
    // {
    //     inc=(m3dk(inf,mykcm,x,pin,g,d,n)-m3dk(0,mykcm,x,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     mykcm+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m3kcm failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // //std::cout << "finishing m3kcm "<< mysum<< std::endl;
    // return mysum;
}
double m4kcm(double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk, gsl_integration_workspace *wkcm )
{
    // if (kcmres<kcm)
    // {
    //     std::cout << "m4res too short failure" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    
    //if (use_kcm_tables) return do_kcm_interp(4, x, pin, kcm, int(g), int(d), int(n));
    
    //std::cout << "starting m4kcm "<< std::endl;
    // exit(EXIT_FAILURE);
    gsl_function F;
    struct makcm_params params = {x,pin,g,d, particleA, n,4,wdk};

    F.function = &makcm_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m4kcm  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, kcm, epsabs, epsrel, limit, wkcm, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=kcm/kcmres;
    // double mykcm=mystep;
    // double inc=0;
    // for (int l = 1; l < kcmres; l++)
    // {
    //     inc=(m4dk(inf,mykcm,x,pin,g,d,n)-m4dk(0,mykcm,x,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     mykcm+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m4kcm failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // //std::cout << "finishing m4kcm "<< mysum<< std::endl;
    // return mysum;
}
double m5kcm(double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk, gsl_integration_workspace *wkcm)
{
    // if (kcmres<kcm)
    // {
    //     std::cout << "m5res too short failure" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    //std::cout << "starting m5kcm "<< std::endl;
    // exit(EXIT_FAILURE);
    
    //if (use_kcm_tables) return do_kcm_interp(5, x, pin, kcm, int(g), int(d), int(n));
    
    gsl_function F;
    struct makcm_params params = {x,pin,g,d, particleA, n,5,wdk};

    F.function = &makcm_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m5kcm  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, kcm, epsabs, epsrel, limit, wkcm, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=kcm/kcmres;
    // double mykcm=mystep;
    // double inc=0;
    // for (int l = 1; l < kcmres; l++)
    // {
    //     inc=(m5dk(inf,mykcm,x,pin,g,d,n)-m5dk(0,mykcm,x,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     mykcm+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m5kcm failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // //std::cout << "finishing m5kcm "<< mysum<< std::endl;
    // return mysum;
}
double m6kcm(double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk, gsl_integration_workspace *wkcm)
{
    // if (kcmres<kcm)
    // {
    //     std::cout << "m6res too short failure" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    //std::cout << "starting m6kcm "<< std::endl;
    // exit(EXIT_FAILURE);
    
    //if (use_kcm_tables) return do_kcm_interp(6, x, pin, kcm, int(g), int(d), int(n));
    
    gsl_function F;
    struct makcm_params params = {x,pin,g,d, particleA, n,6,wdk};

    F.function = &makcm_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m6kcm  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, kcm, epsabs, epsrel, limit, wkcm, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=kcm/kcmres;
    // double mykcm=mystep;
    // double inc=0;
    // for (int l = 1; l < kcmres; l++)
    // {
    //     inc=(m6dk(inf,mykcm,x,pin,g,d,n)-m6dk(0,mykcm,x,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     mykcm+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m6kcm failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // //std::cout << "finishing m6kcm " << mysum << std::endl;
    // return mysum;
}
double m7kcm(double kcm,double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk, gsl_integration_workspace *wkcm)
{
    // if (kcmres<kcm)
    // {
    //     std::cout << "m7res too short failure" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
    if (x>1) return 0;
    if (kcm==0) return 0;
    if (x==0) return 0;
    // exit(EXIT_FAILURE);
    
    //if (use_kcm_tables) return do_kcm_interp(7, x, pin, kcm, int(g), int(d), int(n));
    
    gsl_function F;
    struct makcm_params params = {x,pin,g,d, particleA, n,7,wdk};

    F.function = &makcm_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m7kcm  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, kcm, epsabs, epsrel, limit, wkcm, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=kcm/kcmres;
    // double mykcm=mystep;
    // double inc=0;
    // for (int l = 1; l < kcmres; l++)
    // {
    //     inc=(m7dk(inf,mykcm,x,pin,g,d,n)-m7dk(0,mykcm,x,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     mykcm+=mystep;
    //     if (inc == inf || inc != inc||mysum==inf||mysum!=mysum)
    //         {
    //             std::cout << "m7kcm failure" <<inc<<"  "<<mysum<< std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // //std::cout << "finishing m7kcm "<< mysum<< std::endl;
    // //  std::cout << "m7kcm" <<inc<<"  "<<mysum<< std::endl;
    // return mysum;
}

struct max_params {double pin;double g;double d; double particleA; double n; double a;gsl_integration_workspace *wdk;gsl_integration_workspace *wkcm;};

double max_f (double y, void * p)
{
    struct max_params * params = (struct max_params *)p;
    double pin = (params->pin);
    double g = (params->g);
    double d = (params->d);
    double particleA = (params->particleA);
    double n = (params->n);
    double a= (params->a);
    gsl_integration_workspace *wdk= (params->wdk);
    gsl_integration_workspace *wkcm= (params->wkcm);
    
    if (y==0) return 0;
    if (y>1) return 0;
    if (y!=y)
    {
        std::cout << "max_f gone bad  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    if (a==1) return m1kcm(pin,y,pin,g,d, particleA, n,wdk,wkcm)-m1kcm(0,y,pin,g,d, particleA, n,wdk,wkcm);
    if (a==2) return m2kcm(pin,y,pin,g,d, particleA, n,wdk,wkcm)-m2kcm(0,y,pin,g,d, particleA, n,wdk,wkcm);
    if (a==3) return m3kcm(pin,y,pin,g,d, particleA, n,wdk,wkcm)-m3kcm(0,y,pin,g,d, particleA, n,wdk,wkcm);
    if (a==4) return m4kcm(pin,y,pin,g,d, particleA, n,wdk,wkcm)-m4kcm(0,y,pin,g,d, particleA, n,wdk,wkcm);
    if (a==5) return m5kcm(pin,y,pin,g,d, particleA, n,wdk,wkcm)-m5kcm(0,y,pin,g,d, particleA, n,wdk,wkcm);
    if (a==6) return m6kcm(pin,y,pin,g,d, particleA, n,wdk,wkcm)-m6kcm(0,y,pin,g,d, particleA, n,wdk,wkcm);
    if (a==7) return m7kcm(pin,y,pin,g,d, particleA, n,wdk,wkcm)-m7kcm(0,y,pin,g,d, particleA, n,wdk,wkcm);
}

//Adapt this to read either quark or gluon tables based on particleA
double do_interp(int iX, double x, double pin, int g, int d, int n, int particleA)
{

  if (max_n[iX-1]<n) max_n[iX-1]=n;

  if (n>17) {
    std::cout << "iX= " << iX << " n= " << n << " is too large for current tabulation! EXITING" << std::endl;
    return 0.;
    exit(1);
  }

  int iA = (n-1)*4;
  if (g == -1 && d == -1) iA += 0;
  if (g == -1 && d == 1) iA += 1;
  if (g == 1 && d == -1) iA += 2;
  if (g == 1 && d == 1) iA += 3;

  if (x > 1) {
    std::cout << " x > 1 in mIx, x = " << x << std::endl;
    exit(1);
    //x=1.; //CHECK
  }

  int ux = int(x/step_x);
  double dux = (x - double(ux)*step_x)/step_x;
  int fux = 1;
  if (ux == SIZE_X-1) fux = 0, dux = 0.;

  int up = int((pin - pin_min)/pin_step);
  double dup = (pin - pin_min - double(up)*pin_step)/pin_step;
  int fup = 1;
  if (up == SIZE_PIN-1) fup = 0, dup = 0.;

  double result = 0.;
    
  // particleType: Q = +1 for fermion, G = -1 for boson
  if (particleA == 1) {
    result += q_mx_table[iX-1][iA][up][ux]*(1.-dup)*(1.-dux);
    result += q_mx_table[iX-1][iA][up+fup][ux]*dup*(1.-dux);
    result += q_mx_table[iX-1][iA][up][ux+fux]*(1.-dup)*dux;
    result += q_mx_table[iX-1][iA][up+fup][ux+fux]*dup*dux;
  }
  else if (particleA == -1) {
    result += g_mx_table[iX-1][iA][up][ux]*(1.-dup)*(1.-dux);
    result += g_mx_table[iX-1][iA][up+fup][ux]*dup*(1.-dux);
    result += g_mx_table[iX-1][iA][up][ux+fux]*(1.-dup)*dux;
    result += g_mx_table[iX-1][iA][up+fup][ux+fux]*dup*dux;
  }
  else {
    cout << "Unrecognized particleA= " << particleA << endl;
    exit(1);
  }

  return result;

}

/*
double do_interp(int iX, double x, double pin, int g, int d, int n)
{
  
  if (n>17) {
    std::cout << "iX= " << iX << " n= " << n << " is too large for current tabulation! EXITING" << std::endl;
    return 0.;
    exit(1);
  }

  int iA = (n-1)*4;
  if (g == -1 && d == -1) iA += 0;
  if (g == -1 && d == 1) iA += 1;
  if (g == 1 && d == -1) iA += 2;
  if (g == 1 && d == 1) iA += 3;

  if (x > 1) {
    std::cout << " x > 1 in mIx, x = " << x << std::endl;
    exit(1);
    //x=1.; //CHECK
  }

  int iz = 0;
  double za[SIZE_PIN*SIZE_X];
  for (int iY=0; iY<SIZE_X; iY++) {
    for (int iP=0; iP<SIZE_PIN; iP++) {
      za[iz] = mx_table[iX-1][iA][iP][iY];
      //std::cout << mx_table[iX-1][iA][iP][iY] << " ";
      iz++;
    }
    //std::cout << std::endl;
  }

  gsl_spline2d_init(spline, pin_vals, x_vals, za, SIZE_PIN, SIZE_X);
  return gsl_spline2d_eval(spline, pin, x, xacc, yacc);

}
*/

double m1x(double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx ) 
{
    if (x==0) return 0;

    //std::cout << " n= " << n << std::endl;

    if (use_tables && pin<=1500. && pin>0.) return do_interp(1, x, pin, int(g), int(d), int(n), int(particleA));

    gsl_function F;
    struct max_params params = {pin,g,d, particleA, n,1,wdk,wkcm};
    if (x!=x)
    {
        std::cout << "m1x x gone bad  " << std::endl;
        exit(EXIT_FAILURE);
    }

    F.function = &max_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m1x  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, x, epsabs, epsrel, limit, wx, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=x/x_res;
    // double myx=mystep;
    // double inc=0;
    // for (int l = 0; l < x_res; l++)
    // {
    //     inc=(m1kcm(pin,myx,pin,g,d,n)-m1kcm(0,myx,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     myx+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m1x failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // return mysum;
}
double m2x(double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx )
{
    if (x==0) return 0;
    
    if (use_tables && pin<=1500. && pin>0.) return do_interp(2, x, pin, int(g), int(d), int(n), int(particleA));
    
    gsl_function F;
    struct max_params params = {pin,g,d, particleA, n,2,wdk,wkcm};

    F.function = &max_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m2x  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, x, epsabs, epsrel, limit, wx, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=x/x_res;
    // double myx=mystep;
    // double inc=0;
    // for (int l = 0; l < x_res; l++)
    // {
    //     inc=(m2kcm(pin,myx,pin,g,d,n)-m2kcm(0,myx,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     myx+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m2x failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // return mysum;
}
double m3x(double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx )
{
    if (x==0) return 0;

    if (use_tables && pin<=1500. && pin>0.) return do_interp(3, x, pin, int(g), int(d), int(n), int(particleA));
  
    gsl_function F;
    struct max_params params = {pin,g,d, particleA, n,3,wdk,wkcm};

    F.function = &max_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m3x  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, x, epsabs, epsrel, limit, wx, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=x/x_res;
    // double myx=mystep;
    // double inc=0;
    // for (int l = 0; l < x_res; l++)
    // {
    //     inc=(m3kcm(pin,myx,pin,g,d,n)-m3kcm(0,myx,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     myx+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m3x failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // return mysum;
}
double m4x(double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx )
{
    if (x==0) return 0;

    if (use_tables && pin<=1500. && pin>0.) return do_interp(4, x, pin, int(g), int(d), int(n), int(particleA));
   
    gsl_function F;
    struct max_params params = {pin,g,d, particleA, n,4,wdk,wkcm};

    F.function = &max_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m4x  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, x, epsabs, epsrel, limit, wx, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=x/x_res;
    // double myx=mystep;
    // double inc=0;
    // for (int l = 0; l < x_res; l++)
    // {
    //     inc=(m4kcm(pin,myx,pin,g,d,n)-m4kcm(0,myx,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     myx+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m4x failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // return mysum;
}
double m5x(double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx )
{
    if (x==0) return 0;
    
    if (use_tables && pin<=1500. && pin>0.) return do_interp(5, x, pin, int(g), int(d), int(n), int(particleA));
    
    gsl_function F;
    struct max_params params = {pin,g,d, particleA, n,5,wdk,wkcm};

    F.function = &max_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m5x  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, x, epsabs, epsrel, limit, wx, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=x/x_res;
    // double myx=mystep;
    // double inc=0;
    // for (int l = 0; l < x_res; l++)
    // {
    //     inc=(m5kcm(pin,myx,pin,g,d,n)-m5kcm(0,myx,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     myx+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m1x failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // return mysum;
}

double m6x(double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx )
{
    if (x==0) return 0;
    
    if (use_tables && pin<=1500. && pin>0.) return do_interp(6, x, pin, int(g), int(d), int(n), int(particleA));
    
    gsl_function F;
    struct max_params params = {pin,g,d, particleA, n,6,wdk,wkcm};

    F.function = &max_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m6x  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, x, epsabs, epsrel, limit, wx, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=x/x_res;
    // double myx=mystep;
    // double inc=0;
    // for (int l = 0; l < x_res; l++)
    // {
    //     inc=(m6kcm(pin,myx,pin,g,d,n)-m6kcm(0,myx,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     myx+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m6x failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // return mysum;
}
double m7x(double x,double pin,double g,double d, double particleA, double n,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx )
{
    if (x==0) return 0;

    if (use_tables && pin<=1500. && pin>0.) return do_interp(7, x, pin, int(g), int(d), int(n), int(particleA));

    gsl_function F;
    struct max_params params = {pin,g,d, particleA, n,7,wdk,wkcm};

    F.function = &max_f;
    F.params = &params;
    double result;
    double abserr=0;
    double epsabs=abs_res_kcmx;
    double epsrel=rel_res_kcmx;
    double limit=100000;
    gsl_set_error_handler_off();
    if(!wdk)
    {
        std::cout << "set the workspace in m7x  "<<std::endl;
        exit(EXIT_FAILURE);
       
    } 
    gsl_integration_qags(&F, 0, x, epsabs, epsrel, limit, wx, &result, &abserr);
    // gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    // std::cout<<result<<std::endl;
    return result;
    // long double mysum=0;
    // double mystep=x/x_res;
    // double myx=mystep;
    // double inc=0;
    // for (int l = 0; l < x_res; l++)
    // {
    //     inc=(m7kcm(pin,myx,pin,g,d,n)-m7kcm(0,myx,pin,g,d,n))*mystep;
    //     mysum+=inc;
    //     myx+=mystep;
    //     if (inc == inf || inc != inc)
    //         {
    //             std::cout << "m7x failure" << std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    // }
    // return mysum;
}


















double m1tot(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA, double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double n=0) //dist is # of un-intgrated variables (0 for PhiKThetaP,1 for PhiKTheta,2 for PhiK,3 for Phi,4 for all 4 left)
{
    //return 0; //DEBUG

    if(!wdk)
    {
        std::cout << "set the workspace in m1tot  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    }
    if (x!=x)
    {
        std::cout << "m1tot x gone bad  " << std::endl;
        exit(EXIT_FAILURE);
    }
    if (dist==0) return m1x(x,pin,g,d, particleA, n,wdk,wkcm,wx)-m1x(0,pin,g,d, particleA, n,wdk,wkcm,wx);
    if (dist==1) return m1kcm(kcm,x,pin,g,d, particleA, n,wdk,wkcm)-m1kcm(0,x,pin,g,d, particleA, n,wdk,wkcm);
    if (dist==2 && dk!=0) return m1dk(dk,kcm,x,pin,g,d, particleA, n,wdk)-m1dk(0,kcm,x,pin,g,d, particleA, n,wdk);
    if (dist==2 && dk==0) return 0.; 
    if (dist==3) return m1phi(phi,dk,kcm,x,pin,g,d, particleA)-m1phi(0,dk,kcm,x,pin,g,d, particleA);
    if (dist==4) return m1base(phi,dk,kcm,x,pin,g,d, particleA);
    std::cout<<"m1tot issue"<<std::endl;
    exit(EXIT_FAILURE);
}
double m2tot(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA, double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double n=0)
{
    //return 0; //DEBUG
    if(!wdk)
    {
        std::cout << "set the workspace in m2tot  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    if (dist==0) return m2x(x,pin,g,d, particleA, n,wdk,wkcm,wx)-m2x(0,pin,g,d, particleA, n,wdk,wkcm,wx);
    if (dist==1) return m2kcm(kcm,x,pin,g,d, particleA, n,wdk,wkcm)-m2kcm(0,x,pin,g,d, particleA, n,wdk,wkcm);
    if (dist==2 && dk!=0.) return m2dk(dk,kcm,x,pin,g,d, particleA, n,wdk)-m2dk(0,kcm,x,pin,g,d, particleA, n,wdk);
    if (dist==2 && dk==0.) return 0.; 
    if (dist==3) return m2phi(phi,dk,kcm,x,pin,g,d, particleA)-m2phi(0,dk,kcm,x,pin,g,d, particleA);
    if (dist==4) return m2base(phi,dk,kcm,x,pin,g,d, particleA);
    std::cout<<"m2tot issue"<<std::endl;
    exit(EXIT_FAILURE);
}
double m3tot(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA, double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double n=0)
{
    //return 0; //DEBUG
    if(!wdk)
    {
        std::cout << "set the workspace in m3tot "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    if (dist==0) return m3x(x,pin,g,d, particleA, n,wdk,wkcm,wx)-m3x(0,pin,g,d, particleA, n,wdk,wkcm,wx);
    if (dist==1) return m3kcm(kcm,x,pin,g,d, particleA, n,wdk,wkcm)-m3kcm(0,x,pin,g,d, particleA, n,wdk,wkcm);
    if (dist==2 && dk!=0.) return m3dk(dk,kcm,x,pin,g,d, particleA, n,wdk)-m3dk(0,kcm,x,pin,g,d, particleA, n,wdk);
    if (dist==2 && dk==0.) return 0.;
    if (dist==3) return m3phi(phi,dk,kcm,x,pin,g,d, particleA)-m3phi(0,dk,kcm,x,pin,g,d, particleA);
    if (dist==4) return m3base(phi,dk,kcm,x,pin,g,d, particleA);
    std::cout<<"m3tot issue"<<std::endl;
    exit(EXIT_FAILURE);
}
double m4tot(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA, double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double n=0)
{
    //return 0; //DEBUG
    if(!wdk)
    {
        std::cout << "set the workspace in m4tot  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    if (dist==0) return m4x(x,pin,g,d, particleA, n,wdk,wkcm,wx)-m4x(0,pin,g,d, particleA, n,wdk,wkcm,wx);
    if (dist==1) return m4kcm(kcm,x,pin,g,d, particleA, n,wdk,wkcm)-m4kcm(0,x,pin,g,d, particleA, n,wdk,wkcm);
    if (dist==2 && dk!=0.) return m4dk(dk,kcm,x,pin,g,d, particleA, n,wdk)-m4dk(0,kcm,x,pin,g,d, particleA, n,wdk);
    if (dist==2 && dk==0.) return 0.;
    if (dist==3) return m4phi(phi,dk,kcm,x,pin,g,d, particleA)-m4phi(0,dk,kcm,x,pin,g,d, particleA);
    if (dist==4) return m4base(phi,dk,kcm,x,pin,g,d, particleA);
    std::cout<<"m4tot issue"<<std::endl;
    exit (EXIT_FAILURE);
}
double m5tot(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA, double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double n=0)
{
    //return 0; //DEBUG
    if(!wdk)
    {
        std::cout << "set the workspace in m5tot  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    if (dist==0) return m5x(x,pin,g,d, particleA, n,wdk,wkcm,wx)-m5x(0,pin,g,d, particleA, n,wdk,wkcm,wx);
    if (dist==1) return m5kcm(kcm,x,pin,g,d, particleA, n,wdk,wkcm)-m5kcm(0,x,pin,g,d, particleA, n,wdk,wkcm);
    if (dist==2 && dk!=0.) return m5dk(dk,kcm,x,pin,g,d, particleA, n,wdk)-m5dk(0,kcm,x,pin,g,d, particleA, n,wdk);
    if (dist==2 && dk==0.) return 0.;
    if (dist==3) return m5phi(phi,dk,kcm,x,pin,g,d, particleA)-m5phi(0,dk,kcm,x,pin,g,d, particleA);
    if (dist==4) return m5base(phi,dk,kcm,x,pin,g,d, particleA);
    std::cout<<"m5tot issue"<<std::endl;
    exit(EXIT_FAILURE);
}
double m6tot(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA, double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double n=0)
{
    //return 0; //DEBUG
    if(!wdk)
    {
        std::cout << "set the workspace in m6tot  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    if (dist==0) return m6x(x,pin,g,d, particleA, n,wdk,wkcm,wx)-m6x(0,pin,g,d, particleA, n,wdk,wkcm,wx);
    if (dist==1) return m6kcm(kcm,x,pin,g,d, particleA, n,wdk,wkcm)-m6kcm(0,x,pin,g,d, particleA, n,wdk,wkcm);
    if (dist==2 && dk!=0.) return m6dk(dk,kcm,x,pin,g,d, particleA, n,wdk)-m6dk(0,kcm,x,pin,g,d, particleA, n,wdk);
    if (dist==2 && dk==0.) return 0.;
    if (dist==3) return m6phi(phi,dk,kcm,x,pin,g,d, particleA)-m6phi(0,dk,kcm,x,pin,g,d, particleA);
    if (dist==4) return m6base(phi,dk,kcm,x,pin,g,d, particleA);
    std::cout<<"m6tot issue"<<std::endl;
    exit(EXIT_FAILURE);
}
double m7tot(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA, double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double n=0)
{
    //return 0; //DEBUG
    if(!wdk)
    {
        std::cout << "set the workspace in m7tot  "<<std::endl;
        exit(EXIT_FAILURE);
       
        
    } 
    if (dist==0) return m7x(x,pin,g,d, particleA, n,wdk,wkcm,wx)-m7x(0,pin,g,d, particleA, n,wdk,wkcm,wx);
    if (dist==1) return m7kcm(kcm,x,pin,g,d, particleA, n,wdk,wkcm)-m7kcm(0,x,pin,g,d, particleA, n,wdk,wkcm);
    if (dist==2 && dk!=0.) return m7dk(dk,kcm,x,pin,g,d, particleA, n,wdk)-m7dk(0,kcm,x,pin,g,d, particleA, n,wdk);
    if (dist==2 && dk==0.) return 0.;
    if (dist==3) return m7phi(phi,dk,kcm,x,pin,g,d, particleA)-m7phi(0,dk,kcm,x,pin,g,d, particleA);
    if (dist==4) return m7base(phi,dk,kcm,x,pin,g,d, particleA);
    std::cout<<"m7tot issue"<<std::endl;
    exit(EXIT_FAILURE);
}
std::vector<double> meval(double phi,double dk,double kcm,double x,double pin,double g,double d, double particleA, double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx)
{//This combines some of the sums in mtp in order to speed up runtime and cleans up the rest
    double eta=pow(4*M_PI,-3)/pow(pin,2);
    if(!wdk)
    {
        std::cout << "set the workspace in meval  "<<std::endl;
        exit(EXIT_FAILURE);
        
    } 
    if (dist <=2)
    {
        double sum1=0;
        double sum2=0;
        double sum3=0;
        double sum4=0;
        double sum5=0;
        double sum6=0;
        double sum7=0;

        //return {eta*sum1,eta*sum2,eta*sum3,eta*sum4,eta*sum5,eta*sum6,eta*sum7}; //DEBUG

        double inc=1000;
        double pinc=1000;
        for (int n=1; abs(pinc) >=ep/20 * abs(sum1)||(abs(inc)>abs(pinc)); n++)
        { 
	    //if (n>15) break; //DEBUG
            pinc=inc;
            inc=m1tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,n);
            sum1+=inc;
	    //std::cout << " n1 = " << n << " inc= " << inc << std::endl;
            if(n>=5&&((inc==0&&sum1==0)||(inc==0&&pinc==0))) break;
        }
	//std::cout << std::setprecision(9) << "Sums= " << sum1 << " " << sum2 << " " << sum3 << " " << sum4 << " " << sum5 << " " << sum6 << " " << sum7 << std::endl;
        //return {eta*sum1,eta*sum2,eta*sum3,eta*sum4,eta*sum5,eta*sum6,eta*sum7};

	inc=1000;
        for (int n=1; abs(pinc) >=ep/20 * abs(sum2)||(abs(inc)>abs(pinc)); n++)
        { 
	    //if (n>15) break; //DEBUG
            pinc=inc;
            inc=m2tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,n);
            sum2+=inc;
	    //std::cout << " n2 = " << n << std::endl;
            if(n>=5&&((inc==0&&sum2==0)||(inc==0&&pinc==0))) break;
        }
        inc=1000;
        for (int n=1; abs(pinc) >=ep/20 * abs(sum3)||(abs(inc)>abs(pinc)); n++)
        { 
	    //if (n>15) break; //DEBUG
            pinc=inc;
            inc=m3tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,n);
            sum3+=inc;
	    //std::cout << " n3 = " << n << std::endl;
            if(n>=5&&((inc==0&&sum3==0)||(inc==0&&pinc==0))) break;
        }
        inc=1000;
        for (int n=1; abs(pinc) >=ep/20 * abs(sum4)||(abs(inc)>abs(pinc)); n++)
        { 
            // std::cout<<n<<"   "<<inc<<std::endl;
	    //if (n>15) break; //DEBUG
            pinc=inc;
            inc=m4tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,n);
            sum4+=inc;
	    //std::cout << " n4 = " << n << std::endl;
            if(n>=5&&((inc==0&&sum4==0)||(inc==0&&pinc==0))) break;
        }
        inc=1000;
        for (int n=1; (abs(pinc) >=ep/20 * abs(sum5))||(abs(inc)>abs(pinc)); n++)
        { 
	    //if (n>15) break; //DEBUG
            pinc=inc;
            inc=m5tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,n);
            sum5+=inc;
	    //std::cout << " n5 = " << n << std::endl;
            if(n>=5&&((inc==0&&sum5==0)||(inc==0&&pinc==0))) break;
        }
        inc=1000;
        for (int n=1; (abs(pinc) >=ep/20 * abs(sum6))||(abs(inc)>abs(pinc)); n++)
        { 
	    //if (n>15) break; //DEBUG
            pinc=inc;
            inc=m6tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,n);
            sum6+=inc;
	    //std::cout << " n6 = " << n << std::endl;
            if(n>=5&&((inc==0&&sum6==0)||(inc==0&&pinc==0))) break;
        }
        inc=1000;
        //int myN=0;
        for (int n=1; (abs(pinc) >=ep/100 * abs(sum7))||(abs(inc)>abs(pinc)); n++)
        { 
            //myN+=1;
	    //if (n>15) break; //DEBUG
            pinc=inc;
            inc=m7tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,n);
	    //std::cout << " n7 = " << n << std::endl;
            sum7+=inc;
            if(n>=5&&((inc==0&&sum7==0)||(inc==0&&pinc==0))) break;
        }

        //cout<<myN<<std::endl;
        return {eta*sum1,eta*sum2,eta*sum3,eta*sum4,eta*sum5,eta*sum6,eta*sum7};
        
    }
    return {eta*m1tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,0),eta*m2tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,0),eta*m3tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,0),eta*m4tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,0),eta*m5tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,0),eta*m6tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,0),eta*m7tot(phi,dk,kcm,x,pin,g,d, particleA, dist,wdk,wkcm,wx,0)};
    
}

std::vector<double> Proc1_qq_qq()
{
    double C1=16*pow(dF*CF,2)/dA;
    double C2=16*dF*CF*(CF-CA/2);
    //return C1*(m1tot(phi,dk,kcm,x,g,d,dist)+m2tot(phi,dk,kcm,x,g,d,dist)+m3tot(phi,dk,kcm,x,g,d,dist)-m6tot(phi,dk,kcm,x,g,d,dist)+m7tot(phi,dk,kcm,x,g,d,dist))+C2*(-m2tot(phi,dk,kcm,x,g,d,dist)+m3tot(phi,dk,kcm,x,g,d,dist)-m6tot(phi,dk,kcm,x,g,d,dist));
    return {C1,C2-C1,C1+C2,0,0,C1+C2,C1};
}
std::vector<double>Proc3_qqb_qqb()
{
    double C1=16*pow(dF*CF,2)/dA;
    double C2=16*dF*CF*(CF-CA/2);
    //return C1*(m1tot(phi,dk,kcm,x,g,d,dist)+m2tot(phi,dk,kcm,x,g,d,dist)+m3tot(phi,dk,kcm,x,g,d,dist)+m5tot(phi,dk,kcm,x,g,d,dist)+m4tot(phi,dk,kcm,x,g,d,dist))+C2*(m2tot(phi,dk,kcm,x,g,d,dist)+m4tot(phi,dk,kcm,x,g,d,dist)+2*m3tot(phi,dk,kcm,x,g,d,dist));
    return {C1,-C1-C2,C1+2*C2,-C1-C2,C1,0,0};
    
}
std::vector<double> Proc3twid_qqb_qqb()
{
    double C1=16*pow(dF*CF,2)/dA;
    double C2=16*dF*CF*(CF-CA/2);
    //return C1*(m3tot(phi,dk,kcm,x,g,d,dist)-m6tot(phi,dk,kcm,x,g,d,dist)+m7tot(phi,dk,kcm,x,g,d,dist)+m5tot(phi,dk,kcm,x,g,d,dist)+m4tot(phi,dk,kcm,x,g,d,dist))+C2*(-m4tot(phi,dk,kcm,x,g,d,dist)+m6tot(phi,dk,kcm,x,g,d,dist));
    return {0,0,C1,C2-C1,C1,C1-C2,C1};
    
}
std::vector<double> Proc4_qqp_qqp()
{
    double C=8*pow(dF*CF,2)/dA;
    //return 8*pow(dF*CF,2)/dA*(2*m1tot(phi,dk,kcm,x,g,d,dist)+2*m2tot(phi,dk,kcm,x,g,d,dist)+m3tot(phi,dk,kcm,x,g,d,dist));
    return {2*C,-2*C,C,0,0,0,0};
    
}
std::vector<double> Proc4twid_qqp_qqp()
{
    double C=8*pow(dF*CF,2)/dA;
    //return 8*pow(dF*CF,2)/dA*(m3tot(phi,dk,kcm,x,g,d,dist)-2*m6tot(phi,dk,kcm,x,g,d,dist)+2*m7tot(phi,dk,kcm,x,g,d,dist));
    return {0,0,C,0,0,2*C,2*C};
    
}
std::vector<double> Proc7_qqb_qpqbp()
{
    double C=8*pow(dF*CF,2)/dA;
    //return 8*pow(dF*CF,2)/dA*(2*m5tot(phi,dk,kcm,x,g,d,dist)+2*m4tot(phi,dk,kcm,x,g,d,dist)+m3tot(phi,dk,kcm,x,g,d,dist));
    return {0,0,C,-2*C,2*C,0,0};
    
}
std::vector<double> Proc8_qqb_gg()
{
    double C1=-8*dF*CF*CF;
    double C2=-8*dF*CF*CA;
    //return 8*dF*CF*(CF*(-m6tot(phi,dk,kcm,x,g,d,dist)-m2tot(phi,dk,kcm,x,g,d,dist)-m3tot(phi,dk,kcm,x,g,d,dist))-CA*(2*m5tot(phi,dk,kcm,x,g,d,dist)+2*m4tot(phi,dk,kcm,x,g,d,dist)+m3tot(phi,dk,kcm,x,g,d,dist)));
    return {0,-C1,C1+C2,-2*C2,2*C2,-C1,0};
    
}
std::vector<double> Proc9_qg_qg()
{
    double C1=8*dF*CF*CF;
    double C2=8*dF*CF*CA;
    //return 8*dF*CF*(CF*(2*m3tot(phi,dk,kcm,x,g,d,dist)+m4tot(phi,dk,kcm,x,g,d,dist)-m6tot(phi,dk,kcm,x,g,d,dist))+CA*(2*m1tot(phi,dk,kcm,x,g,d,dist)+2*m2tot(phi,dk,kcm,x,g,d,dist)+m3tot(phi,dk,kcm,x,g,d,dist)));
    return {2*C2,-2*C2,2*C1+C2,-C1,0,C1,0};
    
}
std::vector<double> Proc9twid_qg_qg()
{
    double C1=8*dF*CF*CF;
    double C2=8*dF*CF*CA;
    //return 8*dF*CF*(CF*(-m2tot(phi,dk,kcm,x,g,d,dist)-m4tot(phi,dk,kcm,x,g,d,dist))+CA*(m3tot(phi,dk,kcm,x,g,d,dist)-2*m6tot(phi,dk,kcm,x,g,d,dist)+2*m7tot(phi,dk,kcm,x,g,d,dist)));
    return {0,C1,C2,C1,0,2*C2,2*C2};
    
}
std::vector<double> Proc11_gg_gg()
{
    double C=16*dA*CA*CA;
    //return 16*dA*CA*CA*(3*m3tot(phi,dk,kcm,x,g,d,dist)+m1tot(phi,dk,kcm,x,g,d,dist)+m2tot(phi,dk,kcm,x,g,d,dist)+m7tot(phi,dk,kcm,x,g,d,dist)-m6tot(phi,dk,kcm,x,g,d,dist)+m5tot(phi,dk,kcm,x,g,d,dist)+m4tot(phi,dk,kcm,x,g,d,dist));
    return {C,-C,3*C,-C,C,C,C};
    
}

double dot(std::vector<double> a,std::vector<double> b)
{
    double total=0;
    for (int i=0;i<7;++i) total+=a[i]*b[i];
    return total;
}

//the following will return multi index lists with probs for each of the subprocesses
std::vector<double> F_QQ(double phi, double dk, double kcm, double x, double pin, double deltime, double dist, gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double charm = 0)
{//<the first return index is the probability,-1 for gluon,1 for quark,2 for QBar for each of the two indices.
    std::vector<double> mQQ=meval(phi,dk,kcm,x,pin,Q,Q, Q, dist,wdk,wkcm,wx);
    std::vector<double> mGG = meval(phi, dk, kcm, x, pin, G, G, Q, dist,wdk,wkcm,wx);
    double c1=1;
    double c2=1;
    if (abs(charm)==4)
    {
        c1=0;
        c2=Nf/(Nf-1);
    }
    double kap=pow(gs,4)*deltime;
    double p1=dot(mQQ,Proc1_qq_qq())*c1;
    double p3=dot(mQQ,Proc3_qqb_qqb())*c1;
    double p4=(Nf-1)*2*dot(mQQ,Proc4_qqp_qqp())*c2;
    double p4t=(Nf-1)*dot(mQQ,Proc4twid_qqp_qqp())*c2;
    double p7=(Nf-1)*dot(mQQ,Proc7_qqb_qpqbp())*c1;
    double p9=dot(mGG,Proc9_qg_qg());
    //cout<<p1<<p3<<p7<<std::endl;

    double prob=kap/vq*(p1+p3+p4+p4t+p7+p9);
    if (prob==0) return {0,0};
    double guess=distcon(generator);
    //cout<<p1<<" "<<p3<<" "<<prob<<" "<<guess<<" hi"<<std::endl;
    if (guess<=p1*kap/vq/prob) return {prob,1};
    if (guess<=(p1+p3)*kap/vq/prob) return {prob,3};
    if (guess<=(p1+p3+p4)*kap/vq/prob) return {prob,4};
    if (guess<=(p1+p3+p4+p4t)*kap/vq/prob) return {prob,-4}; 
    if (guess<=(p1+p3+p4+p4t+p7)*kap/vq/prob) return {prob,7};
    if (guess<=1) return {prob,9};
    std::cout<<"F_QQ problem"<<std::endl;
    exit(EXIT_FAILURE);
}
std::vector<double> F_QG(double phi,double dk,double kcm,double x,double pin,double deltime,double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double charm=0)
{
    std::vector<double> mQG=meval(phi,dk,kcm,x,pin,Q,G, G, dist,wdk,wkcm,wx);
    std::vector<double> mGQ=meval(phi,dk,kcm,x,pin,G,Q, G, dist,wdk,wkcm,wx);
    double c1=1;
    if (abs(charm)==4) c1=0;
    double kap=pow(gs,4)*deltime;
    double p8=2*dot(mQG,Proc8_qqb_gg())*c1;
    double p9t=dot(mGQ,Proc9twid_qg_qg());
    double guess=distcon(generator);
    double prob=kap/vq*(p8+p9t);
    if (prob==0) return {0,0};
    if (guess<=p8/(p8+p9t)) return {prob,8}; 
    if (guess<=1) return {prob,-9};
    std::cout<<"F_QG problem"<<std::endl;
    exit (EXIT_FAILURE);
}
std::vector<double> F_GQ(double phi,double dk,double kcm,double x,double pin,double deltime,double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx)
{
    std::vector<double> mQG=meval(phi,dk,kcm,x,pin,Q,G, Q, dist,wdk,wkcm,wx);
    std::vector<double> mGQ=meval(phi,dk,kcm,x,pin,G,Q, Q, dist,wdk,wkcm,wx);
    double kap=pow(gs,4) * deltime;
    double p8=Nf*dot(mGQ,Proc8_qqb_gg());
    double p9t=Nf*dot(mQG,Proc9twid_qg_qg());
    double prob=kap/vg*(p8+p9t);
    if (prob==0) return {0,0};
    double guess=distcon(generator);
    if (guess<=p8/(p8+p9t)) return {prob,8};
    if (guess<=1) return {prob,-9};
    std::cout<<"F_GQ problem"<<std::endl;
    exit (EXIT_FAILURE);
}

std::vector<double> F_GG(double phi,double dk,double kcm,double x,double pin,double deltime,double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx)
{
    std::vector<double> mQQ=meval(phi,dk,kcm,x,pin,Q,Q, G, dist,wdk,wkcm,wx);
    std::vector<double> mGG=meval(phi,dk,kcm,x,pin,G,G, G, dist,wdk,wkcm,wx);
    double kap=pow(gs,4)*deltime;
    //prdouble(mQQ);
    double p9=2*Nf*dot(mQQ,Proc9_qg_qg());
    double p11=dot(mGG,Proc11_gg_gg());
    double prob=kap/vg*(p9+p11);
    if (prob==0) return {0,0};
    double guess=distcon(generator);
    if (guess<=p9/(p9+p11)) return {prob,9};
    if (guess<=1) return {prob,11};
    std::cout<<"F_GG problem"<<std::endl;
    exit (EXIT_FAILURE);
}
std::vector<double> F_QQB(double phi,double dk,double kcm,double x,double pin,double deltime,double dist,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx,double charm=0)
{
    std::vector<double> mQQ=meval(phi,dk,kcm,x,pin,Q,Q, Q, dist,wdk,wkcm,wx);
    double c1=1;
    double c2=1;
    if (abs(charm)==4)
    {
        c1=0;
        c2=Nf/(Nf-1);
    }
    double kap=pow(gs,4)*deltime;
    double p3t=dot(mQQ,Proc3twid_qqb_qqb())*c1;
    double p4t=(Nf-1)*dot(mQQ,Proc4twid_qqp_qqp())*c2;
    double p7=(Nf-1)*dot(mQQ,Proc7_qqb_qpqbp())*c1;
    double prob=kap/vq*(p3t+p4t+p7);
    if (prob==0) return {0,0};
    double guess=distcon(generator);
    if (guess<=p3t/(p3t+p4t+p7)) return {prob,-3};
    if (guess<=(p3t+p4t)/(p3t+p4t+p7)) return {prob,-4};
    if (guess<=1) return {prob,7};
    std::cout<<"F_QQB problem"<<std::endl;
    exit (EXIT_FAILURE);
}

double totalprob(double phi , double dk, double kcm, double x,double pin,double deltime, double dist,double pinPDG,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx)
{
    double pintype=-1000;

    if (pinPDG>=-4&&pinPDG<=4) pintype=Q;
    if (pinPDG==21) pintype=G;

    //cout << "check" << endl;
    if (pintype==Q)
    {
        //cout<<"its a Quark"<<endl;
        //cout << "check1" << endl;
        double FQQ=F_QQ(phi,dk,kcm,x,pin,deltime,dist,wdk,wkcm,wx,pinPDG)[0];
        double FQG=F_QG(phi,dk,kcm,x,pin,deltime,dist,wdk,wkcm,wx,pinPDG)[0];
        double FQQB=F_QQB(phi,dk,kcm,x,pin,deltime,dist,wdk,wkcm,wx,pinPDG)[0];

        //cout<<FQQ<<" "<<FQG<<" "<<FQQB<<endl;
        double prob=FQQ+FQG+FQQB;
        return prob;
        }
    if (pintype==G)
    {
        //cout << "its Gluon" << endl;
        double FGQ = F_GQ(phi, dk, kcm, x, pin, deltime, dist,wdk,wkcm,wx)[0];
        double FGG = F_GG(phi, dk, kcm, x, pin, deltime, dist,wdk,wkcm,wx)[0];
        double prob=FGQ+FGG;
        // return ExpCorr(prob);
        return prob;
    }
    std::cout<<"totalprob is having issues"<<std::endl;
    exit(EXIT_FAILURE);

    
}

std::vector<double> is_split(double pin, double deltime, double pinPDG, gsl_integration_workspace *wdk, gsl_integration_workspace *wkcm, gsl_integration_workspace *wx)
{
    // decides if a splitting happens and to what (input temp-normalized values)
    // returns (a,b,c,d) a is 1 if splitting 0 if no splitting,and -1 if try again with changed time,b is the mult change in the time,c is 1 if FQQ,2 if FQG,3 if FQQB,4 if FGQ,and 5 if FGG and 0 if nothing,d is the number of the process
    // the last entry is the total probability of splitting on this interval
    // cout<<"starting is_split"<<std::endl;

    double b=1;
    double pintype=-1000;

    if (pinPDG>=-4&&pinPDG<=4) pintype=Q;
    if (pinPDG==21) pintype=G;

    //cout<<"check"<<std::endl;
    if (pintype==Q)
    {
        //cout<<"its a Quark"<<std::endl;
        std::vector<double> FQQ_=F_QQ(0,0,0,1,pin,deltime,0,wdk,wkcm,wx,pinPDG);
        std::vector<double> FQG_=F_QG(0,0,0,1,pin,deltime,0,wdk,wkcm,wx,pinPDG);
        std::vector<double> FQQB_=F_QQB(0,0,0,1,pin,deltime,0,wdk,wkcm,wx,pinPDG);
        //cout<<"check1"<<std::endl;
        double FQQ=FQQ_[0];
        double FQG=FQG_[0];
        double FQQB=FQQB_[0];

        //cout<<FQQ<<" "<<FQG<<" "<<FQQB<<std::endl;
        double prob=FQQ+FQG+FQQB;
        double corr=ExpCorr(prob)/prob;
        if (prob==0) corr=1;
        //Why prob corrected twice??
	prob*=corr;
        FQQ*=corr;
        FQG*=corr;
        FQQB*=corr;
        //prob*=corr;
        if (prob>=.5) return {-1,.5,0,0,prob, pintype};//-1 means to half the global time step and run again

        if (prob<=.05) b=2; //means double the global time step on the next run
        double guess=distcon(generator);
         //cout<<guess<<"  "<<prob<<std::endl;
        
	//Is it invariant under interchange of FQQ FGQ and so on?
        if (guess>=prob) return {0,b,0,0,prob, pintype}; //0 means no splitting this interval
        if (guess<=FQQ) return {1,b,1,FQQ_[1],prob, pintype};
        if (guess<=(FQQ+FQG)) return {1,b,2,FQG_[1],prob, pintype};
        if (guess<=prob) return {1,b,3,FQQB_[1],prob, pintype};
        
    }
    if (pintype==G)
    {
        //std::cout<<"its Gluon"<<std::endl;
	//clock_t all_startClock = clock();
        std::vector<double> FGQ_=F_GQ(0,0,0,1,pin,deltime,0,wdk,wkcm,wx);
	//clock_t all_endClock = clock();
        //double all_time = double((all_endClock - all_startClock)) / CLOCKS_PER_SEC;
	//std::cout << " All_time = " << all_time << std::endl;
        std::vector<double> FGG_=F_GG(0,0,0,1,pin,deltime,0,wdk,wkcm,wx);
        double FGQ=FGQ_[0];
        double FGG=FGG_[0];
        double prob=FGQ+FGG;
        double corr=ExpCorr(prob)/prob;
        if (prob==0) corr=1;
        prob*=corr;
        FGQ*=corr;
        FGG*=corr;
        // prob*=corr;
        if (prob>=.5) return {-1,.5,0,0,prob, pintype};
        if (prob<=.05) b=2;
        double guess=distcon(generator);
        //cout<<guess<<"    "<<prob<<std::endl;

        if (guess>=prob) return {0,b,0,0,prob, pintype}; //0 means no splitting this interval
        if (guess<=FGQ) return {1,b,4,FGQ_[1],prob, pintype};
        if (guess<=prob) return {1,b,5,FGG_[1],prob, pintype};
    }
    std::cout<<"is_split is having issues"<<std::endl;
    exit(EXIT_FAILURE);
}
double myproc(double phi,double dk,double kcm,double x,double pin,double g,double d,double dist,double proc,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx, double pinType)
{
    std::vector<double> mQ = meval(phi,dk,kcm,x,pin,g,d, Q, dist,wdk,wkcm,wx);
    std::vector<double> mG = meval(phi,dk,kcm,x,pin,g,d, G, dist,wdk,wkcm,wx);

    if (proc==1) return dot(mQ,Proc1_qq_qq());
    if (proc==3) return dot(mQ,Proc3_qqb_qqb());
    if (proc==-3) return dot(mQ,Proc3twid_qqb_qqb());
    if (proc==4) return dot(mQ,Proc4_qqp_qqp());
    if (proc==-4) return dot(mQ,Proc4twid_qqp_qqp());
    if (proc==7) return dot(mQ,Proc7_qqb_qpqbp());
    if (proc==8) {
        if (pinType==Q) return dot(mG,Proc8_qqb_gg());
        else return dot(mQ,Proc8_qqb_gg());
    }
    if (proc==9) {
        if (g==Q) return dot(mG,Proc9_qg_qg());
        else return dot(mQ,Proc9_qg_qg());
    }
    if (proc==-9) {
        if (g==Q) return dot(mQ,Proc9twid_qg_qg());
        else return dot(mG,Proc9twid_qg_qg());
    }
    if (proc==11) return dot(mG,Proc11_gg_gg());
    std::cout<<"myproc failed"<<std::endl;
    exit (EXIT_FAILURE);
}

std::vector<double> q_nums(double proc, double pinPDG, double pQGQB) //-1 for gluon,1 for quark,2 for QBar for outgoing particle,input is PDG for incoming and process->pinPDG2,pPDG,kPDG,E4PDG,QG1,QG2
{//convention of the output is (pin,p,k,E4)
//doing this for only light quarks 
    double charge=distint(generator)*2-1;
    double flavor=distintNf(generator);

    //cout<<charge<<"   "<<flavor<<std::endl;
    if (proc==1) return {pinPDG,pinPDG,pinPDG,pinPDG,Q,Q}; //(pin,p,k,E4,Q/G from paper,Q/G from paper)
    if (proc==3) return {pinPDG,pinPDG,-pinPDG,-pinPDG,Q,Q}; 
    if (proc==-3) return {pinPDG,-pinPDG,-pinPDG,pinPDG,Q,Q}; 
    if (proc==4) 
    {
        while (flavor==abs(pinPDG)) flavor=distintNf(generator);
        return {pinPDG,pinPDG,charge*flavor,charge*flavor,Q,Q}; 
    }
    if (proc==-4)
    {
        while (flavor==abs(pinPDG)) flavor=distintNf(generator);
        if (pQGQB==1) return {pinPDG,sgn(pinPDG)*flavor,sgn(pinPDG)*flavor,pinPDG,Q,Q}; 
        if (pQGQB==2) return {pinPDG,-sgn(pinPDG)*flavor,-sgn(pinPDG)*flavor,pinPDG,Q,Q}; 
    }
    if (proc==7)
    {
        while (flavor==abs(pinPDG)) flavor=distintNf(generator);
        if (pQGQB==1) return {pinPDG,sgn(pinPDG)*flavor,-pinPDG,-sgn(pinPDG)*flavor,Q,Q}; 
        if (pQGQB==2) return {pinPDG,-sgn(pinPDG)*flavor,-pinPDG,sgn(pinPDG)*flavor,Q,Q}; 
    }
    if (proc==8)
    {
        if (pinPDG==21) return {21,charge*flavor,21,-charge*flavor,G,Q}; 
        return {pinPDG,21,-pinPDG,21,G,Q}; 
    }
    if (proc==9)
    {
        if (pinPDG==21) return {21,21,charge*flavor,charge*flavor,Q,Q}; 
        return {pinPDG,pinPDG,21,21,Q,Q}; 
    }
    if (proc==-9)
    {
        if (pinPDG==21) return {21,charge*flavor,charge*flavor,21,G,Q}; 
        return {pinPDG,21,21,pinPDG,G,Q}; 
    }
    if (proc==11) return {pinPDG,pinPDG,pinPDG,pinPDG,G,G}; 
    std::cout<<"q_nums failed"<<std::endl;
    exit(EXIT_FAILURE);
}

//this and the next 3 assume a picked process and sample that process's distribution 
double gen_phi(double dk,double kcm,double x,double pin,double g,double d,double proc,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx, double pinType)
{
    double l=0;
    double h=2*M_PI;
    double guess=distcon(generator) * myproc(h,dk,kcm,x,pin,g,d,3,proc,wdk,wkcm,wx, pinType);
    while ((h-l)*2/(h+l)>ep/1000)
    {
        double mid=(h+l) / 2;
        double proc_temp=myproc(mid,dk,kcm,x,pin,g,d,3,proc,wdk,wkcm,wx, pinType);
        if (proc_temp >= guess) h=mid;
        if (proc_temp < guess) l=mid;
    }
    return (h+l)/2;
}
double gen_dk(double kcm,double x,double pin,double g,double d,double proc,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx, double pinType)
{
    double l=0;
    double h=pin*10;
    double guess=distcon(generator) * myproc(0,h,kcm,x,pin,g,d,2,proc,wdk,wkcm,wx, pinType);
    //while (myproc(0,inf,kcm,x,pin,g,d,2,proc)<guess) h*=2;
    while ((h-l)*2/(h+l)>ep/1000)
    {
        double mid=(h+l) / 2;
        double proc_temp=myproc(0,mid,kcm,x,pin,g,d,2,proc,wdk,wkcm,wx, pinType);
        if (proc_temp >= guess) h=mid;
        if (proc_temp < guess) l=mid;
    }
    return (h+l)/2;
}
double gen_kcm(double x,double pin,double g,double d,double proc,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx, double pinType)
{
    double l=0;
    double h=pin;
    double guess=distcon(generator) * myproc(0,0,h,x,pin,g,d,1,proc,wdk,wkcm,wx, pinType);
    while ((h-l)*2/(h+l)>ep/1000)
    {
        double mid=(h+l) / 2;
        double proc_temp=myproc(0,0,mid,x,pin,g,d,1,proc,wdk,wkcm,wx, pinType);
        if (proc_temp >= guess) h=mid;
        if (proc_temp < guess) l=mid;
    }
    return (h+l)/2;
}
double gen_x(double pin,double g,double d,double proc,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx, double pinType)
{
    double l=0;
    double h=1;
    double guess=distcon(generator) * myproc(0,0,0,1,pin,g,d,0,proc,wdk,wkcm,wx, pinType);
    while ((h-l)*2/(h+l)>ep/1000)
    {
        double mid=(h+l) / 2;
        //cout<<h-l<<std::endl;
        double proc_temp=myproc(0,0,0,mid,pin,g,d,0,proc,wdk,wkcm,wx, pinType);
        if (proc_temp >= guess) h=mid;
        if (proc_temp < guess) l=mid;
    }
    return (h+l)/2;
}
std::vector<double> MatVec(double rot_tot[3][3], double rot_z[3][3], double vec[3])
{
    double r[3]={0,0,0},i[3]={0,0,0};
    for(int k=0; k<3;++k) for(int l=0;l<3;++l) i[l]+=rot_z[l][k]*vec[k];
    for(int k=0; k<3;++k) for(int l=0;l<3;++l) r[l]+=rot_tot[l][k]*i[k];
    return {r[0],r[1],r[2]}; 
}
/*
double Dot(double* V1,double* V2)
{
    return V1[0]*V2[0]+V1[1]*V2[1]+V1[2]*V2[2];
}

double* Cross(double* V1,double* V2)
{
    double final[3]{ V1[1]*V2[2]-V1[2]*V2[1],V1[2]*V2[0]-V1[0]*V2[2],V1[0]*V2[1]-V1[1]*V2[0]};
    
}
double*Diff(double* V1,double* V2)
{
    double final[3]{V1[0]1V2[0],V1[1]*V2[1],V1[2]*V2[2]};
    
}
*/
std::vector<double> kinematics(double pinx, double piny, double pinz, double g, double d, double proc,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx, double pinType) // after we calculate that a new particle is necessary,this generates it at that point
{
    double pin=sqrt(pinx*pinx+piny*piny+pinz*pinz);

    double x, kcm, dk, phi;
    double theta, k, p, q;
    double utwid;

    while (true) {

      //std::cout << " Generating x " << std::endl;
      x=gen_x(pin,g,d,proc,wdk,wkcm,wx, pinType);
      //std::cout << " Generating kcm " << std::endl;
      kcm=gen_kcm(x,pin,g,d,proc,wdk,wkcm,wx, pinType);
      //std::cout << " Generating dk " << std::endl;
      dk=gen_dk(kcm,x,pin,g,d,proc,wdk,wkcm,wx, pinType);
      //std::cout << " Generating phi " << std::endl;
      phi=gen_phi(dk,kcm,x,pin,g,d,proc,wdk,wkcm,wx, pinType);

      theta=theta_(kcm,x,pin);
      k=k_(dk,kcm,x,pin);
      p=p_(kcm,x,pin);

      //double AR=abs(R_(theta,pin));
      //double ARB=abs(RBar_(theta,pin));
      q=q_(kcm,x,pin);

    

      utwid=amD2/x*(Cu_(dk,kcm,x,pin)-D_(dk,kcm,x,pin)*cos(phi));

      if (utwid<amD2)
      {
        std::cout<<"utwid still somehow negative,very bad"<<std::endl;
	std::cout << " utwid= " << utwid << " amD2= " << amD2 << std::endl;
	std::cout << "x= " << x << " dk= " << dk << " kcm= " << kcm << " phi= " << phi << " pin= " << pin << " Cu_= " << Cu_(dk,kcm,x,pin) << " D_= " << D_(dk,kcm,x,pin) << std::endl;
        std::cout << "Retrying" << std::endl;
	continue;

	//exit(EXIT_FAILURE);
      }

      break;
    }

    //cout<<"phi="<<phi<<std::endl;

    double phi_k=distcon(generator)*(2 * M_PI);
    //cout<<"phi_k="<<phi_k<<std::endl;
    double E4=pin+k-p;
    //double coinflip=distint(generator)*2-1; //y has an arbitrary sign 
    

    // double px=p*sin(theta);
    // double py=0;
    // double pz=p*cos(theta);
    
    double ux=pinx/sqrt(2*pin*(pin+pinz));//these are for the rotation matrix to get it back to the experiment frame
    double uy=piny/sqrt(2*pin*(pin+pinz));
    double uz=(pin+pinz)/sqrt(2*pin*(pin+pinz));
//cout<<"hello"<<std::endl;
    double ConSqrt=sqrt(k*E4/p/pin-sin(theta/2)*sin(theta/2));
    double kx=-2*cos(phi)*sin(theta/2)*(-pin+p*cos(theta))*p*pin/q/q*ConSqrt+p*sin(theta)/2*(1+(k*k-E4*E4)/q/q);
    //cout<<coinflip<<"  cf"<<std::endl;
    double ky=sin(phi)*2*sin(theta/2)*p*pin/q*ConSqrt;
    double kz=(p*cos(theta)-pin)/2*(1+(k*k-E4*E4)/q/q)+2*cos(phi)*p*p*pin/q/q*ConSqrt*sin(theta/2)*sin(theta);
    //cout<<"k_check1 "<<kx*kx+ky*ky+kz*kz-k*k<<std::endl;
    double E4x=kx-p*sin(theta);
    double E4y=ky;
    double E4z=kz-p*cos(theta)+pin;
    //cout<<"sqrtcheck"<<k*k-p*p+pow(p*pin/q*sin(theta),2)<<std::endl;
    //cout<<"E4check1 "<<E4x*E4x+E4y*E4y+E4z*E4z-E4*E4<<std::endl;
    //cout<<"sumcheck "<<pin+k-p-E4<<std::endl;

    //double pinj[3]={0,0,pin};
    double pj[3]={p*sin(theta),0,p*cos(theta)};
    double kj[3]={kx,ky,kz};
    double E4j[3]={E4x,E4y,E4z};

    //double test1[3] {1,2,3};
    //double test2[3]{4,5,6};
    //double pinjt[3]{0,0,pin};
    //double pjt[3]{p*sin(theta),0,p*cos(theta)};
    //double kjt[3]{kx,ky,kz};
    //double E4jt[3]{E4x,E4y,E4z};

    //cout<<Diff(test1,test2)[0]<<Diff(test1,test2)[1]<<Diff(test1,test2)[2]<<"hi"<<std::endl;
    //cout<<0<<"  "<<p*pin*sin(theta)<<"   "<<0<<"  "<<std::endl;
    //cout<<Cross(pj,Diff(pj,pinj))[1]<<std::endl;
    //cout<<pj[0]<<"  "<<pj[2]<<std::endl;
    //double top=Dot(Cross(pjt,Diff(pjkcm,x,pinjt)),Cross(Diff(pjkcm,x,pinjt),kjt));
    //double bot=sqrt(Dot(Cross(pjt,Diff(pjkcm,x,pinjt)),Cross(pjt,Diff(pjkcm,x,pinjt)))*Dot(Cross(Diff(pjkcm,x,pinjt),kjt),Cross(Diff(pjkcm,x,pinjt),kjt)));
    //cout<<cos(phi)<<"   "<<top/bot<<std::endl;

    double rot_tot[3][3]={ {-1+2*ux*ux,2*ux*uy,2*ux*uz},
                             {2*ux*uy,-1+2*uy*uy,2*uy*uz},
                             {2*ux*uz,2*uy*uz,-1+2*uz*uz} };

    double rot_z[3][3]={{cos(phi_k),-sin(phi_k),0},
                        {sin(phi_k),cos(phi_k),0},
                        {0,0,1}};//make sure to apply this one first!!
    //double pin_x_teskcm,x,pin_y_teskcm,x,pin_z_test;
    //double pintest[3]={pinx,piny,pinz};
    //pin_x_teskcm,x,pin_y_teskcm,x,pin_z_test=MatVec(rot_tot,rot_z,pintest);
    //cout<<pin_x_test<<"  pin test"<<std::endl;
    std::vector<double> kv=MatVec(rot_tot,rot_z,kj);
    std::vector<double> pv=MatVec(rot_tot,rot_z,pj);
    std::vector<double> E4v=MatVec(rot_tot,rot_z,E4j);
    //cout<<(px*pinx+py*piny+pz*pinz)/p/pin-cos(theta)<<std::endl;
    //cout<<px<<std::endl;
    //cout<<"pcheck"<<px*px+py*py+pz*pz-p*p<<std::endl;
    //cout<<"k_check2"<<kx*kx+ky*ky+kz*kz-k*k<<std::endl;
    //cout<<"E4check"<<E4x*E4x+E4y*E4y+E4z*E4z-E4*E4<<std::endl;

    return {kv[0],kv[1],kv[2],pv[0],pv[1],pv[2],E4v[0],E4v[1],E4v[2]}; 
}


//This is the only function you will need to call,as it runs everything. Please temperature normalize every input parameter.
std::vector<double> gen_particles(double pinx, double piny, double pinz, double pinPDG, double deltime,gsl_integration_workspace *wdk,gsl_integration_workspace *wkcm,gsl_integration_workspace *wx) // split_or_nah,timechange,px,py,pz,pPDG,kx,ky,kz,kPDG,E4x,E4y,E4z,E4PDG
{
    
    //cout<<"starting gen_particles"<<std::endl;
    double pin=sqrt(pinx*pinx+piny*piny+pinz*pinz);

    double pQGQB;
    std::vector<double> out = is_split(pin, deltime, pinPDG,wdk,wkcm,wx);
    double split_or_nah=out[0]; double timechange=out[1]; double myF=out[2]; double myproc=out[3]; double pinType=out[4];
    //cout<<"is_split has been performed"<<std::endl;
    //check highest n
    //for (int ii=0; ii<7; ii++) {
      //std::cout << " max[" << ii+1 << "] = " << max_n[ii] << std::endl;
    //}
    //
    if (split_or_nah==-1) return {-1,timechange,0,0,0,0,0,0,0,0,0,0,0,0}; 
    if (split_or_nah==0) return {0,timechange,0,0,0,0,0,0,0,0,0,0,0,0};
    std::cout << " GOT SPLITTING " << std::endl;
    //cout<<myproc<<std::endl;
    if (myF==1||myF==4) pQGQB=1;
    if (myF==2||myF==5) pQGQB=-1;
    if (myF==3) pQGQB=2;
    std::vector<double> out2=q_nums(myproc,pinPDG,pQGQB);//pinPDG2,pPDG,kPDG,E4PDG,QG1,QG2
    double pPDG=out2[1]; double kPDG=out2[2]; double E4PDG=out2[3]; double QG1=out2[4]; double QG2=out2[5];
    //cout<<"q_nums has been performed "<<pinPDG2<<"  "<<pPDG<<"  "<<kPDG<<"  "<<E4PDG<<"  "<<QG1<<"   "<<QG2<<std::endl;
    //cout<<pinx<<"  "<<piny<<"  "<<pinz<<"  "<<QG1<<"  "<<QG2<<"   "<<myproc<<std::endl;
    std::vector<double> out1=kinematics(pinx,piny,pinz,QG1,QG2,myproc,wdk,wkcm,wx, pinType);
    double kx=out1[0];double ky=out1[1]; double kz=out1[2]; double px=out1[3]; double py=out1[4]; double pz=out1[5]; double E4x=out1[6]; double E4y=out1[7]; double E4z=out1[8];
    if (px==0&&py==0&&pz==0) return {0,timechange,0,0,0,0,0,0,0,0,0,0,0,0}; 
    //cout<<"kinematics has been performed"<<std::endl;
    return {1,timechange,px,py,pz,pPDG,kx,ky,kz,kPDG,E4x,E4y,E4z,E4PDG}; 
}

//TODO
//check q_nums
//clean up mathematica notebook
//turned off annoying sums is that ok
//fix canisample qualF=1 p dist
//fix qnums to make sure twiddle is understood
//fix sampling errors
//charm quark 4,5,6,9,10,also Nfs 
