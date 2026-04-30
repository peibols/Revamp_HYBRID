#include <iostream>
#include "gsl/gsl_sf_expint.h"
#include "gsl/gsl_sf_hyperg.h"
#include "gsl/gsl_sf_gamma.h"
#include "gsl/gsl_sf_bessel.h"
#include "gsl/gsl_math.h"
#include "gsl/gsl_integration.h"
#include "gsl/gsl_errno.h"
#include <stdlib.h>

//#include <omp.h>
#include <random>
#include <vector>
#include <complex>
#include <cstdlib>

// #include "gsl/gsl_sf_laguerre.h"
double inf = std::numeric_limits<double>::infinity();
double ep=pow(10.0,-1.2);
//double psi_res=1000000;
// double kcmres=100;
// double x_res=100;

double abs_res_kcmx=pow(10,-1);
double rel_res_kcmx=pow(10,-1);
// double abs_res_dk=pow(10,-3);
// double rel_res_dk=pow(10,-3);
double abs_res_dk=pow(10,-15);
double rel_res_dk=pow(10,-10);


double small_number = std::pow(10,-7);

// 1.00021 1.00024 1.00027 2.89406 -2615.45 13.7553 286.85
// 1.00021 1.00024 1.00027 2.77496 -2450.89 12.9503 268.808
// 1.00021 1.00024 1.00027 2.76304 -2434.42 12.8697 267.003
//1.00021 1.00024 1.00027 2.76185 -2432.78 12.8617 266.822
//  double R_(double t, double pin)
//  {
//      return -pin * (1+cos(t))/2.0;
//  }
//  double RBar_(double t,double pin)
//  {
//      return -pin *(1-cos(t))/2.0;
//  }
//  double q_(double t, double p,double pin)
//  {
//      return sqrt(pow(pin-p*cos(t),2)+pow(p*sin(t),2));
//  }
//  double f_(double t, double p,double pin)
//  {
//      return (p+pin+q_(t,p,pin))/2;
//  }
//  double F_(double t, double p,double pin)
//  {
//      return f_(t,p,pin)+R_(t,pin);
//  }
//  double s_(double t, double p,double pin)
//  {
//      return (q_(t,p,pin)-pin-p)/2;
//  }
//  double S_(double t, double p,double pin)
//  {
//      return s_(t,p,pin) -R_(t,pin);
//  }

double choose(double n,double k)
{
    
    if (k==0 or n==k) return 1;
    if (k==1) return n;
    if (n==-1 && int(k)==k) return pow(-1,k);
    // std::cout<<n<<" start  "<<k<<std::endl;
    // if(k>n&&n>0&&int(n)==n&&int(k)==k) return 0;
    if (-abs(int(n-k+1))==n-k+1&&-abs(int(n+1))!=n+1&&-abs(int(k+1))!=k+1) return 0;
    if (n==0&&k!=0) return sin(M_PI*k)/M_PI/k;
    // std::cout << "choose" << n << "   " << k << std::endl;
    // if (n==.5&&k>0&&k<30&&int(k)==k) return gsl_sf_choose(2*k,k)*pow(-1,k+1)/pow(2,2*k)/(2*k-1);
    // if (n>0&&k>0&&int(n)==n&&int(k)==k) return gsl_sf_choose(n,k);
    // if (n>0&&k>0&&int(k)==k&&int(n)!=n) return sin((n-k)*M_PI)/M_PI/(n-k)/choose(k,n);
    if (abs(n)<170&&abs(k)<170&&abs(n-k+1)<170&&(int(n)!=n||n>=1)&&(int(k)!=k||k>=1)) 
    {
        // std::cout<<gsl_sf_gamma(n+1)/gsl_sf_gamma(k+1)/gsl_sf_gamma(n-k)<<std::endl;
        // std::cout<<n<<" start  "<<k<<std::endl;
        double out=gsl_sf_gamma(n+1)/gsl_sf_gamma(k+1)/gsl_sf_gamma(n-k+1);
        // std::cout << n << "  finished " << k << std::endl;
        return out;
    }
    if (abs(k)>4*abs(n-k)&&abs(n)>4*abs(n-k)&&int(n-k)==n-k) return choose(n,n-k);
    if (int(k)==k &&k>0&&int(n)!=n) return n/k*choose(n-1,k-1);
    if (int(k)==k &&k<0&&k>-10&&n!=-1) return (k+1)/(n+1)*choose(n+1,k+1);

    // 
    

    // if (n>=1&&k>=1)
    // {
    //     if (abs(n-k)<20||n<k) return n/k*choose(n-1,k-1);
    //     // if (abs(n-k)>10) return choose(n,n-k);
    //     // {
    //     //     if (abs(k)<10) return n/(n-k)*choose(n-1,k);
    //     //     if (abs(n)<10) return 
    //     // }
        
    // }

    if (abs(n-k)==n-k&&n>k&&5>n-k) return n/(n-k)*choose(n-1,k);
    if (abs(n)<3*k||(int(k)==k&&k>0)) return n/k*choose(n-1,k-1);
    if (n>0&&k>0&&n>k&&n-k>k) return choose(n,n-k);
    if (n==int(n)&&n>0) return n/(n-k)*choose(n-1,k);
    if (n<0&&k<0) return (k+1)/(n+1)*choose(n+1,k+1);
    //if (n==.5&&k>=0) return (n-k+1)/k*choose(n,k-1);
    //if (abs(n-k)<5) return 


    // if ((n>0&&n<1&&k>0)||(k>0&&int(k)==k&&abs(n)<10)) return (n-k+1)/k*choose(n,k-1);
    
    // if (int(n-k)==n-k&&n>k) return ;

    // if (n<0||k<0)
    // {
    //     return (k+1)/(n+1)*choose(n+1,k+1);
    // }

    // if (int(n)==n&&int(k)!=k&&k==.5) return n/(n-k)*choose(n-1,k);

    // if (k==int(k))
    // {
    //     if (k>0&&n>0&&k>300+n) return choose(n,n-k);
    //     if (k>0) return n/k*choose(n-1,k-1);
    //     if (k<0) return (k+1)/(n+1)*choose(n+1,k+1);
    // }
    std::cout << "choose issue " << n << "   " << k << "  "<<n-k<< std::endl;
    exit(EXIT_FAILURE);
    // if (k>n&&n>=0) return 0;
    // if (n<k && n<0) return pow(-1,k)*choose(-n+k-1,k);
    
    
}

double chooserat(double a, double b, double j)
{
    if ((a>=0||int(a)!=a)&&b<0&&int(b)==b) return 0;
    
    double t1=choose(a,j);
    double t2=choose(b,j);
    if ((a<50&&b<50&&j<50)||j==0||a==0||b==0||(t1==t1&&t2==t2))
    {
        // std::cout<<a<<"   "<<b<<"    "<<j<<"   "<<choose(a,j)/choose(b,j)<<std::endl;
        return t1/t2;
    }
    
    return a/b*chooserat(a-1,b-1,j-1);
    std::cout<<"chooserat issue "<< a<< "   "<<b<<"  "<<j<<std::endl;
    exit (EXIT_FAILURE);

}

double expEn(double n,double x)//This will only work for n>=0
{
    if (x==inf) return 0;

    if (x==0 && n>1) return 1/(n-1);
    if (x==0 && n<=1) return std::pow(10.,50.);
    //if (x==0 && n<=1) return inf;

    if (n==0) return 1/x;
    // if ((n==1||n==2)&&x>500) std::cout<<"ExpEn "<< n<< "   "<<x<<std::endl;
    if (n==1) return gsl_sf_expint_E1_scaled(x);
    if (n==2) return gsl_sf_expint_E2_scaled(x);
    if (x>0&&x<=500) return gsl_sf_expint_En_scaled(n, x);
    if (n > 1 &&n<4) return (1 - x * expEn(n - 1, x)) / (n - 1);
    std::cout<<"hoping I don't have to go down here"<<x<<"   "<<n<<std::endl;
    if (n<0) return 1/x-n/x*expEn(n+1,x);//not accurate for large n
    //if (n==1 && x>100) return .5*(.5*log(1+2/x)+log(1+1/x));
    //std::cout<<"hoping I don't have to go down here"<<x<<"   "<<n<<std::endl;
    //if ((n>50&&x>0)) return 1/(x+n)*(1+n*pow(x+n,-2)+n*(n-2*x)*pow(x+n,-4)+n*(6*x*x-8*n*x+n*n)*pow(x+n,-6));
    if ((x>500&&n>=3)) return .5/(x+n)+.5/(x+n-1);
   std::cout<<"ExpEn issue "<< n<< "   "<<x<<std::endl;
    exit (EXIT_FAILURE);
}

double ExpCorr(double p) //prob that at least 1 thing happens
{
    return 1 - exp(-p);
}

template <typename T> double sgn(T val)
{
    return (T(0) < val) - (val < T(0));
}

double fact(double n)//only works for n>0
{
    // if (n!=0&&n!=1&&n!=2&&n<150) std::cout<<"mygamma " <<n+1<<std::endl;
    if (n<150) return gsl_sf_gamma(n+1);
    if (n==1||n==0) return 1;
    if (n==.5) return sqrt(M_PI)/2.;
    if (n<0)
    {
        std::cout << "fact issue" << std::endl;
        exit(EXIT_FAILURE);
    }
    return n*fact(n-1);
}



double expGamma(double n, double x)
{
    if (n==0 &&x==0) 
    {
        std::cout<<n<<" Gamma00inf  "<<x<<std::endl;
        //return inf;
	return pow(10.,50.);
    }
    if (x==0 && n-1>0) return fact(n-1);
    if (n==1) return 1;
    if (n>1) return (n-1)*expGamma(n-1,x)+pow(x,n-1);
    if (n>=-3 && x!=0.) return expEn(1-n,x)*pow(x,n);
    if (n>=-3 && x==0.) return pow(10.,50.);
    std::cout<<"expGamma issue"<<std::endl;
    exit(EXIT_FAILURE);
}

double Gamma(double n, double x)
{
    if (x==inf) return 0;
    return exp(-x)*expGamma(n,x);
    std::cout<<"Gamma issue"<<std::endl;
    exit(EXIT_FAILURE);
}

double powdiff(double x,double y,double n)
{
    if (abs(x/y)>2||abs(y/x)>2) return (pow(x,n)-pow(y,n))/(x-y);
    double sum=0;
    for (int l = 0; l <= n - 1; l++) sum+=pow(x,l)*pow(y,n-1-l);
    return sum;
}

double xofact(double x,double n,double expcorr=0)
{
    if (n==0) return exp(expcorr);
    if (x<0) return xofact(-x,n,expcorr)*pow(-1,n);
    if (expcorr!=0) 
    {
        // std::cout << "xofactexpcorr " <<x<<"   "<<n<<"  "<<expcorr<<"   "<<x*exp(expcorr/n)<< std::endl;
        return xofact(x*exp(expcorr/n),n);
    }
    // std::cout << "xofactexpcorr " <<x<<"   "<<n<<"  "<<expcorr<<"   "<<x*exp(expcorr/n)<< "  "<<exp(n*log(x)-gsl_sf_lngamma(n+1))<<std::endl;
    //  std::cout << "xofactexpcorr " <<x<<"   "<<n<<"  "<<expcorr<<"   "<<x*exp(expcorr/n)<< std::endl;
    return exp(n*log(x)-gsl_sf_lngamma(n+1));
    //  std::cout << "xofactexpcorr " <<x<<"   "<<n<< std::endl;
    // if (n!=0) return x/n*xofact(x,n-1);
    std::cout << "more xofact cases" << std::endl;
    exit(EXIT_FAILURE);
}

// double xofactsum(double x,double n)
// {
//     double sum=0;
//     double store=1;
//     for (int l = 0; l <= n; l++) 
//     {
//         sum+=store;
//         store*=x/n;
//     }
//     return sum;
// }

// double lag(double k, double a, double x,double expcorr=0)
// {
//     // std::cout << "lag " <<k<<"   "<<a<<"  "<<x<<"   "<< std::endl;
//     if (k==0) return exp(expcorr);
//     if (k==1) return exp(expcorr)*(1+a-x);
//     if (a==-k+2) return (xofact(-x*exp(expcorr/k),k)*(pow(x-k,2)-k)/x/x);
//     if (a==-k+1) return ((x-k)/(x)*xofact(-x*exp(expcorr/k),k));
//     if (a==-k) return (xofact(-x*exp(expcorr/k),k));
//     if (a==-k-1) 
//     {

//         double p1= xofactsum(x,k)*exp(expcorr)*pow(-1, k);
//         // double p1=exp(expcorr)*pow(-1, k);
//         // for (int i = 1; i <= k; ++i) p1 += xofact(x*exp(expcorr/i),i)*pow(-1, k);
//         if (p1 == inf || p1 != p1)
//         {
//             std::cout << "lag failure " <<k<<"   "<<a<<"  "<<x<<"   "<< expcorr<<"   "<< p1<< std::endl;
//             exit(EXIT_FAILURE);
//         }
//         return p1;
//     }
//     if (a==-k-2)
//     {
//         double p1=xofactsum(x,k)*exp(expcorr)*pow(-1, k)*(k+1);
//         double p2=xofactsum(x,k-1)*exp(expcorr)*pow(-1, k)*(-x);
//         // for (int i = 1; i <= k; ++i) p1 += xofact(x*exp(expcorr/i),i)*pow(-1, k)*(k+1);
//         // for (int i = 1; i <= k-1; ++i) p1 += -x*xofact(x*exp(expcorr/i),i)*pow(-1, k);
// ;
//         if (p1+p2 == inf || p1+p2 != p1+p2)
//         {
//             std::cout << "lag failure " <<k<<"   "<<a<<"  "<<x<<"   "<< p1+p2<< std::endl;
//             exit(EXIT_FAILURE);
//         }
//         return p1+p2;
//     }
//     if (a>=0) return gsl_sf_laguerre_n(k,a,x)*exp(expcorr);
    
//     // if (a<0) return lag(k,a+1,x)-lag(k-1,a+1,x);
//     std::cout << "explag" << "   "<<k<<"  "<<a<<"   "<<x<<std::endl;
//     std::cout << "more explag cases" << std::endl;
//     exit(EXIT_FAILURE);
//     //if (k==1) return exp(-x)*(1+a-x);
// }

double GammaQOn(double n, double x,double y=1, double z=0)
{

    //  std::cout<<n<<"   "<<x<<" gammaqon "<<y<<"  "<<z<<std::endl;
    if (n==0) return expGamma(0,x)*exp(z-x);
    if ((n>100+x&&x>0)||(x<0&&n>-4*x&&n>30)) return 1./n*pow(exp(z/n)*y,n);
    //std::cout << "cases" << std::endl;
    if (x>=0) 
    {
        // std::cout<<log(log(gsl_sf_gamma_inc_Q(n,x)/n)+n*log(y)+z)<<std::endl;
        // std::cout << "cases1" << std::endl;
        double store=log(gsl_sf_gamma_inc_Q(n,x)/n)+n*log(y)+z;
        // std::cout << "cases2" << std::endl;
        if (store==-inf) return 0;
        if (store!=store)
        {
            std::cout << "nanGammaQnan" << std::endl;
    exit(EXIT_FAILURE);
        }
        return exp(store);
    }
    
    //  if (expGamma(n,x)*xofact(y,n,-x+z)==inf) std::cout << "case "<< n <<"  "<< x<<"  "<<x*y<<"  "<<xofact(y*x,n,-x+z)/x<< std::endl;
     if (n>100) return exp(z)*y*(n-1)/n*GammaQOn(n-1,x,y)+xofact(y*x,n,-x+z)/x;
     std::cout << "more GammaQ cases" << std::endl;
    return expGamma(n,x)*xofact(y,n,-x+z)/n;
    std::cout << "more GammaQ cases" << std::endl;
    exit(EXIT_FAILURE);
}

/*
double lagsum_f(double mRof, double a, double mRbl, dou4le eta, double fnl)
{
    double sum=0;
    double inc = 1000;
    double pinc=0;
    //double N=0;
    for (int k = 0; (abs(inc) >= ep * abs(sum));k++)// || (abs(pinc) <= abs(inc)); k++)
    {   
        //pinc=inc;
        inc = pow(mRof,k)*explag(k,a,mRbl)*expEn(k+eta,fnl);
        //cout << k + eta<<"  "<<fnl<<"   "<<expEn(k+eta,fnl)<<endl;
        sum += inc;
        //N+=1;
        if(k==5&&inc==0&&pinc==0&&sum==0) break;
    }
    //cout<<N<<endl;
    return sum;
}

double lagsum_s(double a, double mRbnl, double fl, double mRopin,double pinof)
{
    //cout<<f<<"     "<<pin<<endl;
    double sum = 0;
    double inc = 1000;
    double pinc= 0;
    //double N=0;
    for (int m = 0; (abs(inc) >= ep * abs(sum));m++)// || (abs(pinc) <= abs(inc)); m++)
    {
        //pinc=inc;
        inc=0;
        for (int k = 0; k<=m; k++) inc+=explag(k,1,mRbnl)*pow(mRopin,k)*pow(m+1-k,a);
        //the above exponential has the Rnl factor
        //cout << "The terms are " <<explag(m,1+a,mRbnl) <<"     "<< pow(mRopin,m)<<"    "<<pow(pin/f,m)<<"     "<<expEn(m+3+a,fl)<<"   "<<m<<"   "<<inc<<endl;
        inc*=expEn(m+3+a,fl)*pow(pinof,m);
        sum += inc;
        //N += 1;
        if (m == 5 && inc == 0 && pinc ==0 && sum == 0) break;
    }
    //cout << N <<"    s"<< std::endl;
    //cout<<"full monty "<<sum<<endl;
    return sum;
}
*/
double h(double y, double z)
{
    if (y==M_PI) return M_PI;
    if (fabs(z-1.)<pow(10.,-7.)) return M_PI;
    if (z>=1) return 2*atan(sqrt((z+1)/(z-1))*tan(y/2.));
    std::cout << "h has z<1" << std::endl;
    exit(EXIT_FAILURE);
    // if (z<1) return -2*atanh(sqrt((z-1)/(z+1))*tan(y/2.))

}

double g(double y, double z)
{
    if (z<1)
    {
        std::cout << "g has z<1" << std::endl;
        exit(EXIT_FAILURE);
    }
    return 2*M_PI*floor(y/2./M_PI+.5)+h(y,z);
}


double ABo2pi(double n, double A, double B, double x) //indef integral for the phi integrals, n is the power of the integral
{
    if (B>A &&abs(1-A/B)<pow(10,-10)) B=A;
    //if (n==-2 && x<=M_PI) return (2*A*atan((A+B)/sqrt(A*A-B*B)*tan(x/2))/pow(A*A-B*B,1.5)+B*sin(x)/((A*A-B*B)*(A-B*cos(x))))/(2*M_PI);
    if (n==-2) return (g(x,A/B)*A/pow(A*A-B*B,1.5)+B*sin(x)/((A*A-B*B)*(A-B*cos(x))))/(2*M_PI);

    //if (n==-1 && x<=M_PI) return (atan((A+B)/sqrt(A*A-B*B)*tan(x/2))/sqrt(A*A-B*B))/M_PI;
    if (n==-1) return g(x,A/B)/sqrt(A*A-B*B)/(2*M_PI);

    if (n==0) return x/(2*M_PI);

    if (n==1) return (A*x-B*sin(x))/(2*M_PI);
    
    if (n==2) return (A*A*x-2*A*B*sin(x)+B*B*(x/2.+sin(2*x)/4.))/(2.*M_PI);
   std::cout << "ABo2pi issue" << std::endl;
    exit (EXIT_FAILURE);
}

double ZifZ(double x, double y) //zero if zero
{
    if (x==0) return 0;
    if (x!=0) return x*y;
   std::cout << "ZifZ issue" << std::endl;
    exit(EXIT_FAILURE);
}

double poch(double a,double b)
{
    //  std::cout << a<<" | "<<b << std::endl;
    return gsl_sf_poch(a,b);
    // std::cout << a<<" || "<<b << std::endl;
    // if (b==0) return 1;
    // if (a<=0) return a/(a+b)*poch(a+1,b);
    // return gsl_sf_gamma(a+b)/gsl_sf_gamma(a);
}

double F21(double a,double b,double c,double z)
{
      //std::cout << "F21 "<<a<<"  "<<b<<"   "<<c<<"    "<<z<< std::endl;
    if (z==1) 
    {
        // std::cout << "hi" << std::endl;
        // if (c>a+b) return gsl_sf_gamma(c)*gsl_sf_gamma(c-a-b)/gsl_sf_gamma(c-a)/gsl_sf_gamma(c-b);
        // if (c>a+b) return choose(c-1,a)/choose(c-b-1,a);
        //  std::cout << "2F1 z==1 " << c-1 << "  " << c-b-1 << "   " << a <<"   "<<chooserat(c-1,c-b-1,a)<< std::endl;
        if (c>a+b) return chooserat(c-1,c-b-1,a);
        std::cout << "2F1 issue " << a << "  " << b << "   " << c << "    " << z << std::endl;
        exit(EXIT_FAILURE);
    }
    if (z==0||a==0||b==0) return 1;
    if (b==-1&&c!=0) return 1-a*z/c;
    if (a==c) return pow(1-z,-b);

    // if (b<0&&int(b)==b&&-b<50) 
    // {
    //     double start=0;
    //     for (int j = 0; j<=-b; j++) start+=pow(-z,j)*choose(-b,j)*chooserat(a+j-1,c+j-1,j);
    //     return start;
    // }
    //  std::cout << "continuing" << std::endl;
    //if (abs(a)<1000&&abs(b)<1000&&abs(c)<1000) 
    if (abs(z)>1)
    {
        std::cout << "2F1 issue2 " << a << "  " << b << "   " << c << "    " << z << std::endl;
        exit(EXIT_FAILURE);
    }


    
    if ((z<.99&&b>-2)||(a>0&&b>0&&c>0&&b<2)) 
    {
         //std::cout << "using 2f1 " << a << "  " << b << "   " << c << "    " << z <<"  "<<gsl_sf_hyperg_2F1(a, b, c, z)<<std::endl;
        double out=gsl_sf_hyperg_2F1(a, b, c, z);
        //  std::cout << "finished 2f1 "<< std::endl;
        return out;
    }
    // if (c>a+b&&z>.999) return F21(a,b,c,1);
    
    if (b<0&&c==1) return choose(c-a-b-1,-b)/choose(-b+c-1,-b)*F21(a,b,a+b-c+1,1-z);

     if (b<0&&a>0&&c>a&&z<.99&&b>-300) return pow(1-z,c-a-b)*F21(c-a,c-b,c,z);
    // std::cout << "continuing" << std::endl;
    // std::cout << "2F1 z!=1 " << c-1 <<  "  " << c-b-1 << "   " << a << std::endl;
    // std::cout << "2F1 z!=12 " << c - 1 << "  " << b - 1 << "   " << c-a<< std::endl;
    if (c-a-b!=int(c-a-b)&&c>a+b) 
    {
        // std::cout << "special" << std::endl;
        // std::cout << a<<"  "<< b<<"   "<<c<<"  "<<z<<"  "<< chooserat(c-1,c-b-1,a)<<"  "<<F21(a,b,a+b+1-c,1-z)<<" || "<<chooserat(c-1,b-1,c-a)<<" "<<pow(1-z,c-a-b)<<"   "<<F21(c-a,c-b,1+c-a-b,1-z)<< std::endl;
        
        if (z>.5&&b<0) 
        {
            // std::cout<<"using"<<std::endl;
            return ZifZ(chooserat(c-1,c-b-1,a),F21(a,b,a+b+1-c,1-z))+ZifZ(chooserat(c-1,b-1,c-a),pow((1-z)/z,c-b)*F21(1-b,c-b,1+c-a-b,(z-1)/z)*pow(1-z,-a));
        }
        return chooserat(c-1,c-b-1,a)*F21(a,b,a+b+1-c,1-z)+chooserat(c-1,b-1,c-a)*pow(1-z,c-a-b)*F21(c-a,c-b,1+c-a-b,1-z);
    }
    //if (b>a&&a>0) return F21(b,a,c,z);
    // if (b<4&&int(b)==b)return ((2 - 2*b + c + (-1 - a + b)*z)/((b - 1)*(z - 1)))*F21(a, b - 1, c, z) + ((-1 + b - c)/((b - 1)*(z - 1)))*F21(a, b - 2, c, z);
    // std::cout << "F21 long way2 " << a << "  " << b << "   " << c << "    " << z << std::endl;
    // if (b<-1300&&c==a+1) return pow(1-z,c-a-b)*F21(c-a,c-b,c,z);
    
    // std::cout << "F21 long way1 "<<a<<"  "<<b<<"   "<<c<<"    "<<z<< std::endl;
        double sum=0;
        double inc = 1000;
        double cr=1;
        //double N=0;
        for (int j = 0; (abs(inc)*10000. >= ep * abs(sum));j++)
        {   
            // std::cout<<alpha-1<<"   "<<j<<std::endl;
            //inc = choose(a+j-1,j)*choose(b+j-1,j)/choose(c+j-1,j)*pow(z,j);
            //chooserat(a + j - 1, c + j - 1,j)
            inc =cr*choose(b + j - 1, j) * pow(z, j);
            cr*=(a+j)/(c+j);
            //cout << k + eta<<"  "<<fnl<<"   "<<expEn(k+eta,fnl)<<endl;
            sum += inc;
              //std::cout << abs(inc)<<"    "<<ep * abs(sum)<<"  "<<"  ||| "<<choose(b + j - 1, j)<<std::endl;
            //N+=1;
            if (abs(inc)==inf||inc!=inc) 
            {
                std::cout << "F21case2 "<<j<<"   "<<a<<"  "<<b<<"    "<<c<<"    "<<z<<std::endl;
                exit(EXIT_FAILURE);
            }
            if(j==5&&inc==0&&sum==0) break;
        }
        return sum;


    // if ((c-a==int(c-a)&&c>a&&a>0)&&(b<0||int(b)==b)) 
    // {
    //     // std::cout << (c-a==int(c-a)&&c>a)<<"  "<<(c-b!=int(b-a))<<"   "<<(c<b)<<"    "<< std::endl;
    //     //std::cout << "flipping" << std::endl;
    //     return F21(b,a,c,z);
    // }
    

    // if (c==1+b&&a<0)
    // {
    //     // double alpha=b;
    //     // double m=-a;
    //     double sum=0;
    //     double inc = 1000;
    //     //double N=0;
    //     for (int j = 0; (abs(inc) >= ep * abs(sum));j++)
    //     {   
    //         // std::cout<<alpha-1<<"   "<<j<<std::endl;
    //         inc = pow(1-z,j-a+1)/(a-j-1)*pow(-1,j)*choose(b-1,j);
    //         //cout << k + eta<<"  "<<fnl<<"   "<<expEn(k+eta,fnl)<<endl;
    //         sum += inc;
    //         // std::cout << z-1<<"    "<<m+j+1<<"  "<<pow(z-1,m+j+1)<<"   "<<inc<<std::endl;
    //         //N+=1;
    //         if (abs(inc)==inf||inc!=inc) 
    //         {
    //             std::cout << "F21case2 "<<j<<"   "<<a<<"  "<<b<<"    "<<c<<"    "<<z<<"    "<<pow(1-z,j-a+1)/(a-j-1)<<"  ||||| "<<pow(-1,j)*choose(b-1,j)<<"    " <<inc<< std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    //         if(j==5&&inc==0&&sum==0) break;
    //     }
    //     // std::cout<<"finished "<<sum<<std::endl;
    //     return b/pow(z,b)*(sum+1./(choose(b-a,b)*b));
        
    // }


    // if (c==1+b&&a>=1&&a<=4&&b>=0)
    // {
    //     // double alpha=b;
    //     // double m=-a;
    //     double sum=0;
    //     double inc = 1000;
    //     //double N=0;
    //     for (int j = 0; (abs(inc) >= ep * abs(sum));j++)
    //     {   
    //         // std::cout<<alpha-1<<"   "<<j<<std::endl;
    //         inc = choose(a+b-1,a)/choose(a+b+j-1,a)*choose(a+j-1,a-1)*pow(z,j);
    //         //cout << k + eta<<"  "<<fnl<<"   "<<expEn(k+eta,fnl)<<endl;
    //         sum += inc;
    //         // std::cout << z-1<<"    "<<m+j+1<<"  "<<pow(z-1,m+j+1)<<"   "<<inc<<std::endl;
    //         //N+=1;
    //         if (abs(inc)==inf||inc!=inc) 
    //         {
    //             std::cout << "F21case2 "<<j<<"   "<<a<<"  "<<b<<"    "<<c<<"    "<<z<<"    "<<pow(1-z,j-a+1)/(a-j-1)<<"  ||||| "<<pow(-1,j)*choose(b-1,j)<<"    " <<inc<< std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    //         if(j==5&&inc==0&&sum==0) break;
    //     }
    //     // std::cout<<"finished "<<sum<<std::endl;
    //     return sum;
        
    // }


    // if (c-b==int(c-b)&&c-b!=1)
    // {
    //     double sum=0;
    //     double inc = 1000;
    //     //double N=0;
    //     double d=c-b;
    //     for (int j = 0; (j<=d-1);j++)
    //     {   
    //         // std::cout<<b<<"   "<<c<<"  "<<j<<std::endl;
    //         inc = pow(-1,j)*choose(d-1,j)*choose(c-1,d-1)*b/(j+b)*F21(a,j+b,j+b+1,z);
    //         //cout << k + eta<<"  "<<fnl<<"   "<<expEn(k+eta,fnl)<<endl;
    //         sum += inc;
    //         if (abs(inc)==inf||inc!=inc) 
    //         {
    //             std::cout << "F21case3 "<<a<<"  "<<j+b<<"    "<<j+b+1<<"    "<<z<<"    "<<sum<<"    " <<inc<< std::endl;
    //             exit(EXIT_FAILURE);
    //         }
    //     }
    //     //  std::cout<<"finished "<<sum<<std::endl;
    //     return sum;
        
    // }
    // std::cout<<"fin2ished "<<gsl_sf_hyperg_2F1(a,b,c,z)<<std::endl;
    std::cout << "F21 issue "<<a<<"  "<<b<<"   "<<c<<"    "<<z<< std::endl;
    exit(EXIT_FAILURE);
    // return gsl_sf_hyperg_2F1(a,b,c,z);
}

double Beta(double x,double a,double b)
{   
    if (x==0) return 0;
    // if (x<.01) return pow(x,a)/a*F21(a,1-b,a+1,x);


    if (x*b<.01) return pow(x,a)*(1./a-((-1+b)*x)/(1+a));
    if (a*log10(x)-log10(a)<-300) return pow(x,a)/a;



    //std::cout<<"running beta "<< x<<"   "<<a  <<"   "<<b<<"  "<<gsl_sf_beta_inc(a,b,x) * gsl_sf_beta(a,b)<<std::endl;
    return gsl_sf_beta_inc(a,b,x) * gsl_sf_beta(a,b);
}


double F1(double a,double b1,double b2,double c,double z1,double z2)
{
    if (abs(z1)>1||abs(z2)>1)
    {
        if (abs(z2-1)<pow(10,-10)) return F1(a,b1,b2,c,z1,1);
        std::cout << "F1 issue " <<z1<<"   "<<z2<< std::endl;
        exit(EXIT_FAILURE);
    }
      
    if (z1==0||z2==0) return F21(a,b1,c,z1)*F21(a,b2,c,z2);
    //double mine=0;
    // std::cout << "F1 "<<z2<<"    "<<(z2==1)<<" "<<z2-1<<std::endl;
    if (abs(z2-1)<pow(10,-15)) 
    {
        // std::cout << "F1 z2=1 "<<z2<< std::endl;
        //mine=F21(a,b2,c,1)*F21(a,b1,c-b2,z1);
        return F21(a,b2,c,1)*F21(a,b1,c-b2,z1);
    }
    if (abs(z1-1)<pow(10,-15)) 
    {
        if (F21(a,b2,c-b1,z2)<0||F21(a,b2,c-b1,z2)==inf)
        {
        std::cout << "F21 in F1 issue " <<F21(a,b2,c-b1,z2)<<"  "<<a<<"  "<<b2<<"   "<<c-b1<<"    "<<z2<< std::endl;
        exit(EXIT_FAILURE);
        }
        //mine=F21(a,b1,c,1)*F21(a,b2,c-b1,z2);
        return F21(a,b1,c,1)*F21(a,b2,c-b1,z2);
    }
    
    // std::cout << "F1 "<<a<<"  "<<b1<<"   "<<b2<<"    "<<c<<"    "<<z1<<"    "<<z2<<"    "<<mine<< std::endl;
    double sum=0;
    double inc = 1000;
    //double N=0;
    for (int r = 0; (abs(inc)*100000/ep >=abs(sum));r++)
    {   
        //std::cout << "F1 "<<a<<"  "<<b1<<"   "<<b2<<"    "<<c<<"    "<<z1<<"    "<<z2<<"    "<<r << std::endl;
        // double mine2=F21(a+r,b2+r,c+2*r,z2);
        // std::cout << "got it "<<std::endl;
        // double mine =F21(a+r,b2+r,c+2*r,z2);
        // std::cout << "got it2 "<<std::endl;
            // inc = F21(a+r,b1+r,c+2*r,z1)*F21(a+r,b2+r,c+2*r,z2)*pow(z1*z2,r)*poch(a,r)*poch(b1,r)*poch(b2,r)*poch(c-a,r)/poch(c+r-1,r)/poch(c,2*r)/fact(r);
            inc = F21(a+r,b1+r,c+2*r,z1)*F21(a+r,b2+r,c+2*r,z2)*pow(z1*z2,r)*poch(b1,r)*poch(b2,r)*poch(c-a,r)/fact(r)*chooserat(a+r-1,c+2*r-2,r)/choose(c+2*r-1,c-1)/fact(2*r);
            // std::cout<<poch(a,r)/poch(c+r-1,r)/poch(c,2*r)/(chooserat(a+r-1,c+2*r-2,r)/choose(c+2*r-1,c-1)/fact(2*r))<<std::endl;
            //check this function????????????????????
            if (inc!=inc||inc==inf) 
            {
                std::cout << "F1inc  "<<F21(a+r,b1+r,c+2*r,z1)<<"  "<<F21(a+r,b2+r,c+2*r,z2)<<"  "<<a+r<<"   "<<b1+r<<"   "<<c+2*r<<"   "<<z1<<std::endl;
                exit(EXIT_FAILURE);
            }
            //cout << k + eta<<"  "<<fnl<<"   "<<expEn(k+eta,fnl)<<endl;
            sum += inc;
            //N+=1;
        if(r==5&&inc==0&&sum==0) break;
    }
    // std::cout<<sum<<std::endl;
    if (abs(sum)==inf||sum!=sum) 
    {
        std::cout << "F1 fail"<<a<<"  "<<b1<<"   "<<b2<<"    "<<c<<"    "<<z1<<"    "<<z2<<"    " <<sum<<"  "<<inc<< std::endl;
        exit(EXIT_FAILURE);
    }
    return sum;

}


// double powcossum(double P,double M, double psi,double j,double b, double cs=0,double sn=0)
// {
//     // std::cout << "starting powcossum "<<P<<"  "<<M<<"   "<<psi<<"  "<<j<<"  "<<b<<"  "<<cs<<"   "<<sn<< std::endl;
//     //int (1-M/(P+b)*cospsi)^j*cos^cs*sin^sn
//     if (abs(P+M+b)<=M) 
//     {
//         std::cout << "powcossum R out of range" <<P<<"  "<<b<<"   "<<M<< std::endl;
//         exit(EXIT_FAILURE);
//     }
//     if (sn !=0)
//     {
//         std::cout << "powcossum sn out of range" << std::endl;
//         exit(EXIT_FAILURE);
//     }

//     double sum = 0;
//     double inc = 1000;
//     double pinc = 1000;
//     for (int l = 0; abs(pinc) >= ep / 20 * abs(sum) || (abs(inc) > abs(pinc)); l++)
//     {
//         pinc = inc;
//         // double hyper=0;
//         // std::cout << "hyper"<<(3+cs+l)/2 <<std::endl;
//         // std::cout << "hyper" <<gsl_sf_gamma((3+cs+l)/2)<< std::endl;
//         // if (cos(psi)*cos(psi)==1) hyper=gsl_sf_gamma((3+cs+l)/2)*sqrt(M_PI)/gsl_sf_gamma((2+cs+l)/2);
//         // if (cos(psi)*cos(psi)<1) hyper=gsl_sf_hyperg_2F1(.5,(1+cs+l)/2,(3+cs+l)/2,cos(psi)*cos(psi));
//         // //std::cout << "got here "<<j<<"  "<<l<<std::endl;
//         // //std::cout << "got2 here "<<choose(-j,l)<<std::endl;
//         // //std::cout << "got3 here "<<hyper<<std::endl;
//         // inc = -pow(-M/(P+b),l)*choose(j,l)*pow(cos(psi),1+cs+l)*hyper/(1+cs+l);//??????????? fix this M

//         inc =-choose(j,l)*pow(-2*M/(P+M+b),l)*Beta(cos(psi/2)*cos(psi/2),l+.5,.5);
//         // std::cout << "powcosinc "<<inc<<"  "<<j<<"    "<<l<<"   "<<choose(j,l)<<"     "<<pow(-2*M/(P+M+b),l)<<"    "<<Beta(cos(psi/2)*cos(psi/2),l+.5,.5)<<std::endl;
//         // std::cout << "got here again "<<choose(-1,1)<<std::endl;
//         sum += inc;
//         if (l == 5 && inc == 0 && sum == 0)
//             break;

//         if (inc == inf || inc != inc)
//         {
//             std::cout << "powcossum failure" << std::endl;
//             exit(EXIT_FAILURE);
//         }
//     }


//     // std::cout << "powcossum result"<< sum<<"  "<<j<< std::endl;
//     return sum;

// }

double coscos(double x, double r, double l)
{
    if (r==l&&l==0) return x-M_PI;
    if (l==r) return (x-M_PI)/2+sin(2*r*x)/4/r;
    return (l*sin(l*x)*cos(r*x)-r*cos(l*x)*sin(r*x))/(l*l-r*r);
}

double BesselI_scaled(double i, double x, double psi)
{
    //int exp(-x*cos(psi)-x)*cos(i*psi) from pi to 0
    if (psi==0) return -gsl_sf_bessel_In_scaled(i,-x);
    //if (psi==0) return -gsl_sf_bessel_In_scaled(i,-x)*M_PI;
    //if (psi==M_PI) return 0;
    if (fabs(psi-M_PI)<small_number) return 0;
    double inc=0;
    double sum=0;

    for (int l = 1; abs(inc)*10000000000/ep >=abs(sum); l++)
        {
            inc=2*gsl_sf_bessel_In_scaled(l,-x)*coscos(psi,i,l)/M_PI;
            //inc=2*gsl_sf_bessel_In_scaled(l,-x)*coscos(psi,i,l);
            sum += inc;
            // std::cout<<inc<<std::endl;
            if (l>i && inc == 0 && sum == 0)break;
        }
        //  std::cout<<i<<"   "<<x<<"  "<<psi<<"  "<<gsl_sf_bessel_I0_scaled(x)*coscos(psi,i,0)<<std::endl;
    return sum+gsl_sf_bessel_I0_scaled(x)*coscos(psi,i,0)/M_PI;
    //return sum+gsl_sf_bessel_I0_scaled(x)*coscos(psi,i,0);
}

struct expK_params { double P;double M;double j;double b; double n;};

double expK_f (double x, void * p)
{
    struct expK_params * params = (struct expK_params *)p;
    double P = (params->P);
    double M = (params->M);
    double n = (params->n);
    double j = (params->j);
    double b = (params->b);
    //std::cout << "P= " << P << "M= " << M << std::endl;
    double K=P+M*cos(x);

    return  exp(-n*(K))*pow(K+b,j)/M_PI;
}

double expKint(double P,double M, double psi,double j,double b, double n,gsl_integration_workspace *w,double a=inf,double cs=0,double sn=0 )
{
    //if (psi==0) return 0;
    //std::cout << "starting expkint" << std::endl;
    //if (psi==M_PI) return 0;
    if (fabs(psi-M_PI)<small_number) return 0;
    if (M==0) 
    {
        double cons=1;
        if (a!=inf) cons=1./(P+a);
        // std::cout << "starting expkintshort" << std::endl;
        // double mine=exp(-n*P)*pow(P+b,j)*cons*sgn(pow(cos(psi),cs-1))/(1+sn)*pow(sin(psi),1+sn)*F21(.5*(1-cs),.5*(1+sn),.5*(3+sn),sin(psi)*sin(psi))/M_PI;
        // std::cout << "starting expkintshort" <<psi << std::endl;
        //return exp(-n*P)*pow(P+b,j)*cons*sgn(pow(cos(psi),cs-1))/(1+sn)*pow(sin(psi),1+sn)*F21(.5*(1-cs),.5*(1+sn),.5*(3+sn),sin(psi)*sin(psi))/M_PI;
        return exp(-n*P)*pow(P+b,j)*cons*sgn(pow(cos(psi),cs-1))*sin(psi)*F21(.5*(1-cs),.5,1.5,sin(psi)*sin(psi))/M_PI;
    }

    double inc=0;
    //////
    // if (psi!=0)
    // {
    
    // long double mysum=0;
    //     double mypsi=M_PI;
    //     double mystep=(psi-M_PI)/psi_res;
    //     double K=P+M*cos(mypsi);
    //     for (int l = 0; l < psi_res; l++)
    //     {
    //         K=P+M*cos(mypsi);
    //         double cons=1;
    //         if (a!=inf) cons=1./(K+a);
    //         inc=exp(-n*K)*pow(K+b,j)*pow(cos(mypsi),cs)*pow(sin(mypsi),sn)*cons*mystep;
    //         mysum+=inc;
    //         mypsi+=mystep;
    //         //  std::cout<<inc<<std::endl;
    //         if (abs(inc) == inf || inc != inc||abs(mysum)==inf||mysum!=mysum)
    //             {
    //                 // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //                 // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //                 // std::cout << "expKintincfailr "<<std::endl;
    //                 std::cout << "expKint failure " <<mypsi<<"  "<<M<<"  "<<P<<"  "<<K<<"  "<<psi<< std::endl;
    //                 exit(EXIT_FAILURE);
    //             }
    //     }
    //     //std::cout << mystep << std::endl;
    //     //std::cout << inc << std::endl;
    //     return mysum/M_PI;
    // }

    //////





    if (sn > 0)
    {
        return expKint(P,M,psi,j,b,n,w,a,cs,sn-2) - expKint(P,M,psi,j,b,n,w,a,cs+2,sn-2);
        // return 2*expKint(P,M,psi,j,b,n,cs+1,sn-2,a) - expKint(P,M,psi,j,b,n,cs+2,sn-2,a);
    }
    if (cs>0) return (-(P+b)*expKint(P,M,psi,j,b,n,w,a,cs-1,sn)+expKint(P,M,psi,j+1,b,n,w,a,cs-1,sn))/M;

    //int dpsi _0^psi e^(-nK)cos^cs*(K+b)^j
    if (sn != 0||cs!=0)
    {
        std::cout << "expKint sn out of range" << std::endl;
        exit(EXIT_FAILURE);
    }
    if (a!=inf)
    {
        if (j<0 && a!=b) return expKint(P,M,psi,j,b,n,w)/(a-b)-expKint(P,M,psi,j+1,b,n,w,a)/(a-b);
        if (j<0 && a==b) return expKint(P,M,psi,j-1,b,n,w);
        if (j==0) return expKint(P,M,psi,-1,a,n,w);
	return expKint(P,M,psi,j-1,b,n,w)+(b-a)*expKint(P,M,psi,j-1,b,n,w,a);
    }

    //    std::cout << "Bessel" << -gsl_sf_bessel_I0_scaled(M*n)*exp(-n*P+n*M)<< std::endl;

    if (j==3) return (((P+b)*(P+b)+1.5*M*M)*(P+b)*BesselI_scaled(0,n*M,psi)+(3*(P+b)*(P+b)+.75*M*M)*M*BesselI_scaled(1,n*M,psi)+.25*M*M*M*BesselI_scaled(3,n*M,psi)+1.5*(P+b)*M*M*BesselI_scaled(2,n*M,psi))*exp(-n*P+n*M);
    if (j==2) return (((P+b)*(P+b)+M*M/2)*BesselI_scaled(0,n*M,psi)+2*(P+b)*M*BesselI_scaled(1,n*M,psi)+M*M/2*BesselI_scaled(2,n*M,psi))*exp(-n*P+n*M);
    if (j==1) return ((P+b)*BesselI_scaled(0,n*M,psi)+M*BesselI_scaled(1,n*M,psi))*exp(-n*P+n*M);
    if (j==0) return BesselI_scaled(0,n*M,psi)*exp(-n*P+n*M);
    //if (j==3) return (((P+b)*(P+b)+1.5*M*M)*(P+b)*BesselI_scaled(0,n*M,psi)+(3*(P+b)*(P+b)+.75*M*M)*M*BesselI_scaled(1,n*M,psi)+.25*M*M*M*BesselI_scaled(3,n*M,psi)+1.5*(P+b)*M*M*BesselI_scaled(2,n*M,psi))*exp(-n*P+n*M)/M_PI;
    //if (j==2) return (((P+b)*(P+b)+M*M/2)*BesselI_scaled(0,n*M,psi)+2*(P+b)*M*BesselI_scaled(1,n*M,psi)+M*M/2*BesselI_scaled(2,n*M,psi))*exp(-n*P+n*M)/M_PI;
    //if (j==1) return ((P+b)*BesselI_scaled(0,n*M,psi)+M*BesselI_scaled(1,n*M,psi))*exp(-n*P+n*M)/M_PI;
    //if (j==0) return BesselI_scaled(0,n*M,psi)*exp(-n*P+n*M)/M_PI;

    gsl_function F;
    struct expK_params params = {P, M, j, b, n};

    F.function = &expK_f;
    F.params = &params;
    double result = 0;
    double abserr = 0;
    double epsabs = abs_res_dk;
    double epsrel = rel_res_dk;
    double limit = 100000;
    gsl_set_error_handler_off();
    if(!w)
    {
        std::cout << "set the workspace in expkint  "<<std::endl;
        exit(EXIT_FAILURE);
        
    } 

    gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    if (result!=result||abs(result)==inf) 
    {
        std::cout << "expKint failure1  "<< result<<std::endl;
        exit(EXIT_FAILURE);
    }
    // std::cout<<result<<std::endl;
    return result;

    double sum=0;

    // std::cout << "Bessels" << j<< std::endl;
    // if (j==0||j==-1||j==1) 
    if (j==-1||j==-2) 
    {
        double al=(P+b)/sqrt(pow(P+b,2)-M*M);
        for (int r = 1; abs(inc)*100000000/ep >= abs(sum); r++)
        {
            inc=0;
            // std::cout << "Besselsr" << r<<"  "<<M*n<< std::endl;
            // if (j==2&&r!=1&&r!=2) inc =2*gsl_sf_bessel_In_scaled(r,M*n)*pow(-1,r)*((P+b)*(P+b)*sin(r*psi)/r+2*M*(P+b)/(r*r-1)*(r*cos(psi)*sin(r*psi)-cos(r*psi)*sin(psi))+M*M/2/r/(r*r-4)*(sin(2*psi)*cos(r*psi)*(-2*r)+sin(r*psi)*(cos(2*psi)*r*r+r*r-4)));
            // if (j==1&&r!=1) inc= 2*gsl_sf_bessel_In_scaled(r,M*n)*pow(-1,r)*((-sin(psi)*cos(r*psi)*r*M+sin(r*psi)*((P+b)*(r*r-1)+M*r*r*cos(psi)))/(r*(r*r-1)));
            // std::cout << "Bessel" << r<<"  "<<M*n << "   "<<inc<<"  "<<sum<<"  "<<psi<<std::endl;
            // if (j==0) inc=2*gsl_sf_bessel_In_scaled(r,M*n)*pow(-1,r)*sin(r*psi)/r;
            // double inc1=0;
            // double sum1=0;
            if (j==-1) inc=pow((al-1)/(1+al),.5*r)*pow(-1,r)*2*BesselI_scaled(r,M*n,psi);
            if (j==-2) inc=pow((al-1)/(1+al),.5*r)*pow(-1,r)*2*BesselI_scaled(r,M*n,psi)*(al+r);
            // inc+=gsl_sf_bessel_In_scaled(r,M*n)*pow((al-1)/(al+1),.5*r)*(sin(2*r*psi)+2*r*(psi-M_PI))/r;
            

            // for (int l = 1; abs(inc1) >= ep / 10000000000000000 * abs(sum1); l++)
            // {
            //     if (psi==0) l=r;
            //     if (j==-1) inc1=4*gsl_sf_bessel_In_scaled(r,M*n)*pow(-1,r+l)*pow((al-1)/(al+1),.5*l)*coscos(psi,r,l);
            //     if (j==-2) inc1=4*gsl_sf_bessel_In_scaled(r,M*n)*pow(-1,r+l)*(al+l)*pow((al-1)/(al+1),.5*l)*coscos(psi,r,l);
            //     sum1 += inc1;
            //     if (l == 5 && inc1 == 0 && sum1 == 0)break;
            //     if (psi==0) break;
            // }
            // inc+=sum1;
            sum+=inc;
            if (r == 5 && inc == 0 && sum == 0)break;
            //return -gsl_sf_bessel_I0_scaled(M*n)*exp(-n*P+n*M);

        }
        //  std::cout << "Besself" << j<<"  "<<M*n << std::endl;
        // double I0=gsl_sf_bessel_I0_scaled(M*n);
        // double I1=gsl_sf_bessel_I1_scaled(M*n);
        // double I2=gsl_sf_bessel_In_scaled(2,M*n);


        // if (j==2) return (sum+I0*????????????)*exp(-n*P+n*M)/M_PI;
        
        if (j==-1)return (sum+BesselI_scaled(0,M*n,psi))*exp(-n*P+n*M)/sqrt(pow(P+b,2)-M*M);
        if (j==-2)return (sum+al*BesselI_scaled(0,M*n,psi))*exp(-n*P+n*M)/(pow(P+b,2)-M*M);
        //if (j==-1)return (sum+BesselI_scaled(0,M*n,psi))*exp(-n*P+n*M)/M_PI/sqrt(pow(P+b,2)-M*M);
        //if (j==-2)return (sum+al*BesselI_scaled(0,M*n,psi))*exp(-n*P+n*M)/M_PI/(pow(P+b,2)-M*M);

    
    }
    //double l = 0;

std::cout << "reimpliment expKint" << j<<"  "<<b<< std::endl;
exit(EXIT_FAILURE);

    

// 







    // //std::cout << "expkint next step" << std::endl;
    
    // // double prevlag=inf;
    // for (int r = 0; abs(inc)*10000/ep  >= abs(sum); r++)
    // {
    //     inc=0;
    //     // std::cout<<l<<"  "<<inc<<std::endl;
    //     // double hyper=0;
    //     // if (cos(psi)*cos(psi)==1) hyper=gsl_sf_gamma((3+cs+l)/2)*sqrt(M_PI)/gsl_sf_gamma((2+cs+l)/2);
    //     // if (cos(psi)*cos(psi)<1) hyper=gsl_sf_hyperg_2F1(.5,(1+cs+l)/2,(3+cs+l)/2,cos(psi)*cos(psi));
    //     // if (prevlag!=inf&&((j==0)||j==1||j==2))
    //     // {
    //     //     double x=-n*(P+M+b);
    //     //     // if (j==-1) prevlag=prevlag*expGamma(l+1,x)/expGamma(l,x)/(-l);
    //     //     if (j==0) prevlag=prevlag*(-x/l);
    //     //     if (j==1) prevlag=prevlag*(x-l)/(x-l+1)*(-x/l);
    //     //     if (j==2) prevlag=prevlag*((x-l)*(x-l)-l)/((x-l+1)*(x-l+1)-l+1)*(-x/l);
    //     // }
    //     // else prevlag=lag(l,j-l,n*(P+M+b),-n*M*(1+cos(psi)));

    //     // inc = -pow(2,cs)*pow(-2*M/(P+M+b), l) *Beta(cos(psi/2.)*cos(psi/2.),.5+cs+l,.5)*pow(P+M+b,j)*prevlag;

    //     // std::cout << "expKintinc " << inc << " "<<j<<"  "<<l<<"    "<<pow(2,cs)<<"   "<<xofact(2*M*n, l)<<"    "<<exp(-n*M*(1+cos(psi)))<<" step "<<"    "<<pow(-2*M/(P+M+b), l)<<"   "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P-M))<<std::endl;
        
    //     ////////////////////
    //     // l=int(M*n*(1-cos(psi)))+r;
    //     // // inc = pow(2,cs+1)*xofact(2*M*n, l)*exp(-n*M*(1+cos(psi)))*sin(psi/2)*pow(P+b-M,j)*F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P-M));
    //     // if (psi!=0) inc += 1./(l+.5)*xofact(M*n*(1-cos(psi)),l,-n*(M+P))*sin(psi/2)*pow(P+b+M,j)*F1(.5+l,.5,-j,1.5+l,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M));
    //     // l=int(M*n*(1-cos(psi)))-1-r;
    //     // if (l>=0&&psi!=0) inc+=1./(l+.5)*xofact(M*n*(1-cos(psi)),l,-n*(M+P))*sin(psi/2)*pow(P+b+M,j)*F1(.5+l,.5,-j,1.5+l,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M));

    //     // l=int(2*M*n)+r;
    //     // inc-=1./(l+.5)*xofact(2*M*n,l,-n*(M+P))*pow(P+b+M,j)*F1(.5+l,.5,-j,1.5+l,1,2*M/(P+b+M));

    //     // l=int(2*M*n)-1-r;
    //     // if (l>=0) inc-=1./(l+.5)*xofact(2*M*n,l,-n*(M+P))*pow(P+b+M,j)*F1(.5+l,.5,-j,1.5+l,1,2*M/(P+b+M));
    //     /////////////////////

    //     // std::cout << "resum this around the center of the poisson" << std::endl;
    //     // exit(EXIT_FAILURE);
    //     // std::cout<<inc<<std::endl;
    //     // if (F1(.5+l,.5,-j,1.5+l,1,2*M/(P+b+M))<0)
    //     // {
    //     //     std::cout << "F1 neg failure" <<F1(.5+l,.5,-j,1.5+l,1,2*M/(P+b+M))<<"   "<<.5+l<<"  "<<.5<<"  "<<-j<<" "<<1.5+l<<"  "<<1<<"  "<<2*M/(P+b+M)   << std::endl;
    //     //     exit(EXIT_FAILURE);
    //     // }
        
    //     // inc = 2*xofact(-2*M*n,l,n*M*(cos(psi)+1))*sin(psi/2)*pow(P+b+M,j)*F1(.5,.5-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M));
    //     // inc-=2*xofact(-2*M*n,l)*pow(P+b+M,j)*F1(.5,.5-l,-j,1.5,1,M*2/(P+b+M));
    //     //std::cout << "F1test " << .5+l << " "<<.5<<"  "<<-j<<"    "<<1.5+l<<"   "<<1<<"    "<<2*M/(P+b+M)<<"   "<<F1(.5+l,.5,-j,1.5+l,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M))<<"   "<<std::endl;
    //           // inc = -pow(-M*cos(psi)/(P+b),l)*pow(cos(psi),1+cs)*explag(l,j-l,n*(P+b))*pow(P+b,j)*hyper/(1+cs+l)*exp(n*b);
    //           // exp(n*b)/fact(l)*pow(-n,l)*pow(P+b,j+l)*powcossum(P,M,psi,j+l,b,cs,sn);
    //     if (inc == inf||inc!=inc||(inc==0&&psi!=M_PI&&l<1000))
    //     {
    //         std::cout << "expKintinc " << inc << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<M*n*(1-cos(psi))<<"  "<<xofact(2*M*n,l,-n*(M+P))<<"    "<<1./(l+.5)*sin(psi/2)*pow(P+b+M,j)<<" step "<<"    "<<F1(.5+l,.5,-j,1.5+l,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M))<<"   "<<std::endl;
    //         // std::cout << "expKintinc " << inc << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<xofact(2*M*n, l)<<"   "<<exp(-2*n*M)<<"    "<<pow(P+b+M,j)<<" step "<<"    "<<F1(.5+l,.5,-j,1.5+l,1,2*M/(P+b+M))<<"   "<<std::endl;
    //         std::cout << "expKint failure" << std::endl;
    //         exit(EXIT_FAILURE);
    //     }
    //     sum += inc;
        
    //     if (l == 5 && inc == 0 && sum == 0)
    //         break;
    // }
    //  if (sum == inf||sum!=sum)
    //     {
    //         std::cout << "expKintinc " << inc << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<M*n*(1-cos(psi))<<"  "<<xofact(2*M*n,l,-n*(M+P))<<"    "<<1./(l+.5)*sin(psi/2)*pow(P+b+M,j)<<" step "<<"    "<<F1(.5+l,.5,-j,1.5+l,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M))<<"   "<<std::endl;
    //         // std::cout << "expKintinc " << inc << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<xofact(2*M*n, l)<<"   "<<exp(-2*n*M)<<"    "<<pow(P+b+M,j)<<" step "<<"    "<<F1(.5+l,.5,-j,1.5+l,1,2*M/(P+b+M))<<"   "<<std::endl;
    //         std::cout << "expKintsum failure" << std::endl;
    //         exit(EXIT_FAILURE);
    //     }
    // // std::cout << "finishing expKint " << sum/M_PI<<"  "<<std::endl;
    // return sum/M_PI;





    // double sum = 0;
    // double inc = 1000;
    // double pinc = 1000;
    // for (int l = 0; abs(pinc) >= ep / 20 * abs(sum) || (abs(inc) > abs(pinc)); l++)
    // {
    //     pinc = inc;
    //     double hyper=0;

    //     if (cos(psi)*cos(psi)==1) hyper=gsl_sf_gamma((3+cs+l)/2)*sqrt(M_PI)/gsl_sf_gamma((2+cs+l)/2);
    //     if (cos(psi)*cos(psi)<1) hyper=gsl_sf_hyperg_2F1(.5,(1+cs+l)/2,(3+cs+l)/2,cos(psi)*cos(psi));

    //     inc = pow(-1,b+1)*pow(R,-b-l)*gsl_sf_choose(-b,l)*gsl_sf_hyperg_1F1(-l,-b-l+1,-R*nM)*pow(cos(psi),1+cs+l)*hyper/(1+cs+l);
    //     sum += inc;
    //     if (l == 5 && inc == 0 && sum == 0)
    //         break;
    // }
    // return sum;

    // if (j==-1)
    // {
    //     double sum = 0;
    //     double inc = 1000;
    //     double pinc = 1000;
    //     for (int l = 0; abs(pinc) >= ep / 20 * abs(sum) || (abs(inc) > abs(pinc)); l++)
    //     {
    //         pinc = inc;
    //         double hyper=0;
    //         if (cos(psi)*cos(psi)==1) hyper=sqrt(M_PI)*gsl_sf_gamma(.5*(l+3))/gsl_sf_gamma(.5*(l+2));
    //         if (cos(psi)*cos(psi)<1) hyper=gsl_sf_hyperg_2F1(.5,.5*(l+1),.5*(l+3),cos(psi)*cos(psi));

    //         inc = -1/M_PI/nP*pow(nM/nP,l)*gsl_sf_gamma_inc_Q(1+l,nP)*pow(cos(psi),l+1)*hyper*exp(-(nP-nM*cos(psi)));
    //         sum += inc;
    //         if (l == 5 && inc == 0 && sum == 0)
    //             break;
    //     }
    //     return sum;
    // }
    // if (j>=0)
    // {
    //     double sum = 0;
    //     double inc = 1000;
    //     double pinc = 1000;
    //     for (int l = 0; abs(pinc) >= ep / 20 * abs(sum) || (abs(inc) > abs(pinc)); l++)
    //     {
    //         pinc = inc;
    //         double hyper=0;
    //         if (cos(psi)*cos(psi)==1) hyper=sqrt(M_PI)*gsl_sf_gamma(.5*(l+j+3))/gsl_sf_gamma(.5*(l+j+2));
    //         if (cos(psi)*cos(psi)<1) hyper=gsl_sf_hyperg_2F1(.5,.5*(l+j+1),.5*(l+j+3),cos(psi)*cos(psi));
    //         inc = -1/M_PI*pow(nM,l)/fact(l)*exp(-nP)*pow(cos(psi),l+j+1)*hyper*exp(-(nP-nM*cos(psi)))/(l+j+1);
    //         sum += inc;
    //         if (l == 5 && inc == 0 && sum == 0)
    //             break;
    //     }
    //     return sum;

    // }
    

    std::cout<<"expKint not complete"<<std::endl;
    exit(EXIT_FAILURE);
}

double BesselItwid(double j, double A, double B,double a=0)
{
    double inc=1000;
    long double sum = 0;
    for (int k = 0; abs(inc) * 100000000000./ep >= abs(sum); k++)
    {
        inc=GammaQOn(2*k+j,A,B/2/A,a+gsl_sf_lnchoose(2*k+j,k));
        sum += inc;
        // 
        //std::cout << inc<<"  "<<sum<<"  "<<2*k+j<<"  "<<A <<"  "<<exp((2*k+j)*log(B/2/A)+gsl_sf_lnchoose(2*k+j,k))<<"  "<<gsl_sf_gamma_inc_Q(2*k+j,A)<< std::endl;
        // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"  "<<a<<"  "<<psi<<"  "<<psis<<"  "<<(P-M+a)<<"  "<<sum1<<std::endl;//F1(.5,.5-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M))<<std::endl;//"    "<<xofact(M*n,l)<<"   "<<expEn(i-l,n*(P+M+a))<<"    "<<powcossum(P,M,psi,j,b,cs+l,sn)<<"   "<<exp(-n*(a))<<"   "<<exp(-n*(P+M))<<std::endl;

        if (k == 5 && inc == 0 && sum == 0) break;

    }
    if (abs(sum) == inf || sum != sum)
        {
            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            // std::cout << "expintKincfail2 " << inc << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<sum1<<"  "<<GammaQOn(l,n*(P+M+a),(1-cos(psis))*M/(P+M+a),n*a)<<" |||| "<<(1-cos(psis))*M/(P+M+a)<<"    "<<n*(P+M+a)<<" step "<<n*a<<"    "<<"   "<<psim<<std::endl;
        
            
            
            std::cout << "BesselItwid failure "<< sum<< "   "<<j<<std::endl;
            exit(EXIT_FAILURE);
        }
    // std::cout <<" done  "<<sum << std::endl;
    return sum;
}
double ExpEncos(double A,double B,double r,double psi,double a=0)
{
    double inc=1000;
    long double sum = 0;
    if (psi==0) return -BesselItwid(r,A,B,a)*pow(-1,r);
    for (int k = 1; abs(inc)*100000000000./ep  >= abs(sum)||k<=r; k++)
    {
        inc=2*BesselItwid(k,A,B,a)*pow(-1,k)*coscos(psi,r,k);
        sum += inc;
        // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"  "<<a<<"  "<<psi<<"  "<<psis<<"  "<<(P-M+a)<<"  "<<sum1<<std::endl;//F1(.5,.5-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M))<<std::endl;//"    "<<xofact(M*n,l)<<"   "<<expEn(i-l,n*(P+M+a))<<"    "<<powcossum(P,M,psi,j,b,cs+l,sn)<<"   "<<exp(-n*(a))<<"   "<<exp(-n*(P+M))<<std::endl;

        if (k == 5 && inc == 0 && sum == 0) break;

    }

    if (abs(sum) == inf || sum != sum)
        {
            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            // std::cout << "expintKincfail2 " << inc << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<sum1<<"  "<<GammaQOn(l,n*(P+M+a),(1-cos(psis))*M/(P+M+a),n*a)<<" |||| "<<(1-cos(psis))*M/(P+M+a)<<"    "<<n*(P+M+a)<<" step "<<n*a<<"    "<<"   "<<psim<<std::endl;
        
            
            
            std::cout << "expencos failure" << std::endl;
            exit(EXIT_FAILURE);
        }
        // std::cout<<(sum+BesselItwid(0,A,B)*coscos(psi,r,0))/M_PI<<std::endl;
    return (sum+BesselItwid(0,A,B,a)*coscos(psi,r,0))/M_PI;

}

struct expint_params { double P;double M; double i;double j;double b; double a;double n;};

double expint_f (double x, void * p)
{
    struct expint_params * params = (struct expint_params *)p;
    double P = (params->P);
    double M = (params->M);
    double n = (params->n);
    double i = (params->i);
    double j = (params->j);
    double a = (params->a);
    double b = (params->b);
    double K=P+M*cos(x);

    return  expEn(i,n*(K+a))*pow(K+b,j)*exp(-n*(K))/M_PI;
}

double expintK(double P,double M, double psi,double i,double j,double b, double a,double n, gsl_integration_workspace *w )
{
    gsl_function F;
    struct expint_params params = {P, M,i, j, b, a, n};

    F.function = &expint_f;
    F.params = &params;
    double result = 0;
    double abserr = 0;
    double epsabs = abs_res_dk;
    double epsrel = rel_res_dk;
    double limit = 100000;
    gsl_set_error_handler_off();
    if(!w)
    {
        std::cout << "set the workspace in expintk  "<<std::endl;
        exit(EXIT_FAILURE);
    } 
    // std::cout<<result<<std::endl;
    gsl_integration_qags(&F, M_PI, psi, epsabs, epsrel, limit, w, &result, &abserr);
    if (result!=result||abs(result)==inf) 
    {
        std::cout << "expintK failure1  "<< result<<"  "<<P<<"  "<<M<<"   "<<i<<"   "<<j<<"   "<<b<<"     "<<a<<"    "<<n<<std::endl;
        exit(EXIT_FAILURE);
    }
    return result;

    // if (b<0)
    // {
    //     std::cout << "b<0" << std::endl;
    //     exit(EXIT_FAILURE);
    // }
    //if (psi==M_PI) return 0;
    if (fabs(psi-M_PI)<small_number) return 0;

    if (i!=1)
    {
        std::cout << "expintK i failure" << std::endl;
        exit(EXIT_FAILURE);
    }
    double inc=1000;
    //std::cout << "starting expintk " <<j<< "   "<<n*(P+a)<<"   "<<n*M<<"   "<<psi<<"   "<<n*a<<"   "<<std::endl;
    //////
    if (P+a>M)
    {
        if (j==0) 
        {
            // std::cout<<ExpEncos(n*(P+a),n*M,0,psi,n*a)<<"  "<<exp(n*a)<<std::endl;
            double mine=ExpEncos(n*(P+a),n*M,0,psi,n*a);
            //std::cout<<ExpEncos(n*(P+a),n*M,0,psi,n*a)<<" finished expintk j=0 "<<std::endl;
            return mine;
        }
        if (j==1) return (P+b)*ExpEncos(n*(P+a),n*M,0,psi,n*a)+M*ExpEncos(n*(P+a),n*M,1,psi,n*a);
        if (j==-1)
        {
            double al=(P+b)/sqrt(pow(P+b,2)-M*M);
            long double sum = 0;
            for (int k = 1; 100000000.*abs(inc)/ep >= abs(sum); k++)
            {
                inc=2*pow(-1,k)*ExpEncos(n*(P+a),n*M,k,psi,n*a)*pow((al-1)/(1+al),.5*k);
                sum += inc;
                 //std::cout << "expintkinc " << inc << " "<<k<<"  "<<sum<<std::endl;//F1(.5,.5-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M))<<std::endl;//"    "<<xofact(M*n,l)<<"   "<<expEn(i-l,n*(P+M+a))<<"    "<<powcossum(P,M,psi,j,b,cs+l,sn)<<"   "<<exp(-n*(a))<<"   "<<exp(-n*(P+M))<<std::endl;

                if (k == 5 && inc == 0 && sum == 0) break;

            }
            return (sum+ExpEncos(n*(P+a),n*M,0,psi,n*a))/sqrt((P+b)*(P+b)-M*M);
        }

        
        std::cout << "bad j "<< j << std::endl;
        exit(EXIT_FAILURE);
        
    }

    if (b<0&&j==1) return expintK(P,M,psi,i,j,0, a,n,w)+b*expintK(P,M,psi,i,0,0, a,n,w);


    

    
    
    //////

    ;

    //////
    
    // long double mysum=0;
    //     double mypsi=M_PI;
    //     double mystep=(psi-M_PI)/psi_res;
    //     double K=P+M*cos(mypsi);
    //     for (int l = 0; l < psi_res; l++)
    //     {
    //         K=P+M*cos(mypsi);
    //         inc=expEn(i,n*(K+a))*pow(K+b,j)*exp(-n*(K))*mystep;
    //         mysum+=inc;
    //         mypsi+=mystep;
    //         // std::cout<<inc<<std::endl;
    //         if (inc == inf || inc != inc)
    //             {
    //                 // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //                 // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
    //                 // std::cout << "expKintincfailr "<<std::endl;
    //                 std::cout << "expintK failure" << std::endl;
    //                 exit(EXIT_FAILURE);
    //             }
    //     }
    //     //std::cout << mystep << std::endl;
    //     //std::cout << inc << std::endl;
    //     return mysum/M_PI;

    //////
    

    double psis=M_PI;
    double psim=0;
//int E_i(K+a)*(K+b)^j
    //if (psi==0) return 0;

    long double sum = 0;
    double sum1=0;
    
    for (int l = 0; abs(inc)*10000./ep >= abs(sum); l++)
    {
        inc=0;
        // double hyper = 0;
        // inc = xofact(M*n,l)*expEn(i-l,n*(P+M+a))*pow(P+M+b,j)*powcossum(P,M,psi,j,b,cs+l,sn)*exp(-n*(M*(1+cos(psi))+a-b));//*exp(n*a-n*b);//*pow(M,l-j);
        //inc = pow(-(P+b-M)/(a-b),l)*pow(P+b-M,j)*exp(n*(P-M*cos(psi)-a))*choose(l-i,l)*fact(-i)*GammaQ(l-i+1,n*(a-b))*pow(n*(a-b),i-1)*pow(2,cs+1)*sin(psi/2)*F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M));//*exp(n*a-n*b);//*pow(M,l-j);
        // inc = pow(2*M/(P+M+a),l)*exp(n*(P-M*cos(psi)+a))*choose(l-i,l)*fact(-i)*pow(n*(P+M+a),i-1)*pow(P+b-M,j)*GammaQ(1+l-i,n*(P+M+a))*pow(2,cs+1)*sin(psi/2)*F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M));//*exp(n*a-n*b);//*pow(M,l-j);
        // inc = xofact(2*M/(P+M+a),l)*exp(n*(-M*cos(psi)-M))*pow(n*(P+M+a),i-1)*pow(P+b-M,j)*expGamma(1+l-i,n*(P+M+a))*2*sin(psi/2.)*F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M));//*exp(n*a-n*b);//*pow(M,l-j);
        
        
        
        // ???????????????inc = xofact(2*M/(P+M+a),l)*exp(n*(-M*cos(psi)-M))*pow(n*(P+M+a),i-1)*pow(P+b-M,j)*expGamma(1+l-i,n*(P+M+a))*2*sin(psi/2.)*F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M));//*exp(n*a-n*b);//*pow(M,l-j);
        // inc = pow(2,cs+1)*xofact(2*M*n, l)*exp(-n*M*(1+cos(psi)))*sin(psi/2)*pow(P+b-M,j)*F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P-M));
        if (P+a<M)
        {
            // std::cout << "start"<<std::endl;
            psis=acos(-(a+P)/M);
            // psis=acos(-(a+P)/M);
            // //std::cout<<-(a+P)/M<<std::endl;
            // if (P+M*cos(psi)+a>=0)
            // {
            //     //std::cout << M*(1-cos(psi))/(P+b+M) << std::endl;
            //     // std::cout << "help" << std::endl;
            //     if (psi!=0) inc+= 1./(l+.5)*pow((1-cos(psi))*M/(P+M+a),l)*sin(psi/2)*F1(.5+l,.5,-j,1.5+l,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M))*GammaQOn(l,n*(P+M+a))*exp(n*a);
                
            //     inc-=1./(l+.5)*sin(psis/2)*F1(.5+l,.5,-j,1.5+l,sin(psis/2)*sin(psis/2),M*(1-cos(psis))/(P+b+M))*GammaQOn(l,n*(P+M+a))*exp(n*a);
                
            // }
            psim=std::max(psi,psis);
            sum1 = 0;
            //std::cout << P+M*cos(psis)<<std::endl;
            double inc1 = 0;

            /////
            // if (l==0)
            // {
            // long double mysum=0;
            // double mypsi=M_PI;
            // double mystep=(M_PI-psim)/psi_res;
            // double K=P+M*cos(mypsi);
            // for (int l = 0; l < psi_res; l++)
            // {
            //     K=P+M*cos(mypsi);
            //     inc1=-expEn(i,n*(K+a))*pow(K+b,j)*exp(-n*(K))*mystep;
            //     mysum+=inc1;
            //     mypsi-=mystep;
            //     // std::cout<<inc<<std::endl;
            //     if (inc == inf || inc != inc)
            //         {
            //             // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            //             // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            //             // std::cout << "expKintincfailr "<<std::endl;
            //             std::cout << "expintK failure" << std::endl;
            //             exit(EXIT_FAILURE);
            //         }
            // }
            // //std::cout << mystep << std::endl;
            // //std::cout << inc << std::endl;
            // inc+=mysum/pow(P+b+M,j);
            // }
            ///

            ////////////////////////////////////
            double EM=0.5772156649;
            //  std::cout << "start "<<2*M/(P+M+b)*sin(psim/2)*sin(psim/2)<<"  "<<2*M/(P+M+b)<<std::endl;
                if (l==0) inc+=-(EM+log(n*(M-a-P)))*exp(n*a)*2*(F1(.5,.5,-j,1.5,sin(psim/2)*sin(psim/2),2*M/(P+M+b)*sin(psim/2)*sin(psim/2))*sin(psim/2)-F1(.5,.5,-j,1.5,1,2*M/(P+M+b)));
            //  std::cout << "finish"<<std::endl;
             

            if (l>0)
            {
                //   inc+=pow(2*M/(M-P-a),l)*2/l*(pow(cos(psim/2)*cos(psim/2),.5+l)*pow((P+b-M)/(P+M*cos(psim)+b)*sin(psim/2)*sin(psim/2),j)*F1(1,1+l+j,-j,1.5,sin(psim/2)*sin(psim/2),(P+b-M)/(P+M*cos(psim)+b)*sin(psim/2)*sin(psim/2))*sin(psim/2)-F1(.5,.5-l,-j,1.5,1,2*M/(P+M+b)))*exp(n*a);
                // std::cout << inc<<"   " << std::endl;
                for (int r = 0; 1000/ep*abs(inc1) >=abs(sum1); r++)
                {
                    // std::cout <<r<<"  "<<l<<std::endl;
                    inc1=0;
                    // std::cout<<2*M/(P+M+b)*sin(psim/2)*sin(psim/2)<<std::endl;
                    //  inc1=-choose(j,r)*pow(-2*M/(P+M+b),r)*F21(.5+l,.5-r,1.5+l,cos(psim/2)*cos(psim/2));
                    // if (2*M/(P+M+b)*sin(psim/2)*sin(psim/2)>1) std::cout<<2*M/(P+M+b)<<std::endl;
                    //  if (r<=l) inc1-=pow(n,l)*xofact((M-P-a),l-r)/l*xofact(-2*M,r,n*a)*2*(sin(psim/2)*F1(.5,.5-r,-j,1.5,sin(psim/2)*sin(psim/2),2*M/(P+M+b)*sin(psim/2)*sin(psim/2))-F1(.5,.5-r,-j,1.5,1,2*M/(P+M+b)));//issues?????????/
                    int w=0;
                    double mine=0;
                    /*DANI
		    for (w = 0; w<=l; w++)///get rid of this one
                    {
                        mine=pow(-1,w)*exp((l-w)*log(-a-P+M)+w*log(2*M*cos(psim/2)*cos(psim/2))-gsl_sf_lnfact(w)-gsl_sf_lnfact(l-w)+n*a+l*log(n));
                        inc1-=-1./(.5+w)*mine*F21(.5+w,.5-r,1.5+w,cos(psim/2)*cos(psim/2));//Beta(cos(psim/2)*cos(psim/2),.5+w,.5+r);
                        if (inc1 == inf || inc1 != inc1)
                        {
                            //Beta(cos(psim/2)*cos(psim/2),.5+l,.5+r)
                            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
                            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
                            std::cout << "expintKincfailr1 " << inc1 << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<sum1<<"  "<<w<<"    "<<mine<<"  "<<(l-w)*log(-a-P+M)+w*log(2*M*cos(psim/2)*cos(psim/2))-gsl_sf_lnfact(w)-gsl_sf_lnfact(l-w)+n*a+l*log(n)<<" step "<<"    " <<"   "<<F21(.5+w,.5-r,1.5+w,cos(psim/2)*cos(psim/2))<<std::endl;
                            std::cout << "expintKr failure" << std::endl;
                            exit(EXIT_FAILURE);
                        }
                        //   std::cout << inc1<<"   " << w<<std::endl;
                    }
		    *///DANI

                    
                       
                       
                       
                     inc1+=-pow(2*M*cos(psim/2)*cos(psim/2)/(M-a-P)*exp(n*a/l),l)/(.5+l)*F21(.5+l,.5-r,l+1.5,cos(psim/2)*cos(psim/2));//Beta(cos(psi/2)*cos(psi/2),.5+l,.5+r)

                     if (inc1 == inf || inc1 != inc1)
                    {
                        //Beta(cos(psim/2)*cos(psim/2),.5+l,.5+r)
                        // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
                        // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
                        std::cout << "expintKincfailr2 " << inc1 << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<sum1<<"  "<<pow(M*(1+cos(psim))/(P+b-M),r)<<"    "<<xofact(-a-P+M,l)<<"  "<<xofact(-M*(1+cos(psim)),l)<<" step "<<"    " <<"   "<<F21(.5,.5+r+l,1.5+r+l,cos(psim/2)*cos(psim/2))<<std::endl;
                        std::cout << "expintKr failure" << std::endl;
                        exit(EXIT_FAILURE);
                    }

                    inc1*=choose(j,r)/l*pow(-2*M/(P+b+M),r)*cos(psim/2);

                    sum1+=inc1;

                    //change to 1+cos
                      //std::cout << sum1<<"   " <<l << std::endl;
                    // if (r>1000) std::cout << "expintkincsubr " << inc1 <<" "<< pinc1<<" "<<ep / 20 *sum1<<"  "<<r<<"    "<< l<<"   "<<pow(-2*M/(P+M+b),r)<<"    "<<Beta(cos(psim/2)*cos(psim/2),.5+l,.5+r)<<std::endl;
                    

                    if (r == 5 && inc1 == 0 && sum1 == 0) break;
                    if (inc1 == inf || inc1 != inc1)
                    {
                        //Beta(cos(psim/2)*cos(psim/2),.5+l,.5+r)
                        // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
                        // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
                        std::cout << "expintKincfailr3 " << inc1 << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<sum1<<"  "<<pow(M*(1+cos(psim))/(P+b-M),r)<<"    "<<xofact(-a-P+M,l)<<"  "<<xofact(-M*(1+cos(psim)),l)<<" step "<<"    " <<"   "<<F21(.5,.5+r+l,1.5+r+l,cos(psim/2)*cos(psim/2))<<std::endl;
                        std::cout << "expintKr failure" << std::endl;
                        exit(EXIT_FAILURE);
                    }
                }
                 inc+=sum1;
            }
            ///////////////////////////////

            ///////////////////////////////////////////////////////
            // double pinc1 = 1000;
            // for (int r = 0; 10000000/ep*abs(pinc1) >=abs(sum1) || (abs(inc1) > abs(pinc1)); r++)
            // {
            //     pinc1= inc1;
            //     inc1=-choose(j,r)*pow(-2*M/(P+M+b),r)*F21(.5+l,.5-r,1.5+l,cos(psim/2)*cos(psim/2));
            //     sum1+=inc1;
            //     // if (r>1000) std::cout << "expintkincsubr " << inc1 <<" "<< pinc1<<" "<<ep / 20 *sum1<<"  "<<r<<"    "<< l<<"   "<<pow(-2*M/(P+M+b),r)<<"    "<<Beta(cos(psim/2)*cos(psim/2),.5+l,.5+r)<<std::endl;
                

            // if (r == 5 && inc1 == 0 && sum1 == 0) break;
            // if (inc1 == inf || inc1 != inc1)
            // {
            //     //Beta(cos(psim/2)*cos(psim/2),.5+l,.5+r)
            //     // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            //     // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            //     std::cout << "expintKincfailr " << inc << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<sum1<<"  "<<pow(-2*M/(P+M+b),r)<<"    "<<" step "<<"    "<<exp(n*(P+M*cos(psi))+a) <<"   "<<GammaQOn(l,n*(P+M+a))<<std::endl;
            //     std::cout << "expintKr failure" << std::endl;
            //     exit(EXIT_FAILURE);
            // }
            // }
            // inc+=sum1*GammaQOn(l,n*(P-M+a),-2*M*cos(psim/2)*cos(psim/2)/(P-M+a),n*a)/(l+.5)*cos(psim/2);

            //////////////////////////////////////////////////

                // std::cout << inc<<std::endl;
                // inc+=sum1;
                
            // }
            

            // 
            // 
            // if (P-M)
            // for (int r = 0; r <= l; r++) sum1+=-GammaQOn(r,n*(P-M+a))*choose(j,l-r)*pow(-(P-M+b)/(P-M+a),r);
            // inc=sum1*F21(.5,.5+l,1.5+l,cos(psim/2)*cos(psim/2))*pow(M*(1+cos(psim))/(P-M+b),l)*pow(P-M+b,j)*cos(psim/2)/(l+.5);

        if (inc == inf || inc != inc||sum/M_PI*exp(n*(a))==inf)
        {
            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            std::cout << "expintKincfail1 " << inc << " "<<j<<"  "<<l<<"   "<<psi<<" |||||   "<<sum1<<"  "<<GammaQOn(l,n*(P-M+a),-2*M*cos(psim/2)*cos(psim/2)/(P-M+a),n*a)<<"  "<<pow(-2*M/(P-M+a),l)<<"    "<<a<<" step "<<"    "<<"   "<<psim<<std::endl;
        
            
            
            std::cout << "expintK failure" << std::endl;
            exit(EXIT_FAILURE);
        }
    
            

            
        }


        if (psi<psis)
        {
            if (psi!=0) inc+= 1./(l+.5)*sin(psi/2)*F1(.5+l,.5,-j,1.5+l,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M))*GammaQOn(l,n*(P+M+a),(1-cos(psi))*M/(P+M+a),n*a);
            inc-= 1./(l+.5)*sin(psis/2)*F1(.5+l,.5,-j,1.5+l,sin(psis/2)*sin(psis/2),M*(1-cos(psis))/(P+b+M))*GammaQOn(l,n*(P+M+a),(1-cos(psis))*M/(P+M+a),n*a);
        }
        if (!((psi < psis)||(P+a<M)))
        {
            std::cout << "expintK ratio failure"<<psi<<"  "<<psis <<"  "<<P+a <<"  "<<M <<std::endl;
            exit(EXIT_FAILURE);
        }
        
        //  std::cout << inc<<"   " << sum<<std::endl;
        sum += inc;
        // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"  "<<a<<"  "<<psi<<"  "<<psis<<"  "<<(P-M+a)<<"  "<<sum1<<std::endl;//F1(.5,.5-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b+M))<<std::endl;//"    "<<xofact(M*n,l)<<"   "<<expEn(i-l,n*(P+M+a))<<"    "<<powcossum(P,M,psi,j,b,cs+l,sn)<<"   "<<exp(-n*(a))<<"   "<<exp(-n*(P+M))<<std::endl;

        if (l == 5 && inc == 0 && sum == 0)
            break;
        if (inc == inf || inc != inc||sum/M_PI==inf)
        {
            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<< xofact(2*M*n,l)<<"   "<<exp(-n*M*(1+cos(psi)))<<"    "<<expEn(i-l,n*(P+M+a))<<"    "<<F1(.5,.5-cs-l,-j,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5-cs,-j-l,1.5,sin(psi/2)*sin(psi/2),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            // std::cout << "expintkinc " << inc << " "<<j<<"  "<<l<<"    "<<2*M/(P+M+a)<<"  "<<M/(P+b-M)<<"   "<<exp(n*(P-M*cos(psi)+a)) <<"  "<<pow(2*M/(P+M+a),l)<<"   "<<choose(l-i,l)*fact(-i)<<"    "<<pow(n*(P+M+a),i-1)*pow(P+b-M,j)<<"    ||| "<<i<<"   "<<GammaQ(1+l-i,n*(P+M+a))<<"   "<<F1(.5,.5-l,-j,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<F1(.5,.5,-j-l,1.5,sin(psi/2.)*sin(psi/2.),M*(1-cos(psi))/(P+b-M))<<"   "<<psi<<"  break"<<std::endl;
            std::cout << "expintKincfail2 " << inc << " "<<j<<"  "<<l<<"   "<<psi<<"    "<<sum1<<"  "<<GammaQOn(l,n*(P+M+a),(1-cos(psis))*M/(P+M+a),n*a)<<" |||| "<<(1-cos(psis))*M/(P+M+a)<<"    "<<n*(P+M+a)<<" step "<<n*a<<"    "<<"   "<<psim<<std::endl;
        
            
            
            std::cout << "expintK failure" << std::endl;
            exit(EXIT_FAILURE);
        }
    }
    //GO BACK TO DOING THE P+M+b-M(1+cos(psi))
    // std::cout << "finish"<<std::endl;
    //std::cout << "finishing expintK " << sum/M_PI<< std::endl;

    return sum/M_PI*pow(P+b+M,j);
}
