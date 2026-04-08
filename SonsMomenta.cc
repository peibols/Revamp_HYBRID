#include <cmath>
#include <iostream>
#include <assert.h>
#include <algorithm>
#include <vector>

#include "vector_operators.h"

using std::vector;
using namespace std;

void quenched_sons(vector<double> p, vector<double> qp, vector<double> &d1, vector<double> &d2);
double normalise(vector<double> &p);
vector<double> vec_prod(vector<double> a, vector<double> b);
std::vector<double> vec_abs(std::vector<double> p);

void quenched_sons(vector<double> p, vector<double> qp, vector<double> &d1, vector<double> &d2)
{
	//Normalise 3-momentum
	double qmod=normalise(qp);
	double modmom=normalise(p);
	double modthis=normalise(d1);
	double modoson=normalise(d2);
	//Define transverse axis for rotation and normalise
	vector<double> axis = vec_prod(p,qp);
	double modaxis=normalise(axis);
	//Find angle in plane
	double angle=0.;
	if (p[0]*qp[0]+p[1]*qp[1]+p[2]*qp[2]>=1.) angle=0.;
	else angle=acos(p[0]*qp[0]+p[1]*qp[1]+p[2]*qp[2]);
	//if (angle!=0.) cout << " angle= " << angle << " modaxis= " << modaxis << endl;
	//Perform Rodrigues rotation
	double thisscal=axis[0]*d1[0]+axis[1]*d1[1]+axis[2]*d1[2];
	vector<double> use = d1*cos(angle) + vec_prod(axis,d1)*sin(angle) + axis*thisscal*(1.-cos(angle));
	double oscal=axis[0]*d2[0]+axis[1]*d2[1]+axis[2]*d2[2];
        vector<double> ouse = d2*cos(angle) + vec_prod(axis,d2)*sin(angle) + axis*oscal*(1.-cos(angle));
	//Update momenta
	double lamp=qmod/modmom;
	double lambda=qp[3]/p[3];
	for (unsigned int i = 0; i<3; i++)
	{
		d1[i] = use[i]*modthis*lamp;
		d2[i] = ouse[i]*modoson*lamp;
	}
	d1[3]*=lambda;
	d2[3]*=lambda;
}

double normalise(vector<double> &p)
{
	double norm=sqrt(pow(p[0],2.)+pow(p[1],2.)+pow(p[2],2.));
	if (norm==0.) return norm;
	p[0]/=norm;
	p[1]/=norm;
	p[2]/=norm;
	return norm;
}

vector<double> vec_prod(vector<double> a, vector<double> b)
{
	assert(a.size()==b.size());

	vector<double> rot;
	rot.reserve(a.size());

	rot.push_back(a[1]*b[2]-a[2]*b[1]);
	rot.push_back(a[2]*b[0]-a[0]*b[2]);
	rot.push_back(a[0]*b[1]-a[1]*b[0]);
	rot.push_back(0.);

	return rot;
}

vector<double> vec_abs(vector<double> p)
{
        for (unsigned int i=0; i<p.size(); i++)
        {
                if (p[i]<0.) p[i]*=-1.;
        }
        return p;
}
