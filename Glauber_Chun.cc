#include <fstream>
#include <iostream>
#include <vector>
#include <cmath>
#include <sstream>
#include <string>
#include <assert.h>
#include <stdlib.h>

#include "Random.h"

//#include "global.h"

using std::vector;
using namespace std;

vector<pair<double, double> > binpos; 

int read_nuclear_ipsat(int nhyd, string cent);
void gxy_ipsat(double &x, double &y, int n);

int read_nuclear_ipsat(int nhyd, string cent)
{
        ostringstream filename;
        //filename << "hydros_" << cent.c_str() << "/job-" << nhyd << "/NcollList.dat";
        filename << "./NcollList.dat";
	ifstream initial;
        initial.open(filename.str().c_str());
	assert(!initial.fail());
        if (initial.fail()) {
          cout << "Initial open fail " << filename.str() << endl;
          exit(1);
        }

	string s, sub;
	//First line: some crap
	getline(initial,s);
        //Rest of lines: x, y; 
	double x, y;
	binpos.clear();
	do {
	  getline(initial,s);
	  if (initial.eof()) break;
	  istringstream iss(s);
	  iss >> x >> y;
	  pair<double, double> temp = make_pair(x, y);
          binpos.push_back(temp);
	} while(true);
	
      return int(binpos.size());

}

void gxy_ipsat(double &x, double &y, int index) {

  x = binpos[index].first;
  y = binpos[index].second;

}
