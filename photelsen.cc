#include "Pythia8/Pythia.h"
#include "Pythia8/Event.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <assert.h>

using namespace std;
using namespace Pythia8;

//Args: N, seed, sqrts (in GeV)
int main(int argc, char** argv) {

assert(argc==7);

double trigger_pt=atof(argv[4]);
double trigger_eta=atof(argv[5]);
int trigger_id=atoi(argv[6]);

char outFile[100];
sprintf(outFile,"tree_%s.txt",argv[2]);
ofstream geneal;
geneal.open (outFile);

Pythia pythia;
Info info = pythia.info;
Event& event = pythia.event;

//#Events
int N=atoi(argv[1]);

//Seed
ostringstream seedstring;
double see=733+atoi(argv[2]);
seedstring << "Random:seed = " << see; 
pythia.readString(seedstring.str().c_str());

ostringstream pythiaset;
pythiaset << "setup_pythia.cmnd";
pythia.readFile(pythiaset.str()); 

int finals, counter, m1, m2, use, qg, mam, ini;
double ultra[10000][7];

//Initialize
pythia.init();

int count=0;
do {

        bool have_trig=0;

	if (!pythia.next()) continue;
	for (int i = 0; i < pythia.event.size(); ++i) {
	  if (pythia.event[i].isFinal() && pythia.event[i].id()==trigger_id) {
	    if (fabs(pythia.event[i].y())<=trigger_eta && pythia.event[i].pT()>trigger_pt)
	    have_trig=1;
	    break;
          }
	}

        if (!have_trig) continue;

	if (pythia.event.size()>950) {
		cout << " size= " << pythia.event.size();
	}
	finals = 0;
	counter = 0;
	for (unsigned e=0;e<10000; e++) {
		for (unsigned f=0;f<7; f++) {
			ultra[e][f]=0.;
		}
	}

	for (int i = 0; i < pythia.event.size(); ++i) {
	        if (pythia.event[i].status()==-23) {
		        int ide = pythia.event[i].id();
			double px=pythia.event[i].px();
			double py=pythia.event[i].py();
			double pz=pythia.event[i].pz();
			double pe=pythia.event[i].e();
			geneal << -1000 << " " << px << " " << py << " " << pz << " " << pe << " " << 0 << " " << 0 << " " << 0 << " " << 0 << " "  << 0 << " "  << 0 << " "  << ide << "\n";

		}
		if (pythia.event[i].isFinal()) {
			if (pythia.event[i].mother2() == 0) continue;
			if (pythia.event.size()>10000) continue;
			finals+=1;
			use = i;
			for (unsigned k = 0; k<30; k++){
				m1 = pythia.event[use].mother1();
				m2 = pythia.event[use].mother2();
				if (m1!=m2) break;
				use = m1;
			}
			ultra[i][0]=pythia.event[i].px();
			ultra[i][1]=pythia.event[i].py();
			ultra[i][2]=pythia.event[i].pz();
			ultra[i][3]=pythia.event[i].e();
			ultra[i][4]=m1;
		}
	}
	donantli:
	for (unsigned l = 0; l<10000; l++) {
		for (unsigned m=0; m<10000; m++) {
			if (ultra[m][4]==ultra[l][4] && m!=l && ultra[m][5]==0 && ultra[m][4]!=0) {
				counter +=1;
				use = int(ultra[m][4]);
				ultra[use][0]=ultra[m][0]+ultra[l][0];
				ultra[use][1]=ultra[m][1]+ultra[l][1];
				ultra[use][2]=ultra[m][2]+ultra[l][2];
				ultra[use][3]=ultra[m][3]+ultra[l][3];
				int early = use;
				for ( unsigned y = 0; y<30; y++) {
					m1 = pythia.event[use].mother1();
					m2 = pythia.event[use].mother2();
					if (m1!=m2) break;
					use = m1;
				}
				if (pythia.event[m1].status()!=-41) {
					ultra[early][4]=m1;
				}
				else ultra[early][4]=3;
				ultra[m][5]=1;
				ultra[l][5]=1;
				if (counter < finals) goto donantli;
			}
		}
	}
	ultra[3][0]=0.;
	ultra[3][1]=0.;
	ultra[3][2]=0.;
	ultra[3][3]=0.;
	int col, acol, ide;
	for (unsigned c=3; c<10000; c++) {
		if (pythia.event[ultra[c][4]].status()==-41) ultra[c][4]=3;
		if ((ultra[c][4]!=0 || c==3) && pythia.event[c].status()!=-41) {
			ultra[c][6]=abs(sqrt(abs(pow(ultra[c][3],2.)-pow(ultra[c][0],2.)-pow(ultra[c][1],2.)-pow(ultra[c][2],2.)-pow(pythia.event[c].m(),2.))));
			ide=abs(pythia.event[c].id());
			if (ide >=1 && ide <=6) qg=1;
			if (ide==21) qg=2;
			if (ide==22) {
				qg=3;
			}
			ini=0;
			use=c;
			for (unsigned w=0; w<50; w++) {
				if (pythia.event[use].status()==-23) {
					ini=use;
					ide=abs(pythia.event[use].id());
					if (ide >=1 && ide <=6) ini*=1.;
					if (ide==21) ini*=-1.;
					break;
				}
				mam=pythia.event[use].mother1();
				use=mam;
			}

			ide=pythia.event[c].id();
			col=pythia.event[c].col();
			acol=pythia.event[c].acol();
			geneal << c << " " << ultra[c][0] << " " << ultra[c][1] << " " << ultra[c][2] << " " << ultra[c][3] << " " << ultra[c][6] << " " << int(ultra[c][4]) << " " << qg << " " << ini << " " << col << " " << acol << " " << ide << "\n";
		}
	}
	for (int c=0; c<pythia.event.size(); c++) {
		if (pythia.event[c].status()==63) {
			ide=pythia.event[c].id();
			col=pythia.event[c].col();
			acol=pythia.event[c].acol();
			double px=pythia.event[c].px();
			double py=pythia.event[c].py();
			double pz=pythia.event[c].pz();
			double pe=pythia.event[c].e();
			double pm=pythia.event[c].m();
			double pq=sqrtpos(pe*pe-px*px-py*py-pz*pz-pm*pm);
			geneal << c << " " << px << " " << py << " " << pz << " " << pe << " " << pq << " " << 0 << " " << 0 << " " << 0 << " "  << col << " "  << acol << " "  << ide << "\n";
		}
	}
	double weight = info.weight();
	double cross = info.sigmaGen();
	geneal << "NEXT" << " 0. " << " 0. " << " 0. " << "0.00" << count << " 0. " << " " << weight << " " << cross << " 0. " << " 0. " << " 0. " << " 0. " << "\n" ;
	count+=1;
} while (count<N);

geneal << "END" << " " << std::setprecision(9) << info.sigmaGen() << " " << info.weightSum() << endl;
geneal.close();
pythia.stat();

return 0;

} //End Program





