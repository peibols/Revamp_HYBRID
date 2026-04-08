// RES HYDRO QUENCHING

#include <sstream>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <time.h>

using namespace std;
int ir;
double TA[4000], step, bmin, bmax;
double norm=1.;
double tau0, deltat, deltax, deltay;
int maxt, maxx, maxy, it, ix, iy;
double hydrot[200][101][101], hydroe[200][101][101], hydrox[200][101][101], hydroy[200][101][101];
double dt, dx, dy;
double rando(int*);
double gTAA(double, double, double);
double gWN(double, double, double);
double gT(double, double, double);
double gVy(double, double, double);
double gVx(double, double, double);
double gE(double, double, double);
double getGrid(double, double, double, double *, double *, double *, int *, int *, int *);
double gQ(double);
double fuckedN=0.;
double totalN=0.;

double alpha;
double kappa;
int qmethod; //1 Coll, 2 Radiative, 3 Strong
int tmethod;
double rpower;
double tquench=0.6;
//ARGUMENTS: File, Alpha, Kappa, N, Model, Tmethod, Cent, Rpower, Rand, sqrts
int main (int argc, char** argv)
{
	assert(argc==11);
	char outFile[200], inpFile[200];
	alpha=atof(argv[2]);
	kappa=atof(argv[3]);
	int N=atoi(argv[4]);
	qmethod=atoi(argv[5]);
	tmethod=atoi(argv[6]);
	std::string Cent=argv[7];
	rpower=atof(argv[8]);
	cout << " rpower= " << rpower << " frand= " << atoi(argv[9]) << endl;
	ir=atoi(argv[9]);
	sprintf(outFile,"./%s.out",argv[1]);
	sprintf(inpFile,"./tree_%s.txt",argv[1]);
	
	//Output
	ofstream check;
	check.open (outFile);
	//Input
	std::ifstream infile1(inpFile);
	assert(!infile1.fail());

	double P[1000][12], num, px, py, pz, E, q, mom, qge, ini, col, acol, ide;
	double x,y,z,b;
	double pexe, peye, peze;
	double gxy(double*, double*, double*);
	double gdE(double, double, double, double, double *, double *, double *, double *,  double *, double *);
	double getRes(double,double,double,double,double,double,double,double,double,double,double);
	int count=0;
	double ene=5.;

	ostringstream sqrts;
        sqrts << argv[10];

	ostringstream hydFile;
        hydFile << "./";
        if (Cent=="0-5") {          //0-5%
                bmin=0.;
                bmax=3.5;
                //if (sqrts.str()=="020") bmin=0., bmax=2.74;	//Doga values
                if (sqrts.str()=="020") bmin=0., bmax=3.3;	//Hirano values
        }
        if (Cent=="5-10") {          //5-10%
                bmin=3.5;
                bmax=4.94;
                //if (sqrts.str()=="020") bmin=2.74, bmax=4.;	//Doga values
		if (sqrts.str()=="020") bmin=3.3, bmax=4.7;	//Hirano values
        }
        if (Cent=="10-20") {          //10-20%
                bmin=4.94;
                bmax=6.98;
                if (sqrts.str()=="020") bmin=4.7, bmax=6.7;
        }
        if (Cent=="20-30") {          //20-30%
                bmin=6.98;
                bmax=8.55;
		if (sqrts.str()=="020") bmin=6.7, bmax=8.2;
        }
        if (Cent=="30-40") {          //30-40%
                bmin=8.55;
                bmax=9.88;
		if (sqrts.str()=="020") bmin=8.2, bmax=9.4;
        }
        if (Cent=="40-50") {          //40-50%
                bmin=9.88;
                bmax=11.04;
		if (sqrts.str()=="020") bmin=9.4, bmax=10.6;
        }
        if (Cent=="50-60") {          //50-60%
                bmin=11.04;
                bmax=12.09;
		if (sqrts.str()=="020") bmin=10.6, bmax=11.6;
        }
        if (Cent=="60-70") {          //60-70%
                bmin=12.09;
                bmax=13.05;
		if (sqrts.str()=="020") bmin=11.6, bmax=12.5;
        }
	hydFile << "hydroinfoPlaintxtHuichaoFormat.dat";
        std::ifstream infile3(hydFile.str().c_str());
        assert(!infile3.fail());

	ostringstream glauFile;
        glauFile << "./TAb2LL.dat";
        std::ifstream infile2(glauFile.str().c_str());
        assert(!infile2.fail());

        cout << " bmin= " << bmin << " bmax= " << bmax << endl;

	//Reading T(t)
	double b2;
	for (unsigned a=0; a<4000; a++)
	{
		infile2 >> b2 >> TA[a];
		if (a == 1) step = b2;
	}

	//Reading HydroTemp
	double enedat, tdat, vxdat, vydat;
	double tou, hor, ver;
	maxt=179;
	maxx=100;
	maxy=100;
	tau0=0.6;
	deltat=0.1;
	deltax=0.3;
	deltay=0.3;
	tou=0.6;
	while (infile3 >> hor >> ver >> tou >> enedat >> tdat >> vxdat >> vydat)//horizontal,vetical, time, energy density, temperature in fluid rest, vel in fluid restframe of plasma
	{

		it = int((tou+deltat/2.-tau0)/deltat);//discretizing the plasma
          	ix = int((hor+deltax*maxx/2.+deltax/2.)/deltax);//discretizing the plasma
           	iy = int((ver+deltay*maxy/2.+deltay/2.)/deltay);//discretizing the plasma

		hydrot[it][ix][iy]=tdat;
		hydroe[it][ix][iy]=enedat;
		hydrox[it][ix][iy]=vxdat;
		hydroy[it][ix][iy]=vydat;
	}
	int event=0;
	clock_t startClock = clock();


	double weight, cross;
	//Event Loop
	do
	{
		//Getting initial x, y and b
		gxy(&x,&y,&b);//impact parameter is b
		//cout << " X= " << x << " Y= " << y << " B= " << b << "\n";
		double xcre=x;//creation x
		double ycre=y;//creation y

		//Getting tree info
		int FinId[500];
		for (unsigned a=0; a<1000; a++) //set a bunch of stuff to 0
		{
			P[a][0]=0.;
			P[a][1]=0.;
			P[a][2]=0.;
			P[a][3]=0.;
			P[a][4]=0.;
			P[a][5]=0.;
			P[a][6]=0.;
			P[a][7]=0.;
			P[a][8]=0.;
			P[a][11]=0.;
		}

		std::string snum;
		//Store information from tree and write remnants
		while (infile1 >> snum >> px >> py >> pz >> E >> q >> mom >> qge >> ini >> col >> acol >> ide)
		{
			if (snum!="NEXT")
			{
				num=strtod(snum.c_str(), NULL);
				if (num==-1000) {
				        int initag;
					if (int(ide)==21) initag=-2000;
					else if (abs(int(ide))<=6) initag=-1000;
					else initag=-3000;
					check << px << " " << py << " " << pz << " " << E  << " " << px << " " << py << " " << pz << " " << E << " " << int(col) << " " << int(acol) << " " << int(initag) << "\n";
					continue;
				}
				P[int(num)][0]=px;
				P[int(num)][1]=py;
				P[int(num)][2]=pz;
				P[int(num)][3]=E;
				P[int(num)][4]=q;	//Virtuality
				P[int(num)][5]=mom;
				//cout<<mom<<endl;
				P[int(num)][6]=qge;	//Quark or Gluon
				P[int(num)][7]=ini;
				P[int(num)][9]=col;
				P[int(num)][10]=acol;
				P[int(num)][11]=ide;
				//Remnants
				if (int(q+0.9)==0 && E!=0. && mom==0.)
				{
					check << px << " " << py << " " << pz << " " << E  << " " << px << " " << py << " " << pz << " " << E << " " << int(col) << " " << int(acol) << " " << int(ide) << "\n";
				}
			}
			else
			{
				weight=mom;
				cross=qge;
				event+=1;
				break;
			}
		}

		int k=0;
		//Assign final status; as I understand this, it counts the final particles and puts them in finid
		for (unsigned a=4; a<1000; a++)
		{
			if (P[int(a)][5]!=0.) //if the mom is not empty
			{
				int canfin=1; //index this variable
				for (unsigned b=4; b<1000; b++)
				{
					if (int(P[int(b)][5])==int(a))
					{
						canfin=0;
						break;
					}
				}
				if (canfin==1) FinId[k]=int(a), k+=1; //k is number of final particles and FinId has all of the particles with out moms in the set of particles
			}
		}
//---------------------test round to get q in test (mine)-----------------------------------------
		int fam1[30];
		double lifetimes[1000][11];// 0 is the single res time, 1 is the creation time, 2 is the finishing time,
		//3 is the total resolve time, 4 is a store, 5 is the inital x pos, 6 is inital y pos,
		//7 is inital zpos, 8 is final pos x, 9 is final pos y, 10 is final z pos
		for (unsigned zb=0; zb<1000; zb++)
		{
			for (unsigned za=0; za<8; za++)
			{
				lifetimes[zb][za]=0.;
			}
		}

		//First pass propagation
		for (int cr=0; cr<k; cr++)
		{
			int ind=FinId[cr];
			//cout << " ind= " << ind << endl;
			int fr=0;
			double acum;
			fam1[fr]=ind;
			//cout << " InFinal= " << ind<< endl;
			for (unsigned a=0; a<50; a++)
			{
				int mom=int(P[ind][5]);
				if (P[mom][8]==1. && mom!=3)
				{
					x=lifetimes[mom][8];		// X
					y=lifetimes[mom][9];		// Y
					z=lifetimes[mom][10];		// Z
					acum=lifetimes[mom][2];		// t
					//cout << " ind= " << ind << " z= " << z  << endl;
					lifetimes[ind][5]=x;//myline
					lifetimes[ind][6]=y;//myline
					lifetimes[ind][7]=z;//myline
					break;
				}
				if (mom==3)
				{
					x=xcre;
					y=ycre;
					z=0.;
					acum=0.;
					lifetimes[ind][5]=x;//myline
					lifetimes[ind][6]=y;//myline
					lifetimes[ind][7]=z;//myline
					break;
				}
				if (P[mom][8]==0.)
				{
					fr+=1;
					fam1[fr]=mom;
					ind=mom;
				}
			}

			for (int w=fr; w>-1; w--)
			{
				pexe=P[fam1[w]][0];
				peye=P[fam1[w]][1];
				peze=P[fam1[w]][2];
				double we=P[fam1[w]][3];
				double previous=acum;
				double tauf=0.2*2.*P[fam1[w]][3]/pow(P[fam1[w]][4],2.);
				lifetimes[fam1[w]][1]=previous; //my special line
				if (fam1[w]==FinId[cr]) tauf=10000.;
				acum=previous+tauf;
				lifetimes[fam1[w]][2]=acum;//myline
				lifetimes[fam1[w]][5]=x;
                		lifetimes[fam1[w]][6]=y;
                		lifetimes[fam1[w]][7]=z;
				x+=P[fam1[w]][0]/we*(acum-previous);
				y+=P[fam1[w]][1]/we*(acum-previous);
				z+=P[fam1[w]][2]/we*(acum-previous);
				lifetimes[fam1[w]][8]=x;
				lifetimes[fam1[w]][9]=y;
				lifetimes[fam1[w]][10]=z;
				P[fam1[w]][8]=1.;		//Parton done
			}
		}


//----------------------------------------------------------------------------------------------------------------------myline from here
		int traceback[1000][2];//myline, this is a map from the parents to the kids, -42 means unfilled, excpet for 0 and 3
		int brothers[1000];//myline, this is a map from a sister to a sister
		int effectivemom[1000];//after all the resolution business is done, who is the real mom
		double letitlive[1000];//this is the actual time the particle should be done for
		int effectivedaughter[1000];//just a map to one of the effective immediate daughters of a particle, I dont care which
		for (int za=0; za<1000; za++)
		{
			effectivemom[za]=int(P[za][5]);
			letitlive[za]=0.;
			effectivedaughter[za]=0;
			traceback[za][0]=-42;
			traceback[za][1]=-42;
			brothers[za]=-42;
		}
		for (int mymoms=0; mymoms<1000; mymoms++)//myline
		{

				for (int mysons=0; mysons<1000; mysons++)//myline
				{
					//if (int(P[mysons][5])==mymoms && traceback[mymoms][0]!=-42 && traceback[mymoms][1]!=-42)
					//{
					//	cout<<"have less kids";
					//}
					if (int(P[mysons][5])==mymoms && traceback[mymoms][0]==-42 && int(P[mysons][5])!=0 &&int(P[mysons][5])!=3)
					{
						traceback[mymoms][0]=mysons;
						//cout<<"traceback["<<mymoms<<"][0]="<<mysons<<endl;
					}
					else if (int(P[mysons][5])==mymoms && traceback[mymoms][0]!=-42 && int(P[mysons][5])!=0 &&int(P[mysons][5])!=3)
					{
						traceback[mymoms][1]=mysons;
						brothers[mysons]=traceback[mymoms][0];
						brothers[traceback[mymoms][0]]=mysons;
						//cout<<"brothers["<<mysons<<"]="<<traceback[mymoms][0]<<endl;
						//cout<<"brothers["<<traceback[mymoms][0]<<"]="<<mysons<<endl;
						//cout<<"traceback["<<mymoms<<"][1]="<<mysons<<endl;
					}
				}
		}

		for (int zb=0; zb<1000; zb++)//this finishes filling the various times of the particles
		{
			if (brothers[zb]!=-42 && P[zb][5]!=3. && P[zb][5]!=0.) //&& (abs(P[zb][11])<=6. || P[zb][11]==21))
			{
				int mommy=int(P[zb][5]);
				int sissy=brothers[zb];
				//cout << " sissy= " << sissy << " other sissy= " << zb << " mommy= " << mommy << " lifetime= " 
				//<< lifetimes[zb][1] << " ini z= " << lifetimes[sissy][7] << " vs " << lifetimes[zb][7] << endl;
				//double store=.6;
				double store=getRes(P[mommy][3],P[mommy][0],P[mommy][1],P[mommy][2],lifetimes[zb][5],lifetimes[zb][6],lifetimes[zb][7],P[sissy][0]/P[sissy][3]-P[zb][0]/P[zb][3],P[sissy][1]/P[sissy][3]-P[zb][1]/P[zb][3],P[sissy][2]/P[sissy][3]-P[zb][2]/P[zb][3],lifetimes[zb][1]);
				//cout << " Resol= " << store << "\n";
				lifetimes[zb][0]=store;
				lifetimes[zb][3]=lifetimes[zb][0]+lifetimes[zb][1];
				lifetimes[zb][4]=lifetimes[zb][3];
				//cout <<","<<lifetimes[zb][3];
			}
		}

		int changes=0;
		while (changes==0)
		{
			changes=1;
			for (int zb=0; zb<1000; zb++)
			{
				if (lifetimes[zb][4]<lifetimes[effectivemom[zb]][4] && effectivemom[zb]!=0 && effectivemom[zb]!=3)//enforce the relationship with the mother daughter relationship
				{
					changes=0;
					lifetimes[effectivemom[zb]][4]=lifetimes[zb][4];
					//cout<<"Daughter Mother Change "<<zb<<" "<<int(effectivemom[zb])<<endl;
				}
				if (lifetimes[zb][4]!=lifetimes[brothers[zb]][4] && effectivemom[zb]!=0 && effectivemom[zb]!=3 && zb!=3 &&zb!=0)//sisters that resolve together stay together
				{
					//cout<<"Sister Sister Change"<<zb<<" "<<brothers[zb]<<endl;//lifetimes[zb][4]-lifetimes[brothers[zb]][4];
					changes=0;
					lifetimes[brothers[zb]][4]=min(lifetimes[brothers[zb]][4],lifetimes[zb][4]);
					lifetimes[zb][4]=min(lifetimes[brothers[zb]][4],lifetimes[zb][4]);
				}
			}
		}
		changes=0;
		while (changes==0)
		{
			changes=1;
			for (int zb=0; zb<1000; zb++)
			{
				if (lifetimes[zb][4]==lifetimes[effectivemom[zb]][4] && effectivemom[effectivemom[zb]]!=0 && effectivemom[effectivemom[zb]]!=3 &&zb!=3 &&zb!=0 && effectivemom[zb]!=0 && effectivemom[zb]!=3)//enforce the relationship with the mother daughter relationship
				{
					changes=0;
					effectivemom[zb]=effectivemom[effectivemom[zb]];
				}
			}
		}
//for (int zb=0; zb<60; zb++)
//			{
//				cout<<zb<<" "<<lifetimes[zb][4]<<" " <<effectivemom[zb]<<P[zb][5]<<endl;
//			}

		int a=0;
		for (int zb=1; zb<1000; zb++)//get a daughter
		{
			a=effectivemom[zb];
			//cout<<"effmom"<<zb<<"="<<a;
			effectivedaughter[a]=zb;
		}
		for (int zb=0; zb<1000; zb++)
		{
			letitlive[zb]=lifetimes[effectivedaughter[zb]][4]-lifetimes[zb][4];
			P[zb][5]=double(effectivemom[zb]);
			P[zb][8]=0.;
		}
		//start of the final loop
		int fam[30];
		double Q[1000][13];//currently only using Q[0-8 inclusive]
/*
for (int zb=0; zb<1000; zb++)
		{
			cout<< P[zb][0] << " " << P[zb][1] << " " << P[zb][2] << " " << P[zb][3]  << " " << Q[zb][0] << " " << Q[zb][1] << " " << Q[zb][2] << " " << Q[zb][3] << " " << 0. << " " << int(P[zb][6]) << " " << int(P[zb][7]) << " " << int(P[zb][9]) << " " << int(P[zb][10]) << " " << int(P[zb][11]) << "\n";
			cout<< "lifetimes"<<lifetimes[zb][0] << " " << lifetimes[zb][1] << " " << lifetimes[zb][2] << " " << lifetimes[zb][3]  << " " << lifetimes[zb][4] << " " << lifetimes[zb][5] << " " << lifetimes[zb][6]  << " " << lifetimes[zb][7]  << " " << 0. << " " << lifetimes[zb][8]  << " " << lifetimes[zb][9]  << " " << lifetimes[zb][10] << "\n";

		}
*/
		//Energy Loss Loop
		for (int cr=0; cr<k; cr++) //loop over the final particles
		{
			int ind=FinId[cr]; //get the particle number of the finished particle
			int fr=0;
			double frac;
			double acum;
			int mymomisko=0;
			fam[fr]=ind;//set the first of the family line (fam)
			for (unsigned a=0; a<50; a++)
			{
				int mom=int(P[ind][5]);
				if (P[mom][8]==1. && mom!=3) //If I have a good mom that is already energy lossed, then set the daughter's
				{
				        cout << " Getting daughter " << ind << " lambda= " << Q[mom][3]/P[mom][3] << " from mom= " << mom << endl;
					frac=Q[mom][3]/P[mom][3];
					Q[ind][0]=frac*P[ind][0];		//Px is the quenched momentum from the mom
					Q[ind][1]=frac*P[ind][1];		//Py ...
					Q[ind][2]=frac*P[ind][2];		//Pz ...
					Q[ind][3]=frac*P[ind][3];		//E is the reduced energy from the mom due to energy loss of the mom's family
					Q[ind][4]=lifetimes[ind][5]+(lifetimes[ind][4]-lifetimes[ind][1])*(P[ind][0]/P[ind][3]);		// X myline
					Q[ind][5]=lifetimes[ind][6]+(lifetimes[ind][4]-lifetimes[ind][1])*(P[ind][1]/P[ind][3]);		// Y myline
					Q[ind][6]=lifetimes[ind][7]+(lifetimes[ind][4]-lifetimes[ind][1])*(P[ind][2]/P[ind][3]);		// Z myline
					//cout<<"mylife"<<lifetimes[ind][4]-lifetimes[ind][1];
					if (Q[mom][3]==0.) mymomisko=1; //If my mom has no energy, just set the position
					acum=Q[mom][7];//+lifetimes[ind][4]-lifetimes[ind][1];		// t ...
					break;
				}
				if (mom==3) //this sets the initial parameters for the particle's line
				{
					Q[ind][0]=P[ind][0];
					Q[ind][1]=P[ind][1];
					Q[ind][2]=P[ind][2];
					Q[ind][3]=P[ind][3];
					cout << "Initial particle= " << ind << endl;
					x=xcre;
					y=ycre;
					z=0.;
					acum=0.;
					break;
				}
				if (P[mom][8]==0.) //if the mom is not done, then add to the family tree
				{
					fr+=1;//fr is the number of particles in the chain up to the last mom that was finished
					fam[fr]=mom;
					cout << "adding to family tree= " << mom << endl;
					ind=mom;
				}

			}
			for (int w=fr; w>-1; w--)
			{
				//make sure the invisible particles dont come to this part
				if (mymomisko==1)
				{
					for (int ku=w; ku>-1; ku--)
					{
                    				Q[fam[ku]][3]=0.;
                    				P[fam[ku]][8]=1.;
                    			}
                			break;
				}
				pexe=Q[fam[w]][0];
				peye=Q[fam[w]][1];
				peze=Q[fam[w]][2];
				x=lifetimes[fam[w]][5]+(lifetimes[fam[w]][4]-lifetimes[fam[w]][1])*(P[fam[w]][0]/P[fam[w]][3]);		// X myline
				y=lifetimes[fam[w]][6]+(lifetimes[fam[w]][4]-lifetimes[fam[w]][1])*(P[fam[w]][1]/P[fam[w]][3]);		// Y myline
				z=lifetimes[fam[w]][7]+(lifetimes[fam[w]][4]-lifetimes[fam[w]][1])*(P[fam[w]][2]/P[fam[w]][3]);		// Z myline
				//acum=lifetimes[fam[w]][1]+lifetimes[fam[w]][4]-lifetimes[fam[w]][1];
				double we=Q[fam[w]][3];
				double previous=acum;
				double tauf=letitlive[fam[w]];//myline
				if (fam[w]==FinId[cr]) tauf=10000.;
				acum=previous+tauf;
				double dE;
				if ((abs(P[fam[w]][11])<=6. || P[fam[w]][11]==21.))//we only apply loss to color particles
				//if ((abs(P[fam[w]][11])<=6. || P[fam[w]][11]==21.) && letitlive[fam[w]]>0.)//we only apply loss to color particles
				{
					dE=gdE(we,P[fam[w]][6],acum-previous,previous,&x,&y,&z,&pexe,&peye,&peze);
					Q[fam[w]][0]=pexe;
                    			Q[fam[w]][1]=peye;
                    			Q[fam[w]][2]=peze;
		                        cout << "Energy loss of fam[w]= " << fam[w] << " lambda= " << max(we-dE,0.)/we << " and letitlive= " << letitlive[fam[w]] << " and tauf= " << tauf << endl;
		                        if (we-dE<=0.)  //what to do when we finish a particle
					{
                    				for (int ko=w; ko>-1; ko--)
                    				{
                    					Q[fam[ko]][3]=0.;//the energy is now 0
                    					P[fam[ko]][8]=1.;//the particle is done
                    				}
                    				break;
                    			}
				}
				else //other particles just get position updates
				{
					dE=0.;
					x+=Q[fam[w]][0]/we*(acum-previous);
					y+=Q[fam[w]][1]/we*(acum-previous);
					z+=Q[fam[w]][2]/we*(acum-previous);
				}

				Q[fam[w]][3]=we-dE;//the new step of energy with the energy loss
				Q[fam[w]][4]=x;// the new position
				Q[fam[w]][5]=y;//...
				Q[fam[w]][6]=z;//...
				Q[fam[w]][7]=acum;//the new time
				P[fam[w]][8]=1.;		//Parton done
				frac=Q[fam[w]][3]/P[fam[w]][3];
				if (w>0)
				{
                			Q[fam[w-1]][0]=P[fam[w-1]][0]*frac;
                			Q[fam[w-1]][1]=P[fam[w-1]][1]*frac;
                			Q[fam[w-1]][2]=P[fam[w-1]][2]*frac;
                			Q[fam[w-1]][3]=P[fam[w-1]][3]*frac;
                		}

			}

			check << P[fam[0]][0] << " " << P[fam[0]][1] << " " << P[fam[0]][2] << " " << P[fam[0]][3]  << " " << Q[fam[0]][0] << " " << Q[fam[0]][1] << " " << Q[fam[0]][2] << " " << Q[fam[0]][3] << " " << int(P[fam[0]][9]) << " " << int(P[fam[0]][10]) << " " << int(P[fam[0]][11]) << "\n";

		}
//----------------------------------------------------------------------------

		if (k!=0) check << "NEXT" << " " << " 0. " << " " << " 0. " << " " << xcre << " " << ycre << " " << b << " " << " 0. " << " " << count << " " << " 0. " << " " << cross << " " << weight << "\n";
	        for (unsigned int i=0; i<1000; i++) {
                  if (P[i][0]==0.) continue;
		  cout << "num= " << i << " lambda= " << Q[i][3]/P[i][3] << endl;
		}

		count+=1;
                int mult=int(double(count)/(ene*100.));
                if (mult>0)
                {
                        cout << " " << count << "\n";
                        ene+=5.;
                }
	
	//End Event Loop
	} while (count < N);

	std::string snum;	
	double sigmaGen, weightSum;
	infile1 >> snum >> sigmaGen >> weightSum;
	if (snum!="END") { cout << " Problem reading file! "; exit(0); }
	check << snum << std::setprecision(9) << " " << sigmaGen << " " << weightSum;

	clock_t endClock = clock();
	cout << " Time= ";
	printf("%ld", (endClock - startClock) / CLOCKS_PER_SEC);
	cout  << " #FuckedN= " << fuckedN/totalN << " ";
	cout << " #Events= " << event << " ";

	check.close();
//End Program
}
//------------------gdE-----------------------------//
double gdE(double E, double M, double t, double t0, double *x, double *y, double *z, double *pexe, double *peye, double *peze) {

	double ti, tih, tot;
	double proper, properh;
	double xp, yp, zp;
	double xph, yph, zph;
	double Temp, Temph;
	double Ev, CF;
	double Tc;
	if (tmethod==0) Tc=170./1000.;//critical energy 0
	if (tmethod==1) Tc=145./1000.;//critical energy 1
	double stepe=0.01;
	double wx, wy, wz;
	double vx, vy, vz;
	double lore, vscalw, v2;
	int marker=0;
	double i=0.;
	double quench=0.;
	double gede=0.;
	double dist=0.;
	double distlab=0.;
	if (int(M)==1) CF=1.; //quark
	if (int(M)==2)
	{ //gluon
		if (qmethod==3) CF=pow(9./4.,1./3.);
		else CF=pow(9./4.,1.);
	}
	Ev=E;
	tot=t0+t;//final time
	int noproblem=1;
	int saidbef=0, saidaft=0;
	//cout << "\n New Rung--- Ei= " << E << "\n";
	//cout << " gamma= " << gamma << "\n";
	//cout << " xini= " << *x << " yini= "<< *y << " zini= " << *z << " tauini= " << sqrt(t0*t0-pow(*z,2.))<< "\n";
	//cout << " Initial Q2= " << E*E-pow(*pexe,2.)-pow(*peye,2.)-pow(*peze,2.) << "\n";
	do {
		ti=t0+i*stepe;//time at this substep
		tih=ti+stepe;//time at the next substep for derivatives I guess
		if (tih>tot) //if we are done with the time we are given
		{
			marker=1;
			tih=tot;
			stepe=tih-ti;//if we run out of time, we dont get to go over
		}
		wx=*pexe/Ev;//momentum over E at this step (reduced momentum)
		wy=*peye/Ev;//momentum over E at this step
		wz=*peze/Ev;//momentum over E at this step

		xp=*x;//initial x pos
		yp=*y;//initial y pos
		zp=*z;//initial z pos
		proper=sqrt(ti*ti-zp*zp);//some proper time in the z direction

		xph=*x+wx*stepe;
        yph=*y+wy*stepe;
        zph=*z+wz*stepe;
        properh=sqrt(tih*tih-zph*zph);

        *x=xph;//update positions
        *y=yph;
        *z=zph;
        if (proper!=proper) proper=0.;
        if (properh!=properh) properh=0.;
		if (proper>=tquench && properh>=tquench) {
			//extract fluid velocity with LAB coordinates
			vx=gVx(proper,xp,yp);//get fluid velo
			vy=gVy(proper,xp,yp);//get fluid velo
			double frap=0.;// this is something
			vz=*z/tih;
			frap=atanh(vz);//rapidity!
			vx/=cosh(frap);//adjust the rest
			vy/=cosh(frap);
			v2=pow(vx,2.)+pow(vy,2.)+pow(vz,2.);
			vscalw=vx*wx+vy*wy+vz*wz;
			if (v2>=1.) {
				v2=0.9999999;
				cout << " V2 > 1 WTFFF !!! \n";
			}
			lore=1./sqrt(1.-v2);//gamma factor
			double infi=stepe*sqrt(wx*wx+wy*wy+wz*wz+lore*lore*(v2-2.*vscalw+vscalw*vscalw));
			if (infi!=infi) infi=0.;
			dist+=infi;
			distlab+=stepe;
			Temp=gT(proper,xp,yp);
			Temph=gT(properh,xph,yph);
			if (Temph<Tc) {//if we are done making it lose energy
				double cospart=xph*wx+yph*wy;
				if (cospart>0.) {
					Temp=0.;
					Temph=0.;
					*x=xph+wx*(tot-tih);//new positions
					*y=yph+wy*(tot-tih);
					*z=zph+wz*(tot-tih);
				}
			}
			if (Temph==0.) marker=1;//if we went through the last loop, we are done 
			if (Ev<=0.) marker=1;// We cant subtract energy from 0 energy
		}
		if (Ev>0.) {
			if (proper>=tquench && properh>=tquench && (Temp+Temph)/2.>Tc+0.00000001) {
		if (kappa!=0.) { //Kappa condition for broadening		
				//Broadening
				double e1x=wy*vz-wz*vy;
				double e1y=wz*vx-wx*vz;
				double e1z=wx*vy-wy*vx;
				double Ne1=sqrt(pow(e1x,2.)+pow(e1y,2.)+pow(e1z,2.));
				if (Ne1!=Ne1) cout << " PROBLEMMMM WITH NE1 ";
				e1x=e1x/Ne1;
				e1y=e1y/Ne1;
				e1z=e1z/Ne1;

				double Nw=sqrt(pow(wx,2.)+pow(wy,2.)+pow(wz,2.));
				if (Nw!=Nw) cout << " PROBLEMMMM WITH NW ";
				double lx=(wy*e1z-wz*e1y)/Nw;
				double ly=(wz*e1x-wx*e1z)/Nw;
				double lz=(wx*e1y-wy*e1x)/Nw;

				//cout << " Modul L= " << lx*lx+ly*ly+lz*lz << " \n";
				double uscalW=lore*(1.-vx*wx-vy*wy-vz*wz);
				double uscall=lore*(-vx*lx-vy*ly-vz*lz);
				double W2=1.-wx*wx-wy*wy-wz*wz;
				//if (wx>1. || wy>1. || wz>1.) cout << " at i= " << i << " W2 loker= " << W2 << " with Ev= " << Ev << " wx= " << wx << " wy= " << wy << " wz= " << wz << "\n";
				double Wpt=1.-lore*W2/uscalW;
				double Wpx=wx-vx*lore*W2/uscalW;
				double Wpy=wy-vy*lore*W2/uscalW;
				double Wpz=wz-vz*lore*W2/uscalW;

				double NN=1.+W2*pow(uscall,2.)/(-pow(uscalW,2.)+W2);
				double Nalpha=-uscall*uscalW/(pow(uscalW,2.)-W2);
				if (v2>1.-W2) {
					//cout << " vfluid= " << v2 << " vparton= " << 1.-W2 << " ";
					//cout << " nume/deno= " << W2*pow(uscall,2.)/(-pow(uscalW,2.)+W2) << " ";
					//cout << " W2= " << W2 << "\n";
				}
				if (sqrt(NN)!=sqrt(NN)) {
					noproblem=0;
					//cout << " Cos= " << (vx*wx+vy*wy+wz*vz)/sqrt(v2)/sqrt(1.-W2) << " ";
					//cout <<  " deno N= " << -pow(uscalW,2.)+W2 << " Ev= " << Ev << "\n";
					//cout << " nume/deno= " << W2*pow(uscall,2.)/(-pow(uscalW,2.)+W2) << " ";
					//cout << " vfluid= " << v2 << " vparton= " << 1.-W2 << " W2= " << W2 << " lore= " << lore << "\n";
				}
				else {
					double e2t=Nalpha*Wpt/sqrt(NN);
					double e2x=(lx+Nalpha*Wpx)/sqrt(NN);
					double e2y=(ly+Nalpha*Wpy)/sqrt(NN);
					double e2z=(lz+Nalpha*Wpz)/sqrt(NN);

					double Ef=Ev*lore*(1.-vscalw);
					double wf2=1.-W2/pow(lore*(1.-vscalw),2.);
					double DelQ2=kappa*pow((Temp+Temph)/2.,3.)*lore*(1.-vscalw)*stepe*5.;
					double qfac=0.;
					int bucle=0;
					if (Ef>(Temp+Temph)/2.) {
						do {
							bucle+=1;
							if (kappa!=0.) qfac=gQ(DelQ2);
							if (bucle>1000) cout << " Ef= " << Ef << " ";
						} while (qfac>Ef*sqrt(wf2) && bucle<1000);
					}
					double qbeta=sqrt(1.-qfac*qfac/Ef/Ef/wf2)-1.;
					if (qbeta!=qbeta) {
						qbeta=-1.;
						qfac=Ef*sqrt(wf2);
						if (qfac!=qfac) cout << " QFAC LOKKKKKKERRRR ";
					}
					double qphi=2.*3.141592654*rando(&ir);

					double et=sin(qphi)*e2t;
					double ex=cos(qphi)*e1x+sin(qphi)*e2x;
					double ey=cos(qphi)*e1y+sin(qphi)*e2y;
					double ez=cos(qphi)*e1z+sin(qphi)*e2z;

					double Wtt=(1.-(uscalW)*lore)/lore/(1.-vscalw);
					double Wtx=(wx-(uscalW)*vx*lore)/lore/(1.-vscalw);
					double Wty=(wy-(uscalW)*vy*lore)/lore/(1.-vscalw);
					double Wtz=(wz-(uscalW)*vz*lore)/lore/(1.-vscalw);

					if (W2<0. && saidbef==0) {
						//cout << " Q² bef= " << Ev*Ev-pow(*pexe,2.)-pow(*peye,2.)-pow(*peze,2.);
						saidbef=1;
					}
					//cout << " Echange= " << qbeta*Ef*Wpt+qfac*et << " ";
					if (qbeta*Ef*Wtt+qfac*et!=qbeta*Ef*Wtt+qfac*et) cout << " Ener WRONG ";
					if (qbeta*Ef*Wtx+qfac*ex!=qbeta*Ef*Wtx+qfac*ex) cout << " Px WRONG";
					if (qbeta*Ef*Wty+qfac*ey!=qbeta*Ef*Wty+qfac*ey) cout << " Py WRONG";
					if (qbeta*Ef*Wtz+qfac*ez!=qbeta*Ef*Wtz+qfac*ez) cout << " Pz WRONG";
					Ev+=qbeta*Ef*Wtt+qfac*et;
					if (Ev<=0.) {
						//cout << " at i= " << i << " Ev= " << Ev << " ";
						Ev=0.;
						marker=1;
					}
					*pexe+=qbeta*Ef*Wtx+qfac*ex;
					*peye+=qbeta*Ef*Wty+qfac*ey;
					*peze+=qbeta*Ef*Wtz+qfac*ez;
					if (W2<0. && saidaft==0) {
						//cout << " Q² aft= " << Ev*Ev-pow(*pexe,2.)-pow(*peye,2.)-pow(*peze,2.) << "\n";
						saidaft=1;
					}
				}
		} //Kappa condition for broadening
				//Coll
				if (qmethod==1 && Ev>0.) {
					double intpiece=alpha*CF*stepe*5.*pow((Temp+Temph)/2.,2.);
					quench=(Ev-intpiece)/Ev;
					Ev-=intpiece;//energy loss
					if (Ev<0.) Ev=0.;
				}
				//Radiative
				if (qmethod==2 && Ev>0.) {
					double intpiece=5.*alpha*CF*stepe*dist*5.*pow((Temp+Temph)/2.,3.);
					quench=(Ev-intpiece)/Ev;
					Ev-=intpiece;//energy loss
					if (Ev<0.) Ev=0.;
				}
				//Strong
				if (qmethod==3 && Ev>0.) {
					double Efs=E*lore*(1.-vscalw);//inital energy of particle (not quenched) (lore is gamma factor)
					double tstop=0.2*pow(Efs,1./3.)/(2.*pow((Temp+Temph)/2.,4./3.)*alpha)/CF; //stopping distance
					double beta=tstop/dist;
					if (beta>1.) {
						double intpiece=Efs*stepe*4./(3.141592)*(1./(beta*tstop*sqrt(beta*beta-1.)));
						quench=(Ev-intpiece)/Ev;
						Ev-=intpiece; //energy of loss in timestep
						if (intpiece<0.) cout << " WTF INTPIECE= " << intpiece << " ";
						if (quench>1.) cout << " QuenchFactor= " << quench << "  and IntPiece= " << intpiece << " and Ev= " << Ev << " ";
					}
					else {
						Ev=0.;
						quench=0.;
						marker=1;
					}
					if (Ev<0.) {
						Ev=0.;
						quench=0.;
						marker=1;
					}
				}
				if (quench>1.) cout << " QuenchFactor= " << quench << " ";
				if (Ev<=0.) quench=0.;
				*pexe*=quench;//momentums quenched
				*peye*=quench;
				*peze*=quench;
			}
		}
		else {
			gede=E;
			Ev=0.;
			marker=1;
		}
		i+=1.;
	} while (marker!=1);
	//cout << " Final Q2= " << Ev*Ev-pow(*pexe,2.)-pow(*peye,2.)-pow(*peze,2.) << "\n";
	//cout << " Final 4mom= " << Ev << " " << *pexe << " " << *peye << " " << *peze << " quench= " << quench << "\n";
	gede=E-Ev;//energy chance
	if (noproblem!=1) fuckedN+=1.;
	totalN+=1.;
	//cout << " xfin= " << *x << " yfin= " << *y << " zfin= " << *z <<" Temp fin= " << (Temp+Temph)/2. << "\n";
	//cout << " DistLab= " << distlab << " DistFluid= " << dist << "\n";
	//cout << " dE= " << gede << " ";
	return gede;
}

//-----------------resolve time (mine)----------------------//
//calculates resolve time of a particle only considering the particle, its parent, and its sister
double getRes(double mome,double momPx,double momPy,double momPz, double momx, double momy, double momz, double DPx,double DPy,double DPz, double t0)
//the names suck. DP is acutally dvelocity now
{ //momx,momy, and momz are the finishing postions, or the starting postions of its daughters
	double myscale=1/3.14/rpower;
	double ti;
	double proper;
	//double xp, yp, zp;
	double xph, yph, zph;
	double Temp;
	double Tc;
	if (tmethod==0) Tc=170./1000.;//critical energy 0
        if (tmethod==1) Tc=145./1000.;//critical energy 1
	double wx, wy, wz;
	//double vx, vy, vz;
	//double v2;
	ti=t0;

        if (rpower<0.000001) return 100000000000;

        //REVIEW this function...

	double stepe=.01;

	wx=momPx/mome;//momentum over E at this step (velocities)
	wy=momPy/mome;//momentum over E at this step
	wz=momPz/mome;//momentum over E at this step
	proper=sqrt(t0*t0-momz*momz);
	//cout << t0 << " ,  " << momz << ", ";
	//cout<<" Ini proper= " << proper << endl;

	stepe=0.1;
	double xprime=0., yprime=0., zprime=0.;
        int dontadv=0;
	for (unsigned i=0; i<=100000000; i++)
	{
		if (i==1000000)
		{
			cout << "went too far";
		}

		double propermom=sqrt(ti*ti-momz*momz);

		if (propermom>=tquench) {
			//vx=gVx(propermom,momx,momy);//get fluid velo
	                //vy=gVy(propermom,momx,momy);//get fluid velo
			//double frap=0.;// this is something
			//vz=momz/ti;
			//frap=atanh(vz);//rapidity!
			//vx/=cosh(frap);//adjust the rest
			//vy/=cosh(frap);
			//v2=pow(vx,2.)+pow(vy,2.)+pow(vz,2.);
			//double myv=sqrt(v2);
			Temp=gT(propermom,momx,momy);
			//double mygamma=1/sqrt(1-v2); //check the lorentz transforms
			//double xprime=((1+(mygamma-1)*(vx/myv)*(vx/myv))*(DPx)+(mygamma-1)*(vx/myv)*(vy/myv)*(DPy)+(mygamma-1)*(vx/myv)*(vz/myv)*(DPz))*i*stepe;
			//double yprime=((mygamma-1)*(vx/myv)*(vy/myv)*(DPx)+(1+(mygamma-1)*(vy/myv)*(vy/myv))*(DPy)+(mygamma-1)*(vy/myv)*(vz/myv)*(DPz))*i*stepe;
			//double zprime=((mygamma-1)*(vx/myv)*(vz/myv)*(DPx)+(mygamma-1)*(vz/myv)*(vy/myv)*(DPy)+(1+(mygamma-1)*(vx/myv)*(vz/myv))*(DPz))*i*stepe;
			//double xprime=DPx*double(i)*stepe;
			//double yprime=DPy*double(i)*stepe;
			//double zprime=DPz*double(i)*stepe;

			//if (sqrt(pow(Q[sister][0]-Q[particle][0],2)+pow(Q[sister][1]-Q[brothers[particle]][1],2)+pow(Q[sister][2]-Q[particle][2],2))*i*stepe>=myscale/Temp)
			if (sqrt(pow(xprime,2)+pow(yprime,2)+pow(zprime,2))>=myscale/(Temp/0.2) || Temp<=Tc)
			{
				if (stepe==0.1) {
                                        if (int(i)==0) {
                                                //cout << " dist= " << sqrt(pow(xprime,2)+pow(yprime,2)+pow(zprime,2)) << " temp= " << Temp << endl;
                                                return ti-t0;
                                        } else dontadv=1;
                                } else return ti-t0;
			}
		}

		if (dontadv==0) {
                        xph=momx+wx*stepe;
                        yph=momy+wy*stepe;
                        zph=momz+wz*stepe;
                        momx=xph;
                        momy=yph;
                        momz=zph;
                        xprime+=DPx*stepe;
                        yprime+=DPy*stepe;
                        zprime+=DPz*stepe;
                        ti+=stepe;
                } else {
                        momx-=wx*stepe;
                        momy-=wy*stepe;
                        momz-=wz*stepe;
                        xprime-=DPx*stepe;
                        yprime-=DPy*stepe;
                        zprime-=DPz*stepe;
                        ti-=stepe;
                        stepe=0.01;
                        dontadv=0;
                }
	}
	return 100000000000;
	//taum= xdiff/sqrt(pow(Q[fam[w]][0]-Q[brothers[fam[w]]][0],2)+pow(Q[fam[w]][1]-Q[brothers[fam[w]]][1],2)+pow(Q[fam[w]][2]-Q[brothers[fam[w]]][2],2)); //myline
}

//------------------getGrid-------------------------//
double getGrid(double tau, double x, double y, double *dt, double *dx, double *dy, int *it, int *ix, int *iy) {

	*it=int((tau-tau0)/deltat);
	*dt=(tau-tau0-double(*it)*deltat)/deltat;

	if (y>0.) {
		*iy = int(y/deltay)+maxy/2;
		*dy = (y - double(*iy-maxy/2.)*deltay)/deltay;
	}
	else {
		*iy = int(y/deltay)+maxy/2-1;
		*dy = (y - double(*iy-maxy/2.)*deltay)/deltay;
	}

	if (x>0.) {
		*ix = int(x/deltax)+maxx/2;
		*dx = (x - double(*ix-maxx/2.)*deltax)/deltax;
	}
	else {
		*ix = int(x/deltax)+maxx/2-1;
		*dx = (x - double(*ix-maxx/2.)*deltax)/deltax;
	}

	return 0;
}

//------------------gT------------------------------//
double gT(double tau, double x, double y) {
// temp in fluid rest frame
	double gete=0.;
	double tau1=18.5;

	if (tau>tau1) {
		//cout << " No data for tau > " << tau1 << " fm ";
		return gete;
	}

	getGrid(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);
	//cout << " " << dt << " " << dx << " " << dy << " " << dh << " " << it << " " << ix << " " << iy << " " << ih << "\n";
	if (ix<0 || ix>maxx || iy<0 || iy>maxy) return gete;
	gete=hydrot[it][ix][iy]*(1.-dt)*(1.-dx)*(1.-dy);
	gete+=hydrot[it+1][ix][iy]*(1.-dx)*(1.-dy)*dt;
	gete+=hydrot[it][ix+1][iy]*(1.-dt)*(1.-dy)*dx;
	gete+=hydrot[it][ix][iy+1]*(1.-dx)*(1.-dt)*dy;
	gete+=hydrot[it+1][ix+1][iy]*dx*(1.-dy)*dt;
	gete+=hydrot[it][ix+1][iy+1]*(1.-dt)*dy*dx;
	gete+=hydrot[it+1][ix][iy+1]*(1.-dx)*dt*dy;
	gete+=hydrot[it+1][ix+1][iy+1]*dx*dt*dy;

	return gete*0.2;
}



//------------------gE------------------------------//
double gE(double tau, double x, double y) {

	double gener=0.;
	double tau1=18.5;

	if (tau>tau1) {
		//cout << " No data for tau > " << tau1 << " fm ";
		return gener;
	}

	getGrid(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);
	//cout << " " << dt << " " << dx << " " << dy << " " << dh << " " << it << " " << ix << " " << iy << " " << ih << "\n";
	if (ix<0 || ix>maxx || iy<0 || iy>maxy) return gener;
	gener=hydroe[it][ix][iy]*(1.-dt)*(1.-dx)*(1.-dy);
	gener+=hydroe[it+1][ix][iy]*(1.-dx)*(1.-dy)*dt;
	gener+=hydroe[it][ix+1][iy]*(1.-dt)*(1.-dy)*dx;
	gener+=hydroe[it][ix][iy+1]*(1.-dx)*(1.-dt)*dy;
	gener+=hydroe[it+1][ix+1][iy]*dx*(1.-dy)*dt;
	gener+=hydroe[it][ix+1][iy+1]*(1.-dt)*dy*dx;
	gener+=hydroe[it+1][ix][iy+1]*(1.-dx)*dt*dy;
	gener+=hydroe[it+1][ix+1][iy+1]*dx*dt*dy;

	return gener;
}

//------------------gVx------------------------------//
double gVx(double tau, double x, double y) {

        double gvelx=0.;
        double tau1=18.5;

        if (tau>tau1) {
                //cout << " No data for tau > " << tau1 << " fm ";
                return gvelx;
        }

        getGrid(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);
	//cout << " " << dt << " " << dx << " " << dy << " " << dh << " " << it << " " << ix << " " << iy << " " << ih << "\n";
	if (ix<0 || ix>maxx || iy<0 || iy>maxy) return gvelx;
	gvelx=hydrox[it][ix][iy]*(1.-dt)*(1.-dx)*(1.-dy);
	gvelx+=hydrox[it+1][ix][iy]*(1.-dx)*(1.-dy)*dt;
	gvelx+=hydrox[it][ix+1][iy]*(1.-dt)*(1.-dy)*dx;
	gvelx+=hydrox[it][ix][iy+1]*(1.-dx)*(1.-dt)*dy;
	gvelx+=hydrox[it+1][ix+1][iy]*dx*(1.-dy)*dt;
	gvelx+=hydrox[it][ix+1][iy+1]*(1.-dt)*dy*dx;
	gvelx+=hydrox[it+1][ix][iy+1]*(1.-dx)*dt*dy;
	gvelx+=hydrox[it+1][ix+1][iy+1]*dx*dt*dy;

        return gvelx;
}

//------------------gVy------------------------------//
double gVy(double tau, double x, double y) {

        double gvely=0.;
        double tau1=18.5;

        if (tau>tau1) {
                //cout << " No data for tau > " << tau1 << " fm ";
                return gvely;
        }

        getGrid(tau,x,y,&dt,&dx,&dy,&it,&ix,&iy);
	//cout << " " << dt << " " << dx << " " << dy << " " << dh << " " << it << " " << ix << " " << iy << " " << ih << "\n";
	if (ix<0 || ix>maxx || iy<0 || iy>maxy) return gvely;
	gvely=hydroy[it][ix][iy]*(1.-dt)*(1.-dx)*(1.-dy);
	gvely+=hydroy[it+1][ix][iy]*(1.-dx)*(1.-dy)*dt;
	gvely+=hydroy[it][ix+1][iy]*(1.-dt)*(1.-dy)*dx;
	gvely+=hydroy[it][ix][iy+1]*(1.-dx)*(1.-dt)*dy;
	gvely+=hydroy[it+1][ix+1][iy]*dx*(1.-dy)*dt;
	gvely+=hydroy[it][ix+1][iy+1]*(1.-dt)*dy*dx;
	gvely+=hydroy[it+1][ix][iy+1]*(1.-dx)*dt*dy;
	gvely+=hydroy[it+1][ix+1][iy+1]*dx*dt*dy;

        return gvely;
}

//------------------gQ------------------------------//
double gQ(double Del) {

	double gaussq;
	double qfac;
	double gaussmax=sqrt(2./Del/exp(1.));

	qhatelsen:
	qfac=5.*sqrt(Del/2.)*rando(&ir);
	gaussq=2.*qfac/Del*exp(-qfac*qfac/Del)/gaussmax;
	double nrand=rando(&ir);
	if (nrand>gaussq) goto qhatelsen;

return qfac;
}

//------------------gxy-----------------------------//
double gxy(double *x, double *y, double *b) {

	double rho,phi;
	double P;

        naiguels:
	*b=sqrt((bmax*bmax-bmin*bmin)*rando(&ir)+bmin*bmin);
        norm=1.;
        norm=gTAA(0.,0.,bmin);
	rho=sqrt(150.*rando(&ir));
	phi=2.*3.141592654*rando(&ir);
	*x=rho*cos(phi);
	*y=rho*sin(phi);
	P=rando(&ir);
	if(P>gTAA(*x,*y,*b)) goto naiguels;
	return 0;
}

//-----------------gTAA------------------------------//
double gTAA(double x, double y, double b)
{
	int il, irr;
	double rho2, use;

	rho2=pow(x+b/2.,2.)+y*y;
	il=int(rho2/step);
	rho2=pow(x-b/2.,2.)+y*y;
	irr=int(rho2/step);
	use=0.;
	if(il<4000 && irr<4000) {
		use=TA[il]*TA[irr]/norm;
	}
	return use;
}

//--------------Random Number-----------------------//
double rando (int *ir)
{
	double da=16807.;
	double db=2147483647.;
	double dc=2147483648.;
	double usran;

	*ir=int(fabs(fmod(da*(*ir),db)+0.5));
	usran=double(*ir)/dc;
	return usran;
}
