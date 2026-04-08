#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>

// uncomment the line below to activate all additional printings
//#define DO_PRINT

using namespace std;
int ir;
double rando(int*);
double factorial(int);

//Args: File, Seed, N
int main (int argc, char** argv) {
	printf(" inpargc=%d inpFile=%s\n",argc,argv[1]);
	assert(argc==4);
	char inpFile[100], outFile[100];
	sprintf(inpFile,"./sor_%s.out",argv[1]);
	sprintf(outFile,"./wake_%s.out",argv[1]);
	printf(" inp=%s  out=%s\n",inpFile, outFile);
	ir=1346+atoi(argv[2]);

#ifdef DO_PRINT
	//BackFiles
	ofstream phifile;
        phifile.open ("/home/peibols/ADSPYTHIA/Backreaction/Dists/Back/Phidebug.txt");
        ofstream rapfile;
        rapfile.open ("/home/peibols/ADSPYTHIA/Backreaction/Dists/Back/Rapdebug.txt");
        ofstream ptfile;
        ptfile.open ("/home/peibols/ADSPYTHIA/Backreaction/Dists/Back/Ptdebug.txt");
        ofstream nfile;
        nfile.open ("/home/peibols/ADSPYTHIA/Backreaction/Dists/Back/Ndebug.txt");
	ofstream raphi;
	raphi.open ("/home/peibols/ADSPYTHIA/Backreaction/Dists/Back/Raphi.txt");
	ofstream cons;
        cons.open ("/home/peibols/ADSPYTHIA/Backreaction/Dists/Back/Cons.txt");
#endif
	ofstream savedback;
	savedback.open (outFile);

	std::ifstream infile1(inpFile);
	assert(!infile1.fail());

	int Ntot=atoi(argv[3]);
	int count=0;

	double transcut=0.;
	double basesig=0.65;
	double msigma=0.65;
	int Nrun=800000;
	double maxptsq=pow(3.5,2.);
	double maxrap=2.5;
	double rbin=2.*maxrap/200.;
	double tole=0.25;
	double masspi=0.1396;
	double masspro=0.938;
	double masstra[2];
	masstra[0]=masspi;
	masstra[1]=masspro;
        double normcoop[2];
        normcoop[0]=30.;
        normcoop[1]=30.;
        double maxcooper[2];
        maxcooper[0]=0.;
        maxcooper[1]=0.;

	int toomuch=0;
	int totalan=0;
	int totalsave=0;

	double passdist[4][20][1];
	for (unsigned a=0; a<4; a++) {
		for (unsigned b=0; b<20; b++) {
			for (unsigned c=0; c<1; c++) {
				passdist[a][b][c]=0.;
			}
		}
	}

	//Back Histos Ini
	double phirap[201][201][2];
	for (unsigned a=0; a<201; a++) {
		for (unsigned b=0; b<201; b++) {
			for (unsigned c=0; c<2; c++) {
				phirap[a][b][c]=0.;
			}
		}
	}
	double phisto[201][2];
	for (unsigned a=0; a<201; a++) {
		for (unsigned b=0; b<2; b++) {
			phisto[a][b]=0.;
		}
	}
	double risto[201][2];
	for (unsigned a=0; a<201; a++) {
		for (unsigned b=0; b<2; b++) {
			risto[a][b]=0.;
		}
	}
	double ptisto[200][2];
	for (unsigned a=0; a<200; a++) {
		for (unsigned b=0; b<2; b++) {
			ptisto[a][b]=0.;
		}
	}
	double nmed=0.;
	double ntwomed=0.;
	double nisto[100];
	for (unsigned a=0; a<100; a++) {
		nisto[a]=0.;
	}
	double Yield[2];
	Yield[0]=0.;
	Yield[1]=0.;

	//Event Loop
	do {
		//Give Fastjet the vacuum and quenched momenta
		double Px, Py, Pz, E, Qx, Qy, Qz, Qe, ide, col, acol;
		std::string sQx,sQy,sQz,sQe;
		std::string sPx;
		int k=0;
		double lostx[500], losty[500], lostz[500], loste[500];
		for (unsigned a=0; a<500; a++) {
			lostx[a]=0.;
			losty[a]=0.;
			lostz[a]=0.;
			loste[a]=0.;
		}
		while (infile1 >> sPx >> Py >> Pz >> E >> sQx >> sQy >> sQz >> sQe >> col >> acol >> ide)
                {
                        Qx = strtod(sQx.c_str(), NULL);
                        Qy = strtod(sQy.c_str(), NULL);
                        Qz = strtod(sQz.c_str(), NULL);
                        Qe = strtod(sQe.c_str(), NULL);
                        if (Qx!=Qx) Qx=0.;
                        if (Qy!=Qy) Qy=0.;
                        if (Qz!=Qz) Qz=0.;
                        if (Qe!=Qe) Qe=0.;
                        if (sPx == "NEXT") {
                                break;
                        }
                        else {
				Px=strtod(sPx.c_str(), NULL);
                                if (Qe==0.) {
                                        Qx=0.;
                                        Qy=0.;
                                        Qz=0.;
                                }
                                lostx[k]=Px-Qx;
                                losty[k]=Py-Qy;
                                lostz[k]=Pz-Qz;
                                loste[k]=E-Qe;
                                k+=1;
                        }
                }
		//cout << " k= " << k << "\n";
		//Backreaction Montecarlo
		double bomx[2000], bomy[2000], bomz[2000], bome[2000], boml[2000], bomm[2000], pxrand, pyrand;
		int kk=0;
		for (int wa=0; wa<k; wa++) {
			double ptlost=sqrt(pow(lostx[wa],2.)+pow(losty[wa],2.));
                        double m2lost=pow(loste[wa],2.)-pow(ptlost,2.)-pow(lostz[wa],2.);
			double raplost=1./2.*log((sqrt(pow(ptlost,2.)+pow(lostz[wa],2.))+lostz[wa])/(sqrt(pow(ptlost,2.)+pow(lostz[wa],2.))-lostz[wa]));
                        //double raplost=1./2.*log((loste[wa]+lostz[wa])/(loste[wa]-lostz[wa]));
                        //if (raplost!=raplost) cout << " Raplost prob= " << raplost << " ptlost= " << ptlost << " pzlost= " << lostz[wa] << endl;
			if (m2lost<0.) {
                        	//cout << " Negative Mass Squared!! mass2= " << m2lost << "\n";
				m2lost=0.;
                        }
                        //double mtlost=sqrt(pow(ptlost,2.)+m2lost);
			double mtlost=loste[wa]/cosh(raplost);
			if (loste[wa]!=0. && loste[wa]<100000000. && mtlost>transcut) {
				totalan+=1;
				double Temp;
				double eneback=0.;
				double pxback=0.;
				double pyback=0.;
				double pzback=0.;
				double phidif;
				double cooper=0.;
				double randian=1.;
				double ptrand=0.;
				double mtrand=0.;
				double phirand=0.;
				double raprand=0.;
				//cout << " M2lost= " << sqrt(m2lost) << "\n";
                                double aa=3.*ptlost/mtlost;
				if (mtlost!=mtlost) cout << " percal pues \n";
				if (aa!=aa) cout << " pt lost= " << ptlost << "\n";
				//cout << " LostE= " << loste[wa] << " LostPx= " << lostx[wa] << " LostPy= " << losty[wa] << " LostPz= " << lostz[wa] << " a= " << wa << "\n";
				//cout << " Mlost/Ptlost= " << sqrt(m2lost)/ptlost << "\n";
				//cout << " Rap= " << raplost << " PseudoRap= " << pseudoraplost << "\n";
				int ntracks=0;
				int dx, dy, dz, de, Npart;
				int spe;
				double difx, dify, difz, dife;
				double pass, newpass, passrand, sign, ppxback, ppyback, ppzback, peneback;
				int runi, encallao, you;
				int numenc=0;
				//Nrun=4000;
				double clocklim=0.1;
				int tooclock=0;
				thisis:
				clock_t startClock = clock();
				kk=0;
				eneback=0.;
				pxback=0.;
				pyback=0.;
				pzback=0.;
				ntracks=0;
				cooper=0.;
				randian=1.;
				do {
					if (rando(&ir)<=0.05) spe=1;
                                        else spe=0;
					do {
						ptrand=sqrt(maxptsq*rando(&ir));
						if (spe==0) {
							Temp=0.211501*pow(ptrand,0.275362);
                                        		if (Temp>0.4) Temp=0.4;
                                        		if (Temp<0.19) Temp=0.19;
						}
						else {
							Temp=0.33*pow(ptrand,0.3);
                                        		if (Temp>0.4) Temp=0.4;
                                        		if (Temp<0.15) Temp=0.15;
						}
						phirand=2.*3.141592*rando(&ir);
						raprand=maxrap*(-1.+2.*rando(&ir));
						pxrand=ptrand*cos(phirand);
						pyrand=ptrand*sin(phirand);
						mtrand=sqrt(pow(ptrand,2.)+pow(masstra[spe],2.));
						phidif=acos((lostx[wa]*pxrand+losty[wa]*pyrand)/(ptlost*ptrand));
						cooper=exp(-mtrand/Temp*cosh(raprand))*mtrand/pow(Temp,5.)*cosh(raprand)*(ptrand*aa*cos(phidif)+mtrand*cosh(raprand))/normcoop[spe];
						//if (cooper<0.) cout << " NegaCooo phidif= " << phidif << " ";
						if (cooper!=cooper) cooper=0.;
						if (fabs(cooper)>maxcooper[spe]) maxcooper[spe]=fabs(cooper);
						if (fabs(cooper)>1.) {
							//cout << " coooooOOOOoooOooOOper= " << cooper << " at ptrand= " << ptrand << "\n";
							normcoop[spe]*=(fabs(cooper)+0.0001);
							cooper=0.;
						}
						randian=rando(&ir);
					} while (fabs(cooper)<randian);
					//if (count==9266) cout << " EneBack= " << eneback << "\n";
					bomx[kk]=pxrand;
					bomy[kk]=pyrand;
					//bome[kk]=sqrt(pow(masstra[spe],2.)+pow(ptrand,2.))*cosh(raprand);
                                        //bomz[kk]=bome[kk]*tanh(raprand);
					bomz[kk]=ptrand*sinh(raprand+raplost);
					bome[kk]=sqrt(pow(ptrand*cosh(raprand+raplost),2.)+pow(masstra[spe],2.));
					bomm[kk]=masstra[spe];
					if (cooper>0.) {
						boml[kk]=0.;
						eneback+=bome[kk];
						ntracks+=1;
					}
					else {
						boml[kk]=1.;
						eneback-=bome[kk];
						ntracks-=1;
					}
					kk+=1;
					cooper=0.;
					randian=1.;
					//cout  << " Pt of feed= " << ptrand << " and Pz= " << ptrand*sinh(raprand) << "\n";
				} while (eneback<loste[wa]);
				//cout << " Ntracks= " << ntracks << " alltracks= " << kk << "\n";
				//Correction
				pxback=0.;
				pyback=0.;
				pzback=0.;
				eneback=0.;
				for (int b=0; b<kk; b++) {
					//cout << " sign= " << pow(-1.,boml[wa]) << "\n";
					pxback+=pow(-1.,boml[b])*bomx[b];
					pyback+=pow(-1.,boml[b])*bomy[b];
					pzback+=pow(-1.,boml[b])*bomz[b];
					eneback+=pow(-1.,boml[b])*bome[b];
				}
				difx=fabs(lostx[wa]-pxback);		// fabsDist Values!
				dify=fabs(losty[wa]-pyback);
				difz=fabs(lostz[wa]-pzback);
				dife=fabs(loste[wa]-eneback);
				msigma=basesig;
                                msigma=sqrt(pow(difx,2.)+pow(dify,2.)+pow(difz,2.)+pow(dife,2.))/sqrt(log(2.));
				//cout << " DPx= " << difx << " DPy= " << dify << " DPz= " << difz << " DE= " << dife << "\n";
				pass=exp(-pow(difx/msigma,2.)/2.)*exp(-pow(dify/msigma,2.)/2.)*exp(-pow(difz/msigma,2.)/2.)*exp(-pow(dife/msigma,2.)/2.);
				//cout << " Pass0= " << pass << "\n";
				runi=0;
				encallao=0;
				do {
				//Choose random particle
					you=int(double(kk)*rando(&ir));
					if (you==kk) cout << " rando dio 111111111 WTF ";
					do {
						ptrand=sqrt(maxptsq*rando(&ir));
                                                if (bomm[you]<0.5) {
							spe=0;
							Temp=0.211501*pow(ptrand,0.275362);
                                        		if (Temp>0.4) Temp=0.4;
                                        		if (Temp<0.19) Temp=0.19;
                                                }
                                                else {
							spe=1;
							Temp=0.33*pow(ptrand,0.3);
                                                        if (Temp>0.4) Temp=0.4;
                                                        if (Temp<0.15) Temp=0.15;
                                                }
						phirand=2.*3.141592*rando(&ir);
						raprand=maxrap*(-1.+2.*rando(&ir));
						pxrand=ptrand*cos(phirand);
						pyrand=ptrand*sin(phirand);
						mtrand=sqrt(pow(ptrand,2.)+pow(bomm[you],2.));
						phidif=acos((lostx[wa]*pxrand+losty[wa]*pyrand)/(ptlost*ptrand));
						cooper=exp(-mtrand/Temp*cosh(raprand))*mtrand/pow(Temp,5.)*cosh(raprand)*(ptrand*aa*cos(phidif)+mtrand*cosh(raprand))/normcoop[spe];
						if (cooper!=cooper) cooper=0.;
						if (fabs(cooper)>maxcooper[spe]) maxcooper[spe]=fabs(cooper);
						if (fabs(cooper)>1.) {
							//cout << " coooooOOOOoooOooOOper= " << cooper << " at ptrand= " << ptrand << "\n";
							normcoop[spe]*=(fabs(cooper)+0.0001);
							cooper=0.;
						}
						randian=rando(&ir);
					} while (pow(-1.,boml[you])*cooper<randian);
					//} while (fabs(cooper)<randian);
					pass=exp(-pow(difx/msigma,2.)/2.)*exp(-pow(dify/msigma,2.)/2.)*exp(-pow(difz/msigma,2.)/2.)*exp(-pow(dife/msigma,2.)/2.);
					sign=0.;
					if (cooper>0.) sign=1.;
					else sign=-1.;
					ppxback=pxback-pow(-1.,boml[you])*bomx[you]+sign*pxrand;
					ppyback=pyback-pow(-1.,boml[you])*bomy[you]+sign*pyrand;
					//peneback=eneback-pow(-1.,boml[you])*bome[you]+sign*sqrt(pow(ptrand,2.)+pow(masstra[spe],2.))*cosh(raprand);
                                        //ppzback=pzback-pow(-1.,boml[you])*bomz[you]+sign*sqrt(pow(ptrand,2.)+pow(masstra[spe],2.))*sinh(raprand);
					ppzback=pzback-pow(-1.,boml[you])*bomz[you]+sign*ptrand*sinh(raprand+raplost);
					peneback=eneback-pow(-1.,boml[you])*bome[you]+sign*sqrt(pow(ptrand*cosh(raprand+raplost),2.)+pow(masstra[spe],2.));
					difx=fabs(lostx[wa]-ppxback);
					dify=fabs(losty[wa]-ppyback);
					difz=fabs(lostz[wa]-ppzback);
					dife=fabs(loste[wa]-peneback);
					newpass=exp(-pow(difx/msigma,2.)/2.)*exp(-pow(dify/msigma,2.)/2.)*exp(-pow(difz/msigma,2.)/2.)*exp(-pow(dife/msigma,2.)/2.);
					if (newpass>pass) {
						bomx[you]=pxrand;
						bomy[you]=pyrand;
						//bome[you]=sqrt(pow(ptrand,2.)+pow(masstra[spe],2.))*cosh(raprand);
                                                //bomz[you]=bome[you]*tanh(raprand);
						bomz[you]=ptrand*sinh(raprand+raplost);
						bome[you]=sqrt(pow(ptrand*cosh(raprand+raplost),2.)+pow(masstra[spe],2.));
						if (cooper>0.) boml[you]=0.;
						else boml[you]=1.;
						pass=newpass;
						pxback=ppxback;
						pyback=ppyback;
						pzback=ppzback;
						eneback=peneback;
						runi+=1;
						encallao=0;
						//msigma=basesig;
						//if (count>3260) cout << " 1 Accepted= " << runi << "\n";
					}
					else {
						/*passrand=rando(&ir);
						if (newpass/pass > passrand) {
							bomx[you]=pxrand;
							bomy[you]=pyrand;
							//bome[you]=sqrt(pow(ptrand,2.)+pow(masstra[spe],2.))*cosh(raprand);
                                                	//bomz[you]=bome[you]*tanh(raprand);
							bomz[you]=ptrand*sinh(raprand+raplost);
							bome[you]=sqrt(pow(ptrand*cosh(raprand+raplost),2.)+pow(masstra[spe],2.));;
							if (cooper>0.) boml[you]=0.;
							else boml[you]=1.;
							pass=newpass;
							pxback=ppxback;
							pyback=ppyback;
							pzback=ppzback;
							eneback=peneback;
							runi+=1;
							encallao=0;
							//msigma=basesig;
							//if (count>3260) cout << " 2 Accepted= " << runi << "\n";
						}
						else */
						encallao+=1;
					}
					difx=fabs(lostx[wa]-pxback);
					dify=fabs(losty[wa]-pyback);
					difz=fabs(lostz[wa]-pzback);
					dife=fabs(loste[wa]-eneback);
					
					if (encallao>50000) {
						numenc+=1;
						if (numenc>5) {
							toomuch+=1;
							break;
						}
						//cout << " ENCALLA DPx= " << difx << " DPy= " << dify << " DPz= " << difz << " DE= " << dife << "\n";
						//cout << " runi encalla= " << runi << "\n";
						goto thisis;
					}
					clock_t endClock = clock();
					if (double((endClock - startClock)) / CLOCKS_PER_SEC>clocklim) {
						//cout << " CLOCK !!! \n";
						clocklim+=0.02;
						tooclock+=1;

                                          double remass2=dife*dife-difx*difx-dify*dify-difz*difz;

                                          //cout << " Got to= " << " DifX= " << difx << " DifY= " << dify << " DifZ= " << difz << " DifE= " << dife << endl;
                                          if (fabs(difx)<3.*tole && fabs(dify)<3.*tole && fabs(difz)<8.*tole && remass2>0.) {
                                            //Add extra particle with remnant momentum
                                            //cout << " OUTSIDE TOLERANCE \n";
                                            //cout << " remnant mass= " << remass2 << " " << sqrt(remass2) << endl;

  	 	 			    difx=(lostx[wa]-pxback);
				 	    dify=(losty[wa]-pyback);
					    difz=(lostz[wa]-pzback);
					    dife=(loste[wa]-eneback);
                                            
					    if (fabs(remass2-masstra[0]*masstra[0])<fabs(remass2-masstra[1]*masstra[1])) spe=0;
                                            else spe=1;
                                            //cout << " spe= " << spe << endl;

                                            //if (nr.rando()<=0.05) spe=1;    //is proton
                                            //else spe=0;                     //is pion

                                            double stat;
                                            if (dife<0.) stat=-1., boml[kk]=1;
                                            else stat=1., boml[kk]=0;

                                            vector<double> p;
                                            double tdife=sqrt(difx*difx+dify*dify+difz*difz+masstra[spe]*masstra[spe]);
                                            p.push_back(stat*difx), p.push_back(stat*dify), p.push_back(stat*difz), p.push_back(tdife);
         				    bomx[kk]=p[0], bomy[kk]=p[1], bomz[kk]=p[2], bome[kk]=p[3];                    
					    
					    pxback+=pow(-1.,boml[kk])*bomx[kk];
					    pyback+=pow(-1.,boml[kk])*bomy[kk];
					    pzback+=pow(-1.,boml[kk])*bomz[kk];
					    eneback+=pow(-1.,boml[kk])*bome[kk];
					    
  	 	 			    difx=fabs(lostx[wa]-pxback);
				 	    dify=fabs(losty[wa]-pyback);
					    difz=fabs(lostz[wa]-pzback);
					    dife=fabs(loste[wa]-eneback);
					    kk+=1;
                                            break;
                                          }

					  if (tooclock>100) break;
					  goto thisis;
	                                }
					//cout << " Time lapsed= " << (endClock - startClock) / CLOCKS_PER_SEC << " with RapLost= " << raplost << "\n";
					cooper=0.;
					randian=1.;
				} while (runi<Nrun && (difx>tole || dify>tole || difz>tole || dife>tole));
				//cout << " DPx= " << difx << " DPy= " << dify << " DPz= " << difz << " DE= " << dife << " runi= " << runi  << "\n";
				
				if (difx!=0.) {

  	 	 			    difx=(lostx[wa]-pxback);
				 	    dify=(losty[wa]-pyback);
					    difz=(lostz[wa]-pzback);
					    dife=(loste[wa]-eneback);
                          	  
				  double remass2=dife*dife-difx*difx-dify*dify-difz*difz;

                                  if (fabs(remass2-masstra[0]*masstra[0])<fabs(remass2-masstra[1]*masstra[1])) spe=0;
                                  else spe=1;

                                  double stat;
                                  if (dife<0.) stat=-1., boml[kk]=1;
                                  else stat=1., boml[kk]=0;
                                           
				  vector<double> p;
                                  double tdife=sqrt(difx*difx+dify*dify+difz*difz+masstra[spe]*masstra[spe]);
                                  p.push_back(stat*difx), p.push_back(stat*dify), p.push_back(stat*difz), p.push_back(tdife);
         			  bomx[kk]=p[0], bomy[kk]=p[1], bomz[kk]=p[2], bome[kk]=p[3];                    
					    
				  pxback+=pow(-1.,boml[kk])*bomx[kk];
			          pyback+=pow(-1.,boml[kk])*bomy[kk];
				  pzback+=pow(-1.,boml[kk])*bomz[kk];
				  eneback+=pow(-1.,boml[kk])*bome[kk];
					    
  	 	 	          difx=fabs(lostx[wa]-pxback);
				  dify=fabs(losty[wa]-pyback);
				  difz=fabs(lostz[wa]-pzback);
				  dife=fabs(loste[wa]-eneback);
					
				  kk+=1;

                                }

				//cout << difx << " " << dify << " " << difz << " " << dife << endl;
			        Npart=0;
				for (int c=0; c<kk; c++) {
                                        ptrand=sqrt(pow(bomx[c],2.)+pow(bomy[c],2.));
                                        mtrand=sqrt(pow(ptrand,2.)+pow(bomm[c],2.));
					phirand=acos((bomy[c]*losty[wa]+bomx[c]*lostx[wa])/ptrand/ptlost);
                                        if (phirand<0.) phirand+=2.*3.141592;
                                        //raprand=asinh(bomz[c]/mtrand);
					raprand=asinh(bomz[c]/ptrand);
                                        int nphi=int(phirand/(2.*3.141592)*200.);
                                        int nrap=0;
                                        raprand-=raplost;
                                        if (raprand>0.) nrap=100+int((raprand+rbin/2.)/rbin);
                                        if (raprand<0.) nrap=100+int((raprand-rbin/2.)/rbin);
                                        int npt=int(200.*(ptrand/sqrt(maxptsq)));
					spe=0;
					if (bomm[c]>0.5) spe=1;
                                        if (npt<200) ptisto[npt][spe]+=pow(-1.,boml[c]);
                                        if (nrap<201 && nrap>-1) risto[nrap][spe]+=pow(-1.,boml[c]);
                                        if (nphi<201) phisto[nphi][spe]+=pow(-1.,boml[c]);
                                        if (nrap<201 && nphi<201 && nrap>-1) {
						phirap[nphi][nrap][spe]+=pow(-1.,boml[c]);
					}
					Yield[spe]+=int(pow(-1.,boml[c]));
					Npart+=int(pow(-1.,boml[c]));
					savedback << bomx[c] << " " << bomy[c] << " " << bomz[c] << " " << bome[c] << " " << int(pow(-1.,boml[c])) << "\n";
				}
				if (lostx[wa]-pxback>0.) dx=10+int(10.*(lostx[wa]-pxback+0.25)/5.);
				else dx=10+int(10.*(lostx[wa]-pxback-0.25)/5.);
				if (losty[wa]-pyback>0.) dy=10+int(10.*(losty[wa]-pyback+0.25)/5.);
				else dy=10+int(10.*(losty[wa]-pyback-0.25)/5.);
				if (lostz[wa]-pzback>0.) dz=10+int(10.*(lostz[wa]-pzback+0.25)/5.);
				else dz=10+int(10.*(lostz[wa]-pzback-0.25)/5.);
				if (loste[wa]-eneback>0.) de=10+int(10.*(loste[wa]-eneback+0.25)/5.);
				else de=10+int(10.*(loste[wa]-eneback-0.25)/5.);
				if (dx<0) dx=0;
				if (dx>19) dx=19;
				if (dy<0) dy=0;
				if (dy>19) dy=19;
				if (dz<0) dz=0;
				if (dz>19) dz=19;
				if (de<0) de=0;
				if (de>19) de=19;
				passdist[0][dx][0]+=1.;
				passdist[1][dy][0]+=1.;
				passdist[2][dz][0]+=1.;
				passdist[3][de][0]+=1.;

				nmed+=Npart;
				ntwomed+=pow(Npart,2.);
				double aveN=mtlost/3./0.35;
                                if (fabs(Npart-aveN)<50.) {
                                        nisto[49+int(Npart-aveN)]+=1.;
                                }
			}
		}
		if (k!=0) savedback << " 0.000123 " << " 0. " << " 0. " << " 0. " << count << "\n";

		//cout << " NEXT EVENT ";
		count+=1;
		//cout << " " << count << " \n";
	//End Event Loop
	} while (count < Ntot);
	totalsave=Ntot;
	//cout << " MaxCooper= " << maxcooper << "\n";

	//Dists Histo
#ifdef DO_PRINT
	for (unsigned a=0; a<200; a++) {
		for (unsigned b=0; b<200; b++) {
			raphi << 2.*3.141592/400.+double(a)/200.*2.*3.141592 << " " << -maxrap+double(b)*rbin << " " << phirap[a][b][0] << " " << phirap[a][b][1] <<  "\n";
		}
	}
	for (unsigned a=0; a<200; a++) {
                phifile << 2.*3.141592/400.+double(a)/200.*2.*3.141592 << " " << phisto[a][0]/double(totalsave) << " " << sqrt(phisto[a][0])/double(totalsave) << " " << phisto[a][1]/double(totalsave) << " " << sqrt(phisto[a][1])/double(totalsave) << "\n";
        }
        for (unsigned a=0; a<200; a++) {
                rapfile << -maxrap+double(a)*rbin << " " << risto[a][0]/double(totalsave) << " " << sqrt(risto[a][0])/double(totalsave) << " " << risto[a][1]/double(totalsave) << sqrt(risto[a][1])/double(totalsave) << "\n";
        }
        for (unsigned a=0; a<200; a++) {
                ptfile << double(a)/200.*sqrt(maxptsq)+sqrt(maxptsq)/400. << " " << ptisto[a][0]/double(totalsave) << " " << sqrt(ptisto[a][0])/double(totalsave) << " " << ptisto[a][1]/double(totalsave) << " " << sqrt(ptisto[a][1])/double(totalsave) << "\n";
        }
        for (unsigned a=0; a<100; a++) {
                nfile << -49+double(a) << " " << nisto[a]/double(count) << " " << sqrt(nisto[a])/count << "\n";
        }
        double nerr=sqrt(ntwomed+double(count)*pow(nmed/double(count),2.)-2.*pow(nmed,2.)/double(count))/sqrt(double(count)*(double(count)-1.));
        cout << " <N>= " << nmed/double(count) << " " << nerr << "\n";
	cout << " <Pion Yield>= " << Yield[0]/double(totalsave) << " <Proton Yield>= " << Yield[1]/double(totalsave) << " <Proton/Pion Yield> " << Yield[1]/Yield[0] << "\n";

	raphi.close();
        phifile.close();
        rapfile.close();
        ptfile.close();
        nfile.close();

	//EM Conservation Histo
	for (unsigned b=0; b<20; b++) {
                cons << -5.+double(b)*0.5 << " " << passdist[0][b][0] << " " << passdist[1][b][0] << " " << passdist[2][b][0] << " " << passdist[3][b][0] << "\n";
        }
	cons.close();

#endif

//End Main
}

//--------------Random Number-----------------------//
double rando (int *ir)
{
	double da=16807.;
	double db=2147483647.;
	double dc=2147483648.;
	double usran;

	*ir=(int)(fabs(fmod(da*(*ir),db)+0.5));
	usran=double(*ir)/dc;
	return usran;
}

//Factorial
double factorial(int n)
{
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}
