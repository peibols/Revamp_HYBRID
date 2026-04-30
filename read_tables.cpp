#include <iostream>
#include <utility>
#include <tuple>
#include <fstream>
#include <sstream>
#include <gsl/gsl_math.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "read_tables.hpp"

//const gsl_interp2d_type *TmX = gsl_interp2d_bilinear; //SO MUCH FASTER THAN BICUBIC
//const gsl_interp2d_type *TmX = gsl_interp2d_bicubic;
//gsl_interp_accel *xacc;
//gsl_interp_accel *yacc;
//gsl_spline2d *spline;

extern const int SIZE_PIN = 501;
extern const int SIZE_X = 101;

double ****q_mx_table;
double ****g_mx_table;
double pin_vals[SIZE_PIN];
double pin_step, pin_min, pin_max;

double x_vals[SIZE_X];
double step_x;

bool use_tables;

void read_tables(std::string tables_path) {

  //Allocate tables quark
  int nX=7, nA=68, nP=SIZE_PIN, nY=SIZE_X;
  q_mx_table = (double ****)malloc(nX * sizeof(double ***));
  if (q_mx_table==NULL)
  {
    fprintf(stderr, "out of memory\n");
    exit(0);
  }
  for (int i=0; i<nX; i++)
  {
    q_mx_table[i]=(double ***)malloc(nA * sizeof(double **));
    if (q_mx_table[i]==NULL)
    {
      fprintf(stderr, "out of memory\n");
      exit(0);
    }
  }
  for (int i=0; i<nX; i++)
  {
    for (int j=0; j<nA; j++)
    {
      q_mx_table[i][j]=(double **)malloc(nP * sizeof(double *));
      if (q_mx_table[i][j]==NULL)
      {
        fprintf(stderr, "out of memory\n");
        exit(0);
      }
    }
  }
  for (int i=0; i<nX; i++)
  {
    for (int j=0; j<nA; j++)
    {
      for (int k=0; k<nP; k++)
      {
        q_mx_table[i][j][k]=(double *)malloc(nY * sizeof(double));
        if (q_mx_table[i][j][k]==NULL)
        {
          fprintf(stderr, "out of memory\n");
          exit(0);
        }
      }
    }
  }
  
  //Allocate tables gluon
  g_mx_table = (double ****)malloc(nX * sizeof(double ***));
  if (g_mx_table==NULL)
  {
    fprintf(stderr, "out of memory\n");
    exit(0);
  }
  for (int i=0; i<nX; i++)
  {
    g_mx_table[i]=(double ***)malloc(nA * sizeof(double **));
    if (g_mx_table[i]==NULL)
    {
      fprintf(stderr, "out of memory\n");
      exit(0);
    }
  }
  for (int i=0; i<nX; i++)
  {
    for (int j=0; j<nA; j++)
    {
      g_mx_table[i][j]=(double **)malloc(nP * sizeof(double *));
      if (g_mx_table[i][j]==NULL)
      {
        fprintf(stderr, "out of memory\n");
        exit(0);
      }
    }
  }
  for (int i=0; i<nX; i++)
  {
    for (int j=0; j<nA; j++)
    {
      for (int k=0; k<nP; k++)
      {
        g_mx_table[i][j][k]=(double *)malloc(nY * sizeof(double));
        if (g_mx_table[i][j][k]==NULL)
        {
          fprintf(stderr, "out of memory\n");
          exit(0);
        }
      }
    }
  }

  //Read quark tables
  int nbins_pin, nbins_x;
  for (int iX=0; iX<nX; iX++) {
  //for (int iX=0; iX<1; iX++) { //DEBUG
    int n = 1;
    for (int iA=0; iA<nA; iA++) {
      int res = (iA % 4);
      int g, d;
      if (res == 0) g = -1, d = -1;
      else if (res == 1) g = -1, d = 1;
      else if (res == 2) g = 1, d = -1;
      else if (res == 3) g = 1, d = 1;
      else {
        std::cout << "WTF res= " << res << std::endl;
	exit(1);
      }
      std::ostringstream fsq;
      fsq << tables_path.c_str() << "quark_tables/m" << iX+1 << "x_g_" << int(g) << "_d_" << int(d) << "_n_" << int(n) << ".dat";
      std::ifstream qinfile(fsq.str().c_str());
      if (qinfile.fail()) {
        std::cout << "No file = " << fsq.str().c_str() << std::endl;
	exit(1);
      }
      std::ostringstream fsg;
      fsg << tables_path.c_str() << "gluon_tables/m" << iX+1 << "x_g_" << int(g) << "_d_" << int(d) << "_n_" << int(n) << ".dat";
      std::ifstream ginfile(fsg.str().c_str());
      if (ginfile.fail()) {
        std::cout << "No file = " << fsg.str().c_str() << std::endl;
	exit(1);
      }

      std::cout << "Reading file= " << fsq.str() << " and " << fsg.str() << std::endl;
      std::string sq,sg;
      getline(qinfile,sq); //Header
      getline(ginfile,sg); //Header
      if (iX==0 && iA==0) {
        std::istringstream ifq(sq);
        std::istringstream ifg(sg);
	std::string scrap;
        ifq >> scrap >> pin_min >> scrap >> pin_max >> scrap >> nbins_pin >> scrap >> nbins_x;
        ifg >> scrap >> pin_min >> scrap >> pin_max >> scrap >> nbins_pin >> scrap >> nbins_x;
      }

      for (int iY=0; iY<nY; iY++) {
        getline(qinfile,sq);
        std::istringstream ifq(sq);
        getline(ginfile,sg);
        std::istringstream ifg(sg);
        for (int iP=0; iP<nP; iP++) ifq >> q_mx_table[iX][iA][iP][iY];
        for (int iP=0; iP<nP; iP++) ifg >> g_mx_table[iX][iA][iP][iY];
      }
      qinfile.close();
      ginfile.close();

      //break; //DEBUG
      if (res==3) n++;

    }
    
    //break; //DEBUG
  }
  std::cout << "Finished reading mX tables" << std::endl;

  //spline = gsl_spline2d_alloc(TmX, nP, nY);
  //xacc = gsl_interp_accel_alloc();
  //yacc = gsl_interp_accel_alloc();

  pin_step = (pin_max - pin_min)/double(nbins_pin);
  for (int iP=0; iP<nP; iP++) pin_vals[iP] = pin_min + pin_step*double(iP);
  //for (int iP=0; iP<nP; iP++) std::cout << " pin val iP= " << iP << " = " << pin_vals[iP] << std::endl; 

  step_x = 1./double(nbins_x);
  for (int iY=0; iY<nY; iY++) x_vals[iY] = step_x*double(iY);
  //for (int iY=0; iY<nY; iY++) std::cout << " x_vals iY= " << iY << " = " << x_vals[iY] << std::endl;

  return;

}
