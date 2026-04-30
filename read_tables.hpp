#ifndef read_tables_H
#define read_tables_H

#include <iostream>
#include <utility>
#include <fstream>
#include <sstream>
//#include <gsl/gsl_math.h>
//#include <gsl/gsl_interp2d.h>
//#include <gsl/gsl_spline2d.h>

//const gsl_interp2d_type *TmX = gsl_interp2d_bilinear; //SO MUCH FASTER THAN BICUBIC
//const gsl_interp2d_type *TmX = gsl_interp2d_bicubic;
//gsl_interp_accel *xacc;
//gsl_interp_accel *yacc;
//gsl_spline2d *spline;

extern const int SIZE_PIN;
extern const int SIZE_X;

extern double ****q_mx_table;
extern double ****g_mx_table;
extern double pin_vals[];
extern double pin_step, pin_min, pin_max;

extern double x_vals[];
extern double step_x;

extern bool use_tables;

void read_tables(std::string tables_path);

#endif
