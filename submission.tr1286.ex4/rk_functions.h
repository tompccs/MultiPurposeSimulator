/*
	(c) Tom Robbins 2012

*/

#import "lib.h"
#include <stdio.h>
#include <math.h>
#include <float.h>

#define GRAVITATIONAL_CONSTANT 6.673E-11
// #define VERBOSE_DEBUG

extern int body_count;

double simple_2d_orbit_functions(double * vars_in, int function_ref);
int simple_2d_orbit_runge_kutta_4th(double * vars_in, double * vars_out, double step);
double free_2d_orbit_functions(double * vars_in, int function_ref);
int free_2d_orbit_runge_kutta_4th(double * vars_in, double * vars_out, double step);
double free_3d_orbit_functions(double * vars_in, int function_ref);
int free_3d_orbit_runge_kutta_4th(double * vars_in, double * vars_out, double step);
