/*
	(c) Tom Robbins 2012

*/
#ifndef LIB_INCLUDED
#define LIB_INCLUDED

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define CONTINUE_ITERATING 0
#define STOP_ITERATING 1

#define FLAG_EULER 1
#define FLAG_FELIX 2
#define FLAG_CENTRAL 4
#define FLAG_VARY_DRAG 8

static double * variable_pool = NULL;

void iterate_to_file(int(*iter_func)(double *, double *, double), int variable_count, double * starting_values, double independent_variable_step, char ** variable_labels, FILE * fout);
int process_flags(int argc, char ** args, int flagnc, char ** flagns);
int process_numeric_args(int argc, char ** args, double * processed_args);
void runge_kutta_4th(double(*func)(double *, int), double * vars_in, double * vars_out, int var_count, int const_count, double step);
void set_up_runge_kutta_4th(int variable_count, int constant_count);
void free_runge_kutta_4th();

#endif
