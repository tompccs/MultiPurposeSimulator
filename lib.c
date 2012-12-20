/*
    (c) Tom Robbins 2012

*/

#include "lib.h"
#define USE_VAR_POOL
//#define PRINT_KVALS

void iterate_to_file(int(*iter_func)(double *, double *, double), int variable_count, double * starting_values, double independent_variable_step, char ** variable_labels, FILE * fout){
    if(fout == NULL){
        return;
    }

    int i;
    for(i=0;i<variable_count;i++){
        fprintf(fout, "%s,", variable_labels[i]);
    }
    fprintf(fout, "\n");

    // copy the starting values in case they need to be used elsewhere
    double * variables = malloc(sizeof(double) * variable_count);
    for(i=0;i<variable_count;i++){
        variables[i] = starting_values[i];
    }

    double * previous = malloc(sizeof(double) * variable_count);
    double * next = variables;

    // iter_func is called with the following parameters:
    // iter_func(double * in_variables, double * out_variables, double step)
    do{
        for(i=0;i<variable_count;i++){
            fprintf(fout, "%lf,", next[i]);
        }

        // swap the pointers (ah yes, the old switcheroo)
        double * temp = previous;
        previous = next;
        next = temp;

        fprintf(fout, "\n");
    }while((iter_func(previous, next, independent_variable_step) == CONTINUE_ITERATING));
}

void runge_kutta_4th(double(*func)(double *, int), double * vars_in, double * vars_out, int var_count, int const_count, double step){
    /*
    The elements in the arrays vars_in and vars_out represent the following:
    0               --> independent variable
    1               --> dependent variables
    .                   ""
    .                   ""
    (var_count - 1)     ""
    var_count       --> constants
    .                   ""
    .                   ""
    (const_count -1)    ""
    
    */

    #ifdef PRINT_KVALS
    char * labels[] = {"x", "y", "xvel", "yvel"};
    #endif

    // set up variables using the variable pool:
    #ifdef USE_VAR_POOL
    double * temp = variable_pool;
    double * k[] = {&(variable_pool[const_count + var_count]), &(variable_pool[const_count + var_count + (var_count - 1)]),
         &(variable_pool[const_count + var_count + (var_count - 1) * 2]), &(variable_pool[const_count + var_count + (var_count - 1) * 3])};
    #else
    double * temp = malloc(sizeof(double) * (const_count + var_count));
    double * k[] = { malloc(sizeof(double) * (var_count - 1)), malloc(sizeof(double) * (var_count - 1)), malloc(sizeof(double) * (var_count - 1)),
        malloc(sizeof(double) * (var_count - 1)) };
    #endif

    // copy the constants first:
    int i, j;
    for(i=var_count;i<const_count + var_count;i++){
        vars_out[i] = vars_in[i];
        temp[i] = vars_in[i];
    }

    // initialise the k-array. Access like k[k-index - 1][variable]
    // for each variable, calculate the k values using the corresponding function:

    // k1
    for(i=1;i<var_count;i++){
        // set the k-values:
        k[0][i-1] = func(vars_in, i);
    }

    // k2
    temp[0] = vars_in[0] + step/2;
    for(j=1;j<var_count;j++){
        temp[j] = vars_in[j] + step * k[0][j-1] / 2;
    }

    for(i=1;i<var_count;i++){
        k[1][i-1] = func(temp, i);
    }

    // k3
    for(j=1;j<var_count;j++){
        temp[j] = vars_in[j] + step * k[1][j-1] / 2;
    }

    for(i=1;i<var_count;i++){
        k[2][i-1] = func(temp, i);
    }

    // k4 + calculate new variables
    temp[0] = vars_in[0] + step;
    for(j=1;j<var_count;j++){
        temp[j] = vars_in[j] + step * k[2][j-1];
    }
    
    #ifdef PRINT_KVALS
    printf("iteration at t=%lf:\n", vars_in[0]);
    #endif
    for(i=1;i<var_count;i++){
        k[3][i-1] = func(temp, i);
        vars_out[i] = vars_in[i] + step/6*(k[0][i-1] + 2*k[1][i-1] + 2*k[2][i-1] + k[3][i-1]);
        #ifdef PRINT_KVALS
        for(j=0;j<4;j++){
            printf("     k%d%s = %lf\n", j + 1, labels[i-1], k[j][i-1]);
        }
        printf("     %s = %lf\n", labels[i-1], vars_out[i]);
        #endif
    }
    #ifdef PRINT_KVALS
    printf("\n");
    #endif

    // always assume the 0th variable is the independent one.
    vars_out[0] = vars_in[0] + step;

    #ifndef USE_VAR_POOL
    free(temp);
    for(i=0;i<4;i++){
        free(k[i]);
    }
    #endif

}

void set_up_runge_kutta_4th(int variable_count, int constant_count){
    // this is an optimisation so that only one malloc call needs to be
    // made per simulation.

    // for each variable, we need 4 k-values, as well as a temporary copy.
    variable_pool = malloc(sizeof(double) * (4 * (variable_count - 1) + constant_count + variable_count));
}

void free_runge_kutta_4th(){
    // you must always call this function between simulations!!!
    free(variable_pool);
}

int process_flags(int argc, char ** args, int flagnc, char ** flagns){
    int flags = 0, i, j;
    for(i=0; i<argc; i++) {
        for(j=0; j<flagnc; j++){
            if(!strcmp(flagns[j], args[i])){
                flags += (1 << j);
                break;
            }
        }
    }
    return flags;
}

int process_numeric_args(int argc, char ** args, double * processed_args){
    int i, j = 0;
    for(i=0;i<argc;i++){
        double number;
        if (sscanf(args[i], "%lf", &number) >= 1){
            processed_args[j] = number;
            j ++;
        }
    }
    return j;
}

