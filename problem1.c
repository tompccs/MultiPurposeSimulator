/*
    (c) Tom Robbins 2012

*/

// #define VERBOSE_DEBUG

#define FLAG_ARRAY {"--orbit", "--simple", "--free", "--2D", "--3D", "--help", "--stdout"}
#define FLAG_ARRAY_SIZE 7

#define FLAG_ORBIT 1
#define FLAG_SIMPLE 2
#define FLAG_FREE 4
#define FLAG_2D 8
#define FLAG_3D 16
#define FLAG_HELP 32
#define FLAG_STDOUT 64

#define GRAVITATIONAL_CONSTANT 6.673E-11

#include "lib.h"

void help(){
    printf("Usage:\n");
    printf("    simulator [filepath] --flag1 --flag2... var1 var2 var3...\n\n");

    printf("Available options:\n");
    printf("    --orbit\n");
    printf("        Simulates an orbit.\n");
    printf("    --simple\n");
    printf("        In the orbit case, simulates one static body being \n        orbited by another.\n");
    printf("    --free\n");
    printf("        In the orbit case, simulates two or more bodies.\n");
    printf("    --2D\n");
    printf("        Simulates in 2 dimensions.\n");
    printf("    --3D\n");
    printf("        Simulates in 3 dimensions. Results in undefined \n        behaviour where --2D has also been set.\n");
    printf("    --stdout\n");
    printf("        Writes to the standard out (ie, the command line) \n        rather than a physical file. In this case a filepath \n        does not need to be set. \n        Useful for visulaising simulations 'live'.\n");
    printf("    --help\n");
    printf("        Displays this message.\n\n");

    printf("Example:\n");
    printf("To run a simulation of a simple 2D orbit with the following \nparameters and write the file './out.csv': \n");
    printf("    start time = 0\n");
    printf("    start xpos = 10E3\n");
    printf("    start ypos = 0\n");
    printf("    start xvel = 0\n");
    printf("    start yvel = 1000\n");
    printf("    object mass = 1E9\n");
    printf("    time limit = 1000\n");
    printf("    time step = 0.1\n");
    printf("you would run the command:\n");
    printf("    simulator ./out.csv --orbit --simple --2D 0 10E3 0 0 1000 1E9 1000 0.1\n");
    printf("Where there are multiple bodies, the numeric arguments are entered \nas follows (where there are N variables and M constants per body):\n");
    printf("    <body1var1>...<body1varN> <body1const1>...<body1constM>\n");
}

double simple_2d_orbit_functions(double * vars_in, int function_ref){
    /*
    Our variables are as follows:
    0 --> time
    1 --> x-position of satellite
    2 --> y-position of satellite
    3 --> x-velocity of satellite
    4 --> y-velocity of satellite
    CONSTANTS:
    5 --> mass of object
    6 --> time limit
    USEFUL THINGS TO PRINT:
    7 --> total energy / unit mass of satellite
    */

    // Orbit function with one static object being orbited by another in 2 dimensions
    switch(function_ref){
    case 1:
        return vars_in[3];
    case 2:
        return vars_in[4];
    case 3:
        return - GRAVITATIONAL_CONSTANT * vars_in[5]*vars_in[1]/pow(pow(vars_in[1],2)+pow(vars_in[2],2), 1.5);
    case 4:
        return - GRAVITATIONAL_CONSTANT * vars_in[5]*vars_in[2]/pow(pow(vars_in[1],2)+pow(vars_in[2],2), 1.5);
    default:
        return vars_in[function_ref];
    }
}

int simple_2d_orbit_runge_kutta_4th(double * vars_in, double * vars_out, double step){
    runge_kutta_4th(&simple_2d_orbit_functions, vars_in, vars_out, 5, 2, step);

    // we also want to know the total energy per unit mass of the satellite. So whack that into a variable.
    vars_out[7] = .5 * (pow(vars_out[3], 2) + pow(vars_out[4], 2)) - 
        GRAVITATIONAL_CONSTANT * vars_out[5] / sqrt(pow(vars_out[1],2) + pow(vars_out[2],2));

    // check the terminating condition. In this case, a time limit.
    if(vars_out[0] >= vars_in[6]){
        return STOP_ITERATING;
    }else{
        return CONTINUE_ITERATING;
    }
}

static int body_count = 1;

double free_2d_orbit_functions(double * vars_in, int function_ref){
    /*
    Our variables and constants are as follows:
    0 --> time
    1 --> xpos of body 0
    2 --> ypos of body 0
    3 --> xvel of body 0
    4 --> yvel of body 0
    5 --> mass of body 0
    .
    .
    (body count * 5 + 1)
      --> time limit
    */

    if (function_ref > body_count * 5){
        return 0;
    }

    int i, actual_function_ref = function_ref % 5;
    int a = function_ref - actual_function_ref; // array offset for this body.
    int body_index = a / 5;

    #ifdef VERBOSE_DEBUG
    printf("fr = %d, afr = %d, a = %d, bi = %d\n", function_ref, actual_function_ref, a, body_index);
    #endif

    double xdiff, ydiff, total_acceleration;
    switch(actual_function_ref){
    case 1:
        #ifdef VERBOSE_DEBUG
        printf("body %d x-velocity = %lf (index %d)\n", body_index, vars_in[a + 3], a+3);
        #endif
        return vars_in[a + 3];
    case 2:
        #ifdef VERBOSE_DEBUG
        printf("body %d y-velocity = %lf (index %d)\n", body_index, vars_in[a + 4], a+4);
        #endif
        return vars_in[a + 4];
    case 3:
        // x-acceleration
        total_acceleration = 0;
        for(i=0;i<body_count;i++){
            if(i != body_index){
                xdiff = vars_in[a + 1] - vars_in[i*5 + 1];
                ydiff = vars_in[a + 2] - vars_in[i*5 + 2];

                // check for collision
                if(xdiff <= DBL_EPSILON || ydiff <= DBL_EPSILON){
                    return 0;
                }

                #ifdef VERBOSE_DEBUG
                printf("mass of body %d = %lf\n", i, vars_in[i * 5 + 5]);
                printf("body %d-->%d xdiff = %lf, ydiff = %lf\n", body_index, i, xdiff, ydiff);
                #endif
                total_acceleration -=
                    GRAVITATIONAL_CONSTANT * vars_in[i * 5 + 5] * xdiff / pow(pow(xdiff, 2) + pow(ydiff, 2), 1.5);
            }
        }
        #ifdef VERBOSE_DEBUG
        printf("body %d x-acceleration = %lf\n", body_index, total_acceleration);
        #endif
        return total_acceleration;
    case 4:
        // y-acceleration
        total_acceleration = 0;
        for(i=0;i<body_count;i++){
            if(i != body_index){
                xdiff = vars_in[a + 1] - vars_in[i*5 + 1];
                ydiff = vars_in[a + 2] - vars_in[i*5 + 2];

                // check for collision
                if(xdiff <= DBL_EPSILON || ydiff <= DBL_EPSILON){
                    return 0;
                }

                #ifdef VERBOSE_DEBUG
                printf("body %d-->%d xdiff = %lf, ydiff = %lf\n", body_index, i, xdiff, ydiff);
                #endif
                total_acceleration -=
                    GRAVITATIONAL_CONSTANT * vars_in[i * 5 + 5] * ydiff / pow(pow(xdiff, 2) + pow(ydiff, 2), 1.5);
            }
        }
        #ifdef VERBOSE_DEBUG
        printf("body %d y-acceleration = %lf\n", body_index, total_acceleration);
        #endif
        return total_acceleration;
    default:
        return 0;
    }

}

int free_2d_orbit_runge_kutta_4th(double * vars_in, double * vars_out, double step){
    // this is a 2D simulation, so the number of variables is 5 * body_count + 1
    // (4 variables: xpos, ypos, xvel, yvel and 1 constant: mass of body)
    runge_kutta_4th(&free_2d_orbit_functions, vars_in, vars_out, 5 * body_count + 1, 1, step);

    // check the terminating condition. In this case, a time limit.
    if(vars_out[0] >= vars_in[5 * body_count + 1]){
        return STOP_ITERATING;
    }else{
        return CONTINUE_ITERATING;
    }
}

int main(int argc, char ** args){
    char *flag_array[] = FLAG_ARRAY;
    int flags = process_flags(argc, args, FLAG_ARRAY_SIZE, flag_array);
    double * numeric_args = malloc(sizeof(double) * (argc - 1));
    int numeric_arg_count = process_numeric_args(argc, args, numeric_args);

    if(flags & FLAG_HELP){
        help();
        return 0;
    }

    // first argument should always be a file unless --stdout
    FILE * fout;
    if(flags & FLAG_STDOUT){
        // write to stdout. This is useful if we want to pipe the data somewhere (ie, for visualisation)
        fout = stdout;
    }else{
        fout = fopen(args[1], "w");
        if(fout == NULL){
            printf("Could not open file at %s for writing. No such directory or permission denied.\n", args[1]);
            return 1;
        }
    }

    if(flags & FLAG_ORBIT){
        if(flags & FLAG_SIMPLE){
            if(flags & FLAG_2D){
                if(numeric_arg_count == 8){
                    char *labels[8] = {"time", "xpos", "ypos", "xvel", "yvel", "object_mass", "time_limit", "total_energy"};
                    set_up_runge_kutta_4th(5, 2);
                    iterate_to_file(&simple_2d_orbit_runge_kutta_4th, 8, numeric_args, numeric_args[7], labels, fout);
                    free_runge_kutta_4th();
                }else{
                    printf("Invalid number of numerical arguments. Need 8, %d given.\n", numeric_arg_count);
                }
            }else{
                printf("Simple 3D simulation has not been implemented.\n");
            }
        }else if(flags & FLAG_FREE){
            if(flags & FLAG_2D){
                // make sure a valid number of args has been entered
                if(numeric_arg_count >= 8 && (numeric_arg_count - 3) % 5 == 0){
                    body_count = (numeric_arg_count - 3) / 5;
                    char ** labels = malloc(sizeof(char *) * (body_count * 5 + 2));
                    labels[0] = "time";
                    int i;
                    for(i=0;i<body_count;i++){
                        labels[i * 5 + 1] = malloc(8);
                        sprintf(labels[i * 5 + 1], "%d.xpos", i);
                        labels[i * 5 + 2] = malloc(8);
                        sprintf(labels[i * 5 + 2], "%d.ypos", i);
                        labels[i * 5 + 3] = malloc(8);
                        sprintf(labels[i * 5 + 3], "%d.xvel", i);
                        labels[i * 5 + 4] = malloc(8);
                        sprintf(labels[i * 5 + 4], "%d.yvel", i);
                        labels[i * 5 + 5] = malloc(8);
                        sprintf(labels[i * 5 + 5], "%d.mass", i);
                    }
                    labels[body_count * 5 + 1] = "time_limit";
                    set_up_runge_kutta_4th(5 * body_count + 1, 1);
                    iterate_to_file(&free_2d_orbit_runge_kutta_4th, 5 * body_count + 2, numeric_args, numeric_args[numeric_arg_count - 1], labels, fout);
                    free_runge_kutta_4th();
                }else{
                    printf("Invalid number of numerical arguments. Need at least 8 with 5 arguments for each body. %d given\n", numeric_arg_count);
                    help();
                    return 1;
                }
            }else if(flags & FLAG_3D){
                
            }
        }else{
            printf("Please specify either --simple or --free\n");
            help();
            return 1;
        }
    }

    return 0;
}
