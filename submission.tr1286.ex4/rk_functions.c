/*
    (c) Tom Robbins 2012

*/


#include "rk_functions.h"

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
        if(fabs(vars_in[1]) <= DBL_EPSILON && fabs(vars_in[2]) <= DBL_EPSILON){
            return 0;
        }else{
            return vars_in[3];
        }
    case 2:
        if(fabs(vars_in[1]) <= DBL_EPSILON && fabs(vars_in[2]) <= DBL_EPSILON){
            return 0;
        }else{
            return vars_in[4];
        }
    case 3:
        if(fabs(vars_in[1]) <= DBL_EPSILON){
            return 0;
        }else{
            return - GRAVITATIONAL_CONSTANT * vars_in[5]*vars_in[1]/pow(pow(vars_in[1],2)+pow(vars_in[2],2), 1.5);   
        }
    case 4:
        if(fabs(vars_in[2]) <= DBL_EPSILON){
            return 0;
        }else{
            return - GRAVITATIONAL_CONSTANT * vars_in[5]*vars_in[2]/pow(pow(vars_in[1],2)+pow(vars_in[2],2), 1.5);
        }
    default:
        return 0;
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
                if(fabs(xdiff) <= DBL_EPSILON){
                    continue;
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
                if(fabs(ydiff) <= DBL_EPSILON){
                    continue;
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

double free_3d_orbit_functions(double * vars_in, int function_ref){
    /*
    Our variables and constants are as follows:
    0 --> time
    1 --> xpos of body 0
    2 --> ypos of body 0
    3 --> zpos of body 0
    4 --> xvel of body 0
    5 --> yvel of body 0
    6 --> zvel of body 0
    7 --> mass of body 0
    .
    .
    (body count * 7 + 1)
      --> time limit
    */

    if (function_ref > body_count * 7){
        return 0;
    }

    int i, actual_function_ref = function_ref % 7;
    int a = function_ref - actual_function_ref; // array offset for this body.
    int body_index = a / 7;

    #ifdef VERBOSE_DEBUG
    printf("fr = %d, afr = %d, a = %d, bi = %d\n", function_ref, actual_function_ref, a, body_index);
    #endif

    double xdiff, ydiff, zdiff, total_acceleration;
    switch(actual_function_ref){
    case 1:
        #ifdef VERBOSE_DEBUG
        printf("body %d x-velocity = %lf (index %d)\n", body_index, vars_in[a + 4], a+4);
        #endif
        return vars_in[a + 4];
    case 2:
        #ifdef VERBOSE_DEBUG
        printf("body %d y-velocity = %lf (index %d)\n", body_index, vars_in[a + 5], a+5);
        #endif
        return vars_in[a + 5];
    case 3:
        #ifdef VERBOSE_DEBUG
        printf("body %d z-velocity = %lf (index %d)\n", body_index, vars_in[a + 6], a+6);
        #endif
        return vars_in[a + 6];
    case 4:
        // x-acceleration
        total_acceleration = 0;
        for(i=0;i<body_count;i++){
            if(i != body_index){
                xdiff = vars_in[a + 1] - vars_in[i*7 + 1];
                ydiff = vars_in[a + 2] - vars_in[i*7 + 2];
                zdiff = vars_in[a + 3] - vars_in[i*7 + 3];

                // check for collision
                if(fabs(xdiff) <= DBL_EPSILON){
                    continue;
                }

                #ifdef VERBOSE_DEBUG
                printf("mass of body %d = %lf\n", i, vars_in[i * 7 + 7]);
                printf("body %d-->%d xdiff = %lf, ydiff = %lf, zdiff = %lf\n", body_index, i, xdiff, ydiff, zdiff);
                #endif
                total_acceleration -=
                    GRAVITATIONAL_CONSTANT * vars_in[i * 7 + 7] * xdiff / pow(pow(xdiff, 2) + pow(ydiff, 2) + pow(zdiff, 2), 1.5);
            }
        }
        #ifdef VERBOSE_DEBUG
        printf("body %d x-acceleration = %lf\n", body_index, total_acceleration);
        #endif
        return total_acceleration;
    case 5:
        // y-acceleration
        total_acceleration = 0;
        for(i=0;i<body_count;i++){
            if(i != body_index){
                xdiff = vars_in[a + 1] - vars_in[i*7 + 1];
                ydiff = vars_in[a + 2] - vars_in[i*7 + 2];
                zdiff = vars_in[a + 3] - vars_in[i*7 + 3];

                // check for collision
                if(fabs(ydiff) <= DBL_EPSILON){
                    continue;
                }

                #ifdef VERBOSE_DEBUG
                printf("body %d-->%d xdiff = %lf, ydiff = %lf, zdiff = %lf\n", body_index, i, xdiff, ydiff, zdiff);
                #endif
                total_acceleration -=
                    GRAVITATIONAL_CONSTANT * vars_in[i * 7 + 7] * ydiff / pow(pow(xdiff, 2) + pow(ydiff, 2) + pow(zdiff, 2), 1.5);
            }
        }
        #ifdef VERBOSE_DEBUG
        printf("body %d y-acceleration = %lf\n", body_index, total_acceleration);
        #endif
        return total_acceleration;
    case 6:
        // z-acceleration
        total_acceleration = 0;
        for(i=0;i<body_count;i++){
            if(i != body_index){
                xdiff = vars_in[a + 1] - vars_in[i*7 + 1];
                ydiff = vars_in[a + 2] - vars_in[i*7 + 2];
                zdiff = vars_in[a + 3] - vars_in[i*7 + 3];

                // check for collision
                if(fabs(zdiff) <= DBL_EPSILON){
                    continue;
                }

                #ifdef VERBOSE_DEBUG
                printf("body %d-->%d xdiff = %lf, ydiff = %lf, zdiff = %lf\n", body_index, i, xdiff, ydiff, zdiff);
                #endif
                total_acceleration -=
                    GRAVITATIONAL_CONSTANT * vars_in[i * 7 + 7] * zdiff / pow(pow(xdiff, 2) + pow(ydiff, 2) + pow(zdiff, 2), 1.5);
            }
        }
        #ifdef VERBOSE_DEBUG
        printf("body %d z-acceleration = %lf\n", body_index, total_acceleration);
        #endif
        return total_acceleration;
    default:
        return 0;
    }
}

int free_3d_orbit_runge_kutta_4th(double * vars_in, double * vars_out, double step){
    // this is a 3D simulation, so the number of variables is 7 * body_count + 1
    // (6 variables: xpos, ypos, zpos, xvel, yvel, zvel and 1 constant/body: mass of body)
    runge_kutta_4th(&free_3d_orbit_functions, vars_in, vars_out, 7 * body_count + 1, 1, step);
    // check the terminating condition. In this case, a time limit.
    if(vars_out[0] >= vars_in[7 * body_count + 1]){
        return STOP_ITERATING;
    }else{
        return CONTINUE_ITERATING;
    }
}
