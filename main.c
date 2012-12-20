/*
    (c) Tom Robbins 2012

*/

#include "main.h"

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
    printf("    --resume\n");
    printf("        Resumes an existing simulation. The file specified \n        is appended to rather than overwritten (unless --stdout\n        is specified). The only numerical argument required is then the time step.\n");
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

int main(int argc, char ** args){
    char *flag_array[] = FLAG_ARRAY;
    int flags = process_flags(argc, args, FLAG_ARRAY_SIZE, flag_array);
    double * numeric_args = malloc(sizeof(double) * (argc - 1));
    int numeric_arg_count = process_numeric_args(argc, args, numeric_args);

    if(flags & FLAG_HELP){
        help();
        return 0;
    }

    if(flags & FLAG_RESUME){
        printf("Unimplemented. Sorry.\n");
        return 1;
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
                printf("Simple 3D simulation has not been implemented. Please specify --2D.\n");
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
                // make sure a valid number of args have been entered
                if(numeric_arg_count >= 10 && (numeric_arg_count - 3) % 7 == 0){
                    body_count = (numeric_arg_count - 3) / 7;
                    char ** labels = malloc(sizeof(char *) * (body_count * 7 + 2));
                    labels[0] = "time";
                    int i;
                    for(i=0;i<body_count;i++){
                        labels[i * 7 + 1] = malloc(8);
                        sprintf(labels[i * 7 + 1], "%d.xpos", i);
                        labels[i * 7 + 2] = malloc(8);
                        sprintf(labels[i * 7 + 2], "%d.ypos", i);
                        labels[i * 7 + 3] = malloc(8);
                        sprintf(labels[i * 7 + 3], "%d.zpos", i);
                        labels[i * 7 + 4] = malloc(8);
                        sprintf(labels[i * 7 + 4], "%d.xvel", i);
                        labels[i * 7 + 5] = malloc(8);
                        sprintf(labels[i * 7 + 5], "%d.yvel", i);
                        labels[i * 7 + 6] = malloc(8);
                        sprintf(labels[i * 7 + 6], "%d.zvel", i);
                        labels[i * 7 + 7] = malloc(8);
                        sprintf(labels[i * 7 + 7], "%d.mass", i);
                    }
                    labels[body_count * 7 + 1] = "time_limit";
                    set_up_runge_kutta_4th(7 * body_count + 1, 1);
                    iterate_to_file(&free_3d_orbit_runge_kutta_4th, 7 * body_count + 2, numeric_args, numeric_args[numeric_arg_count - 1], labels, fout);
                    free_runge_kutta_4th();
                }else{
                    printf("Invalid number of numerical arguments. Need at least 10 with 7 arguments for each body. %d given\n", numeric_arg_count);
                    help();
                    return 1;
                }
            }else{
                printf("Please specify either --2D or --3D\n");
                help();
                return 1;
            }
        }else{
            printf("Please specify either --simple or --free\n");
            help();
            return 1;
        }
    }

    return 0;
}
