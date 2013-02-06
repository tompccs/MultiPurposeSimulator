/*
	(c) Tom Robbins 2012

*/


#import "lib.h"
#include "rk_functions.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#define FLAG_ARRAY {"--orbit", "--simple", "--free", "--2D", "--3D", "--help", "--stdout", "--resume"}
#define FLAG_ARRAY_SIZE 7

#define FLAG_ORBIT 1
#define FLAG_SIMPLE 2
#define FLAG_FREE 4
#define FLAG_2D 8
#define FLAG_3D 16
#define FLAG_HELP 32
#define FLAG_STDOUT 64
#define FLAG_RESUME 128

int body_count;

int main(int argc, char ** args);
void help();
