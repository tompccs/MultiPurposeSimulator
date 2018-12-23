MultiPurposeSimulator
=====================

A multi-purpose physics simulator. A collection of functions which can be mixed and matched to run a variety of physics simulations. Out of the box it works as an (up to) 3D n-body gravity simulator by 4th order Runge-Kutta method.

About
=====

The program is made up of 3 .c files (each with its own header .h file):

<dl>
  <dt>main.c</dt>
  <dd>Point of entry. Sets up simulations according to arguments and contains help text.</dd>

  <dt>lib.c</dt>
  <dd>A general purpouse library for solving differential equations. In theory, any method for solving DEs can be implemented by writing a function and passing a function pointer to iterate_to_file(). 4th order Runge-Kutta (runge_kutta_4th()) is implemented in this way. It also contains some helper functions for parsing command line arguments.</dd>

  <dt>rk_functions.c</dt>
  <dd>For each type of simulation, there are two functions. One is passed into iterate_to_file(), in turn calling runge_kutta_4th(), passing into it a function pointer to its corresponding set of functions, typically one for each variable in the system of differential equations. In theory, these functions can be used in solving their	system of equations by other methods, such as Gauss' or higher-order RK.</dd>
</dl>

There is also a python script for visualising the data. It reads from the standard input and takes two arguments. Argument 1 is the scale factor (expressed as half the number of distance units displayed) and argument 2 is a file to output the data to exactly as it is read from the standard input. You may also specify --noflush to read in a simulation (perhaps from an existing file) before drawing it.

To Compile
==========
```
gcc ./lib.h ./rk_functions.h ./main.h ./lib.c ./rk_functions.c ./main.c -lm -o ./simulator
```

Example Commands
================

```
./simulator ./out.csv --simple --orbit --2D 0 3E7 0 0 10E3 0 1E26 1E6 10
```
Simulates a simple 2D orbit with one body orbiting another static body and writes the result to ./out.csv

```
./simulator --stdout --orbit --free --2D starttime: 0s earth: 0 0 0 0 6E24 moon: 384.4E6 0 0 1E3 7.3E22 probe: 7E6 0 10.35E3 2.75E3 20E3kg timelim: 1E9s step: 10s | ./visual.py 1E9 ./out.csv
```
Simulates 3 bodies in 2D corresponding to the given initial conditions. Note that the annotations (such as 'starttime:', 'earth:', 'moon:', 'timelim:' and 'step:') around the numerical arguments are not required and will be ignored by the program (bodies are indexed numerically in the output) but they are useful for keeping track of which variable is which. The result is visualised and written to ./out.csv 

```
./simulator --stdout --3D --free --orbit timestart: 0 earth: pos 0, 0, 0 vel 0, 0, 0 mass 5.97E24 moon: pos 0, .5E8 -3.84E8 vel 1E3, 0, 0 7.3477E22 timelimit: 5E6 600 | ./visual.py 4E9 ./myout.csv
```
Simulates 2 bodies in 3D and visualises them with a 2D top-down projection, also writing result to ./myout.csv

More Info
=========
I'll be putting some more information on the wiki
