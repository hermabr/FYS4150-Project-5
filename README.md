# FYS4150 Project 5 - Schrödinger equation

# Introduction
This repo contains code for simulating the two-dimensional time-dependent Schrödinger equation.
To do so, the Crank-Nicolson method for PDE's has been implemented. We apply this to a 
double-slit-in-a-box setup and variations on this. To better illustrate the results, there is 
a Python script for generating animations and various plots.

The C++ script uses config1.in and config2.in with at double slit setup, and config3.in for
a simple, double and triple slit setup. Running will produce in all five binary files 
containing information about the system at each timestep.

# Usage

## Python

```

```

## C++

Compiling

```bash
make
```

Running

```

```

# Structure

```
.
├── include - header files
├── makefile - the makefile used to compile the project
├── output - all the output generated by the programs
   ├── animations - all the animations generated by the Python script
   ├── plots - all the plots generated by the Python script
   ├── data - all the binary files generated by the C++ program. Filenames indicate configuration and number of slits
├── README.md - README-file
├── config1.in - The first configuration of the system to run. Using double slits.
├── config2.in - The second configuration of the system to run. Using double slits.
├── config2.in - The third configuration of the system to run. Using simple, double and triple slits
└── src - the code
   ├── main.cpp - the main script for running c++ code

```

# Code documentation

For more in depth description of files see the docstring in the python file for python and the header files for c++.

# Notation and indices

The indexing in the C++ code differs from the notation in the report. Here are the major differences:

* We use i to index the y-coordinates and j to index the x-coordinates
* All indexing starts at 0.