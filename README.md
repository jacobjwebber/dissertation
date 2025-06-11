# README
8 years after i did it i have decided to open source my research on 3d room acoustics computer models using GPU.

This contains the following directories.

## 2d_example
This contains a very early prototype that was used for project preparation
to sho viability. 

## matlab_scripts
This contains the script written by Brian Hamilton which implements the 
new boundary conditions with the bounding box method.

## structured
This contains code from Craig Webb's PhD thesis that uses the old boundary
conditions and the structured method.

## semi-structured
This contains all the new work including

### lib
All the files reusable between experiments

### memory-capacity
A source file for conducting the memory capacity experiment and a makfile 
that compiles it with the necessary lib files

### stencil-performance
A source file for conducting the stencil performance experiment and a makfile 
that compiles it with the necessary lib files

### test
A test file and makefile

### results
A script that runs the memory-capacity experiment and writes the results to
a graph (using gnuplot) and a csv file.

To compile any of the experiments `cd` to either directory and use `make`.
`make clean` deletes object and binary files.
