Two-sided Jacobi SVD
====================

This project implements the two-sided Jacobi algorithm for singular value
decomposition of a real square matrix.

For a matrix A, the program computes

    A = U D V^T

where D is diagonal with non-negative entries, and U and V are orthogonal
matrices.

Files
-----

main.cc       : tests the SVD implementation and runs timing measurements
matrix.h      : vector and matrix declarations
matrix.cc     : vector and matrix implementations
svd.h         : declaration of the SVD class
svd.cc        : implementation of the two-sided Jacobi SVD algorithm
Makefile      : builds and runs the project
timing.gpi    : gnuplot script for the timing plot

Build and run
-------------

Run

    make

This compiles the program, runs it, and creates

    Out.txt
    out.times.data
    out.times.svg

Out.txt contains the numerical tests. The file out.times.data contains
timing measurements, and out.times.svg contains a plot of the runtime
together with a fitted cubic model.

Clean
-----

Run

    make clean

to remove generated files.

Tests
-----

The program tests the decomposition by checking

    U D V^T ≈ A
    U^T U ≈ I
    V^T V ≈ I

It also checks that D is diagonal by computing the norm of the off-diagonal
elements.

Timing
------

The timing test measures the runtime of the SVD algorithm for different
matrix sizes. The gnuplot script fits the data to a cubic model

    t(N) = a N^3

This is expected because one Jacobi rotation costs O(N), and one sweep goes
through O(N^2) pairs of indices. Therefore one sweep costs O(N^3).

Use of previous code and AI assistance
--------------------------------------

The files matrix.h and matrix.cc are reused from previous homework code from
the course and adapted for this project.

The rough project structure, including the choice of files, the testing setup,
the Makefile structure, and the timing/gnuplot setup, was made by me.

The files svd.h and svd.cc were written with AI assistance. I used AI to help
translate the two-sided Jacobi SVD algorithm into C++ code, while I tested,
compiled, modified, and checked the implementation myself.

The file main.cc and the Makefile were written by me with some AI assistance
for testing ideas and cleanup.