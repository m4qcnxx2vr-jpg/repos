Two-sided Jacobi SVD
====================

The task is to implement the two-sided Jacobi algorithm for the singular
value decomposition (SVD) of a real square matrix A:

    A = U D V^T

where D is diagonal with non-negative entries, and U and V are orthogonal
matrices.

This project does three things:

1. Implements the two-sided Jacobi SVD algorithm and tests it against the
   above properties.
2. Times the implementation for matrices of increasing size and fits the
   measured runtime to a cubic model, to check that it behaves as the
   theoretical O(N^3) cost per sweep predicts.
3. Compares the SVD timing against the plain (one-sided) Jacobi eigenvalue
   algorithm, to see what the extra work needed to handle general
   (non-symmetric) matrices actually costs in practice.

Files
-----

main.cc       : tests the SVD implementation and runs timing measurements
                for both SVD and the plain Jacobi eigenvalue algorithm
matrix.h      : vector and matrix declarations
matrix.cc     : vector and matrix implementations
svd.h         : declaration of the SVD class
svd.cc        : implementation of the two-sided Jacobi SVD algorithm
jacobi.h      : declaration of the Jacobi eigenvalue class
jacobi.cc     : implementation of the cyclic Jacobi eigenvalue algorithm
Makefile      : builds and runs the project
timing.gpi    : gnuplot script for the SVD timing plot
svd_vs_jacobi.gpi : gnuplot script comparing SVD and Jacobi timing

Build and run
-------------

Run

    make

This compiles the program, runs it, and creates

    Out.txt
    out.times.data
    out.times.svg
    svd_vs_jacobi_timing.svg

Out.txt contains the numerical tests. The file out.times.data contains
timing measurements, and out.times.svg contains a plot of the runtime
together with a fitted cubic model. Additionally we times the Jacobi eigenvalue algorithm at the same matrix sizes
and produces svd_vs_jacobi_timing.svg, comparing it against SVD.

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

SVD vs. Jacobi comparison
--------------------------

The two-sided SVD algorithm is built directly on top of the plain Jacobi
eigenvalue algorithm: for each pair (p,q), SVD performs one extra rotation
(to symmetrize the p,q block) before applying the same rotation that plain
Jacobi uses to zero out the off-diagonal entry. So SVD does roughly double
the rotation work per pair, every sweep, compared to Jacobi.

To check this, the same test matrices were timed with both algorithms and
fitted to t(N) = a N^3. Both fits land close to cubic, as expected from the
O(N^3) per-sweep cost, and the SVD curve sits consistently above the Jacobi
curve with a larger fitted constant a, which matches the extra rotation SVD
performs at every pair. In short: the comparison shows numerically what the
algorithm tells us in theory, that supporting general (non-symmetric)
matrices via SVD costs roughly a constant factor more than the symmetric-only
Jacobi algorithm it is built from, not a higher power of N.

Use of previous code and AI assistance
--------------------------------------

I gave a rough structure for the project, including the choice of files, the
testing setup, the Makefile structure, and the timing/gnuplot setup. AI wrote
the files on top of that structure.
