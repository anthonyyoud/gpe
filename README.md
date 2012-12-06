Description
===========
Modular Gross-Pitaevskii equation solver in 3D parallelised using MPI.  The
time stepping scheme can be either explicit second-order Euler, explicit
fourth-order Runge-Kutta, or explicit fourth-order adaptive
Runge-Kutta-Fehlberg using the algorithm proposed in Numerical Recipes (section
16.2, page 708).  The spatial discretisation is via either second-order or
fourth-order accurate centred finite differences.  The boundary conditions can
be either periodic or reflective.  A variety of initial conditions are possible
including vortex lines, vortex rings, and rarefaction pulses.  Any combination
of these is possible by simply multiplying the desired initial conditions
together.

The code uses MPI (Message Passing Interface) as the method of parallelisation
in solving the GP equation.

__Note:__ Github thinks this is a Prolog project because of all the `.pro`
files.  These are IDL procedure files; IDL can be used for visualisation of the
results.  There is not a single line of Prolog anywhere here!

Directory structure
===================
* doc - Contains documentation.  Currently only the manual exists here.
* examples - Some example simulations already set up and ready to compile.
* idl - IDL routines for visualisation are kept in this directory.
* scripts - Various helper scripts are kept in this directory.
* src - The source for the code itself.

Manual
======
For instructions on how to use the code see the manual in the doc directory.

Before compiling for the first time
===================================
* In the `src` directory, copy `Makefile.dist` to `Makefile` and edit to suit.
* In each of the example directories, copy each `*.in.dist` file to `*.in`.
  If you want to make changes to the example parameters, edit the `*.in` files,
and leave the `*.in.dist` files unchanged.
