.. _e2i_main:

Inequivalent Impurities Solver Routines 
============================================

The module :f:mod:`E2I_MAIN` contains the procedures that initializes,
launch and finalize multiple **EDIpack** solvers for independent inequivalent quantum impurity
problems. 

The functionalities of the main `EDIpack` functions are here naturally
extended to the case of multiple impurities:  `ed_init` setups and allocates all the
internal variables and memory used in the code,  which remain
available to the user until  `ed_finalize` is called.  

The `ed_solve` functions aims to diagonalize the
impurity problems, evaluate the dynamical response functions and local
observables using different parallelization schemes according to the
value of the input variable :f:var:`mpi_lanc`: serial on the inequivalent
impurities and parallel in the diagonalization, parallel on the
inequivalent impurities and serial in the diagonalization.   

.. f:automodule::   e2i_main

