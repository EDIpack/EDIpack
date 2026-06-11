.. _main:

Impurity Problem Solver Routines 
====================================

.. preferred-crossrefs::
   :ed_solve: f/ed_main/ed_solve

The module :f:mod:`ED_MAIN` contains the procedures that initializes,
launch and finalize the **EDIpack** solver for the quantum impurity
problem.


The initialization, :f:func:`ed_init` setups and allocates all the
internal variables and memory used in the code,  which remain
available to the user until :f:func:`ed_finalize` is called.  

The key procedure is :f:func:`ed_solve` which aims to diagonalize the
impurity problem, evaluate the dynamical response functions and local
observables making them available to user through input/output
procedures. 

.. f:automodule::   ed_main


   
