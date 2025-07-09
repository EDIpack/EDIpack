.. _edipack:

The EDIpack Fortran Library
#################################################################################

Here we give an overview of the structure of the |edipack| library,
with a detailed description of the relevant modules and procedures.


:ref:`Library Frontend <edipackmodule>`
=================================================================

The :f:mod:`EDIPACK` module represents the main user interface (or Fortran API) 
of the |edipack| library. This modules gives access to all the
available procedures and variables as needed to solve quantum impurity problems.  

.. toctree::
   :caption: Library Frontend
   :maxdepth: 2
   :hidden:
      
   edipack/01_edipackmodule
      

:ref:`Core Solver Routines <main>`
=================================================================
     
The module :f:mod:`ED_MAIN` wraps the |edipack| algorithm into
three functions, one for each step: initialization, solution and finalization. 

.. toctree::
   :caption: Core Solver Routines
   :maxdepth: 2
   :hidden:

   edipack/02_main


General Environment
=================================================================

This part of the library includes a set of modules which contain
variables or procedures serving the rest of the code. 
The set includes:

:f:mod:`ED_INPUT_VARS` which defines global input variables and input read method.

:f:mod:`ED_VARS_GLOBAL` containing global shared variables and classes,such as the key class  :f:type:`effective_bath` which gathers the different bath components according to value of :f:var:`bath_type` and :f:var:`ed_mode`.

:f:mod:`ED_AUX_FUNX`  defining procedures which are used throughout the code.

:f:mod:`ED_SETUP`  which handles global memory (de)allocation and the evaluation of the dimensions of the symmetry sectors.

:f:mod:`ED_PARSE_UMATRIX`  which describes how a custom interaction Hamiltonian can be initialized.



.. toctree::
   :maxdepth: 2
   :glob:

   edipack/03_general/*


Classes
=================================================================   

This is a group of modules implementing different classes which are
used throughout the software to deal with special needs.
The group includes:

The  :f:mod:`ED_SPARSE_MATRIX` class, which contains the implementation of a
suitable Compact Row Stored (CRS) sparse matrix, used to store the sparse Hamiltonians in each symmetry
sector. 

The :f:mod:`ED_EIGENSPACE` class, which implements an ordered
single linked list to efficiently store the lower part of the energy spectrum.


The :f:mod:`ED_GFMATRIX` class which stores the
results of the dynamical Lanczos calculation of dynamical correlation
functions in terms of weights and poles of the Kallen-Lehmann representation.  

.. toctree::
   :maxdepth: 2
   :glob:

   edipack/04_classes/*


  


:ref:`Symmetry Sectors <sectors>`
=================================================================

The :f:mod:`ED_SECTOR` module implements the construction of the symmetry
sectors for the three sets of quantum numbers :math:`\vec{Q}` considered in |edipack|, that we
recall are:

* :math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` for which
  either the number of total or orbital spin up and down electrons is
  conserved: **NORMAL** 

* :math:`\vec{Q}=S_z`  with conserved total magnetization:  **SUPERConducting**  

* :math:`\vec{Q}=N_{\rm tot}`  where spin degrees freedom is not fully
  conserved:  **NON-SU(2)**
  


.. toctree::
   :caption: Symmetry Sectors
   :maxdepth: 2
   :hidden:

   edipack/05_sectors



:ref:`Quantum Impurity Bath  <bath>`
=================================================================

In |edipack| the bath is handled using a Reverse Communication
Strategy. All the procedures designed to define or handle the
discretized bath as well as those to evaluate suitable functions of
the bath are grouped in set of modules.   

    
.. toctree::
   :caption: Bath
   :maxdepth: 2
   :hidden:

   edipack/06_bath
   

:ref:`Hamiltonian <hamiltonian>`
=================================================================

This part of the |edipack| code implements the setup of the
sector Hamiltonian in each  operational modes,  corresponding to the
choice of one of the  symmetries implemented in the code selected by the variable :f:var:`ed_mode` =  :code:`normal, superc, nosu2`. See
:f:mod:`ED_SECTOR` for more info about the symmetries implemented in
the code.

Any of three different mode is implemented in a distinct class of
modules performing all the  operations required for the construction
of the Hamiltonian and the definition of the corresponding matrix
vector product, required by the Arpack/Lanczos. 

.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:
      
   edipack/07_hamiltonian


   

:ref:`Exact Diagonalization <diag>`
=================================================================

This part of the |edipack| code implements the exact
diagonalization of the general, single-site, multi-orbital quantum
impurity problem in each of the specific symmetry implemented in the
code. The operational modes are selected by the variable
:f:var:`ed_mode` =  :code:`normal, superc, nosu2` (see
:f:mod:`ED_SECTOR` for more info about the symmetries implemented
in the code).

As above, we implemented the three different channels in
distinct class of modules. 

.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:

   edipack/08_diag
   



:ref:`Green's Functions  <greensfunctions>`
=================================================================
     
This part of the |edipack| code implements the calculation of the
impurity interacting Green's functions, self-energy functions and
impurity susceptibilities. Calculations are performed in  each operational mode,  corresponding to the
choice of the specific symmetry implemented in the code, i.e. which
quantum numbers are to be conserved. The operational modes are selected
by the variable :f:var:`ed_mode` =  :code:`normal, superc, nosu2` (see
:f:mod:`ED_SECTOR` for more info about the symmetries implemented in
the code).


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:
      
   edipack/09_greensfunctions


:ref:`Susceptibilities  <chifunctions>`
=================================================================
     
This part of the |edipack| code implements the calculation of the
impurity susceptibilities in different physical channels: spin,
charge, pair and excitonic. The calculations are performed for
:f:var:`ed_mode` =  :code:`normal`.  


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:
      
   edipack/10_chifunctions


   

:ref:`Observables <observables>`
=================================================================

This part of the |edipack| code implements the calculation of the
impurity observables and static correlations, such as density,
internal energy or double occupation. Calculations are performed in
each operational mode,  corresponding to the choice of the specific
symmetry implemented in the code, i.e. which quantum numbers are to be
conserved. The operational modes are selected by the variable
:f:var:`ed_mode` =  :code:`normal, superc, nosu2`. The observables
are printed in plain text files and are accessible to the user
through the routines lited in :f:mod:`ED_IO`  (see
:f:mod:`ED_SECTOR` for more info about the symmetries implemented in
the code).


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:

   edipack/11_observables


   
:ref:`Reduced Density Matrix <rdm>`
=================================================================

This part of the |edipack| code implements the calculation of the
impurity Reduced Density Matrix (RDM). Calculations are performed in
each operational mode,  corresponding to the choice of the specific
symmetry implemented in the code, i.e. which quantum numbers are to be
conserved. The operational modes are selected by the variable
:f:var:`ed_mode` =  :code:`normal, superc, nosu2`.
The evaluation of the RDM exploit an efficient sparse algorithm to
avoid tracing over the exponentially many bath configuration.
The RDM is saved to file and made available to the user 
through the routines listed in :f:mod:`ED_IO`. 


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:

   edipack/12_rdm


   

:ref:`Input/Output <ed_io>`
=================================================================

This module provides access to the results of the exact
diagonalization. All quantities such as dynamical response functions,
self-energy components  or impurity observables  can be retrieved  
using specific functions. Additionally we provide procedure to perform
on-the-fly re-calculation of the impurity Green's functions and
self-energy on a given arbitrary set of points in the complex
frequency domain.    

.. toctree::
   :caption: Input/Output
   :maxdepth: 2
   :hidden:
      
   edipack/13_io




:ref:`Bath Optimization <fit>`
=================================================================

.. preferred-crossrefs::
   :ed_chi2_fitgf: f/ed_bath_fit/ed_chi2_fitgf

In  this module we provide to the user a generic function
:f:func:`ed_chi2_fitgf` performing the minimization of a user provided
Weiss field against the corresponding model of non-interacting
Anderson Green's function with the aim of updating the user bath parameters.


.. toctree::
   :caption: Bath Optimization
   :maxdepth: 2
   :hidden:
      
   edipack/14_fit




