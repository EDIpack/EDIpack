Library structure
#################################################################################

The  |edipack|  library is made of 3 distinct parts. 
Each part builds into a library (static or dynamic) with its own
internal dependencies.  Alongside the library, specific module files
are generated and stored in the same location during installation.  

The organization of the library structure is as follows: 

**EDIpack**: :ref:`edipack`
=======================================
The core library, implementing the parallel Lanczos-based exact diagonalization of a quantum
impurity problem and dynamical correlation functions
calculations.
It builds into the static libary `libedipack.a` and it is
accessible via the Fortran module :f:mod:`EDIPACK`


**EDIpack2ineq**: :ref:`edipack2ineq`
=======================================
This is a sub-library implementing the inequivalent impurities
extension for |edipack|, tackling the solution of systems represented by several
independent quantum impurity problems. Either serial or parallel execution
schemes with respect to the inequivalent impurities are provided, as
well as automatic deployment of the inequivalent files. 
This |edipack2ineq| builds into the same static library
`libedipack.a`, yet it is accessible only through the module
:f:mod:`EDIPACK2INEQ`.
The naming convention of the |edipack| procedures is preserved. 


**EDIpack C-bindings**: :ref:`edipack_cbindings`
=================================================
A module implementing the Fortran-C interface for |edipack|,
including |edipack2ineq| extension where included.
The  language interoperability leverages over the implicit Fortran
:code:`ISO_C_BINDING` features and is achieved using C-types to wrap the
relevant `EDIpack` procedures.
The C-bindings module serves as a key tool to setup additional API for 
|edipack|, as for instance the Python API (see EDIpack2py_ ).
This module builds into dynamic libraries
`libedipack_cbindings.so/.dylib`. 
Dependency: `edipack` and, if required, `edipack2ineq`.  



.. _EDIpack2py: https://github.com/edipack/EDIpack2py
