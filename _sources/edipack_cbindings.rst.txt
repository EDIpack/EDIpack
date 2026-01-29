.. _edipack_cbindings:

EDIpack Fortran-C interface
=================================================================

EDIpack offers Fortran-C interoperability making use of the
:f:mod:`ISO_C_BINDING` Fortran module.

The :f:mod:`EDIPACK_C` module consists of a series of
C-type procedures which wrap one-by-one the relevant
functions as required to perform a full calculation.  

At compile time, if :code:`make` or :code:`make edipack_cbindings` is issued,
a dynamic library called :code:`libedipack_cbindings.so` will be created and
put in the same folder as the static Fortran library.
An API to another language will need to interface to this library. 

A Python API called EDIpack2py_ is provided and mantained.
An experimental Julia API called EDipack2jl_ is in development.

A :code:`C++` header file called :code:`edipack_cbindings.h` is installed into the
include directory at build time. The list of variables and functions therein provided
is documented in the following. 
Caution needs to be applied for functions requiring array input parameters: since arrays 
have to be passed as raw pointers, additional
arrays containing the dimensions of the former need to be passed.

An example :code:`C++` program using EDIpack is provided  in the examples_ folder 
of the GitHub repository.
A more exotic usecase is available here_, come fosse antani.

.. _EDIpack2py: https://github.com/edipack/EDIpack2py
.. _EDIpack2jl: https://github.com/edipack/EDIpack2jl
.. _examples: https://github.com/EDIpack/EDIpack/tree/master/examples/cpp
.. _here: https://github.com/lcrippa/prematurata_la_dmft


.. toctree::
   :maxdepth: 2
   :glob:

   c_bindings/*

