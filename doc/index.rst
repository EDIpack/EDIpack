Welcome
################


.. sidebar:: EDIpack |PROJECT_VERSION|

      .. image:: _static/pictures/logo_github.png
         :width: 75%
         :align: center
         :target: https://github.com/EDIpack/EDIpack

**EDIpack** is an Exact Diagonalization solver for the solution of generic 
Quantum Impurity problems, exploiting MPI distributed memory parallelization.

List of Features
"""""""""""""""""""
* Support for zero and low temperature calculations
* Support for multiple star bath geometries
* Support for superconductive systems
* Support for systems where spin degrees freedom is not fully conserved
* One and two-particle response function on imaginary frequency, real frequency and imaginary time axes
* Support for reduce density matrices
* Support for Holstein phonon modes
* Support for inequivalent impurity sites through a Real-Space DMFT module
* APIs for FORTRAN, C++, Python and Julia
* Interoperability with the TRIQS software suite through a fully-featured thin compatibility layer
* Published in `SciPost Physics Codebases <https://scipost.org/SciPostPhysCodeb.58>`_

Authors
""""""""""

**EDIpack** is collective effort by different authors:

* `Adriano Amaricci`_ (leading author)
  
* `Lorenzo Crippa`_
  
* `Samuele Giuli`_

* `Gabriele Bellomia`_

* `Igor Krivenko`_
  
* Massimo Capone

.. _Adriano Amaricci: https://github.com/aamaricci
.. _Lorenzo Crippa: https://github.com/lcrippa    
.. _Samuele Giuli: https://github.com/SamueleGiuli
.. _Gabriele Bellomia: https://github.com/beddalumia
.. _Igor Krivenko: https://github.com/krivenko



Installation
"""""""""""""

:doc:`dependencies`
     Software requirements to install the |edipack| library.
     
:doc:`installation`
     Build, install and configure the library in the OS.
     

Usage
""""""""

:doc:`quickstart`
     A quick start guide with two simple examples.

:doc:`examples`
     Further examples showcasing some potentialities of the software. 



Fortran library
"""""""""""""""""

:doc:`edipack`
     A detailed presentation  of the library with a thorough 
     description of the relevant modules, data types and procedures.

:doc:`edipack2ineq`
     The inequivalent impurities extension of |edipack|


Compatibility
"""""""""""""""""

:doc:`edipack_cbindings`
     The Fortran-C interface for |edipack| and |edipack2ineq|
     
:doc:`edipack_python`
     Link to the documentation of the EDIpack python API
     EDIack2py, as well as the TRIQS compatibilty layer
     

Source navigation
""""""""""""""""""""

:ref:`Visual representation <browsecode>` of the interdependent modules
and files, with links to the source.



.. Hidden TOCs

.. toctree::
   :caption: Installation
   :maxdepth: 2
   :hidden:

   dependencies
   installation

.. toctree::
   :caption: Usage
   :maxdepth: 2
   :hidden:

   quickstart
   examples
   

.. toctree::
   :caption: Core Library
   :maxdepth: 2
   :hidden:
      
   edipack


.. toctree::
   :caption: EDIpack2ineq
   :maxdepth: 2
   :hidden:

   edipack2ineq

   
.. toctree::
   :caption: Compatibility
   :maxdepth: 2
   :hidden:

   edipack_cbindings
   
   edipack_python
   

.. toctree::
   :caption: Source navigation
   :maxdepth: 2
   :hidden:

   browsecode

.. toctree::
   :caption: External Links

   EDIpack on GitHub <https://github.com/edipack/EDIpack>   
   EDIpack2py on GitHub <https://github.com/edipack/EDIpack2py>
   EDIpack2Triqs on GitHub <https://github.com/edipack/edipack2triqs>
   SciFortran on GitHub <https://github.com/SciFortran/SciFortran>





   

