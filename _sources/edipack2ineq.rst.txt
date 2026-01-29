.. _edipack2ineq:

EDIpack for inequivalent impurities
#################################################################################

|edipack2ineq| is a complementary sub-library of |edipack|, which
extends its functionalities to systems with several 
independent impurities. These maye emerge from RealSpace-DMFT treatment of
unit cells with different inequivalent atoms or large super-cells with
broken translational symmetries (heterostructures, disordered
systems, etc.).

As for its scope, the structure of this sub-library closely mimics that
of |edipack| adding, where needed, new procedures adapted to the 
inequivalent impurities case. A standard fortran interface is then used to group
all the procedures (from `EDIpack` and `EDIpack2ineq`) under the same
name.      






:ref:`Library Frontend <edipack2ineqmodule>`
=================================================================

The :f:mod:`EDIPACK2INEQ` module represents the main user interface
(or Fortran API). This modules gives access to all the available
procedures and variables as needed to solve quantum impurity problems.


.. toctree::
   :caption: Library Frontend
   :maxdepth: 2
   :hidden:

   edipack2ineq/01_edipack2ineqmodule
      

   
:ref:`Core Solver Routines <e2i_main>`
=================================================================

The module :f:mod:`E2I_MAIN` extends wrapping of the main algorithms into
three functions, one for each step: initialization,
solution and finalization. 

.. toctree::
   :caption: Core Solver Routines
   :maxdepth: 2
   :hidden:

   edipack2ineq/02_main


General Environment
=================================================================

This part of the library includes a set of global variables and
procedures  for `EDIpack2ineq`.  This includes:

:f:mod:`E2I_VARS_GLOBAL`: contains global shared
variables storing site dependent quantities in the memory

Finally, :f:mod:`E2I_AUX_FUNX`  defines procedures which are used throughout the code.

.. toctree::
   :maxdepth: 2
   :glob:

   edipack2ineq/03_general/*


   

:ref:`Quantum Impurity Bath  <e2i_bath>`
=================================================================

The module :f:mod:`E2I_BATH` contains some extensions to procedures
enabling user side control of the bath setup and the implementation of
conventional symmetry operations. 

      
.. toctree::
   :caption: Bath
   :maxdepth: 2
   :hidden:

   edipack2ineq/04_bath
   



   

:ref:`Input/Output <e2i_io>`
=================================================================

The module :f:mod:`E2I_IO` provides to extended input and output
procedures returning observables or dynamical functions per
inequivalent site. 

.. toctree::
   :caption: Input/Output
   :maxdepth: 2
   :hidden:

   edipack2ineq/05_io






   

:ref:`Bath Optimization <e2i_fit>`
=================================================================

.. preferred-crossrefs::
   :ed_chi2_fitgf: f/e2i_bath_fit/ed_chi2_fitgf

The module :f:mod:`E2I_BATH_FIT` contains the extensions of the
generic function :f:func:`ed_chi2_fitgf` to perform 
multiple independent minimizations for all inequivalent sites.  


.. toctree::
   :caption: Bath Optimization
   :maxdepth: 2
   :hidden:

   edipack2ineq/06_fit

 






