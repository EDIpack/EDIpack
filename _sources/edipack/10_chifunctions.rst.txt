.. _chifunctions:

:ref:`Impurity Susceptibilities <imp_chifunctions>`
---------------------------------------------------------------------------

The :f:mod:`ED_CHI_FUNCTIONS` wraps the  methods to evaluate the
dynamical impurity response functions. Because
these are often used as a proxy to test instability in some channel 
(e.g. spin, charge, etc.) the susceptibilities are accessible only for :f:var:`ed_mode` =
:code:`normal`. 

The main method in this module is called directly in the :f:mod:`ED_MAIN` Fortran API. 

.. toctree::
   :maxdepth: 1
   :hidden:

   10_chifunctions/impurity_chifunctions
   

:ref:`Spin channel <chi_spin>`
---------------------------------------------------------------------------

This module implements the calculations of impurity spin dynamical
response functions

.. toctree::
   :maxdepth: 1
   :hidden:
      
   normal/04_chi_spin
   

:ref:`Charge channel <chi_dens>`
---------------------------------------------------------------------------

This module implements the calculations of impurity charge dynamical
response functions

.. toctree::
   :maxdepth: 1
   :hidden:
      
   normal/04_chi_dens



:ref:`Pair channel <chi_pair>`
---------------------------------------------------------------------------

This module implements the calculations of impurity Cooper pair dynamical
response functions

.. toctree::
   :maxdepth: 1
   :hidden:
      
   normal/04_chi_pair




:ref:`Exciton channel <chi_exct>`
---------------------------------------------------------------------------

This module implements the calculations of impurity excitonic dynamical
response functions

.. toctree::
   :maxdepth: 2
   :hidden:
      
   normal/04_chi_exct
