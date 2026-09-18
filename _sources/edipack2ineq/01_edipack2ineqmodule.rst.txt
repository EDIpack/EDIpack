.. _edipack2ineqmodule:

EDIpack2ineq Fortran Module
======================================

:f:mod:`EDIPACK2INEQ` is the top module of the |EDIpack2ineq|
library. It provides access to all the functions extending |edipack|
to the case of inequivalent sites, realizing an addition to the Fortran API.

The user needs to invoke use of this module alongside the original
library one, as:

.. code-block:: fortran

   program test
     USE EDIPACK
     USE EDIPACK2INEQ
     ...

.. preferred-crossrefs::
   :e2i_bath/set_hgeneral: f/e2i_bath_replica/set_hgeneral
   :e2i_bath/set_hreplica: f/e2i_bath_replica/set_hreplica
   :e2i_bath/spin_symmetrize_bath:  f/e2i_bath_user/spin_symmetrize_bath
   :e2i_bath/orb_symmetrize_bath: f/e2i_bath_user/orb_symmetrize_bath
   :e2i_bath/orb_equality_bath: f/e2i_bath_user/orb_equality_bath
   :e2i_bath/ph_symmetrize_bath: f/e2i_bath_user/ph_symmetrize_bath
   :e2i_bath/ph_trans_bath: f/e2i_bath_user/ph_trans_bath
   :e2i_bath/break_symmetry_bath: f/e2i_bath_user/break_symmetry_bath
   :e2i_bath/enforce_normal_bath: f/e2i_bath_user/enforce_normal_bath
   :e2i_bath/save_array_as_bath: f/e2i_bath_user/save_array_as_bath

   		   
.. f:automodule::   edipack2ineq
