MPI all-2-all transposition
==============================================

This module :f:mod:`ED_HAMILTONIAN_NORMAL_COMMON` defines several
variables shared across the Hamiltonian setup in the :f:var:`ed_mode`
= :code:`normal` mode. It also contains the procedure :f:func:`vector_transpose_mpi` implementing the  MPI :code:`Allv-2-Allv`
parallel transposition of a matrix. This is the key function of the massively parallel execution of matrix-vector products discussed in `j.cpc.2021.108261`_.

.. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261

.. f:automodule::  ed_hamiltonian_normal_common
   :members: vector_transpose_mpi
