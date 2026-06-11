Hamiltonian Setup
============================


In :f:mod:`ED_HAMILTONIAN_NORMAL` we  setup and build the sector
Hamiltonian, which are then used elsewhere to obtain the low part of
the energy spectrum or to construct the dynamical response functions. 

The symmetry sector electronic Hamiltonian is:

.. math::

   H_e = \vec{H}_\downarrow \otimes \vec{1}_\uparrow + \vec{1}_\downarrow \otimes
   \vec{H}_\uparrow + H_d + H_{nd} 

considering the electron-phonon coupling one has:

.. math::

   H = 1_{ph} \otimes H_e + H_{ph} \otimes 1_e + H_{ph-ph}\otimes H_{e-ph}
   

where :math:`H_\sigma` have different shape according to the value of
:f:var:`ed_total_ud`. See  `j.cpc.2021.108261`_ for further
information about this.  

In this operational mode the vectors :math:`\vec{v}` are stored as
matrices :math:`v_{ i_\uparrow i_\downarrow}`. In parallel mode
anumber of :math:`Q_\downarrow=\frac{D_\downarrow}{N_{threads}}` of columns are assigned to each thread. 

The matrices :math:`H_\sigma`, :math:`H_d`, having a small memory footprint, are entirely stored on each thread. The :math:`H_{nd}` is
instead splits across the threads and applied using the :code:`allgather_mpi`  algorithm (see `j.cpc.2021.108261`_ ). 

The matrix-vector operation proceeds as follow:
First the diagonal part :math:`H_d` is
applied to the vector. This step is local in the memory of each
thread. 
The :math:`H_\uparrow` is then applied  to  each column of the vector
:math:`\vec{v}`. This operation is also local on each thread. Next the matrix :math:`H_\downarrow` is applied to the
rows of the vector. In parallel mode this required parallel
transposition of the vector matrix, which allows to perform this step
locally on each thread. Finally, is non-vanishing, the last term
:math:`H_{nd}` is applied. In parallel mode, this step requires to
perform a :f:func:`allgather_vector_mpi` which represents the
bottleneck of the parallel execution. 



.. _j.cpc.2021.108261: https://doi.org/10.1016/j.cpc.2021.108261

.. f:automodule::  ed_hamiltonian_normal

