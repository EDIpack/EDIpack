.. _hamiltonian:

Hamiltonian Setup
-----------------------------------------------

.. f:automodule::   ed_hamiltonian
   :hide-output: True

The :f:mod:`ED_HAMILTONIAN` module wraps all the different
Hamiltonian setup procedures for each operational modes, as described below.


Normal mode
-----------------------------------------------

This set of modules implements the Hamiltonian setup for each symmetry
sector assuming :math:`\vec{Q}=\left[\vec{N}_\uparrow,\vec{N}_\downarrow \right]`.
Where :math:`\vec{N}_\sigma=N_\sigma` if the total number of electrons
with spin :math:`\sigma` is conserved (:f:var:`ed_total_ud` = T ) or
:math:`\vec{N}_\sigma=[ N_{1\sigma},\dots,N_{N_{orb}\sigma} ]` if the
number of electrons in the orbital :math:`\alpha=1,\dots,N_{orb}` and
spin :math:`\sigma` is conserved (:f:var:`ed_total_ud` = F). 

This case corresponds to the normal phase in presence of spin
conservation, possibly reduced to :math:`U(1)` in presence of long
range magnetic order along :math:`z` quantization axis of the spin
operator.   

.. toctree::
   :maxdepth: 2
   :glob:

   normal/02_hamiltonian
   normal/02_hamiltonian_common
   normal/02_hamiltonian_sparse
   normal/02_hamiltonian_direct


Superconductive mode
-----------------------------------------------

This set of modules implements the Hamiltonian setup for each symmetry
sector  assuming 
:math:`\vec{Q}\equiv S_z=N_\uparrow-N_\downarrow`.

This case corresponds to the superconductive phase with :math:`s-`
wave spin singlet pairing.


.. toctree::
   :maxdepth: 2
   :glob:

   superc/02_hamiltonian
   superc/02_hamiltonian_sparse
   superc/02_hamiltonian_direct



   
Non-SU(2) mode
-----------------------------------------------

This set of modules implements the Hamiltonian setup for each symmetry
sector  assuming 
:math:`\vec{Q}\equiv N_{tot}=N_\uparrow+N_\downarrow`.

This case corresponds to the normal phase in the absence of spin
conservation, as for instance in presence of Spin-Orbit coupling.  


.. toctree::
   :maxdepth: 2
   :glob:

   nonsu2/02_hamiltonian
   nonsu2/02_hamiltonian_sparse
   nonsu2/02_hamiltonian_direct
