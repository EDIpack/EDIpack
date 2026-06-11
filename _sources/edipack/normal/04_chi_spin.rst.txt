.. _chi_spin:
   
Spin Susceptibility
==================================

In :f:mod:`ed_chi_spin` we evaluate the impurity spin-spin
susceptibility, defined as:

.. math::

   \chi^{z}_{ab}(\omega) = \frac{1}{\cal
   Z}\sum_n e^{-\beta E_n} (\langle n | S^z_a [\omega-H]^{-1} S^z_b  | n \rangle +
   \langle n | S^z_b [\omega+H]^{-1} S^z_a  | n \rangle )

where :math:`S^z_a` is the z-component of the spin operator of the
orbital :math:`a` and :math:`\omega \in {\mathbb C}`. As for the
Green's functions, the susceptibility is evaluated using the dynamical
Lanczos method: a) the partial tridiagonalization of the 
sector Hamiltonian :math:`H` with quantum numbers
:math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` on the Krylov
basis of :math:`S^z_a|n\rangle` is obtained; b) the resulting
tridiagonal matrix is further diagonalized to obtained excitations
amplitudes or **weights**  :math:`\langle m | S^z_a | n \rangle` for
any state :math:`| m \rangle` in the spectrum (*without knowing the
state itself* ) and the excitations energies :math:`\delta E = E_m -
E_n` or **poles**; c) an controlled approximation to the
Kallen-Lehmann sum is constructed for  :math:`a,b=1,\dots,N_{\rm
orb}`. 


.. note::

   A more general susceptibility function for the other components of
   the spin operators :math:`S_x, S_y` should be implemented. 
   
.. f:automodule::  ed_chi_spin
