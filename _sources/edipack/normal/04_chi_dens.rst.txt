.. _chi_dens:

Charge Susceptibility
============================


In :f:mod:`ed_chi_dens` we evaluate the impurity density-density 
susceptibility, defined as:

.. math::

   \chi^n_{ab}(\omega)  = \frac{1}{\cal
   Z}\sum_m e^{-\beta E_m} (\langle m | n_a [\omega-H]^{-1} n_b  | m \rangle +
   \langle m | n_b [\omega+H]^{-1} n_a  | m \rangle)

where :math:`n_a` is the density operator of the
orbital :math:`a` and :math:`\omega \in {\mathbb C}`. As for the
Green's functions, the susceptibility is evaluated using the dynamical
Lanczos method: a) the partial tridiagonalization of the 
sector Hamiltonian :math:`H` with quantum numbers
:math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` on the Krylov
basis of :math:`n_a|m\rangle` is obtained; b) the resulting
tridiagonal matrix is further diagonalized to obtained excitations
amplitudes or **weights**  :math:`\langle p | n_a | m \rangle` for
any state :math:`| p \rangle` in the spectrum (*without knowing the
state itself* ) and the excitations energies :math:`\delta E = E_p -
E_m` or **poles**; c) an controlled approximation to the
Kallen-Lehmann sum is constructed for  :math:`a,b=1,\dots,N_{\rm
orb}`. 



.. f:automodule::  ed_chi_dens

