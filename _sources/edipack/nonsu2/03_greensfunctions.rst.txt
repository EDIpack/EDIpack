.. _gf_nonsu2:

Green's Function calculation: non-SU(2)
=====================================================================

In this module :f:mod:`ed_gf_nonsu2` the interacting impurity Green's
functions  :math:`\hat{G}(z)` are evaluated for :f:var:`ed_mode` =
:code:`nonsu2`.

Briefly, for any eigenstate :math:`|n\rangle` in the
:f:var:`state_list` contributing to the low energy part of the Hamiltonian spectrum the normal Green's functions:

.. math::

   G_{ab\sigma\sigma'}(z) = \frac{1}{\cal Z}\sum_n e^{\beta E_n}\langle n| c^\dagger_{a\sigma} [z-H]^{-1} c_{b\sigma'} |n
   \rangle + \langle n | c_{a\sigma} [z-H]^{-1} c^\dagger_{b\sigma'} | n \rangle

are evaluated using dynamical Lanczos method: a) the partial tridiagonalization of the
sector Hamiltonian :math:`H` with quantum numbers
:math:`\vec{Q}=N_{\rm tot} = N_\uparrow+N\downarrow` on the Krylov basis of :math:`c_{a\sigma}|n\rangle`
or  :math:`c^\dagger_{a\sigma}|n\rangle` is obtained; b) the resulting
tridiagonal matrix is further diagonalized to obtained excitations
amplitudes or **weights**  :math:`\langle m | c_{a\sigma} | n \rangle` or :math:`\langle m |
c^\dagger_{a\sigma} | n \rangle` for any state :math:`| m \rangle` in the
spectrum (*without knowing the state itself* ) and the excitations
energies :math:`\delta E = E_m - E_n` or **poles**; c) an controlled
approximation to the  Kallen-Lehmann sum is constructed for any
value of :math:`z\in{\mathbb C}` and :math:`a,b=1,\dots,N_{\rm orb}`,
:math:`\sigma=\uparrow,\downarrow`. 


While the Green's functions are evaluated in a given set of Matsubara
:f:var:`impgmats` and Real-axis points
:f:var:`impgreal`, the weights and the poles
obtained in this calculation are stored in a dedicated data
structure :f:var:`gfmatrix` for a fast recalculation on any given
intervals of frequencies in the complex plane.


Finally, the self-energy functions are constructed using impurity
Dyson equation :math:`\hat{\Sigma}(z) = \hat{G}^{-1}_0(z) -
\hat{G}^{-1}(z)`. 


.. f:automodule::  ed_gf_nonsu2

