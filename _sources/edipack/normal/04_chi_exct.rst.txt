.. _chi_exct:

Exciton Susceptibility
============================

In :f:mod:`ed_chi_exct` we evaluate the impurity exciton-exciton  
susceptibility, defined as:

.. math::

   \chi^{X}_{ab}(\omega) = \frac{1}{\cal
   Z}\sum_m e^{-\beta E_m} (\langle m | X^\dagger_{ab} [\omega-H]^{-1} X_{ab}  | m \rangle
   + \langle m | X_{ab} [\omega+H]^{-1} X^\dagger_{ab}  | m \rangle )

where :math:`X_{ab}=S_{ab},T^x_{ab},T^z_{ab}` are, respectively, the singlet and
triplet :math:`x,z`  exciton operators: 

.. math::

   S_{ab}     & = \sum_{rs} c^\dagger_{ar} \sigma^0_{rs} c_{bs}\\
   T^x_{ab} & = \sum_{rs} c^\dagger_{ar} \sigma^x_{rs} c_{bs}\\
   T^z_{ab} & = \sum_{rs} c^\dagger_{ar} \sigma^z_{rs} c_{bs}


and :math:`\omega \in {\mathbb C}`. As for the
Green's functions, the susceptibility is evaluated using the dynamical
Lanczos method: a) the partial tridiagonalization of the 
sector Hamiltonian :math:`H` with quantum numbers
:math:`\vec{Q}=[\vec{N}_\uparrow,\vec{N}_\downarrow]` on the Krylov
basis of :math:`X_{ab}|m\rangle` is obtained; b) the resulting
tridiagonal matrix is further diagonalized to obtained excitations
amplitudes or **weights**  :math:`\langle p | X_{ab} | m \rangle` for
any state :math:`| p \rangle` in the spectrum (*without knowing the
state itself* ) and the excitations energies :math:`\delta E = E_p -
E_m` or **poles**; c) an controlled approximation to the
Kallen-Lehmann sum is constructed for  :math:`a,b=1,\dots,N_{\rm
orb}`. 






.. f:automodule::  ed_chi_exct

