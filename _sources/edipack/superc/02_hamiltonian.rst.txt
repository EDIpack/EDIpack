Hamiltonian Setup
============================


The :f:mod:`ED_HAMILTONIAN_SUPERC` module contains the implementation of
the  superconducting sector Hamiltonian constructor.
The Hamiltonian is then used elsewhere to obtain
the low part of the energy spectrum or construct the dynamical response
functions. 

The symmetry sector electronic Hamiltonian has the form:

.. math::

   H_e = H_{\rm imp} + H_{\rm int}+ H_{\rm bath} + H_{\rm hyb}

where :math:`H_{\rm imp}` describes all the impurity Hamiltonian
terms

.. math::

      H_{\rm imp}  = \sum_{ab,\sigma} \left[ H^{\rm imp}_{ab,\sigma} -\mu\delta_{ab}\right]
      c^{\dagger}_{a\sigma}c_{b\sigma} + \sum_{a} P_{a} c_{a\uparrow} c_{a\downarrow} + H.c.

where :math:`P_a` is the external pair field coupled to the orbital
resolved pair amplitudes.  

The term :math:`H_{\rm int}` describes the local electron-electron
interaction in the generic Hubbard-Kanamori form with tunable
parameters:

.. math::

      H_{\rm int}  = \sum_{a} U_a n_{a\uparrow}n_{a\downarrow} +
      U'\sum_{a<b,\sigma} n_{a\sigma}n_{b\bar{\sigma}} +
      (U'-J)\sum_{a<b,\sigma} n_{a\sigma}n_{b\sigma} +
      J_x \sum_{ab} c^{\dagger}_{a\uparrow}c^{\dagger}_{b\uparrow}c_{a\downarrow}c_{b\uparrow} +
      J_p \sum_{ab}c^{\dagger}_{a\uparrow}c^{\dagger}_{a\downarrow}c_{b\downarrow}c_{b\uparrow}

where :math:`U_a=` :f:var:`uloc`,   :math:`U'=` :f:var:`ust`,
:math:`J=` :f:var:`jh`, :math:`J_x=` :f:var:`jx` and :math:`J_p=`
:f:var:`jp`. 

The :math:`H_{\rm bath}` and  :math:`H_{\rm hyb}` describe, respectively, the bath terms of the
Hamiltonian and the hopping between the impurity and the bath levels.

.. math::

      H_{\rm bath} & = \sum_p \sum_{ab\sigma} h^p_{ab\sigma}p^{\dagger}_{a\sigma}p_{b\sigma} +
      \sum_p \sum_{a} \Delta^p_a p_{a\downarrow}p_{a\uparrow} + H.c.\\\\
      H_{\rm hyb} & = \sum_p \sum_{a\sigma} V^p_{a\sigma}
      c^{\dagger}_{a\sigma} p_{p\sigma}  + H.c.


In addition we include electron-phonon coupling terms of the form:

.. math::

   H = H_{ph} \otimes 1_e + H_{ph-ph}\otimes H_{e-ph}
   

.. f:automodule::  ed_hamiltonian_superc
