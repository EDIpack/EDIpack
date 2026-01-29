Hamiltonian Setup
============================


The :f:mod:`ED_HAMILTONIAN_NONSU2` module contains the implementation of
the  superconducting sector Hamiltonian constructor.
The Hamiltonian is then used elsewhere to obtain
the low part of the energy spectrum or construct the dynamical response
functions. 

The sector electron Hamiltonian has the form:

.. math::

   H_e = H_{\rm imp} + H_{\rm int}+ H_{\rm bath} + H_{\rm hyb}

where :math:`H_{\rm imp}` describes all the impurity Hamiltonian
terms

.. math::

      H_{\rm imp}  = \sum_{ab,\sigma} \left[ H^{\rm imp}_{ab,\sigma\sigma`} -\mu\delta_{ab}\delta_{\sigma\sigma`}\right]
      c^{\dagger}_{a\sigma}c_{b\sigma'} + \sum_{i=0,x,y,z}\sum_{a}
      c^\dagger_{a\sigma} T^i_{a,\sigma\sigma'}  c_{a\sigma'} + H.c.

where :math:`T^i_a` are the excitonic field components, each proportional to
the 4 generators of the :math:`SU(2)`, coupled to the exciton amplitudes.  

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

      H_{\rm bath} & = \sum_p \sum_{ab\sigma\sigma'}
      h^p_{ab\sigma\sigma'}p^{\dagger}_{a\sigma}p_{b\sigma'} + H.c.\\\\      
      H_{\rm hyb} & = \sum_p \sum_{a\sigma} V^p_{a\sigma}
      c^{\dagger}_{a\sigma} p_{p\sigma}  + H.c.


   
.. note::

   In the current implementation the :code:`nonsu2` mode does not
   support (yet?) the presence of electron-phonon coupling terms. 



.. f:automodule::  ed_hamiltonian_nonsu2

