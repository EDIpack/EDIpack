.. _rdm:

:ref:`Impurity Reduced Density Matrix <rdm>`
---------------------------------------------------------------------------

The :f:mod:`ED_RDM`  provides a single interface to the evaluation of
the Reduced Density Matrix (RDM) :math:`\rho_{imp}={\rm
Tr}_{bath}\rho` for any value of :f:var:`ed_mode`.
This is used in the :f:mod:`ED_MAIN` Fortran API. 

.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:
      
   12_rdm/impurity_rdm
   

:ref:`Normal mode <rdm_normal>`
---------------------------------------------------------------------------

This module implements the  evaluation of impurity RDM
assuming the conserved quantum numbers are
:math:`\vec{Q}=\left[\vec{N}_\uparrow,\vec{N}_\downarrow \right]`.
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
   :hidden:
      

   normal/06_rdm


:ref:`Superconductive mode <rdm_superc>`
---------------------------------------------------------------------------

This module implements the  evaluation of impurity RDM 
assuming  :math:`\vec{Q}\equiv S_z=N_\uparrow-N_\downarrow`.

This case corresponds to the superconductive phase with :math:`s-`
wave pairing.


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:
      

   superc/06_rdm



:ref:`Non-SU(2) mode <rdm_nonsu2>`
---------------------------------------------------------------------------

This module implements the  evaluation of impurity RDM 
assuming 
:math:`\vec{Q}\equiv N_{tot}=N_\uparrow+N_\downarrow`.

This case corresponds to the normal phase in the absence of spin
conservation, as for instance in presence of Spin-Orbit coupling.  


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:
      

   nonsu2/06_rdm

