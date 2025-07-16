.. _greensfunctions:

:ref:`Impurity Green's functions <imp_greensfunctions>`
---------------------------------------------------------------------------

The :f:mod:`ED_GREENS_FUNCTIONS` wraps the different Green's functions
calculation methods available in the code in a single procedure.  
This is used  in the :f:mod:`ED_MAIN` Fortran API. 

.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:

   09_greensfunctions/impurity_greensfunctions
   

:ref:`Normal mode <gf_normal>`
---------------------------------------------------------------------------

This set of modules implements the calculations of impurity dynamical
response functions, e.g. the Green's functions and different
susceptibilities,  assuming :math:`\vec{Q}=\left[\vec{N}_\uparrow,\vec{N}_\downarrow \right]`.
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
      
   normal/03_greensfunctions
   

:ref:`Superconductive mode <gf_superc>`
---------------------------------------------------------------------------

This set of modules implements the calculations of impurity dynamical
response functions, e.g. the Green's functions,  assuming 
:math:`\vec{Q}\equiv S_z=N_\uparrow-N_\downarrow`.

This case corresponds to the superconductive phase with :math:`s-`
wave pairing.


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:
      
   superc/03_greensfunctions



:ref:`Non-SU(2) mode <gf_nonsu2>`
---------------------------------------------------------------------------

This set of modules implements the calculations of impurity dynamical
response functions, e.g. the Green's functions,  assuming 
:math:`\vec{Q}\equiv N_{tot}=N_\uparrow+N_\downarrow`.

This case corresponds to the normal phase in the absence of spin
conservation, as for instance in presence of Spin-Orbit coupling.  


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:
      
   nonsu2/03_greensfunctions

