.. _observables:

:ref:`Impurity Observables <imp_observables>`
---------------------------------------------------------------------------

The :f:mod:`ED_OBSERVABLES`  provides a single interface  to all the different
observables and static correlation calculation procedures available in the code.
This is used in the :f:mod:`ED_MAIN` Fortran API. 

.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:

   11_observables/impurity_observables
   

:ref:`Normal mode <observables_normal>`
---------------------------------------------------------------------------

This set of modules implements the  evaluation of impurity observables
and other static correlations  assuming :math:`\vec{Q}=\left[\vec{N}_\uparrow,\vec{N}_\downarrow \right]`.
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


   normal/05_observables


:ref:`Superconductive mode <observables_superc>`
---------------------------------------------------------------------------

This set of modules implements the  evaluation of impurity observables
and other static correlations  assuming 
:math:`\vec{Q}\equiv S_z=N_\uparrow-N_\downarrow`.

This case corresponds to the superconductive phase with :math:`s-`
wave pairing.


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:


   superc/05_observables



:ref:`Non-SU(2) mode <observables_nonsu2>`
---------------------------------------------------------------------------

This set of modules implements the  evaluation of impurity observables
and other static correlations   assuming 
:math:`\vec{Q}\equiv N_{tot}=N_\uparrow+N_\downarrow`.

This case corresponds to the normal phase in the absence of spin
conservation, as for instance in presence of Spin-Orbit coupling.  


.. toctree::
   :maxdepth: 2
   :glob:
   :hidden:


   nonsu2/05_observables

