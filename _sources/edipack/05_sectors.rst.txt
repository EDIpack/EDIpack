.. _sectors:

Quantum Numbers Sectors
============================

The direct diagonalization of the Hamiltonian for a finite system
becomes quickly impossible due to the exponential increase of its
size :math:`D^{Ns}`, where :math:`D` is the dimension of the local
Hilbert space (:math:`D=4` for spin-:math:`1/2` electrons).

A way to partially soften such severe behavior is to take into account
conserved quantities. For any quantity :math:`{\cal Q}` such that
:math:`[H,{\cal Q}]=0`,  we can separate the structure of the
Hamiltonian matrix into different blocks, each corresponding to a
discrete value :math:`q` of the spectrum of :math:`Q`.
Analyzing one-by-one such symmetry sector it becomes possible to
construct at least part of the energy spectrum.


The :mod:`ED_SECTOR` represents a key part of the ED code. Here the
symmetry sectors construction  essentially consist in the determination of the
injective map :math:`{\cal M}:{\cal S}(\vec{Q})\rightarrow{\cal F}`
relating the states of the sectors :math:`|i\rangle\in{\cal
S}(\vec{Q})` to the corresponding ones :math:`|I\rangle\in{\cal F}` in
the Fock space.

The map is constructed iterating over the states of the Fock spaces
while verified the quantum number conditions. If verified for a given
state :math:`|I\rangle` a new state is found and it is added to the
list. The following code snippet illustrates the case for a system
conserving the total magnetization :math:`S_z=N_\uparrow-N_\downarrow`

.. code-block:: fortran

       dim=0
       ! Iterate over down spin configurations
       do idw=0,2**Ns-1
          ! use bit operation to get :math:`N_\downarrow`
          ndw_= popcnt(idw)
	  ! Iterate over up spin configurations
          do iup=0,2**Ns-1
             ! use bit operation to get :math:`N_\uparrow`
             nup_ = popcnt(iup)
	     ! Get :math:`S_z=N_\uparrow-N_\downarrow`
             sz_  = nup_ - ndw_
	     ! if the evaluated magnetization is the required one:
	     ! increase counter and store the Fock state into the map
             if(sz_ == self%Sz)then
                dim=dim+1
                self%H(1)%map(dim) = iup + idw*2**Ns
             endif
          enddo
       enddo

   
The dimension of the symmetry sector, i.e. the size of  the map
:math:`{\cal M}(\vec{Q})`, is determined by a simple combinatorial
algorithm  (:math:`N` being the total number of electronic levels):


.. list-table:: Symmetry Sectors Dimensions
   :widths: 10 10 80
   :header-rows: 1

   * - :f:var:`ed_mode`
     - Quantum Numbers
     - Sector Dimension
       
   * - :code:`normal`
     - :math:`[\vec{N}_\uparrow,\vec{N}_\downarrow]`
     - :math:`\prod_{\alpha}\binom{N}{N_{\alpha\uparrow}}\binom{N}{N_{\alpha\downarrow}}`
       
   * - :code:`superc`
     - :math:`S_z=N_\uparrow-N_\downarrow`
     - :math:`\sum_i 2^{N-S_z-2i}\binom{N}{N-S_z-2i}\binom{S_z+2i}{i}`

   * - :code:`nonsu2`
     - :math:`N_{tot}=N_\uparrow+N_\downarrow`
     - :math:`\binom{2N}{N_{tot}}`




.. f:automodule::   ed_sector
   :members: build_sector,map_allocate,get_sector,get_quantumnumbers,get_dim



