Stored Hamiltonian :math:`H\times\vec{v}`  
==============================================

The module :f:mod:`ED_HAMILTONIAN_SUPERC_STORED_HXV` constructs the
the different contributions to the sector Hamiltonian, storing them into different
:f:var:`sparse_matrix` instances: :f:var:`sph0` for the electronic
part and three others for the phononic and electron-phonon terms. 

The main output of this module are the matrix vector products performed
using the stored sparse matrices :math:`\vec{w} = H\times \vec{v}`.



.. f:automodule::  ed_hamiltonian_superc_stored_hxv

