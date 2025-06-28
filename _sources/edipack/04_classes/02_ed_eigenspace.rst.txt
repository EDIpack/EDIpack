.. _eigenspace:

EigenSpace 
=======================

The  module :f:mod:`ED_EIGENSPACE` contains a class to store the
energy spectrum of the quantum impurity problem in
:f:var:`sparse_espace`. Each instance of the class is a linked list of
:f:var:`sparse_estate`, each containing: eigenvectors, eigenvalues, symmetry sector index and
twin states, i.e. states in sectors with exchanged quantum numbers which do
have same energy. The module provides a number of procedures
implementing natural operations over the :f:var:`sparse_espace`, such
as add or removing :f:var:`sparse_estate`, returning eigenvector,
eigenvalue or minimum/maximum of the energy. 

.. f:automodule::   ed_eigenspace
   :members: state_list,sparse_estate, sparse_espace, es_insert_state, es_add_state, es_return_dvector, es_return_cvector, es_return_sector, es_return_energy, es_delete_espace

