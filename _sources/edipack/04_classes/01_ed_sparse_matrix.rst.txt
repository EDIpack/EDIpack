.. _sparse_matrix:

Sparse Matrices 
=======================

Implements sparse matrices data structures as dedicated CSR matrices
:f:type:`sparse_matrix_csr`.
Each instance of :f:type:`sparse_matrix_csr` corresponds to a rank-1
array of a tuple of dynamically reallocated arrays,
:f:type:`sparse_row_csr`: 
one array store the columns indices and one array store the values of
the non-zero elements of a given matrix. 

This structure guarantees :math:`O(1)` access.
If MPI parallelization is enabled the structure includes an
automatic rows split balanced among the different threads and a
further subdivision of the columns/values in local and non-local
blocks.   



.. note::
   A more modern, object oriented, class for sparse matrices is
   available in SciFortran_

.. _SciFortran: https://github.com/SciFortran/SciFortran/tree/master/src/SF_SPARSE


.. f:automodule::   ed_sparse_matrix
   :members: sparse_row_csr,sparse_matrix_csr, sp_init_matrix, sp_insert_element,sp_dump_matrix, sp_insert_element,sp_set_mpi_matrix




