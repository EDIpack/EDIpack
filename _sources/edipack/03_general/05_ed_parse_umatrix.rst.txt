.. _parse_umatrix:

Parse Umatrix
=======================

EDIpack offers the user two ways of providing the interaction Hamiltonian.

The default way is to use the in-built variables :f:var:`ULOC`, :f:var:`UST`,
:f:var:`JH`, :f:var:`JX`, :f:var:`JP`, which implement a Hubbard-Kanamori Hamiltonian for a maximum of 5 orbitals.
This is achieved by setting the :f:var:`ED_USE_KANAMORI` flag to :code:`T` in the input file.

Alternatively, the user can provide a plain text file containing the second-quantized
two-body operators. This is achieved by setting the :f:var:`ED_READ_UMATRIX` flag to :code:`T`
in the input file. Note that :f:var:`ED_READ_UMATRIX` and :f:var:`ED_USE_KANAMORI` cannot be :code:`T` at the same time.
The filename (without extension) has to be provide via the :f:var:`UMATRIX_FILE` 
variable (default name :code:`umatrix` ). The actual file in the execution folder will need to have
the :code:`.restart` suffix. When doing a real-space DMFT simulation, the prefix will need to be
:code:`_ineqXXXX.restart` , where :code:`XXXX` is a 4-digit incresing number corresponding to the index of 
the inequivalent site.

After the interaction terms have been set in either way, a properly formatted file 
called :f:var:`UMATRIX_FILE` :code:`.used` 
will be generated containing the list of all two-body operators.

Each umatrix file needs to have the following format:

.. code-block:: text

    NORB BANDS
    i1 j1 k1 l1 U_i1j1k1l1
    i2 j2 k2 l2 U_i2j2k2l2
     ...


:f:var:`NORB` is the number of orbitals. Empty lines and lines starting with :code:`#,%,!` are ignored.
The indices :code:`i,j,k,l` are of the form :code:`o s` where :code:`o` is an integer number for the 
orbital index (starting from 1, maximum :f:var:`NORB` ) and :code:`s` is :code:`u` for :math:`\uparrow`
spin and :code:`d` for :math:`\downarrow` spin. Each line of the file represents an operator of this form

.. math::
    \frac{1}{2}\sum_{i,j,k,l} c^{\dagger}_i c^{\dagger}_j U_{ijkl} c_l c_k
    
Note the inversion of the :code:`l` and :code:`k` operators with respect to :math:`U_{ijkl}`. 
The coefficient is given by the following formula

.. math::
    U_{ijkl} = \int \mathrm{d}x \int \mathrm{d}y \phi_i^{*}(x) \phi_j^{*}(y) U(x, y) \phi_k(x) \phi_l(y).

Note that the :math:`\frac{1}{2}` prefactor is applied internally and does not need to be included by the user.

What follows is an example for a Hubbard-Kanamori interaction corresponding to the default input variables
:f:var:`ULOC` = :code:`100`, :f:var:`UST` = :code:`10`, :f:var:`JH` = :f:var:`JX` = :f:var:`JP` = :code:`1` .

.. code-block:: text
   
   2 BANDS
   
   #ULOC
   1 d 1 u 1 d 1 u 100.0
   1 u 1 d 1 u 1 d 100.0
   2 d 2 u 2 d 2 u 100.0
   2 u 2 d 2 u 2 d 100.0
   
   #UST
   1 d 2 u 1 d 2 u 10.0
   1 u 2 d 1 u 2 d 10.0
   2 d 1 u 2 d 1 u 10.0
   2 u 1 d 2 u 1 d 10.0
   
   #UST-JH 
   1 u 2 u 1 u 2 u 9.0
   1 d 2 d 1 d 2 d 9.0
   2 d 1 d 2 d 1 d 9.0
   2 u 1 u 2 u 1 u 9.0
   
   #JX
   1 d 2 u 2 d 1 u 1.0
   1 u 2 d 2 u 1 d 1.0
   2 d 1 u 1 d 2 u 1.0
   2 u 1 d 1 u 2 d 1.0
   
   #JP
   1 d 1 u 2 d 2 u 1.0
   1 u 1 d 2 u 2 d 1.0
   2 d 2 u 1 d 1 u 1.0
   2 u 2 d 1 u 1 d 1.0


.. f:automodule::   ed_parse_umatrix
   :members: read_umatrix_file, save_umatrix_file, reset_umatrix, add_twobody_operator
