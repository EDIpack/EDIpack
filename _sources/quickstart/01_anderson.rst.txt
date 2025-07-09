.. _01_anderson:

Solving a Single Impurity Anderson model
=============================================

.. preferred-crossrefs::
   :ed_solve: f/ed_main/ed_solve

In this section we use the methods in |edipack|  to solve a simple
example of Anderson quantum impurity problem.

Looking forward for a DMFT application, here we consider the Bethe
lattice DOS :math:`\rho(x)=\frac{1}{2D}\sqrt{D^2-x^2}` and build the
corresponding non-interacting  Green's function :math:`G_0(z) =
\int_{\mathbb{R}}\rho(\epsilon)\left[ z -\epsilon \right]^{-1}`.

We construct a discretized bath by fitting such function on the
Matsubara frequencies  :math:`G_0(i\omega_n)` using the methods in
:ref:`fit`. Finally we input this bath into the :f:func:`ed_solve`
solver of  |edipack| in presence of local interaction on the impurity.  

Source code
------------------------------

The initialization of the code is:

.. code-block:: fortran
		
   program lancED
      USE EDIPACK
      USE SCIFOR
      implicit none
      ! Bethe lattice half-bandwidth = energy unit
      real(8),parameter                       :: D=1d0
      ! Bath size and Nso=Nspin*Norb (here =1)
      integer                                 :: Nb,Nso
      ! User bath, allocatable see below
      real(8),allocatable                     :: Bath(:)
      ! Non-interacting Bethe lattice Green's function (naming
      ! convention will be clear in the following section)
      complex(8),allocatable,dimension(:,:,:) :: Weiss

      !> READ THE input using EDIpack procedure: 
      call ed_read_input('inputED.conf')

      
where we load both the |edipack| and SciFortran_ libraries through
their main module :f:mod:`EDIPACK` and :f:mod:`SCIFOR`. We also define
some local variables and proceed with reading the input file
:code:`"inputED.conf"` using the function
:f:func:`ed_read_input`.

Next we construct the non-interacting Green's
function, using procedures available in SciFortran_:

.. _SciFortran: https://github.com/SciFortran/SciFortran


.. code-block:: fortran
		
   !> Get the Bethe lattice non-interacting Matsubara GF as a guess for the bath 
   allocate(Weiss(Nso,Nso,Lmats))
   call bethe_guess_g0(Weiss(1,1,:),D,beta,hloc=0d0)



Then we initialize the solver. This step requires the user to pass
the user bath as a rank-1 double precision array of a given size. The
correct size of the bath array is evaluated internally by the
|edipack| code through the function
:f:func:`ed_get_bath_dimension`.

.. code-block:: fortran

   !> Init solver: 
   Nb=ed_get_bath_dimension()
   allocate(bath(Nb))
   call ed_init_solver(bath)


Upon initialization the bath is guessed from a flat distribution
centered around zero and with half-width :f:var:`ed_hw_bath`. Here we
update the bath optimizing it against the non-interacting Bethe
lattice Green's function:

.. code-block:: fortran

   !> Fit the bath against G0 guess: the outcome is a bath discretizing the Bethe DOS.
   call ed_chi2_fitgf(Weiss,bath,ispin=1,iorb=1)


We are now ready to solve the quantum impurity problem for a given set
of parameters specified in the input file (see below)


.. code-block:: fortran

   !> Solve SIAM with this given bath
   call ed_solve(bath)


.. raw:: html

   <hr>


Results
------------------------------

Here is a snapshot of the results obtained for :math:`U=1.0, 10.0`.

   
.. image:: 01_anderson_fig.svg
   :class: with-border
   :width: 800px
	     
In the top panel we report the impurity spectral functions :math:`-\Im
G^{im}(\omega)/\pi` compared to the  Bethe density of states (filled
curve). In the bottom panel we show the real part of the
impurity real-axis self-energy functions :math:`\Sigma(\omega)`
near the Fermi level :math:`\omega=0`. The linear fit :math:`y =
A\omega` gives a direct estimate of the derivatives :math:`A\simeq
\tfrac{\partial\Re\Sigma}{\partial\omega}_{|_{\omega\rightarrow 0}}` and thus of the quasi-particle
renormalization constants :math:`Z=\left( 1 - A \right)^{-1}` as
reported in the legend.

As a direct comparison we report also the values of :math:`Z` 
estimated from the Matsubara axis using the relation 
:math:`\frac{\Im\Sigma(i\omega_n}{\omega_n}_{|_{\omega_n\rightarrow
0}}= \frac{1}{\pi}\int_{\mathbb R}d\epsilon
\frac{\Re\Sigma(\epsilon)}{\epsilon^2}=
\frac{\partial\Re\Sigma}{\partial\omega}_{|_{\omega\rightarrow
0}}`.

We obtained:  :math:`Z=0.74` and :math:`Z=0.002`, respectively, for :math:`U=1.0`
and  :math:`U=10.0`.





   

   
.. raw:: html

   <hr>


Here is an example of input file used in actual calculations:

.. code-block::

   NORB=1                                        !Number of impurity orbitals (max 5).
   NBATH=9                                       !Number of bath sites
   ULOC=1d0				       
   ED_MODE=normal                                !Flag to set ED type
   BATH_TYPE=normal                              !flag to set bath type
   BETA=1000.000000000                           !Inverse temperature, at T=0 is used as a IR cut-off.
   ED_VERBOSE=3                                  !Verbosity level: 0=almost nothing --> 5:all. Really: all
   LMATS=4096                                    !Number of Matsubara frequencies.
   LFIT=2048
   EPS=1.000000000E-02                           !Broadening on the real-axis.
   CG_FTOL=1.000000000E-10                       !Conjugate-Gradient tolerance.
   CG_NITER=2048                                 !Max. number of Conjugate-Gradient iterations.
   ED_TWIN=T
   LANC_NGFITER=500                              !Number of Lanczos iteration in GF determination. Number of momenta.
   ED_HW_BATH=1.000000000                        !half-bandwidth for the bath initialization: flat in -ed_hw_bath:ed_hw_bath


   
