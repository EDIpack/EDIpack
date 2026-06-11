.. _fit:


:math:`\chi^2` Fit
---------------------------------------------------------------------------


In :f:mod:`ED_BATH_FIT` we perform an optimisation of the user bath by minimizing
the projection of a user supplied function with respect to the
corresponding quantum impurity problem, bath dependent, one
:math:`\min_{\vec{b}}\chi^2(\vec{b})` with:

.. math::

   \chi^2(\vec{b}) = || F(z) - F^{And}(z;\vec{b}) ||
   
:math:`F(z)`  is either the Weiss field :math:`{\cal G}_0` or
the hybridization function :math:`\Delta` as evaluated by the user. 
:math:`F^{And}(z;\vec{b})` is, respectively, the quantum impurity
non-interacting Green's function
:math:`G^{And}_0(z;\vec{b})=[z+\mu-h_0-\Delta(z;\vec{b})]^{-1}`  or the
hybrization function :math:`\Delta(z,c)=\sum_p
V_p[z-h^p]^{-1}V_p`. Finally :math:`\vec{b}` is the array containing the
discretized bath parameters. 

The minimization with respect to :math:`\vec{b}` is performed using conjugate gradient algorithm. 

.. preferred-crossrefs::
   :sf_linalg/operator(.x.): f/sf_sparse_array_csc/operator(.x.)

.. f:automodule::   ed_bath_fit
