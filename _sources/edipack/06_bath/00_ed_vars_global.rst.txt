Bath class definition
=======================


In :f:mod:`ED_VARS_GLOBAL`  a suitable data structure
:f:var:`effective_bath` is defined. This type  encapsulates the
discretized effective bath of the quantum impurity problem. A
global instance of the bath :f:var:`dmft_bath` is then shared through
the code to store the current state of the user provided bath.  



.. f:automodule::   ed_vars_global
   :members: effective_bath_component,effective_bath,dmft_bath



