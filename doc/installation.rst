.. _edipack_install:

Install
##############################

|edipack| is available in the form of a static Fortran library
(`libedipack.a`) and the related Fortran module :f:mod:`EDIPACK`.
The release version includes additional modules to extend the
software functionalities: an inequivalent impurities extension |edipack2ineq|
and a shared dynamical library implementing the Fortran-C interface. 

A standard installation method is available through CMake.

An alternative approach is provided via Anaconda. 


From source
========================

Building
------------------------------

We assume that SciFortran_ and MPI_ have been correctly installed
and are available in the system. See related documentation. Note that
the installation of |edipack| closely follows the SciFortran_ template.

.. _SciFortran: https://github.com/SciFortran/SciFortran
.. _MPI: https://github.com/open-mpi/ompi


Clone the repo:

.. code-block:: bash
		
   git clone https://github.com/edipack/EDIpack EDIpack



Optionally define the fortran compiler:

.. code-block:: bash
		
   export FC=mpif90/gfortran/ifort


From the repository directory (`cd EDIpack`) make a standard
out-of-source CMake compilation:

GNU Make
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Using GNU `make` is the default CMake workflow, with widest version
support (CMake > 3.0). Note that parallel `make` execution is tested
and working.

.. code-block:: bash
		
   mkdir build 
   cd build  
   cmake .. 
   make -j


..
    Ninja
    ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

    Using `ninja` if a fortran-capable version of `ninja
    <https://ninja-build.org>`_ is available in your system (and CMake can
    take advantage of it), you can use it to build the library at lightning, multi-threaded, speed. 

    .. code-block:: bash
		    
       mkdir build    
       cd build  
       cmake -GNinja ..  
       ninja

The `CMake` compilation can be customized using the following
additional variables:   

* ``-DLONG_PREFIX``
  
  * **Scope:** Set install prefix. If ``True``, installs to
    ``~/opt/EDIpack/TAG/PLAT/BRANCH``; otherwise to
    ``~/opt/EDIpack/custom``.
  * **Value:** ``yes``/``True`` or ``no``/``False``
    (default: ``True``).

* ``-DCMAKE_INSTALL_PREFIX``

  * **Scope:** Specify a custom installation path. Overrides
    ``-DLONG_PREFIX``. No default subdirectories are added.
  * **Value:** User-defined path.

* ``-DUSE_MPI``

  * **Scope:** Enable MPI support at compile time.
  * **Value:** ``yes``/``True`` or ``no``/``False``
    (default: ``True``).

* ``-DWITH_INEQ``

  * **Scope:** Include the inequivalent impurity extension
    (``EDIPACK2INEQ``).
  * **Value:** ``yes``/``True`` or ``no``/``False``
    (default: ``True``).

* ``-DCMPLX_NORMAL``

  * **Scope:** Enable complex eigenvectors for the
    :f:var:`ED_MODE` = ``NORMAL`` case.
  * **Value:** ``yes``/``True`` or ``no``/``False``
    (default: ``False``).
    
* ``-DVERBOSE``

  * **Scope:** Enable verbose CMake output. Superseded by
    ``make VERBOSE=yes/no``.
  * **Value:** ``yes``/``True`` or ``no``/``False``
    (default: ``True``).

* ``-DBUILD_TYPE``

  * **Scope:** Set build mode flags.
  * **Value:** ``RELEASE``, ``DEBUG``, ``TESTING``,
    ``AGGRESSIVE`` (default: ``RELEASE``).


.. warning::
   
   :code:`BUILD_TYPE=AGGRESSIVE`  includes many deep level debug options which might not compile on some systems or breakdown compilation at linking step.



The default target builds either the main library and the C-bindings. A
specific building for each library is also available specifying the
required target. For user convenience at the end CMake configuration
step the following recap is printed:

.. code-block:: bash

   *Build edipack [Default]:  
   $ make -j [all/edipack, default=all]
   
   *Build C-bindings: 
   $ make edipack_cbindings
      
   *Install: 
   $ make [all/edipack/edipack_cbindings, default=all] install
   
   *Uninstall: 
   $ make uninstall
   
   *Build documenation: 
   $ make doc
   
   *Build and Runtest: 
   $ make test







   
Installing
------------------------------

System-wide installation is completed after the build step using

.. code-block:: bash

   make install

  
Please follow the instructions on the screen to complete installation on your environment.  
The library can be loaded using one of the following, automatically generated, files :  

*  A generated `environment module`_ , installed to
   `~/.modules.d/edipack/<PREFIX/default if LONG_PREFIX=T/F>`
  
* A generated `bash` script installed to `<PREFIX>/bin/configvars.sh`, to be sourced for permanent loading.

*  A generated `pkg-config`_ file to, installed to `~/.pkg-config.d/edipack.pc`  

.. _environment module: https://github.com/cea-hpc/modules
.. _pkg-config: https://github.com/freedesktop/pkg-config

For ease of use a specific and automatically generated recap message is printed after installation. 




Uninstalling
------------------------------




CMake does not officially provide uninstall procedures in the
generated makefiles. Hence EDIpack supplies a homebrew
method to remove the generated files by calling (from the relevant
build folder):

.. code-block:: bash
		
   make uninstall






   


From Anaconda
==============================

We provide Linux and MacOS packages for the `Anaconda <https://www.anaconda.com/>`_ 
distribution. To install the module, the virtual environment of choice should include
python 3.10 or later.

Once a command-line tool such as `conda <https://www.anaconda.com/>`_ or 
`mamba <https://mamba.readthedocs.io/en/latest/>`_ is installed, an environment 
using one of the available python version can be created, and then the EDIpack 
package can be installed:

.. code-block:: shell

   conda create -n edipack
   conda activate edipack
   conda install -c conda-forge -c edipack edipack


this installs a bundle of the :code:`scifor` and :code:`edipack` libraries. In order to compile a
fortran program linking the libraries, we provide  :code:`.pc` files which are readable via 
:code:`pkg-config`. If not present, the :code:`compilers` and :code:`pkg-config` conda
packages need to be installed

.. code-block:: shell

   conda install compilers
   conda install pkg-config
   
The inclusion and linking flag can then be obtained via 

.. code-block:: shell

   pkg-config --cflags edipack scifor
   pkg-config --libs   edipack scifor

