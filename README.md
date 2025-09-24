# EDIpack: a generic and interoperable Lanczos-based  exact diagonalization solver for Quantum Impurity problems 

[![TestSuite](https://img.shields.io/github/actions/workflow/status/edipack/EDIpack/PushWorkflow.yml?label=TestSuite&logo=Fortran&style=flat-square)](https://github.com/edipack/EDIpack/actions/workflows/PushWorkflow.yml) 
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://edipack.github.io/EDIpack/)
[![Anaconda-Server Badge](https://anaconda.org/edipack/edipack/badges/version.svg)](https://anaconda.org/edipack/edipack)

<!-- TO BE SETUP ASAP
[![Coverage]()]()
[![api docs](https://img.shields.io/static/v1?label=API&message=documentation&color=734f96&logo=read-the-docs&logoColor=white&style=flat-square)](https://qcmplab.github.io/DMFT_ED)
-->

This is the latest version of [EDIpack](https://github.com/edipack/EDIpack): a  Lanczos based method 
for the solution of generic Quantum Impurity problems,  exploiting distributed memory MPI parallelisation.
This version, aims to solve single-site, multi-orbital models, in either  *normal*, *superconducting* (s-wave) or *Spin-non-conserving* (e.g. with Spin-Orbit Coupling or in-plane magnetization) phases, including electron-phonons coupling. The code works at zero and low temperatures.   
See the EDIpack associated [publications](#reference)  

### Install 
*EDIpack* is available in the form of a static Fortran library (`libedipack.a`) and the related Fortran module `EDIPACK`.
The release version includes additional modules to extend the software functionalities: i) an inequivalent impurities extension `Edipack2ineq`
and ii) a shared dynamical library `edipack_cbindings.so` implementing the Fortran-C interface. 

A standard installation from source is available through `CMake`, via the standard out-of-source method. 

An alternative approach is provided via `Anaconda`. 

Detailed information can be found at [edipack.github.io/EDIpack/installation](https://edipack.github.io/EDIpack/installation.html)


### Documentation
All the informations about the structure of the library and its use can be found in the documenation at [edipack.github.io/EDIpack/](https://edipack.github.io/EDIpack/)  


### Use
In [Quickstart](https://edipack.github.io/EDIpack/quickstart/02_dmft.html) we illustrate the use and the capabilities of EDIpack as a solver for Dynamical Mean-Field Theory calculation. 



### Reference
A full overview of EDIpack can be found in the related publication which can be found at [SciPostPhysCodeb.58](https://scipost.org/10.21468/SciPostPhysCodeb.58). The paper source is also available at [EDIpack/EDIpackManuscript](https://github.com/EDIpack/EDIpackManuscript).  

Useful informations about the distributed memory parallel algorithms underlying the functioning of EDIpack and their benchmarks are presented in [j.cpc.2021.108261](https://doi.org/10.1016/j.cpc.2021.108261) (also freely available in the arXiv).   

**Should you use EDIpack or any of the derived packages, please consider citing these papers**:

```
@article{Crippa2025SPCa,
	author = {Lorenzo Crippa and Igor Krivenko and Samuele Giuli and Gabriele Bellomia and Alexander Kowalski and Francesco Petocchi and Alberto Scazzola and Markus Wallerberger and Giacomo Mazza and Luca de Medici and Giorgio Sangiovanni and Massimo Capone and Adriano Amaricci},
	date-added = {2025-09-24 16:21:06 +0200},
	date-modified = {2025-09-24 16:21:10 +0200},
	doi = {10.21468/SciPostPhysCodeb.58},
	journal = {SciPost Phys. Codebases},
	pages = {58},
	publisher = {SciPost},
	title = {{Next-generation EDIpack: A Lanczos-based package for quantum impurity models featuring general broken-symmetry phases, flexible bath topologies and multi-platform interoperability}},
	url = {https://scipost.org/10.21468/SciPostPhysCodeb.58},
	year = {2025},
	bdsk-url-1 = {https://scipost.org/10.21468/SciPostPhysCodeb.58},
	bdsk-url-2 = {https://doi.org/10.21468/SciPostPhysCodeb.58}}




@article{Crippa2025SPC,
	author = {Lorenzo Crippa and Igor Krivenko and Samuele Giuli and Gabriele Bellomia and Alexander Kowalski and Francesco Petocchi and Alberto Scazzola and Markus Wallerberger and Giacomo Mazza and Luca de Medici and Giorgio Sangiovanni and Massimo Capone and Adriano Amaricci},
	date-added = {2025-09-24 16:21:06 +0200},
	date-modified = {2025-09-24 16:21:10 +0200},
	doi = {10.21468/SciPostPhysCodeb.58-r5.3},
	journal = {SciPost Phys. Codebases},
	pages = {58-r5.3},
	publisher = {SciPost},
	title = {{Codebase release 5.3 for EDIpack}},
	url = {https://scipost.org/10.21468/SciPostPhysCodeb.58-r5.3},
	year = {2025},
	bdsk-url-1 = {https://scipost.org/10.21468/SciPostPhysCodeb.58-r5.3},
	bdsk-url-2 = {https://doi.org/10.21468/SciPostPhysCodeb.58-r5.3}}

```

and 

```
@article{amaricci2022CPC,
	author = {A. Amaricci and L. Crippa and A. Scazzola and F. Petocchi and G. Mazza and L. {de Medici} and M. Capone},
	doi = {https://doi.org/10.1016/j.cpc.2021.108261},
	issn = {0010-4655},
	journal = {Computer Physics Communications},
	keywords = {Exact diagonalization, Quantum impurity models, Strongly correlated electrons, Dynamical mean-field theory},
	pages = {108261},
	title = {EDIpack: A parallel exact diagonalization package for quantum impurity problems},
	url = {https://www.sciencedirect.com/science/article/pii/S0010465521003738},
	volume = {273},
	year = {2022},
	bdsk-url-1 = {https://www.sciencedirect.com/science/article/pii/S0010465521003738},
	bdsk-url-2 = {https://doi.org/10.1016/j.cpc.2021.108261}}
```


### Issues
If you encounter bugs, difficulties or have any other query please [file an issue](https://github.com/edipack/EDIpack/issues/new/choose) or send an email to [edipack@iom.cnr.it](mailto:edipack@iom.cnr.it).          

### Authors
EDIpack authors:   
[Adriano Amaricci](https://github.com/aamaricci)   
[Lorenzo Crippa](https://github.com/lcrippa)    
[Samuele Giuli](https://github.com/SamueleGiuli)    
[Gabriele Bellomia](https://github.com/beddalumia)    
[Igor Krivenko](https://github.com/krivenko)  
Alberto Scazzola   
Luca de Medici   
[Giacomo Mazza](https://github.com/GiacMazza)  
Francesco Petocchi  
[Massimo Capone](https://github.com/massimocapone)

Other important authors contributed to the development of the EDIpack ecosystem:   
[Alexander Kowalski](https://github.com/alexkowalski)  
[Markus Wallerberger](https://github.com/mwallerb)   
[Giorgio Sangiovanni](https://github.com/sangiova) 
