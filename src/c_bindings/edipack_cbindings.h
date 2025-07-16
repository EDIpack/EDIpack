#ifndef EDIPACK_CBINDINGS
#define EDIPACK_CBINDINGS

#include <stdint.h>
#include <stdbool.h>
#include <complex.h>

#ifdef __cplusplus
extern "C" {
#endif


/**
*This flag is `1` if real-space DMFT support enabled, `0` otherwhise
*/ 
extern int has_ineq;

/*!
\rst
Inverse temperature. 
Interfaces to :f:var:`beta`.
\endrst
*/ 
extern double beta;

/*!
\rst
Convergence threshold. 
Interfaces to :f:var:`dmft_error`.
\endrst
*/ 
extern double dmft_error;

/*!
\rst
Twin sector simmetry flag. 
Interfaces to :f:var:`ed_twin`.
\endrst
*/ 
extern bool ed_twin;

/*!
\rst
Flag setting orbital-resolved sectorization of Fock space Hamiltonian. 
Interfaces to :f:var:`ed_total_ud`.
\endrst
*/ 
extern bool ed_total_ud;

/*!
\rst
Real axis broadening. 
Interfaces to :f:var:`eps`.
\endrst
*/ 
extern double eps;

/*!
\rst
Hund's coupling. 
Interfaces to :f:var:`Jh`.
\endrst
*/ 
extern double Jh;

/*!
\rst
Pair-hopping interaction strength. 
Interfaces to :f:var:`Jp`.
\endrst
*/ 
extern double Jp;

/*!
\rst
Spin-exchange interaction strength. 
Interfaces to :f:var:`Jx`.
\endrst
*/ 
extern double Jx;

/*!
\rst
Number of frequencies used in :math:`\chi^{2}` fit. 
Interfaces to :f:var:`Lfit`.
\endrst
*/ 
extern int Lfit;

/*!
\rst
Number of Matsubara frequencies. 
Interfaces to :f:var:`Lmats`.
\endrst
*/ 
extern int Lmats;

/*!
\rst
Number of real frequencies. 
Interfaces to :f:var:`Lreal`.
\endrst
*/ 
extern int Lreal;

/*!
\rst
Number of points in Probability Distribution Function lattice for phonons. 
Interfaces to :f:var:`Lpos`.
\endrst
*/ 
extern int Lpos;

/*!
\rst
Number of imaginary-time points. 
Interfaces to :f:var:`Ltau`.
\endrst
*/ 
extern int Ltau;

/*!
\rst
Number of bath sites. 
Interfaces to :f:var:`Nbath`.
\endrst
*/ 
extern int Nbath;

/*!
\rst
Number of DMFT loops. 
Interfaces to :f:var:`Nloop`.
\endrst
*/ 
extern int Nloop;

/*!
\rst
Number of orbitals. 
Interfaces to :f:var:`Norb`.
\endrst
*/ 
extern int Norb;

/*!
\rst
Number of phonons. 
Interfaces to :f:var:`Nph`.
\endrst
*/ 
extern int Nph;

/*!
\rst
Occupation value for fixed-density calculations. 
Interfaces to :f:var:`nread`.
\endrst
*/ 
extern double nread;

/*!
\rst
Number of spins. 
Interfaces to :f:var:`Nspin`.
\endrst
*/ 
extern int Nspin;

/*!
\rst
Number of iterations under convergence threshod. 
Interfaces to :f:var:`Nsuccess`.
\endrst
*/ 
extern int Nsuccess;

/*!
\rst
Symmetry-breaking field. 
Interfaces to :f:var:`sb_field`.
\endrst
*/ 
extern double sb_field;

/*!
\rst
Hubbard interaction strength. 
Interfaces to :f:var:`Uloc`.
\endrst
*/ 
extern double Uloc[5];

/*!
\rst
Density-density Kanamori interaction strength. 
Interfaces to :f:var:`Ust`.
\endrst
*/ 
extern double Ust;

/*!
\rst
Upper bound on real axis frequencies. 
Interfaces to :f:var:`wfin`.
\endrst
*/ 
extern double wfin;

/*!
\rst
Upper bound on real axis frequencies. 
Interfaces to :f:var:`wini`.
\endrst
*/ 
extern double wini;

/*!
\rst
Maximum of the x-range for the local lattice probability distribution function. 
Interfaces to :f:var:`xmax`. 
\endrst
*/ 
extern double xmax;

/*!
\rst
Minimum of the x-range for the local lattice probability distribution function. 
Interfaces to :f:var:`xmin`. 
\endrst
*/ 
extern double xmin;


/*!
\rst
Chemical potential. 
Interfaces to :f:var:`xmu`.
\endrst
*/ 
extern double xmu;




/*!
\rst
This function reads the input. Interfaces to :f:func:`ed_read_input`
\endrst
*
* @param instr: the input file name
*/
void read_input(char *instr);


/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/ed_aux_funx/ed_set_Hloc`.
Rank-2 array variant for single-site DMFT
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
*/
void ed_set_Hloc_single_N2(std::complex<double> *Hloc, 
                           int64_t *d
                           );

/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/ed_aux_funx/ed_set_Hloc`.
Rank-4 array variant for single-site DMFT.
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
*/       
void ed_set_Hloc_single_N4(std::complex<double> *Hloc, 
                           int64_t *d
                           );
                           
                           
/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/e2i_aux_funx/ed_set_Hloc`.
Rank-2 array variant for real-space DMFT.
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
* @param Nlat: number of inequivalent sites
*/                      
void ed_set_Hloc_lattice_N2(std::complex<double> *Hloc,
                            int64_t *d,
                            int Nlat
                            );

/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/e2i_aux_funx/ed_set_Hloc`.
Rank-3 array variant for real-space DMFT.
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
* @param Nlat: number of inequivalent sites
*/ 
void ed_set_Hloc_lattice_N3(std::complex<double> *Hloc,
                            int64_t *d,
                            int Nlat
                            );

/*!
\rst
This function sets the local Hamiltonian.
Interfaces to :f:func:`f/e2i_aux_funx/ed_set_Hloc`.
Rank-5 array variant for real-space DMFT.
\endrst
*
* @param Hloc: the local Hamiltonian
* @param d: array of dimensions of the local Hamiltonian
* @param Nlat: number of inequivalent sites
*/ 
void ed_set_Hloc_lattice_N5(std::complex<double> *Hloc, 
                            int64_t *d, 
                            int Nlat
                            );



/*!
\rst
This function gets the dimension of the user-accessible bath array.
Interfaces to :f:func:`get_bath_dimension`.
\endrst
*
* @return the dimension of the bath array
*/ 
int  get_bath_dimension_direct(void);

/*!
\rst
This function gets the dimension of the user-accessible bath array.
Specific for replica and general baths.
Interfaces to :f:func:`get_bath_dimension`.
\endrst
*
* @param Nsym: the number of components of the replica matrix linear combination
* @return the dimension of the bath array
*/ 
int  get_bath_dimension_symmetries(int Nsym);

/*!
\rst
This function updates a variable (usually the chemical potential) trying to 
achieve a desired density set by :f:var:`nread`.
Interfaces to :f:func:`ed_search-variable`.
\endrst
*
* @param var: the variable to be adjusted
* @param ntmp: the density value at a given iteration
* @param converged: the DMFT loop convergence status.
*/ 
void search_variable(double *var,
                    double *ntmp, 
                    int64_t *converged
                    );


/*!
\rst
This function sets the replica bath H. 
Interfaces to :f:func:`f/ed_bath_replica/set_hreplica`.
Rank-3 array case for single-site DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the replica H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the array linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                                 
void init_Hreplica_symmetries_d3(std::complex<double> *Hvec, 
                                 int64_t *d_hvec, 
                                 double *lambdavec, 
                                 int64_t *d_lambdavec
                                 );

/*!
\rst
This function sets the replica bath H. 
Interfaces to :f:func:`f/ed_bath_replica/set_hreplica`.
Rank-5 array case for single-site DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the replica H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the array linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/ 
void init_Hreplica_symmetries_d5(std::complex<double> *Hvec, 
                                 int64_t *d_hvec, 
                                 double *lambdavec, 
                                 int64_t *d_lambdavec
                                 );
                                                                  
/*!
\rst
This function sets the general bath H. 
Interfaces to :f:func:`f/ed_bath_replica/set_hgeneral`.
Rank-3 array case for single-site DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the general H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the array linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                               
void init_Hgeneral_symmetries_d3(std::complex<double> *Hvec, 
                                 int64_t *d_hvec, 
                                 double *lambdavec, 
                                 int64_t *d_lambdavec
                                 );                              
                                 
/*!
\rst
This function sets the general bath H. 
Interfaces to :f:func:`f/ed_bath_replica/set_hgeneral`.
Rank-5 array case for single-site DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the general H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the array linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                             
void init_Hgeneral_symmetries_d5(std::complex<double> *Hvec, 
                                 int64_t *d_hvec, 
                                 double *lambdavec, 
                                 int64_t *d_lambdavec
                                 );

/*!
\rst
This function sets the replica bath H. 
Interfaces to :f:func:`f/e2i_bath_replica/set_hreplica`.
Rank-3 array case for real-space DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the replica H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the array linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                                            
void init_Hreplica_symmetries_lattice_d3(std::complex<double> *Hvec, 
                                         int64_t *d_hvec, 
                                         double *lambdavec, 
                                         int64_t *d_lambdavec
                                         );

/*!
\rst
This function sets the replica bath H. 
Interfaces to :f:func:`f/e2i_bath_replica/set_hreplica`.
Rank-5 array case for real-space DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the replica H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the array linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/   
void init_Hreplica_symmetries_lattice_d5(std::complex<double> *Hvec, 
                                         int64_t *d_hvec, 
                                         double *lambdavec, 
                                         int64_t *d_lambdavec
                                         );

/*!                                         
\rst
This function sets the general bath H. 
Interfaces to :f:func:`f/e2i_bath_replica/set_hgeneral`.
Rank-3 array case for real-space DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the general H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the array linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/                                       
void init_Hgeneral_symmetries_lattice_d3(std::complex<double> *Hvec, 
                                         int64_t *d_hvec, 
                                         double *lambdavec, 
                                         int64_t *d_lambdavec
                                         );

                                                                                 
/*!
\rst
This function sets the general bath H. 
Interfaces to :f:func:`f/e2i_bath_replica/set_hgeneral`.
Rank-5 array case for real-space DMFT.
\endrst
*
* @param Hvec: array of matrices summing up to the general H
* @param d_hvec: dimensions of the array of matrices
* @param lambdavec: array of coefficients of the array linear combination
* @param d_lambdavec: dimensions of the array of coefficients
*/ 
void init_Hgeneral_symmetries_lattice_d5(std::complex<double> *Hvec, 
                                         int64_t *d_hvec, 
                                         double *lambdavec, 
                                         int64_t *d_lambdavec
                                         );

/*!                                         
\rst
This function breaks the bath symmetry. 
Interfaces to :f:func:`f/ed_bath_user/break_symmetry_bath`. 
Single-site DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param field: symmetry breaking field 
* @param sgn: sign of the symmetry breaking field 
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/   
void break_symmetry_bath_site(double *bath, 
                              int64_t *dim_bath, 
                              double field, 
                              double sgn, 
                              int sav
                              );
                              
/*!                                         
\rst
This function breaks the bath symmetry. 
Interfaces to :f:func:`f/e2i_bath_user/break_symmetry_bath`. 
Real-space DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath and sign array
* @param field: symmetry breaking field 
* @param sgn: array of signs of the symmetry breaking field 
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/                               
void break_symmetry_bath_ineq(double *bath, 
                              int64_t *dim_bath, 
                              double field, 
                              double *sgn, 
                              int sav
                              );

/*!                                         
\rst
This function enforces a paramagnetic bath
Interfaces to :f:func:`f/ed_bath_user/spin_symmetrize_bath`. 
Single-site DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath and sign array
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/ 
void spin_symmetrize_bath_site(double *bath, 
                               int64_t *dim_bath, 
                               int sav
                               );
                               
/*!                                         
\rst
This function enforces a paramagnetic bath
Interfaces to :f:func:`f/e2i_bath_user/spin_symmetrize_bath`. 
Real-space DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/ 
void spin_symmetrize_bath_ineq(double *bath, 
                               int64_t *dim_bath, 
                               int sav
                               );

/*!                                         
\rst
This function enforces an orbital-symmetric bath. 
Interfaces to :f:func:`f/ed_bath_user/orb_symmetrize_bath`. 
Single-site DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param orb1: first orbital to symmetrize
* @param orb2: second orbital to symmetrize
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/ 
void orb_symmetrize_bath_site(double *bath, 
                              int64_t *dim_bath, 
                              int orb1, 
                              int orb2, 
                              int sav
                              );
                              
/*!                                         
\rst
This function enforces an orbital-symmetric bath. 
Interfaces to :f:func:`f/e2i_bath_user/orb_symmetrize_bath`. 
Real-space DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param orb1: first orbital to symmetrize
* @param orb2: second orbital to symmetrize
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/                               
void orb_symmetrize_bath_ineq(double *bath, 
                              int64_t *dim_bath, 
                              int orb1, 
                              int orb2, 
                              int sav
                              );

/*!                                         
\rst
This function enforces an orbita-symmetric bath. 
Interfaces to :f:func:`f/ed_bath_user/orb_symmetrize_bath`. 
Single-site DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param indx: index of the orbital to copy
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/ 
void orb_equality_bath_site(double *bath, 
                            int64_t *dim_bath, 
                            int indx, 
                            int sav
                            );
                            
/*!                                         
\rst
This function enforces an orbita-symmetric bath. 
Interfaces to :f:func:`f/e2i_bath_user/orb_symmetrize_bath`. 
Real-space DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param indx: index of the orbital to copy
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/ 
void orb_equality_bath_ineq(double *bath, 
                            int64_t *dim_bath, 
                            int indx, 
                            int sav
                            );

/*!                                         
\rst
This function enforces a particle-hole symmetric bath.
Interfaces to :f:func:`f/ed_bath_user/ph_symmetrize_bath`. 
Single-site DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/ 
void ph_symmetrize_bath_site(double *bath, 
                             int64_t *dim_bath, 
                             int sav
                             );
                             
/*!                                         
\rst
This function enforces a particle-hole symmetric bath.
Interfaces to :f:func:`f/e2i_bath_user/ph_symmetrize_bath`. 
Real-space DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param sav: integer flag to save the output bath. `1` if true, `0` if false.
*/                              
void ph_symmetrize_bath_ineq(double *bath, 
                             int64_t *dim_bath, 
                             int sav
                             );

/*!                                         
\rst
This function saves the bath array in a properly formatted file.
Interfaces to :f:func:`f/ed_bath_user/save_array_as_bath`. 
Single-site DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
*/
void save_array_as_bath_site(double *bath, 
                             int64_t *dim_bath
                             );
                             
/*!                                         
\rst
This function saves the bath array in a properly formatted file.
Interfaces to :f:func:`f/e2i_bath_user/save_array_as_bath`. 
Real-space DMFT variant.
\endrst
*
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
*/                          
void save_array_as_bath_ineq(double *bath, 
                             int64_t *dim_bath
                             );




/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/ed_bath_fit/ed_chi2_fitgf`. 
Single-site DMFT variant for :f:var:`ED_MODE` = :code:`NORMAL/NONSU2` and rank-3 arrays.
\endrst
*
* @param g: function to fit
* @param dim_g: dimensions of the function to fit
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
* @param iorb: orbital index. If `0`, the fortran function is called without `iorb`
* @param fmpi: integer flag to enable mpi. `1` for `True`, `0` for `False`
*/ 
void chi2_fitgf_single_normal_n3(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin, 
                                int iorb, 
                                int fmpi
                                );

/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/ed_bath_fit/ed_chi2_fitgf`. 
Single-site DMFT variant for :f:var:`ED_MODE` = :code:`NORMAL/NONSU2` and rank-5 arrays.
\endrst
*
* @param g: function to fit
* @param dim_g: dimensions of the function to fit
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
* @param iorb: orbital index. If `0`, the fortran function is called without `iorb`
* @param fmpi: integer flag to enable mpi. `1` for `True`, `0` for `False`
*/                                 
void chi2_fitgf_single_normal_n5(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin, 
                                int iorb, 
                                int fmpi
                                );

/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/ed_bath_fit/ed_chi2_fitgf`. 
Single-site DMFT variant for :f:var:`ED_MODE` = :code:`SUPERC` and rank-3 arrays.
\endrst
*
* @param g: function to fit, normal component
* @param dim_g: dimensions of the function to fit, normal component
* @param f: function to fit, anomalous component
* @param dim_f: dimensions of the function to fit, anomalous component
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
* @param iorb: orbital index. If `0`, the fortran function is called without `iorb`
* @param fmpi: integer flag to enable mpi. `1` for `True`, `0` for `False`
*/                                 
void chi2_fitgf_single_superc_n3(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin, 
                                 int iorb, 
                                 int fmpi
                                 );
                                 
                                 
/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/ed_bath_fit/ed_chi2_fitgf`. 
Single-site DMFT variant for :f:var:`ED_MODE` = :code:`SUPERC` and rank-5 arrays.
\endrst
*
* @param g: function to fit, normal component
* @param dim_g: dimensions of the function to fit, normal component
* @param f: function to fit, anomalous component
* @param dim_f: dimensions of the function to fit, anomalous component
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
* @param iorb: orbital index. If `0`, the fortran function is called without `iorb`
* @param fmpi: integer flag to enable mpi. `1` for `True`, `0` for `False`
*/                                   
void chi2_fitgf_single_superc_n5(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin, 
                                 int iorb, 
                                 int fmpi
                                 );

/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/e2i_bath_fit/ed_chi2_fitgf`. 
Real-space DMFT variant for :f:var:`ED_MODE` = :code:`NORMAL/NONSU2` and rank-3 arrays.
\endrst
*
* @param g: function to fit
* @param dim_g: dimensions of the function to fit
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
*/ 
void chi2_fitgf_lattice_normal_n3(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin
                                );
                              
/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/e2i_bath_fit/ed_chi2_fitgf`. 
Real-space DMFT variant for :f:var:`ED_MODE` = :code:`NORMAL/NONSU2` and rank-4 arrays.
\endrst
*
* @param g: function to fit
* @param dim_g: dimensions of the function to fit
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
*/                               
void chi2_fitgf_lattice_normal_n4(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin
                                );
                                
/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/e2i_bath_fit/ed_chi2_fitgf`. 
Real-space DMFT variant for :f:var:`ED_MODE` = :code:`NORMAL/NONSU2` and rank-6 arrays.
\endrst
*
* @param g: function to fit
* @param dim_g: dimensions of the function to fit
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
*/                                 
void chi2_fitgf_lattice_normal_n6(std::complex<double> *g, 
                                int64_t *dim_g, 
                                double *bath, 
                                int64_t *dim_bath, 
                                int ispin
                                );


/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/e2i_bath_fit/ed_chi2_fitgf`. 
Real-space DMFT variant for :f:var:`ED_MODE` = :code:`SUPERC` and rank-3 arrays.
\endrst
*
* @param g: function to fit
* @param dim_g: dimensions of the function to fit
* @param f: function to fit, anomalous component
* @param dim_f: dimensions of the function to fit, anomalous component
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
*/ 
void chi2_fitgf_lattice_superc_n3(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin
                                 );
                                 
/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/e2i_bath_fit/ed_chi2_fitgf`. 
Real-space DMFT variant for :f:var:`ED_MODE` = :code:`SUPERC` and rank-4 arrays.
\endrst
*
* @param g: function to fit
* @param dim_g: dimensions of the function to fit
* @param f: function to fit, anomalous component
* @param dim_f: dimensions of the function to fit, anomalous component
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
*/                                  
void chi2_fitgf_lattice_superc_n4(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin
                                 );
                                 
                                 
/*!                                         
\rst
This function fits the Weiss field or hybridization function.
Interfaces to :f:func:`f/e2i_bath_fit/ed_chi2_fitgf`. 
Real-space DMFT variant for :f:var:`ED_MODE` = :code:`SUPERC` and rank-6 arrays.
\endrst
*
* @param g: function to fit
* @param dim_g: dimensions of the function to fit
* @param f: function to fit, anomalous component
* @param dim_f: dimensions of the function to fit, anomalous component
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param ispin: spin index
*/                                  
void chi2_fitgf_lattice_superc_n6(std::complex<double> *g, 
                                 int64_t *dim_g, 
                                 std::complex<double> *f, 
                                 int64_t *dim_f, 
                                 double *bath, 
                                 int64_t *dim_bath, 
                                 int ispin
                                 );
                                 

/*!                                         
\rst
This function obtains the density.
Interfaces to :f:func:`f/ed_io/ed_get_dens`. 
Single-site DMFT version
\endrst
*
* @param self: orbitally-resolved density array
*/                                 
void ed_get_dens_n1(std::complex<double> *self);

/*!                                         
\rst
This function obtains the density.
Interfaces to :f:func:`f/e2i_io/ed_get_dens`. 
Real-space DMFT version
\endrst
*
* @param self: density array
* @param Nlat: number of inequivalent sites
*/   
void ed_get_dens_n2(std::complex<double> *self, int Nlat);

/*!                                         
\rst
This function obtains the magnetization.
Interfaces to :f:func:`f/ed_io/ed_get_mag`. 
Single-site DMFT version
\endrst
*
* @param self: magnetization array
*/  
void ed_get_mag_n2(std::complex<double> *self);

/*!                                         
\rst
This function obtains the magnetization.
Interfaces to :f:func:`f/e2i_io/ed_get_mag`. 
Real-space DMFT version
\endrst
*
* @param self: magnetization array
* @param Nlat: number of inequivalent sites
*/  
void ed_get_mag_n3(std::complex<double> *self, int Nlat);


/*!                                         
\rst
This function obtains the double-occupation.
Interfaces to :f:func:`f/ed_io/ed_get_docc`. 
Single-site DMFT version
\endrst
*
* @param self: double-occupation array
*/  
void ed_get_docc_n1(std::complex<double> *self);

/*!                                         
\rst
This function obtains the double occupation.
Interfaces to :f:func:`f/e2i_io/ed_get_docc`. 
Real-space DMFT version
\endrst
*
* @param self: double-occupation array
* @param Nlat: number of inequivalent sites
*/  
void ed_get_docc_n2(std::complex<double> *self, int Nlat);

/*!                                         
\rst
This function obtains the superconductive order parameter modulus.
Interfaces to :f:func:`f/ed_io/ed_get_phi`. 
Single-site DMFT version
\endrst
*
* @param self: superconductive order parameter modulus array
*/  
void ed_get_phisc_n2(std::complex<double> *self);

/*!                                         
\rst
This function obtains the superconductive order parameter modulus.
Interfaces to :f:func:`f/e2i_io/ed_get_phi`. 
Real-space DMFT version
\endrst
*
* @param self: superconductive order parameter modulus array
* @param Nlat: number of inequivalent sites
*/  
void ed_get_phisc_n3(std::complex<double> *self, int Nlat);



/*!                                         
\rst
This function obtains the superconductive order parameter argument.
Interfaces to :f:func:`f/ed_io/ed_get_argphi`. 
Single-site DMFT version
\endrst
*
* @param self: superconductive order parameter argument array
*/  
void ed_get_argsc_n2(std::complex<double> *self);

/*!                                         
\rst
This function obtains the superconductive order parameter argument.
Interfaces to :f:func:`f/e2i_io/ed_get_argphi`. 
Real-space DMFT version
\endrst
*
* @param self: superconductive order parameter argument array
* @param Nlat: number of inequivalent sites
*/  
void ed_get_argsc_n3(std::complex<double> *self, int Nlat);

  
  
/*!                                         
\rst
This function obtains the local energy.
Interfaces to :f:func:`f/ed_io/ed_get_eimp`. 
Single-site DMFT version
\endrst
*
* @param self: energy array
*/  
void ed_get_eimp_n1(std::complex<double> *self);

/*!                                         
\rst
This function obtains the local energy.
Interfaces to :f:func:`f/e2i_io/ed_get_eimp`. 
Real-space DMFT version
\endrst
*
* @param self: energy array
* @param Nlat: number of inequivalent sites
*/  
void ed_get_eimp_n2(std::complex<double> *self, int Nlat);


/*!
\rst
This function obtains the self-energy.
Interfaces to :f:func:`f/ed_io/ed_get_sigma` for a rank-3 array.
Single-site DMFT variant.
\endrst
* @param self: the self-energy array
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/   
void get_sigma_site_n3(std::complex<double> *self, 
                       int axis, 
                       int typ, 
                       std::complex<double> *zeta, 
                       int dz, 
                       int zflag
                       );
                       
/*!
\rst
This function obtains the self-energy.
Interfaces to :f:func:`f/ed_io/ed_get_sigma` for a rank-5 array.
Single-site DMFT variant.
\endrst
* @param self: the self-energy array
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/                         
void get_sigma_site_n5(std::complex<double> *self, 
                       int axis, 
                       int typ, 
                       std::complex<double> *zeta, 
                       int dz, 
                       int zflag
                       );

/*!
\rst
This function obtains the self-energy.
Interfaces to :f:func:`f/e2i_io/ed_get_sigma` for a rank-3 array.
Real-space DMFT variant.
\endrst
* @param self: the self-energy array
* @param Nineq: numer of inequivalent sites
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/                         
void get_sigma_lattice_n3(std::complex<double> *self, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );
                          
/*!
\rst
This function obtains the self-energy.
Interfaces to :f:func:`f/e2i_io/ed_get_sigma` for a rank-4 array.
Real-space DMFT variant.
\endrst
* @param self: the self-energy array
* @param Nineq: numer of inequivalent sites
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/                           
void get_sigma_lattice_n4(std::complex<double> *self, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );
                          
/*!
\rst
This function obtains the self-energy.
Interfaces to :f:func:`f/e2i_io/ed_get_sigma` for a rank-6 array.
Real-space DMFT variant.
\endrst
* @param self: the self-energy array
* @param Nineq: numer of inequivalent sites
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/                            
void get_sigma_lattice_n6(std::complex<double> *self, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );

/*!
\rst
This function obtains the impurity Green's function.
Interfaces to :f:func:`f/ed_io/ed_get_gimp` for a rank-3 array.
Single-site DMFT variant.
\endrst
* @param gimp: the impurity Green's function array
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/   
void get_gimp_site_n3(std::complex<double> *gimp, 
                       int axis, 
                       int typ, 
                       std::complex<double> *zeta, 
                       int dz, 
                       int zflag
                       );
                       
                       
/*!
\rst
This function obtains the impurity Green's function.
Interfaces to :f:func:`f/ed_io/ed_get_gimp` for a rank-5 array.
Single-site DMFT variant.
\endrst
* @param gimp: the impurity Green's function array
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/                          
void get_gimp_site_n5(std::complex<double> *gimp, 
                       int axis, 
                       int typ, 
                       std::complex<double> *zeta, 
                       int dz, 
                       int zflag
                       );

/*!
\rst
This function obtains the impurity Green's function.
Interfaces to :f:func:`f/e2i_io/ed_get_gimp` for a rank-3 array.
Real-space DMFT variant.
\endrst
* @param gimp: the impurity Green's function array
* @param Nineq: numer of inequivalent sites
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/                          
void get_gimp_lattice_n3(std::complex<double> *gimp, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );
                          
/*!
\rst
This function obtains the impurity Green's function.
Interfaces to :f:func:`f/e2i_io/ed_get_gimp` for a rank-4 array.
Real-space DMFT variant.
\endrst
* @param gimp: the impurity Green's function array
* @param Nineq: numer of inequivalent sites
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/                             
void get_gimp_lattice_n4(std::complex<double> *gimp, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );
                          
/*!
\rst
This function obtains the impurity Green's function.
Interfaces to :f:func:`f/e2i_io/ed_get_gimp` for a rank-6 array.
Real-space DMFT variant.
\endrst
* @param gimp: the impurity Green's function array
* @param Nineq: numer of inequivalent sites
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
* @param zeta: frequency array
* @param dz: dimension of the frequency array
* @param zflag: flag to set (`1`) or not (`0`) calculation with internal frequency array
*/                             
void get_gimp_lattice_n6(std::complex<double> *gimp, 
                          int Nineq, 
                          int axis, 
                          int typ, 
                          std::complex<double> *zeta, 
                          int dz, 
                          int zflag
                          );


/*!
\rst
This function obtains the Weiss field.
Interfaces to :f:func:`f/ed_io/ed_get_g0and` for a rank-3 array.
Single-site DMFT variant.
\endrst
* @param warray: Array of frequencies
* @param dim_warray: Dimension of the array of frequencies
* @param bath: User-accessible bath array
* @param dim_bath: Dimension of the bath array
* @param G0and: the Weiss field array
* @param dim_g0and: the dimensions of the Weiss field array
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
*
*/ 
void get_g0and_n3(std::complex<double> *warray, 
                  int64_t *dim_warray, 
                  double *bath, 
                  int dim_bath, 
                  std::complex<double> *G0and, 
                  int64_t *dim_g0and, 
                  char *axis, 
                  char *typ
                  );
/*!
\rst
This function obtains the Weiss field.
Interfaces to :f:func:`f/ed_io/ed_get_g0and` for a rank-5 array.
Single-site DMFT variant.
\endrst
* @param warray: Array of frequencies
* @param dim_warray: Dimension of the array of frequencies
* @param bath: User-accessible bath array
* @param dim_bath: Dimension of the bath array
* @param G0and: the Weiss field array
* @param dim_g0and: the dimensions of the Weiss field array
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
*
*/                
void get_g0and_n5(std::complex<double> *warray,
                  int64_t *dim_warray, 
                  double *bath, 
                  int dim_bath, 
                  std::complex<double> *G0and, 
                  int64_t *dim_g0and, 
                  char *axis = (char *)"m", 
                  char *typ = (char *)"n"
                  );

/*!
\rst
This function obtains the hybridization function.
Interfaces to :f:func:`f/ed_io/ed_get_delta` for a rank-3 array.
Single-site DMFT variant.
\endrst
* @param warray: Array of frequencies
* @param dim_warray: Dimension of the array of frequencies
* @param bath: User-accessible bath array
* @param dim_bath: Dimension of the bath array
* @param Delta: the Weiss field array
* @param dim_delta: the dimensions of the Weiss field array
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
*
*/                   
void get_delta_n3(std::complex<double> *warray, 
                  int64_t *dim_warray, 
                  double *bath, 
                  int dim_bath, 
                  std::complex<double> *Delta, 
                  int64_t *dim_delta, 
                  char *axis, 
                  char *typ
                  );
                  
/*!
\rst
This function obtains the hybridization function.
Interfaces to :f:func:`f/ed_io/ed_get_delta` for a rank-5 array.
Single-site DMFT variant.
\endrst
* @param warray: Array of frequencies
* @param dim_warray: Dimension of the array of frequencies
* @param bath: User-accessible bath array
* @param dim_bath: Dimension of the bath array
* @param Delta: the Weiss field array
* @param dim_delta: the dimensions of the Weiss field array
* @param axis: integer flag for axis: `1` = `r`, otherwise `m`
* @param typ: integer flag for type: `1` = `a`, otherwise `n`
*
*/   
void get_delta_n5(std::complex<double> *warray, 
                  int64_t *dim_warray, 
                  double *bath, 
                  int dim_bath, 
                  std::complex<double> *Delta, 
                  int64_t *dim_delta, 
                  char *axis, 
                  char *typ
                  );


/*!
\rst
This function obtains the spin susceptibility.
Interfaces to :f:func:`f/ed_io/ed_spinchi` or :f:func:`f/e2i_io/ed_spinchi`.
\endrst
* @param self: susceptibility array
* @param zeta: frequency array
* @param dim_zeta: dimension of the frequency array
* @param zetaflag: flag to set (`1`) or not (`0`) calculation with internal frequency/time  array
* @param axis: flag for axis
* @param Nsites: number of inequivalent sites (`1` for single-impurity)
* @param latticeflag: flag to set (`1`) or not (`0`) calculation with real-space DMFT
*/ 
void ed_get_spinchi(std::complex<double> *self, 
                    std::complex<double> *zeta, 
                    int dim_zeta, 
                    int zetaflag, 
                    char *axis, 
                    int Nsites, 
                    int latticeflag
                    );
                    
/*!
\rst
This function obtains the charge susceptibility.
Interfaces to :f:func:`f/ed_io/ed_denschi` or :f:func:`f/e2i_io/ed_denschi`.
\endrst
* @param self: susceptibility array
* @param zeta: frequency array
* @param dim_zeta: dimension of the frequency array
* @param zetaflag: flag to set (`1`) or not (`0`) calculation with internal frequency/time  array
* @param axis: flag for axis
* @param Nsites: number of inequivalent sites (`1` for single-impurity)
* @param latticeflag: flag to set (`1`) or not (`0`) calculation with real-space DMFT
*/                     
void ed_get_denschi(std::complex<double> *self, 
                    std::complex<double> *zeta, 
                    int dim_zeta, 
                    int zetaflag, 
                    char *axis, 
                    int Nsites, 
                    int latticeflag
                    );
                    
                    
/*!
\rst
This function obtains the pairing susceptibility.
Interfaces to :f:func:`f/ed_io/ed_pairchi` or :f:func:`f/e2i_io/ed_pairchi`.
\endrst
* @param self: susceptibility array
* @param zeta: frequency array
* @param dim_zeta: dimension of the frequency array
* @param zetaflag: flag to set (`1`) or not (`0`) calculation with internal frequency/time array
* @param axis: flag for axis
* @param Nsites: number of inequivalent sites (`1` for single-impurity)
* @param latticeflag: flag to set (`1`) or not (`0`) calculation with real-space DMFT
*/                       
void ed_get_pairchi(std::complex<double> *self, 
                    std::complex<double> *zeta, 
                    int dim_zeta, 
                    int zetaflag, 
                    char *axis, 
                    int Nsites, 
                    int latticeflag
                    );
                    
/*!
\rst
This function obtains the excitonic susceptibility.
Interfaces to :f:func:`f/ed_io/ed_exctchi` or :f:func:`f/e2i_io/ed_exctchi`.
\endrst
* @param self: susceptibility array
* @param zeta: frequency array
* @param dim_zeta: dimension of the frequency array
* @param zetaflag: flag to set (`1`) or not (`0`) calculation with internal frequency/time array
* @param axis: flag for axis
* @param Nsites: number of inequivalent sites (`1` for single-impurity)
* @param latticeflag: flag to set (`1`) or not (`0`) calculation with real-space DMFT
*/                      
void ed_get_exctchi(std::complex<double> *self, 
                    std::complex<double> *zeta, 
                    int dim_zeta, 
                    int zetaflag, 
                    char *axis, 
                    int Nsites, 
                    int latticeflag
                    );

/*!
\rst
This function obtains the impurity reduced density matrix.
Interfaces to :f:func:`f/ed_io/ed_get_impurity_rdm`.
\endrst
* @param rdm: impurity reduced density matrix
* @param doprint: flag to print to file (`1`) or not (`0`) the RDM
*/                      
void ed_get_impurity_rdm(std::complex<double> *rdm, 
                    int doprint
                    );             

/*!
\rst
This function initializes the solver.
Interfaces to :f:func:`f/ed_main/ed_init_solver`.
Single-site DMFT variant.
\endrst
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
*/ 
void init_solver_site(double *bath, int64_t *dim_bath);

/*!
\rst
This function initializes the solver for an isolated impurity.
Interfaces to :f:func:`f/ed_main/ed_init_solver`.
Single-site DMFT variant.
\endrst
*/ 
void init_solver_site_nobath();

/*!
\rst
This function initializes the solver.
Interfaces to :f:func:`f/ed_main/ed_init_solver`.
Real-space DMFT variant.
\endrst
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
*/ 
void init_solver_ineq(double *bath, int64_t *dim_bath);


/*!
\rst
This function initializes the solver for an isolated impurity.
Interfaces to :f:func:`f/ed_main/ed_init_solver`.
Real-space DMFT variant.
\endrst
* @param Nlat: number of inequivalent sites
*/ 
void init_solver_ineq_nobath(int Nlat);

/*!
\rst
Solves the impurity problem.
Interfaces to :f:func:`f/ed_main/solve`.
Single-site DMFT variant.
\endrst
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param flag_gf: Green's function calculation flag (`1`=`True`, `0`=`False`)
* @param flag_mpi: mpi setting flag (`1`=`True`, `0`=`False`)
*/ 
void solve_site(double *bath, 
                int64_t *dim_bath, 
                int flag_gf, 
                int flag_mpi
                );
                
/*!
\rst
This function solves the impurity problem for the case of
an isolated impurity.
Interfaces to :f:func:`f/ed_main/solve`.
Single-site DMFT variant.
\endrst
* @param flag_gf: Green's function calculation flag (`1`=`True`, `0`=`False`)
* @param flag_mpi: mpi setting flag (`1`=`True`, `0`=`False`)
*/ 
void solve_site_nobath(
                int flag_gf, 
                int flag_mpi
                );
                
/*!
\rst
This function solves the impurity problem.
Interfaces to :f:func:`f/e2i_main/solve`.
Real-space DMFT variant.
\endrst
* @param bath: user-accessible bath array
* @param dim_bath: dimensions of the bath array
* @param flag_gf: Green's function calculation flag (`1`=`True`, `0`=`False`)
* @param mpi_lanc: parallelization setting flag (`1`=`True`, `0`=`False`)
*/                 
void solve_ineq(double *bath, 
                int64_t *dim_bath, 
                int flag_gf,                
                int mpi_lanc 
                );

/*!
\rst
This function solves the impurity problem for the case of
an isolated impuirty.
Interfaces to :f:func:`f/e2i_main/solve`.
Real-space DMFT variant.
\endrst
* @param Nlat: number of inequivalent sites
* @param flag_gf: Green's function calculation flag (`1`=`True`, `0`=`False`)
* @param mpi_lanc: parallelization setting flag (`1`=`True`, `0`=`False`)
*/                 
void solve_ineq_nobath(int Nlat, 
                       int flag_gf,                
                       int mpi_lanc 
                       );

/*!
\rst
This function finalizes the solver.
Interfaces to :f:func:`f/ed_main/finalize_solver` and :f:func:`f/e2i_main/finalize_solver`.
\endrst
* @param Nineq: number of inequivalent sites. Set to `0` for single-site DMFT.
*/  
void finalize_solver (int Nineq);

/*!
\rst
This function resets the interaction coefficients and user-provided terms
Interfaces to :f:func:`f/ed_parse_umatrix/reset_umatrix`.
\endrst
*/ 
void reset_umatrix();

/*!
\rst
This function sets the a two-body interaction term
Interfaces to :f:func:`f/ed_parse_umatrix/add_twobody_operator`.
\endrst
* @param o1: orbital index of first creation operator
* @param s1: spin index of first creation operator
* @param o2: orbital index of second creation operator
* @param s2: spin index of second creation operator
* @param o3: orbital index of first annihilation operator
* @param s3: spin index of first annihilation operator
* @param o4: orbital index of second annihilation operator
* @param s4: spin index of second annihilation operator
* @param U: interaction coefficient
*/ 
void add_twobody_operator(int o1,
                          int s1,
                          int o2,
                          int s2,
                          int o3,
                          int s3,
                          int o4,
                          int s4,
                          double U
                     );

#ifdef __cplusplus
}
#endif
#endif
