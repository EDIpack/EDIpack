MODULE EDIPACK
  !:synopsis: EDIpack library frontend
  USE ED_INPUT_VARS  , only: &
       ed_read_input , &
       ed_update_input,&
       Norb          , &
       Nspin         , &
       Nbath         , &
       Nloop         , &
       Nph           , &
       Uloc          , &
       Ust           , &
       Jh            , &
       Jx            , &
       Jp            , &
       xmu           , &
       beta          , &
       g_ph          , &
       w0_ph         , &
       eps           , &
       wini          , &
       wfin          , &
       xmin          , &
       xmax          , &
       Nsuccess      , &
       dmft_error    , &
       sb_field      , &
       cg_Scheme     , &
       nread         , &
       Lmats         , &
       Lreal         , &
       Ltau          , &
       Lpos          , &
       Hfile         , &
       HLOCfile      , &
       LOGfile       , &
       ed_mode       , &
       ed_verbose    , &
       ed_read_umatrix,&
       ed_use_kanamori,&
       ed_sparse_H    ,&
       ed_hw_bath     ,&
       ed_input_file  ,&
       lanc_nstates_total ,&
       bath_type 


  USE ED_BATH, only:                                                     &
       ed_read_dmft_bath               => read_dmft_bath               , &
       ed_get_g0and                                                    , &
       ed_get_delta                                                    , &
       ed_allocate_Hreplica            => allocate_Hreplica            , &
       ed_deallocate_Hreplica          => deallocate_Hreplica          , &
       ed_set_linit_Hreplica           => set_linit_Hreplica           , &
       ed_set_hsym_Hreplica            => set_hsym_Hreplica            , &
       ed_build_Hreplica               => build_Hreplica               , &
       ed_print_Hreplica               => print_Hreplica               , &
       ed_set_Hreplica                 => set_Hreplica                 , &
       ed_set_Hgeneral                 => set_Hgeneral                 , &
       ed_Hreplica_mask                => Hreplica_mask                , &
       ed_Hgeneral_mask                => Hgeneral_mask                , &       
       ed_get_bath_dimension           => get_bath_dimension           , &
       ed_get_g0and                    => ed_get_g0and                 , &
       ed_get_delta                    => ed_get_delta                 , &
       ed_spin_symmetrize_bath         => spin_symmetrize_bath         , &
       ed_orb_symmetrize_bath          => orb_symmetrize_bath          , &
       ed_orb_equality_bath            => orb_equality_bath            , &
       ed_ph_symmetrize_bath           => ph_symmetrize_bath           , &
       ed_ph_trans_bath                => ph_trans_bath                , &
       ed_break_symmetry_bath          => break_symmetry_bath          , &
       ed_enforce_normal_bath          => enforce_normal_bath          , &
       ed_save_array_as_bath           => save_array_as_bath       


  USE ED_AUX_FUNX, only:                        &
       ed_set_Hloc                                                     , &
       ed_set_suffix                                                   , &
       ed_reset_suffix                                                 , &
       ed_read_ImpGMatrix                 =>   read_ImpGMatrix         , &
       ed_read_ImpDMatrix                 =>   read_ImpDMatrix         , &
       ed_read_spinChimatrix              =>   read_spinChimatrix      , &
       ed_read_densChimatrix              =>   read_densChimatrix      , &
       ed_read_pairChimatrix              =>   read_pairChimatrix      , &
       ed_read_exctChimatrix              =>   read_exctChimatrix      , &
       ed_search_variable                                              , &
       ed_search_chemical_potential       => search_chemical_potential


  USE ED_IO, only: &
       ed_get_gimp            , &
       ed_get_dimp            , &
       ed_get_sigma           , &
       ed_get_g0imp           , &
       ed_get_spinChi         , &
       ed_get_densChi         , &
       ed_get_pairChi         , &
       ed_get_exctChi         , &       
       ed_get_dens            , &
       ed_get_phi             , &
       ed_get_argphi          , &
       ed_get_exct            , &
       ed_get_mag             , &
       ed_get_docc            , &
       ed_get_eimp            , &
       ed_get_epot            , &
       ed_get_eint            , &
       ed_get_ehartree        , &
       ed_get_eknot           , &
       ed_get_doubles         , &
       ed_get_evals           , &
       ed_get_imp_info        , &
       ed_get_nsectors        , &
       ed_get_neigen_sector   , &
       ed_set_neigen_sector   , &
       ed_get_impurity_rdm    , &
       ed_get_sp_dm           


  USE ED_RDM, only: &
       ed_get_reduced_rdm   => get_reduced_rdm


  USE ED_PARSE_UMATRIX, only: &
       ed_add_twobody_operator => add_twobody_operator, &
       ed_reset_umatrix => reset_umatrix

  USE ED_GREENS_FUNCTIONS, only: &
       ed_build_impG   => get_impG ,&
       ed_build_impF   => get_impF ,&
       ed_build_impD   => get_impD ,&
       ed_build_Sigma  => get_Sigma ,&
       ed_build_Self   => get_Self


  USE ED_CHI_FUNCTIONS, only: &
       ed_build_spinChi  => get_spinChi ,&
       ed_build_densChi  => get_densChi ,&
       ed_build_pairChi  => get_pairChi ,&
       ed_build_exctChi  => get_exctChi

  
  USE ED_MAIN, only:    &
       ed_init_solver , &
       ed_solve       , &
       ed_finalize_solver


  
  USE ED_BATH_FIT,  only: ed_chi2_fitgf


END MODULE EDIPACK

