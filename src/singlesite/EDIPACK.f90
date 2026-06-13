MODULE EDIPACK
  !:synopsis: EDIpack library frontend
  USE ED_INPUT_VARS , only: &
       ed_read_input      , &
       ed_update_input    , &
       bath_type          , &
       beta               , &
       cg_Scheme          , &
       dmft_error         , &
       ed_hw_bath         , &
       ed_input_file      , &
       ed_mode            , &
       ed_read_umatrix    , &
       ed_sparse_H        , &
       ed_use_kanamori    , &
       ed_verbose         , &
       eps                , &
       g_ph               , &
       Hfile              , &
       HLOCfile           , &
       Jh                 , &
       Jp                 , &
       Jx                 , &
       lanc_nstates_total , &
       Lmats              , &
       LOGfile            , &
       Lpos               , &
       Lreal              , &
       Ltau               , &
       Ltimes             , &
       Nbath              , &
       Nloop              , &
       Norb               , &
       Nph                , &
       nread              , &
       Nspin              , &
       Nsuccess           , &
       sb_field           , &
       Uloc               , &
       Ust                , &
       Tmax               , &
       w0_ph              , &
       wfin               , &
       wini               , &
       xmax               , &
       xmin               , &
       xmu



  USE ED_BATH, only:                                                     &
       ed_allocate_Hreplica            => allocate_Hreplica            , &
       ed_break_symmetry_bath          => break_symmetry_bath          , &
       ed_build_Hreplica               => build_Hreplica               , &
       ed_deallocate_Hreplica          => deallocate_Hreplica          , &
       ed_enforce_normal_bath          => enforce_normal_bath          , &
       ed_get_bath_dimension           => get_bath_dimension           , &
       ed_get_delta                                                    , &
       ed_get_delta                    => ed_get_delta                 , &
       ed_get_g0and                                                    , &
       ed_get_g0and                    => ed_get_g0and                 , &
       ed_Hgeneral_mask                => Hgeneral_mask                , &       
       ed_Hreplica_mask                => Hreplica_mask                , &
       ed_orb_equality_bath            => orb_equality_bath            , &
       ed_orb_symmetrize_bath          => orb_symmetrize_bath          , &
       ed_ph_symmetrize_bath           => ph_symmetrize_bath           , &
       ed_ph_trans_bath                => ph_trans_bath                , &
       ed_print_Hreplica               => print_Hreplica               , &
       ed_read_dmft_bath               => read_dmft_bath               , &
       ed_save_array_as_bath           => save_array_as_bath           , &        
       ed_set_Hgeneral                 => set_Hgeneral                 , &
       ed_set_Hreplica                 => set_Hreplica                 , &
       ed_set_hsym_Hreplica            => set_hsym_Hreplica            , &
       ed_set_linit_Hreplica           => set_linit_Hreplica           , &
       ed_spin_symmetrize_bath         => spin_symmetrize_bath


  USE ED_AUX_FUNX, only: &
       ed_read_densChimatrix           =>   read_densChimatrix      , &
       ed_read_exctChimatrix           =>   read_exctChimatrix      , &
       ed_read_ImpDMatrix              =>   read_ImpDMatrix         , &
       ed_read_ImpGMatrix              =>   read_ImpGMatrix         , &
       ed_read_pairChimatrix           =>   read_pairChimatrix      , &
       ed_read_spinChimatrix           =>   read_spinChimatrix      , &
       ed_reset_suffix                                              , &
       ed_search_chemical_potential    => search_chemical_potential , &
       ed_search_variable                                           , &
       ed_set_A_ph                                                  , &
       ed_set_G_ph                                                  , &
       ed_set_Hloc                                                  , &
       ed_set_suffix


  USE ED_IO, only: &
       ed_get_argphi          , &
       ed_get_dens            , &
       ed_get_densChi         , &
       ed_get_denmat          , &
       ed_get_dimp            , &
       ed_get_docc            , &
       ed_get_doubles         , &
       ed_get_ehartree        , &
       ed_get_eimp            , &
       ed_get_eint            , &
       ed_get_eknot           , &
       ed_get_ephon           , &
       ed_get_epot            , &
       ed_get_evals           , &
       ed_get_exct            , &
       ed_get_exctChi         , &       
       ed_get_g0imp           , &
       ed_get_gimp            , &
       ed_get_imp_info        , &
       ed_get_impurity_rdm    , &
       ed_get_Kcdg            , &
       ed_get_Kn              , &
       ed_get_mag             , &
       ed_get_neigen_sector   , &
       ed_get_Normcdg         , &
       ed_get_Normn           , &
       ed_get_nsectors        , &
       ed_get_Pcdg            , &
       ed_get_pairChi         , &
       ed_get_Pn              , &
       ed_get_phi             , &
       ed_get_phon            , &
       ed_get_Scdg            , &
       ed_get_sigma           , &
       ed_get_Sn              , &
       ed_get_sp_dm           , &
       ed_get_spinChi         , &
       ed_set_neigen_sector


  USE ED_RDM, only: &
       ed_get_reduced_rdm   => get_reduced_rdm


  USE ED_PARSE_UMATRIX, only: &
       ed_add_twobody_operator => add_twobody_operator, &
       ed_reset_umatrix => reset_umatrix

  USE ED_GREENS_FUNCTIONS, only: &
       ed_build_impD   => get_impD  , &  
       ed_build_impF   => get_impF  , &
       ed_build_impG   => get_impG  , &
       ed_build_Self   => get_Self  , &
       ed_build_Sigma  => get_Sigma


  USE ED_CHI_FUNCTIONS, only: &
       ed_build_densChi  => get_densChi , &
       ed_build_exctChi  => get_exctChi , &
       ed_build_pairChi  => get_pairChi , &
       ed_build_spinChi  => get_spinChi

  USE ED_KRYLOV, only: &
       ed_build_krylov_state_complexity => krylov_state_complexity

  
  USE ED_MAIN, only:    &
       ed_finalize_solver , &
       ed_init_solver     , &
       ed_solve


  
  USE ED_BATH_FIT,  only: &
       ed_chi2_fitgf


END MODULE EDIPACK
