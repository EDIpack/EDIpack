MODULE EDIPACK2INEQ
  !:synopsis: EDIpack real-space DMFT extension
  USE EDIPACK
  USE E2I_AUX_FUNX, only:                        &
       ed_set_Hloc


  USE E2I_BATH, only:                                                    &
       ed_set_Hreplica                 => set_Hreplica                 , &
       ed_set_Hgeneral                 => set_Hgeneral                 , &
       ed_spin_symmetrize_bath         => spin_symmetrize_bath         , &
       ed_orb_symmetrize_bath          => orb_symmetrize_bath          , &
       ed_orb_equality_bath            => orb_equality_bath            , &
       ed_ph_symmetrize_bath           => ph_symmetrize_bath           , &
       ed_ph_trans_bath                => ph_trans_bath                , &
       ed_break_symmetry_bath          => break_symmetry_bath          , &
       ed_enforce_normal_bath          => enforce_normal_bath          , &
       ed_save_array_as_bath           => save_array_as_bath       




  USE E2I_IO, only: &
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
       ed_get_arg             , &
       ed_get_mag             , &
       ed_get_docc            , &
       ed_get_eimp            , &
       ed_get_epot            , &
       ed_get_eint            , &
       ed_get_ehartree        , &
       ed_get_eknot           , &
       ed_get_doubles!          , &
  ! ed_get_sp_dm  


  USE E2I_MAIN, only:    &
       ed_init_solver , &
       ed_solve       , &
       ed_finalize_solver

  USE E2I_BATH_FIT,  only: ed_chi2_fitgf


END MODULE EDIPACK2INEQ

