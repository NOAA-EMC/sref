





MODULE module_state_description

  INTEGER, PARAMETER :: passiveqv = 0
  INTEGER, PARAMETER :: kesslerscheme = 1
  INTEGER, PARAMETER :: linscheme = 2
  INTEGER, PARAMETER :: wsm3scheme = 3
  INTEGER, PARAMETER :: wsm5scheme = 4
  INTEGER, PARAMETER :: etampnew = 5
  INTEGER, PARAMETER :: wsm6scheme = 6
  INTEGER, PARAMETER :: gsfcgcescheme = 7
  INTEGER, PARAMETER :: thompson = 8
  INTEGER, PARAMETER :: milbrandt2mom = 9
  INTEGER, PARAMETER :: morr_two_moment = 10
  INTEGER, PARAMETER :: cammgmpscheme = 11
  INTEGER, PARAMETER :: sbu_ylinscheme = 13
  INTEGER, PARAMETER :: wdm5scheme = 14
  INTEGER, PARAMETER :: wdm6scheme = 16
  INTEGER, PARAMETER :: nssl_2mom = 17
  INTEGER, PARAMETER :: nssl_2momccn = 18
  INTEGER, PARAMETER :: nssl_1mom = 19
  INTEGER, PARAMETER :: nssl_1momlfo = 21
  INTEGER, PARAMETER :: etampold = 95
  INTEGER, PARAMETER :: nodfimoist = -1
  INTEGER, PARAMETER :: passiveqv_dfi = 0
  INTEGER, PARAMETER :: kesslerscheme_dfi = 1
  INTEGER, PARAMETER :: linscheme_dfi = 2
  INTEGER, PARAMETER :: wsm3scheme_dfi = 3
  INTEGER, PARAMETER :: wsm5scheme_dfi = 4
  INTEGER, PARAMETER :: etampnew_dfi = 5
  INTEGER, PARAMETER :: wsm6scheme_dfi = 6
  INTEGER, PARAMETER :: gsfcgcescheme_dfi = 7
  INTEGER, PARAMETER :: thompson_dfi = 8
  INTEGER, PARAMETER :: milbrandt2mom_dfi = 9
  INTEGER, PARAMETER :: morr_two_moment_dfi = 10
  INTEGER, PARAMETER :: wdm5scheme_dfi = 14
  INTEGER, PARAMETER :: wdm6scheme_dfi = 16
  INTEGER, PARAMETER :: nssl_2mom_dfi = 17
  INTEGER, PARAMETER :: nssl_2mom_dficcn = 18
  INTEGER, PARAMETER :: nssl_1mom_dfi = 19
  INTEGER, PARAMETER :: nssl_1momlfo_dfi = 21
  INTEGER, PARAMETER :: etampold_dfi = 95
  INTEGER, PARAMETER :: noprogn = 0
  INTEGER, PARAMETER :: progndrop = 1
  INTEGER, PARAMETER :: rrtmscheme = 1
  INTEGER, PARAMETER :: camlwscheme = 3
  INTEGER, PARAMETER :: rrtmg_lwscheme = 4
  INTEGER, PARAMETER :: goddardlwscheme = 5
  INTEGER, PARAMETER :: flglwscheme = 7
  INTEGER, PARAMETER :: gfdllwscheme = 99
  INTEGER, PARAMETER :: heldsuarez = 31
  INTEGER, PARAMETER :: swradscheme = 1
  INTEGER, PARAMETER :: gsfcswscheme = 2
  INTEGER, PARAMETER :: camswscheme = 3
  INTEGER, PARAMETER :: rrtmg_swscheme = 4
  INTEGER, PARAMETER :: goddardswscheme = 5
  INTEGER, PARAMETER :: flgswscheme = 7
  INTEGER, PARAMETER :: gfdlswscheme = 99
  INTEGER, PARAMETER :: sfclayscheme = 1
  INTEGER, PARAMETER :: myjsfcscheme = 2
  INTEGER, PARAMETER :: gfssfcscheme = 3
  INTEGER, PARAMETER :: qnsesfcscheme = 4
  INTEGER, PARAMETER :: mynnsfcscheme = 5
  INTEGER, PARAMETER :: pxsfcscheme = 7
  INTEGER, PARAMETER :: temfsfcscheme = 10
  INTEGER, PARAMETER :: sfclayrevscheme = 11
  INTEGER, PARAMETER :: idealscmsfcscheme = 89
  INTEGER, PARAMETER :: noahucmscheme = 1
  INTEGER, PARAMETER :: bepscheme = 2
  INTEGER, PARAMETER :: bep_bemscheme = 3
  INTEGER, PARAMETER :: slabscheme = 1
  INTEGER, PARAMETER :: lsmscheme = 2
  INTEGER, PARAMETER :: ruclsmscheme = 3
  INTEGER, PARAMETER :: noahmpscheme = 4
  INTEGER, PARAMETER :: clmscheme = 5
  INTEGER, PARAMETER :: pxlsmscheme = 7
  INTEGER, PARAMETER :: ssibscheme = 8
  INTEGER, PARAMETER :: ysuscheme = 1
  INTEGER, PARAMETER :: myjpblscheme = 2
  INTEGER, PARAMETER :: gfsscheme = 3
  INTEGER, PARAMETER :: qnsepblscheme = 4
  INTEGER, PARAMETER :: qnsepbl09scheme = 94
  INTEGER, PARAMETER :: mynnpblscheme2 = 5
  INTEGER, PARAMETER :: mynnpblscheme3 = 6
  INTEGER, PARAMETER :: mynn_tkebudget = 1
  INTEGER, PARAMETER :: acmpblscheme = 7
  INTEGER, PARAMETER :: boulacscheme = 8
  INTEGER, PARAMETER :: camuwpblscheme = 9
  INTEGER, PARAMETER :: mrfscheme = 99
  INTEGER, PARAMETER :: temfpblscheme = 10
  INTEGER, PARAMETER :: gbmpblscheme = 12
  INTEGER, PARAMETER :: kfetascheme = 1
  INTEGER, PARAMETER :: bmjscheme = 2
  INTEGER, PARAMETER :: gdscheme = 93
  INTEGER, PARAMETER :: sasscheme = 84
  INTEGER, PARAMETER :: meso_sas = 85
  INTEGER, PARAMETER :: osasscheme = 4
  INTEGER, PARAMETER :: g3scheme = 5
  INTEGER, PARAMETER :: gfscheme = 3
  INTEGER, PARAMETER :: camzmscheme = 7
  INTEGER, PARAMETER :: g3tave = 1
  INTEGER, PARAMETER :: tiedtkescheme = 6
  INTEGER, PARAMETER :: nsasscheme = 14
  INTEGER, PARAMETER :: kfscheme = 99
  INTEGER, PARAMETER :: g3shcuscheme = 1
  INTEGER, PARAMETER :: camuwshcuscheme = 2
  INTEGER, PARAMETER :: grimsshcuscheme = 3
  INTEGER, PARAMETER :: fogsettling0 = 0
  INTEGER, PARAMETER :: fogsettling1 = 1
  INTEGER, PARAMETER :: fogsettling2 = 2
  INTEGER, PARAMETER :: psufddagd = 1
  INTEGER, PARAMETER :: psusfddagd = 1
  INTEGER, PARAMETER :: spnudging = 2
  INTEGER, PARAMETER :: slopeopt = 1
  INTEGER, PARAMETER :: gwdopt = 1
  INTEGER, PARAMETER :: omlscheme = 1
  INTEGER, PARAMETER :: pwp3dscheme = 2
  INTEGER, PARAMETER :: scmopt = 1
  INTEGER, PARAMETER :: prec_acc = 1
  INTEGER, PARAMETER :: bucketropt = 1
  INTEGER, PARAMETER :: restofwrf = 0
  INTEGER, PARAMETER :: original_mom = 1
  INTEGER, PARAMETER :: weno_mom = 3
  INTEGER, PARAMETER :: original = 0
  INTEGER, PARAMETER :: positivedef = 1
  INTEGER, PARAMETER :: monotonic = 2
  INTEGER, PARAMETER :: weno_scalar = 3
  INTEGER, PARAMETER :: wenopd_scalar = 4
  INTEGER, PARAMETER :: maxmin_output = 1
  INTEGER, PARAMETER :: nwp_output = 1
  INTEGER, PARAMETER :: dfi_setup = 0
  INTEGER, PARAMETER :: dfi_bck = 1
  INTEGER, PARAMETER :: dfi_fwd = 2
  INTEGER, PARAMETER :: dfi_fst = 3
  INTEGER, PARAMETER :: dfi_startfwd = 4
  INTEGER, PARAMETER :: dfi_startbck = 5
  INTEGER, PARAMETER :: dfi_nodfi = 0
  INTEGER, PARAMETER :: dfi_dfl = 1
  INTEGER, PARAMETER :: dfi_ddfi = 2
  INTEGER, PARAMETER :: dfi_tdfi = 3
  INTEGER, PARAMETER :: realonly = 1
  INTEGER, PARAMETER :: tconly = 2
  INTEGER, PARAMETER :: reg_interp = 0
  INTEGER, PARAMETER :: flat_p_interp = 1
  INTEGER, PARAMETER :: notenddiag = 0
  INTEGER, PARAMETER :: usetenddiag = 1
  INTEGER, PARAMETER :: no_trajectory = 0
  INTEGER, PARAMETER :: um_trajectory = 1
  INTEGER, PARAMETER :: albsi_zero = 0
  INTEGER, PARAMETER :: albsi_one = 1
  INTEGER, PARAMETER :: albsi_two = 2
  INTEGER, PARAMETER :: snowsi_zero = 0
  INTEGER, PARAMETER :: snowsi_one = 1
  INTEGER, PARAMETER :: icedepth_zero = 0
  INTEGER, PARAMETER :: icedepth_one = 1
  INTEGER, PARAMETER :: notseries = 0
  INTEGER, PARAMETER :: tseries = 1
  INTEGER, PARAMETER :: ltng_none = 0
  INTEGER, PARAMETER :: ltng_crm_pr92w = 1
  INTEGER, PARAMETER :: ltng_crm_pr92z = 2
  INTEGER, PARAMETER :: ltng_cpm_pr92z = 11
  INTEGER, PARAMETER :: io_intio = 1
  INTEGER, PARAMETER :: io_netcdf = 2
  INTEGER, PARAMETER :: io_hdf = 3
  INTEGER, PARAMETER :: io_phdf5 = 4
  INTEGER, PARAMETER :: io_grib1 = 5
  INTEGER, PARAMETER :: io_mcel = 6
  INTEGER, PARAMETER :: io_esmf = 7
  INTEGER, PARAMETER :: io_yyy = 8
  INTEGER, PARAMETER :: io_zzz = 9
  INTEGER, PARAMETER :: io_grib2 = 10
  INTEGER, PARAMETER :: io_pnetcdf = 11
  INTEGER, PARAMETER :: no_wrfhydro = 0
  INTEGER, PARAMETER :: wrfhydro = 1
  INTEGER, PARAMETER :: fire_sfire = 2
  INTEGER, PARAMETER :: noavgflxem = 0
  INTEGER, PARAMETER :: avgflxem = 1
  INTEGER, PARAMETER :: noavgflxcugd = 0
  INTEGER, PARAMETER :: avgflxcugd = 1
  INTEGER, PARAMETER :: no_stoch_force = 0
  INTEGER, PARAMETER :: stoch_backscatter = 1
  INTEGER, PARAMETER :: nosfs = 0
  INTEGER, PARAMETER :: nba1 = 1
  INTEGER, PARAMETER :: nba2 = 2
  INTEGER, PARAMETER :: mout = 1
  INTEGER, PARAMETER :: skip_press_diags = 0
  INTEGER, PARAMETER :: press_diags = 1
  INTEGER, PARAMETER :: no_perturb_bdy = 0
  INTEGER, PARAMETER :: perturb_bdy_stoch_patrn = 1
  INTEGER, PARAMETER :: perturb_bdy_user_patrn = 2
  INTEGER, PARAMETER :: tracer_test1 = 2

  INTEGER, PARAMETER :: PARAM_qv = 1
  INTEGER :: P_qv = 1
  LOGICAL :: F_qv = .FALSE.
  INTEGER, PARAMETER :: PARAM_qc = 2
  INTEGER :: P_qc = 1
  LOGICAL :: F_qc = .FALSE.
  INTEGER, PARAMETER :: PARAM_qr = 3
  INTEGER :: P_qr = 1
  LOGICAL :: F_qr = .FALSE.
  INTEGER, PARAMETER :: PARAM_qi = 4
  INTEGER :: P_qi = 1
  LOGICAL :: F_qi = .FALSE.
  INTEGER, PARAMETER :: PARAM_qs = 5
  INTEGER :: P_qs = 1
  LOGICAL :: F_qs = .FALSE.
  INTEGER, PARAMETER :: PARAM_qg = 6
  INTEGER :: P_qg = 1
  LOGICAL :: F_qg = .FALSE.
  INTEGER, PARAMETER :: PARAM_qh = 7
  INTEGER :: P_qh = 1
  LOGICAL :: F_qh = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_moist = 8
  INTEGER :: NUM_moist = 1
  INTEGER, PARAMETER :: PARAM_dfi_qv = 1
  INTEGER :: P_dfi_qv = 1
  LOGICAL :: F_dfi_qv = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qc = 2
  INTEGER :: P_dfi_qc = 1
  LOGICAL :: F_dfi_qc = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qr = 3
  INTEGER :: P_dfi_qr = 1
  LOGICAL :: F_dfi_qr = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qi = 4
  INTEGER :: P_dfi_qi = 1
  LOGICAL :: F_dfi_qi = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qs = 5
  INTEGER :: P_dfi_qs = 1
  LOGICAL :: F_dfi_qs = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qg = 6
  INTEGER :: P_dfi_qg = 1
  LOGICAL :: F_dfi_qg = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qh = 7
  INTEGER :: P_dfi_qh = 1
  LOGICAL :: F_dfi_qh = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_dfi_moist = 8
  INTEGER :: NUM_dfi_moist = 1
  INTEGER, PARAMETER :: PARAM_qndrop = 1
  INTEGER :: P_qndrop = 1
  LOGICAL :: F_qndrop = .FALSE.
  INTEGER, PARAMETER :: PARAM_qni = 2
  INTEGER :: P_qni = 1
  LOGICAL :: F_qni = .FALSE.
  INTEGER, PARAMETER :: PARAM_qt = 3
  INTEGER :: P_qt = 1
  LOGICAL :: F_qt = .FALSE.
  INTEGER, PARAMETER :: PARAM_qns = 4
  INTEGER :: P_qns = 1
  LOGICAL :: F_qns = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnr = 5
  INTEGER :: P_qnr = 1
  LOGICAL :: F_qnr = .FALSE.
  INTEGER, PARAMETER :: PARAM_qng = 6
  INTEGER :: P_qng = 1
  LOGICAL :: F_qng = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnh = 7
  INTEGER :: P_qnh = 1
  LOGICAL :: F_qnh = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnn = 8
  INTEGER :: P_qnn = 1
  LOGICAL :: F_qnn = .FALSE.
  INTEGER, PARAMETER :: PARAM_qnc = 9
  INTEGER :: P_qnc = 1
  LOGICAL :: F_qnc = .FALSE.
  INTEGER, PARAMETER :: PARAM_qvolg = 10
  INTEGER :: P_qvolg = 1
  LOGICAL :: F_qvolg = .FALSE.
  INTEGER, PARAMETER :: PARAM_qke_adv = 11
  INTEGER :: P_qke_adv = 1
  LOGICAL :: F_qke_adv = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_scalar = 12
  INTEGER :: NUM_scalar = 1
  INTEGER, PARAMETER :: PARAM_dfi_qndrop = 1
  INTEGER :: P_dfi_qndrop = 1
  LOGICAL :: F_dfi_qndrop = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qni = 2
  INTEGER :: P_dfi_qni = 1
  LOGICAL :: F_dfi_qni = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qt = 3
  INTEGER :: P_dfi_qt = 1
  LOGICAL :: F_dfi_qt = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qns = 4
  INTEGER :: P_dfi_qns = 1
  LOGICAL :: F_dfi_qns = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qnr = 5
  INTEGER :: P_dfi_qnr = 1
  LOGICAL :: F_dfi_qnr = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qng = 6
  INTEGER :: P_dfi_qng = 1
  LOGICAL :: F_dfi_qng = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qnh = 7
  INTEGER :: P_dfi_qnh = 1
  LOGICAL :: F_dfi_qnh = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qnn = 8
  INTEGER :: P_dfi_qnn = 1
  LOGICAL :: F_dfi_qnn = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qnc = 9
  INTEGER :: P_dfi_qnc = 1
  LOGICAL :: F_dfi_qnc = .FALSE.
  INTEGER, PARAMETER :: PARAM_dfi_qvolg = 10
  INTEGER :: P_dfi_qvolg = 1
  LOGICAL :: F_dfi_qvolg = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_dfi_scalar = 11
  INTEGER :: NUM_dfi_scalar = 1
  INTEGER, PARAMETER :: PARAM_ocarbon = 1
  INTEGER :: P_ocarbon = 1
  LOGICAL :: F_ocarbon = .FALSE.
  INTEGER, PARAMETER :: PARAM_seasalt = 2
  INTEGER :: P_seasalt = 1
  LOGICAL :: F_seasalt = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust = 3
  INTEGER :: P_dust = 1
  LOGICAL :: F_dust = .FALSE.
  INTEGER, PARAMETER :: PARAM_bcarbon = 4
  INTEGER :: P_bcarbon = 1
  LOGICAL :: F_bcarbon = .FALSE.
  INTEGER, PARAMETER :: PARAM_sulfate = 5
  INTEGER :: P_sulfate = 1
  LOGICAL :: F_sulfate = .FALSE.
  INTEGER, PARAMETER :: PARAM_upperaer = 6
  INTEGER :: P_upperaer = 1
  LOGICAL :: F_upperaer = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_aerod = 7
  INTEGER :: NUM_aerod = 1
  INTEGER, PARAMETER :: PARAM_mth01 = 1
  INTEGER :: P_mth01 = 1
  LOGICAL :: F_mth01 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth02 = 2
  INTEGER :: P_mth02 = 1
  LOGICAL :: F_mth02 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth03 = 3
  INTEGER :: P_mth03 = 1
  LOGICAL :: F_mth03 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth04 = 4
  INTEGER :: P_mth04 = 1
  LOGICAL :: F_mth04 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth05 = 5
  INTEGER :: P_mth05 = 1
  LOGICAL :: F_mth05 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth06 = 6
  INTEGER :: P_mth06 = 1
  LOGICAL :: F_mth06 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth07 = 7
  INTEGER :: P_mth07 = 1
  LOGICAL :: F_mth07 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth08 = 8
  INTEGER :: P_mth08 = 1
  LOGICAL :: F_mth08 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth09 = 9
  INTEGER :: P_mth09 = 1
  LOGICAL :: F_mth09 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth10 = 10
  INTEGER :: P_mth10 = 1
  LOGICAL :: F_mth10 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth11 = 11
  INTEGER :: P_mth11 = 1
  LOGICAL :: F_mth11 = .FALSE.
  INTEGER, PARAMETER :: PARAM_mth12 = 12
  INTEGER :: P_mth12 = 1
  LOGICAL :: F_mth12 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_ozmixm = 13
  INTEGER :: NUM_ozmixm = 1
  INTEGER, PARAMETER :: PARAM_sul = 1
  INTEGER :: P_sul = 1
  LOGICAL :: F_sul = .FALSE.
  INTEGER, PARAMETER :: PARAM_sslt = 2
  INTEGER :: P_sslt = 1
  LOGICAL :: F_sslt = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust1 = 3
  INTEGER :: P_dust1 = 1
  LOGICAL :: F_dust1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust2 = 4
  INTEGER :: P_dust2 = 1
  LOGICAL :: F_dust2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust3 = 5
  INTEGER :: P_dust3 = 1
  LOGICAL :: F_dust3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_dust4 = 6
  INTEGER :: P_dust4 = 1
  LOGICAL :: F_dust4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_ocpho = 7
  INTEGER :: P_ocpho = 1
  LOGICAL :: F_ocpho = .FALSE.
  INTEGER, PARAMETER :: PARAM_bcpho = 8
  INTEGER :: P_bcpho = 1
  LOGICAL :: F_bcpho = .FALSE.
  INTEGER, PARAMETER :: PARAM_ocphi = 9
  INTEGER :: P_ocphi = 1
  LOGICAL :: F_ocphi = .FALSE.
  INTEGER, PARAMETER :: PARAM_bcphi = 10
  INTEGER :: P_bcphi = 1
  LOGICAL :: F_bcphi = .FALSE.
  INTEGER, PARAMETER :: PARAM_bg = 11
  INTEGER :: P_bg = 1
  LOGICAL :: F_bg = .FALSE.
  INTEGER, PARAMETER :: PARAM_volc = 12
  INTEGER :: P_volc = 1
  LOGICAL :: F_volc = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_aerosolc = 13
  INTEGER :: NUM_aerosolc = 1
  INTEGER, PARAMETER :: PARAM_u_ndg_new = 1
  INTEGER :: P_u_ndg_new = 1
  LOGICAL :: F_u_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_v_ndg_new = 2
  INTEGER :: P_v_ndg_new = 1
  LOGICAL :: F_v_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_t_ndg_new = 3
  INTEGER :: P_t_ndg_new = 1
  LOGICAL :: F_t_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_q_ndg_new = 4
  INTEGER :: P_q_ndg_new = 1
  LOGICAL :: F_q_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_ph_ndg_new = 5
  INTEGER :: P_ph_ndg_new = 1
  LOGICAL :: F_ph_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_u_ndg_old = 6
  INTEGER :: P_u_ndg_old = 1
  LOGICAL :: F_u_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_v_ndg_old = 7
  INTEGER :: P_v_ndg_old = 1
  LOGICAL :: F_v_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_t_ndg_old = 8
  INTEGER :: P_t_ndg_old = 1
  LOGICAL :: F_t_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_q_ndg_old = 9
  INTEGER :: P_q_ndg_old = 1
  LOGICAL :: F_q_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_ph_ndg_old = 10
  INTEGER :: P_ph_ndg_old = 1
  LOGICAL :: F_ph_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_fdda3d = 11
  INTEGER :: NUM_fdda3d = 1
  INTEGER, PARAMETER :: PARAM_mu_ndg_new = 1
  INTEGER :: P_mu_ndg_new = 1
  LOGICAL :: F_mu_ndg_new = .FALSE.
  INTEGER, PARAMETER :: PARAM_mu_ndg_old = 2
  INTEGER :: P_mu_ndg_old = 1
  LOGICAL :: F_mu_ndg_old = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_fdda2d = 3
  INTEGER :: NUM_fdda2d = 1
  INTEGER, PARAMETER :: PARAM_advh_qv = 1
  INTEGER :: P_advh_qv = 1
  LOGICAL :: F_advh_qv = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_advh_t = 2
  INTEGER :: NUM_advh_t = 1
  INTEGER, PARAMETER :: PARAM_advz_qv = 1
  INTEGER :: P_advz_qv = 1
  LOGICAL :: F_advz_qv = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_advz_t = 2
  INTEGER :: NUM_advz_t = 1
  INTEGER, PARAMETER :: PARAM_m11 = 1
  INTEGER :: P_m11 = 1
  LOGICAL :: F_m11 = .FALSE.
  INTEGER, PARAMETER :: PARAM_m22 = 2
  INTEGER :: P_m22 = 1
  LOGICAL :: F_m22 = .FALSE.
  INTEGER, PARAMETER :: PARAM_m33 = 3
  INTEGER :: P_m33 = 1
  LOGICAL :: F_m33 = .FALSE.
  INTEGER, PARAMETER :: PARAM_m12 = 4
  INTEGER :: P_m12 = 1
  LOGICAL :: F_m12 = .FALSE.
  INTEGER, PARAMETER :: PARAM_m13 = 5
  INTEGER :: P_m13 = 1
  LOGICAL :: F_m13 = .FALSE.
  INTEGER, PARAMETER :: PARAM_m23 = 6
  INTEGER :: P_m23 = 1
  LOGICAL :: F_m23 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_nba_mij = 7
  INTEGER :: NUM_nba_mij = 1
  INTEGER, PARAMETER :: PARAM_r12 = 1
  INTEGER :: P_r12 = 1
  LOGICAL :: F_r12 = .FALSE.
  INTEGER, PARAMETER :: PARAM_r13 = 2
  INTEGER :: P_r13 = 1
  LOGICAL :: F_r13 = .FALSE.
  INTEGER, PARAMETER :: PARAM_r23 = 3
  INTEGER :: P_r23 = 1
  LOGICAL :: F_r23 = .FALSE.
  INTEGER, PARAMETER :: PARAM_smnsmn = 4
  INTEGER :: P_smnsmn = 1
  LOGICAL :: F_smnsmn = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_nba_rij = 5
  INTEGER :: NUM_nba_rij = 1
  INTEGER, PARAMETER :: PARAM_NUM_chem = 1
  INTEGER :: NUM_chem = 1
  INTEGER, PARAMETER :: PARAM_tr17_1 = 1
  INTEGER :: P_tr17_1 = 1
  LOGICAL :: F_tr17_1 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_2 = 2
  INTEGER :: P_tr17_2 = 1
  LOGICAL :: F_tr17_2 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_3 = 3
  INTEGER :: P_tr17_3 = 1
  LOGICAL :: F_tr17_3 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_4 = 4
  INTEGER :: P_tr17_4 = 1
  LOGICAL :: F_tr17_4 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_5 = 5
  INTEGER :: P_tr17_5 = 1
  LOGICAL :: F_tr17_5 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_6 = 6
  INTEGER :: P_tr17_6 = 1
  LOGICAL :: F_tr17_6 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_7 = 7
  INTEGER :: P_tr17_7 = 1
  LOGICAL :: F_tr17_7 = .FALSE.
  INTEGER, PARAMETER :: PARAM_tr17_8 = 8
  INTEGER :: P_tr17_8 = 1
  LOGICAL :: F_tr17_8 = .FALSE.
  INTEGER, PARAMETER :: PARAM_NUM_tracer = 9
  INTEGER :: NUM_tracer = 1
  INTEGER, PARAMETER :: P_XSB = 1
  INTEGER, PARAMETER :: P_XEB = 2
  INTEGER, PARAMETER :: P_YSB = 3
  INTEGER, PARAMETER :: P_YEB = 4
  INTEGER, PARAMETER :: NUM_TIME_LEVELS = 2
  INTEGER , PARAMETER :: PARAM_FIRST_SCALAR = 2
CONTAINS
SUBROUTINE init_module_state_description
END SUBROUTINE init_module_state_description
END MODULE module_state_description
