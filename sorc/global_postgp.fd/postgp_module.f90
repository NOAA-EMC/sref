!-------------------------------------------------------------------------------
module postgp_module
!$$$  Module documentation block
!
! Module:    postgp_module   Shared data for postgp
!   Prgmmr: Iredell          Org: W/PX23     Date: 1999-01-18
!
! Abstract: This module provides shared data for the postgp program.
!   
! Program history log:
!   1999-01-18  Mark Iredell
!   2005-02-02  Mark Iredell  add new flux fields
!
! Public Variables:
!   nflxs           Integer parameter number of input flux scalar fields
!   iflxs_*         Integer parameter index for input flux scalar field *
!   jpdsflxs        Integer (1:24,nflxs) parameter PDS type and level ids
!                   for input flux scalar fields
!   nflxv           Integer parameter number of input flux vector fields
!   iflxu_*         Integer parameter index for input flux vector X field *
!   iflxv_*         Integer parameter index for input flux vector Y field *
!   jpdsflxu        Integer (1:24,nflxv) parameter PDS type and level ids
!                   for input flux vector X fields
!   jpdsflxv        Integer (1:24,nflxv) parameter PDS type and level ids
!                   for input flux vector Y fields
!   nfoxs           Integer parameter number of output flux scalar fields
!   ifoxs_*         Integer parameter index for output flux scalar field *
!   jpdsfoxs        Integer (1:24,nfoxs) parameter PDS type and level ids
!                   for output flux scalar fields
!   nfoxv           Integer parameter number of output flux vector fields
!   ifoxu_*         Integer parameter index for output flux vector X field *
!   ifoxv_*         Integer parameter index for output flux vector Y field *
!   jpdsfoxu        Integer (1:24,nfoxv) parameter PDS type and level ids
!                   for output flux vector X fields
!   jpdsfoxv        Integer (1:24,nfoxv) parameter PDS type and level ids
!                   for output flux vector Y fields
!   nsuns           Integer parameter number of output sundry scalar fields
!   isuns_*         Integer parameter index for output sundry scalar field *
!   jpdssuns        Integer (1:24,nsuns) parameter PDS type and level ids
!                   for output sundry scalar fields
!   nsunv           Integer parameter number of output sundry vector fields
!   isunu_*         Integer parameter index for output sundry vector X field *
!   isunv_*         Integer parameter index for output sundry vector Y field *
!   jpdssunu        Integer (1:24,nsunv) parameter PDS type and level ids
!                   for output sundry vector X fields
!   jpdssunv        Integer (1:24,nsunv) parameter PDS type and level ids
!                   for output sundry vector Y fields
!   npos            Integer parameter number of output p-level scalar fields
!   ipos_*          Integer parameter index for output p-level scalar field *
!   jpdspos         Integer (1:24,npos) parameter PDS type and level ids
!                   for output p-level scalar fields
!   npov            Integer parameter number of output p-level vector fields
!   ipou_*          Integer parameter index for output p-level vector X field *
!   ipov_*          Integer parameter index for output p-level vector Y field *
!   jpdspou         Integer (1:24,npov) parameter PDS type and level ids
!                   for output p-level vector X fields
!   jpdspov         Integer (1:24,npov) parameter PDS type and level ids
!                   for output p-level vector Y fields
!   npts            Integer parameter number of output p-layer scalar fields
!   ipts_*          Integer parameter index for output p-layer scalar field *
!   jpdspts         Integer (1:24,npts) parameter PDS type and level ids
!                   for output p-layer scalar fields
!   nptv            Integer parameter number of output p-layer vector fields
!   iptu_*          Integer parameter index for output p-layer vector X field *
!   iptv_*          Integer parameter index for output p-layer vector Y field *
!   jpdsptu         Integer (1:24,nptv) parameter PDS type and level ids
!                   for output p-layer vector X fields
!   jpdsptv         Integer (1:24,nptv) parameter PDS type and level ids
!                   for output p-layer vector Y fields
!   nzzs            Integer parameter number of output z-level scalar fields
!   izzs_*          Integer parameter index for output z-level scalar field *
!   jpdszzs         Integer (1:24,nzzs) parameter PDS type and level ids
!                   for output z-level scalar fields
!   nzzv            Integer parameter number of output z-level vector fields
!   izzu_*          Integer parameter index for output z-level vector X field *
!   izzv_*          Integer parameter index for output z-level vector Y field *
!   jpdszzu         Integer (1:24,nzzv) parameter PDS type and level ids
!                   for output z-level vector X fields
!   jpdszzv         Integer (1:24,nzzv) parameter PDS type and level ids
!                   for output z-level vector Y fields
!   nths            Integer parameter number of output theta-level scalar fields
!   iths_*          Integer parameter index for output theta-level scalar field *
!   jpdsths         Integer (1:24,nths) parameter PDS type and level ids
!                   for output theta-level scalar fields
!   nthv            Integer parameter number of output theta-level vector fields
!   ithu_*          Integer parameter index for output theta-level vector X field *
!   ithv_*          Integer parameter index for output theta-level vector Y field *
!   jpdsthu         Integer (1:24,nthv) parameter PDS type and level ids
!                   for output theta-level vector X fields
!   jpdsthv         Integer (1:24,nthv) parameter PDS type and level ids
!                   for output theta-level vector Y fields
!   npvs            Integer parameter number of output pv-level scalar fields
!   ipvs_*          Integer parameter index for output pv-level scalar field *
!   jpdspvs         Integer (1:24,npvs) parameter PDS type and level ids
!                   for output pv-level scalar fields
!   npvv            Integer parameter number of output pv-level vector fields
!   ipvu_*          Integer parameter index for output pv-level vector X field *
!   ipvv_*          Integer parameter index for output pv-level vector Y field *
!   jpdspvu         Integer (1:24,npvv) parameter PDS type and level ids
!                   for output pv-level vector X fields
!   jpdspvv         Integer (1:24,npvv) parameter PDS type and level ids
!                   for output pv-level vector Y fields
!
! Attributes:
!   Language: Fortran 90
!
!$$$
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  input flux fields (scalar)
  integer,parameter:: nflxs                   =62
  integer,parameter:: iflxs_shtfl_sfc         =01
  integer,parameter:: iflxs_lhtfl_sfc         =02
  integer,parameter:: iflxs_tmp_sfc           =03
  integer,parameter:: iflxs_soilw_dlr_0_10    =04
  integer,parameter:: iflxs_soilw_dlr_10_40   =05
  integer,parameter:: iflxs_soilw_dlr_40_100  =06
  integer,parameter:: iflxs_soilw_dlr_100_200 =07
  integer,parameter:: iflxs_soilw_dlr_10_200  =08
  integer,parameter:: iflxs_tmp_dlr_0_10      =09
  integer,parameter:: iflxs_tmp_dlr_10_40     =10
  integer,parameter:: iflxs_tmp_dlr_40_100    =11
  integer,parameter:: iflxs_tmp_dlr_100_200   =12
  integer,parameter:: iflxs_tmp_dlr_10_200    =13
  integer,parameter:: iflxs_soill_dlr_0_10    =14
  integer,parameter:: iflxs_soill_dlr_10_40   =15
  integer,parameter:: iflxs_soill_dlr_40_100  =16
  integer,parameter:: iflxs_soill_dlr_100_200 =17
  integer,parameter:: iflxs_cnwat_sfc         =18
  integer,parameter:: iflxs_snod_sfc          =19
  integer,parameter:: iflxs_weasd_sfc         =20
  integer,parameter:: iflxs_dlwrf_sfc         =21
  integer,parameter:: iflxs_ulwrf_sfc         =22
  integer,parameter:: iflxs_ulwrf_toa         =23
  integer,parameter:: iflxs_uswrf_toa         =24
  integer,parameter:: iflxs_uswrf_sfc         =25
  integer,parameter:: iflxs_dswrf_sfc         =26
  integer,parameter:: iflxs_duvb_sfc          =27
  integer,parameter:: iflxs_cduvb_sfc         =28
  integer,parameter:: iflxs_tcdc_hcl          =29
  integer,parameter:: iflxs_pres_hct          =30
  integer,parameter:: iflxs_pres_hcb          =31
  integer,parameter:: iflxs_tmp_hct           =32
  integer,parameter:: iflxs_tcdc_mcl          =33
  integer,parameter:: iflxs_pres_mct          =34
  integer,parameter:: iflxs_pres_mcb          =35
  integer,parameter:: iflxs_tmp_mct           =36
  integer,parameter:: iflxs_tcdc_lcl          =37
  integer,parameter:: iflxs_pres_lct          =38
  integer,parameter:: iflxs_pres_lcb          =39
  integer,parameter:: iflxs_tmp_lct           =40
  integer,parameter:: iflxs_prate_sfc         =41
  integer,parameter:: iflxs_cprat_sfc         =42
  integer,parameter:: iflxs_gflux_sfc         =43
  integer,parameter:: iflxs_land_sfc          =44
  integer,parameter:: iflxs_icec_sfc          =45
  integer,parameter:: iflxs_icetk_sfc         =46
  integer,parameter:: iflxs_tmp_hag_2         =47
  integer,parameter:: iflxs_spfh_hag_2        =48
  integer,parameter:: iflxs_pres_sfc          =49
  integer,parameter:: iflxs_tmax_hag_2        =50
  integer,parameter:: iflxs_tmin_hag_2        =51
  integer,parameter:: iflxs_watr_sfc          =52
  integer,parameter:: iflxs_pevpr_sfc         =53
  integer,parameter:: iflxs_cwork_clm         =54
  integer,parameter:: iflxs_hpbl_sfc          =55
  integer,parameter:: iflxs_pwat_clm          =56
  integer,parameter:: iflxs_albdo_sfc         =57
  integer,parameter:: iflxs_tcdc_clm          =58
  integer,parameter:: iflxs_tcdc_cvc          =59
  integer,parameter:: iflxs_pres_cvt          =60
  integer,parameter:: iflxs_pres_cvb          =61
  integer,parameter:: iflxs_tcdc_blc          =62
  integer jpdsflxs(1:24,nflxs)
  data jpdsflxs(19,iflxs_shtfl_sfc         ) /002/
  data jpdsflxs(19,iflxs_lhtfl_sfc         ) /002/
  data jpdsflxs(19,iflxs_tmp_sfc           ) /002/
  data jpdsflxs(19,iflxs_soilw_dlr_0_10    ) /002/
  data jpdsflxs(19,iflxs_soilw_dlr_10_40   ) /002/
  data jpdsflxs(19,iflxs_soilw_dlr_40_100  ) /002/
  data jpdsflxs(19,iflxs_soilw_dlr_100_200 ) /002/
  data jpdsflxs(19,iflxs_soilw_dlr_10_200  ) /002/
  data jpdsflxs(19,iflxs_tmp_dlr_0_10      ) /002/
  data jpdsflxs(19,iflxs_tmp_dlr_10_40     ) /002/
  data jpdsflxs(19,iflxs_tmp_dlr_40_100    ) /002/
  data jpdsflxs(19,iflxs_tmp_dlr_100_200   ) /002/
  data jpdsflxs(19,iflxs_tmp_dlr_10_200    ) /002/
  data jpdsflxs(19,iflxs_soill_dlr_0_10    ) /130/
  data jpdsflxs(19,iflxs_soill_dlr_10_40   ) /130/
  data jpdsflxs(19,iflxs_soill_dlr_40_100  ) /130/
  data jpdsflxs(19,iflxs_soill_dlr_100_200 ) /130/
  data jpdsflxs(19,iflxs_cnwat_sfc         ) /002/
  data jpdsflxs(19,iflxs_snod_sfc          ) /002/
  data jpdsflxs(19,iflxs_weasd_sfc         ) /002/
  data jpdsflxs(19,iflxs_dlwrf_sfc         ) /002/
  data jpdsflxs(19,iflxs_ulwrf_sfc         ) /002/
  data jpdsflxs(19,iflxs_ulwrf_toa         ) /002/
  data jpdsflxs(19,iflxs_uswrf_toa         ) /002/
  data jpdsflxs(19,iflxs_uswrf_sfc         ) /002/
  data jpdsflxs(19,iflxs_dswrf_sfc         ) /002/
  data jpdsflxs(19,iflxs_duvb_sfc          ) /129/
  data jpdsflxs(19,iflxs_cduvb_sfc         ) /129/
  data jpdsflxs(19,iflxs_tcdc_hcl          ) /002/
  data jpdsflxs(19,iflxs_pres_hct          ) /002/
  data jpdsflxs(19,iflxs_pres_hcb          ) /002/
  data jpdsflxs(19,iflxs_tmp_hct           ) /002/
  data jpdsflxs(19,iflxs_tcdc_mcl          ) /002/
  data jpdsflxs(19,iflxs_pres_mct          ) /002/
  data jpdsflxs(19,iflxs_pres_mcb          ) /002/
  data jpdsflxs(19,iflxs_tmp_mct           ) /002/
  data jpdsflxs(19,iflxs_tcdc_lcl          ) /002/
  data jpdsflxs(19,iflxs_pres_lct          ) /002/
  data jpdsflxs(19,iflxs_pres_lcb          ) /002/
  data jpdsflxs(19,iflxs_tmp_lct           ) /002/
  data jpdsflxs(19,iflxs_prate_sfc         ) /002/
  data jpdsflxs(19,iflxs_cprat_sfc         ) /002/
  data jpdsflxs(19,iflxs_gflux_sfc         ) /002/
  data jpdsflxs(19,iflxs_land_sfc          ) /002/
  data jpdsflxs(19,iflxs_icec_sfc          ) /002/
  data jpdsflxs(19,iflxs_icetk_sfc         ) /002/
  data jpdsflxs(19,iflxs_tmp_hag_2         ) /002/
  data jpdsflxs(19,iflxs_spfh_hag_2        ) /002/
  data jpdsflxs(19,iflxs_pres_sfc          ) /002/
  data jpdsflxs(19,iflxs_tmax_hag_2        ) /002/
  data jpdsflxs(19,iflxs_tmin_hag_2        ) /002/
  data jpdsflxs(19,iflxs_watr_sfc          ) /002/
  data jpdsflxs(19,iflxs_pevpr_sfc         ) /002/
  data jpdsflxs(19,iflxs_cwork_clm         ) /002/
  data jpdsflxs(19,iflxs_hpbl_sfc          ) /002/
  data jpdsflxs(19,iflxs_pwat_clm          ) /002/
  data jpdsflxs(19,iflxs_albdo_sfc         ) /002/
  data jpdsflxs(19,iflxs_tcdc_clm          ) /002/
  data jpdsflxs(19,iflxs_tcdc_cvc          ) /002/
  data jpdsflxs(19,iflxs_pres_cvt          ) /002/
  data jpdsflxs(19,iflxs_pres_cvb          ) /002/
  data jpdsflxs(19,iflxs_tcdc_blc          ) /002/
  data jpdsflxs(5:7,iflxs_shtfl_sfc         ) /122,001,000/
  data jpdsflxs(5:7,iflxs_lhtfl_sfc         ) /121,001,000/
  data jpdsflxs(5:7,iflxs_tmp_sfc           ) /011,001,000/
  data jpdsflxs(5:7,iflxs_soilw_dlr_0_10    ) /144,112,010/
  data jpdsflxs(5:7,iflxs_soilw_dlr_10_40   ) /144,112,2600/
  data jpdsflxs(5:7,iflxs_soilw_dlr_40_100  ) /144,112,10340/
  data jpdsflxs(5:7,iflxs_soilw_dlr_100_200 ) /144,112,25800/
  data jpdsflxs(5:7,iflxs_soilw_dlr_10_200  ) /144,112,2760/
  data jpdsflxs(5:7,iflxs_tmp_dlr_0_10      ) /011,112,010/
  data jpdsflxs(5:7,iflxs_tmp_dlr_10_40     ) /011,112,2600/
  data jpdsflxs(5:7,iflxs_tmp_dlr_40_100    ) /011,112,10340/
  data jpdsflxs(5:7,iflxs_tmp_dlr_100_200   ) /011,112,25800/
  data jpdsflxs(5:7,iflxs_tmp_dlr_10_200    ) /011,112,2760/
  data jpdsflxs(5:7,iflxs_soill_dlr_0_10    ) /160,112,010/
  data jpdsflxs(5:7,iflxs_soill_dlr_10_40   ) /160,112,2600/
  data jpdsflxs(5:7,iflxs_soill_dlr_40_100  ) /160,112,10340/
  data jpdsflxs(5:7,iflxs_soill_dlr_100_200 ) /160,112,25800/
  data jpdsflxs(5:7,iflxs_cnwat_sfc         ) /223,001,000/
  data jpdsflxs(5:7,iflxs_snod_sfc          ) /066,001,000/
  data jpdsflxs(5:7,iflxs_weasd_sfc         ) /065,001,000/
  data jpdsflxs(5:7,iflxs_dlwrf_sfc         ) /205,001,000/
  data jpdsflxs(5:7,iflxs_ulwrf_sfc         ) /212,001,000/
  data jpdsflxs(5:7,iflxs_ulwrf_toa         ) /212,008,000/
  data jpdsflxs(5:7,iflxs_uswrf_toa         ) /211,008,000/
  data jpdsflxs(5:7,iflxs_uswrf_sfc         ) /211,001,000/
  data jpdsflxs(5:7,iflxs_dswrf_sfc         ) /204,001,000/
  data jpdsflxs(5:7,iflxs_duvb_sfc          ) /200,001,000/
  data jpdsflxs(5:7,iflxs_cduvb_sfc         ) /201,001,000/
  data jpdsflxs(5:7,iflxs_tcdc_hcl          ) /071,234,000/
  data jpdsflxs(5:7,iflxs_pres_hct          ) /001,233,000/
  data jpdsflxs(5:7,iflxs_pres_hcb          ) /001,232,000/
  data jpdsflxs(5:7,iflxs_tmp_hct           ) /011,233,000/
  data jpdsflxs(5:7,iflxs_tcdc_mcl          ) /071,224,000/
  data jpdsflxs(5:7,iflxs_pres_mct          ) /001,223,000/
  data jpdsflxs(5:7,iflxs_pres_mcb          ) /001,222,000/
  data jpdsflxs(5:7,iflxs_tmp_mct           ) /011,223,000/
  data jpdsflxs(5:7,iflxs_tcdc_lcl          ) /071,214,000/
  data jpdsflxs(5:7,iflxs_pres_lct          ) /001,213,000/
  data jpdsflxs(5:7,iflxs_pres_lcb          ) /001,212,000/
  data jpdsflxs(5:7,iflxs_tmp_lct           ) /011,213,000/
  data jpdsflxs(5:7,iflxs_prate_sfc         ) /059,001,000/
  data jpdsflxs(5:7,iflxs_cprat_sfc         ) /214,001,000/
  data jpdsflxs(5:7,iflxs_gflux_sfc         ) /155,001,000/
  data jpdsflxs(5:7,iflxs_land_sfc          ) /081,001,000/
  data jpdsflxs(5:7,iflxs_icec_sfc          ) /091,001,000/
  data jpdsflxs(5:7,iflxs_icetk_sfc         ) /092,001,000/
  data jpdsflxs(5:7,iflxs_tmp_hag_2         ) /011,105,002/
  data jpdsflxs(5:7,iflxs_spfh_hag_2        ) /051,105,002/
  data jpdsflxs(5:7,iflxs_pres_sfc          ) /001,001,000/
  data jpdsflxs(5:7,iflxs_tmax_hag_2        ) /015,105,002/
  data jpdsflxs(5:7,iflxs_tmin_hag_2        ) /016,105,002/
  data jpdsflxs(5:7,iflxs_watr_sfc          ) /090,001,000/
  data jpdsflxs(5:7,iflxs_pevpr_sfc         ) /145,001,000/
  data jpdsflxs(5:7,iflxs_cwork_clm         ) /146,200,000/
  data jpdsflxs(5:7,iflxs_hpbl_sfc          ) /221,001,000/
  data jpdsflxs(5:7,iflxs_pwat_clm          ) /054,200,000/
  data jpdsflxs(5:7,iflxs_albdo_sfc         ) /084,001,000/
  data jpdsflxs(5:7,iflxs_tcdc_clm          ) /071,200,000/
  data jpdsflxs(5:7,iflxs_tcdc_cvc          ) /071,244,000/
  data jpdsflxs(5:7,iflxs_pres_cvt          ) /001,243,000/
  data jpdsflxs(5:7,iflxs_pres_cvb          ) /001,242,000/
  data jpdsflxs(5:7,iflxs_tcdc_blc          ) /071,211,000/
!  input flux fields (vector)
  integer,parameter:: nflxv                   = 3
  integer,parameter:: iflxu_uflx_sfc          =01
  integer,parameter:: iflxv_vflx_sfc          =01
  integer,parameter:: iflxu_ugrd_hag_10       =02
  integer,parameter:: iflxv_vgrd_hag_10       =02
  integer,parameter:: iflxu_ugwd_sfc          =03
  integer,parameter:: iflxv_vgwd_sfc          =03
  integer jpdsflxu(1:24,nflxv),jpdsflxv(1:24,nflxv)
  data jpdsflxu(19,:) /nflxv*-1/
  data jpdsflxv(19,:) /nflxv*-1/
  data jpdsflxu(5:7,iflxu_uflx_sfc          ) /124,001,000/
  data jpdsflxv(5:7,iflxv_vflx_sfc          ) /125,001,000/
  data jpdsflxu(5:7,iflxu_ugrd_hag_10       ) /033,105,010/
  data jpdsflxv(5:7,iflxv_vgrd_hag_10       ) /034,105,010/
  data jpdsflxu(5:7,iflxu_ugwd_sfc          ) /147,001,000/
  data jpdsflxv(5:7,iflxv_vgwd_sfc          ) /148,001,000/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  output flux fields (scalar)
  integer,parameter:: nfoxs                   =67
  integer,parameter:: ifoxs_shtfl_sfc         =01
  integer,parameter:: ifoxs_lhtfl_sfc         =02
  integer,parameter:: ifoxs_tmp_sfc           =03
  integer,parameter:: ifoxs_soilw_dlr_0_10    =04
  integer,parameter:: ifoxs_soilw_dlr_10_40   =05
  integer,parameter:: ifoxs_soilw_dlr_40_100  =06
  integer,parameter:: ifoxs_soilw_dlr_100_200 =07
  integer,parameter:: ifoxs_soilw_dlr_10_200  =08
  integer,parameter:: ifoxs_tmp_dlr_0_10      =09
  integer,parameter:: ifoxs_tmp_dlr_10_40     =10
  integer,parameter:: ifoxs_tmp_dlr_40_100    =11
  integer,parameter:: ifoxs_tmp_dlr_100_200   =12
  integer,parameter:: ifoxs_tmp_dlr_10_200    =13
  integer,parameter:: ifoxs_soill_dlr_0_10    =14
  integer,parameter:: ifoxs_soill_dlr_10_40   =15
  integer,parameter:: ifoxs_soill_dlr_40_100  =16
  integer,parameter:: ifoxs_soill_dlr_100_200 =17
  integer,parameter:: ifoxs_cnwat_sfc         =18
  integer,parameter:: ifoxs_snod_sfc          =19
  integer,parameter:: ifoxs_weasd_sfc         =20
  integer,parameter:: ifoxs_dlwrf_sfc         =21
  integer,parameter:: ifoxs_ulwrf_sfc         =22
  integer,parameter:: ifoxs_ulwrf_toa         =23
  integer,parameter:: ifoxs_uswrf_toa         =24
  integer,parameter:: ifoxs_uswrf_sfc         =25
  integer,parameter:: ifoxs_dswrf_sfc         =26
  integer,parameter:: ifoxs_duvb_sfc          =27
  integer,parameter:: ifoxs_cduvb_sfc         =28
  integer,parameter:: ifoxs_tcdc_hcl          =29
  integer,parameter:: ifoxs_pres_hct          =30
  integer,parameter:: ifoxs_pres_hcb          =31
  integer,parameter:: ifoxs_tmp_hct           =32
  integer,parameter:: ifoxs_tcdc_mcl          =33
  integer,parameter:: ifoxs_pres_mct          =34
  integer,parameter:: ifoxs_pres_mcb          =35
  integer,parameter:: ifoxs_tmp_mct           =36
  integer,parameter:: ifoxs_tcdc_lcl          =37
  integer,parameter:: ifoxs_pres_lct          =38
  integer,parameter:: ifoxs_pres_lcb          =39
  integer,parameter:: ifoxs_tmp_lct           =40
  integer,parameter:: ifoxs_prate_sfc         =41
  integer,parameter:: ifoxs_cprat_sfc         =42
  integer,parameter:: ifoxs_gflux_sfc         =43
  integer,parameter:: ifoxs_land_sfc          =44
  integer,parameter:: ifoxs_icec_sfc          =45
  integer,parameter:: ifoxs_icetk_sfc         =46
  integer,parameter:: ifoxs_tmp_hag_2         =47
  integer,parameter:: ifoxs_spfh_hag_2        =48
  integer,parameter:: ifoxs_tmax_hag_2        =49
  integer,parameter:: ifoxs_tmin_hag_2        =50
  integer,parameter:: ifoxs_watr_sfc          =51
  integer,parameter:: ifoxs_pevpr_sfc         =52
  integer,parameter:: ifoxs_cwork_clm         =53
  integer,parameter:: ifoxs_hpbl_sfc          =54
  integer,parameter:: ifoxs_albdo_sfc         =55
  integer,parameter:: ifoxs_tcdc_clm          =56
  integer,parameter:: ifoxs_tcdc_cvc          =57
  integer,parameter:: ifoxs_pres_cvt          =58
  integer,parameter:: ifoxs_pres_cvb          =59
  integer,parameter:: ifoxs_tcdc_blc          =60
  integer,parameter:: ifoxs_apcp_sfc          =61
  integer,parameter:: ifoxs_acpcp_sfc         =62
  integer,parameter:: ifoxs_crain_sfc         =63
  integer,parameter:: ifoxs_cfrzr_sfc         =64
  integer,parameter:: ifoxs_cicep_sfc         =65
  integer,parameter:: ifoxs_csnow_sfc         =66
  integer,parameter:: ifoxs_rh_hag_2          =67
  integer jpdsfoxs(1:24,nfoxs)
  data jpdsfoxs(19,ifoxs_shtfl_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_lhtfl_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_tmp_sfc           ) /002/
  data jpdsfoxs(19,ifoxs_soilw_dlr_0_10    ) /002/
  data jpdsfoxs(19,ifoxs_soilw_dlr_10_40   ) /002/
  data jpdsfoxs(19,ifoxs_soilw_dlr_40_100  ) /002/
  data jpdsfoxs(19,ifoxs_soilw_dlr_100_200 ) /002/
  data jpdsfoxs(19,ifoxs_soilw_dlr_10_200  ) /002/
  data jpdsfoxs(19,ifoxs_tmp_dlr_0_10      ) /002/
  data jpdsfoxs(19,ifoxs_tmp_dlr_10_40     ) /002/
  data jpdsfoxs(19,ifoxs_tmp_dlr_40_100    ) /002/
  data jpdsfoxs(19,ifoxs_tmp_dlr_100_200   ) /002/
  data jpdsfoxs(19,ifoxs_tmp_dlr_10_200    ) /002/
  data jpdsfoxs(19,ifoxs_soill_dlr_0_10    ) /130/
  data jpdsfoxs(19,ifoxs_soill_dlr_10_40   ) /130/
  data jpdsfoxs(19,ifoxs_soill_dlr_40_100  ) /130/
  data jpdsfoxs(19,ifoxs_soill_dlr_100_200 ) /130/
  data jpdsfoxs(19,ifoxs_cnwat_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_snod_sfc          ) /002/
  data jpdsfoxs(19,ifoxs_weasd_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_dlwrf_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_ulwrf_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_ulwrf_toa         ) /002/
  data jpdsfoxs(19,ifoxs_uswrf_toa         ) /002/
  data jpdsfoxs(19,ifoxs_uswrf_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_dswrf_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_duvb_sfc          ) /129/
  data jpdsfoxs(19,ifoxs_cduvb_sfc         ) /129/
  data jpdsfoxs(19,ifoxs_tcdc_hcl          ) /002/
  data jpdsfoxs(19,ifoxs_pres_hct          ) /002/
  data jpdsfoxs(19,ifoxs_pres_hcb          ) /002/
  data jpdsfoxs(19,ifoxs_tmp_hct           ) /002/
  data jpdsfoxs(19,ifoxs_tcdc_mcl          ) /002/
  data jpdsfoxs(19,ifoxs_pres_mct          ) /002/
  data jpdsfoxs(19,ifoxs_pres_mcb          ) /002/
  data jpdsfoxs(19,ifoxs_tmp_mct           ) /002/
  data jpdsfoxs(19,ifoxs_tcdc_lcl          ) /002/
  data jpdsfoxs(19,ifoxs_pres_lct          ) /002/
  data jpdsfoxs(19,ifoxs_pres_lcb          ) /002/
  data jpdsfoxs(19,ifoxs_tmp_lct           ) /002/
  data jpdsfoxs(19,ifoxs_prate_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_cprat_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_gflux_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_land_sfc          ) /002/
  data jpdsfoxs(19,ifoxs_icec_sfc          ) /002/
  data jpdsfoxs(19,ifoxs_icetk_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_tmp_hag_2         ) /002/
  data jpdsfoxs(19,ifoxs_spfh_hag_2        ) /002/
  data jpdsfoxs(19,ifoxs_tmax_hag_2        ) /002/
  data jpdsfoxs(19,ifoxs_tmin_hag_2        ) /002/
  data jpdsfoxs(19,ifoxs_watr_sfc          ) /002/
  data jpdsfoxs(19,ifoxs_pevpr_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_cwork_clm         ) /002/
  data jpdsfoxs(19,ifoxs_hpbl_sfc          ) /002/
  data jpdsfoxs(19,ifoxs_albdo_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_tcdc_clm          ) /002/
  data jpdsfoxs(19,ifoxs_tcdc_cvc          ) /002/
  data jpdsfoxs(19,ifoxs_pres_cvt          ) /002/
  data jpdsfoxs(19,ifoxs_pres_cvb          ) /002/
  data jpdsfoxs(19,ifoxs_tcdc_blc          ) /002/
  data jpdsfoxs(19,ifoxs_apcp_sfc          ) /002/
  data jpdsfoxs(19,ifoxs_acpcp_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_crain_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_cfrzr_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_cicep_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_csnow_sfc         ) /002/
  data jpdsfoxs(19,ifoxs_rh_hag_2          ) /002/
  data jpdsfoxs(5:7,ifoxs_shtfl_sfc         ) /122,001,000/
  data jpdsfoxs(5:7,ifoxs_lhtfl_sfc         ) /121,001,000/
  data jpdsfoxs(5:7,ifoxs_tmp_sfc           ) /011,001,000/
  data jpdsfoxs(5:7,ifoxs_soilw_dlr_0_10    ) /144,112,010/
  data jpdsfoxs(5:7,ifoxs_soilw_dlr_10_40   ) /144,112,2600/
  data jpdsfoxs(5:7,ifoxs_soilw_dlr_40_100  ) /144,112,10340/
  data jpdsfoxs(5:7,ifoxs_soilw_dlr_100_200 ) /144,112,25800/
  data jpdsfoxs(5:7,ifoxs_soilw_dlr_10_200  ) /144,112,2760/
  data jpdsfoxs(5:7,ifoxs_tmp_dlr_0_10      ) /011,112,010/
  data jpdsfoxs(5:7,ifoxs_tmp_dlr_10_40     ) /011,112,2600/
  data jpdsfoxs(5:7,ifoxs_tmp_dlr_40_100    ) /011,112,10340/
  data jpdsfoxs(5:7,ifoxs_tmp_dlr_100_200   ) /011,112,25800/
  data jpdsfoxs(5:7,ifoxs_tmp_dlr_10_200    ) /011,112,2760/
  data jpdsfoxs(5:7,ifoxs_soill_dlr_0_10    ) /160,112,010/
  data jpdsfoxs(5:7,ifoxs_soill_dlr_10_40   ) /160,112,2600/
  data jpdsfoxs(5:7,ifoxs_soill_dlr_40_100  ) /160,112,10340/
  data jpdsfoxs(5:7,ifoxs_soill_dlr_100_200 ) /160,112,25800/
  data jpdsfoxs(5:7,ifoxs_cnwat_sfc         ) /223,001,000/
  data jpdsfoxs(5:7,ifoxs_snod_sfc          ) /066,001,000/
  data jpdsfoxs(5:7,ifoxs_weasd_sfc         ) /065,001,000/
  data jpdsfoxs(5:7,ifoxs_dlwrf_sfc         ) /205,001,000/
  data jpdsfoxs(5:7,ifoxs_ulwrf_sfc         ) /212,001,000/
  data jpdsfoxs(5:7,ifoxs_ulwrf_toa         ) /212,008,000/
  data jpdsfoxs(5:7,ifoxs_uswrf_toa         ) /211,008,000/
  data jpdsfoxs(5:7,ifoxs_uswrf_sfc         ) /211,001,000/
  data jpdsfoxs(5:7,ifoxs_dswrf_sfc         ) /204,001,000/
  data jpdsfoxs(5:7,ifoxs_duvb_sfc          ) /200,001,000/
  data jpdsfoxs(5:7,ifoxs_cduvb_sfc         ) /201,001,000/
  data jpdsfoxs(5:7,ifoxs_tcdc_hcl          ) /071,234,000/
  data jpdsfoxs(5:7,ifoxs_pres_hct          ) /001,233,000/
  data jpdsfoxs(5:7,ifoxs_pres_hcb          ) /001,232,000/
  data jpdsfoxs(5:7,ifoxs_tmp_hct           ) /011,233,000/
  data jpdsfoxs(5:7,ifoxs_tcdc_mcl          ) /071,224,000/
  data jpdsfoxs(5:7,ifoxs_pres_mct          ) /001,223,000/
  data jpdsfoxs(5:7,ifoxs_pres_mcb          ) /001,222,000/
  data jpdsfoxs(5:7,ifoxs_tmp_mct           ) /011,223,000/
  data jpdsfoxs(5:7,ifoxs_tcdc_lcl          ) /071,214,000/
  data jpdsfoxs(5:7,ifoxs_pres_lct          ) /001,213,000/
  data jpdsfoxs(5:7,ifoxs_pres_lcb          ) /001,212,000/
  data jpdsfoxs(5:7,ifoxs_tmp_lct           ) /011,213,000/
  data jpdsfoxs(5:7,ifoxs_prate_sfc         ) /059,001,000/
  data jpdsfoxs(5:7,ifoxs_cprat_sfc         ) /214,001,000/
  data jpdsfoxs(5:7,ifoxs_gflux_sfc         ) /155,001,000/
  data jpdsfoxs(5:7,ifoxs_land_sfc          ) /081,001,000/
  data jpdsfoxs(5:7,ifoxs_icec_sfc          ) /091,001,000/
  data jpdsfoxs(5:7,ifoxs_icetk_sfc         ) /092,001,000/
  data jpdsfoxs(5:7,ifoxs_tmp_hag_2         ) /011,105,002/
  data jpdsfoxs(5:7,ifoxs_spfh_hag_2        ) /051,105,002/
  data jpdsfoxs(5:7,ifoxs_tmax_hag_2        ) /015,105,002/
  data jpdsfoxs(5:7,ifoxs_tmin_hag_2        ) /016,105,002/
  data jpdsfoxs(5:7,ifoxs_watr_sfc          ) /090,001,000/
  data jpdsfoxs(5:7,ifoxs_pevpr_sfc         ) /145,001,000/
  data jpdsfoxs(5:7,ifoxs_cwork_clm         ) /146,200,000/
  data jpdsfoxs(5:7,ifoxs_hpbl_sfc          ) /221,001,000/
  data jpdsfoxs(5:7,ifoxs_albdo_sfc         ) /084,001,000/
  data jpdsfoxs(5:7,ifoxs_tcdc_clm          ) /071,200,000/
  data jpdsfoxs(5:7,ifoxs_tcdc_cvc          ) /071,244,000/
  data jpdsfoxs(5:7,ifoxs_pres_cvt          ) /001,243,000/
  data jpdsfoxs(5:7,ifoxs_pres_cvb          ) /001,242,000/
  data jpdsfoxs(5:7,ifoxs_tcdc_blc          ) /071,211,000/
  data jpdsfoxs(5:7,ifoxs_apcp_sfc          ) /061,001,000/
  data jpdsfoxs(5:7,ifoxs_acpcp_sfc         ) /063,001,000/
  data jpdsfoxs(5:7,ifoxs_crain_sfc         ) /140,001,000/
  data jpdsfoxs(5:7,ifoxs_cfrzr_sfc         ) /141,001,000/
  data jpdsfoxs(5:7,ifoxs_cicep_sfc         ) /142,001,000/
  data jpdsfoxs(5:7,ifoxs_csnow_sfc         ) /143,001,000/
  data jpdsfoxs(5:7,ifoxs_rh_hag_2          ) /052,105,002/
!  output flux fields (vector)
  integer,parameter:: nfoxv                   = 3
  integer,parameter:: ifoxu_uflx_sfc          =01
  integer,parameter:: ifoxv_vflx_sfc          =01
  integer,parameter:: ifoxu_ugrd_hag_10       =02
  integer,parameter:: ifoxv_vgrd_hag_10       =02
  integer,parameter:: ifoxu_ugwd_sfc          =03
  integer,parameter:: ifoxv_vgwd_sfc          =03
  integer jpdsfoxu(1:24,nfoxv),jpdsfoxv(1:24,nfoxv)
  data jpdsfoxu(19,:) /nfoxv*2/
  data jpdsfoxv(19,:) /nfoxv*2/
  data jpdsfoxu(5:7,ifoxu_uflx_sfc          ) /124,001,000/
  data jpdsfoxv(5:7,ifoxv_vflx_sfc          ) /125,001,000/
  data jpdsfoxu(5:7,ifoxu_ugrd_hag_10       ) /033,105,010/
  data jpdsfoxv(5:7,ifoxv_vgrd_hag_10       ) /034,105,010/
  data jpdsfoxu(5:7,ifoxu_ugwd_sfc          ) /147,001,000/
  data jpdsfoxv(5:7,ifoxv_vgwd_sfc          ) /148,001,000/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  output sundry fields (scalar)
  integer,parameter:: nsuns                   =32
  integer,parameter:: isuns_pres_sfc          =01
  integer,parameter:: isuns_pwat_clm          =02
  integer,parameter:: isuns_rh_clm            =03
  integer,parameter:: isuns_hgt_trp           =04
  integer,parameter:: isuns_tmp_trp           =05
  integer,parameter:: isuns_pres_trp          =06
  integer,parameter:: isuns_vwsh_trp          =07
  integer,parameter:: isuns_lftx_sfc          =08
  integer,parameter:: isuns_cape_sfc          =09
  integer,parameter:: isuns_cin_sfc           =10
  integer,parameter:: isuns_blftx_sfc         =11
  integer,parameter:: isuns_cape_plg_180_0    =12
  integer,parameter:: isuns_cin_plg_180_0     =13
  integer,parameter:: isuns_hgt_mwl           =14
  integer,parameter:: isuns_tmp_mwl           =15
  integer,parameter:: isuns_pres_mwl          =16
  integer,parameter:: isuns_hgt_sfc           =17
  integer,parameter:: isuns_prmsl_msl         =18
  integer,parameter:: isuns_rh_slr_044_100    =19
  integer,parameter:: isuns_rh_slr_072_094    =20
  integer,parameter:: isuns_rh_slr_044_072    =21
  integer,parameter:: isuns_rh_slr_033_100    =22
  integer,parameter:: isuns_pot_sig_9950      =23
  integer,parameter:: isuns_tmp_sig_9950      =24
  integer,parameter:: isuns_vvel_sig_9950     =25
  integer,parameter:: isuns_rh_sig_9950       =26
  integer,parameter:: isuns_tozne_clm         =27
  integer,parameter:: isuns_cwat_clm          =28
  integer,parameter:: isuns_hgt_zdeg          =29
  integer,parameter:: isuns_rh_zdeg           =30
  integer,parameter:: isuns_hgt_htfl          =31
  integer,parameter:: isuns_rh_htfl           =32
  integer jpdssuns(1:24,nsuns)
  data jpdssuns(19,:) /nsuns*2/
  data jpdssuns(5:7,isuns_pres_sfc          ) /001,001,000/
  data jpdssuns(5:7,isuns_pwat_clm          ) /054,200,000/
  data jpdssuns(5:7,isuns_rh_clm            ) /052,200,000/
  data jpdssuns(5:7,isuns_hgt_trp           ) /007,007,000/
  data jpdssuns(5:7,isuns_tmp_trp           ) /011,007,000/
  data jpdssuns(5:7,isuns_pres_trp          ) /001,007,000/
  data jpdssuns(5:7,isuns_vwsh_trp          ) /136,007,000/
  data jpdssuns(5:7,isuns_lftx_sfc          ) /131,001,000/
  data jpdssuns(5:7,isuns_cape_sfc          ) /157,001,000/
  data jpdssuns(5:7,isuns_cin_sfc           ) /156,001,000/
  data jpdssuns(5:7,isuns_blftx_sfc         ) /132,001,000/
  data jpdssuns(5:7,isuns_cape_plg_180_0    ) /157,116,46080/
  data jpdssuns(5:7,isuns_cin_plg_180_0     ) /156,116,46080/
  data jpdssuns(5:7,isuns_hgt_mwl           ) /007,006,000/
  data jpdssuns(5:7,isuns_tmp_mwl           ) /011,006,000/
  data jpdssuns(5:7,isuns_pres_mwl          ) /001,006,000/
  data jpdssuns(5:7,isuns_hgt_sfc           ) /007,001,000/
  data jpdssuns(5:7,isuns_prmsl_msl         ) /002,102,000/
  data jpdssuns(5:7,isuns_rh_slr_044_100    ) /052,108,11364/
  data jpdssuns(5:7,isuns_rh_slr_072_094    ) /052,108,18526/
  data jpdssuns(5:7,isuns_rh_slr_044_072    ) /052,108,11336/
  data jpdssuns(5:7,isuns_rh_slr_033_100    ) /052,108,8548/
  data jpdssuns(5:7,isuns_pot_sig_9950      ) /013,107,9950/
  data jpdssuns(5:7,isuns_tmp_sig_9950      ) /011,107,9950/
  data jpdssuns(5:7,isuns_vvel_sig_9950     ) /039,107,9950/
  data jpdssuns(5:7,isuns_rh_sig_9950       ) /052,107,9950/
  data jpdssuns(5:7,isuns_tozne_clm         ) /010,200,000/
  data jpdssuns(5:7,isuns_cwat_clm          ) /076,200,000/
  data jpdssuns(5:7,isuns_hgt_zdeg          ) /007,004,000/
  data jpdssuns(5:7,isuns_rh_zdeg           ) /052,004,000/
  data jpdssuns(5:7,isuns_hgt_htfl          ) /007,204,000/
  data jpdssuns(5:7,isuns_rh_htfl           ) /052,204,000/
!  output sundry fields (vector)
  integer,parameter:: nsunv                   = 3
  integer,parameter:: isunu_ugrd_trp          =01
  integer,parameter:: isunv_vgrd_trp          =01
  integer,parameter:: isunu_ugrd_mwl          =02
  integer,parameter:: isunv_vgrd_mwl          =02
  integer,parameter:: isunu_ugrd_sig_9950     =03
  integer,parameter:: isunv_vgrd_sig_9950     =03
  integer jpdssunu(1:24,nsunv),jpdssunv(1:24,nsunv)
  data jpdssunu(19,:) /nsunv*2/
  data jpdssunv(19,:) /nsunv*2/
  data jpdssunu(5:7,isunu_ugrd_trp          ) /033,007,000/
  data jpdssunv(5:7,isunv_vgrd_trp          ) /034,007,000/
  data jpdssunu(5:7,isunu_ugrd_mwl          ) /033,006,000/
  data jpdssunv(5:7,isunv_vgrd_mwl          ) /034,006,000/
  data jpdssunu(5:7,isunu_ugrd_sig_9950     ) /033,107,9950/
  data jpdssunv(5:7,isunv_vgrd_sig_9950     ) /034,107,9950/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  output pressure-level fields (scalar)
  integer,parameter:: npos                    = 9
  integer,parameter:: ipos_hgt_prs            =01
  integer,parameter:: ipos_tmp_prs            =02
  integer,parameter:: ipos_vvel_prs           =03
  integer,parameter:: ipos_rh_prs             =04
  integer,parameter:: ipos_absv_prs           =05
  integer,parameter:: ipos_spfh_prs           =06
  integer,parameter:: ipos_o3mr_prs           =07
  integer,parameter:: ipos_clwmr_prs          =08
  integer,parameter:: ipos_fwavh_prs          =09
  integer jpdspos(1:24,npos)
  data jpdspos(19,:) /npos*2/
  data jpdspos(5:6,ipos_hgt_prs             ) /007,100/
  data jpdspos(5:6,ipos_tmp_prs             ) /011,100/
  data jpdspos(5:6,ipos_vvel_prs            ) /039,100/
  data jpdspos(5:6,ipos_rh_prs              ) /052,100/
  data jpdspos(5:6,ipos_absv_prs            ) /041,100/
  data jpdspos(5:6,ipos_spfh_prs            ) /051,100/
  data jpdspos(5:6,ipos_o3mr_prs            ) /154,100/
  data jpdspos(5:6,ipos_clwmr_prs           ) /153,100/
  data jpdspos(5:6,ipos_fwavh_prs           ) /222,100/
!  output pressure-level fields (vector)
  integer,parameter:: npov                    = 1
  integer,parameter:: ipou_ugrd_prs           =01
  integer,parameter:: ipov_vgrd_prs           =01
  integer jpdspou(1:24,npov),jpdspov(1:24,npov)
  data jpdspou(19,:) /npov*2/
  data jpdspov(19,:) /npov*2/
  data jpdspou(5:6,ipou_ugrd_prs            ) /033,100/
  data jpdspov(5:6,ipov_vgrd_prs            ) /034,100/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  output pressure-thickness fields (scalar)
  integer,parameter:: npts                    = 3
  integer,parameter:: ipts_tmp_plg            =01
  integer,parameter:: ipts_rh_plg             =02
  integer,parameter:: ipts_spfh_plg           =03
  integer jpdspts(1:24,npts)
  data jpdspts(19,:) /npts*2/
  data jpdspts(5:6,ipts_tmp_plg             ) /011,116/
  data jpdspts(5:6,ipts_rh_plg              ) /052,116/
  data jpdspts(5:6,ipts_spfh_plg            ) /051,116/
!  output pressure-thickness fields (vector)
  integer,parameter:: nptv                    = 1
  integer,parameter:: iptu_ugrd_plg           =01
  integer,parameter:: iptv_vgrd_plg           =01
  integer jpdsptu(1:24,nptv),jpdsptv(1:24,nptv)
  data jpdsptu(19,:) /nptv*2/
  data jpdsptv(19,:) /nptv*2/
  data jpdsptu(5:6,iptu_ugrd_plg            ) /033,116/
  data jpdsptv(5:6,iptv_vgrd_plg            ) /034,116/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  output height-level fields (scalar)
  integer,parameter:: nzzs                    = 2
  integer,parameter:: izzs_tmp_hml            =01
  integer,parameter:: izzs_rh_hml             =02
  integer jpdszzs(1:24,nzzs)
  data jpdszzs(19,:) /nzzs*2/
  data jpdszzs(5:6,izzs_tmp_hml             ) /011,103/
  data jpdszzs(5:6,izzs_rh_hml              ) /052,103/
!  output height-level fields (vector)
  integer,parameter:: nzzv                    = 1
  integer,parameter:: izzu_ugrd_hml           =01
  integer,parameter:: izzv_vgrd_hml           =01
  integer jpdszzu(1:24,nzzv),jpdszzv(1:24,nzzv)
  data jpdszzu(19,:) /nzzv*2/
  data jpdszzv(19,:) /nzzv*2/
  data jpdszzu(5:6,izzu_ugrd_hml            ) /033,103/
  data jpdszzv(5:6,izzv_vgrd_hml            ) /034,103/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  output isentropic-level fields (scalar)
  integer,parameter:: nths                    = 3
  integer,parameter:: iths_mntsf_thel         =01
  integer,parameter:: iths_tmp_thel           =02
  integer,parameter:: iths_pvort_thel         =03
  integer jpdsths(1:24,nths)
  data jpdsths(19,:) /nths*2/
  data jpdsths(5:6,iths_mntsf_thel          ) /037,113/
  data jpdsths(5:6,iths_tmp_thel            ) /011,113/
  data jpdsths(5:6,iths_pvort_thel          ) /004,113/
!  output isentropic-level fields (vector)
  integer,parameter:: nthv                    = 1
  integer,parameter:: ithu_ugrd_thel          =01
  integer,parameter:: ithv_vgrd_thel          =01
  integer jpdsthu(1:24,nthv),jpdsthv(1:24,nthv)
  data jpdsthu(19,:) /nthv*2/
  data jpdsthv(19,:) /nthv*2/
  data jpdsthu(5:6,ithu_ugrd_thel           ) /033,117/
  data jpdsthv(5:6,ithv_vgrd_thel           ) /034,117/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  output potential-vorticity-level fields (scalar)
  integer,parameter:: npvs                    = 4
  integer,parameter:: ipvs_hgt_pvl            =01
  integer,parameter:: ipvs_tmp_pvl            =02
  integer,parameter:: ipvs_pres_pvl           =03
  integer,parameter:: ipvs_vwsh_pvl           =04
  integer jpdspvs(1:24,npvs)
  data jpdspvs(19,:) /npvs*2/
  data jpdspvs(5:6,ipvs_hgt_pvl             ) /007,117/
  data jpdspvs(5:6,ipvs_tmp_pvl             ) /011,117/
  data jpdspvs(5:6,ipvs_pres_pvl            ) /001,117/
  data jpdspvs(5:6,ipvs_vwsh_pvl            ) /136,117/
!  output potential-vorticity-level fields (vector)
  integer,parameter:: npvv                    = 1
  integer,parameter:: ipvu_ugrd_pvl           =01
  integer,parameter:: ipvv_vgrd_pvl           =01
  integer jpdspvu(1:24,npvv),jpdspvv(1:24,npvv)
  data jpdspvu(19,:) /npvv*2/
  data jpdspvv(19,:) /npvv*2/
  data jpdspvu(5:6,ipvu_ugrd_pvl            ) /033,117/
  data jpdspvv(5:6,ipvv_vgrd_pvl            ) /034,117/
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end module
