
! !module: gfs_dynamics_initialize_slg_mod 
!          --- initialize module of the
!              gridded component of the gfs dynamics system.
!              gfs dynamics related initilization 
!
! !description: gfs dynamics gridded component initialize module.
!
! !revision history:
!
!  november 2004  weiyu yang     initial code.
!  january 2006  s. moorthi      update to the new gfs version
!  august 2006   h. juang        add option to run generalized coordinates
!  december 2006 s. moorthi      gfsio included
!  january 2007 h. juang         change for dynamics only
!  May     2008 j. wang          change for gfs wrt grid component
!  Oct 04  2009 sarah lu         init xlon, xlat, lats_nodes_a_fix
!  Oct 05  2009 sarah lu         grid_gr unfolded from 2D to 3D
!  Oct 16  2009 sarah lu         initialize gfs_dyn_tracer
!  november 2009 j. wang         grid_gr_dfi for digital filter
!  Feb 05  2010 j. wang          add option to read in  restart file
!  Aug 19  2010 S. Moorthi       Updated for T574 + added num_reduce to namelist
!  Aug 25  2010 sarah lu         add option to compute tracer global sum
!  Sep 08  2010 J. Wang          changed gfsio file to nemsio file
!  Nov 01  2010 H. Juang         add non-iteration dimensional-split
!                                semi-Lagrangian dynamics (ndslfv)
!                                with mass_dp, process_split options.
!  Dec 16  2010 J. Wang          changed to nemsio library
!  Feb 20  2011 H. Juang         implement into nems for mass_dp and ndsl
!  Feb 28  2011 Sarah Lu         add thermodyn_id and sfcpress_id
!  Apr 06  2012 H. Juang         add idea
!  Sep 20  2012 J. Wang          add sigio option
!  Feb 04  2013 W. Yang          modified for the slg version.
!
!
! !interface:
!
      module gfs_dynamics_initialize_slg_mod
!
!!uses:
!
      use gfs_dynamics_getcf_mod
      use gfs_dyn_machine, only : kind_io4
      use nemsio_module , only : nemsio_init
!
      use gfs_dyn_write_state, only : buff_mult_pieceg
      use gfs_dyn_layout1, only : ipt_lats_node_a, lats_node_a_max,lon_dim_a
      use gfs_dyn_resol_def, only : adiabatic, thermodyn_id, sfcpress_id
      use namelist_dynamics_def, only : fhrot,fhini,num_reduce,nemsio_in, lingg_a, redgg_a, semilag
      use gfs_dyn_tracer_config, only: gfs_dyn_tracer, tracer_config_init,gfs_dyn_tracer
      use gfs_dyn_io_header, only: z_r,z
      use gfs_dyn_coordinate_def
#ifndef IBM
      USE omp_lib
#endif

      implicit none

      contains

      subroutine gfs_dynamics_initialize_slg(gis_dyn, rc)

! this subroutine set up the internal state variables,
! allocate internal state arrays for initializing the gfs system.
!----------------------------------------------------------------
!
      implicit none
!
      type(gfs_dynamics_internal_state), pointer, intent(inout) :: gis_dyn
      integer,                                    intent(out)   :: rc

      integer 		:: ierr, npe_single_member, n1, n2

      integer 		:: j, l, n, ilat, locl, ikey, nrank_all
!!
! idea add
! idea-related changes
! Introduced arrays cvd00 and cvd00m (global mean temperature,h2o,o3,cld,
!  O and O2) similar to coef00 and coef00m and added diffusion coeffs.
      REAL(KIND=kind_mpi) ,allocatable:: cvd00m(:,:) ! temp, h2o,o3,water,o,o2
      REAL(KIND=kind_evod),allocatable:: cvd00 (:,:) ! temp, h2o,o3,water,o,o2
      REAL(KIND=kind_evod),allocatable:: visc(:),cond(:),diff(:),plyr(:)
      real p0
      integer indlsev,jbasev,i,k,kk
      integer num_parthds
      character(20) cfile, cfile2

      indlsev(n,l) = jbasev + (n-l)/2 + 1

! set up gfs internal state dimension and values for dynamics etc
!-------------------------------------------------------------------
      me     = gis_dyn%me
      if (me == 0)                                                      &
      write(0,*)'in initial,nbefore allocate lonsperlat,',allocated(gis_dyn%lonsperlat),'latg=',latg
!
      nodes  = gis_dyn%nodes
      nlunit = gis_dyn%nam_gfs_dyn%nlunit
      npe_single_member = gis_dyn%npe_single_member

      semilag = gis_dyn%SLG_FLAG
      call compns_dynamics(gis_dyn%deltim, gis_dyn%iret, gis_dyn%ntrac,	&
                           gis_dyn%nxpt,   gis_dyn%nypt, gis_dyn%jintmx,&
                           gis_dyn%jcap,                              	&
                           gis_dyn%levs,   gis_dyn%levr, 		&
                           gis_dyn%lonf,   gis_dyn%latg,          	&
                           gis_dyn%ntoz,   gis_dyn%ntcw, gis_dyn%ncld, 	&
                           gis_dyn%spectral_loop,               	&
                           me, gis_dyn%thermodyn_id,gis_dyn%sfcpress_id,&
                           gis_dyn%nam_gfs_dyn%nlunit, 			&
                           gis_dyn%nam_gfs_dyn%gfs_dyn_namelist,        &
                           gis_dyn%ndfi)                           ! jw
!
      call get_tracer_const(gis_dyn%ntrac,me,gis_dyn%nam_gfs_dyn%nlunit)
!
! met+chem tracer specification (Sarah Lu)
!
!      call tracer_config_init( gis_dyn%gfs_dyn_tracer, gis_dyn%ntrac,   &
      call tracer_config_init( gis_dyn%ntrac,   &
                               gis_dyn%ntoz, gis_dyn%ntcw,              &
                               gis_dyn%ncld,  me )
!      gfs_dyn_tracer = gis_dyn%gfs_dyn_tracer
      if( me == 0) then
       write(0,*)'LU_TRC, exit tracer_config_init in dyn'
       write(0,*)'LU_TRC, ntrac=     ',gfs_dyn_tracer%ntrac,ntrac
       write(0,*)'LU_TRC, ntrac_met =',gfs_dyn_tracer%ntrac_met
       write(0,*)'LU_TRC, ntrac_chem=',gfs_dyn_tracer%ntrac_chem
       do n = 1, gfs_dyn_tracer%ntrac
        write(0,*)'LU_TRC, tracer_vname=',gfs_dyn_tracer%vname(n, :)
       enddo
      endif
!
      ntrac   = gis_dyn%ntrac
      nxpt    = gis_dyn%nxpt
      nypt    = gis_dyn%nypt
      jintmx  = gis_dyn%jintmx
      jcap    = gis_dyn%jcap
      levs    = gis_dyn%levs
      levr    = gis_dyn%levr
      lonf    = gis_dyn%lonf
      latg    = gis_dyn%latg
      ntoz    = gis_dyn%ntoz
      ntcw    = gis_dyn%ntcw
      ncld    = gis_dyn%ncld
      thermodyn_id = gis_dyn%thermodyn_id
      sfcpress_id  = gis_dyn%sfcpress_id
      nemsio_in    = gis_dyn%nemsio_in
      lon_dim_a = lonf + 2

      if (gis_dyn%nam_gfs_dyn%total_member <= 1) then
        ens_nam=' '
      else
        write(ens_nam,'("_",i2.2)') gis_dyn%nam_gfs_dyn%member_id
      endif
!
      levh   = ntrac*levs
      latgd  = latg+ 2*jintmx 
      jcap1  = jcap+1 
      jcap2  = jcap+2 
      latg2  = latg/2 
      levm1  = levs-1 
      levp1  = levs+1 
      lonfx  = lonf + 1 + 2*nxpt+1 
      lnt    = jcap2*jcap1/2 
      lnuv   = jcap2*jcap1 
      lnt2   = 2*lnt 
      lnt22  = 2*lnt+1 
      lnte   = (jcap2/2)*((jcap2/2)+1)-1 
      lnto   = (jcap2/2)*((jcap2/2)+1)-(jcap2/2) 
      lnted  = lnte 
      lntod  = lnto 

!jw      ngrids_gg       = 2+levs*(4+ntrac)
      ngrids_gg       = 2+levs*(5+ntrac)
      gis_dyn%lnt2    = lnt2

      allocate(lat1s_a(0:jcap))
      allocate(lon_dims_a(latgd))
      allocate(lon_dims_ext(latgd))

      allocate(colrad_a(latg2))
      allocate(wgt_a(latg2))
      allocate(wgtcs_a(latg2))
      allocate(rcs2_a(latg2))
      allocate(sinlat_a(latg))
      allocate(coslat_a(latg))
      coslat_a = 0

      allocate(am(levs,levs))
      allocate(bm(levs,levs))
      allocate(cm(levs,levs))
      allocate(dm(levs,levs,jcap1))
      allocate(tor(levs))
      allocate(si(levp1))
      allocate(sik(levp1))
      allocate(sl(levs))
      allocate(slk(levs))
      allocate(del(levs))
      allocate(rdel2(levs))
      allocate(ci(levp1))
      allocate(cl(levs))
      allocate(tov(levs))
      allocate(sv(levs))
      allocate(tor_slg(LEVS))
      allocate(y_ecm(LEVS,LEVS))
      allocate(t_ecm(LEVS,LEVS))
      allocate(am_slg(LEVS,LEVS))
      allocate(bm_slg(LEVS,LEVS))
      allocate(sv_ecm(LEVS))
      allocate(sv_slg(LEVS))
      allocate(D_slg_m(levs,levs,jcap1))

      allocate(ak5(levp1))
      allocate(bk5(levp1))
      allocate(ck5(levp1)) 
      allocate(thref(levp1))
      allocate(ck(levs))
      allocate(dbk(levs))
      allocate(bkl(levs))
      allocate(amhyb(levs,levs))
      allocate(bmhyb(levs,levs))
      allocate(smhyb(levs,levs))
      allocate(hmhyb(levs,levs))
      allocate(svhyb(levs))
      allocate(tor_hyb(levs))
      allocate(d_hyb_m(levs,levs,jcap1))
      allocate(dm205_hyb(jcap1,levs,levs))

      allocate(spdmax(levs))

      if (me == 0)                                                       &
      write(0,*)'before allocate lonsperlat,',allocated(gis_dyn%lonsperlat),'latg=',latg
      allocate(gis_dyn%lonsperlat(latg))

      IF(redgg_a) THEN
          If(me == 0) PRINT *,' run with reduced quardratic grid '
          IF(lingg_a) THEN
              CALL set_lonsgg_redgg_lin (gis_dyn%lonsperlat, latg, me)
          ELSE
              CALL set_lonsgg_redgg_quad(gis_dyn%lonsperlat, latg, me)
          END IF
      ELSE ! next, for full grid.
          IF(me == 0) PRINT *,' run with full gaussian grid '

          IF(lingg_a) THEN
              CALL set_lonsgg_fullgg_lin (gis_dyn%lonsperlat, latg, me)
          ELSE
              CALL set_lonsgg_fullgg_quad(gis_dyn%lonsperlat, latg, me)
          END IF
      END IF
!
! spectral location
      P_GZ  = 0*LEVS+0*LEVH+1  !      GZE/O(LNTE/OD,2),
      P_ZEM = 0*LEVS+0*LEVH+2  !     ZEME/O(LNTE/OD,2,LEVS),
      P_DIM = 1*LEVS+0*LEVH+2  !     DIME/O(LNTE/OD,2,LEVS),
      P_TEM = 2*LEVS+0*LEVH+2  !     TEME/O(LNTE/OD,2,LEVS),
      P_QM  = 3*LEVS+0*LEVH+2  !      QME/O(LNTE/OD,2),
      P_ZE  = 3*LEVS+0*LEVH+3  !      ZEE/O(LNTE/OD,2,LEVS),
      P_DI  = 4*LEVS+0*LEVH+3  !      DIE/O(LNTE/OD,2,LEVS),
      P_TE  = 5*LEVS+0*LEVH+3  !      TEE/O(LNTE/OD,2,LEVS),
      P_Q   = 6*LEVS+0*LEVH+3  !       QE/O(LNTE/OD,2),
      P_DLAM= 6*LEVS+0*LEVH+4  !  DPDLAME/O(LNTE/OD,2),
      P_DPHI= 6*LEVS+0*LEVH+5  !  DPDPHIE/O(LNTE/OD,2),
      P_ULN = 6*LEVS+0*LEVH+6  !     ULNE/O(LNTE/OD,2,LEVS),
      P_VLN = 7*LEVS+0*LEVH+6  !     VLNE/O(LNTE/OD,2,LEVS),
      P_W   = 8*LEVS+0*LEVH+6  !       WE/O(LNTE/OD,2,LEVS),
      P_X   = 9*LEVS+0*LEVH+6  !       XE/O(LNTE/OD,2,LEVS),
      P_Y   =10*LEVS+0*LEVH+6  !       YE/O(LNTE/OD,2,LEVS),
      P_ZQ  =11*LEVS+0*LEVH+6  !      ZQE/O(LNTE/OD,2)
      P_RT  =11*LEVS+0*LEVH+7  !      RTE/O(LNTE/OD,2,LEVH),
      P_RM  =11*LEVS+1*LEVH+7  !      RME/O(LNTE/OD,2,LEVH),
      P_RQ  =11*LEVS+2*LEVH+7  !      RQE/O(LNTE/OD,2,LEVH),

      lotls  = 11*LEVS+3*LEVH+6

      g_gz   = 1
      g_uum  = g_gz  + 1        !  for grid point
      g_vvm  = g_uum + levs     !  for grid point
      g_ttm  = g_vvm + levs     !  for grid point
      g_qm   = g_ttm + levs     !  for grid point
      g_uu   = g_qm  + 1        !  for grid point
      g_vv   = g_uu  + levs     !  for grid point
      g_tt   = g_vv  + levs     !  for grid point
      g_q    = g_tt  + levs     !  for grid point

      g_rm   = g_ttm + levs     !  for grid point
      g_dpm  = g_rm  + levh     !  for grid point

      g_rq   = g_tt  + levs     !  for grid point
      g_dp   = g_rq  + levh     !  for grid point

      g_u    = g_q   + 1        !  for grid point
      g_v    = g_u   + levs     !  for grid point
      g_t    = g_v   + levs     !  for grid point
      g_rt   = g_t   + levs     !  for grid point
      g_dpn  = g_rt  + levh     !  for grid point
      g_zq   = g_dpn + levs     !  for grid point

      g_p    = g_zq  + 1   	!  for grid point 
      g_dpdt = g_p   + levs   	!  for grid point 
      g_zz   = g_dpdt+ levs     !  for grid point

      g_uup   = g_zz  + levs    !  for grid point
      g_vvp   = g_uup + levs    !  for grid point
      g_ttp   = g_vvp + levs    !  for grid point
      g_rqp   = g_ttp + levs    !  for grid point
      g_dpp   = g_rqp + levh    !  for grid point
      g_zqp   = g_dpp + levs    !  for grid point

      lotgr  = g_zqp
      lotgr6 = 4*levs+1*levh+1
!c
        lots = 6*levs+1*levh+5
        lots_slg = 8*levs+1*levh+4
        lotd = 6*levs+2*levh+0
        lota = 5*levs+1*levh+2
      lotp = 4*levs
!
        ksz     =1
        ksd     =ksz+levs
        kst     =ksd+levs
        ksr     =kst+levs
        ksdp    =ksr+levh
        ksq     =ksdp+levs
        ksplam  =ksq+1
        kspphi  =ksplam+1
        ksu     =kspphi+1
        ksv     =ksu+levs
        kzslam  =ksv+levs
        kzsphi  =kzslam+1
!
        kau     =1
        kav     =kau+levs
        kat     =kav+levs
        kar     =kat+levs
        kadp    =kar+levh
        kaps    =kadp+levs
        kazs    =kaps+1
        kap2    =kazs+1
!
      kdpphi  =1
      kzzphi  =kdpphi+levs
      kdplam  =kzzphi+levs
      kzzlam  =kdplam+levs
!
        kdtphi  =1
        kdrphi  =kdtphi+levs
        kdtlam  =kdrphi+levh
        kdrlam  =kdtlam+levs
        kdulam  =kdrlam+levh
        kdvlam  =kdulam+levs
        kduphi  =kdvlam+levs
        kdvphi  =kduphi+levs
!
! point to internal state
        gis_dyn%p_zem   = p_zem            !     zeme/o(lnte/od,2,levs),
        gis_dyn%p_dim   = p_dim            !     dime/o(lnte/od,2,levs),
        gis_dyn%p_tem   = p_tem            !     teme/o(lnte/od,2,levs),
        gis_dyn%p_rm    = p_rm             !      rme/o(lnte/od,2,levh),
        gis_dyn%p_qm    = p_qm             !      qme/o(lnte/od,2),
        gis_dyn%p_ze    = p_ze             !      zee/o(lnte/od,2,levs),
        gis_dyn%p_di    = p_di             !      die/o(lnte/od,2,levs),
        gis_dyn%p_te    = p_te             !      tee/o(lnte/od,2,levs),
        gis_dyn%p_rq    = p_rq             !      rqe/o(lnte/od,2,levh),
        gis_dyn%p_q     = p_q              !       qe/o(lnte/od,2),
        gis_dyn%p_dlam  = p_dlam           !  dpdlame/o(lnte/od,2),
        gis_dyn%p_dphi  = p_dphi           !  dpdphie/o(lnte/od,2),
        gis_dyn%p_uln   = p_uln            !     ulne/o(lnte/od,2,levs),
        gis_dyn%p_vln   = p_vln            !     vlne/o(lnte/od,2,levs),
        gis_dyn%p_w     = p_w              !       we/o(lnte/od,2,levs),
        gis_dyn%p_x     = p_x              !       xe/o(lnte/od,2,levs),
        gis_dyn%p_y     = p_y              !       ye/o(lnte/od,2,levs),
        gis_dyn%p_rt    = p_rt             !      rte/o(lnte/od,2,levh),
        gis_dyn%p_zq    = p_zq             !      zqe/o(lnte/od,2)
        gis_dyn%p_gz    = p_gz             !      gze/o(lnte/od,2),
!
      gis_dyn%g_gz    = g_gz
      gis_dyn%g_uum   = g_uum
      gis_dyn%g_vvm   = g_vvm
      gis_dyn%g_ttm   = g_ttm
      gis_dyn%g_rm    = g_rm
      gis_dyn%g_dpm   = g_dpm
      gis_dyn%g_qm    = g_qm
      gis_dyn%g_uu    = g_uu
      gis_dyn%g_vv    = g_vv
      gis_dyn%g_tt    = g_tt
      gis_dyn%g_rq    = g_rq
      gis_dyn%g_dp    = g_dp
      gis_dyn%g_q     = g_q
      gis_dyn%g_u     = g_u
      gis_dyn%g_v     = g_v
      gis_dyn%g_t     = g_t
      gis_dyn%g_rt    = g_rt
      gis_dyn%g_dpn   = g_dpn
      gis_dyn%g_zq    = g_zq
      gis_dyn%g_p     = g_p
      gis_dyn%g_dpdt  = g_dpdt
      gis_dyn%g_zz    = g_zz
      gis_dyn%g_uup   = g_uup
      gis_dyn%g_vvp   = g_vvp
      gis_dyn%g_ttp   = g_ttp
      gis_dyn%g_rqp    = g_rqp
      gis_dyn%g_dpp    = g_dpp
      gis_dyn%g_zqp    = g_zqp
!
      gis_dyn%lotls = lotls
      gis_dyn%lotgr = lotgr
      gis_dyn%lots = lots 
      gis_dyn%lots_slg = lots_slg
      gis_dyn%lotd = lotd
      gis_dyn%lota = lota
      gis_dyn%lotp = lotp
!
      allocate(gis_dyn%tee1(levs))

!     print *,' finish dimension in gfs_dynamics_initialize '

!
!@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
!!
!!      create io communicator and comp communicator
!!
      nodes_comp=nodes
!c
      call f_hpminit(me,"evod")  !jjt hpm stuff
!c
      call f_hpmstart(25,"get_ls_gftlons")
!c
!!
      call synchro
      call init_countperf(latg)
!$$$      time0=timer()
!jfe  call countperf(0,15,0.)
!
      if (me == 0) then
        print 100, jcap,levs
100   format (' smf ',i3,i3,' created august 2000 ev od ri ')
        print*,'number of threads is',num_parthds()
        print*,'number of mpi procs is',nodes
      endif
!
      gis_dyn%cons0    =    0.0d0     !constant
      gis_dyn%cons0p5  =    0.5d0     !constant
      gis_dyn%cons1200 = 1200.d0      !constant
      gis_dyn%cons3600 = 3600.d0      !constant
!
      ls_dim = (jcap1-1)/nodes+1
!!
!cxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
!
      if (me == 0)                                                       &
       write(0,*)'befpre allocate ls_nodes,',allocated(gis_dyn%ls_nodes),&
       'ls_dim=', ls_dim,'nodes=',nodes
!
      allocate (      gis_dyn%ls_node (ls_dim*3) )
      allocate (      gis_dyn%ls_nodes(ls_dim,nodes) )
      allocate (  gis_dyn%max_ls_nodes(nodes) )
!
      allocate (  gis_dyn%lats_nodes_a_fix(nodes))     ! added for mGrid
!
      allocate (  gis_dyn%lats_nodes_a(nodes) )
      allocate ( gis_dyn%global_lats_a(latg) )
!
      allocate (   gis_dyn%lats_nodes_ext(nodes) )
      allocate ( gis_dyn%global_lats_ext(latg+2*jintmx+2*nypt*(nodes-1)) )

! For creating the ESMF interface state with the GFS
! internal parallel structure.   Weiyu.
!---------------------------------------------------
      ALLOCATE(gis_dyn%TRIE_LS_SIZE      (npe_single_member))
      ALLOCATE(gis_dyn%TRIO_LS_SIZE      (npe_single_member))
      ALLOCATE(gis_dyn%TRIEO_LS_SIZE     (npe_single_member))
      ALLOCATE(gis_dyn%LS_MAX_NODE_GLOBAL(npe_single_member))
      ALLOCATE(gis_dyn%LS_NODE_GLOBAL    (LS_DIM*3, npe_single_member))

      gis_dyn%LS_NODE_GLOBAL     = 0
      gis_dyn%LS_MAX_NODE_GLOBAL = 0
      gis_dyn%TRIEO_TOTAL_SIZE   = 0

      DO i = 1, npe_single_member
          CALL GET_LS_NODE(i-1, gis_dyn%LS_NODE_GLOBAL(1, i),               &
                            gis_dyn%LS_MAX_NODE_GLOBAL(i), gis_dyn%IPRINT)
          gis_dyn%TRIE_LS_SIZE(i) = 0
          gis_dyn%TRIO_LS_SIZE(i) = 0
          DO LOCL = 1, gis_dyn%LS_MAX_NODE_GLOBAL(i)
              gis_dyn%LS_NODE_GLOBAL(LOCL+  LS_DIM, i)   = gis_dyn%TRIE_LS_SIZE(i)
              gis_dyn%LS_NODE_GLOBAL(LOCL+  2*LS_DIM, i) = gis_dyn%TRIO_LS_SIZE(i)

              L = gis_dyn%LS_NODE_GLOBAL(LOCL, i)

              gis_dyn%TRIE_LS_SIZE(i) = gis_dyn%TRIE_LS_SIZE(i) + (JCAP+3-L)/2
              gis_dyn%TRIO_LS_SIZE(i) = gis_dyn%TRIO_LS_SIZE(i) + (JCAP+2-L)/2
          END DO
          gis_dyn%TRIEO_LS_SIZE(i) = gis_dyn%TRIE_LS_SIZE(i)  + gis_dyn%TRIO_LS_SIZE(i) + 3
          gis_dyn%TRIEO_TOTAL_SIZE = gis_dyn%TRIEO_TOTAL_SIZE + gis_dyn%TRIEO_LS_SIZE(i)
      END DO


!---------------------------------------------------
!
      gis_dyn%iprint = 0
      call get_ls_node( me, gis_dyn%ls_node, ls_max_node, gis_dyn%iprint )
!
!
      len_trie_ls=0
      len_trio_ls=0
      do locl=1,ls_max_node
           gis_dyn%ls_node(locl+  ls_dim)=len_trie_ls
          gis_dyn%ls_node(locl+2*ls_dim)=len_trio_ls
         l=gis_dyn%ls_node(locl)
         len_trie_ls=len_trie_ls+(jcap+3-l)/2
         len_trio_ls=len_trio_ls+(jcap+2-l)/2
      enddo
      if (me == 0) print *,'ls_node=',gis_dyn%ls_node(1:ls_dim),'2dim=',  &
         gis_dyn%ls_node(ls_dim+1:2*ls_dim),'3dim=',  &
         gis_dyn%ls_node(2*ls_dim+1:3*ls_dim)
!c
!c
      allocate (       gis_dyn%epse  (len_trie_ls) )
      allocate (       gis_dyn%epso  (len_trio_ls) )
      allocate (       gis_dyn%epsedn(len_trie_ls) )
      allocate (       gis_dyn%epsodn(len_trio_ls) )
!c
      allocate (      gis_dyn%snnp1ev(len_trie_ls) )
      allocate (      gis_dyn%snnp1od(len_trio_ls) )
!c
      allocate (       gis_dyn%ndexev(len_trie_ls) )
      allocate (       gis_dyn%ndexod(len_trio_ls) )
!c
      allocate (      gis_dyn%plnev_a(len_trie_ls,latg2) )
      allocate (      gis_dyn%plnod_a(len_trio_ls,latg2) )
      allocate (      gis_dyn%pddev_a(len_trie_ls,latg2) )
      allocate (      gis_dyn%pddod_a(len_trio_ls,latg2) )
      allocate (      gis_dyn%plnew_a(len_trie_ls,latg2) )
      allocate (      gis_dyn%plnow_a(len_trio_ls,latg2) )
!c
      gis_dyn%maxstp=36

 
      if(me == 0) 							&
        print*,'from compns_dynamics : iret=',gis_dyn%iret		&
       ,' nsout=',nsout,' nsres=',nsres
      if(gis_dyn%iret/=0) then
        if(me == 0) print *,' incompatible namelist - aborted in main'
        call mpi_quit(13)
      endif
!!
      gis_dyn%lats_nodes_ext = 0
      call getcon_dynamics(gis_dyn%n3,gis_dyn%n4,			&
           gis_dyn%ls_node,gis_dyn%ls_nodes,gis_dyn%max_ls_nodes,       &
           gis_dyn%lats_nodes_a,gis_dyn%global_lats_a,                  &
           gis_dyn%lonsperlat,gis_dyn%lats_node_a_max,                  &
           gis_dyn%lats_nodes_ext,gis_dyn%global_lats_ext,              &
           gis_dyn%epse,gis_dyn%epso,gis_dyn%epsedn,gis_dyn%epsodn,     &
           gis_dyn%snnp1ev,gis_dyn%snnp1od,				&
           gis_dyn%ndexev,gis_dyn%ndexod,  				&
           gis_dyn%plnev_a,gis_dyn%plnod_a,				&
           gis_dyn%pddev_a,gis_dyn%pddod_a, 				&
           gis_dyn%plnew_a,gis_dyn%plnow_a,gis_dyn%colat1)
!
      gis_dyn%lats_node_a=gis_dyn%lats_nodes_a(me+1)
      gis_dyn%ipt_lats_node_a=ipt_lats_node_a
      if (me == 0)                                                       &
      write(0,*)'after getcon_dynamics,lats_node_a=',gis_dyn%lats_node_a &
       ,'ipt_lats_node_a=',gis_dyn%ipt_lats_node_a,'ngptc=',ngptc
!
      print *,'ls_nodes=',gis_dyn%ls_nodes(1:ls_dim,me+1)
!!
      gis_dyn%nblck=lonf/ngptc+1

!
! initialize coord def (xlon,xlat) and lats_nodes_a_fix
!
      gis_dyn%lats_nodes_a_fix(:) = gis_dyn%lats_node_a_max
      allocate ( gis_dyn%XLON(lonf,gis_dyn%lats_node_a) )
      allocate ( gis_dyn%XLAT(lonf,gis_dyn%lats_node_a) )

      call gfs_dyn_lonlat_para(gis_dyn%global_lats_a,        &
              gis_dyn%xlon, gis_dyn%xlat, gis_dyn%lonsperlat)

      call countperf(0,18,0.)
!!
      call countperf(0,15,0.)
      allocate (   gis_dyn%trie_ls (len_trie_ls,2,lotls) )
      allocate (   gis_dyn%trio_ls (len_trio_ls,2,lotls) )
      allocate (   gis_dyn%grid_gr (lonf,lats_node_a_max,lotgr    ) )
      allocate (   gis_dyn%grid_gr6(lonf,lats_node_a_max,lotgr6*2 ) )
      allocate (   gis_dyn%pwat    (lonf,lats_node_a) )
      allocate (   gis_dyn%ptot    (lonf,lats_node_a) )
      allocate (   gis_dyn%ptrc    (lonf,lats_node_a,ntrac) )         !glbsum
      allocate (   z(lnt2) )
      allocate (   z_r(lnt2) )
!c
      allocate (   gis_dyn%syn_ls_a(4*ls_dim,gis_dyn%lots,latg2) )
      allocate (   gis_dyn%dyn_ls_a(4*ls_dim,gis_dyn%lotd,latg2) )
!c
      allocate (   gis_dyn%syn_gr_a_1(lonfx*gis_dyn%lots,lats_dim_ext) )
      allocate (   gis_dyn%syn_gr_a_2(lonfx*gis_dyn%lots,lats_dim_ext) )
      allocate (   gis_dyn%pyn_gr_a_1(lonfx*gis_dyn%lotp,lats_dim_ext) )
      allocate (   gis_dyn%pyn_gr_a_2(lonfx*gis_dyn%lotp,lats_dim_ext) )
      allocate (   gis_dyn%dyn_gr_a_1(lonfx*gis_dyn%lotd,lats_dim_ext) )
      allocate (   gis_dyn%dyn_gr_a_2(lonfx*gis_dyn%lotd,lats_dim_ext) )
      allocate (   gis_dyn%anl_gr_a_1(lonfx*gis_dyn%lota,lats_dim_ext) )
      allocate (   gis_dyn%anl_gr_a_2(lonfx*gis_dyn%lota,lats_dim_ext) )
!!
!** allocate digital filter vars
      gis_dyn%grid_gr_dfi%z_imp=0;gis_dyn%grid_gr_dfi%ps_imp=0
      gis_dyn%grid_gr_dfi%z_imp=0;gis_dyn%grid_gr_dfi%ps_imp=0
      gis_dyn%grid_gr_dfi%temp_imp=0;gis_dyn%grid_gr_dfi%u_imp=0
      gis_dyn%grid_gr_dfi%v_imp=0;gis_dyn%grid_gr_dfi%tracer_imp=0
      gis_dyn%grid_gr_dfi%p_imp=0;gis_dyn%grid_gr_dfi%dp_imp=0
      gis_dyn%grid_gr_dfi%dpdt_imp=0
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%z_import == 1) then
        gis_dyn%grid_gr_dfi%z_imp=1
        allocate( gis_dyn%grid_gr_dfi%hs(lonf,lats_node_a_max,1))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%ps_import == 1) then
        gis_dyn%grid_gr_dfi%ps_imp=1
        allocate( gis_dyn%grid_gr_dfi%ps(lonf,lats_node_a_max,1))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%temp_import == 1) then
        gis_dyn%grid_gr_dfi%temp_imp=1
        allocate( gis_dyn%grid_gr_dfi%t(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%u_import == 1) then
        gis_dyn%grid_gr_dfi%u_imp=1
        allocate( gis_dyn%grid_gr_dfi%u(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%v_import == 1) then
        gis_dyn%grid_gr_dfi%v_imp=1
        allocate( gis_dyn%grid_gr_dfi%v(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%tracer_import == 1) then
        gis_dyn%grid_gr_dfi%tracer_imp=1
        allocate( gis_dyn%grid_gr_dfi%tracer(lonf,lats_node_a_max,ntrac*levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%p_import == 1) then
        gis_dyn%grid_gr_dfi%p_imp=1
        allocate( gis_dyn%grid_gr_dfi%p(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%dp_import == 1) then
        gis_dyn%grid_gr_dfi%dp_imp=1
        allocate( gis_dyn%grid_gr_dfi%dp(lonf,lats_node_a_max,levs))
      endif
      if(gis_dyn%ndfi>0.and. gis_dyn%esmf_sta_list%dpdt_import == 1) then
        gis_dyn%grid_gr_dfi%dpdt_imp=1
        allocate( gis_dyn%grid_gr_dfi%dpdt(lonf,lats_node_a_max,levs))
      endif
!
!## allocate output vars
      allocate(buff_mult_pieceg(lonf,lats_node_a_max,ngrids_gg))
      buff_mult_pieceg = 0.
      adiabatic = gis_dyn%adiabatic
!##

!!
      allocate (      gis_dyn%fhour_idate(1,5) )
!      write(0,*)'after allocate fhour_idate'
!
      if (me == 0) then
        print*, ' lats_dim_a=', lats_dim_a, ' lats_node_a=', lats_node_a
        print*, ' lats_dim_ext=', lats_dim_ext,                           &
                ' lats_node_ext=', lats_node_ext
      endif
!c
      gis_dyn%grid_gr  = 0.0
      gis_dyn%grid_gr6 = 0.0
      gis_dyn%ptot     = 0.0
      gis_dyn%pwat     = 0.0
      gis_dyn%ptrc     = 0.0                            !glbsum

        ilat=lats_node_a
      call countperf(1,15,0.)
!!
!c......................................................................
!c
      call countperf(0,15,0.)
      call f_hpmstop(25)
!c
!      write(0,*) 'number of latitudes ext. :',lats_node_ext,              &
!                  lats_dim_ext,lats_node_a
!!
      call countperf(1,15,0.)
!!
!      print *,' sig_ini=',gis_dyn%nam_gfs_dyn%sig_ini,			  &
!              ' sig_ini2=',gis_dyn%nam_gfs_dyn%sig_ini2 
      call countperf(0,18,0.)
      gis_dyn%pdryini = 0.0

      if (me == 0) then
        print *,' grid_ini=',gis_dyn%nam_gfs_dyn%grid_ini,'fhrot=',fhrot,    &
        'fhini=',fhini,'restart_run=',gis_dyn%restart_run
      endif

      cfile  = gis_dyn%nam_gfs_dyn%sig_ini
      cfile2 = gis_dyn%nam_gfs_dyn%sig_ini2
      n1     = 11
      n2     = 12
!
      CALL input_fields_slg(n1, n2, gis_dyn%pdryini,                       &
        gis_dyn%trie_ls, gis_dyn%trio_ls,                                  &
        gis_dyn%ls_node, gis_dyn%ls_nodes, gis_dyn%max_ls_nodes,           &
        gis_dyn%snnp1ev, gis_dyn%snnp1od,                                  &
        gis_dyn%epse, gis_dyn%epso,                                        &
        cfile, cfile2, gis_dyn%restart_run)
!
      IF(.NOT. gis_dyn%restart_run) THEN
          gis_dyn%start_step    = .true.
          gis_dyn%reset_step    = .false.
          gis_dyn%restart_step  = .false.
      ELSE
        PRINT *,'restart,filename=', TRIM(cfile), TRIM(cfile2)
        gis_dyn% start_step    = .false.
        gis_dyn% reset_step    = .false.
        gis_dyn% restart_step  = .true.
      END IF
!!
      call f_hpmstart(26,"step1")
!c
!!
      call countperf(1,18,0.)
!!
      gis_dyn%zhour        = fhour
      rc = 0
!
!
      end subroutine gfs_dynamics_initialize_slg
!


!***********************************************************************
      subroutine set_lonsgg_fullgg_lin(lonsperlatx,latx,me)
      use gfs_dyn_resol_def, only : jcap
      implicit none

!sela fix lonsperlatx for all resolutions except t878


      integer lonsperlatx(latx),latx,me
!
      integer lonsperlatx_42(64)
      integer lonsperlatx_62(64)
      integer lonsperlatx_126(190)
      integer lonsperlatx_170(256)
      integer lonsperlatx_190(288)
      integer lonsperlatx_254(384)
      integer lonsperlatx_382(384)
      integer lonsperlatx_510(576)  ! lin. grid.
      integer lonsperlatx_574(576)  ! lin. grid.
      INTEGER LONSPERLATX_878(880)

      data lonsperlatx_42/32*128,32*0/
      data lonsperlatx_62/32*128,32*0/
      data lonsperlatx_126/95*384,95*0/
      data lonsperlatx_170/128*512,128*0/
      data lonsperlatx_190/144*576,144*0/
      data lonsperlatx_254/192*768,192*0/
      data lonsperlatx_382/192*768,192*0/
      data lonsperlatx_510/288*1152,288*0/  ! lin. grid.
      data lonsperlatx_574/288*1152,288*0/  ! lin. grid.
      data lonsperlatx_878/440*1760,440*0/

      integer i

      if (jcap == 42)  lonsperlatx = lonsperlatx_42
      if (jcap == 62)  lonsperlatx = lonsperlatx_62
      if (jcap == 126) lonsperlatx = lonsperlatx_126
      if (jcap == 170) lonsperlatx = lonsperlatx_170
      if (jcap == 190) lonsperlatx = lonsperlatx_190
      if (jcap == 254) lonsperlatx = lonsperlatx_254
      if (jcap == 382) lonsperlatx = lonsperlatx_382
      if (jcap == 510) lonsperlatx = lonsperlatx_510
      if (jcap == 574) lonsperlatx = lonsperlatx_574
      if (jcap == 878) lonsperlatx = lonsperlatx_878
!
      if (jcap /= 42  .and. jcap /= 62  .and. jcap /= 126 .and.        &
          jcap /= 170 .and. jcap /= 190 .and. jcap /= 254 .and.        &
          jcap /= 382 .and. jcap /= 510 .and. jcap /= 574 .and.        &
          jcap /= 878) then
         print*,' Resolution not supported - lonsperlar/lonsperlatx    &
         &data is needed in read_lonsgg '
         stop 55
      endif
!
      if (me == 0) then
        print*,' LINGG in set_lonsgg_fullgg_lin jcap,lonsperlatx =     &
&      ',jcap,lonsperlatx
        print*,' in set_lonsgg_fullgg min=',minval(lonsperlatx)
        print*,' in set_lonsgg_fullgg max=',maxval(lonsperlatx)
      endif
      return
      end subroutine set_lonsgg_fullgg_lin




      subroutine set_lonsgg_fullgg_quad(lonsperlatx,latx,me)
      use gfs_dyn_resol_def, only : jcap
      implicit none
      integer lonsperlatx(latx),latx,me

      integer lonsperlatx_42(64)
      integer lonsperlatx_62(94)
      integer lonsperlatx_126(190)
      integer lonsperlatx_170(256)
      integer lonsperlatx_190(288)
      integer lonsperlatx_254(384)
      integer lonsperlatx_382(576)
      integer lonsperlatx_510(766)
      integer lonsperlatx_574(880)

      data lonsperlatx_42/32*128,32*0/
      data lonsperlatx_62/47*192,47*0/
      data lonsperlatx_126/95*384,95*0/
      data lonsperlatx_170/128*512,128*0/
      data lonsperlatx_190/144*576,144*0/
      data lonsperlatx_254/192*768,192*0/
      data lonsperlatx_382/288*1152,288*0/
      data lonsperlatx_510/383*1536,383*0/
      data lonsperlatx_574/440*1536,440*0/

      integer i

      if (jcap == 42)  lonsperlatx = lonsperlatx_42
      if (jcap == 62)  lonsperlatx = lonsperlatx_62
      if (jcap == 126) lonsperlatx = lonsperlatx_126
      if (jcap == 170) lonsperlatx = lonsperlatx_170
      if (jcap == 190) lonsperlatx = lonsperlatx_190
      if (jcap == 254) lonsperlatx = lonsperlatx_254
      if (jcap == 382) lonsperlatx = lonsperlatx_382
      if (jcap == 510) lonsperlatx = lonsperlatx_510
      if (jcap == 574) lonsperlatx = lonsperlatx_574
!
      if (jcap /= 42  .and. jcap /= 62  .and. jcap /= 126 .and.         &
          jcap /= 170 .and. jcap /= 190 .and. jcap /= 254 .and.         &
          jcap /= 382 .and. jcap /= 510 .and. jcap /= 574) then          
         print*,' Resolution not supported - lonsperlar/lonsperlatx     &
&         data is needed in read_lonsgg '
         stop 55
      endif
!
      if (me == 0) then
        print*,' in set_lonsgg_fullgg jcap,lonsperlatx=',jcap,          &
                 lonsperlatx
        print*,' in set_lonsgg_fullgg min,max of lonsperlatx = ',       &
                 minval(lonsperlatx), maxval(lonsperlatx)
      endif
      return
      end subroutine set_lonsgg_fullgg_quad



      subroutine set_lonsgg_redgg_lin(lonsperlatx,latx,me)
      use gfs_dyn_resol_def, only : jcap
      implicit none
      integer lonsperlatx(latx),latx,me
      integer lonsperlatx_t382(384)
      integer lonsperlatx_t510_1152_576(576)
      integer lonsperlatx_t510_1024_512(512)
      INTEGER LONSPERLATX_t574_1152_576(576)
      INTEGER LONSPERLATX_t768_1536_768(768)
      INTEGER LONSPERLATX_t878_1760_880(880)
      INTEGER LONSPERLATX_t1022_2048_1024(1024)
      INTEGER LONSPERLATX_t1148_2304_1152(1152)
      INTEGER LONSPERLATX_t1500_3072_1536(1536)
      INTEGER LONSPERLATX_t1534_3072_1536(1536)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      data lonsperlatx_t382      /                                   &
          64,   64,   64,   64,   64,   64,   72,   72,   80,   90,  &
          96,  110,  110,  120,  120,  128,  144,  144,  154,  160,  &
         168,  180,  180,  180,  192,  192,  210,  220,  220,  240,  &
         240,  240,  240,  252,  256,  280,  280,  280,  288,  288,  &
         288,  308,  308,  320,  320,  320,  330,  360,  360,  360,  &
         360,  360,  360,  384,  384,  384,  384,  420,  420,  420,  &
         440,  440,  440,  440,  440,  440,  462,  462,  462,  480,  &
         480,  480,  480,  480,  480,  504,  504,  504,  504,  512,  &
         512,  560,  560,  560,  560,  560,  560,  576,  576,  576,  &
         576,  576,  576,  576,  576,  616,  616,  616,  616,  616,  &
         616,  640,  640,  640,  640,  640,  640,  640,  640,  640,  &
         640,  660,  660,  660,  720,  720,  720,  720,  720,  720,  &
         720,  720,  720,  720,  720,  720,  720,  720,  720,  720,  &
         720,  720,  720,  720,  720,  720,  720,  720,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  192*0/

!sela  aug 21 2007  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     data lonsperlatx_t510_1152_576      /
!    .   64,  64,  64,  64,  64,  64,  64,  70,  80,  84,
!    .   88,  96, 110, 110, 120, 126, 132, 140, 144, 154,
!    .  160, 168, 176, 180, 192, 192, 198, 210, 220, 220,
!    .  224, 240, 240, 252, 252, 256, 264, 280, 280, 280,
!    .  288, 308, 308, 308, 320, 320, 330, 336, 352, 352,
!    .  352, 360, 384, 384, 384, 384, 396, 396, 420, 420,
!    .  420, 420, 440, 440, 440, 448, 448, 462, 462, 480,
!    .  480, 480, 504, 504, 504, 504, 512, 528, 528, 528,
!    .  560, 560, 560, 560, 560, 560, 576, 576, 616, 616,
!    .  616, 616, 616, 616, 616, 616, 630, 630, 640, 640,
!    .  660, 660, 660, 660, 672, 672, 704, 704, 704, 704,
!    .  704, 704, 720, 720, 720, 768, 768, 768, 768, 768,
!    .  768, 768, 768, 768, 768, 792, 792, 792, 792, 792,
!    .  840, 840, 840, 840, 840, 840, 840, 840, 840, 840,
!    .  880, 880, 880, 880, 880, 880, 880, 880, 880, 880,
!    .  896, 896, 896, 896, 924, 924, 924, 924, 924, 924,
!    .  960, 960, 960, 960, 960, 960, 960, 960, 960, 960,
!    .  990, 990, 990, 990, 990, 990, 990, 990,1008,1008,
!    . 1008,1008,1008,1008,1024,1024,1024,1024,1024,1024,
!    . 1056,1056,1056,1056,1056,1056,1056,1056,1056,1056,
!    . 1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,
!    . 1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,
!    . 1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,
!    . 1120,1152,1152,1152,1152,1152,1152,1152,1152,1152,
!    . 1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,
!    . 1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,
!    . 1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,
!    . 1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,
!    . 1152,1152,1152,1152,1152,1152,1152,1152, 288*   0/
!sela  aug 21 2007  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Added here by Moorthi from T382 Quadratic grid on 12/02/2010
!
      data lonsperlatx_t510_1152_576      /                          &
         64,  64,  64,  64,  64,  64,  64,  70,  80,  84,            & 
         88,  96, 110, 110, 120, 126, 132, 140, 144, 154,            &
        160, 168, 176, 180, 192, 192, 198, 210, 220, 220,            &
        224, 240, 240, 252, 252, 256, 264, 280, 280, 280,            &
        288, 308, 308, 308, 320, 320, 330, 336, 352, 352,            &
        352, 360, 384, 384, 384, 384, 396, 396, 420, 420,            &
        420, 420, 440, 440, 440, 448, 448, 462, 462, 480,            &
        480, 480, 504, 504, 504, 504, 512, 528, 528, 528,            &
        560, 560, 560, 560, 560, 560, 576, 576, 616, 616,            &
        616, 616, 616, 616, 616, 616, 630, 630, 640, 640,            &
        660, 660, 660, 660, 672, 672, 704, 704, 704, 704,            &
        704, 704, 720, 720, 720, 768, 768, 768, 768, 768,            &
        768, 768, 768, 768, 768, 792, 792, 792, 792, 792,            &
        840, 840, 840, 840, 840, 840, 840, 840, 840, 840,            &
        880, 880, 880, 880, 880, 880, 880, 880, 880, 880,            &
        896, 896, 896, 896, 924, 924, 924, 924, 924, 924,            &
        960, 960, 960, 960, 960, 960, 960, 960, 960, 960,            &
        990, 990, 990, 990, 990, 990, 990, 990,1008,1008,            &
       1008,1008,1008,1008,1024,1024,1024,1024,1024,1024,            &
       1056,1056,1056,1056,1056,1056,1056,1056,1056,1056,            &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,            &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,            &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,            &
       1120,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152, 288*   0/
!
!sela  aug 21 2007  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     IBM FFT RESTRICTED
!     DATA T574_LAT576/
!     DATA LONSPERLATX_t574_1152_576      /
!    .    64,   64,   64,   64,   64,   64,   64,   72,   80,   80,
!    .    90,   96,  110,  110,  120,  120,  128,  144,  144,  154,
!    .   160,  160,  180,  180,  180,  192,  192,  210,  220,  220,
!    .   240,  240,  240,  240,  252,  252,  256,  280,  280,  288,
!    .   288,  288,  308,  308,  320,  320,  320,  360,  360,  360,
!    .   360,  360,  360,  384,  384,  384,  384,  420,  420,  420,
!    .   440,  440,  440,  440,  440,  462,  462,  462,  480,  480,
!    .   480,  480,  480,  504,  504,  504,  504,  512,  512,  560,
!    .   560,  560,  560,  560,  576,  576,  576,  576,  576,  576,
!    .   616,  616,  616,  616,  640,  640,  640,  640,  640,  640,
!    .   640,  640,  660,  720,  720,  720,  720,  720,  720,  720,
!    .   720,  720,  720,  720,  720,  720,  720,  768,  768,  768,
!    .   768,  768,  768,  768,  768,  768,  768,  840,  840,  840,
!    .   840,  840,  840,  840,  840,  840,  880,  880,  880,  880,
!    .   880,  880,  880,  880,  880,  880,  880,  880,  924,  924,
!    .   924,  924,  924,  924,  924,  924,  924,  960,  960,  960,
!    .   960,  960,  960,  960,  960,  960,  960,  960,  960,  960,
!    .   960,  960,  960,  990,  990,  990, 1008, 1008, 1008, 1008,
!    .  1008, 1008, 1008, 1008, 1024, 1024, 1024, 1024, 1024, 1024,
!    .  1024, 1024, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120,
!    .  1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120,
!    .  1120, 1120, 1120, 1120, 1152, 1152, 1152, 1152, 1152, 1152,
!    .  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,
!    .  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,
!    .  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,
!    .  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,
!    .  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,
!    .  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,
!    .  1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,  288*0/
!sela  aug 21 2007  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
!     Added here by Moorthi from T382 Quadratic grid on 12/02/2010
!     IBM FFT RESTRICTED
!     DATA T574_LAT576/
      DATA LONSPERLATX_t574_1152_576      /                          &
         64,  64,  64,  64,  64,  64,  64,  70,  80,  84,            &
         88,  96, 110, 110, 120, 126, 132, 140, 144, 154,            &
        160, 168, 176, 180, 192, 192, 198, 210, 220, 220,            &
        224, 240, 240, 252, 252, 256, 264, 280, 280, 280,            &
        288, 308, 308, 308, 320, 320, 330, 336, 352, 352,            &
        352, 360, 384, 384, 384, 384, 396, 396, 420, 420,            &
        420, 420, 440, 440, 440, 448, 448, 462, 462, 480,            &
        480, 480, 504, 504, 504, 504, 512, 528, 528, 528,            &
        560, 560, 560, 560, 560, 560, 576, 576, 616, 616,            &
        616, 616, 616, 616, 616, 616, 630, 630, 640, 640,            &
        660, 660, 660, 660, 672, 672, 704, 704, 704, 704,            &
        704, 704, 720, 720, 720, 768, 768, 768, 768, 768,            &
        768, 768, 768, 768, 768, 792, 792, 792, 792, 792,            &
        840, 840, 840, 840, 840, 840, 840, 840, 840, 840,            &
        880, 880, 880, 880, 880, 880, 880, 880, 880, 880,            &
        896, 896, 896, 896, 924, 924, 924, 924, 924, 924,            &
        960, 960, 960, 960, 960, 960, 960, 960, 960, 960,            &
        990, 990, 990, 990, 990, 990, 990, 990,1008,1008,            &
       1008,1008,1008,1008,1024,1024,1024,1024,1024,1024,            &
       1056,1056,1056,1056,1056,1056,1056,1056,1056,1056,            &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,            &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,            &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,            &
       1120,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152, 288*   0/
!
!     FFT RESTRICTED
!     DATA T  768 _LAT  768/
!sela768  INTEGER LONSPERLAT(LATG)
      data LONSPERLATX_t768_1536_768  /                              &
         210,  210,  210,  210,  210,  210,  210,  210,  210,  210,  &
         210,  210,  210,  210,  210,  210,  210,  210,  210,  210,  &
         210,  210,  210,  210,  210,  210,  210,  210,  210,  210,  &
         220,  224,  240,  240,  252,  252,  256,  264,  280,  280,  &
         288,  288,  308,  308,  308,  320,  320,  330,  336,  352,  &
         352,  352,  360,  384,  384,  384,  384,  396,  396,  420,  &
         420,  420,  440,  440,  440,  440,  448,  462,  462,  480,  &
         480,  480,  504,  504,  504,  504,  512,  528,  528,  528,  &
         560,  560,  560,  560,  560,  576,  576,  576,  616,  616,  &
         616,  616,  616,  616,  616,  630,  630,  640,  640,  660,  &
         660,  660,  672,  672,  704,  704,  704,  704,  704,  704,  &
         720,  720,  720,  768,  768,  768,  768,  768,  768,  768,  &
         768,  770,  792,  792,  792,  792,  840,  840,  840,  840,  &
         840,  840,  840,  840,  880,  880,  880,  880,  880,  880,  &
         880,  880,  896,  896,  896,  924,  924,  924,  924,  924,  &
         960,  960,  960,  960,  960,  960,  960,  990,  990,  990,  &
         990,  990,  990,  990, 1008, 1008, 1008, 1024, 1024, 1024,  &
        1056, 1056, 1056, 1056, 1056, 1056, 1056, 1120, 1120, 1120,  &
        1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120,  &
        1120, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1232,  &
        1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232,  &
        1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1260, 1260,  &
        1260, 1260, 1260, 1260, 1260, 1260, 1280, 1280, 1280, 1280,  &
        1280, 1320, 1320, 1320, 1320, 1320, 1320, 1320, 1320, 1320,  &
        1320, 1320, 1320, 1344, 1344, 1344, 1344, 1344, 1344, 1344,  &
        1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386,  &
        1386, 1386, 1386, 1386, 1386, 1408, 1408, 1408, 1408, 1408,  &
        1408, 1408, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440,  &
        1440, 1440, 1440, 1440, 1440, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536,  384*0/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!          Feb 3 2009                 !!!!!!
!     IBM FFT RESTRICTED
!     DATA T  878 _LAT  880/
!     DATA LONSPERLATx_t878_1760_880       /
!    .    18,   20,   28,   36,   42,   56,   56,   64,   70,   80,
!    .    84,   96,  110,  110,  112,  120,  126,  132,  140,  144,
!    .   154,  160,  168,  176,  180,  192,  192,  198,  210,  210,
!    .   220,  224,  240,  240,  252,  252,  256,  264,  280,  280,
!    .   288,  288,  308,  308,  308,  320,  320,  330,  336,  352,
!    .   352,  352,  360,  384,  384,  384,  384,  396,  396,  420,
!    .   420,  420,  440,  440,  440,  440,  448,  462,  462,  480,
!    .   480,  480,  504,  504,  504,  504,  512,  528,  528,  528,
!    .   560,  560,  560,  560,  560,  576,  576,  576,  616,  616,
!    .   616,  616,  616,  616,  630,  630,  630,  640,  660,  660,
!    .   660,  660,  672,  672,  704,  704,  704,  704,  704,  720,
!    .   720,  720,  768,  768,  768,  768,  768,  768,  768,  768,
!    .   770,  792,  792,  792,  840,  840,  840,  840,  840,  840,
!    .   840,  840,  840,  880,  880,  880,  880,  880,  880,  880,
!    .   896,  896,  896,  924,  924,  924,  924,  924,  960,  960,
!    .   960,  960,  960,  960,  960,  990,  990,  990,  990,  990,
!    .   990, 1008, 1008, 1008, 1024, 1024, 1024, 1056, 1056, 1056,
!    .  1056, 1056, 1056, 1120, 1120, 1120, 1120, 1120, 1120, 1120,
!    .  1120, 1120, 1120, 1120, 1120, 1120, 1152, 1152, 1152, 1152,
!    .  1152, 1152, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232,
!    .  1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1260,
!    .  1260, 1260, 1260, 1260, 1260, 1280, 1280, 1280, 1280, 1280,
!    .  1320, 1320, 1320, 1320, 1320, 1320, 1320, 1320, 1320, 1344,
!    .  1344, 1344, 1344, 1344, 1344, 1386, 1386, 1386, 1386, 1386,
!    .  1386, 1386, 1386, 1386, 1386, 1386, 1408, 1408, 1408, 1408,
!    .  1408, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440,
!    .  1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,
!    .  1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,
!    .  1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1540, 1584,
!    .  1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584,
!    .  1584, 1584, 1584, 1584, 1680, 1680, 1680, 1680, 1680, 1680,
!    .  1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,
!    .  1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,
!    .  1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,
!    .  1680, 1680, 1680, 1680, 1760, 1760, 1760, 1760, 1760, 1760,
!    .  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,
!    .  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,
!    .  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,
!    .  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,
!    .  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,
!    .  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,
!    .  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,
!    .  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,
!    .  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,
!    .   440*0/

!
!     Added here by Moorthi from T574 Quadratic grid on 12/02/2010
!     IBM FFT RESTRICTED
!     DATA T  878 _LAT  880/
      DATA LONSPERLATx_t878_1760_880       /                         &
          18,   28,   32,   42,   48,   56,   64,   72,   80,   84,  &
          90,  110,  110,  110,  120,  126,  132,  140,  144,  154,  &
         160,  168,  176,  176,  192,  192,  198,  210,  210,  220,  &
         224,  240,  240,  252,  252,  256,  264,  280,  280,  288,  &
         288,  308,  308,  308,  320,  320,  330,  330,  352,  352,  &
         352,  360,  384,  384,  384,  384,  396,  396,  420,  420,  &
         420,  420,  440,  440,  440,  448,  462,  462,  462,  480,  &
         480,  480,  504,  504,  504,  504,  512,  528,  528,  528,  &
         560,  560,  560,  560,  560,  576,  576,  576,  616,  616,  &
         616,  616,  616,  616,  630,  630,  630,  640,  660,  660,  &
         660,  660,  672,  672,  704,  704,  704,  704,  704,  720,  &
         720,  720,  768,  768,  768,  768,  768,  768,  768,  768,  &
         770,  792,  792,  792,  792,  840,  840,  840,  840,  840,  &
         840,  840,  840,  880,  880,  880,  880,  880,  880,  880,  &
         896,  896,  896,  896,  924,  924,  924,  924,  924,  960,  &
         960,  960,  960,  960,  960,  960,  990,  990,  990,  990,  &
         990, 1008, 1008, 1008, 1008, 1024, 1024, 1024, 1056, 1056,  &
        1056, 1056, 1056, 1056, 1120, 1120, 1120, 1120, 1120, 1120,  &
        1120, 1120, 1120, 1120, 1120, 1120, 1120, 1152, 1152, 1152,  &
        1152, 1152, 1152, 1152, 1232, 1232, 1232, 1232, 1232, 1232,  &
        1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232,  &
        1232, 1260, 1260, 1260, 1260, 1260, 1260, 1260, 1280, 1280,  &
        1280, 1280, 1280, 1320, 1320, 1320, 1320, 1320, 1320, 1320,  &
        1320, 1320, 1344, 1344, 1344, 1344, 1344, 1344, 1386, 1386,  &
        1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1408,  &
        1408, 1408, 1408, 1408, 1440, 1440, 1440, 1440, 1440, 1440,  &
        1440, 1440, 1440, 1440, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1584, 1584, 1584, 1584, 1584, 1584,  &
        1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584,  &
        1584, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
         440*0/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     FFT RESTRICTED with even number of ons from ecmwf
!     DATA T1022_LAT1024/
      DATA LONSPERLATx_t1022_2048_1024      /                        &
          18,   32,   32,   40,   56,   56,   60,   60,   72,   72,  &
          90,   90,   90,   96,   96,  110,  110,  120,  128,  128,  &
         144,  144,  154,  160,  160,  180,  180,  180,  192,  192,  &
         210,  220,  220,  240,  240,  240,  240,  252,  252,  256,  &
         280,  280,  288,  288,  288,  308,  320,  320,  320,  320,  &
         360,  360,  360,  360,  360,  360,  384,  384,  384,  384,  &
         420,  420,  420,  440,  440,  440,  440,  462,  462,  462,  &
         480,  480,  480,  480,  480,  504,  504,  504,  512,  512,  &
         560,  560,  560,  560,  560,  576,  576,  576,  576,  576,  &
         576,  616,  616,  616,  640,  640,  640,  640,  640,  640,  &
         640,  660,  720,  720,  720,  720,  720,  720,  720,  720,  &
         720,  720,  720,  720,  768,  768,  768,  768,  768,  768,  &
         768,  768,  840,  840,  840,  840,  840,  840,  840,  880,  &
         880,  880,  880,  880,  880,  880,  880,  880,  880,  924,  &
         924,  924,  924,  924,  924,  960,  960,  960,  960,  960,  &
         960,  960,  960,  960,  960,  960,  990,  990, 1008, 1008,  &
        1008, 1008, 1008, 1024, 1024, 1024, 1024, 1120, 1120, 1120,  &
        1120, 1120, 1120, 1120, 1120, 1120, 1120, 1152, 1152, 1152,  &
        1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,  &
        1152, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232,  &
        1280, 1280, 1280, 1280, 1280, 1280, 1280, 1280, 1280, 1280,  &
        1280, 1280, 1280, 1280, 1280, 1280, 1320, 1320, 1320, 1386,  &
        1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386,  &
        1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440,  &
        1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848,  &
        1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848,  &
        1848, 1848, 1848, 1920, 1920, 1920, 1920, 1920, 1920, 1920,  &
        1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920,  &
        1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920,  &
        1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920,  &
        1920, 1920, 1920, 1920, 1920, 1920, 1920, 1980, 1980, 1980,  &
        1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 2016, 2016,  &
        2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016,  &
        2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016,  &
        2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016,  &
        2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2048, 2048,  512*0/

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!sela 1148 data
      DATA LONSPERLATx_t1148_2304_1152 /  &
          18,   22,   30,   40,   44,   56,   60,   66,   72,   80,  &
          88,   96,  110,  110,  112,  120,  126,  132,  140,  154,  &
         154,  160,  168,  176,  180,  192,  192,  198,  210,  220,  &
         220,  224,  240,  240,  252,  252,  256,  264,  280,  280,  &
         288,  308,  308,  308,  308,  320,  330,  330,  336,  352,  &
         352,  360,  360,  384,  384,  384,  396,  396,  420,  420,  &
         420,  420,  440,  440,  440,  448,  448,  462,  462,  480,  &
         480,  480,  504,  504,  504,  504,  512,  528,  528,  560,  &
         560,  560,  560,  560,  576,  576,  576,  616,  616,  616,  &
         616,  616,  616,  616,  630,  630,  640,  660,  660,  660,  &
         660,  672,  672,  704,  704,  704,  704,  704,  720,  720,  &
         720,  768,  768,  768,  768,  768,  768,  768,  770,  792,  &
         792,  792,  840,  840,  840,  840,  840,  840,  840,  840,  &
         840,  880,  880,  880,  880,  880,  880,  896,  896,  896,  &
         924,  924,  924,  924,  924,  960,  960,  960,  960,  960,  &
         960,  990,  990,  990,  990,  990, 1008, 1008, 1008, 1024,  &
        1024, 1024, 1056, 1056, 1056, 1056, 1056, 1056, 1120, 1120,  &
        1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1152,  &
        1152, 1152, 1152, 1152, 1152, 1232, 1232, 1232, 1232, 1232,  &
        1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232,  &
        1260, 1260, 1260, 1260, 1260, 1280, 1280, 1280, 1280, 1320,  &
        1320, 1320, 1320, 1320, 1320, 1320, 1344, 1344, 1344, 1344,  &
        1344, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1408,  &
        1408, 1408, 1408, 1440, 1440, 1440, 1440, 1440, 1440, 1440,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1792, 1792, 1792, 1792, 1792, 1792, 1792, 1792, 1848, 1848,  &
        1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848,  &
        1848, 1848, 1848, 1920, 1920, 1920, 1920, 1920, 1920, 1920,  &
        1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920,  &
        1920, 1920, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980,  &
        1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980,  &
        2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016,  &
        2016, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2048, 2048, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112,  &
        2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112,  &
        2112, 2112, 2112, 2112, 2112, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304,  576*0/

!T1500
!  This reduced grid was derived from Henry Juang's reduced grid method
!  with numreduce=4

      DATA LONSPERLATx_t1500_3072_1536/                              &
     &    18,   28,   36,   42,   56,   60,   64,   72,   80,   88,  &
     &    96,  110,  110,  112,  120,  126,  132,  140,  154,  154,  &
     &   160,  168,  176,  180,  192,  192,  198,  210,  210,  220,  &
     &   224,  240,  240,  252,  252,  256,  264,  280,  280,  288,  &
     &   288,  308,  308,  308,  320,  320,  330,  336,  352,  352,  &
     &   352,  360,  384,  384,  384,  384,  396,  396,  420,  420,  &
     &   420,  420,  440,  440,  440,  448,  462,  462,  480,  480,  &
     &   480,  504,  504,  504,  504,  512,  528,  528,  528,  560,  &
     &   560,  560,  560,  560,  576,  576,  576,  616,  616,  616,  &
     &   616,  616,  616,  630,  630,  640,  640,  660,  660,  660,  &
     &   672,  672,  704,  704,  704,  704,  704,  720,  720,  720,  &
     &   768,  768,  768,  768,  768,  768,  768,  768,  792,  792,  &
     &   792,  792,  840,  840,  840,  840,  840,  840,  840,  840,  &
     &   880,  880,  880,  880,  880,  880,  896,  896,  896,  924,  &
     &   924,  924,  924,  924,  960,  960,  960,  960,  960,  960,  &
     &   990,  990,  990,  990,  990, 1008, 1008, 1008, 1024, 1024,  &
     &  1024, 1056, 1056, 1056, 1056, 1056, 1120, 1120, 1120, 1120,  &
     &  1120, 1120, 1120, 1120, 1120, 1120, 1120, 1152, 1152, 1152,  &
     &  1152, 1152, 1152, 1232, 1232, 1232, 1232, 1232, 1232, 1232,  &
     &  1232, 1232, 1232, 1232, 1232, 1232, 1232, 1260, 1260, 1260,  &
     &  1260, 1260, 1280, 1280, 1280, 1320, 1320, 1320, 1320, 1320,  &
     &  1320, 1320, 1344, 1344, 1344, 1344, 1344, 1386, 1386, 1386,  &
     &  1386, 1386, 1386, 1386, 1408, 1408, 1408, 1408, 1440, 1440,  &
     &  1440, 1440, 1440, 1440, 1536, 1536, 1536, 1536, 1536, 1536,  &
     &  1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
     &  1536, 1536, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584,  &
     &  1584, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
     &  1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1760,  &
     &  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
     &  1760, 1760, 1760, 1760, 1760, 1792, 1792, 1792, 1792, 1792,  &
     &  1792, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848,  &
     &  1848, 1848, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920,  &
     &  1920, 1920, 1920, 1920, 1920, 1920, 1920, 1980, 1980, 1980,  &
     &  1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980,  &
     &  2016, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2048, 2048,  &
     &  2048, 2048, 2048, 2048, 2112, 2112, 2112, 2112, 2112, 2112,  &
     &  2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2240,  &
     &  2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
     &  2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
     &  2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2304, 2304,  &
     &  2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
     &  2304, 2304, 2304, 2304, 2310, 2310, 2464, 2464, 2464, 2464,  &
     &  2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,  &
     &  2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,  &
     &  2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,  &
     &  2464, 2464, 2464, 2464, 2464, 2464, 2520, 2520, 2520, 2520,  &
     &  2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520,  &
     &  2520, 2520, 2560, 2560, 2560, 2560, 2560, 2560, 2560, 2560,  &
     &  2560, 2560, 2560, 2560, 2640, 2640, 2640, 2640, 2640, 2640,  &
     &  2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
     &  2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
     &  2688, 2688, 2688, 2688, 2688, 2688, 2688, 2688, 2688, 2688,  &
     &  2688, 2688, 2688, 2688, 2688, 2688, 2772, 2772, 2772, 2772,  &
     &  2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772,  &
     &  2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772,  &
     &  2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772, 2816, 2816,  &
     &  2816, 2816, 2816, 2816, 2816, 2816, 2816, 2816, 2816, 2816,  &
     &  2816, 2816, 2816, 2816, 2816, 2816, 2816, 2880, 2880, 2880,  &
     &  2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880,  &
     &  2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880,  &
     &  2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,  &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 768*0/


!T1534
!  This reduced grid was derived from Henry Juang's reduced grid method
!  with numreduce=4

      DATA LONSPERLATx_t1534_3072_1536/                               &
     &    18,   28,   36,   44,   56,   60,   66,   72,   80,   88,   &
     &    96,  110,  110,  120,  126,  128,  140,  144,  154,  160,   &
     &   168,  168,  176,  192,  192,  198,  210,  210,  220,  224,   &
     &   240,  240,  252,  252,  256,  264,  280,  280,  280,  288,   &
     &   308,  308,  308,  320,  320,  330,  336,  352,  352,  352,   &
     &   360,  384,  384,  384,  384,  396,  420,  420,  420,  420,   &
     &   440,  440,  440,  448,  462,  462,  462,  480,  480,  504,   &
     &   504,  504,  504,  512,  528,  528,  528,  560,  560,  560,   &
     &   560,  560,  576,  576,  576,  616,  616,  616,  616,  616,   &
     &   616,  630,  630,  640,  640,  660,  660,  660,  672,  672,   &
     &   704,  704,  704,  704,  704,  720,  720,  768,  768,  768,   &
     &   768,  768,  768,  768,  768,  770,  792,  792,  792,  840,   &
     &   840,  840,  840,  840,  840,  840,  840,  880,  880,  880,   &
     &   880,  880,  880,  896,  896,  896,  924,  924,  924,  924,   &
     &   924,  960,  960,  960,  960,  960,  990,  990,  990,  990,   &
     &   990, 1008, 1008, 1008, 1024, 1024, 1024, 1056, 1056, 1056,   &
     &  1056, 1056, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120,   &
     &  1120, 1120, 1120, 1152, 1152, 1152, 1152, 1152, 1232, 1232,   &
     &  1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232,   &
     &  1232, 1232, 1260, 1260, 1260, 1260, 1260, 1280, 1280, 1280,   &
     &  1320, 1320, 1320, 1320, 1320, 1320, 1320, 1344, 1344, 1344,   &
     &  1344, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1408,   &
     &  1408, 1408, 1408, 1440, 1440, 1440, 1440, 1440, 1536, 1536,   &
     &  1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,   &
     &  1536, 1536, 1536, 1536, 1536, 1540, 1584, 1584, 1584, 1584,   &
     &  1584, 1584, 1584, 1584, 1680, 1680, 1680, 1680, 1680, 1680,   &
     &  1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,   &
     &  1680, 1680, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,   &
     &  1760, 1760, 1760, 1760, 1760, 1760, 1760, 1792, 1792, 1792,   &
     &  1792, 1792, 1792, 1848, 1848, 1848, 1848, 1848, 1848, 1848,   &
     &  1848, 1848, 1848, 1848, 1920, 1920, 1920, 1920, 1920, 1920,   &
     &  1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1980, 1980,   &
     &  1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980,   &
     &  1980, 2016, 2016, 2016, 2016, 2016, 2016, 2016, 2048, 2048,   &
     &  2048, 2048, 2048, 2048, 2048, 2112, 2112, 2112, 2112, 2112,   &
     &  2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2240, 2240,   &
     &  2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,   &
     &  2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,   &
     &  2240, 2240, 2240, 2240, 2240, 2240, 2240, 2304, 2304, 2304,   &
     &  2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,   &
     &  2304, 2304, 2310, 2464, 2464, 2464, 2464, 2464, 2464, 2464,   &
     &  2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,   &
     &  2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,   &
     &  2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,   &
     &  2464, 2464, 2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520,   &
     &  2520, 2520, 2520, 2520, 2520, 2520, 2520, 2560, 2560, 2560,   &
     &  2560, 2560, 2560, 2560, 2560, 2560, 2560, 2560, 2640, 2640,   &
     &  2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,   &
     &  2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,   &
     &  2640, 2688, 2688, 2688, 2688, 2688, 2688, 2688, 2688, 2688,   &
     &  2688, 2688, 2688, 2688, 2688, 2688, 2772, 2772, 2772, 2772,   &
     &  2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772,   &
     &  2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772, 2772,   &
     &  2772, 2772, 2772, 2772, 2772, 2816, 2816, 2816, 2816, 2816,   &
     &  2816, 2816, 2816, 2816, 2816, 2816, 2816, 2816, 2816, 2816,   &
     &  2816, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880,   &
     &  2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880,   &
     &  2880, 2880, 2880, 2880, 2880, 2880, 2880, 2880, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072,   &
     &  3072, 3072, 3072, 3072, 3072, 3072, 3072, 3072, 768*0/



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      integer i
      if (jcap == 382) then
         lonsperlatx = lonsperlatx_t382
      endif
      if (jcap == 510 .and. latx .eq. 576) then
         lonsperlatx = lonsperlatx_t510_1152_576
      endif
      if (jcap == 510 .and. latx .eq. 512) then
         lonsperlatx = lonsperlatx_t510_1024_512
      endif
      if (jcap == 574) then
         lonsperlatx = LONSPERLATX_t574_1152_576
      endif

      if (jcap == 766) then
         lonsperlatx = LONSPERLATX_t768_1536_768
      endif

      if (jcap == 878) then
         lonsperlatx = LONSPERLATX_t878_1760_880
      endif

      if (jcap == 1022) then
         lonsperlatx = LONSPERLATX_t1022_2048_1024  
      endif

      if (jcap == 1148) then
         lonsperlatx = LONSPERLATX_t1148_2304_1152 
      endif

      if (jcap .eq. 1500) then
         lonsperlatx = LONSPERLATX_t1500_3072_1536
      endif
      if (jcap .eq. 1534) then
         lonsperlatx = LONSPERLATX_t1534_3072_1536
      endif

      if (jcap /= 382 .and. jcap /= 510 .and. jcap /= 574 .and.    &
          jcap /= 766 .and. jcap /= 878 .and. jcap /= 1022 .and.   &
          jcap /= 1148 .and. jcap /= 1500 .and. jcap /= 1534) then  
         print*,' resolution not supported - lonsperlar/lonsperlatx &
     & in set_lonsgg_redgg_lin data is needed in read_lonsgg '
         stop 55
      endif
      if (me == 0) then
        print*,' set_lonsgg_redgg_lin jcap = ',jcap
        print*,' set_lonsgg_redgg_lin lonsperlatx = ',lonsperlatx
        print*,'min,max of lonsperlatx = ',minval(lonsperlatx),          &
                maxval(lonsperlatx)
      endif
      return
      end subroutine set_lonsgg_redgg_lin

      subroutine set_lonsgg_redgg_quad(lonsperlatx,latx,me)
      use gfs_dyn_resol_def, only : jcap
      implicit none
      integer lonsperlatx(latx),latx,me
!
      integer lonsperlatx_42(64)
      integer lonsperlatx_62(94)
      integer lonsperlatx_126(190)
      integer lonsperlatx_170(256)
      integer lonsperlatx_190(288)
      integer lonsperlatx_254(384)
      integer lonsperlatx_382(576)
      integer lonsperlatx_510(766)
      integer lonsperlatx_574(880)
      integer lonsperlatx_766(1150)
      integer lonsperlatx_878(1320)

      data lonsperlatx_42/32*128,32*0/

      data lonsperlatx_62/                                           &
        30,  30,  30,  40,  48,  56,  60,  72,  72,  80,  90,  90,   &
        96, 110, 110, 120, 120, 128, 144, 144, 144, 144, 154, 160,   &
       160, 168, 168, 180, 180, 180, 180, 180, 180, 192, 192, 192,   &
       192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 47*0/
!
      data lonsperlatx_126      /                                    &
          30,   30,   36,   48,   56,   60,   72,   72,   80,   90,  &
          96,  110,  110,  120,  120,  128,  144,  144,  154,  160,  &
         160,  180,  180,  180,  192,  192,  210,  210,  220,  220,  &
         240,  240,  240,  240,  240,  252,  256,  280,  280,  280,  &
         280,  288,  288,  288,  288,  308,  308,  308,  320,  320,  &
         320,  320,  330,  330,  360,  360,  360,  360,  360,  360,  &
         360,  360,  360,  360,  360,  360,  384,  384,  384,  384,  &
         384,  384,  384,  384,  384,  384,  384,  384,  384,  384,  &
         384,  384,  384,  384,  384,  384,  384,  384,  384,  384,  &
         384,  384,  384,  384,  384, 95*0 /
!
      data lonsperlatx_170 /                                         &
         48,  48,  48,  48,  48,  56,  60,  72,  72,  80,  90,  96,  &
        110, 110, 120, 120, 128, 144, 144, 144, 154, 160, 168, 180,  &
        180, 180, 192, 210, 210, 220, 220, 240, 240, 240, 240, 240,  &
        252, 256, 280, 280, 280, 288, 288, 288, 308, 308, 320, 320,  &
        320, 320, 330, 360, 360, 360, 360, 360, 360, 360, 384, 384,  &
        384, 384, 384, 384, 420, 420, 420, 440, 440, 440, 440, 440,  &
        440, 440, 440, 440, 462, 462, 462, 462, 462, 480, 480, 480,  &
        480, 480, 480, 480, 480, 480, 480, 480, 504, 504, 504, 504,  &
        504, 504, 504, 504, 504, 512, 512, 512, 512, 512, 512, 512,  &
        512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512, 512,  &
        512, 512, 512, 512, 512, 512, 512, 512, 128*0 /
!
      data lonsperlatx_190 /                                         &
        64,  64,  64,  64,  64,  64,  64,  70,  80,  84,             &
        88, 110, 110, 110, 120, 126, 132, 140, 144, 154,             &
       160, 168, 176, 176, 192, 192, 198, 210, 210, 220,             &
       220, 240, 240, 240, 252, 252, 256, 264, 280, 280,             &
       280, 288, 308, 308, 308, 320, 320, 320, 330, 336,             &
       352, 352, 352, 352, 360, 384, 384, 384, 384, 384,             &
       396, 396, 420, 420, 420, 420, 420, 440, 440, 440,             &
       440, 440, 448, 448, 462, 462, 462, 480, 480, 480,             &
       480, 480, 504, 504, 504, 504, 504, 504, 504, 512,             &
       512, 528, 528, 528, 528, 528, 528, 560, 560, 560,             &
       560, 560, 560, 560, 560, 560, 560, 560, 560, 560,             &
       560, 576, 576, 576, 576, 576, 576, 576, 576, 576,             &
       576, 576, 576, 576, 576, 576, 576, 576, 576, 576,             &
       576, 576, 576, 576, 576, 576, 576, 576, 576, 576,             &
       576, 576, 576, 576, 144*   0/
!
      data lonsperlatx_254      /                                    &
          64,   64,   64,   64,   64,   64,   72,   72,   80,   90,  & 
          96,  110,  110,  120,  120,  128,  144,  144,  154,  160,  &
         168,  180,  180,  180,  192,  192,  210,  220,  220,  240,  &
         240,  240,  240,  252,  256,  280,  280,  280,  288,  288,  &
         288,  308,  308,  320,  320,  320,  330,  360,  360,  360,  &
         360,  360,  360,  384,  384,  384,  384,  420,  420,  420,  &
         440,  440,  440,  440,  440,  440,  462,  462,  462,  480,  &
         480,  480,  480,  480,  480,  504,  504,  504,  504,  512,  &
         512,  560,  560,  560,  560,  560,  560,  576,  576,  576,  &
         576,  576,  576,  576,  576,  616,  616,  616,  616,  616,  &
         616,  640,  640,  640,  640,  640,  640,  640,  640,  640,  &
         640,  660,  660,  660,  720,  720,  720,  720,  720,  720,  &
         720,  720,  720,  720,  720,  720,  720,  720,  720,  720,  &
         720,  720,  720,  720,  720,  720,  720,  720,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  768,  768,  768,  768,  768,  768,  768,  768,  &
         768,  768,  192*0/
!
      data lonsperlatx_382      /                                    &
         64,  64,  64,  64,  64,  64,  64,  70,  80,  84,            &
         88,  96, 110, 110, 120, 126, 132, 140, 144, 154,            &
        160, 168, 176, 180, 192, 192, 198, 210, 220, 220,            &
        224, 240, 240, 252, 252, 256, 264, 280, 280, 280,            &
        288, 308, 308, 308, 320, 320, 330, 336, 352, 352,            &
        352, 360, 384, 384, 384, 384, 396, 396, 420, 420,            &
        420, 420, 440, 440, 440, 448, 448, 462, 462, 480,            &
        480, 480, 504, 504, 504, 504, 512, 528, 528, 528,            &
        560, 560, 560, 560, 560, 560, 576, 576, 616, 616,            &
        616, 616, 616, 616, 616, 616, 630, 630, 640, 640,            &
        660, 660, 660, 660, 672, 672, 704, 704, 704, 704,            &
        704, 704, 720, 720, 720, 768, 768, 768, 768, 768,            &
        768, 768, 768, 768, 768, 792, 792, 792, 792, 792,            &
        840, 840, 840, 840, 840, 840, 840, 840, 840, 840,            &
        880, 880, 880, 880, 880, 880, 880, 880, 880, 880,            &
        896, 896, 896, 896, 924, 924, 924, 924, 924, 924,            &
        960, 960, 960, 960, 960, 960, 960, 960, 960, 960,            &
        990, 990, 990, 990, 990, 990, 990, 990,1008,1008,            &
       1008,1008,1008,1008,1024,1024,1024,1024,1024,1024,            &
       1056,1056,1056,1056,1056,1056,1056,1056,1056,1056,            &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,            &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,            &
       1120,1120,1120,1120,1120,1120,1120,1120,1120,1120,            &
       1120,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152,1152,1152,            &
       1152,1152,1152,1152,1152,1152,1152,1152, 288*   0/
!
      data lonsperlatx_510      /                                    &
          64,   64,   64,   64,   64,   64,   72,   72,   80,   90,  &
          96,  110,  110,  120,  120,  128,  144,  144,  154,  160,  &
         168,  180,  180,  180,  192,  210,  210,  220,  220,  240,  &
         240,  240,  240,  252,  256,  280,  280,  288,  288,  288,  &
         308,  308,  320,  320,  320,  330,  360,  360,  360,  360,  &
         360,  384,  384,  384,  384,  420,  420,  440,  440,  440,  &
         440,  440,  440,  462,  462,  462,  480,  480,  480,  480,  &
         504,  504,  504,  504,  512,  512,  560,  560,  560,  560,  &
         576,  576,  576,  576,  576,  576,  616,  616,  616,  616,  &
         640,  640,  640,  640,  640,  640,  640,  660,  720,  720,  &
         720,  720,  720,  720,  720,  720,  720,  720,  720,  720,  &
         720,  768,  768,  768,  768,  768,  768,  768,  768,  840,  &
         840,  840,  840,  840,  840,  840,  840,  880,  880,  880,  &
         880,  880,  880,  880,  880,  880,  880,  924,  924,  924,  &
         924,  924,  924,  924,  960,  960,  960,  960,  960,  960,  &
         960,  960,  960,  960,  960,  990,  990,  990, 1008, 1008,  &
        1008, 1008, 1008, 1024, 1024, 1024, 1024, 1024, 1120, 1120,  &
        1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120, 1120,  &
        1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152, 1152,  &
        1152, 1152, 1152, 1152, 1152, 1152, 1232, 1232, 1232, 1232,  &
        1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1260, 1260,  &
        1260, 1260, 1260, 1260, 1260, 1260, 1260, 1260, 1260, 1260,  &
        1280, 1280, 1280, 1280, 1280, 1280, 1280, 1280, 1280, 1320,  &
        1320, 1320, 1320, 1386, 1386, 1386, 1386, 1386, 1386, 1386,  &
        1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386,  &
        1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440,  &
        1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440,  &
        1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440, 1440,  &
        1440, 1440, 1440, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536,  383*0/
!
      data lonsperlatx_574      /                                    &
          18,   28,   32,   42,   48,   56,   64,   72,   80,   84,  &
          90,  110,  110,  110,  120,  126,  132,  140,  144,  154,  &
         160,  168,  176,  176,  192,  192,  198,  210,  210,  220,  &
         224,  240,  240,  252,  252,  256,  264,  280,  280,  288,  &
         288,  308,  308,  308,  320,  320,  330,  330,  352,  352,  &
         352,  360,  384,  384,  384,  384,  396,  396,  420,  420,  &
         420,  420,  440,  440,  440,  448,  462,  462,  462,  480,  &
         480,  480,  504,  504,  504,  504,  512,  528,  528,  528,  &
         560,  560,  560,  560,  560,  576,  576,  576,  616,  616,  &
         616,  616,  616,  616,  630,  630,  630,  640,  660,  660,  &
         660,  660,  672,  672,  704,  704,  704,  704,  704,  720,  &
         720,  720,  768,  768,  768,  768,  768,  768,  768,  768,  &
         770,  792,  792,  792,  792,  840,  840,  840,  840,  840,  &
         840,  840,  840,  880,  880,  880,  880,  880,  880,  880,  &
         896,  896,  896,  896,  924,  924,  924,  924,  924,  960,  &
         960,  960,  960,  960,  960,  960,  990,  990,  990,  990,  &
         990, 1008, 1008, 1008, 1008, 1024, 1024, 1024, 1056, 1056,  &
        1056, 1056, 1056, 1056, 1120, 1120, 1120, 1120, 1120, 1120,  &
        1120, 1120, 1120, 1120, 1120, 1120, 1120, 1152, 1152, 1152,  &
        1152, 1152, 1152, 1152, 1232, 1232, 1232, 1232, 1232, 1232,  &
        1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232,  &
        1232, 1260, 1260, 1260, 1260, 1260, 1260, 1260, 1280, 1280,  &
        1280, 1280, 1280, 1320, 1320, 1320, 1320, 1320, 1320, 1320,  &
        1320, 1320, 1344, 1344, 1344, 1344, 1344, 1344, 1386, 1386,  &
        1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1386, 1408,  &
        1408, 1408, 1408, 1408, 1440, 1440, 1440, 1440, 1440, 1440,  &
        1440, 1440, 1440, 1440, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1584, 1584, 1584, 1584, 1584, 1584,  &
        1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584,  &
        1584, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
         440*0/

      DATA lonsperlatx_766      /                                    &
          18,   28,   32,   42,   48,   56,   64,   72,   80,   84,  &
          96,  110,  110,  120,  120,  126,  132,  140,  154,  154,  &
         168,  168,  176,  180,  192,  198,  210,  210,  220,  224,  &
         240,  240,  240,  252,  252,  264,  280,  280,  280,  288,  &
         308,  308,  308,  320,  320,  330,  336,  352,  352,  352,  &
         360,  384,  384,  384,  384,  396,  396,  420,  420,  420,  &
         440,  440,  440,  440,  448,  462,  462,  480,  480,  480,  &
         504,  504,  504,  504,  512,  528,  528,  560,  560,  560,  &
         560,  560,  560,  576,  576,  616,  616,  616,  616,  616,  &
         616,  630,  630,  630,  640,  660,  660,  660,  660,  672,  &
         672,  704,  704,  704,  704,  704,  720,  720,  768,  768,  &
         768,  768,  768,  768,  768,  768,  770,  792,  792,  792,  &
         840,  840,  840,  840,  840,  840,  840,  840,  880,  880,  &
         880,  880,  880,  880,  896,  896,  896,  896,  924,  924,  &
         924,  924,  960,  960,  960,  960,  960,  960,  990,  990,  &
         990,  990,  990,  990, 1008, 1008, 1008, 1024, 1024, 1056,  &
        1056, 1056, 1056, 1056, 1056, 1120, 1120, 1120, 1120, 1120,  &
        1120, 1120, 1120, 1120, 1120, 1120, 1152, 1152, 1152, 1152,  &
        1152, 1152, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232,  &
        1232, 1232, 1232, 1232, 1232, 1232, 1232, 1260, 1260, 1260,  &
        1260, 1260, 1280, 1280, 1280, 1320, 1320, 1320, 1320, 1320,  &
        1320, 1320, 1320, 1344, 1344, 1344, 1344, 1344, 1386, 1386,  &
        1386, 1386, 1386, 1386, 1386, 1386, 1408, 1408, 1408, 1408,  &
        1440, 1440, 1440, 1440, 1440, 1440, 1440, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1584, 1584, 1584, 1584,  &
        1584, 1584, 1584, 1584, 1584, 1584, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1792, 1792, 1792, 1792,  &
        1792, 1792, 1792, 1848, 1848, 1848, 1848, 1848, 1848, 1848,  &
        1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1920, 1920,  &
        1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920,  &
        1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1980, 1980,  &
        1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980,  &
        1980, 1980, 1980, 1980, 1980, 2016, 2016, 2016, 2016, 2016,  &
        2016, 2016, 2016, 2016, 2016, 2016, 2048, 2048, 2048, 2048,  &
        2048, 2048, 2048, 2048, 2048, 2048, 2112, 2112, 2112, 2112,  &
        2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112,  &
        2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304,  575*0/

      data lonsperlatx_878      /                                    &
          18,   28,   32,   42,   48,   56,   64,   72,   80,   84,  &
          96,  110,  110,  120,  120,  126,  132,  140,  154,  154,  &
         160,  168,  176,  180,  192,  192,  210,  210,  220,  220,  &
         240,  240,  240,  252,  252,  264,  280,  280,  280,  288,  &
         308,  308,  308,  320,  320,  330,  330,  336,  352,  352,  &
         360,  384,  384,  384,  384,  396,  396,  420,  420,  420,  &
         420,  440,  440,  440,  448,  462,  462,  480,  480,  480,  &
         504,  504,  504,  504,  512,  528,  528,  560,  560,  560,  &
         560,  560,  560,  576,  576,  616,  616,  616,  616,  616,  &
         616,  630,  630,  630,  640,  660,  660,  660,  672,  672,  &
         704,  704,  704,  704,  704,  704,  720,  720,  768,  768,  &
         768,  768,  768,  768,  768,  768,  792,  792,  792,  792,  &
         840,  840,  840,  840,  840,  840,  840,  840,  880,  880,  &
         880,  880,  880,  880,  896,  896,  896,  924,  924,  924,  &
         924,  924,  960,  960,  960,  960,  960,  960,  990,  990,  &
         990,  990,  990, 1008, 1008, 1008, 1024, 1024, 1056, 1056,  &
        1056, 1056, 1056, 1056, 1120, 1120, 1120, 1120, 1120, 1120,  &
        1120, 1120, 1120, 1120, 1152, 1152, 1152, 1152, 1152, 1152,  &
        1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232, 1232,  &
        1232, 1232, 1232, 1232, 1260, 1260, 1260, 1260, 1260, 1280,  &
        1280, 1280, 1280, 1320, 1320, 1320, 1320, 1320, 1320, 1320,  &
        1344, 1344, 1344, 1344, 1386, 1386, 1386, 1386, 1386, 1386,  &
        1386, 1386, 1408, 1408, 1408, 1408, 1440, 1440, 1440, 1440,  &
        1440, 1440, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536, 1536,  &
        1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584, 1584,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680,  &
        1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1680, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760, 1760,  &
        1760, 1760, 1760, 1760, 1760, 1792, 1792, 1792, 1792, 1792,  &
        1792, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848, 1848,  &
        1848, 1848, 1848, 1848, 1920, 1920, 1920, 1920, 1920, 1920,  &
        1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920, 1920,  &
        1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980, 1980,  &
        1980, 1980, 1980, 1980, 2016, 2016, 2016, 2016, 2016, 2016,  &
        2016, 2016, 2016, 2048, 2048, 2048, 2048, 2048, 2048, 2048,  &
        2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112, 2112,  &
        2112, 2112, 2112, 2112, 2112, 2112, 2112, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240, 2240,  &
        2240, 2240, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304, 2304,  &
        2304, 2304, 2310, 2464, 2464, 2464, 2464, 2464, 2464, 2464,  &
        2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,  &
        2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,  &
        2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,  &
        2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464,  &
        2464, 2464, 2464, 2464, 2464, 2464, 2464, 2464, 2520, 2520,  &
        2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520,  &
        2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520, 2520,  &
        2520, 2520, 2520, 2520, 2520, 2520, 2560, 2560, 2560, 2560,  &
        2560, 2560, 2560, 2560, 2560, 2560, 2560, 2560, 2560, 2560,  &
        2560, 2560, 2560, 2560, 2560, 2560, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640, 2640,  &
        660*0/

      integer i
      if (jcap == 42) then
         lonsperlatx = lonsperlatx_42
      endif
      if (jcap == 62) then
         lonsperlatx = lonsperlatx_62
      endif
      if (jcap == 126) then
         lonsperlatx = lonsperlatx_126
      endif
      if (jcap == 170) then
         lonsperlatx = lonsperlatx_170
      endif
      if (jcap == 190) then
         lonsperlatx = lonsperlatx_190
      endif
      if (jcap == 254) then
         lonsperlatx = lonsperlatx_254
      endif
      if (jcap == 382) then
         lonsperlatx = lonsperlatx_382
      endif
      if (jcap == 510) then
         lonsperlatx = lonsperlatx_510
      endif
      if (jcap == 574) then
         lonsperlatx = lonsperlatx_574
      endif
      if (jcap == 766) then
         lonsperlatx = lonsperlatx_766
      endif
      if (jcap == 878) then
         lonsperlatx = lonsperlatx_878
      endif
!
      if (jcap /= 42  .and. jcap /= 62  .and. jcap /= 126 .and.  &
          jcap /= 170 .and. jcap /= 190 .and. jcap /= 254 .and.  &
          jcap /= 382 .and. jcap /= 510 .and. jcap /= 574 .and.  &
          jcap /= 766 .and. jcap /= 878) then
         print*,' Resolution not supported - lonsperlar/lonsperlatx    &
&         data is needed in read_lonsgg '
         stop 55
      endif
!
      if (me == 0) then
        print*,'  in lonsgg_redgg jcap, lonsperlatx = ',jcap,          &
                  lonsperlatx
        print*,'min,max of lonsperlatx = ',minval(lonsperlatx),        &
                  maxval(lonsperlatx)
      endif
      return
      end subroutine set_lonsgg_redgg_quad

      end module gfs_dynamics_initialize_slg_mod
