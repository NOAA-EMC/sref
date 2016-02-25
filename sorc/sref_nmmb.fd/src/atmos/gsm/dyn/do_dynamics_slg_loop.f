      SUBROUTINE do_dynamics_slg_loop(deltim,kdt,PHOUR,
     &                 TRIE_LS,TRIO_LS,GRID_GR,grid_gr_dfi, 
     &                 GRID_GR6,
     &                 LS_NODE,LS_NODES,MAX_LS_NODES,
     &                 LATS_NODES_A,GLOBAL_LATS_A,
     &                 LONSPERLAT,
     &                 EPSE,EPSO,EPSEDN,EPSODN,
     &                 SNNP1EV,SNNP1OD,NDEXEV,NDEXOD,
     &                 PLNEV_A,PLNOD_A,PDDEV_A,PDDOD_A,
     &                 PLNEW_A,PLNOW_A,
     &                 LSOUT,COLAT1,CFHOUR1,SPS,
     &                 SYN_GR_A_1, SYN_GR_A_2,
     &                 ANL_GR_A_1, ANL_GR_A_2,
     &                 start_step,restart_step,reset_step,end_step,
     &                 dfiend_step, pdryini,nblck,ZHOUR,
     &                 pwat,ptot,ptrc,nfcstdate7)
!!
      use gfs_dyn_machine     , only : kind_evod, kind_grid
      use gfs_dyn_resol_def   , only : latg,latg2,levh,levs,levp1,
     &                                 ntrac,ncld,lonf,
     &                                 lotgr,lotgr6,lotls,
     &                                 p_di,p_dim,p_q,p_qm,p_rm,p_rq,
     &                                 p_rt,p_te,p_tem,p_uln,p_vln,
     &                                 p_w,p_x,p_y,p_ze,p_zem,p_zq,
     &                                 adiabatic, lonfx, p_dpn, p_dp, 
     &                                 p_dpm, lots, lots_slg, lota

      use gfs_dyn_layout1     , only : len_trie_ls,len_trio_ls,
     &                                 ls_dim,ls_max_node,
     &                                 me,me_l_0,nodes,lats_dim_a,
     &                                 ipt_lats_node_a,lats_node_a,
     &                                 lats_node_a_max, lats_dim_ext
      use gfs_dyn_vert_def,     only : am,bm,si,sl,sv,tov
      use gfs_dyn_date_def    , only : fhour,idate,shour,spdmax
      use namelist_dynamics_def,only : ens_nam, mass_dp,
     &                                 gen_coord_hybrid, gg_tracers,
     &                                 hybrid, igen, explicit,
     &                                 nsres, fhdfi, ldfi_spect,
     &                                 sl_epsln, nsout, filta
      use gfs_dyn_mpi_def     , only : kind_mpi,
     &                                 mc_comp,mpi_r_mpi
!      USE module_gfs_mpi_def,     ONLY: num_pes_fcst
      USE gfs_dyn_dfi_mod,        ONLY: gfs_dfi_grid_gr
      USE gfs_dyn_coordinate_def, ONLY: ak5, bk5
      USE gfs_dyn_gg_def,         ONLY: wgt_a
      USE bfilt_def,              ONLY: bfilte, bfilto
      USE gfs_dyn_tracer_config,  ONLY: glbsum
      USE do_dynamics_mod,        ONLY: do_dynamics_spectdfini_slg,
     &                                  do_dynamics_syn2gridn,
     &                                  do_dynamics_gridfilter,
     &                                  do_dynamics_spectn2c,
     &                                  do_dynamics_gridn2c,
     &                                  do_dynamics_gridn2m,
     &                                  do_dynamics_gridp2n,
     &                                  do_dynamics_gridn2anl,
     &                                  do_dynamics_spectc2n,
     &                                  do_dynamics_spectn2m

      IMPLICIT NONE
!!     
      integer lat
      INTEGER, SAVE :: ndfih, kdtdfi
      CHARACTER(16)             :: CFHOUR1
      INTEGER,INTENT(IN)        :: LONSPERLAT(LATG)
!!     
      INTEGER,INTENT(IN)                   :: nfcstdate7(7)
      REAL(KIND=KIND_EVOD),  INTENT(in)    :: deltim,PHOUR
      REAL(KIND=KIND_EVOD),INTENT(INOUT)   :: ZHOUR

      TYPE(gfs_dfi_grid_gr), INTENT(inout) :: grid_gr_dfi
      REAL(KIND=KIND_GRID) GRID_GR(lonf*lats_node_a_max,lotgr)
      REAL(KIND=KIND_GRID) GRID_GR6(lonf*lats_node_a_max,lotgr6 * 2)

      REAL(KIND = kind_evod) :: trie_ls_rqt(len_trie_ls,2,levs)
      REAL(KIND = kind_evod) :: trio_ls_rqt(len_trio_ls,2,levs)
      REAL(KIND = kind_evod) :: trie_ls_sfc(len_trie_ls,2)
      REAL(KIND = kind_evod) :: trio_ls_sfc(len_trio_ls,2)
      REAL(KIND = kind_grid) :: typdel(levs)
      REAL(KIND = kind_grid) :: SYN_GR_A_1(LONFX*LOTS, LATS_DIM_EXT)
      REAL(KIND = kind_grid) :: SYN_GR_A_2(LONFX*LOTS, LATS_DIM_EXT)
      REAL(KIND = kind_grid) :: ANL_GR_A_1(LONFX*LOTA, LATS_DIM_EXT)
      REAL(KIND = kind_grid) :: ANL_GR_A_2(LONFX*LOTA, LATS_DIM_EXT)

      REAL(KIND = kind_grid) :: ptotj(lats_node_a), pwatj(lats_node_a)
      REAL(KIND = kind_grid) :: ptotg(latg), pwatg(latg), typical_pgr
      REAL(KIND = kind_grid) :: sumwa, sumto, ptotp, pwatp, pdryg, 
     &                          pdryini, pcorr, filtb

      REAL(KIND = kind_grid) :: ptrc(lonf,lats_node_a,ntrac)                !glbsum
      REAL(KIND = kind_grid) :: ptrcj(lats_node_a,ntrac)                    !glbsum
      REAL(KIND = kind_grid) :: tmpj(lats_node_a)                           !glbsum
      REAL(KIND = kind_grid) :: ptrcg(latg),sumtrc(ntrac),ptrcp(ntrac)      !glbsum

      REAL(KIND = kind_grid) :: ptot(lonf,lats_node_a)
      REAL(KIND = kind_grid) :: pwat(lonf,lats_node_a)
      REAL(KIND = kind_grid), SAVE :: dt

      integer ifirst
      data ifirst /1/
      save ifirst
!
      real, allocatable   :: gzie_ln(:,:),gzio_ln(:,:),factor_b2t_ref(:)
      save gzie_ln,gzio_ln,factor_b2t_ref

      INTEGER NBLCK
      REAL(KIND=KIND_EVOD) TRIE_LS(LEN_TRIE_LS,2,lotls)
      REAL(KIND=KIND_EVOD) TRIO_LS(LEN_TRIO_LS,2,lotls)
!!
      integer              ls_node(ls_dim,3)
!!
!     ls_node(1,1) ... ls_node(ls_max_node,1) : values of L
!     ls_node(1,2) ... ls_node(ls_max_node,2) : values of jbasev
!     ls_node(1,3) ... ls_node(ls_max_node,3) : values of jbasod
!!
      INTEGER              LS_NODES(LS_DIM,NODES)
      INTEGER          MAX_LS_NODES(NODES)
      INTEGER               LATS_NODES_A(NODES)
      INTEGER              GLOBAL_LATS_A(LATG)
!
      REAL(KIND=KIND_EVOD)    EPSE(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)    EPSO(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD)  EPSEDN(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD)  EPSODN(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD) SNNP1EV(LEN_TRIE_LS)
      REAL(KIND=KIND_EVOD) SNNP1OD(LEN_TRIO_LS)
      INTEGER               NDEXEV(LEN_TRIE_LS)
      INTEGER               NDEXOD(LEN_TRIO_LS)
      REAL(KIND=KIND_EVOD)   PLNEV_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNOD_A(LEN_TRIO_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PDDEV_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PDDOD_A(LEN_TRIO_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNEW_A(LEN_TRIE_LS,LATG2)
      REAL(KIND=KIND_EVOD)   PLNOW_A(LEN_TRIO_LS,LATG2)
!
      INTEGER kdt, lan, IERR, I, J, K, L, LOCL, N, item, jtem, 
     &           ltem, ktem, lons_lat
      real(kind=kind_evod)  batah, colat1
      REAL(KIND=kind_mpi)  coef00m(LEVS,ntrac)! temp. ozone clwater  
      REAL(KIND=kind_evod) coef00(LEVS,ntrac) ! temp. ozone clwater  
      REAL(KIND=kind_evod), parameter :: CONS0P5 = 0.5, CONS1=1.0, 
     &                                   CONS2=2.0
      INTEGER              INDLSEV,JBASEV
      INTEGER              INDLSOD,JBASOD
      integer iprint
      include 'function2'

      LOGICAL LSOUT, SPS, wrt_g3d
      LOGICAL, SAVE :: fwd_step
      LOGICAL :: start_step, reset_step, end_step, restart_step, 
     &           dfiend_step
      LOGICAL, PARAMETER :: ladj = .false.
!
!
! timings
      real(kind=kind_evod) global_times_a(latg,nodes)
      real*8               rtc ,timer1,timer2
!
      shour = kdt * deltim
      fhour = shour/3600.
      filtb = (cons1-filta)*cons0p5
!
!----------------------------------------------------------
!**********************************************************ME<NUM_PES_FCST*
!      if (me < num_pes_fcst) then
!----------------------------------------------------------
!
! -----------------------------------
        ndfih=nint(fhdfi*3600./deltim)
        if( start_step ) then
! --------------------------------------------------------------
! if the first step, from internal state, prepare syn for nonlinear
! -------- this section is called once only ---------
! --------------------------------------------------------------
          print *,' start step from internal state (grid and spectral)'

          fwd_step = .true.
          dt   = deltim*0.5

          call get_cd_hyb(dt)

! dfi.
!         ldfi_spect = .false.
         kdtdfi=kdt+ndfih
         if( me==0 ) print *,' call spectdfini with ndfih=',ndfih
         IF(ndfih /= 0) THEN
             call do_dynamics_spectdfini_slg(-ndfih-1,ndfih,trie_ls,
     &           trio_ls)
         END IF

         start_step = .false.
! -------------------------------------------------------
        else if( reset_step ) then
! --------------------------------------------------------------
! if it is reset step to reset internal state to be import state
! -------- this section is called once only ---------
! --------------------------------------------------------------
!         print *,' reset internal values by import for all '

          fwd_step = .true.
          batah = 1.0 + sl_epsln       ! Moorthi
          dt   = deltim*0.5
          call get_cd_hyb(dt)

          if (ldfi_spect) then
            if( me==0 ) print *,' finialize spectdfini '
            call do_dynamics_spectdfini_slg
     &          (ndfih+1,ndfih,trie_ls,trio_ls)
            call do_dynamics_spectc2n(trie_ls,trio_ls)
            call do_dynamics_spectn2m(trie_ls,trio_ls)
            ldfi_spect=.false.
          endif

          reset_step = .false.
! -------------------------------------------------------
        else if( restart_step ) then
! --------------------------------------------------------------
          if(me==0)print *,'in restart step'
          fwd_step = .false.
          dt   = deltim
          call get_cd_hyb(dt)
!
! -------------------------------------------------------
        else    ! end start_step, begin not start_step
! ------------------------------------------------------
! start linear computation
! -----------------------------------------------------------
!
!******************************************************************ADIABATIC
          if( .not. adiabatic ) then    ! logical variable in gfs_dyn_resol_def
! --------------------------------
! do after physics, not from input
! -----------------------------------------------------------
! update after physics

! move data from physics to n+1
          call do_dynamics_gridp2n(grid_gr,global_lats_a,lonsperlat)
!
          call do_dynamics_gridn2anl(grid_gr,anl_gr_a_2,
     &                               global_lats_a,lonsperlat)
!
! transform values in grid to spectral
!
          call grid_to_spect(anl_gr_a_1,anl_gr_a_2,
     &       trie_ls,trio_ls,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,lonsperlat,
     &       epse,epso,plnew_a,plnow_a)

!
!
!----------------------------------------------------------
!
!     print *,' ----- do bfilter ------ '

          do k=1,levs
            item = p_rt + k - 1
            jtem = p_rq + k - 1
            ktem = p_w  + k - 1
            ltem = p_ze + k - 1
            do i=1,len_trie_ls
              trie_ls_rqt(i,1,k) = bfilte(i)*
     &                            (trie_ls(i,1,item)-trie_ls(i,1,jtem))
              trie_ls_rqt(i,2,k) = bfilte(i)*
     &                            (trie_ls(i,2,item)-trie_ls(i,2,jtem))
              trie_ls(i,1,ktem)  = trie_ls(i,1,ltem) +
     &                 bfilte(i) *(trie_ls(i,1,ktem)-trie_ls(i,1,ltem))
              trie_ls(i,2,ktem)  = trie_ls(i,2,ltem) +
     &                 bfilte(i) *(trie_ls(i,2,ktem)-trie_ls(i,2,ltem))
            enddo
!!
            do i=1,len_trio_ls
              trio_ls_rqt(i,1,k) = bfilto(i)*
     &                            (trio_ls(i,1,item)-trio_ls(i,1,jtem))
              trio_ls_rqt(i,2,k) = bfilto(i)*
     &                            (trio_ls(i,2,item)-trio_ls(i,2,jtem))
              trio_ls(i,1,ktem)  = trio_ls(i,1,ltem) +
     &                 bfilto(i) *(trio_ls(i,1,ktem)-trio_ls(i,1,ltem))
              trio_ls(i,2,ktem)  = trio_ls(i,2,ltem) +
     &                 bfilto(i) *(trio_ls(i,2,ktem)-trio_ls(i,2,ltem))
            enddo
          enddo
!!
          do k=1,levh
            item = p_rt+k-1
            jtem = p_rq+k-1
            do i=1,len_trie_ls
              trie_ls(i,1,item) =  trie_ls(i,1,jtem) +
     &                bfilte(i) * (trie_ls(i,1,item)-trie_ls(i,1,jtem))
              trie_ls(i,2,item) =  trie_ls(i,2,jtem) +
     &                bfilte(i) * (trie_ls(i,2,item)-trie_ls(i,2,jtem))
            enddo
            do i=1,len_trio_ls
              trio_ls(i,1,item) =  trio_ls(i,1,jtem) +
     &                bfilto(i) * (trio_ls(i,1,item)-trio_ls(i,1,jtem))
              trio_ls(i,2,item) =  trio_ls(i,2,jtem) +
     &                bfilto(i) * (trio_ls(i,2,item)-trio_ls(i,2,jtem))
            enddo
          enddo
!!
!----------------------------------------------------------------------
!         print *,' ----- do pdry correction ------ '

          if(hybrid)then
            typical_pgr=85.
            do k=1,levp1
              si(levs+2-k) = ak5(k)/typical_pgr + bk5(k)
            enddo
          endif

          do k=1,levs
            typdel(k) = si(k)-si(k+1)
          enddo

!----------------------------------------------------------------------
! adjust moisture changes to the total mass to conserve dry mass
!
          do lan=1,lats_node_a
            lat      = global_lats_a(ipt_lats_node_a-1+lan)
            lons_lat = lonsperlat(lat)
            ptotp    = 0.
            pwatp    = 0.
            do i=1,lons_lat
               ptotp = ptotp + ptot(i,lan)
               pwatp = pwatp + pwat(i,lan)
            enddo
            pwatj(lan) = pwatp / (2.*lons_lat)
            ptotj(lan) = ptotp / (2.*lons_lat)
          enddo
          call excha(lats_nodes_a,global_lats_a,ptotj,pwatj,ptotg,pwatg)
          sumwa = 0.
          sumto = 0.
          do lat=1,latg
            sumto = sumto + wgt_a(min(lat,latg-lat+1))*ptotg(lat)
            sumwa = sumwa + wgt_a(min(lat,latg-lat+1))*pwatg(lat)
          enddo

          pdryg = sumto - sumwa
          if(pdryini <= 0.) pdryini = pdryg

          if ( glbsum ) then                                             !glbsum
            do lan=1,lats_node_a
              lat      = global_lats_a(ipt_lats_node_a-1+lan)
              lons_lat = lonsperlat(lat)
              ptrcp(:) = 0.
              do n = 1, ntrac
                do i=1,lons_lat
                  ptrcp(n)   = ptrcp(n) + ptrc(i,lan,n)
                enddo
                ptrcj(lan,n) = ptrcp(n) / (2.*lons_lat)
              enddo
            enddo
            do n = 1, ntrac
              sumtrc(n) = 0.
              tmpj(:)   = ptrcj(:,n)
              call excha(lats_nodes_a,global_lats_a,ptotj,tmpj,
     &                                             ptotg,ptrcg)
              do lat=1,latg
               sumtrc(n)=sumtrc(n)+wgt_a(min(lat,latg-lat+1))*ptrcg(lat)
              enddo
            enddo
          endif                                                          !glbsum

          pcorr = (pdryini-pdryg)/sumto*sqrt(2.)
!
          if (me == 0) write(0,*)'pcorr pdryini pdryg ',pcorr,pdryini,
     &                            pdryg

          if (glbsum .and. me == me_l_0) then                            !glbsum
            write(70,111) kdt,fhour,idate                                !glbsum
            write(71,*)   kdt,(sumtrc(n),n=4,ntrac)                      !glbsum
          endif                                                          !glbsum
111       format ('kdt, fhour, idate=',i6,1x,f10.3,2x,4(i4,2x))          !glbsum


!*********************************************************************LADJ

          if (ladj) then
!
            trie_ls_sfc = 0.0
            trio_ls_sfc = 0.0

            do k=1,levs
              do i=1,len_trie_ls
                trie_ls_sfc(i,1) = trie_ls_sfc(i,1)
     &                           + typdel(k)*trie_ls_rqt(i,1,k)
                trie_ls_sfc(i,2) = trie_ls_sfc(i,2)
     &                           + typdel(k)*trie_ls_rqt(i,2,k)
              enddo
              do i=1,len_trio_ls
                trio_ls_sfc(i,1) = trio_ls_sfc(i,1)
     &                           + typdel(k)*trio_ls_rqt(i,1,k)
                trio_ls_sfc(i,2) = trio_ls_sfc(i,2)
     &                           + typdel(k)*trio_ls_rqt(i,2,k)
              enddo
            enddo

            do i=1,len_trie_ls
              trie_ls(i,1,p_q) = trie_ls_sfc(i,1)
              trie_ls(i,2,p_q) = trie_ls_sfc(i,2)
            enddo
            do i=1,len_trio_ls
              trio_ls(i,1,p_q) = trio_ls_sfc(i,1)
              trio_ls(i,2,p_q) = trio_ls_sfc(i,2)
            enddo

            if (me == me_l_0) then
              trie_ls(1,1,p_q) = trie_ls(1,1,p_q) + pcorr
            endif
!!
!!
            do k=1,levs
              item = p_x +k-1
              jtem = p_di+k-1
              ktem = p_y +k-1
              ltem = p_te+k-1
              do i=1,len_trie_ls
               trie_ls(i,1,item) =  trie_ls(i,1,jtem) +
     &                 bfilte(i) * (trie_ls(i,1,item)-trie_ls(i,1,jtem))
               trie_ls(i,2,item) =  trie_ls(i,2,jtem) +
     &                 bfilte(i) * (trie_ls(i,2,item)-trie_ls(i,2,jtem))
               trie_ls(i,1,ktem) =  trie_ls(i,1,ltem) +
     &                 bfilte(i) * (trie_ls(i,1,ktem)-trie_ls(i,1,ltem))
               trie_ls(i,2,ktem) =  trie_ls(i,2,ltem) +
     &                 bfilte(i) * (trie_ls(i,2,ktem)-trie_ls(i,2,ltem))
              enddo
              do i=1,len_trio_ls
               trio_ls(i,1,item) =  trio_ls(i,1,jtem) +
     &                 bfilto(i) * (trio_ls(i,1,item)-trio_ls(i,1,jtem))
               trio_ls(i,2,item) =  trio_ls(i,2,jtem) +
     &                 bfilto(i) * (trio_ls(i,2,item)-trio_ls(i,2,jtem))
               trio_ls(i,1,ktem) =  trio_ls(i,1,ltem) +
     &                 bfilto(i) * (trio_ls(i,1,ktem)-trio_ls(i,1,ltem))
               trio_ls(i,2,ktem) =  trio_ls(i,2,ltem) +
     &                 bfilto(i) * (trio_ls(i,2,ktem)-trio_ls(i,2,ltem))
              enddo
            enddo

!---------------------------------------------------------
!$OMP parallel do private(locl)
              do locl=1,ls_max_node
                call impadje_hyb(trie_ls(1,1,p_x),trie_ls(1,1,p_y),
     &                           trie_ls(1,1,p_zq),trie_ls(1,1,p_di),
     &                           trie_ls(1,1,p_te),trie_ls(1,1,p_q),
     &                           dt ,
     &                           trie_ls(1,1,p_uln),
     &                           trie_ls(1,1,p_vln),
     &                           snnp1ev,ndexev,ls_node,locl)
!!
                call impadjo_hyb(trio_ls(1,1,p_x),trio_ls(1,1,p_y),
     &                           trio_ls(1,1,p_zq),trio_ls(1,1,p_di),
     &                           trio_ls(1,1,p_te),trio_ls(1,1,p_q),
     &                           dt ,
     &                           trio_ls(1,1,p_uln),
     &                           trio_ls(1,1,p_vln),
     &                           snnp1od,ndexod,ls_node,locl)
              enddo

            call countperf(1,9,0.)

!---------------------------------------------------------

! -----------------------------------------------------------
          else  ! fin impadj,    following is with no impadj
! -----------------------------------------------------------

            if (me.eq.me_l_0) then
              trie_ls(1,1,p_zq)=trie_ls(1,1,p_zq)+pcorr
            endif
!!
            do k=1,levs
              do i=1,len_trie_ls
                trie_ls(i,1,p_zq) = trie_ls(i,1,p_zq)
     &                            + typdel(k)*trie_ls_rqt(i,1,k)
                trie_ls(i,2,p_zq) = trie_ls(i,2,p_zq)
     &                            + typdel(k)*trie_ls_rqt(i,2,k)
              enddo
              do i=1,len_trio_ls
                trio_ls(i,1,p_zq) = trio_ls(i,1,p_zq)
     &                            + typdel(k)*trio_ls_rqt(i,1,k)
                trio_ls(i,2,p_zq) = trio_ls(i,2,p_zq)
     &                            + typdel(k)*trio_ls_rqt(i,2,k)
              enddo
            enddo
!!
            do k=1,levs
              item = p_x +k-1
              jtem = p_di+k-1
              ktem = p_y +k-1
              ltem = p_te+k-1
              do i=1,len_trie_ls
               trie_ls(i,1,item) =  trie_ls(i,1,jtem) +
     &                 bfilte(i) * (trie_ls(i,1,item)-trie_ls(i,1,jtem))
               trie_ls(i,2,item) =  trie_ls(i,2,jtem) +
     &                 bfilte(i) * (trie_ls(i,2,item)-trie_ls(i,2,jtem))
               trie_ls(i,1,ktem) =  trie_ls(i,1,ltem) +
     &                 bfilte(i) * (trie_ls(i,1,ktem)-trie_ls(i,1,ltem))
               trie_ls(i,2,ktem) =  trie_ls(i,2,ltem) +
     &                 bfilte(i) * (trie_ls(i,2,ktem)-trie_ls(i,2,ltem))
         trie_ls(i,1,p_dpn+k-1)=trie_ls(i,1,p_dp+k-1)+bfilte(i)*
     &  (trie_ls(i,1,p_dpn+k-1)-trie_ls(i,1,p_dp+k-1))
         trie_ls(i,2,p_dpn+k-1)=trie_ls(i,2,p_dp+k-1)+bfilte(i)*
     &  (trie_ls(i,2,p_dpn+k-1)-trie_ls(i,2,p_dp+k-1))
              enddo
              do i=1,len_trio_ls
               trio_ls(i,1,item) =  trio_ls(i,1,jtem) +
     &                 bfilto(i) * (trio_ls(i,1,item)-trio_ls(i,1,jtem))
               trio_ls(i,2,item) =  trio_ls(i,2,jtem) +
     &                 bfilto(i) * (trio_ls(i,2,item)-trio_ls(i,2,jtem))
               trio_ls(i,1,ktem) =  trio_ls(i,1,ltem) +
     &                 bfilto(i) * (trio_ls(i,1,ktem)-trio_ls(i,1,ltem))
               trio_ls(i,2,ktem) =  trio_ls(i,2,ltem) +
     &                 bfilto(i) * (trio_ls(i,2,ktem)-trio_ls(i,2,ltem))
         trio_ls(i,1,p_dpn+k-1)=trio_ls(i,1,p_dp+k-1)+bfilto(i)*
     &  (trio_ls(i,1,p_dpn+k-1)-trio_ls(i,1,p_dp+k-1))
         trio_ls(i,2,p_dpn+k-1)=trio_ls(i,2,p_dp+k-1)+bfilto(i)*
     &  (trio_ls(i,2,p_dpn+k-1)-trio_ls(i,2,p_dp+k-1))
              enddo
            enddo
! ----------------------------------------
          endif   ! fin no impadj
!*********************************************************************LADJ
! ----------------------------------------
       endif    ! fin no adiabatic
!******************************************************************ADIABATIC
          if ( fwd_step ) then
            fwd_step = .false.
            dt       = deltim
            if( gen_coord_hybrid ) then
              if( mass_dp ) then

                call get_cd_hyb_gcdp(dt)
              else
                call get_cd_hyb_gc(dt)
              endif
            else if( hybrid ) then
              call get_cd_hyb(dt)
            else
              call get_cd_sig(am,bm,dt,tov,sv)
            endif
          endif
!
! ------------------------------------------------
        endif   ! end not start_step
! ------------------------------------------------
!      endif     ! only for fcst node
!**********************************************************ME<NUM_PES_FCST*
!--------------------------------------------
! =====================================================================
!--------------------------------------------
      IF(ldfi_spect) THEN
        if( me==0 ) print *,' do spectdfini at kdt=',kdt
        call do_dynamics_spectdfini_slg
     &      (kdt-kdtdfi,ndfih,trie_ls,trio_ls)
        if( kdt-kdtdfi == ndfih ) reset_step=.true.
      END IF
!
! =====================================================================
      IF(.not.restart_step) THEN
!
!      print *,'in slg loop,lsout=',lsout,'kdt=',kdt,
!    &   'nsres=',nsres
!--------------------------------------------
        IF (lsout.and.kdt /= 0) THEN
!--------------------------------------------
!!
          CALL f_hpmstart(32,"wrtout_dynamics")
!
          CALL countperf(0,18,0.)
!
          CALL WRTOUT_dynamics(PHOUR,FHOUR,ZHOUR,IDATE,
     &                TRIE_LS,TRIO_LS,grid_gr,
     &                SL,SI,
     &                ls_node,LS_NODES,MAX_LS_NODES,
     &                lats_nodes_a,global_lats_a,lonsperlat,nblck,
     &                COLAT1,CFHOUR1,snnp1ev,snnp1od,
     &                epsedn,epsodn,plnev_a,plnod_a,
     &                epse  ,epso  ,plnew_a,plnow_a,
     &                pdryini,'SIG.F')
!
          CALL f_hpmstop(32)
!
          CALL countperf(1,18,0.)
!
! ----------------------------------
        ENDIF ! if ls_out
! ----------------------------------
!
        IF (mod(kdt,nsres) == 0 .and. kdt /= 0) THEN
!!
          CALL wrt_restart_dynamics(TRIE_LS,TRIO_LS,grid_gr,
     &         SI,fhour,idate,
     &         igen,pdryini,
     &         ls_node,ls_nodes,max_ls_nodes,
     &         global_lats_a,lonsperlat,lats_nodes_a,
     &         epse,epso,plnew_a,plnow_a,
     &         ens_nam,kdt,nfcstdate7)

        ENDIF
!
!-- end of restart step
      ELSE
          restart_step=.false.
      ENDIF

! =====================================================================
!
! if the last step,
! --------------------------
      if( end_step ) then
        return
      else if(dfiend_step .or. reset_step) then
!
!  dfi end step, return to dfi routine
!------------------------------------
          RETURN
      end if
!
!
      kdt = kdt + 1

!*********************************************** following is the tdostep.
      SHOUR = SHOUR + deltim

!      if (me < num_pes_fcst) then
!

!Start do_tstep slg part.
!------------------------

!         batah = 0.
!         batah = 1.                   ! Commented by Moorthi 11/23/2010
          batah = 1.0 + sl_epsln       ! Moorthi

          if(ifirst == 1) then
            allocate ( factor_b2t_ref(levs), gzie_ln(len_trie_ls,2),
     &                 gzio_ln(len_trio_ls,2) )
            call get_cd_hyb_slg(deltim,batah)

            k = 0
            CALL deldifs(
     .                TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     X                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     X                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),    ! hmhj
     X                TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     X                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     X                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),    ! hmhj
     X                deltim,SL,LS_NODE,coef00,k,hybrid,                ! hmhj
     &                gen_coord_hybrid)

            ifirst=0
          endif
          global_times_a = 0.
          timer1         = rtc()
          call gloopa_hyb_slg
     &      (deltim,trie_ls,trio_ls,gzie_ln,gzio_ln,
     &       ls_node,ls_nodes,max_ls_nodes,
     &       lats_nodes_a,global_lats_a,
     &       lonsperlat,
     &       epse,epso,epsedn,epsodn,
     &       snnp1ev,snnp1od,ndexev,ndexod,
     &       plnev_a,plnod_a,pddev_a,pddod_a,plnew_a,plnow_a,
     &       global_times_a,kdt,batah,lsout)
          timer2 = rtc()

!$omp parallel do private(locl)
          do locl=1,ls_max_node
            call sicdife_hyb_slg(trie_ls(1,1,p_x  ), trie_ls(1,1,p_y ),
     x                         trie_ls(1,1,p_zq ), deltim/2.,
     x                         trie_ls(1,1,p_uln), trie_ls(1,1,p_vln),
     x                         ls_node,snnp1ev,ndexev,locl,batah)
            call sicdifo_hyb_slg(trio_ls(1,1,p_x  ), trio_ls(1,1,p_y ),
     x                         trio_ls(1,1,p_zq ), deltim/2.,
     x                         trio_ls(1,1,p_uln), trio_ls(1,1,p_vln),
     x                         ls_node,snnp1od,ndexod,locl,batah)
          enddo
          do j=1,len_trie_ls
            trie_ls(j,1,p_zq ) = trie_ls(j,1,p_zq)-gzie_ln(j,1)
            trie_ls(j,2,p_zq ) = trie_ls(j,2,p_zq)-gzie_ln(j,2)
          enddo
          do j=1,len_trio_ls
            trio_ls(j,1,p_zq ) = trio_ls(j,1,p_zq)-gzio_ln(j,1)
            trio_ls(j,2,p_zq ) = trio_ls(j,2,p_zq)-gzio_ln(j,2)
          enddo
!save n-1 values for diffusion, not really part of samilag scheme
          do j=1,len_trie_ls
            trie_ls(j,1,p_qm ) = trie_ls(j,1,p_zq)
            trie_ls(j,2,p_qm ) = trie_ls(j,2,p_zq)
          enddo
          do j=1,len_trio_ls
            trio_ls(j,1,p_qm ) = trio_ls(j,1,p_zq)
            trio_ls(j,2,p_qm ) = trio_ls(j,2,p_zq)
          enddo

          do k=1,levs
            do j=1,len_trie_ls
              trie_ls(j,1,p_tem+k-1) = trie_ls(j,1,p_y+k-1)
              trie_ls(j,2,p_tem+k-1) = trie_ls(j,2,p_y+k-1)
            enddo
          enddo
          do k=1,levs
            do j=1,len_trio_ls
              trio_ls(j,1,p_tem+k-1) = trio_ls(j,1,p_y+k-1)
              trio_ls(j,2,p_tem+k-1) = trio_ls(j,2,p_y+k-1)
            enddo
          enddo
!--------------------------------------------------------
          coef00(:,:) = 0.0
          IF ( ME .EQ. ME_L_0 ) THEN
            DO LOCL=1,LS_MAX_NODE
              l      = ls_node(locl,1)
              jbasev = ls_node(locl,2)
              IF ( L  ==  0 ) THEN
                N = 0
! 1 Corresponds to temperature,  2 corresponds to ozon, 3 to clwater
                DO K=1,LEVS
                  coef00(K,1) = TRIE_LS(INDLSEV(N,L),1,P_Y +K-1)
                ENDDO
              ENDIF
            END DO
          END IF

          coef00m = coef00
          CALL MPI_BCAST(coef00m,levs*ntrac,MPI_R_MPI,ME_L_0,MC_COMP,
     &                                                       IERR)
          coef00=coef00m
          call updown(sl,coef00(1,1))
c
          if (gg_tracers) then
!$omp parallel do shared(TRIE_LS,TRIO_LS)
!$omp+shared(deltim,SL,LS_NODE,coef00,hybrid)
            do k=1,levs
              CALL deldifs_tracers(
     .                TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     X                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     X                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),
     X                TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     X                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     X                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),
     X                deltim,SL,LS_NODE,coef00,k,hybrid,
     &                gen_coord_hybrid)
            enddo
          else
!
!$omp parallel do shared(TRIE_LS,TRIO_LS)
!$omp+shared(deltim,SL,LS_NODE,coef00,hybrid)
            do k=1,levs
              CALL deldifs(
     &                TRIE_LS(1,1,P_RT+k-1), TRIE_LS(1,1,P_W+k-1),
     &                TRIE_LS(1,1,P_QM    ), TRIE_LS(1,1,P_X+k-1),
     &                TRIE_LS(1,1,P_Y +k-1), TRIE_LS(1,1,P_TEM+k-1),
     &                TRIO_LS(1,1,P_RT+k-1), TRIO_LS(1,1,P_W+k-1),
     &                TRIO_LS(1,1,P_QM    ), TRIO_LS(1,1,P_X+k-1),
     &                TRIO_LS(1,1,P_Y +k-1), TRIO_LS(1,1,P_TEM+k-1),
     &                deltim,SL,LS_NODE,coef00,k,hybrid,
     &                gen_coord_hybrid)
            enddo
          endif
!--------------------------------------------------------
          do j=1,len_trie_ls
            trie_ls(j,1,p_q ) = trie_ls(j,1,p_zq)
            trie_ls(j,2,p_q ) = trie_ls(j,2,p_zq)
          enddo
          do j=1,len_trio_ls
            trio_ls(j,1,p_q ) = trio_ls(j,1,p_zq)
            trio_ls(j,2,p_q ) = trio_ls(j,2,p_zq)
          enddo
!         if (iprint .eq. 1) print*,' me = beg gloopb ',me
          timer1 = rtc()

!!$omp parallel do shared(trie_ls,ndexev,trio_ls,ndexod)
!!$omp+shared(sl,spdmax,deltim,ls_node)
!         do k=1,levs
!sela       call damp_speed(trie_ls(1,1,p_x+k-1), trie_ls(1,1,p_w +k-1),
!selax                  trie_ls(1,1,p_y+k-1), trie_ls(1,1,p_rt+k-1),
!selax                  ndexev,
!selax                  trio_ls(1,1,p_x+k-1), trio_ls(1,1,p_w +k-1),
!selax                  trio_ls(1,1,p_y+k-1), trio_ls(1,1,p_rt+k-1),
!selax                  ndexod,
!selax                  sl,spdmax(k),deltim,ls_node)
!         enddo

          do k=1,levs
            do j=1,len_trie_ls
              trie_ls(j,1,p_di+k-1) = trie_ls(j,1,p_x+k-1)
              trie_ls(j,2,p_di+k-1) = trie_ls(j,2,p_x+k-1)
              trie_ls(j,1,p_ze+k-1) = trie_ls(j,1,p_w+k-1)
              trie_ls(j,2,p_ze+k-1) = trie_ls(j,2,p_w+k-1)
              trie_ls(j,1,p_te+k-1) = trie_ls(j,1,p_y+k-1)
              trie_ls(j,2,p_te+k-1) = trie_ls(j,2,p_y+k-1)
            enddo
          enddo
          do k=1,levs
            do j=1,len_trio_ls
              trio_ls(j,1,p_di+k-1) = trio_ls(j,1,p_x+k-1)
              trio_ls(j,2,p_di+k-1) = trio_ls(j,2,p_x+k-1)
              trio_ls(j,1,p_ze+k-1) = trio_ls(j,1,p_w+k-1)
              trio_ls(j,2,p_ze+k-1) = trio_ls(j,2,p_w+k-1)
              trio_ls(j,1,p_te+k-1) = trio_ls(j,1,p_y+k-1)
              trio_ls(j,2,p_te+k-1) = trio_ls(j,2,p_y+k-1)
            enddo
          enddo
          if(.not. gg_tracers)then
            do k=1,levh
              do j=1,len_trie_ls
                trie_ls(j,1,p_rq+k-1) = trie_ls(j,1,p_rt+k-1)
                trie_ls(j,2,p_rq+k-1) = trie_ls(j,2,p_rt+k-1)
              enddo
            enddo
            do k=1,levh
              do j=1,len_trio_ls
                trio_ls(j,1,p_rq+k-1) = trio_ls(j,1,p_rt+k-1)
                trio_ls(j,2,p_rq+k-1) = trio_ls(j,2,p_rt+k-1)
              enddo
            enddo
          endif !  if(.not.gg_tracers)

      END SUBROUTINE do_dynamics_slg_loop
