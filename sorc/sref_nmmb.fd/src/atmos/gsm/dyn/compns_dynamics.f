!-----------------------------------------------------------------------
      subroutine compns_dynamics (deltim,iret,
     &  ntrac,nxpt,nypt,jintmx,jcap,
     &  levs,levr,lonf,latg, ntoz,
     &  ntcw,ncld, spectral_loop, me,
     &  thermodyn_id,sfcpress_id,
     &  nlunit, gfs_dyn_namelist,ndfi)
!$$$  Subprogram Documentation Block
!
! Subprogram:  compns     Check and compute namelist frequencies
!   Prgmmr: Iredell       Org: NP23          Date: 1999-01-26
!
! Abstract: This subprogram checks global spectral model namelist
!           frequencies in hour units for validity.  If they are valid,
!           then the frequencies are computed in timestep units.
!           The following rules are applied:
!             1. the timestep must be positive;
!             2. the output frequency must be positive and
!                a multiple of the timestep to within tolerance;
!             3. the shortwave frequency must be positive and
!                a multiple of the timestep to within tolerance;
!             4. the longwave frequency must be positive and
!                a multiple of the timestep to within tolerance and
!                a multiple of the shortwave frequency;
!             5. the zeroing frequency must be positive and
!                a multiple of the timestep to within tolerance and
!                a multiple of the output frequency;
!             6. the restart frequency must be positive and
!                a multiple of the timestep to within tolerance and
!                a multiple of the longwave frequency and
!                a multiple of the zeroing frequency;
!             7. the initialization window must be non-negative and
!                a multiple of the timestep to within tolerance and
!                a multiple of the longwave frequency and
!                no longer than the restart frequency;
!             8. the cycling frequency must be non-negative and
!                a multiple of the timestep to within tolerance and
!                a multiple of the longwave frequency.
!             9. the difgital filter  must be non-negative and
!
! Program History Log:
!   1999-01-26  Iredell
!   2007-02-01  H.-M. H. Juang modify to be used for dynamics only
!   2009-11-09  Jun Wang       added ndfi 
!   2010-09-08  Jun Wang       change gfsio to nemsio
!   2011-02-11  Henry Juang    add codes to fit mass_dp and ndslfv
!   2011-02-28  Sarah Lu       add thermodyn_id and sfcpress_id
!   2012-04-06  Henry Juang    add idea for lsidea
!   2012-10-05  Jun Wang       add sigio_out
!   2012-10-18  S. Morthi      add hdif_fac, hdif_fac2, slrd0
!   2013-02-02  Henry Juang    revised reduced grid and add x number
!   2013-04-02  Jun Wang       add dfilevs for digital filer
!
! Usage:    call compns(deltim,
!    &                  fhout,fhres,
!    &                  nsout,nsres,
!    &                  iret)
!   Input Arguments:
!     tol      - real error tolerance allowed for input frequencies
!                (e.g. 0.01 for 1% of timestep maximum error allowed)
!     deltim   - real timestep in seconds
!     fhout    - real output frequency in hours
!     fhres    - real restart frequency in hours
!   Output Arguments:
!     nsout    - integer output frequency in timesteps
!     nsres    - integer restart frequency in timesteps
!     iret     - integer return code (0 if successful or
!                between 1 and 8 for which rule above was broken)
!
! Attributes:
!   Language: Fortran 90
!
!$$$

      
      use namelist_dynamics_def
      use pmgrid              , only : wgt_slg, quamon
      use gfs_dyn_mpi_def, only : liope
      implicit none

      real tol
 
      character (len=*), intent(in) :: gfs_dyn_namelist
      integer, intent(in)           :: me, nlunit
      real,intent(inout)            :: deltim
      integer,intent(out)           :: iret
      integer ntrac,nxpt,nypt,jintmx,jcap,levs,lonf,latg
      integer levr,lsoil,nmtvr,lonr,latr,ndfi,k
      integer ntoz,ntcw,ncld,spectral_loop,member_num
      integer thermodyn_id, sfcpress_id
      real    tfiltc,tem

c - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!
      namelist /nam_dyn/FHMAX,FHOUT,FHRES,FHROT,FHDFI,DELTIM,IGEN,
     & NGPTC,shuff_lats_a,reshuff_lats_a,
     & nxpt,nypt,jintmx,jcap,levs,lonf,latg,levr,
     & ntrac,ntoz,ntcw,ncld,nsout,tfiltc,
     & nemsio_in,nemsio_out,liope,ref_temp,lsidea,
     & explicit,hybrid,gen_coord_hybrid,process_split, 
     & spectral_loop,ndslfv,mass_dp,semi_implicit_temp_profile,
     & reduced_grid,hdif_fac,hdif_fac2,slrd0,
     & thermodyn_id, sfcpress_id,gg_tracers,phigs_d,
     & zflxtvd,num_reduce,sigio_out,ldfi_spect,redgg_a,lingg_a,
     & sl_epsln,settls_dep3ds,settls_dep3dg,
     & herm_x,herm_y,herm_z,lin_xyz,wgt_cub_lin_xyz,lin_xy,
     & cont_eq_opt1,quamon,wgtm,time_extrap_etadot,iter_one_no_interp,
     & opt1_3d_qcubic,dfilevs

!
      num_reduce=0	! to turn on juang reduced grid, set it 4
      fhmax=0
      fhout=0
      fhres=0
      fhrot=0
      fhdfi=0
      dfilevs=levs
      deltim=0
      igen=0
      tfiltc  = 0.85
      NGPTC=lonf

      IF(semilag) THEN
          redgg_a = .true.
          lingg_a = .true.
          sl_epsln   = 0.1
          phigs_d    = 60.0
          gg_tracers = .false.

          herm_x     = .true.
          herm_y     = .true.
          herm_z     = .true.
          lin_xyz    = .false.
          lin_xy     = .false.
          wgt_cub_lin_xyz = .false.
!         levwgt(1) = 34
!         levwgt(2) = 43
          levwgt(1) = 15
          levwgt(2) = 35
          wgtm(1)   = 0.0
          wgtm(2)   = 1.0
          quamon     = .false.
          cont_eq_opt1 = .false.
          time_extrap_etadot = .false.
          iter_one_no_interp = .false.
          opt1_3d_qcubic = .false.

          settls_dep3ds = .false.
          settls_dep3dg = .false.
      END IF
!
      hdif_fac    = 1.0
      hdif_fac2   = 1.0
      slrd0       = 0.002 ! sigma level above which Raleigh damping applied
!
      shuff_lats_a     = .true.
      reshuff_lats_a   = .false.
!
      reduced_grid     = .true.
!
      explicit         = .false.
      hybrid           = .false.
      gen_coord_hybrid = .false.                                     !hmhj
      mass_dp          = .false.                                     !hmhj
      semi_implicit_temp_profile    = .false.                        !hmhj
      process_split    = .false.                                     !hmhj
      liope            = .true.
!
      thermodyn_id     = 1
      sfcpress_id      = 1
!
      zflxtvd          = .false.
      ndslfv           = .false. ! non_iteration semi_Lagrangian finite volume
      ldfi_spect       = .false. ! digital filter on spectral coefficients
!
      nemsio_in         = .true.
      nemsio_out        = .true.
      sigio_out         = .false.
!
      lsidea            = .false. ! idea add (WAM) logical variablei 
!
      ref_temp          = 300.0
!
      nsout            = 0
      levr             = 0
      spectral_loop    = 2	! 1 for one-loop or 2 for two-loop
!
      print *,' nlunit=',nlunit,' gfs_dyn_namelist=',gfs_dyn_namelist
!$$$      read(5,nam_dyn)
      open(unit=nlunit,file=gfs_dyn_namelist)
      rewind (nlunit)
      read(nlunit,nam_dyn)

      IF(semilag) THEN
          if (herm_x .or. herm_y .or. herm_z) then
             if (me == 0)
     &  print*,'hermite interp turn off lin_xyz,wgt_cub_lin_xyz,quamon'
     &       ,',lin_xy'
             lin_xyz         = .false.
             wgt_cub_lin_xyz = .false.
             lin_xy          = .false.
             quamon          = .false.
             cont_eq_opt1    = .false.
          endif

          if (settls_dep3dg .neqv. settls_dep3ds) then
             if (me == 0)
     &  print*,'**error: settls_dep3dg should be same as settls_dep3ds'
           print*,' setting both to true '
             settls_dep3dg  = .true.
             settls_dep3ds  = .true.
          endif

!
!  Allocate and set up weighting function for two time level semi-Lagangian
!  Here levels ae from top to bottom
!
          allocate (wgt_slg(levs))
          do k = 1,levwgt(1)
              wgt_slg(k) = wgtm(1)
          enddo
          tem = (wgtm(2) - wgtm(1)) / (levwgt(2) - levwgt(1) + 1)
          do k = levwgt(1)+1,levwgt(2)
              wgt_slg(k) = wgtm(1) + tem * (k-levwgt(1))
          enddo
          do k = levwgt(2)+1,levs
              wgt_slg(k) = wgtm(2)
          enddo
!         wgt_slg(:) = 1.0
          if (me == 0) print *,' wgt_slg=',wgt_slg
      END IF ! end of semilag.
!
      if( lsidea ) then       ! idea add (WAM) case needs larger ref_temp
        if (levs > 100 .and. ref_temp < 400.0) ref_temp = 2500.0
      endif

      if (me.eq.0) write(6,nam_dyn)
      filta = tfiltc
      phigs = phigs_d * (atan(1.0)/45.0)
      if (me == 0 .and. semilag) then
        write(0,*)' phigs_d=',phigs_d,' phigs=',phigs
      endif
!
      if (levr == 0) then
        levr = levs
      endif

!
!sela - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      tol=0.01
!  Check rule 1.
      if(deltim.le.0) then
        iret=1
        return
      endif
!      write(0,*)'lver=',levr,'deltim=',deltim,'nsout=',nsout,'fhout=',
!     & fhout,'fhres=',fhres,'gen_coord_hybrid=',gen_coord_hybrid, 
!     & 'ntoz=',ntoz
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsout and check rule 2.
      if(nsout.gt.0) fhout=nsout*deltim/3600.
      nsout=nint(fhout*3600./deltim)
      if(nsout.le.0.or.abs(nsout-fhout*3600./deltim).gt.tol) then
        iret=2
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute nsres and check rule 6.
      nsres=nint(fhres*3600./deltim)
      if(nsres.le.0.or.abs(nsres-fhres*3600./deltim).gt.tol) then
        iret=6
        return
      endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  Compute ndfi and check rule 7.
      if(fhdfi.eq.0.) then
        ndfi=0
        ldfi_spect=.false.
      else
        ndfi=nint(2*fhdfi*3600./deltim)
        if(ndfi.le.0.or.abs(ndfi-2*fhdfi*3600./deltim).gt.tol.or.
     &     ndfi.gt.nsres) then
           print *,'ndfi=',ndfi,'is not equal to',2*fhdfi*3600./deltim
          iret=7
          return
        endif
      endif
!
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!!
      IF (NGPTC.GT.lonf) THEN
         NGPTC=lonf
         WRITE(0,*) "NGPTC IS TOO BIG, RESET NGPTC TO lonf",NGPTC
      ENDIF
      IF (ME.EQ.0)   WRITE(0,*) "NGPTC IS SET TO NGPTC :",NGPTC
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  All checks are successful.
!     print *,' done compns_dynamics '
      iret=0
c
      return
      end subroutine compns_dynamics
