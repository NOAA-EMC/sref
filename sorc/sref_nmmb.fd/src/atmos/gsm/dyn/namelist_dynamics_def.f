      module namelist_dynamics_def

!
! program lot
! 06 Apr 2012:    Henry Juang add some options for NDSL
! 05 Oct 2012:    Jun Wang    add sigio_out
! 02 Apr 2014:    Jun Wang    add dfilevs
!
      use gfs_dyn_machine
      implicit none
      
      integer nsres,nsout,igen,ngptc,num_reduce,levwgt(2)
      integer dfilevs
      real(kind=kind_evod) fhrot,fhmax,fhout,fhres,fhini,fhdfi
      real(kind=kind_evod) filta,ref_temp,sl_epsln
      real(kind=kind_evod) hdif_fac,hdif_fac2,slrd0,wgtm(2)
      REAL(KIND = kind_evod) :: phigs,phigs_d
      logical lsfwd,ldfi_spect, semilag
      logical shuff_lats_a,reshuff_lats_a
      logical,target :: hybrid,gen_coord_hybrid
      logical zflxtvd,explicit,gg_tracers

      logical nemsio_in, nemsio_out, sigio_out
      logical reduced_grid, semi_implicit_temp_profile
      logical mass_dp, process_split
      logical settls_dep3ds,settls_dep3dg
      LOGICAL :: cont_eq_opt1, herm_x, herm_y, herm_z, 
     &           lin_xyz, lin_xy, time_extrap_etadot
      LOGICAL :: redgg_a, wgt_cub_lin_xyz, opt1_3d_qcubic,
     &           iter_one_no_interp, lingg_a

      logical ndslfv
! hmhj idea add
      logical lsidea

      character*20 ens_nam
!
      end module namelist_dynamics_def
