      module d3d_def
      use machine
      implicit none
!
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DT3DT(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DQ3DT(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DU3DT(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: DV3DT(:,:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: upd_mf(:,:,:,:)
     &,                                   dwn_mf(:,:,:,:)
     &,                                   det_mf(:,:,:,:)
     &,                                   dkh(:,:,:,:)
     &,                                   rnp(:,:,:,:)
      REAL(KIND=KIND_RAD) ,ALLOCATABLE :: CLDCOV(:,:,:)
!
      contains
!
      subroutine d3d_init(ngptc,nblck,lonr,lats_node_r,levs,pl_coeff,
     &                    ldiag3d,lggfs3d)
      implicit none
      integer ngptc,nblck,lonr,lats_node_r,levs,pl_coeff
      logical ldiag3d, lggfs3d
!
      if (ldiag3d) then
        allocate (DT3DT(NGPTC,LEVS,6,NBLCK,lats_node_r))
        allocate (DU3DT(NGPTC,LEVS,4,NBLCK,lats_node_r))
        allocate (DV3DT(NGPTC,LEVS,4,NBLCK,lats_node_r))
 !    else
 !      allocate (DT3DT(1,1,6,NBLCK,lats_node_r))
 !      allocate (DU3DT(1,1,4,NBLCK,lats_node_r))
 !      allocate (DV3DT(1,1,4,NBLCK,lats_node_r))
      endif
      if (ldiag3d .or. lggfs3d) then
        allocate (DQ3DT(NGPTC,LEVS,5+pl_coeff,NBLCK,lats_node_r))
        allocate (CLDCOV(levs,LONR,lats_node_r))
        allocate (upd_mf(NGPTC,LEVS,NBLCK,lats_node_r))
        allocate (dwn_mf(NGPTC,LEVS,NBLCK,lats_node_r))
        allocate (det_mf(NGPTC,LEVS,NBLCK,lats_node_r))
 !    else
 !      allocate (DQ3DT(1,1,5+pl_coeff,NBLCK,lats_node_r))
 !      allocate (CLDCOV(1,1,lats_node_r))
 !      allocate (upd_mf(1,1,NBLCK,lats_node_r))
 !      allocate (dwn_mf(1,1,NBLCK,lats_node_r))
 !      allocate (det_mf(1,1,NBLCK,lats_node_r))
      endif
      if (lggfs3d) then
        allocate (dkh(NGPTC,LEVS,NBLCK,lats_node_r))
        allocate (rnp(NGPTC,LEVS,NBLCK,lats_node_r))
 !    else
 !      allocate (dkh(1,1,NBLCK,lats_node_r))
 !      allocate (rnp(1,1,NBLCK,lats_node_r))
      endif

!
      end subroutine d3d_init
!
      subroutine d3d_zero(ldiag3d,lggfs3d)
      implicit none
      logical ldiag3d, lggfs3d
      real, parameter :: zero=0.0
!
       if (ldiag3d) then
         DT3DT  = zero
         DU3DT  = zero
         DV3DT  = zero
       endif
       if (ldiag3d .or. lggfs3d) then
         DQ3DT  = zero
         CLDCOV = zero
         upd_mf = zero
         dwn_mf = zero
         det_mf = zero
       endif
       if (lggfs3d) then
         dkh    = zero
         rnp    = zero
       endif
!
      end subroutine d3d_zero
      end module d3d_def
