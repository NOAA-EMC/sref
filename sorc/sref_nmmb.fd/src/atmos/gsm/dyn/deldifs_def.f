      module deldifs_def
      use gfs_dyn_MACHINE
      implicit none
      save
      REAL(KIND=KIND_EVOD),ALLOCATABLE :: DNE(:),DNO(:),
     . SF(:),RTRD(:),RTHK(:),BKLY(:),CKLY(:)			! hmhj
      end module deldifs_def
