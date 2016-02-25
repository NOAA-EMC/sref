SUBROUTINE wrf_debug( level , str )
  IMPLICIT NONE
  CHARACTER*(*) str
  INTEGER level

  WRITE(0,*)str

  RETURN
END SUBROUTINE wrf_debug

SUBROUTINE wrf_abort
  STOP 'wrf_abort'
END SUBROUTINE wrf_abort

LOGICAL FUNCTION wrf_dm_on_monitor()
  wrf_dm_on_monitor = .true.
END FUNCTION wrf_dm_on_monitor

