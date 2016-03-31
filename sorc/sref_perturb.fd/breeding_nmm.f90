PROGRAM breeding_nmm

  IMPLICIT NONE

  INCLUDE 'netcdf.inc'

  INTEGER :: io_form
  INTEGER, PARAMETER :: io_netcdf = 2 

  INTEGER :: mem, m, m1, m2

  CHARACTER(LEN=7) :: input_fname
  INTEGER :: rcode, istat
  INTEGER :: id_var, id
  INTEGER :: im, jm, lm
  INTEGER :: im_id, jm_id, lm_id
  INTEGER :: i,k,j
  REAL :: rmse_t, en_limit, en_scale

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: pd
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: t,q,u,v

  INTEGER :: iret, iunit

!-------------------------------------------------------------------------------

  io_form = 2

  mem = 3

  DO m=1,mem

    WRITE(input_fname,"(A,I2.2)") "fort.", m+19

    SELECT CASE ( io_form  )
    CASE ( IO_NETCDF )
      CALL check( nf_open(input_fname, NF_NOWRITE, iunit) )
      CALL check( nf_inq_dimid(iunit, 'west_east', im_id) )
      CALL check( nf_inq_dimlen(iunit, im_id, im) )
      CALL check( nf_inq_dimid(iunit, 'south_north', jm_id) )
      CALL check( nf_inq_dimlen(iunit, jm_id, jm) )
      CALL check( nf_inq_dimid(iunit, 'bottom_top', lm_id) )
      CALL check( nf_inq_dimlen(iunit, lm_id, lm) )

      IF ( m == 1 ) THEN
        print *, "im,jm,lm = ",im,jm,lm
        ALLOCATE (PD(im,jm,mem))
        ALLOCATE (T(im,jm,lm,mem))
        ALLOCATE (Q(im,jm,lm,mem))
        ALLOCATE (U(im,jm,lm,mem))
        ALLOCATE (V(im,jm,lm,mem))
      END IF

      CALL check( nf_inq_varid (iunit, "PD", id_var) )
      CALL check( nf_get_var_real(iunit, id_var, PD(:,:,m)) )

      CALL check( nf_inq_varid (iunit, "T", id_var) )
      CALL check( nf_get_var_real(iunit, id_var, T(:,:,:,m)) )

      CALL check( nf_inq_varid (iunit, "Q", id_var) )
      CALL check( nf_get_var_real(iunit, id_var, Q(:,:,:,m)) )

      CALL check( nf_inq_varid (iunit, "U", id_var) )
      CALL check( nf_get_var_real(iunit, id_var, U(:,:,:,m)) )

      CALL check( nf_inq_varid (iunit, "V", id_var) )
      CALL check( nf_get_var_real(iunit, id_var, V(:,:,:,m)) )

      CALL check( nf_close(iunit) )

    CASE DEFAULT
      WRITE(0,*)'unknown io_form ', io_form
      STOP
    END SELECT

    print *, "--------------------------------------------"
    print *, "member = ", m, " id:", iunit
    print *, "PD min, PD max = ",minval(PD(:,:,m)), maxval(PD(:,:,m))
    print *, "T min, T max = ",minval(T(:,:,:,m)), maxval(T(:,:,:,m))
    print *, "Q min, Q max = ",minval(Q(:,:,:,m)), maxval(Q(:,:,:,m))
    print *, "U min, U max = ",minval(U(:,:,:,m)), maxval(U(:,:,:,m))
    print *, "V min, V max = ",minval(V(:,:,:,m)), maxval(V(:,:,:,m))

  END DO

!
! Calculate perturbances
!

  m1=2    ! p memeber
  m2=3    ! n memeber

! calculate RMS errors
  rmse_t = 0.0
! approximately 850mb level
  DO k = 15, 15
    DO j = 1, jm
     DO i = 1, im
       rmse_t = rmse_t + (T(i,j,k,m1) - T(i,j,k,m2))**2
     END DO
   END DO
  END DO
  rmse_t = sqrt( rmse_t / float(im*jm) )
  WRITE(6,*) "rmse_t: ", rmse_t

  en_scale=1.0
! en_limit = en_scale*0.7/rmse_t
  en_limit = en_scale*0.5/rmse_t
  WRITE(6,*) "en_limit: ", en_limit

  DO k = 1, lm
     DO j = 1, jm
        DO i = 1, im
           T(i,j,k,1) = T(i,j,k,1) + en_limit*( T(i,j,k,m1) - T(i,j,k,m2) )
           Q(i,j,k,1) = Q(i,j,k,1) + en_limit*( Q(i,j,k,m1) - Q(i,j,k,m2) )
           U(i,j,k,1) = U(i,j,k,1) + en_limit*( U(i,j,k,m1) - U(i,j,k,m2) )
           V(i,j,k,1) = V(i,j,k,1) + en_limit*( V(i,j,k,m1) - V(i,j,k,m2) )
        END DO
     END DO
  END DO

  DO j = 1 , jm
     DO i = 1 , im
        PD(i,j,1) = PD(i,j,1) + en_limit*( PD(i,j,m1) - PD(i,j,m2) )
     END DO
  END DO

!
! OverWRITE control file with perturbed variables 
!
  print *, "--------------------------------------------"
  print *, "PD min, PD max = ",minval(PD(:,:,1)), maxval(PD(:,:,1))
  print *, "T min, T max = ",minval(T(:,:,:,1)), maxval(T(:,:,:,1))
  print *, "Q min, Q max = ",minval(Q(:,:,:,1)), maxval(Q(:,:,:,1))
  print *, "U min, U max = ",minval(U(:,:,:,1)), maxval(U(:,:,:,1))
  print *, "V min, V max = ",minval(V(:,:,:,1)), maxval(V(:,:,:,1))

  m=1
  WRITE(input_fname,"(A,I2.2)") "fort.", m+19

  SELECT CASE ( io_form  )
  CASE ( IO_NETCDF )

    CALL check( nf_open(input_fname, NF_WRITE, iunit) )

    CALL check( nf_inq_varid (iunit, "PD", id_var) )
    CALL check( nf_put_var_real(iunit, id_var, PD(:,:,1)) )
    CALL check( nf_inq_varid (iunit, "T", id_var) )
    CALL check( nf_put_var_real(iunit, id_var, T(:,:,:,1)) )
    CALL check( nf_inq_varid (iunit, "Q", id_var) )
    CALL check( nf_put_var_real(iunit, id_var, Q(:,:,:,1)) )
    CALL check( nf_inq_varid (iunit, "U", id_var) )
    CALL check( nf_put_var_real(iunit, id_var, U(:,:,:,1)) )
    CALL check( nf_inq_varid (iunit, "V", id_var) )
    CALL check( nf_put_var_real(iunit, id_var, V(:,:,:,1)) )

    CALL check( nf_close(iunit) )

  CASE DEFAULT
    WRITE(0,*)'unknown io_form ', io_form
    STOP
  END SELECT

  print *, "End of breeding_nmm"

END PROGRAM breeding_nmm

SUBROUTINE check(status)
  IMPLICIT NONE
  INCLUDE 'netcdf.inc'
  INTEGER, INTENT (IN) :: status
  
  IF (status /= nf_noerr) THEN 
    PRINT *, trim(nf_strerror(status))
    STOP "Stopped"
  END IF
END SUBROUTINE check  

SUBROUTINE wrf_abort
  STOP
END SUBROUTINE wrf_abort

