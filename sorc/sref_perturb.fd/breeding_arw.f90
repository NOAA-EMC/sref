PROGRAM breeding_arw

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

  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: p
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: t,q,u,v

  INTEGER :: iunit

  INTEGER :: iret

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
!     im = im - 1
!     jm = jm - 1
!     lm = lm - 1

      IF ( m == 1 ) THEN
        print *, "im,jm,lm = ",im,jm,lm
        ALLOCATE (P(im,jm,mem))
        ALLOCATE (T(im,lm,jm,mem))
        ALLOCATE (Q(im,lm,jm,mem))
        ALLOCATE (U(im+1,lm,jm,mem))
        ALLOCATE (V(im,lm,jm+1,mem))
      END IF

      CALL check( nf_inq_varid (iunit, "MU", id_var) )
      CALL check( nf_get_var_real(iunit, id_var, P(:,:,m)) )

      CALL check( nf_inq_varid (iunit, "T", id_var) )
      CALL check( nf_get_var_real(iunit, id_var, T(:,:,:,m)) )

      CALL check( nf_inq_varid (iunit, "QVAPOR", id_var) )
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
    print *, "P min, P max = ",minval(P(:,:,m)), maxval(P(:,:,m))
    do k=1,lm
    print *, "T min, T max = ",K,minval(T(:,K,:,m)), maxval(T(:,K,:,m))
    print *, "Q min, Q max = ",K,minval(Q(:,K,:,m)), maxval(Q(:,K,:,m))
    print *, "U min, U max = ",K,minval(U(:,K,:,m)), maxval(U(:,K,:,m))
    print *, "V min, V max = ",K,minval(V(:,K,:,m)), maxval(V(:,K,:,m))
    end do

  END DO

!
! Calculate perturbances
!

  m1=2    ! p memeber
  m2=3    ! n memeber

! calculate RMS errors
  rmse_t = 0.0
! approximately 850mb level
  DO j = 1, jm
    DO k = 20, 20
     DO i = 1, im
       rmse_t = rmse_t + (T(i,k,j,m1) - T(i,k,j,m2))**2
     END DO
   END DO
  END DO
  rmse_t = sqrt( rmse_t / float(im*jm) )
  WRITE(6,*) "rmse_t: ", rmse_t

  en_scale=1.0
  if ( rmse_t.ne.0.0 ) then
! en_limit = en_scale*0.7/rmse_t
  en_limit = en_scale*0.4/rmse_t
  else
  en_limit = en_scale
  endif
  WRITE(6,*) "en_limit: ", en_limit

  DO j = 1, jm
     DO k = 1, lm
        DO i = 1, im
           T(i,k,j,1) = T(i,k,j,1) + en_limit*( T(i,k,j,m1) - T(i,k,j,m2) )
           Q(i,k,j,1) = Q(i,k,j,1) + en_limit*( Q(i,k,j,m1) - Q(i,k,j,m2) )
           if(Q(i,k,j,1).lt.0.0) Q(i,k,j,1) = 0.000001 !unit is g/kg
           if(Q(i,k,j,1).gt.30.0) Q(i,k,j,1) = 30.0
        END DO
     END DO
  END DO

  DO j = 1, jm
     DO k = 1, lm
        DO i = 1, im+1
           U(i,k,j,1) = U(i,k,j,1) + en_limit*( U(i,k,j,m1) - U(i,k,j,m2) )
        END DO
     END DO
  END DO

  DO j = 1, jm+1
     DO k = 1, lm
        DO i = 1, im
           V(i,k,j,1) = V(i,k,j,1) + en_limit*( V(i,k,j,m1) - V(i,k,j,m2) )
        END DO
     END DO
  END DO


  DO j = 1 , jm
     DO i = 1 , im
        P(i,j,1) = P(i,j,1) + en_limit*( P(i,j,m1) - P(i,j,m2) )
     END DO
  END DO

!
! OverWRITE control file with perturbed variables 
!
  print *, "--------------------------------------------"
  print *, "P min, P max = ",minval(P(:,:,1)), maxval(P(:,:,1))
  print *, "T min, T max = ",minval(T(:,:,:,1)), maxval(T(:,:,:,1))
  print *, "Q min, Q max = ",minval(Q(:,:,:,1)), maxval(Q(:,:,:,1))
  print *, "U min, U max = ",minval(U(:,:,:,1)), maxval(U(:,:,:,1))
  print *, "V min, V max = ",minval(V(:,:,:,1)), maxval(V(:,:,:,1))

  m=1
  WRITE(input_fname,"(A,I2.2)") "fort.", m+19

  SELECT CASE ( io_form  )
  CASE ( IO_NETCDF )

    CALL check( nf_open(input_fname, NF_WRITE, iunit) )

    CALL check( nf_inq_varid (iunit, "MU", id_var) )
    CALL check( nf_put_var_real(iunit, id_var, P(:,:,1)) )
    CALL check( nf_inq_varid (iunit, "T", id_var) )
    CALL check( nf_put_var_real(iunit, id_var, T(:,:,:,1)) )
    CALL check( nf_inq_varid (iunit, "QVAPOR", id_var) )
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

  print *, "End of breeding_arw"

END PROGRAM breeding_arw

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

