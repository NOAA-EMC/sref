PROGRAM breeding

  USE dio

  IMPLICIT NONE

  INCLUDE "mpif.h"

  TYPE(dio_file) :: dfile

  INTEGER :: mem, m, m1, m2

  CHARACTER(LEN=7) :: input_fname
  CHARACTER(LEN=6) :: lvl
  INTEGER :: ios
  INTEGER :: im, jm, lm
  INTEGER :: im_save, jm_save, lm_save
  INTEGER :: i,j,l
  REAL :: rmse_t, en_limit, en_scale

  INTEGER, DIMENSION(3) :: idat
  INTEGER :: ihrst
  INTEGER :: iyear_fcst,imonth_fcst,iday_fcst,ihour_fcst
  REAL, ALLOCATABLE, DIMENSION(:,:,:) :: pd
  REAL, ALLOCATABLE, DIMENSION(:,:,:,:) :: t,q,u,v

  INTEGER :: iret, index, ntasks, rank, ikey
  INTEGER :: io_color, other_color
  INTEGER :: mc_io, mc_other, mpi_comm_all_dup

!-------------------------------------------------------------------------------

  CALL MPI_Init(iret)
  CALL MPI_Comm_size(MPI_COMM_WORLD, ntasks, iret)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, iret)
  print *, "ntasks = ", ntasks, " rank = ", rank, " iret = ", iret
  CALL MPI_Comm_Dup(MPI_COMM_WORLD, mpi_comm_all_dup, iret)
  print *, " mpi_comm_all_dup = ", mpi_comm_all_dup, " iret = ", iret
  CALL MPI_Barrier(mpi_comm_all_dup, iret)

  IF ( rank == 0 ) THEN
    ikey=1
    io_color=1
    other_color=MPI_UNDEFINED
  ELSE
    ikey=2
    io_color=MPI_UNDEFINED
    other_color=2
  END IF

  CALL MPI_Comm_split(mpi_comm_all_dup,io_color,ikey,mc_io,iret)
  CALL MPI_Comm_split(mpi_comm_all_dup,other_color,ikey,mc_other,iret)

  print *, "rank = ", rank, " io_color = ", io_color, " other_color = ",other_color, " mc_io = ", mc_io, " mc_other = ", mc_other


  IF ( rank == 0 ) THEN

  mem = 3

  DO m=1,mem

    WRITE(input_fname,"(A,I2.2)") "fort.", m+19

    CALL dio_open(dfile,trim(input_fname),"READ",iret=ios)

    CALL dio_read(dfile,"IM",im,iret=ios)
    CALL dio_read(dfile,"JM",jm,iret=ios)
    CALL dio_read(dfile,"LM",lm,iret=ios)

    IF ( m == 1 ) THEN
      print *, "im,jm,lm = ",im,jm,lm
      im_save = im
      jm_save = jm
      lm_save = lm
      ALLOCATE (PD(im,jm,mem))
      ALLOCATE (T(im,jm,lm,mem))
      ALLOCATE (Q(im,jm,lm,mem))
      ALLOCATE (U(im,jm,lm,mem))
      ALLOCATE (V(im,jm,lm,mem))
    ELSE
      IF (im/=im_save .or. jm/=jm_save .or. lm/=lm_save) THEN
         WRITE(0,*) " incorrect im or jm or in input or output file "
         WRITE(0,*) "    im,jm,lm in input file     ", im_save,jm_save,lm_save
         WRITE(0,*) "    im,jm,lm in output file    ", im,jm,lm
         CALL MPI_Abort(iret)
         STOP 911
      END IF
    END IF

    ! date check
    IF ( m == 1 ) THEN
      CALL dio_read(dfile,"IDAT",idat,iret=ios)
      CALL dio_read(dfile,"IHRST",ihrst,iret=ios)
    ELSE
      CALL dio_read(dfile,"IYEAR_FCST",iyear_fcst,iret=ios)
      CALL dio_read(dfile,"IMONTH_FCST",imonth_fcst,iret=ios)
      CALL dio_read(dfile,"IDAY_FCST",iday_fcst,iret=ios)
      CALL dio_read(dfile,"IHOUR_FCST",ihour_fcst,iret=ios)

      IF (iyear_fcst/=idat(3) .or. imonth_fcst/=idat(2) .or. iday_fcst/=idat(1) .or. ihour_fcst/=ihrst ) THEN
         WRITE(0,*) " incorrect date/time in input or output file "
         WRITE(0,*) "    idat,ihrst in input file                 ", idat,ihrst
         WRITE(0,*) "    iyear,imonth,iday,ihour in output file   ", iyear_fcst,imonth_fcst,iday_fcst,ihour_fcst
         CALL MPI_Abort(iret)
         STOP 911
      END IF
    END IF

    CALL dio_read(dfile,"PD",PD(:,:,m),iret=ios)
    DO l=1,lm
    write(lvl,"(A,I2.2,A)") "_",l,"_2D"
    CALL dio_read(dfile,"T"//lvl,T(:,:,l,m),iret=ios)
    CALL dio_read(dfile,"Q"//lvl,Q(:,:,l,m),iret=ios)
    CALL dio_read(dfile,"U"//lvl,U(:,:,l,m),iret=ios)
    CALL dio_read(dfile,"V"//lvl,V(:,:,l,m),iret=ios)
    END DO

    CALL dio_close(dfile,iret=ios)

    print *, "--------------------------------------------"
    print *, "member = ", m
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
  DO l = 15, 15
    DO j = 1, jm
     DO i = 1, im
       rmse_t = rmse_t + (T(i,j,l,m1) - T(i,j,l,m2))**2
     END DO
   END DO
  END DO
  rmse_t = sqrt( rmse_t / float(im*jm) )
  WRITE(6,*) "rmse_t: ", rmse_t

  en_scale=1.0
  en_limit = en_scale*0.7/rmse_t
  WRITE(6,*) "en_limit: ", en_limit

  DO l = 1, lm
     DO j = 1, jm
        DO i = 1, im
           T(i,j,l,1) = T(i,j,l,1) + en_limit*( T(i,j,l,m1) - T(i,j,l,m2) )
           Q(i,j,l,1) = Q(i,j,l,1) + en_limit*( Q(i,j,l,m1) - Q(i,j,l,m2) )
           U(i,j,l,1) = U(i,j,l,1) + en_limit*( U(i,j,l,m1) - U(i,j,l,m2) )
           V(i,j,l,1) = V(i,j,l,1) + en_limit*( V(i,j,l,m1) - V(i,j,l,m2) )
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

  CALL dio_open(dfile,trim(input_fname),"READWRITE",mode="OVERWRITE",iret=ios)

  CALL dio_write(dfile,"PD",PD(:,:,1),iret=ios)
  DO l=1,lm
  write(lvl,"(A,I2.2,A)") "_",l,"_2D"
  CALL dio_write(dfile,"T"//lvl,T(:,:,l,m),iret=ios)
  CALL dio_write(dfile,"Q"//lvl,Q(:,:,l,m),iret=ios)
  CALL dio_write(dfile,"U"//lvl,U(:,:,l,m),iret=ios)
  CALL dio_write(dfile,"V"//lvl,V(:,:,l,m),iret=ios)
  END DO

  CALL dio_close(dfile,iret=ios)


  END IF ! ( rank == 0 )

  CALL MPI_Finalize(iret)

  print *, "End of breeding"

END PROGRAM breeding
