PROGRAM breeding_nmb

  USE nemsio_module

  IMPLICIT NONE

#ifdef MPI
  INCLUDE "mpif.h"
#endif

  INTEGER :: iret
!
! breeding stuff
!
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

#ifdef MPI
  INTEGER :: index, ntasks, rank, ikey
  INTEGER :: io_color, other_color
  INTEGER :: mc_io, mc_other, mpi_comm_all_dup
#endif

!
! nemsio stuff
!
  ! meta1
  CHARACTER(LEN=8) :: gtype
  CHARACTER(LEN=8) :: gdatatype
  CHARACTER(LEN=8) :: modelname
  INTEGER :: version
  INTEGER :: nmeta
  INTEGER :: lmeta

  ! meta2
  INTEGER :: nrec,dimx,dimy,dimz,idate(7)
  INTEGER :: nfday,nfhour,nfminute,nfsecond,nfsecondn,nfsecondd
  INTEGER :: nframe,ntrac,nsoil
  LOGICAL :: extrameta

  INTEGER :: fieldsize

  CHARACTER(LEN=8),ALLOCATABLE :: recname(:)
  CHARACTER(LEN=16),ALLOCATABLE :: reclevtyp(:)
  INTEGER,ALLOCATABLE :: reclev(:)

  REAL,ALLOCATABLE :: glat1d(:),glon1d(:),gvlat1d(:),gvlon1d(:)
  REAL,ALLOCATABLE :: dx(:),dy(:),vcoord(:,:,:)
  REAL,ALLOCATABLE :: cpi(:),ri(:)

  INTEGER                      :: nmetavari,nmetavarr,nmetavarl,nmetavarc
  INTEGER                      :: nmetaaryi,nmetaaryr,nmetaaryl,nmetaaryc
  CHARACTER(LEN=8),ALLOCATABLE :: variname(:),varrname(:),varlname(:),varcname(:)
  CHARACTER(LEN=8),ALLOCATABLE :: aryiname(:),aryrname(:),arylname(:),arycname(:)

  INTEGER,ALLOCATABLE :: varival(:),aryilen(:),aryrlen(:),aryival(:,:)
  REAL,ALLOCATABLE :: varrval(:),aryrval(:,:)
  LOGICAL,ALLOCATABLE :: varlval(:)

  CHARACTER(LEN=8) :: vname
  CHARACTER(LEN=16) :: vlevtyp
  INTEGER :: vlev
  INTEGER :: jrec
  REAL,ALLOCATABLE :: data(:)
  REAL :: stime,etime
  REAL(KIND=8) :: timef
  INTEGER :: n

  TYPE(nemsio_gfile) :: gfile
!
!===============================================================================
!
#ifdef MPI
  CALL MPI_Init(iret)
  CALL MPI_Comm_size(MPI_COMM_WORLD, ntasks, iret)
  CALL MPI_Comm_rank(MPI_COMM_WORLD, rank, iret)
  PRINT *, "ntasks = ", ntasks, " rank = ", rank, " iret = ", iret
  CALL MPI_Comm_Dup(MPI_COMM_WORLD, mpi_comm_all_dup, iret)
  PRINT *, " mpi_comm_all_dup = ", mpi_comm_all_dup, " iret = ", iret
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

  PRINT *, "rank = ", rank, " io_color = ", io_color, " other_color = ",other_color, " mc_io = ", mc_io, " mc_other = ", mc_oth
er

  IF ( rank == 0 ) THEN
#endif

  CALL nemsio_init(iret=iret)
  PRINT *,'nemsio_init, iret=',iret
  IF (iret .ne.0 ) THEN
    STOP
  END IF


  mem = 3

  DO m=1,mem

    WRITE(input_fname,"(A,I2.2)") "fort.", m+19

  CALL nemsio_open(gfile,trim(input_fname),'READ',iret=iret)
  PRINT *,'nemsio_open, iret=',iret
  IF (iret .ne.0 ) THEN
    stop
  END IF


  ! meta1 
  CALL nemsio_getfilehead(gfile,iret=iret, &
                          gtype=gtype, &
                          gdatatype=gdatatype, &
                          modelname=modelname, &
                          version=version, &
                          nmeta=nmeta, &
                          lmeta=lmeta)
!  PRINT *, ' meta1 iret = ', iret
!  PRINT *, ' ------------------------- '
!  PRINT *, ' gtype        = ', gtype
!  PRINT *, ' gdatatype    = ', gdatatype
!  PRINT *, ' modelname    = ', modelname
!  PRINT *, ' version      = ', version
!  PRINT *, ' nmeta        = ', nmeta
!  PRINT *, ' lmeta        = ', lmeta
!  PRINT *, ' '
  

  ! meta2
  CALL nemsio_getfilehead(gfile,iret=iret, &
                          nrec=nrec, &
                          idate=idate, &
                          dimx=dimx, &
                          dimy=dimy, &
                          dimz=dimz, &
                          nfday=nfday, &
                          nfhour=nfhour, &
                          nfminute=nfminute, &
                          nfsecondn=nfsecondn, &
                          nfsecondd=nfsecondd, &
                          nframe=nframe, &
                          ntrac=ntrac, &
                          nsoil=nsoil, &
                          extrameta=extrameta)
  PRINT *, ' meta2 iret = ', iret
  PRINT *, ' ------------------------- '
  PRINT *, ' nrec         = ', nrec
  PRINT *, ' idate        = ', idate
  PRINT *, ' dimx         = ', dimx
  PRINT *, ' dimy         = ', dimy
  PRINT *, ' dimz         = ', dimz
  PRINT *, ' nfday        = ', nfday
  PRINT *, ' nfhour       = ', nfhour
  PRINT *, ' nfminute     = ', nfminute
  PRINT *, ' nfsecondn    = ', nfsecondn
  PRINT *, ' nfsecondd    = ', nfsecondd
  PRINT *, ' nframe       = ', nframe
  PRINT *, ' ntrac        = ', ntrac
  PRINT *, ' nsoil        = ', nsoil
  PRINT *, ' extrameta    = ', extrameta
  PRINT *, ' '

  im=dimx
  jm=dimy
  lm=dimz

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
#ifdef MPI
         CALL MPI_Abort(iret)
#endif
         STOP 911
      END IF
    END IF

  fieldsize=dimx*dimy

  go to 109

  allocate(recname(nrec))
  allocate(reclevtyp(nrec))
  allocate(reclev(nrec))
  allocate(vcoord(dimz+1,3,2))
  allocate(glat1d(fieldsize))
  allocate(glon1d(fieldsize))
  allocate(dx(fieldsize))
  allocate(dy(fieldsize))
  allocate(cpi(ntrac+1))
  allocate(ri(ntrac+1))
  CALL nemsio_getfilehead(gfile,iret=iret, &
                          recname=recname, &
                          reclevtyp=reclevtyp, &
                          reclev=reclev &
!                          vcoord=vcoord, &
!                          lat=glat1d, &
!                          lon=glon1d, &
!                          dx=dx, &
!                          dy=dy, &
!                          cpi=cpi, &
!                          ri=ri &
                          )

  PRINT *, ' meta3-meta12 iret = ', iret
  PRINT *, ' ------------------------- '
  DO n=1,nrec
  PRINT *, ' recname, reclevtyp, reclev ',recname(n), reclevtyp(n), reclev(n)
  END DO
!  PRINT *,'vcoord(:,1,1)=',vcoord(:,1,1)
!  PRINT *,'vcoord(:,2,1)=',vcoord(:,2,1)
!  PRINT *,'vcoord(:,3,1)=',vcoord(:,3,1)
!  PRINT *,'vcoord(:,1,2)=',vcoord(:,1,2)
!  PRINT *,'vcoord(:,2,2)=',vcoord(:,2,2)
!  PRINT *,'vcoord(:,3,2)=',vcoord(:,3,2)
!  PRINT *,'glat1d=',maxval(glat1d),minval(glat1d)
!  PRINT *,'glon1d=',maxval(glon1d),minval(glon1d)
!  PRINT *,'dx=',maxval(dx),minval(dx)
!  PRINT *,'dy=',maxval(dy),minval(dy)
  deallocate(recname)
  deallocate(reclevtyp)
  deallocate(reclev)
  deallocate(vcoord)
  deallocate(glat1d)
  deallocate(glon1d)
  deallocate(dx)
  deallocate(dy)
  deallocate(cpi)
  deallocate(ri)

  CALL nemsio_getfilehead(gfile,iret=iret, &
              nmetavari=nmetavari,nmetavarr=nmetavarr,nmetavarl=nmetavarl,nmetavarc=nmetavarc, &
              nmetaaryi=nmetaaryi,nmetaaryr=nmetaaryr,nmetaaryl=nmetaaryl,nmetaaryc=nmetaaryc)

  PRINT *, ' meta3 iret = ', iret
  PRINT *, ' ------------------------- '
  PRINT *, ' nmetavari    ', nmetavari
  PRINT *, ' nmetavarr    ', nmetavarr
  PRINT *, ' nmetavarl    ', nmetavarl
  PRINT *, ' nmetavarc    ', nmetavarc
  PRINT *, ' nmetaaryi    ', nmetaaryi
  PRINT *, ' nmetaaryr    ', nmetaaryr
  PRINT *, ' nmetaaryl    ', nmetaaryl
  PRINT *, ' nmetaaryc    ', nmetaaryc

  allocate(variname(nmetavari),varrname(nmetavarr),varlname(nmetavarl))
  allocate(aryiname(nmetaaryi),aryrname(nmetaaryr))

  allocate(varival(nmetavari),varrval(nmetavarr),varlval(nmetavarl))
  allocate(aryilen(nmetaaryi),aryrlen(nmetaaryr))

  CALL nemsio_getfilehead(gfile,iret=iret,variname=variname,varrname=varrname, &
      varlname=varlname,varival=varival,varrval=varrval,varlval=varlval,         &
      aryiname=aryiname,aryrname=aryrname,aryilen=aryilen,aryrlen=aryrlen)

  DO n=1,nmetavari
  PRINT *,'variname,varival=',variname(n),varival(n)
  END DO
  DO n=1,nmetavarr
  PRINT *,'varrname,varrval=',varrname(n),varrval(n)
  END DO
  DO n=1,nmetavarl
  PRINT *,'varlname,varlval=',varlname(n),varlval(n)
  END DO

  allocate(aryival(maxval(aryilen),nmetaaryi),aryrval(maxval(aryrlen),nmetaaryr))

  CALL nemsio_getfilehead(gfile,iret=iret,aryival=aryival,aryrval=aryrval)

  DO n=1,nmetaaryi
  PRINT *, ' aryiname, aryilen, aryival ', aryiname(n),aryilen(n),aryival(:,n) ! ,' . . . ',aryival(aryilen(n),n)
  END DO
  DO n=1,nmetaaryr
  PRINT *, ' aryrname, aryrlen, aryrval ', aryrname(n),aryrlen(n),aryrval(1,n),' . . . ',aryrval(aryrlen(n),n)
  END DO

  deallocate(variname,varrname,varlname)
  deallocate(aryiname,aryrname)
  deallocate(varival,varrval,varlval)
  deallocate(aryilen,aryrlen)
  deallocate(aryival,aryrval)

109 continue

  allocate(data(fieldsize))
!  stime=timef()
!  DO jrec=1,nrec
!     CALL nemsio_getrechead(gfile,jrec,vname,vlevtyp,vlev,iret=iret)
!     CALL nemsio_readrec(gfile,jrec,data(:),gdatatype='bin4',nframe=0,iret=iret)
!!     PRINT *,jrec,vname,vlevtyp,vlev,maxval(data(:)),minval(data(:))
!  enddo
!  etime=timef()
!  PRINT *,'after reading ',nrec,' records takes ', etime-stime


  CALL nemsio_readrecv(gfile,"dpres","hybrid sig lev",1,data,iret=iret)
  PD(:,:,m) = RESHAPE(data, (/im,jm/))
  PRINT *,minval(PD(:,:,m)),maxval(PD(:,:,m))

  DO l=1,lm
    CALL nemsio_readrecv(gfile,"tmp","mid layer",l,data,iret=iret)
    T(:,:,l,m) = RESHAPE(data, (/im,jm/))
    PRINT *,minval(T(:,:,l,m)),maxval(T(:,:,l,m))
    CALL nemsio_readrecv(gfile,"spfh","mid layer",l,data,iret=iret)
    Q(:,:,l,m) = RESHAPE(data, (/im,jm/))
    PRINT *,minval(Q(:,:,l,m)),maxval(Q(:,:,l,m))
    CALL nemsio_readrecv(gfile,"ugrd","mid layer",l,data,iret=iret)
    U(:,:,l,m) = RESHAPE(data, (/im,jm/))
    PRINT *,minval(U(:,:,l,m)),maxval(U(:,:,l,m))
    CALL nemsio_readrecv(gfile,"vgrd","mid layer",l,data,iret=iret)
    V(:,:,l,m) = RESHAPE(data, (/im,jm/))
    PRINT *,minval(V(:,:,l,m)),maxval(V(:,:,l,m))
  END DO

    print *, "--------------------------------------------"
    print *, "member = ", m
    print *, "PD min, PD max = ",minval(PD(:,:,m)), maxval(PD(:,:,m))
    print *, "T min, T max = ",minval(T(:,:,:,m)), maxval(T(:,:,:,m))
    print *, "Q min, Q max = ",minval(Q(:,:,:,m)), maxval(Q(:,:,:,m))
    print *, "U min, U max = ",minval(U(:,:,:,m)), maxval(U(:,:,:,m))
    print *, "V min, V max = ",minval(V(:,:,:,m)), maxval(V(:,:,:,m))



  CALL nemsio_close(gfile,iret=iret)
  PRINT *,'nemsio_close, iret=',iret
  IF (iret .ne.0 ) THEN
    STOP
  END IF

  deallocate(data)

  END DO

!
! Calculate perturbances
!

  m1=2    ! p memeber
  m2=3    ! n memeber

! calculate RMS errors
  rmse_t = 0.0
! approximately 850mb level
  DO l = 20, 20
    DO j = 1, jm
     DO i = 1, im
       rmse_t = rmse_t + (T(i,j,l,m1) - T(i,j,l,m2))**2
     END DO
   END DO
  END DO
  rmse_t = sqrt( rmse_t / float(im*jm) )
  WRITE(6,*) "rmse_t: ", rmse_t

  en_scale=1.0
! en_limit = en_scale*0.7/rmse_t
  en_limit = en_scale*0.5/rmse_t
  WRITE(6,*) "en_limit: ", en_limit

  DO l = 1, lm
     DO j = 1, jm
        DO i = 1, im
           T(i,j,l,1) = T(i,j,l,1) + en_limit*( T(i,j,l,m1) - T(i,j,l,m2) )
           Q(i,j,l,1) = Q(i,j,l,1) + en_limit*( Q(i,j,l,m1) - Q(i,j,l,m2) )
           U(i,j,l,1) = U(i,j,l,1) + en_limit*( U(i,j,l,m1) - U(i,j,l,m2) )
           V(i,j,l,1) = V(i,j,l,1) + en_limit*( V(i,j,l,m1) - V(i,j,l,m2) )
           if(Q(i,j,l,1).lt.0.0) Q(i,j,l,1) = 0.0000000001 !unit is kg/kg
           if(Q(i,j,l,1).gt.0.03) Q(i,j,l,1) = 0.03
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
  CALL nemsio_open(gfile,trim(input_fname),'RDWR',iret=iret)

  allocate(data(fieldsize))
 
  data = RESHAPE(PD(:,:,m), (/im*jm/))
  CALL nemsio_writerecv(gfile,"dpres","hybrid sig lev",1,data,iret=iret)

  DO l=1,lm
    data = RESHAPE(T(:,:,l,m), (/im*jm/))
    CALL nemsio_writerecv(gfile,"tmp","mid layer",l,data,iret=iret)

    data = RESHAPE(Q(:,:,l,m), (/im*jm/))
    CALL nemsio_writerecv(gfile,"spfh","mid layer",l,data,iret=iret)

    data = RESHAPE(U(:,:,l,m), (/im*jm/))
    CALL nemsio_writerecv(gfile,"ugrd","mid layer",l,data,iret=iret)

    data = RESHAPE(V(:,:,l,m), (/im*jm/))
    CALL nemsio_writerecv(gfile,"vgrd","mid layer",l,data,iret=iret)

  END DO
  deallocate(data)

  CALL nemsio_close(gfile,iret=iret)





  CALL nemsio_finalize()
  PRINT *,'nemsio final'

#ifdef MPI
  END IF ! ( rank == 0 )

  CALL MPI_Finalize(iret)
#endif

  PRINT *, "End of breeding_nmb"

END PROGRAM breeding_nmb
