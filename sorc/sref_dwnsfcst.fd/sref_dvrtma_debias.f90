program dvrtma_debias   
!
! main program: dvrtma_debias    
!
! prgmmr: Bo Cui           org: np/wx20        date: 2006-07-21
!
! abstract: downscale forecast by removing downscaling vector from it 
!
! Change log: 
! 12/09/11, Jun Du:   Added a check point to prevent humidity field from being <0.0
!
! usage:
!
!   input file: grib
!     unit 11 -    : downscaling vector                                                 
!     unit 12 -    : ensemble forecast
!
!   output file: grib
!     unit 51 -    : downscaled forecast  
!
!   parameters
!     fgrid -      : forecast
!     bias  -      : downscaling vector
!     nvar  -      : number of variables

! programs called:
!   baopenr          grib i/o
!   baopenw          grib i/o
!   baclose          grib i/o
!   getgbeh          grib reader
!   getgbe           grib reader
!   putgbe           grib writer

! exit states:
!   cond =   0 - successful run
!   cond =   1 - I/O abort
!
! attributes:
!   language: fortran 90
!
!$$$

implicit none

integer     nvar,ivar,i,icstart
parameter   (nvar=4)

real,       allocatable :: fgrid(:),bias(:)
logical(1), allocatable :: lbms(:),lbmsout(:)
real        dmin,dmax

integer     ifile,cfile,ofile
integer     jpds(200),jgds(200),jens(200),kpds(200),kgds(200),kens(200)
integer     kpdsout(200),kgdsout(200),kensout(200)
integer     pds5(nvar),pds6(nvar),pds7(nvar)

! variables: u10m v10m t2m q2m
 
!data pds5/1, 33, 34, 11/
!data pds6/1,105,105,105/
!data pds7/0, 10, 10,  2/
!data pds5/33, 34, 11, 51/
!data pds6/105,105,105, 105/
!data pds7/10, 10,  2, 2/
data pds5/11, 51, 33, 34/
data pds6/105,105,105, 105/
data pds7/2, 2, 10, 10/

integer     maxgrd,ndata                                
integer     index,j,n,iret,jret             
character*7 cfortnn

!namelist/message/icstart

!read(5,message,end=1020)
!write(6,message)
icstart=0

ifile=11
cfile=12
ofile=51

! index=0, to get index buffer from the grib file not the grib index file
! j=0, to search from beginning,  <0 to read index buffer and skip -1-j messages
! lbms, logical*1 (maxgrd or kf) unpacked bitmap if present

index=0
j=0
iret=0

! set the fort.* of intput files

if(icstart.eq.0) then
  write(cfortnn,'(a5,i2)') 'fort.',ifile
  call baopenr(ifile,cfortnn,iret)
  if (iret.ne.0) then; print*,'there is no rtma downscaling vector'; endif
endif

write(cfortnn,'(a5,i2)') 'fort.',cfile
call baopenr(cfile,cfortnn,iret)
if (iret.ne.0) then; print*,'there is no forecast, stop!'; endif
if (iret.ne.0) goto 1020 

! set the fort.* of output file

write(cfortnn,'(a5,i2)') 'fort.',ofile
call baopenw(ofile,cfortnn,iret)
if (iret.ne.0) then; print*,'there is no downscaled output, stop!'; endif
if (iret.ne.0) goto 1020

! find grib message. input: jpds,jgds and jens.  output: kpds,kgds,kens
! ndata: integer number of bites in the grib message

jpds=-1  
jgds=-1
jens=-1

call getgbeh(cfile,index,j,jpds,jgds,jens,ndata,maxgrd,j,kpds,kgds,kens,iret)
if (iret.ne.0) then; print*,' getgbeh ,cfile,index,iret =',cfile,index,iret; endif
if (iret.ne.0) goto 1020

allocate (fgrid(maxgrd),bias(maxgrd),lbms(maxgrd),lbmsout(maxgrd))

do ivar = 1, nvar  

  index=0
  j=0
  iret=0
  jpds=-1
  jgds=-1
  jens=-1
  kpds=-1
  kgds=-1
  kens=-1
  kpdsout=-1
  kgdsout=-1
  kensout=-1

  bias=0.0            
  fgrid=-9999.9999

  ! read and process variable of input data

  jpds(5)=pds5(ivar)
  jpds(6)=pds6(ivar)
  jpds(7)=pds7(ivar)

  ! get initialized downscaling vector

  if(icstart.eq.1) then
    print *, '----- Cold Start for Bias Correction -----'
    bias=0.0
  else
    print *, '----- Initialized Bias for Current Time -----'
    call getgbe(ifile,index,maxgrd,j,jpds,jgds,jens,ndata,j,kpds,kgds,kens,lbms,bias,iret)
    if (iret.ne.0) then; print*, 'iret ne 0', jpds(5),jpds(6),jpds(7); endif
    if (iret.eq.99) then; print*, 'there is no variable ', jpds(5),jpds(6),jpds(7); endif
    if (iret.ne.0) bias=0.0 

    call grange(maxgrd,lbms,bias,dmin,dmax)
    print*, 'Irec pds5 pds6 pds7 pds8 pds9 pd10 pd11 pd14 e2  e3  ndata   Maximun    Minimum    data'
    print '(i4,8i5,2i4,i8,3f10.2)',ivar,(kpds(i),i=5,11),kpds(14),(kens(i),i=2,3),maxgrd,dmax,dmin,bias(30000)
    print*, '  '
  endif

  ! get operational forecast 

  index=0
  j=0
  iret=0

  print *, '----- Operational Forecast for Current Time ------'
  call getgbe(cfile,index,maxgrd,j,jpds,jgds,jens,ndata,j,kpdsout,kgdsout,kensout,lbmsout,fgrid,iret)

  if (iret.ne.0) then; print*, 'iret ne 0', jpds(5),jpds(6),jpds(7); endif
  if (iret.eq.99) then; print*, 'there is no variable ', jpds(5),jpds(6),jpds(7); endif
  if (iret.ne.0) goto 100

  call grange(maxgrd,lbmsout,fgrid,dmin,dmax)
  print*, 'Irec pds5 pds6 pds7 pds8 pds9 pd10 pd11 pd14 e2  e3  ndata   Maximun    Minimum    data'
  print '(i4,8i5,2i4,i8,3f10.2)',ivar,(kpdsout(i),i=5,11),kpdsout(14),(kensout(i),i=2,3),maxgrd,dmax,dmin,fgrid(30000)
  print*, '  '

  ! apply downscaling approach

  call debias(bias,fgrid,maxgrd)

! Prevent humidity from being less than 0.0
! if(jpds(5).eq.51) then
  if(kpdsout(5).eq.51) then
   do i=1,maxgrd
    if(fgrid(i).lt.0.0) fgrid(i) = 0.0001
   enddo
  endif

  ! save the downscaled output

  print *, '----- Output Downscaled Forecast -----'
  call putgbe(ofile,maxgrd,kpdsout,kgdsout,kensout,lbmsout,fgrid,jret)

  call grange(maxgrd,lbmsout,fgrid,dmin,dmax)
  print*, 'Irec pds5 pds6 pds7 pds8 pds9 pd10 pd11 pd14 e2  e3  ndata   Maximun    Minimum    data'
  print '(i4,8i5,2i4,i8,3f10.2)',ivar,(kpdsout(i),i=5,11),kpdsout(14),(kensout(i),i=2,3),maxgrd,dmax,dmin,fgrid(30000)
  print*, '  '

! end of downscaling process for one forecast lead time

100 continue

enddo

call baclose(ifile,iret)
call baclose(cfile,iret)
call baclose(ofile,iret)

print *,'Downscaling Process Successfully Complete'

stop

1020  continue

print *,'Wrong Data Input, Output or Wrong Message Input'

stop
end

subroutine debias(bias,fgrid,maxgrd)

!     apply the downscaling process    
!
!     parameters
!                  fgrid  ---> ensemble forecast
!                  bias   ---> downscaling vector

implicit none

integer maxgrd,ij
real bias(maxgrd),fgrid(maxgrd)

do ij=1,maxgrd
  if(fgrid(ij).gt.-999.0.and.fgrid(ij).lt.999999.0.and.bias(ij).gt.-999.0.and.bias(ij).lt.999999.0) then
    fgrid(ij)=fgrid(ij)-bias(ij)
  else
    fgrid(ij)=fgrid(ij)
  endif
enddo

return
end

subroutine grange(n,ld,d,dmin,dmax)
!$$$  SUBPROGRAM DOCUMENTATION BLOCK
!
! SUBPROGRAM: GRANGE(N,LD,D,DMIN,DMAX)
!   PRGMMR: YUEJIAN ZHU       ORG:NP23          DATE: 97-03-17
!
! ABSTRACT: THIS SUBROUTINE WILL ALCULATE THE MAXIMUM AND
!           MINIMUM OF A ARRAY
!
! PROGRAM HISTORY LOG:
!   97-03-17   YUEJIAN ZHU (WD20YZ)
!
! USAGE:
!
!   INPUT ARGUMENTS:
!     N        -- INTEGER
!     LD(N)    -- LOGICAL OF DIMENSION N
!     D(N)     -- REAL ARRAY OF DIMENSION N
!
!   OUTPUT ARGUMENTS:
!     DMIN     -- REAL NUMBER ( MINIMUM )
!     DMAX     -- REAL NUMBER ( MAXIMUM )
!
! ATTRIBUTES:
!   LANGUAGE: FORTRAN
!
!$$$
logical(1) ld(n)
real d(n)
real dmin,dmax
integer i,n
dmin=1.e38
dmax=-1.e38
do i=1,n
  if(ld(i)) then
    dmin=min(dmin,d(i))
    dmax=max(dmax,d(i))
  endif
enddo
return
end

