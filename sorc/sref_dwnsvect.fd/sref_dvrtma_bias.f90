program rtma_bias          
!
! main program: rtma_bias             
!
! prgmmr: Bo Cui           org: np/wx20        date: 2007-08-28
! Change log:
! 12/02/2012, Jun Du: Replaced sfc pressure with 2-meter specific humidity field
!                     for SREF; increase precision for q field
!
! abstract: update bias estimation between rtma analysis and NCEP operational analyis
!
! usage:
!
!   input file: grib
!     unit 11 -    : prior bias estimation                                               
!     unit 12 -    : rtma analysis
!     unit 13 -    : ncep operational analysis
!
!   output file: grib
!     unit 51 -    : updated bias estimation pgrba file
!
!   parameters
!     fgrid -      : ensemble forecast
!     agrid -      : rtma 5km analysis data
!     bias  -      : bias estimation
!     dec_w -      : decay averaging weight 
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

integer     nvar,ivar,i,k,icstart
parameter   (nvar=4)

real,       allocatable :: agrid(:),fgrid(:),bias(:)
logical(1), allocatable :: lbms(:),lbmsout(:)
real        dmin,dmax

real        dec_w

integer     ifile,afile,cfile,ofile
integer     jpds(200),jgds(200),jens(200),kpds(200),kgds(200),kens(200)
integer     kpdsout(200),kgdsout(200),kensout(200)
integer     pds5(nvar),pds6(nvar),pds7(nvar)
real        minval(4)

! variables: u10m v10m t2m q2m
 
data pds5/33, 34, 11, 51/
data pds6/105,105,105, 105/
data pds7/10, 10,  2, 2/
data minval/-100,-100,200,0.0000000001/

integer     maxgrd,ndata                                
integer     index,j,n,iret,jret             
character*7 cfortnn

namelist/message/icstart,dec_w
read(5,message,end=1020)
write(6,message)
!icstart=1 
!dec_w=0.30

ifile=11
afile=12
cfile=13
ofile=51

! index=0, to get index buffer from the grib file not the grib index file
! j=0, to search from beginning,  <0 to read index buffer and skip -1-j messages
! lbms, logical*1 (maxgrd or kf) unpacked bitmap if present

index=0
j=0
iret=0
jret=0
jpds=-1  
jgds=-1
jens=-1

! set the fort.* of intput files

write(cfortnn,'(a5,i2)') 'fort.',afile
call baopenr(afile,cfortnn,iret)
if (iret.ne.0) then; print*,'there is no rtma data, stop!'; endif
if (iret.ne.0) goto 1020

write(cfortnn,'(a5,i2)') 'fort.',cfile
call baopenr(cfile,cfortnn,iret)
if (iret.ne.0) then; print*,'there is no NCEP analysis data, stop!'; endif
if (iret.ne.0) goto 1020

if(icstart.eq.0) then
  write(cfortnn,'(a5,i2)') 'fort.',ifile
  call baopenr(ifile,cfortnn,iret)
  if (iret.ne.0) then; print*,'there is no bias estimation data, please check!'; endif
endif

! set the fort.* of output file

write(cfortnn,'(a5,i2)') 'fort.',ofile
call baopenw(ofile,cfortnn,iret)
if (iret.ne.0) then; print*,'there is no output bias data, stop!'; endif
if (iret.ne.0) goto 1020

! find grib message. input: jpds,jgds and jens.  output: kpds,kgds,kens
! ndata: integer number of bites in the grib message

call getgbeh(afile,index,j,jpds,jgds,jens,ndata,maxgrd,k,kpds,kgds,kens,iret)
if (iret.ne.0) then; print*,' getgbeh ,afile,index,iret =',afile,index,iret; endif
if (iret.ne.0) goto 1020

allocate (agrid(maxgrd),fgrid(maxgrd),bias(maxgrd),lbms(maxgrd),lbmsout(maxgrd))

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
  kpdsout=0 
  kgdsout=0 
  kensout=0 

  bias=-9999.9999
  agrid=-9999.9999
  fgrid=-9999.9999

  ! read and process variable of input data

  jpds(5)=pds5(ivar)
  jpds(6)=pds6(ivar)
  jpds(7)=pds7(ivar)

  ! get initialized bias estimation

  if(icstart.eq.1) then
    print *, '----- Cold Start for Bias Estimation -----'
    print*, '  '
    bias=0.0
  else
    print *, '----- Initialized Bias for Current Time -----'
    call getgbe(ifile,index,maxgrd,j,jpds,jgds,jens,ndata,k,kpds,kgds,kens,lbms,bias,iret)
    if (iret.ne.0) then; print*, 'iret ne 0', jpds(5),jpds(6),jpds(7); endif
    if (iret.eq.99) then; print*, 'there is no variable ', jpds(5),jpds(6),jpds(7); endif
    if (iret.ne.0) bias=0.0

    call grange(maxgrd,lbms,bias,dmin,dmax)
    print*, 'Irec pds5 pds6 pds7 pds8 pds9 pd10 pd11 pd14 e2  e3  ndata   Maximun    Minimum    data'
    print '(i4,8i5,2i4,i8,3f10.2)',ivar,(kpds(i),i=5,11),kpds(14),(kens(i),i=2,3),maxgrd,dmax,dmin,bias(30000)
    print*, '  '
  endif

  ! get rtma data

  index=0
  j=0
  iret=0

  print *, '----- RTMA Analysis for Current Time ------'

  call getgbe(afile,index,maxgrd,j,jpds,jgds,jens,ndata,k,kpds,kgds,kens,lbms,agrid,iret)

  if (iret.ne.0) then; print*, 'iret ne 0, jpds5, 6, 7= ', jpds(5),jpds(6),jpds(7); endif
  if (iret.eq.99) then; print*, 'there is no variable ', jpds(5),jpds(6),jpds(7); endif
    call grange(maxgrd,lbms,agrid,dmin,dmax)
    print *,"Compare:",dmin,minval(ivar)
    if(dmin.lt.minval(ivar)) then
    print *,"Ribbons of shame!"
    j=-1-k
    print *,"Try again at j=",j
    call getgbe(afile,index,maxgrd,j,jpds,jgds,jens,ndata,k,kpds,kgds,kens,lbms,agrid,iret)
    call grange(maxgrd,lbms,agrid,dmin,dmax)
    endif

  ! if there is no rtma analyis, save previous data message 

  if (iret.ne.0) then
    kpdsout=kpds
    kgdsout=kgds
    kensout=kens
    lbmsout=lbms
  endif

  if (iret.ne.0) goto 200

  call grange(maxgrd,lbms,agrid,dmin,dmax)
  print*, 'Irec pds5 pds6 pds7 pds8 pds9 pd10 pd11 pd14 e2  e3  ndata   Maximun    Minimum    data'
  print '(i4,8i5,2i4,i8,3f10.2)',ivar,(kpds(i),i=5,11),kpds(14),(kens(i),i=2,3),maxgrd,dmax,dmin,agrid(30000)
  print*, '  '

  ! get operational analysis

  index=0
  j=0
  iret=0

  print *, '----- NCEP operational analysis for current Time ------'
  call getgbe(cfile,index,maxgrd,j,jpds,jgds,jens,ndata,k,kpdsout,kgdsout,kensout,lbmsout,fgrid,iret)

  if (iret.ne.0) then; print*, 'iret ne 0, jpds5, 6, 7= ', jpds(5),jpds(6),jpds(7); endif
  if (iret.eq.99) then; print*, 'there is no variable ', jpds(5),jpds(6),jpds(7); endif

  ! if there is no rtma analyis, save previous data message 

  if (iret.ne.0) then
    kpdsout=kpds
    kgdsout=kgds
    kensout=kens
    lbmsout=lbms
  endif

  if (iret .ne. 0) goto 200

  call grange(maxgrd,lbmsout,fgrid,dmin,dmax)
  print*, 'Irec pds5 pds6 pds7 pds8 pds9 pd10 pd11 pd14 e2  e3  ndata   Maximun    Minimum    data'
  print '(i4,8i5,2i4,i8,3f10.2)',ivar,(kpdsout(i),i=5,11),kpdsout(14),(kensout(i),i=2,3),maxgrd,dmax,dmin,fgrid(30000)
  print*, '  '

  ! apply the decay average
  call decay(bias,fgrid,agrid,maxgrd,dec_w)
  print *,bias(30000),fgrid(30000),agrid(30000),maxgrd,dec_w
  ! save the first moment bias, give the bias a time one day before
  200 continue
 
  ! increase the decimal scale factor to 2

  if(kpdsout(5).eq.11.or.kpdsout(5).eq.33.or.kpdsout(5).eq.34) then
    kpdsout(22)=2
  endif
  if(kpdsout(5).eq.51) then
    kpdsout(22)=6
  endif

!  do 999 n=1,200
!  print *,n,kgdsout(n),kpdsout(n),kensout(200) 
!999 continue
  print *,ofile,maxgrd,kpdsout(200),kgdsout(200),kensout(200),lbmsout(maxgrd),bias(maxgrd),jret
  print *, '----- Output Bias Estimation for Current Time ------'
  call putgbe(ofile,maxgrd,kpdsout,kgdsout,kensout,lbmsout,bias,jret)

  print *,"Did putgbe work?"
  call grange(maxgrd,lbms,bias,dmin,dmax)
  print*, 'Irec pds5 pds6 pds7 pds8 pds9 pd10 pd11 pd14 e2  e3  ndata   Maximun    Minimum    data'
  print '(i4,8i5,2i4,i8,3f10.2)',ivar,(kpdsout(i),i=5,11),kpdsout(14),(kensout(i),i=2,3),maxgrd,dmax,dmin,bias(30000)
  print*, '  '

! end of bias estimation for one forecast lead time

100 continue

enddo

call baclose(ifile,iret)
call baclose(afile,iret)
call baclose(cfile,iret)
call baclose(ofile,iret)

print *,'Bias Estimation Successfully Complete'

stop

1020  continue

print *,'Wrong Data Input, Output or Wrong Message Input'

stop
end

subroutine decay(aveeror,fgrid,agrid,maxgrd,dec_w)

!     apply the decaying average scheme
!
!     parameters
!                  fgrid  ---> ensemble forecast
!                  agrid  ---> analysis data
!                  aveeror---> bias estimation
!                  dec_w  ---> decay weight

implicit none

integer maxgrd,ij
real aveeror(maxgrd),fgrid(maxgrd),agrid(maxgrd)
real dec_w           

do ij=1,maxgrd
  if(fgrid(ij).gt.-999.0.and.fgrid(ij).lt.999999.0.and.agrid(ij).gt.-999.0.and.agrid(ij).lt.999999.0) then
      if(aveeror(ij).gt.-999.0.and.aveeror(ij).lt.999999.0) then
        aveeror(ij)= (1-dec_w)*aveeror(ij)+dec_w*(fgrid(ij)-agrid(ij))
      else
        aveeror(ij)= dec_w*(fgrid(ij)-agrid(ij))
      endif
  else
    if(aveeror(ij).gt.-999.0 .and.aveeror(ij).lt.999999.0) then
      aveeror(ij)= aveeror(ij)                   
    else
      aveeror(ij)= 0.0                                
    endif
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

