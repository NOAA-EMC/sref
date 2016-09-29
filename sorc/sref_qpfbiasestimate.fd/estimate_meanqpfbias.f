      program qpfbias

C$$  MAIN PROGRAM DOCUMENTATION BLOCK
C      
C MAIN PROGRAM: qpfbias
C PRGMMR: Jun Du          ORG: NP21       DATE: Oct. 15, 2010
C
C ABSTRACT: 
c  1. read in current forecast 
c  2. read in past forecast and observed PDFs
c  3. find and write out the accumulated forecast and observed PDFs 
c     during a past period (e.g., 20 days)
c  4. use PDF-matching method to calibrate raw precipitation forecasts
c  5. write out new debiased forecasts
c
c PROGRAM HISTORY LOG:
c  10/15/2010, Jun Du, initial program
c  11/03/2011, Jun Du, modified to the new 2012 NEMS-based SREF (v6.0.0).
c              This particular version is for ensemble mean or a single member
c
c INPUT:
c
c OUTPUT:
c
c SUBPROGRAMS CALLED: 
c  getname
c  grange
c  CPCOEF
c  stpint
c  qtpint (substitute for stpint)
c
C ATTRIBUTES:
c   language: Fortran 77
c   machine: IBM-SP
C
C$$$

C OUTPUT:
C KPDS(25) INTEGER*4
C ARRAY TO CONTAIN THE ELEMENTS OF THE PRODUCT DEFINITION SEC
C (VERSION 1)
C KPDS(1)  - ID OF CENTER
C KPDS(2)  - MODEL IDENTIFICATION (SEE "GRIB" TABLE 1)
C KPDS(3)  - GRID IDENTIFICATION (SEE "GRIB" TABLE 2)
C KPDS(4)  - GDS/BMS FLAG
C            BIT       DEFINITION
C            25        0 - GDS OMITTED
C                      1 - GDS INCLUDED
C            26        0 - BMS OMITTED
C                      1 - BMS INCLUDED
C            NOTE:- LEFTMOST BIT = 1,
C                   RIGHTMOST BIT = 32
C KPDS(5)  - INDICATOR OF PARAMETER (SEE "GRIB" TABLE 5)
C KPDS(6)  - TYPE OF LEVEL (SEE "GRIB" TABLES 6 & 7)
C KPDS(7)  - HEIGHT,PRESSURE,ETC  OF LEVEL
C KPDS(8)  - YEAR INCLUDING CENTURY
C KPDS(9)  - MONTH OF YEAR
C KPDS(10) - DAY OF MONTH
C KPDS(11) - HOUR OF DAY
C KPDS(12) - MINUTE OF HOUR
C KPDS(13) - INDICATOR OF FORECAST TIME UNIT (SEE "GRIB"
C            TABLE 8)
C KPDS(14) - TIME 1               (SEE "GRIB" TABLE 8A)
C KPDS(15) - TIME 2               (SEE "GRIB" TABLE 8A)
C KPDS(16) - TIME RANGE INDICATOR (SEE "GRIB" TABLE 8A)
C KPDS(17) - NUMBER INCLUDED IN AVERAGE
C KPDS(18) - EDITION NR OF GRIB SPECIFICATION
C KPDS(19) - VERSION NR OF PARAMETER TABLE
C
C KPDS(20)  - NR MISSING FROM AVERAGE/ACCUMULATION
C KPDS(21)  - CENTURY OF REFERENCE TIME OF DATA
C KPDS(22)  - UNITS DECIMAL SCALE FACTOR
C KPDS(23)  - SUBCENTER NUMBER
C KPDS(24)  - PDS BYTE 29, FOR NMC ENSEMBLE PRODUCTS
C         128 IF FORECAST FIELD ERROR
C          64 IF BIAS CORRECTED FCST FIELD
C          32 IF SMOOTHED FIELD
C         WARNING: CAN BE COMBINATION OF MORE THAN 1
C KPDS(25)  - PDS BYTE 30, NOT USED
C KPDS(26-35)  - RESERVED
C KPDS(36-N)   - CONSECUTIVE BYTES EXTRACTED FROM PROGRAM
C                DEFINITION SECTION (PDS) OF GRIB MESSAGE
c     kpds(1)=7
c     kpds(2)=130
c     kpds(3)=212
c     kpds(4)=0
c     kpds(5)=variable dependent
c     kpds(6)=variable dependent
c     kpds(7)=variable dependent
c     kpds(8)=iyr-(iyr/100)*100
c     kpds(9)=imon
c     kpds(10)=idy
c     kpds(11)=ihr
c     kpds(12)=0
c     kpds(13)=1
c     kpds(14)=variable dependent
c     kpds(15)=variable dependent
c     kpds(16)=variable dependent
c     kpds(17)=-1
c     kpds(18)=1
c     kpds(19)=2
c     kpds(20)=-1
c     kpds(21)=21
c     kpds(22)=variable dependent (units decimal scale factor/precision)
c     kpds(23)=2  !if ensemble data
c     kpds(23)=0  !if not ensemble data
c     kpds(24)=128
c     kpds(25)=-1
       
c     kgds(20)=255
CCCCCCCCCCCCCCCCCCCC

 
C nfile - 1-21: ens mems
      parameter(mem=1,ntime=29,ncat=6,isample=20) ! 20 days
c     parameter(mem=1,ntime=29,ncat=6,isample=14)  ! 14 days
c     parameter(mem=1,ntime=29,ncat=6,isample=7)  ! 7 days
      parameter(jf=185*129,lon=185,lat=129,ngrid=lon*lat)  !for Grid 212
c     parameter(jf=512*256,lon=185,lat=129,ngrid=lon*lat)

      dimension f(jf),data(ngrid,mem),thresh(ncat),thresh_r(ncat),
     &r(ncat),r_r(ncat),fcst_past(mem,ncat,isample),fcst(ncat),
     &obs_past(ncat,isample),obs(ncat),ppt_new(ngrid,mem),
     &fcst_BC(ngrid)
      dimension kkpds(25)
      dimension kftime(mem,ntime),jpds(25),jgds(22),kpds(25),kgds(22)
      integer yy1,mm1,dd1,hh1,type
      logical lb(jf)
      character*2  date
      character*13  fname1
      character*15  fname2
      character*21 fcst_new1
      namelist/namin/yy1,mm1,dd1,hh1,type

      data thresh/1.0,0.75,0.5,0.25,0.1,0.01/  !in inch (for 24h-apcp)
c     data thresh/0.5,0.40,0.3,0.20,0.1,0.01/  !in inch (for 03h-apcp)
      do i=1,ncat
       thresh(i)=thresh(i)*25.4  !convert to mm if used inch
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C decoding current-day forecast grib data
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c Passing over date information
      read(5,namin,end=1000)
1000  continue
      print*,yy1,mm1,dd1,hh1,type
      read(7,*) ((kftime(i,j),j=1,ntime),i=1,mem)
      do i=1,mem
       print *,(kftime(i,j),j=1,ntime)
      enddo

      read(8,*)
      read(8,*)
      read(8,*)

      if (type.eq.1) then
       do k=2,isample
        read(8+k,*)
        read(8+k,*)
        read(8+k,*)
       enddo
      endif

      do 10 nt=1,ntime

      read(8,*)
      read(8,*)
      read(8,*)

      if (type.eq.1) then
       do k=2,isample
        read(8+k,*)
        read(8+k,*)
        read(8+k,*)
       enddo
      endif

        print *,' '
        print *,'----------------------------------------------- '
        print *,'NEW TIME: nt = ',nt
        print *,' '

        j=0 
        data=0.0
        lugb=9         
        lugi=39     

        do 5 nf=1,mem
         read(8,*) (fcst_past(nf,i,1),i=1,ncat)
         print*, (fcst_past(nf,i,1),i=1,ncat)
         if (type.eq.1) then
          do k=2,isample
           read(8+k,*) (fcst_past(nf,i,k),i=1,ncat)
           print*, (fcst_past(nf,i,k),i=1,ncat)
          enddo
         endif

        jpds=-1
        jgds=-1
C specify date information here:
c       jpds(8)=yy1-(yy1/100)*100
c       jpds(8)=yy1
c       jpds(9)=mm1
        jpds(10)=dd1
        jpds(11)=hh1

        print *,' '
        print *,'----------'
        print *,'New File: nt = ',nt,' nf= ',nf
        print *,' '

       jpds(5)=61
       jpds(6)=1
       jpds(7)=0

       if(kftime(nf,nt).gt.3) then
        jpds(14)=kftime(nf,nt)-3
       else
        jpds(14)=0
       endif
       jpds(15)=kftime(nf,nt)

       print*,'jpds10,jpds11,jpds14,jpds15,jpds5,jpds6,jpds7=',
     & jpds(10),jpds(11),jpds(14),jpds(15),jpds(5),jpds(6),jpds(7)

      lugb=lugb+1
      lugi=lugi+1
      print *,'lugb= ',lugb,' lugi= ',lugi
      call getname(nf,fname1,fname2)
      print *,' fname1= ',fname1
      print *,' fname2= ',fname2
      call baopenr(lugb,fname1,ierr)
      call baopenr(lugi,fname2,ierr)
      call getgb(lugb,lugi,jf,j,jpds,jgds,
     &             kf,k,kpds,kgds,lb,f,iret)
      call baclose(lugb,ierr)
      call baclose(lugi,ierr)
c     print *,'lugb= ',lugb,' lugi= ',lugi
 
      call grange(kf,lb,f,dmin,dmax) 
      print *,'immediately after call to getgb, iret= ',iret
     &,' j= ',j,' k= ',k,' nf= ',nf,' nt= ',nt
     &,' j5= ',jpds(5),' j6= ',jpds(6),' j7= ',jpds(7)
     &,' j8= ',jpds(8),' j9= ',jpds(9),' j10= ',jpds(10)
     &,' j11= ',jpds(11),' j14= ',jpds(14),' j15= ',jpds(15)
      print '(i4,2x,9i5,i8,2g12.4)',
     &k,(kpds(i),i=5,11),kpds(14),kpds(15),kf,dmin,dmax

      if (iret.eq.0) then
       do ipt=1,ngrid
        data(ipt,nf) = f(ipt)
        if(nt.eq.1) then
         do i=1,25
          kkpds(i)=kpds(i)
         enddo
        endif
       enddo
      else
       print *, "Something wrong in decoding ensemble forecasts!!!!"
      endif
5     continue    !files (nf: mems)

      read(8,*) (obs_past(i,1),i=1,ncat)
      print*, (obs_past(i,1),i=1,ncat)
      if (type.eq.1) then
       do k=2,isample
        read(8+k,*) (obs_past(i,k),i=1,ncat)
        print*, (obs_past(i,k),i=1,ncat)
       enddo
      endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C bias estimation
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        write(date,'(i2.2)') (nt-1)*3+3

         fcst_new1='BC_smn.pgrb212.c1.f' // date
        call baopen(10,fcst_new1,ierr)

         do nf=1,mem
C Find and write out the new forecast and observed PDFs using the decaying 
C average method or simply accumulating all points from the past
      do i=1,ncat
       fcst(i)=0.0
       if(nf.eq.1) obs(i)=0.0
      if (type.eq.1) then
c      fcst(i)=(1-w)*fcst_acum(nf,i)+w*fcst_past(nf,i,1)
c      if(nf.eq.1) obs(i)=(1-w)*obs_acum(i)+w*obs_past(i,1)
       do k=1,isample
        fcst(i)=fcst(i)+fcst_past(nf,i,k)
        if(nf.eq.1) obs(i)=obs(i)+obs_past(i,k)
       enddo
      else
       fcst(i)=fcst_past(nf,i,1)
       if(nf.eq.1) obs(i)=obs_past(i,1)
      endif
      enddo

      print*,'nt=',nt, ' nf=',nf
      print*,'thresh=',(thresh(i),i=1,ncat)
      print*,'obs=',(obs(i),i=1,ncat)
      print*,'fcst=',(fcst(i),i=1,ncat)
      if(nf.eq.1) write(50,*) (nt-1)*3+3
      write(50,*) (fcst(i),i=1,ncat)
      if(nf.eq.mem) write(50,*) (obs(i),i=1,ncat)

C Find the coefficient based on fcst vs. obs
      call CPCOEF(thresh,obs,fcst,r,ncat,2)
      print*,'raw coefficients',(r(i),i=1,ncat)

c In some occasions, to prevent heavier rain amount from being 
c reduced given the fact of "normally being too dry" 
C Do not increase:
c     if(nf.ge.1.and.nf.le.5) then
c      r(1)=1.0
c      r(2)=1.0
c      r(3)=1.0
c     endif
C Do not reduce:
c     do i=1,4
c     do i=1,1
c      if(r(i).lt.1.0) r(i)=1.0
c     enddo
c To prevent from being corrected too much:
      do i=1,ncat
       if(r(i).gt.1.2) r(i)=1.2
       if(r(i).lt.0.6) r(i)=0.6
      enddo
c If no correction:
c     do i=1,ncat
c      r(i)=1.0
c     enddo
      print*,'modified coefficients',(r(i),i=1,ncat)

c reverse the order of thresholds
      do ii = 1, ncat
       thresh_r(ii) = thresh(ncat-ii+1)
       r_r(ii) = r(ncat-ii+1)
      enddo

C Use the PDF-matching method to calibrate raw precipitation forecasts
      do ipt=1,ngrid
       ppt_old=data(ipt,nf)
       call qtpint(thresh_r,r_r,ncat,2,ppt_old,bbb,1) !Jim Purser's routine
c      call stpint(thresh_r,r_r,ncat,2,ppt_old,bbb,1,aux,naux) !essl lib
c      ppt_new(ipt,nf)=ppt_old*bbb
       ppt_new(ipt,nf)=data(ipt,nf)*bbb
       if(ipt.eq.(lat+lon)/2) then
        print*,'ppt_old=',ppt_old
        print*,'ratio=',bbb
        print*,'ppt_new=',ppt_new(ipt,nf)
       endif
       if(bbb.gt.1.) then
        print*,'ppt_old=',ppt_old
        print*,'ratio=',bbb
        print*,'ppt_new=',ppt_new(ipt,nf)
       endif
      enddo

C Write out new debiased forecasts
         jpds=-1
         jgds=-1
         jpds(5)=61
         jpds(6)=1
         jpds(7)=0
         jpds(14)=kftime(1,nt)-3  !same as ensemble member 1
         jpds(15)=kftime(1,nt)
         jpds(21)=21
         jpds(22)=4

         do i=1,25
          kpds(i)=kkpds(i)
         enddo
c        kpds(4)=0  !bitmap
         kpds(21)=21
         kpds(22)=4
         kpds(23)=0  !if not ensemble data
         kpds(8)=yy1-(yy1/100)*100
         kpds(9)=mm1
         kpds(10)=dd1
         kpds(11)=hh1
         kpds(12)=0
         kpds(14)=(nt-1)*3
         kpds(15)=(nt-1)*3+3
         kpds(5)=61
         kpds(6)=1
         kpds(7)=0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C write out debiased new ensemble fcsts:
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
           do ipt=1,ngrid
            fcst_BC(ipt)=ppt_new(ipt,nf)
            if(fcst_BC(ipt).lt.0.) fcst_BC(ipt)=0.  
           enddo
           if(nf.eq.1) call putgb(10,jf,kpds,kgds,lb,fcst_BC,iret)
         enddo

           call baclose(10,ierr)

10    continue  !time

9999	END

C************************************************
      subroutine grange(n,ld,d,dmin,dmax)
      logical ld
      dimension ld(n),d(n)
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
 
C***************************************************
      subroutine getname (fnum,fname1,fname2)
      character*13  fname1
      character*15  fname2
      integer    fnum
      if (fnum.eq.1)  fname1='r_gribawips01'
      if (fnum.eq.2)  fname1='r_gribawips02'
      if (fnum.eq.3)  fname1='r_gribawips03'
      if (fnum.eq.4)  fname1='r_gribawips04'
      if (fnum.eq.5)  fname1='r_gribawips05'
      if (fnum.eq.6)  fname1='r_gribawips06'
      if (fnum.eq.7)  fname1='r_gribawips07'
      if (fnum.eq.8)  fname1='r_gribawips08'
      if (fnum.eq.9)  fname1='r_gribawips09'
      if (fnum.eq.10) fname1='r_gribawips10'
      if (fnum.eq.11) fname1='r_gribawips11'
      if (fnum.eq.12) fname1='r_gribawips12'
      if (fnum.eq.13) fname1='r_gribawips13'
      if (fnum.eq.14) fname1='r_gribawips14'
      if (fnum.eq.15) fname1='r_gribawips15'
      if (fnum.eq.16) fname1='r_gribawips16'
      if (fnum.eq.17) fname1='r_gribawips17'
      if (fnum.eq.18) fname1='r_gribawips18'
      if (fnum.eq.19) fname1='r_gribawips19'
      if (fnum.eq.20) fname1='r_gribawips20'
      if (fnum.eq.21) fname1='r_gribawips21'
      if (fnum.eq.22) fname1='r_gribawips22'
      if (fnum.eq.23) fname1='r_gribawips23'
      if (fnum.eq.24) fname1='r_gribawips24'
      if (fnum.eq.25) fname1='r_gribawips25'
      if (fnum.eq.26) fname1='r_gribawips26'
      if (fnum.eq.27) fname1='r_gribawips27'
      if (fnum.eq.28) fname1='r_gribawips28'
      if (fnum.eq.29) fname1='r_gribawips29'
      if (fnum.eq.30) fname1='r_gribawips30'

      if (fnum.eq.1)  fname2='r_gribawips01.i'
      if (fnum.eq.2)  fname2='r_gribawips02.i'
      if (fnum.eq.3)  fname2='r_gribawips03.i'
      if (fnum.eq.4)  fname2='r_gribawips04.i'
      if (fnum.eq.5)  fname2='r_gribawips05.i'
      if (fnum.eq.6)  fname2='r_gribawips06.i'
      if (fnum.eq.7)  fname2='r_gribawips07.i'
      if (fnum.eq.8)  fname2='r_gribawips08.i'
      if (fnum.eq.9)  fname2='r_gribawips09.i'
      if (fnum.eq.10) fname2='r_gribawips10.i'
      if (fnum.eq.11) fname2='r_gribawips11.i'
      if (fnum.eq.12) fname2='r_gribawips12.i'
      if (fnum.eq.13) fname2='r_gribawips13.i'
      if (fnum.eq.14) fname2='r_gribawips14.i'
      if (fnum.eq.15) fname2='r_gribawips15.i'
      if (fnum.eq.16) fname2='r_gribawips16.i'
      if (fnum.eq.17) fname2='r_gribawips17.i'
      if (fnum.eq.18) fname2='r_gribawips18.i'
      if (fnum.eq.19) fname2='r_gribawips19.i'
      if (fnum.eq.20) fname2='r_gribawips20.i'
      if (fnum.eq.21) fname2='r_gribawips21.i'
      if (fnum.eq.22) fname2='r_gribawips22.i'
      if (fnum.eq.23) fname2='r_gribawips23.i'
      if (fnum.eq.24) fname2='r_gribawips24.i'
      if (fnum.eq.25) fname2='r_gribawips25.i'
      if (fnum.eq.26) fname2='r_gribawips26.i'
      if (fnum.eq.27) fname2='r_gribawips27.i'
      if (fnum.eq.28) fname2='r_gribawips28.i'
      if (fnum.eq.29) fname2='r_gribawips29.i'
      if (fnum.eq.30) fname2='r_gribawips30.i'
      return
      end


      subroutine CPCOEF(x,y,z,r,n,ictl)
C
C     This program will calculate the precipitation calibration
C     coefficence by using statistical distributions
C
C     Program: IBM-ASP     By: Yuejian Zhu ( 03/22/2001 )
C              Modified    By: Yuejian Zhu ( 03/21/2004 )
C              Modified    By: Jun Du      ( 10/15/2010 )
C
C     Input:    x(n)-vector for preciptation threshold amount
C               y(n)-vector for observation numbers at x(n)
C               z(n)-vector for forecasting numbers at x(n)
C               n   -length of the vector
C               ictl-to control the interpolation bases
C                    1 - standard
C                    2 - logrithm
C     Output:   r(n)-coefficents/ratio of each threshold amount x(n)
C               which will apply to particular precipitation amount
C               multiply the ratio at this point (linear interpolation)
C
      dimension x(n),y(n),z(n),r(n),a(n),b(n)
      dimension rnti(n-2),rnto(n-2)

C     Safty check of x-axis (dimension y)
C      Repeating to confirm the x-axis is ascending
      do i = 1, n-1
        if (y(i).ge.y(i+1)) then
         y(i+1) = y(i+1) + 1.0
        endif
      enddo

      if ( ictl.eq.1) then

C
C   input: y(n) as abscissas (x-axis)
C   input: x(n) as ordinates (y-axis)
C   input: z(n) as request abscissas values ( x-axis )
C   output: b(n) as cooresponding values of z(n)
C           similar to the thread amount, but shifted

       call qtpint(y,x,n,2,z,r,n) !Jim Purser's routine
c      call stpint(y,x,n,2,z,r,n,aux,naux) !essl lib
c      call stpint(prob_obs,thresh_old,n,2,prob_fcst,thresh_new,n,aux,naux)

       write(*,991) (y(i)/100.0,i=1,n)
       write(*,993) (x(i),i=1,n)
       write(*,992) (z(i)/100.0,i=1,n)
       write(*,994) (r(i),i=1,n)
       write(*,995) (r(i)/x(i),i=1,n)

c      do i = 1, n-2
c       rnti(i) = x(n-i)
c       rnto(i) = r(n-i)/x(n-i)
c      enddo
c      call qtpint(rnti,rnto,n-2,2,x,r,n) !Jim Purser's routine
cc      call stpint(rnti,rnto,n-2,2,x,r,n,aux,naux) !essl lib
c      write(*,993) (r(i),i=1,n)
c      if (r(n).lt.0.0.and.r(n-1).gt.0.0) then
c       r(n) = r(n-1)*r(n-2)
c      endif
c      write(*,993) (r(i),i=1,n)

      else
C
C  Tested both of log and log10, log is better
C
C   input: y(n) as abscissas (x-axis)
C   input: x(n) as ordinates (y-axis) -- using logrithem a(n) instead
C   input: z(n) as request abscissas values ( x-axis )
C   output: b(n) as cooresponding values of z(n)
C           similar to the thread amount, but shifted

       do i = 1, n
        a(i) = alog(x(i))
c       a(i) = log10(x(i))
       enddo

       call qtpint(y,a,n,2,z,b,n) !Jim Purser's routine
c      call stpint(y,a,n,2,z,b,n,aux,naux) !essl lib

       do i = 1, n
        r(i) = exp(b(i))
c       r(i) = exp(b(i)*log(10.0))
        if (r(i).gt.100.0.or.r(i).le.100.0) then
        else
         print *, "i=",i,"  r(i)=",r(i)," problem, use default"
         r(i) = x(i)
        endif
       enddo

       write(*,991) (y(i)/100.0,i=1,n)
       write(*,993) (x(i),i=1,n)
       write(*,992) (z(i)/100.0,i=1,n)
       write(*,994) (r(i),i=1,n)
       write(*,995) (r(i)/x(i),i=1,n)

c      do i = 1, n-2
c       rnti(i) = x(n-i)
c       rnto(i) = exp(b(n-i))/x(n-i)
c       rnto(i) = exp(b(n-i)*log(10.0))/x(n-i)
c      enddo
c      call qtpint(rnti,rnto,n-2,2,x,r,n) !Jim Purser's routine
cc      call stpint(rnti,rnto,n-2,2,x,r,n,aux,naux) !essl lib
c      write(*,993) (r(i),i=1,n)
c      if (r(n).lt.0.0.and.r(n-1).gt.0.0) then
c       r(n) = r(n-1)*r(n-2)
c      endif
c      write(*,993) (r(i),i=1,n)
      endif

      do i = 1, n
       r(i) = r(i)/x(i)    !convert to coefficient
      enddo

 991  format ('input = ',10f7.2,' OBS/100')
 992  format ('input = ',10f7.2,' FST/100')
 993  format ('input = ',10f7.2,' thrd    ')
 994  format ('output= ',10f7.2,' thrd_FST')
 995  format ('ratio = ',10f7.2,' tFST/thrd')
      return
      end

C=============================================================================
          subroutine qtpint(x,y,n,nint,t,s,m)
C=============================================================================
C
C * Jim Purser    *
C * 31st Aug 2012 *
C
C Quick-fix substitute for old STPINT routine, when nint=2 only.
C Linearly interpolate from n source pairs (x,y) to m targets (t,s)
C=============================================================================
        integer n,m,nint
        dimension x(n),y(n)
        dimension t(m),s(m)
C-----------------------------------------------------------------------------
        integer i,ip,j
        real    tj
C=============================================================================
        if(nint/=2)stop 'In qtpint; nint must be equal to 2'
        i=1
        do j=1,m
           tj=t(j)
           do ip=i+1,n-1
              if(x(ip)>=tj)exit
           enddo
           i=ip-1
           s(j)=( (x(ip)-tj)*y(i)+(tj-x(i))*y(ip) )/(x(ip)-x(i))
        enddo
        return
        end
