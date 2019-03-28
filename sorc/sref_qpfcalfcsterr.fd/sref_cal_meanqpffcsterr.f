      program calprecipfcsterror

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     1. Read in areal mask data for selected regions such as CONUS, 
c        EAST, WEST ...
c     2. Decode grib data
c     3. Calculate accumulated PDF for both fcst and obs over a region
c     4. Find the difference in PDF between fcst and obs
c     5. Find forecast error over various categories
c 
c Dependency: threshold (ncat), ensemble size (nens), fcst length (ntime)
c             and grid definition
c data ((kftime(i,j),j=1,ntime),i=1,nfile) 
c
c Log:
c 11/05/2009, Jun Du, initial coding
c 09/05/2010, Jun Du, refined the code
c 11/04/2011, Jun Du, a version for ensemble mean forecast or one member
c 03/28/2019, Jun Du, set fcsts to obs if past fcst data unavailable
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      parameter(nens=1,ntime=29)
      parameter(nfile=nens+1,ncat=6,nvar0=1,nvar=1)
      parameter(jf=185*129,lon=185,lat=129)   !grib212
c     parameter(jf=512*256,lon=185,lat=129)   !grib212
      parameter(numreg=30)
c 1-nens: fcsts; nens+1: obs.

      dimension f(jf),data(lon*lat,nfile,nvar),thresh(ncat),
     &mask(lon*lat),mask_2d(lon,lat),Pr(ncat+1,nfile,nvar),
     &kftime(nfile,ntime)
      dimension jpds(25),jgds(25),kpds(25),kgds(25)
      integer yy1,mm1,dd1,hh1,dataOK
      integer yy2(ntime),mm2(ntime),dd2(ntime),hh2(ntime)
      logical lb(jf)
      character*3  mdl
      character*5  reg
      character*13  fname1
      character*15  fname2
      character fmask(numreg)*3
      data fmask/ 'RFC','ATC','WCA','ECA','NAK','SAK','HWI','NPO',
     & 'SPO','NWC','SWC','NMT','SMT','NFR','SFR','NPL','SPL',
     & 'NMW','SMW','APL','NEC','SEC','NAO','SAO','PRI','MEX',
     & 'GLF','CAR','CAM','NSA'/
      namelist/namin/yy1,mm1,dd1,hh1,mdl,reg,dataOK

c     data thresh/76.2,63.5,50.8,38.1,25.4,19.1,12.7,6.35,2.54,.254/ !in mm
c     data thresh/ 3.0, 2.5, 2.0, 1.5, 1.0,0.75, 0.5,0.25, 0.1,0.01/ !in in
c     data thresh/50.8,44.45,38.1,31.75,25.4,19.05,12.7,6.35,2.54,.254/ !in mm
c     data thresh/ 2.0, 1.75, 1.5, 1.25, 1.0, 0.75, 0.5,0.25, 0.1,0.01/ !in in
c     data thresh/25.4,19.05,12.7,6.35,2.54,.254/  !in mm
      data thresh/ 1.0, 0.75, 0.5,0.25, 0.1,0.01/  !in in
      do i=1,ncat
       thresh(i)=thresh(i)*25.4  !convert to mm if used inch
      enddo

c Passing over date information
      read(5,namin,end=1000)
1000  continue
      print*,yy1,mm1,dd1,hh1
      print*,mdl,reg,dataOK

      open(55,file='mask.data',status='unknown')
      open(90,file='data.out',status='unknown')
      OPEN(888,FILE='regmask',FORM='UNFORMATTED')


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Part I: Read in MASK data for valid precipitation regions
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C read in regional masks for 212 grid
      READ(888) (mask(kk),kk=1,lat*lon)
      point=0.
      kk=0
      do j=1,lat
       do i=1,lon
        kk=kk+1

C  the 13 CONUS regions (see fmask above)
        if(reg.eq.'CONUS') then
         if (mask(kk).ge.10.and.mask(kk).le.22) then
          point=point+1.
          mask(kk)=1
         else
          mask(kk)=0
         endif
         mask_2d(i,j)=mask(kk)
        endif
C  the 7 WEST regions (see fmask above)
        if(reg.eq.'WESTR') then
         if (mask(kk).ge.10.and.mask(kk).le.14) then
          point=point+1.
          mask(kk)=1
         else
          mask(kk)=0
         endif
         mask_2d(i,j)=mask(kk)
        endif
C  the 6 EAST regions (see fmask above)
        if(reg.eq.'EASTR') then
         if (mask(kk).ge.15.and.mask(kk).le.22) then
          point=point+1.
          mask(kk)=1
         else
          mask(kk)=0
         endif
         mask_2d(i,j)=mask(kk)
        endif

       enddo  !i
      enddo   !j

      print*,'point_mask=',point,' out of ',lat*lon
      if(point.eq.0.0) print*,"There is a problem in defining a region!"
c make sure entire mask area digitized into "0 and 1"
      write(55,555) (mask(kk),kk=1,lat*lon)
      write(55,555)
c print out a portion of mask area map
      do j=lat,1,-1
       write(55,555) (mask_2d(i,j),i=53,147)
      enddo
555   format(185i1)
      print*, 'ok till here 1'


CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Part II: Decoding grib data
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c read in obs and fcst dates' info
      read(7,*) (yy2(i),i=1,ntime),(mm2(i),i=1,ntime),
     &(dd2(i),i=1,ntime),(hh2(i),i=1,ntime)
      print*,(yy2(i),i=1,ntime),(mm2(i),i=1,ntime),
     &(dd2(i),i=1,ntime),(hh2(i),i=1,ntime)
      read(7,*) ((kftime(i,j),j=1,ntime),i=1,nfile)
      do i=1,nfile
       print *,(kftime(i,j),j=1,ntime)
      enddo
      rewind(7)

      data=0.0
      do 70 nt=1,ntime


        print *,' '
        print *,'----------------------------------------------- '
        print *,'NEW TIME: nt = ',nt
        print *,' '

        j=0
        lugb=9        !if nvar=1
        lugi=39       !if nvar=1


       do 5  nf =1,nfile 

      if(dataOK.eq.1.and.nf.lt.nfile) then
       print*, 'Past forecast is not available and skip data decoding'
       goto 888
      endif

c       lugb=9        !if nvar>1
c       lugi=39       !if nvar>1
c be careful of year and month when in the turn of time
        jpds=-1
c       jgds=-1
C specify date information here:
       if(nf.le.nens) then
c       jpds(8)=yy1
c       jpds(9)=mm1
        jpds(10)=dd1
        jpds(11)=hh1
       endif
       if(nf.eq.nfile) then
c       jpds(8)=yy2(nt)
c       jpds(9)=mm2(nt)
        jpds(10)=dd2(nt)
        jpds(11)=hh2(nt)
       endif

        print *,' '
        print *,'----------'
        print *,'New File: nt = ',nt,' nf= ',nf
        print *,' '

c define variables to be verified
c nvar=1 (total ppt)
c nvar=2 (stratiform ppt)
c nvar=3 (convective ppt)
c......................
      do 111 nv=nvar0,nvar
         if(nv.eq.1)  then
          jpds(5)=61
          jpds(7)=0
         endif
         if(nv.eq.2)  then
          jpds(5)=62
          jpds(7)=0
         endif
         if(nv.eq.3)  then
          jpds(5)=63
          jpds(7)=0
         endif

        if(jpds(5).eq.61.or.jpds(5).eq.62.or.jpds(5).eq.63) then
         if(kftime(nf,nt).gt.3) then
          jpds(14)=kftime(nf,nt)-3
         else
          jpds(14)=0
         endif
         jpds(15)=kftime(nf,nt)
        else
         jpds(14)=kftime(nf,nt)
         jpds(15)=-1
        endif
       if(nf.eq.nfile) jpds(14)=0
       if(nf.eq.nfile) jpds(15)=3
      print*, jpds(10)
      print*, jpds(15)


      lugb=lugb+1
      lugi=lugi+1
      print *,'lugb= ',lugb,' lugi= ',lugi
      call getname(nf,fname1,fname2)
           print *,' fname1= ',fname1
           print *,' fname2= ',fname2
c     call baopen(lugb,fname1,ierr1)
      call baopenr(lugb,fname1,ierr1)
      if(ierr1.ne.0) then 
        write (*,*) 'fname1 error'
        stop 5555
      end if
c     call baopen(lugi,fname2,ierr2)
      call baopenr(lugi,fname2,ierr2)
      if(ierr2.ne.0) then
       write (*,*) 'fname2 error'
       stop 6666
      end if

      do k1=1,25
       write(*,*) k1, jpds(k1)
      end do

      call getgb(lugb,lugi,jf,j,jpds,jgds,
     &          kf,k,kpds,kgds,lb,f,iret)
      if (iret.ne.0) then
       write(*,*) 'iret=',iret
      else
       do k1=jf-100,jf
c       write(*,*) k1, f(k1)
       end do
      end if

      call baclose(lugb,ierr)
      call baclose(lugi,ierr)

      call grange(kf,lb,f,dmin,dmax)
      print *,'immediately after call to getgb, iret= ',iret
     &,' j= ',j,' k= ',k,' nf= ',nf,' nt= ',nt
     &,' j5= ',jpds(5),' j6= ',jpds(6),' j7= ',jpds(7)
     &,' j8= ',jpds(8),' j9= ',jpds(9),' j10= ',jpds(10)
     &,' j11= ',jpds(11),' j14= ',jpds(14),' j15= ',jpds(15)

c     print *,'immediately after call to getgb,jgds(1-25)='
c    &,(jgds(i),i=1,25)

C Note: Open the following both print statements will cause ensemble 
C mean file reading error in getgb! But open either one of them is ok.
C Very strange for unknown reason -- Jun Du, 2/7/2013
c     print '(i4,2x,9i5,i8,2f12.4)',
c    &k,(kpds(i),i=5,11),kpds(14),kpds(15),kf
      print*, 'dmin,dmax=',dmin,dmax
c Original:
c     print '(i4,2x,9i5,i8,2f12.4)',
c    &k,(kpds(i),i=5,11),kpds(14),kpds(15),kf,dmin,dmax

        if (iret.eq.0) then
         zero=0. !for testing
         do ipt=1,lon*lat
          data(ipt,nf,nv)=f(ipt)   !3-hr apcp
          if(nf.eq.11.and.f(ipt).gt.0.254) zero=zero+1  !checking to see if any precip in raw data
         enddo
        else
         print *, "Something wrong here!!!!"
        endif
        if(nf.eq.11) print*,'zero=',zero,' at nt=',nt

111   continue

888   continue   !jump to here if there is no input fcst data

5     continue
      print*, 'ok till here 2'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Part III: Calculate accumulated PDFs of fcsts and obs over a region
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      print*,'thresh=',(thresh(i),i=1,ncat)
      pr=0.

      do 222 nv=nvar0,nvar
       do 333 nf=1,nfile

        if(dataOK.eq.1.and.nf.lt.nfile) then !use obs for all if no fcst data
         do ipt=1,lat*lon
          data(ipt,nf,nv)=data(ipt,nfile,nv)
         enddo
        endif

        do icat=1,ncat
         do ipt=1,lat*lon
          if(mask(ipt).eq.1) then  !mask starts
           if(data(ipt,nf,nv).ge.thresh(icat)) then
            Pr(icat,nf,nv)=Pr(icat,nf,nv)+1.         !points
c           Pr(icat,nf,nv)=Pr(icat,nf,nv)+1./point   !pdf
           endif
          endif                    !mask ends
         enddo
        enddo

c       print*,(pr(i,nf,nv)*100,i=1,ncat)    !pdf
        print*,(pr(i,nf,nv),i=1,ncat)        !points
333     continue    !nf

222     continue    !variables--nv



CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Part IV: Output
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do nv=nvar0,nvar
       if(nt.eq.1.and.nv.eq.1) then
        write(90,900) mdl,reg,hh1,mm1,dd1,yy1
        write(90,901)
        write(90,902) (thresh(i),i=1,ncat)
       endif
        write(90,*)
        write(90,903) 3*nt
c pdf
        write(90,911) nens,point
        do nf=1,nfile
c        write(90,904) (pr(i,nf,nv)*100,i=1,ncat)    !pdf
c        write(90,904) (pr(i,nf,nv),i=1,ncat)        !points
c To avoid zero precip point by starting from at least 1 point
         write(90,904) (pr(i,nf,nv)+1,i=1,ncat)      !points
        enddo

      enddo   !nv

70      continue    !time--nt
        print*, 'ok till here 3'

900   format(a5,1x,a5,7h 3hapcp,1x,
     &16h Initiating from, 1x,i2,1hz,1x,i2,1h/,i2,1h/,i4)
901   format(19htotal precipitation)
902   format(10hthreshold=,10f7.3,3h mm)
903   format(i2,1hh)
904   format(2x,10f7.1)
905   format(2x,11f6.0)
906   format(2x,11f6.0)
907   format(2x,11f6.0)
908   format(2x,11f6.0)
909   format(2x,11f6.0)
910   format(2x,11f6.0)
911   format(2x,6hPDF(1-,i2,18h mem, obs); point=,f6.0)
912   format(2x,24hAreal Bias over a domain)
913   format(2x,25hAmount Bias over a domain)
914   format(2x,11hTotal Error)
915   format(2x,22f5.2)
916   format(2x,32hNormalized Amount Bias per point)
917   format(2x,15hAccumulated PDF)
918   format(2x,19hCategorilized error)


9999	END

C************************************************
      subroutine grange(n,ld,d,dmin,dmax)
      logical ld
      dimension ld(n),d(n)
      dmax=-1.e38
 
      do i=1,n
c       if(ld(i)) then
          dmin=min(dmin,d(i))
          dmax=max(dmax,d(i))
c       endif
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

      return
      end

