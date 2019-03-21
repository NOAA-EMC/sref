      program fcsterror

C$$  MAIN PROGRAM DOCUMENTATION BLOCK
C      
C MAIN PROGRAM: fcsterror
C PRGMMR: Jun Du          ORG: NP21       DATE: 2004-04-09
C
C ABSTRACT: 
c     1. decode grib data 
c     2. calculate any subgroup ensemble mean, spread and fcst error
c     3. write out subgroup ensemble mean, spread and fcst error for 
c        next program to use
c
c PROGRAM HISTORY LOG:
c Apr 09, 2004, J. Du: initial program
c Sep 09, 2004, J. Du: added grib-output capability
c Feb 27, 2006, J. Du: expanded from 15 mem to 21 mem and output subgroups' 
c                      mean and spread in grib format too
c May 21, 2007, J. Du: write bias etc fiedls of all past days in one single
c                      file for suitability of operational implementation
c Jul 28, 2008, J. Du: increased nmm and arw membership from 3 to 5 and reduced
c                      eta membership from 10 to 6 
c Jun 08, 2009, J. Du: added a flag to indicate if the input fcst and analysis 
c                      files are ok, if not ok (dataOK=1) the fcst error is set 
c                      to zero
c Nov 02, 2011, J. Du: modified to the new 2012 NEMS-based SREF (v6.0.0) by 
c                      reducing grouping from 4 (Eta, RSM, NMM and ARW) to
c                      3 (NMMB, NMM and ARW)
c Mar 10, 2015, J. Du: modified to the new 2015 2-model 26-member SREF (v7.0.0)
c 
c INPUT:
c OUTPUT:
c SUBPROGRAMS CALLED: 
c   mean_spread
c   correlation
c   grange
c   getname
C ATTRIBUTES:
c   language: Fortran 77
c   machine: IBM-SP
C
C$$$
 
      parameter(mem=26,nfile=mem+1,ntime=29,nvmax=150,nday=1)
      parameter(nSREF=2) !number of subgroups you want to divide
C nfile - 1-mem:ens mems; mem+1:NDAS
      parameter(jf=185*129,lon=185,lat=129,ngrid=lon*lat)  !for Grid 212

      integer mem_sta(nSREF),mem_end(nSREF)
      dimension f(jf),data(ngrid,nfile,nvmax),ver(jf,nvmax),
     &fcst(nfile),fmean(ngrid),spread(ngrid),error(ngrid),
     &abserr(ngrid),abserr0(ngrid),
     &fmean0(ngrid),spread0(ngrid),error0(ngrid),x(ngrid),y(ngrid),
     &pastfcst(ngrid,nday),pastsprd(ngrid,nday),pasterro(ngrid,nday),
     &corr_fcst(nday),corr_erro(nday),corr_sprd(nday),corr_spsk(nday)
      dimension jpds5(nvmax),jpds6(nvmax),jpds7(nvmax),kkpds(25,nvmax)
      dimension kftime(nfile,ntime),jpds(25),jgds(22),kpds(25),kgds(22)
      integer yy1,mm1,dd1,hh1,dataOK
      integer yy3(ntime),mm3(ntime),dd3(ntime),hh3(ntime)
      integer pyy(nday),pmm(nday),pdd(nday)
      logical lb(jf)
      character*4  var(nvmax)
      character*4  nam
      character*2  date
      character*2  sys
      character*13  fname1
      character*15  fname2
      character*18 fout1,fout2,fout3
      character*16 finp1,finp3,finp5
      character*18 finp2,finp4,finp6
      namelist/namin/yy1,mm1,dd1,hh1,nvar0,nvar,dataOK
      data (mem_sta(i),i=1,nSREF) /01,14/
      data (mem_end(i),i=1,nSREF) /13,26/

c Passing over date information
      read(5,namin,end=1000)
1000  continue
      print*,yy1,mm1,dd1,hh1,nvar0,nvar,dataOK

      read(7,*) (yy3(i),i=1,ntime),(mm3(i),i=1,ntime),
     &(dd3(i),i=1,ntime),(hh3(i),i=1,ntime)
      print*,(yy3(i),i=1,ntime),(mm3(i),i=1,ntime),
     &(dd3(i),i=1,ntime),(hh3(i),i=1,ntime)
      read(7,*) ((kftime(i,j),j=1,ntime),i=1,nfile)
      do i=1,nfile
       print *,(kftime(i,j),j=1,ntime)
      enddo
      rewind(7)

      read(8,*) (pyy(i),i=1,nday),(pmm(i),i=1,nday),
     &(pdd(i),i=1,nday)
      print*,(pyy(i),i=1,nday)
      print*,(pmm(i),i=1,nday)
      print*,(pdd(i),i=1,nday)

      do 10 nt=1,ntime

        print *,' '
        print *,'----------------------------------------------- '
        print *,'NEW TIME: nt = ',nt
        print *,' '

        j=0 
        data=0.0
        ver=0.0

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c part I: Decoding fcst and analysis
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c define variables to be decoded
       nunit=101
       open (nunit, file='variable.tbl', status='old')
       read (nunit,*)
       read (nunit,*)
       do nv=nvar0,nvar
        read(nunit,102) var(nv),jpds5(nv),jpds6(nv),jpds7(nv)
       enddo
       close(101)
102    format(4x,1A4,1x,3i8)

        do 5 nf=1,nfile
 	lugb=8        !if nvar>1
        lugi=9        !if nvar>1
c be careful of year and month when in the turn of time
        jpds=-1
        jgds=-1
C specify date information here:
       if(nf.le.mem) then
c       jpds(8)=yy1
c       jpds(9)=mm1
        jpds(10)=dd1
        jpds(11)=hh1
       endif
       if(nf.eq.(mem+1)) then
      read(7,*) (yy3(i),i=1,ntime),(mm3(i),i=1,ntime),
     &(dd3(i),i=1,ntime),(hh3(i),i=1,ntime)
      rewind(7)
c       jpds(8)=yy3(nt)
c       jpds(9)=mm3(nt)
        jpds(10)=dd3(nt)
        jpds(11)=hh3(nt)
       endif

        print *,' '
        print *,'----------'
        print *,'New File: nt = ',nt,' nf= ',nf
        print *,' '

       do 111 nv=nvar0,nvar
        jpds(5)=jpds5(nv)
        jpds(6)=jpds6(nv)
        jpds(7)=jpds7(nv)
       print*,jpds(5), ' ',jpds(6),' ',jpds(7)
        jpds(14)=kftime(nf,nt)
        jpds(15)=-1
      print*, jpds(10)
      print*, jpds(15)

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
 
        call grange(kf,lb,f,dmin,dmax) 
        print *,'immediately after call to getgb, iret= ',iret
     &,' j= ',j,' k= ',k,' nf= ',nf,' nt= ',nt
     &,' j5= ',jpds(5),' j6= ',jpds(6),' j7= ',jpds(7)
     &,' j8= ',jpds(8),' j9= ',jpds(9),' j10= ',jpds(10)
     &,' j11= ',jpds(11),' j14= ',jpds(14),' j15= ',jpds(15)
      print '(i4,2x,9i5,i8,2g12.4)',
     &k,(kpds(i),i=5,11),kpds(14),kpds(15),kf,dmin,dmax

      if (iret.eq.0) then
       do ipt=1,lon*lat

       data(ipt,nf,nv) = f(ipt)
       if(nt.eq.1) then
        do i=1,25
         kkpds(i,nv)=kpds(i)
        enddo
       endif

c against regional analysis
       if(nf.eq.(mem+1)) ver(ipt,nv)=f(ipt)
       if(nf.eq.(mem+1).and.ipt.eq.1) print*,'***REG ANL used***'

       enddo
      else
       print *, "Something wrong in fcst and analysis decoding!!!!"
c      exit   !make sure all neccessary data exist
      endif

111     continue    !diff variables (nv)
5	continue    !files (mems etc.)

	print*, 'fcsts & analy decoding is successfully done'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c part II: read in past's ensemble mean fcst, spread and fcst error
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        write(date,'(i2.2)') (nt-1)*3+3
       do nens=1,nSREF  !how many sub-SREFs
        write(sys,'(i2.2)') nens
c readin files
        finp1='yesterdayerro.' // sys
        finp2='yesterdayerro.' // sys // '.i'
        finp3='yesterdayfcst.' // sys
        finp4='yesterdayfcst.' // sys // '.i'
        finp5='yesterdaysprd.' // sys
        finp6='yesterdaysprd.' // sys // '.i'
        lugb_err=40
        lugi_err=41
         call baopenr(lugb_err,finp1,ierr)
         call baopenr(lugi_err,finp2,ierr)
        lugb_fcst=42
        lugi_fcst=43
         call baopenr(lugb_fcst,finp3,ierr)
         call baopenr(lugi_fcst,finp4,ierr)
        lugb_sprd=44
        lugi_sprd=45
         call baopenr(lugb_sprd,finp5,ierr)
         call baopenr(lugi_sprd,finp6,ierr)

c writeout files
        fout1='grb_pasterro' // sys // '.f' // date
        fout2='grb_pastfcst' // sys // '.f' // date
        fout3='grb_pastsprd' // sys // '.f' // date
        call baopen(71,fout1,ierr)
        call baopen(72,fout2,ierr)
        call baopen(73,fout3,ierr)
c     open(90,file='asc_pasterro' // sys // '.f' // date, 
c    &status='unknown')
c     open(91,file='asc_pastfcst' // sys // '.f' // date,
c    &status='unknown')
c     open(92,file='asc_pastsprd' // sys // '.f' // date,
c    &status='unknown')
      open(93,file='asc_pastcorr' // sys // '.f' // date,
     &status='unknown')

        do nv=nvar0,nvar  !for selected variables
         write(nam,'(a4)') var(nv)    !variable name

         jpds=-1
         jgds=-1
         jpds(5)=jpds5(nv)
         jpds(6)=jpds6(nv)
         jpds(7)=jpds7(nv)
         jpds(12)=0
         jpds(14)=kftime(1,nt)   !same as ensemble member 1
c        jpds(15)=-1
         jpds(15)=0
         jpds(21)=21

c read pasterr, pastfcst, pastsprd and pastbias
        do npast=1,nday-1
          jpds(8)=pyy(npast)-(pyy(npast)/100)*100
          jpds(9)=pmm(npast)
          jpds(10)=pdd(npast)
          jpds(11)=hh1
c ----------error-------------------
          call getgb(lugb_err,lugi_err,jf,j,jpds,jgds,kf,k,kpds,
     &               kgds,lb,f,iret)
          call grange(kf,lb,f,dmin,dmax)
          if (nens.eq.1.and.npast.eq.1) then
          print *,'read pasterr, iret=, nt= ',iret,nt,dmin,dmax
          endif
c         print *,'read pasterr, iret=, nt= ',iret,nt,dmin,dmax
c         print *,'reading pastday error, iret= ',iret
c    &,' j= ',j,' k= ',k,' nsub= ',nens,' nt= ',nt
c    &,' j5= ',jpds(5),' j6= ',jpds(6),' j7= ',jpds(7)
c    &,' j8= ',jpds(8),' j9= ',jpds(9),' j10= ',jpds(10)
c    &,' j11= ',jpds(11),' j14= ',jpds(14),' j15= ',jpds(15)
c         print '(i4,2x,9i5,i8,2g12.4)',
c    &k,(kpds(i),i=5,11),kpds(14),kpds(15),kf,dmin,dmax
          if (iret.eq.0) then
           do ipt=1,ngrid
            pasterro(ipt,npast) = f(ipt)
           enddo
          else
           print *,"Something wrong in reading pastday error!"
          endif
c----------------fcst---------------------
          call getgb(lugb_fcst,lugi_fcst,jf,j,jpds,jgds,kf,k,kpds,
     &               kgds,lb,f,iret)
          call grange(kf,lb,f,dmin,dmax)
          if (nens.eq.1.and.npast.eq.1) then
          print *,'read pastfcst, iret=, nt= ',iret,nt,dmin,dmax
          endif
c         print *,'read pastfcst, iret=, nt= ',iret,nt,dmin,dmax
          if (iret.eq.0) then
           do ipt=1,ngrid
            pastfcst(ipt,npast) = f(ipt)
           enddo
          else
           print *,"Something wrong in reading pastday fcst!"
          endif
c----------------sprd---------------------
          call getgb(lugb_sprd,lugi_sprd,jf,j,jpds,jgds,kf,k,kpds,
     &               kgds,lb,f,iret)
          call grange(kf,lb,f,dmin,dmax)
          if (nens.eq.1.and.npast.eq.1) then
          print *,'read pastsprd, iret=, nt= ',iret,nt,dmin,dmax
          endif
c         print *,'read pastsprd, iret=, nt= ',iret,nt,dmin,dmax
          if (iret.eq.0) then
           do ipt=1,ngrid
            pastsprd(ipt,npast) = f(ipt)
           enddo
          else
           print *,"Something wrong in reading pastday spread!"
          endif

        enddo  !npast

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c part III:calculate and write out ensemble mean, spread and fcst error
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        do i=1,25
         kpds(i)=kkpds(i,nv)  !force pds(4) to 128 not 192 to avoid bitmap undefined error
c        kpds(4)=0
        enddo
        kpds(21)=21
        kpds(8)=yy1-(yy1/100)*100
        kpds(9)=mm1
        kpds(10)=dd1
        kpds(11)=hh1
        kpds(12)=0
        kpds(14)=(nt-1)*3+3
        kpds(15)=0
        kpds(5)=jpds5(nv)
        kpds(6)=jpds6(nv)
        kpds(7)=jpds7(nv)
c Output precision
        if(kpds(5).eq.51) then
         kpds(22)=8
        else
         kpds(22)=2
        endif

c calculate mean, spread and fcst error over a subgroup
        do ipt=1,lon*lat
         do nf=mem_sta(nens),mem_end(nens) 
          fcst(nf)=data(ipt,nf,nv)
         enddo
      call mean_spread(nfile,mem_sta(nens),mem_end(nens),fcst,ave,std)
         fmean0(ipt)=ave
         if(dataOK.eq.0) then
          error0(ipt)=ave-ver(ipt,nv)
          abserr0(ipt)=abs(ave-ver(ipt,nv))

C testing
c         if(kpds(5).eq.11.and.kpds(6).eq.105) then
c          if(abserr0(ipt).gt.500) print*,'error in T2m',
c    c' nt=',nt,' nv=',nv,' ipt=',ipt,' group=',nens,
c    c' ave=',ave,' ver=',ver(ipt,nv),'abserr=',abserr0(ipt),
c    c' mem_sta and mem_end',mem_sta(nens),mem_end(nens)
c         endif

         else
          error0(ipt)=0.0
          abserr0(ipt)=0.0
          if (ipt.eq.1) 
      &print *,"error is set to zero due to unavailable input files!"
         endif
         spread0(ipt)=std
        enddo
        print*,'I am here at 1'
        call correlation(fmean0,fmean0,r,lat,lon)
         corr_fcst(1)=r
        call correlation(error0,error0,r,lat,lon)
         corr_erro(1)=r
        call correlation(spread0,spread0,r,lat,lon)
         corr_sprd(1)=r
        call correlation(abserr0,spread0,r,lat,lon)
         corr_spsk(1)=r

c write out mean, spread and fcst error over a subgroup
c       write(90,*) nam
c       write(90,*) (error(ipt),ipt=1,lon*lat)
c       write(91,*) nam
c       write(91,*) (fmean(ipt),ipt=1,lon*lat)
c       write(92,*) nam
c       write(92,*) (spread(ipt),ipt=1,lon*lat)

        print*,'I am here at 2'
        if (nv.ge.2) then   !because the first record is misteriously bad
        call putgb(71,jf,kpds,kgds,lb,error0,iret)
        call putgb(72,jf,kpds,kgds,lb,fmean0,iret)
        call putgb(73,jf,kpds,kgds,lb,spread0,iret)
        endif

        do i=1,nday-1
         do ipt=1,ngrid
          error(ipt)=pasterro(ipt,i)
          abserr(ipt)=abs(pasterro(ipt,i))
          fmean(ipt)=pastfcst(ipt,i)
          spread(ipt)=pastsprd(ipt,i)
         enddo
         kpds(8)=pyy(i)-(pyy(i)/100)*100
         kpds(9)=pmm(i)
         kpds(10)=pdd(i)
        if (nv.ge.2) then   !because the first record is misteriously bad
         call putgb(71,jf,kpds,kgds,lb,error,iret)
         call putgb(72,jf,kpds,kgds,lb,fmean,iret)
         call putgb(73,jf,kpds,kgds,lb,spread,iret)
         call correlation(fmean0,fmean,r,lat,lon)
        endif
          corr_fcst(i+1)=r
         call correlation(error0,error,r,lat,lon)
          corr_erro(i+1)=r
         call correlation(spread0,spread,r,lat,lon)
          corr_sprd(i+1)=r
         call correlation(abserr,spread,r,lat,lon)
          corr_spsk(i+1)=r

        enddo 

      if (iret.eq.0) then
       print *, "successful in grib writing at ntime=",(nt-1)*3+3
      else
       print *, "error in grib writing at ntime=",(nt-1)*3+3
      endif
c      print *,'immediately after call to putgb, iret= ',iret
c    &,' nens= ',nens,' nt= ',(nt-1)*3+3
cc   &,' lb=',lb
c    &,' k1= ',kpds(1),' k2= ',kpds(2),' k3= ',kpds(3)
c    &,' k4= ',kpds(4)
c    &,' k5= ',kpds(5),' k6= ',kpds(6),' k7= ',kpds(7)
c    &,' k8= ',kpds(8),' k9= ',kpds(9),' k10= ',kpds(10)
c    &,' k11= ',kpds(11),' k12= ',kpds(12),' k13= ',kpds(13)
c    &,' k14= ',kpds(14),' k15= ',kpds(15),' k16= ',kpds(16)
c    &,' k17= ',kpds(17),' k18= ',kpds(18),' k19= ',kpds(19)
c    &,' k20= ',kpds(20),' k21= ',kpds(21),' k22= ',kpds(22)
c    &,' k23= ',kpds(23),' k24= ',kpds(24),' k25= ',kpds(25)

c write out correlations among fcst, error, spread, sprd_skill, and bias for monitering
        if (nv.ge.2) then   !because the first record is misteriously bad
        write(93,104) nam
        write(93,103) (corr_fcst(i),i=1,nday)
        write(93,103) (corr_erro(i),i=1,nday)
        write(93,103) (corr_sprd(i),i=1,nday)
        write(93,103) (corr_spsk(i),i=1,nday)
        write(93,*)
        endif
103     format(1x,20f5.1)
104     format(1x,1a4,1x,
     &41h(corr in fcst, error, spread, sprd_skill))

        enddo    !variables
        print*,'I am here at 3'
        if (nv.ge.2) then   !because the first record is misteriously bad
        call baclose(lugb_err,ierr)
        call baclose(lugi_err,ierr)
        call baclose(lugb_fcst,ierr)
        call baclose(lugi_fcst,ierr)
        call baclose(lugb_sprd,ierr)
        call baclose(lugi_sprd,ierr)
        call baclose(71,ierr)
        call baclose(72,ierr)
        call baclose(73,ierr)
        close(90)
        close(91)
        close(92)
        close(93)
        endif
       enddo  !nSREF

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

C******************************************************
      subroutine mean_spread(nmax,mem_sta,mem_end,f,ave,std)
c.............................................................
c  Purpose: computing ensemble mean and spread among x number of members
c           at a given grid point (could be a full set or a subset which
c           is especially useful for multi-model ensemble system)
c  Input:
c	 nmax=maximum number of ensemble members allowed
c	 mem_sta=the number of the first member in a sub-ensemble set
c	 mem_end=the number of the last member in a sub-ensemble set
c	 f=ensemble forecasts data
c  Output:
c	 ave=ensemble mean
c	 std=ensemble spread
c  
c  Author: Jun Du, April, 2004
c.............................................................
      dimension f(nmax)
      real ave,std
      ave=0.
      std=0.
      do i=mem_sta,mem_end
       ave=ave+f(i)/float(mem_end-mem_sta+1)
      enddo
      do i=mem_sta,mem_end
       std=std+(f(i)-ave)**2.0
      enddo
      std=sqrt(std/float(mem_end-mem_sta+1))

      if(abs(ave).gt.1.0E+30) then
       print*,'f=',(f(i),i=mem_sta,mem_end)
      endif
      return
      end

C**********************************************
        subroutine correlation(x,y,r,lat,lon)
        dimension x(lon*lat),y(lon*lat)
        r=0.0
        r1=0.0
        r2=0.0
        xm=0.0
        ym=0.0
C calculate mean values over the selected region
        do ipt=1,lat*lon
         xm=xm+x(ipt)/(lat*lon)
         ym=ym+y(ipt)/(lat*lon)
        enddo
c calculate anamoally correlation coefficient
         do ipt=1,lat*lon
        r1=r1+(x(ipt)-xm)**2
        r2=r2+(y(ipt)-ym)**2
        r=r+(x(ipt)-xm)*(y(ipt)-ym)
         enddo
        if(r1.eq.0.0.and.r2.eq.0.0) r=100.0
        if(r1.eq.0.0.and.r2.ne.0.0) r1=0.01
        if(r2.eq.0.0.and.r1.ne.0.0) r2=0.01
        if(r1.ne.0.0.and.r2.ne.0.0)
     &  r=r/sqrt(r1)/sqrt(r2)*100.
        return
        end

