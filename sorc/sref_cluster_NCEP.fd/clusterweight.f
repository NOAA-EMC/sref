      program main_program

c#############################################################
c This program estimates each cluster's weight out of the whole 
c ensemble and then writes it out in grib1 format.
c
c log:
c 06/03/2012, Jun Du -- initial program for SREFv6.0.0 implementation
c
c#############################################################
 
      parameter(mem=1,iens=26,lead=87,interv=3,nvar0=1,nvar=1)
      parameter(ntime=lead/interv+1,nclus=6)
      
      parameter(jf=185*129,lon=185,lat=129,ngrid=lon*lat)  !for Grid 212
c     parameter(jf=512*256,lon=185,lat=129,ngrid=lon*lat)

      dimension f(jf),weight(ngrid),kftime(mem,ntime),
     &jpds(25),jgds(22),kpds(25),kgds(22)
      integer yy2(ntime),mm2(ntime),dd2(ntime),hh2(ntime)
      integer yy1,mm1,dd1,hh1,iar,inhr
      logical lb(jf)
      character*13  fname1
      character*15  fname2
      character*26 wgts
      character*2  date
      character*2  cyc
      character*1  order
      namelist/namin/yy1,mm1,dd1,hh1,iar

      open(50,file='sref_cluster_info',status='old')

      read(50,*)
      read(50,*) 

c Passing over date information
      read(5,namin,end=1000)
1000  continue
      print*,'yy1,mm1,dd1,hh1,iar=',yy1,mm1,dd1,hh1,iar

      read(7,*) (yy2(i),i=1,ntime),(mm2(i),i=1,ntime),
     &(dd2(i),i=1,ntime),(hh2(i),i=1,ntime)
      print*,(yy2(i),i=1,ntime),(mm2(i),i=1,ntime),
     &(dd2(i),i=1,ntime),(hh2(i),i=1,ntime)
      read(7,*) ((kftime(i,j),j=1,ntime),i=1,mem)

      do 10 nt=1,ntime

        print *,' '
        print *,'----------------------------------------------- '
        print *,'NEW TIME: nt = ',nt
        print *,' '

        j=0 
        data=0.0
 	lugb=9        !if nvar=1
        lugi=29       !if nvar=1

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c part I: To obtain pds, gds info
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        do 5 nf=1,mem
c	lugb=9        !if nvar>1
c       lugi=29       !if nvar>1
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

C for instantaneous fields (not accumulative quantity)
        jpds(14)=kftime(nf,nt)
        jpds(15)=-1

        print *,' '
        print *,'----------'
        print *,'New File: nt = ',nt,' nf= ',nf
        print *,' '

c       do 111 nv=nvar0,nvar
c iar=1 (500H)
         if(iar.eq.1)  then
          jpds(5)=7          !variable
          jpds(6)=100        !level type
          jpds(7)=500        !level
         endif
c iar=2 (slp)
         if(iar.eq.2)  then
          jpds(5)=2          !variable
          jpds(6)=102        !level type
          jpds(7)=0          !level
         endif

        print*, 'j5=',jpds(5)
        print*, 'j6=',jpds(6)
        print*, 'j7=',jpds(7)
        print*, 'j8=',jpds(8)
        print*, 'j9=',jpds(9)
        print*, 'j10=',jpds(10)
        print*, 'j11=',jpds(11)
        print*, 'j12=',jpds(12)
        print*, 'j13=',jpds(13)
        print*, 'j14=',jpds(14)
        print*, 'j15=',jpds(15)

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
     &,' j= ',j,' k= ',k,' ifile= ',ifile,' nt= ',nt
     &,' j1= ',jpds(1),' j2= ',jpds(2),' j3= ',jpds(3),' j4= ',jpds(4)
     &,' j5= ',jpds(5),' j6= ',jpds(6),' j7= ',jpds(7),' j8= ',jpds(8)
     &,' j9= ',jpds(9),' j10= ',jpds(10),' j11= ',jpds(11),' j12= ',jpds(12)
     &,' j13= ',jpds(13),' j14= ',jpds(14),' j15= ',jpds(15)
     &,' j16= ',jpds(16),' j17= ',jpds(17),' j18= ',jpds(18)
     &,' j19= ',jpds(19),' j20= ',jpds(20),' j21= ',jpds(21)
     &,' j22= ',jpds(22),' j23= ',jpds(23),' j24= ',jpds(24)
     &,' j25= ',jpds(25)
      print '(i4,2x,25i5,i8,2g12.4)',k,(kpds(i),i=1,25),kf,dmin,dmax

      if (iret.ne.0) then
       print *, "Something wrong in grib decoding!!!!"
      endif
   
c111     continue    !variable (nv)
5	continue    !files (nf)
	print*, 'ok till here 1 -- decoding'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c part II: Find membership in each of the six clusters and 
c write them out in grib format
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
         kpds(1)=7
         kpds(2)=130
         kpds(3)=212
c        kpds(4)=0  !bitmap Table
         kpds(5)=184
         kpds(6)=200
         kpds(7)=0
         kpds(8)=yy1-(yy1/100)*100
         kpds(9)=mm1
         kpds(10)=dd1
         kpds(11)=hh1
         kpds(12)=0
         kpds(14)=(nt-1)*3
         kpds(15)=0
         kpds(19)=129
         kpds(21)=21
         kpds(22)=4  !precision
c        kpds(23)=2  !if ensemble
         kpds(23)=0  !if not ensemble
         kpds(24)=128
         kpds(25)=-1


       read(50,*)
       read(50,100) inhr
       read(50,*)
       print*,'Processing fcst hour=',nt*3-3
       print*,'Read in fcst hour=',inhr
       print*,'The above two should be the same'
       write(date,'(i2.2)') (nt-1)*3

       do 20 ic=1,nclus
        read(50,200) nn
         do ipt=1,ngrid
c         weight(ipt)=100.0*(float(nn)/float(iens))  !in %
          weight(ipt)=float(nn)/float(iens)
         enddo
        do i=1,nn
         read(50,*)
        enddo
        write(order,'(i1.1)') ic
        write(cyc,'(i2.2)') hh1
        wgts='wgt.t'//cyc//'z.pgrb212.clus' // order // '.f' // date
        call baopen(90+ic,wgts,ierr)
        call putgb(90+ic,lat*lon,kpds,kgds,lb,weight,iret)
c       call putgb(90+ic,jf,kpds,kgds,lb,weight,iret)
        call baclose(90+ic,ierr)
20     continue
10	continue

100    format(1x,i2)
200    format(17x,i2)

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

      if (fnum.eq.1)  fname2='r_gribawips01.i'
      if (fnum.eq.2)  fname2='r_gribawips02.i'
      if (fnum.eq.3)  fname2='r_gribawips03.i'

      return
      end

