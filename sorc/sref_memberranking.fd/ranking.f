      program main_program

c#############################################################
c     1. decode grib data 
c     2. calculate distances between members in terms of magnitude and
c        pattern differences
c     3. find member performance ranking such as best, worst members
c        and assign corresponding weights to each members
c     4. find member clusters (incomplete yet)
c     5. output an information table
c                            ---- Jun Du, 12/09/2009
c log:
c 07/25/2006, Jun Du -- initial idea and programming
c 11/30/2009, Jun Du -- modified and performed a study (Du and Zhou, Oct 2011 issue of MWR)
c 11/05/2011, Jun Du -- further modified and implemented it to SREF operation
c 03/13/2015, Jun Du -- modified to 2-model 26-member SREF.v7.0.0
c
c#############################################################
 
      parameter(mem=26,lead=87,interv=3,nvar0=1,nvar=16)
c     parameter(mem=26,lead=87,interv=3,nvar0=1,nvar=13)
      parameter(ntime=lead/interv+1)
      
      parameter(jf=185*129,lon=185,lat=129,ngrid=lon*lat)  !for Grid 212
c     parameter(jf=512*256,lon=185,lat=129,ngrid=lon*lat)

      dimension f(jf),data(lon*lat,mem,nvar),ff(mem),
     &xmean(lon*lat,nvar),xmedian(lon*lat,nvar),x(lat*lon),
     &var(lon*lat,ntime,nvar),avevar(ntime,nvar),y(lat*lon),
     &dist((mem*mem-mem)/2,ntime,nvar),dist1(mem,ntime,nvar),
     &rho((mem*mem-mem)/2,ntime,nvar),rho1(mem,ntime,nvar),
     &sum_mag((mem*mem-mem)/2,ntime),sum_mag1(2,mem,ntime),
     &sum_pat((mem*mem-mem)/2,ntime),sum_pat1(2,mem,ntime),
     &fcst1(mem,ntime),fcst((mem*mem-mem)/2,ntime),tot(ntime),
     &index_mag((mem*mem-mem)/2,2,ntime),wgt_mag1(mem,ntime),
     &index_pat((mem*mem-mem)/2,2,ntime),tot_mag1(ntime),
     &tot_pat1(ntime),wgt_pat1(mem,ntime),weight(ngrid),
     &value_mag((mem*mem-mem)/2,ntime),value_mag1(2,mem,ntime),
     &value_pat((mem*mem-mem)/2,ntime),value_pat1(2,mem,ntime),
     &alinear(ntime,nvar),ave_alinear(ntime),ave_avevar(ntime),
     &kftime(mem,ntime),jpds(200),jgds(200),kpds(200),kgds(200)
c    &kftime(mem,ntime),jpds(25),jgds(22),kpds(25),kgds(22)
      integer lorder_mag((mem*mem-mem)/2,ntime),
     &lorder_mag1(2,mem,ntime),
     &lorder_pat((mem*mem-mem)/2,ntime),lorder_pat1(2,mem,ntime)
      integer yy,mm,dd,hh,meth
      logical lb(jf)
      character*13  fname1
      character*15  fname2
      character*18 wgts
      character*2  date
      character*2  order
      namelist/namin/yy,mm,dd,hh

      open(90,file='data.out',status='unknown')
      open(95,file='wgt.out',status='unknown')

c Passing over date information
      read(5,namin,end=1000)
1000  continue
      print*,yy,mm,dd,hh

      read(7,*) ((kftime(i,j),j=1,ntime),i=1,mem)
      do i=1,mem
       print *,(kftime(i,j),j=1,ntime)
      enddo

      do 10 nt=1,ntime

        print *,' '
        print *,'----------------------------------------------- '
        print *,'NEW TIME: nt = ',nt
        print *,' '

        j=0 
        data=0.0
c	lugb=9        !if nvar=1
c       lugi=39       !if nvar=1

CCCCCCCCCCCCCCCCCCCCCCCC
c part I: Decoding data:
CCCCCCCCCCCCCCCCCCCCCCCC
        do 5 nf=1,mem
 	lugb=9        !if nvar>1
        lugi=39       !if nvar>1
c be careful of year and month when in the turn of time
        jpds=-1
        jgds=-1
C specify date information here:
       if(nf.le.mem) then
c       jpds(8)=yy
c       jpds(9)=mm
        jpds(10)=dd
        jpds(11)=hh
       endif

C for instantaneous fields (not accumulative quantity)
        jpds(14)=kftime(nf,nt)
        jpds(15)=-1

        print *,' '
        print *,'----------'
        print *,'New File: nt = ',nt,' nf= ',nf
        print *,' '

c define parameters to be used in clustering
c SLP(2), H(7), T(11), u(33), v(34), RH(52), ppt(61)
c......................
      do 111 nv=nvar0,nvar
         if(nv.eq.1)  then
c nvar=1 (slp)
          jpds(5)=2          !variable
          jpds(6)=102        !level type
          jpds(7)=0          !level
         endif

         if(nv.eq.2)  then
c nvar=2 (H)
          jpds(5)=7          
          jpds(6)=100    
          jpds(7)=500          
         endif
         if(nv.eq.3)  then
c nvar=3 (H)
          jpds(5)=7          
          jpds(6)=100    
          jpds(7)=700          
         endif
         if(nv.eq.4)  then
c nvar=4 (H)
          jpds(5)=7          
          jpds(6)=100    
          jpds(7)=850          
         endif

         if(nv.eq.5)  then
c nvar=5 (T)
          jpds(5)=11         
          jpds(6)=100    
          jpds(7)=500         
         endif
         if(nv.eq.6)  then
c nvar=6 (T)
          jpds(5)=11         
          jpds(6)=100    
          jpds(7)=700         
         endif
         if(nv.eq.7)  then
c nvar=7 (T)
          jpds(5)=11         
          jpds(6)=100    
          jpds(7)=850         
         endif

         if(nv.eq.8)  then
c nvar=8 (U)
          jpds(5)=33          
          jpds(6)=100    
          jpds(7)=500        
         endif
         if(nv.eq.9)  then
c nvar=9 (U)
          jpds(5)=33
          jpds(6)=100    
          jpds(7)=700
         endif
         if(nv.eq.10)  then
c nvar=10 (U)
          jpds(5)=33
          jpds(6)=100    
          jpds(7)=850
         endif

         if(nv.eq.11)  then
c nvar=11 (V)
          jpds(5)=34          
          jpds(6)=100    
          jpds(7)=500        
         endif
         if(nv.eq.12)  then
c nvar=12 (V)
          jpds(5)=34          
          jpds(6)=100    
          jpds(7)=700        
         endif
         if(nv.eq.13)  then
c nvar=13 (V)
          jpds(5)=34          
          jpds(6)=100    
          jpds(7)=850        
         endif

         if(nv.eq.14)  then
c nvar=14 (RH)
          jpds(5)=52        
          jpds(6)=100    
          jpds(7)=500        
         endif
         if(nv.eq.15)  then
c nvar=15 (RH)
          jpds(5)=52        
          jpds(6)=100    
          jpds(7)=700        
         endif
         if(nv.eq.16)  then
c nvar=16 (RH)
          jpds(5)=52        
          jpds(6)=100    
          jpds(7)=850        
         endif

         if(nv.eq.17)  then
c nvar=17 (3h-apcp)
          jpds(7)=0        
          jpds(5)=61        
          jpds(14)=kftime(nf,nt)-3
          jpds(15)=kftime(nf,nt)
         endif

        print*, 'j5=',jpds(5)
        print*, 'j6=',jpds(6)
        print*, 'j7=',jpds(7)
        print*, 'j10=',jpds(10)
        print*, 'j11=',jpds(11)
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
     &,' j= ',j,' k= ',k,' nf= ',nf,' nt= ',nt
     &,' j5= ',jpds(5),' j6= ',jpds(6),' j7= ',jpds(7)
     &,' j8= ',jpds(8),' j9= ',jpds(9),' j10= ',jpds(10)
     &,' j11= ',jpds(11),' j14= ',jpds(14),' j15= ',jpds(15)
      print '(i4,2x,9i5,i8,2g12.4)',
     &k,(kpds(i),i=5,11),kpds(14),kpds(15),kf,dmin,dmax

      if (iret.eq.0) then
       do ipt=1,lon*lat
        data(ipt,nf,nv) = f(ipt)
       enddo
      else
       print *, "Something wrong in grib decoding!!!!"
      endif
   
c od -f fort.19
      if(nt.eq.11.and.nv.eq.1) write(19)(f(ipt),ipt=1,lon*lat)

111     continue    !variables (nv)
5	continue    !files (nf)
	print*, 'ok till here 1 -- decoding'

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c part II: find clusters, best member and extreme members
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      do meth=1,2   !meth in calculation

      do 222 nv=nvar0,nvar

      if(meth.eq.1) then  !need to be done only once
        alinear(nt,nv)=0.0
c Find Ensemble Mean and Median
        do 100 ipt=1,lon*lat
          do nf=1,mem      
           ff(nf)=data(ipt,nf,nv)
          enddo
          call mean_median(mem,ff,xm,xmed)
          xmean(ipt,nv)=xm        
          xmedian(ipt,nv)=xmed   

C Measuring linearity by differencing ctl mem and ens mean
          alinear(nt,nv)=alinear(nt,nv)+
     &abs(data(ipt,1,nv)-xmean(ipt,nv))/(lon*lat) 
100     continue

c Find standard deviation among members (ensemble spread)
      avevar(nt,nv)=0.0
      do 50 ipt=1,lon*lat
        var(ipt,nt,nv)=0.0
        do nf=1,mem 
         var(ipt,nt,nv)=var(ipt,nt,nv)+
     &   (data(ipt,nf,nv)-xmean(ipt,nv))**2
        enddo 
        var(ipt,nt,nv)=sqrt(var(ipt,nt,nv)/(mem-1))
        avevar(nt,nv)=avevar(nt,nv)+var(ipt,nt,nv)
50    continue
      avevar(nt,nv)=avevar(nt,nv)/(lon*lat) !ave over area
      endif

c Measuring difference by its magnitude: using absolute value difference as "distance"
c Find distance between member and ensemble mean or median for "best member" 
      do nf=1,mem 
       dist1(nf,nt,nv)=0.0
       do ipt=1,lon*lat
       if(meth.eq.1) dist1(nf,nt,nv)=dist1(nf,nt,nv)+
     &abs(data(ipt,nf,nv)-xmean(ipt,nv))      !using mean
       if(meth.eq.2) dist1(nf,nt,nv)=dist1(nf,nt,nv)+
     &abs(data(ipt,nf,nv)-xmedian(ipt,nv))    !using median
       enddo 
       dist1(nf,nt,nv)=dist1(nf,nt,nv)/(lon*lat)     !raw distance
       dist1(nf,nt,nv)=dist1(nf,nt,nv)/avevar(nt,nv) !normalized by spread
       sum_mag1(meth,nf,nt)=sum_mag1(meth,nf,nt)+
     &dist1(nf,nt,nv)/nvar  !sum over variables
      enddo
c Find distance matrix among all members for finding clusters and extreme members
      if(meth.eq.1) then
       m=0
      do i=1,mem-1   
       do j=i+1,mem
       m=m+1
       index_mag(m,1,nt)=i
       index_mag(m,2,nt)=j
       dist(m,nt,nv)=0.0
       do ipt=1,lon*lat
       dist(m,nt,nv)=dist(m,nt,nv)+abs(data(ipt,i,nv)-data(ipt,j,nv)) 
       enddo 
       dist(m,nt,nv)=dist(m,nt,nv)/(lon*lat)     !raw distance
       dist(m,nt,nv)=dist(m,nt,nv)/avevar(nt,nv) !normalized by spread
       sum_mag(m,nt)=sum_mag(m,nt)+dist(m,nt,nv)/nvar    !sum over variables
       enddo
      enddo
      endif

c Measuring difference by its pattern: using corr. coeff. as a "distance" measure
c Find corr. coeff. between member and ensemble mean or median for best member
      do nf=1,mem   
       do ipt=1,lat*lon
        x(ipt)=data(ipt,nf,nv)
        if(meth.eq.1) y(ipt)=xmean(ipt,nv)     !using mean
        if(meth.eq.2) y(ipt)=xmedian(ipt,nv)   !using median
c       if(meth.eq.2) y(ipt)=x(ipt)            !for code testing
       enddo 
       call correlation(x,y,r,lat,lon)
       rho1(nf,nt,nv)=r
       sum_pat1(meth,nf,nt)=sum_pat1(meth,nf,nt)+
     &rho1(nf,nt,nv)/nvar !sum over variables
      enddo
c Find corr. coeff. matrix among members for finding clusters and extreme members
      if(meth.eq.1) then
       m=0
      do i=1,mem-1   
       do j=i+1,mem
       m=m+1
       index_pat(m,1,nt)=i
       index_pat(m,2,nt)=j
       do ipt=1,lat*lon
        x(ipt)=data(ipt,i,nv)
        y(ipt)=data(ipt,j,nv)
       enddo 
       call correlation(x,y,r,lat,lon)
       rho(m,nt,nv)=r
       sum_pat(m,nt)=sum_pat(m,nt)+rho(m,nt,nv)/nvar !sum over variables
       enddo
      enddo
      endif

222     continue

c Sort out new multi-variate combined-distance in ascending order (magnitude)
      do i=1,mem
       fcst1(i,1)=i
       fcst1(i,2)=sum_mag1(meth,i,nt)
       fcst1(i,3)=0
      enddo
      call sortm(fcst1,mem,3,2)
      do i=1,mem
       lorder_mag1(meth,i,nt)=int(fcst1(i,3))
       value_mag1(meth,i,nt)=fcst1(i,2)
      enddo
      print*,'best member'
c     print*,(fcst1(i,1),i=1,mem)
c     print*,'distance:',(fcst1(i,2),i=1,mem)
      print*,'distance:',(value_mag1(meth,i,nt),i=1,mem)
      print*,'member IDs:',(int(fcst1(i,3)),i=1,mem)

      if(meth.eq.1) then
      print*,'m=',m
      print*,'(mem*mem-mem)/2=',(mem*mem-mem)/2
      do i=1,(mem*mem-mem)/2
       fcst(i,1)=i
       fcst(i,2)=sum_mag(i,nt)
       fcst(i,3)=0
      enddo
      call sortm(fcst,(mem*mem-mem)/2,3,2)
      do i=1,(mem*mem-mem)/2
       lorder_mag(i,nt)=int(fcst(i,3))
       value_mag(i,nt)=fcst(i,2)
      enddo
      print*,'extreme members'
c     print*,(fcst(i,1),i=1,(mem*mem-mem)/2)
c     print*,'distance:',(fcst(i,2),i=1,(mem*mem-mem)/2)
      print*,'distance:',(value_mag(i,nt),i=1,(mem*mem-mem)/2)
      print*,'member IDs:',(int(fcst(i,3)),i=1,(mem*mem-mem)/2)
      do i=1,(mem*mem-mem)/2
      print*,'(',index_mag(int(fcst(i,3)),1,nt),',',
     &index_mag(int(fcst(i,3)),2,nt),')'
      enddo
	print*, 'ok till here 2 -- magnitude distance calculation'
      endif

c Sort out new multi-variate combined-distance in ascending order (pattern)
      do i=1,mem
       fcst1(i,1)=i
       fcst1(i,2)=sum_pat1(meth,i,nt)
       fcst1(i,3)=0
      enddo
      call sortm(fcst1,mem,3,2)
      do i=1,mem
       lorder_pat1(meth,i,nt)=int(fcst1(i,3))
       value_pat1(meth,i,nt)=fcst1(i,2)
      enddo
      print*,'best member'
c     print*,(fcst1(i,1),i=1,mem)
c     print*,'distance:',(fcst1(i,2),i=1,mem)
      print*,'distance:',(value_pat1(meth,i,nt),i=1,mem)
      print*,'member IDs:',(int(fcst1(i,3)),i=1,mem)

      if(meth.eq.1) then
      print*,'m=',m
      print*,'(mem*mem-mem)/2=',(mem*mem-mem)/2
      do i=1,(mem*mem-mem)/2
       fcst(i,1)=i
       fcst(i,2)=sum_pat(i,nt)
       fcst(i,3)=0
      enddo
      call sortm(fcst,(mem*mem-mem)/2,3,2)
      do i=1,(mem*mem-mem)/2
       lorder_pat(i,nt)=int(fcst(i,3))
       value_pat(i,nt)=fcst(i,2)
      enddo
      print*,'extreme members'
c     print*,(fcst(i,1),i=1,(mem*mem-mem)/2)
c     print*,'distance:',(fcst(i,2),i=(mem*mem-mem)/2,1,-1)
      print*,'distance:',(value_pat(i,nt),i=(mem*mem-mem)/2,1,-1)
      print*,'member IDs:',(int(fcst(i,3)),i=(mem*mem-mem)/2,1,-1)
      do i=(mem*mem-mem)/2,1,-1
      print*,'(',index_pat(int(fcst(i,3)),1,nt),',',
     &index_pat(int(fcst(i,3)),2,nt),')'
      enddo
	print*, 'ok till here 3 -- pattern distance calculation'
      endif

C Write out wgts in grib format
      if (meth.eq.1) then   !write out only the first method (magnitude difference with ensemble mean as ref)
C kpds might be mistakenly altered
C        kpds(4)=0  !bitmap
         kpds(5)=184
         kpds(6)=200
         kpds(7)=0
         kpds(8)=yy-(yy/100)*100
         kpds(9)=mm
         kpds(10)=dd
         kpds(11)=hh
         kpds(12)=0
         kpds(14)=(nt-1)*3
         kpds(15)=0
         kpds(19)=129
         kpds(21)=21
         kpds(22)=4  !precision
c        kpds(23)=2  !ensemble
         kpds(23)=0  !not ensemble

        print*,'make sure all kpds are correctly set'
        print*,'kpds(1-25)=',(kpds(i),i=1,25)

        write(date,'(i2.2)') (nt-1)*3
        tot_mag1=0.
        do i=1,mem
         tot_mag1(nt)=tot_mag1(nt)+value_mag1(meth,mem,nt)-
     &value_mag1(meth,i,nt)
        enddo

        do i=1,mem
         write(order,'(i2.2)') lorder_mag1(meth,i,nt)
         wgts='wgt.pgrb212.' // order // '.f' // date
         call baopen(9+i,wgts,ierr)

         do ipt=1,ngrid
          weight(ipt)=(value_mag1(meth,mem,nt)-
     &value_mag1(meth,i,nt))/tot_mag1(nt)
          if(weight(ipt).le.0.001) weight(ipt)=0.001
         enddo

         call putgb(9+i,lat*lon,kpds,kgds,lb,weight,iret)
c        call putgb(9+i,jf,kpds,kgds,lb,weight,iret)
         call baclose(9+i,ierr)
        enddo !mem

      endif  !meth=1

      enddo  !meth in calculation

10	continue

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Part III: Output 
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
       write(90,501)
       write(90,502)
       write(90,503)
       write(90,504)
       write(90,*)
       write(90,500) hh,mm,dd,yy
500   format(15hInitiating from,1x,i2,1hz,1x,i2,1h/,i2,1h/,i4)
501   format(41hmember order:  1   2   3   4   5   6   7 ,
     &55h  8   9  10  11  12  13  14  15  16  17  18  19  20  21,
     &20h  22  23  24  25  26)
502   format(42hmember names: Bc0,Bn1,Bp1,Bn2,Bp2,Bn3,Bp3 ,
     &51hBn4,Bp4,Bn5,Bp5,Bn6,Bp6,Ac0,An1,Ap1,An2,Ap2,An3,Ap3,
     &23hAn4,Ap4,An5,Ap5,An6,Ap6)
503   format(45hNote: B=NMMB model; N=NMM model; A=ARW model.,
     &53h Wgts frm the 1st approach below is suggested to use.)
504   format(52hMethod: Du-Zhou Dynamical Performance-Ranking Method,
     &41h (Du and Zhou, Oct. 2011, Mon. Wea. Rew.))

C Combined linearity and Spread
       ave_alinear=0.
       ave_avevar=0.
       do nt=1,ntime
        do nv=1,nvar
         if(nv.eq.1) then !SLP
          ave_alinear(nt)=ave_alinear(nt)+alinear(nt,nv)/100.
          ave_avevar(nt)=ave_avevar(nt)+avevar(nt,nv)/100.
         else
          ave_alinear(nt)=ave_alinear(nt)+alinear(nt,nv)
          ave_avevar(nt)=ave_avevar(nt)+avevar(nt,nv)
         endif
        enddo
       enddo

c      write(90,510)
510   format(33hnonlinearity and spread over time)
c       write(90,750) (interv*(nt-1),nt=1,ntime)
cc       write(90,700) (alinear(nt,1)/100,nt=1,ntime)
cc       write(90,700) (avevar(nt,1)/100,nt=1,ntime)
cc       write(90,700) (alinear(nt,7),nt=1,ntime)
cc       write(90,700) (avevar(nt,7),nt=1,ntime)
c       write(90,700) (ave_alinear(nt)/nvar,nt=1,ntime)
c       write(90,700) (ave_avevar(nt)/nvar,nt=1,ntime)
c       write(90,*)
700   format(30f4.1)
750   format(30i4)

C Member performance ranking
C--------by magnitude------------------------------
       do meth=1,2  !meth in output
       if(meth.eq.1) write(90,521) meth
       if(meth.eq.2) write(90,522) meth
       write(90,520)
       write(95,520)
521   format(24hDerived from the method ,i1,
     &24h (ens mean as reference))
522   format(24hDerived from the method ,i1,
     &26h (ens median as reference))
520   format(41h    member from best (1) to worst (26) by,
     &39h "magnitude measure" over forecast time)
        write(90,900) 'rnk ',(interv*(nt-1),nt=1,ntime)
        write(95,900) 'rnk ',(interv*(nt-1),nt=2,ntime,2)
        do i=1,mem
         write(90,900) ' ',i,(lorder_mag1(meth,i,nt),nt=1,ntime)
         write(95,900) ' ',i,(lorder_mag1(meth,i,nt),nt=2,ntime,2)
        enddo
        write(90,*)
900   format(1a,i2,30i3)

       write(90,530)
530   format(44h    "distance" for best (1) to worst (26) by,
     &39h "magnitude measure" over forecast time)
        write(90,802) 'rnk',(interv*(nt-1),nt=1,ntime)
        do i=1,mem
         write(90,901) ' ',i,(value_mag1(meth,i,nt),nt=1,ntime)
        enddo
        write(90,*)
802   format(1a,30i5)
901   format(1a,i2,30f5.2)

       write(90,540)
       write(95,540)
540   format(37h    wgt for best (1) to worst (26) by,
     &39h "magnitude measure" over forecast time)
        write(90,803) 'rnk',(interv*(nt-1),nt=1,ntime)
        write(95,803) 'rnk',(interv*(nt-1),nt=2,ntime,2)
        tot_mag1=0.
        tot=0.
        do nt=1,ntime
        do i=1,mem
         tot_mag1(nt)=tot_mag1(nt)+value_mag1(meth,mem,nt)-
     &value_mag1(meth,i,nt)
        enddo
        enddo
        do i=1,mem
         do nt=1,ntime
      wgt_mag1(i,nt)=(value_mag1(meth,mem,nt)-
     &value_mag1(meth,i,nt))/tot_mag1(nt)
         if(wgt_mag1(i,nt).le.0.001) wgt_mag1(i,nt)=0.001
         tot(nt)=tot(nt)+wgt_mag1(i,nt)
         enddo
         write(90,903) ' ',i,(wgt_mag1(i,nt),nt=1,ntime)
         write(95,903) ' ',i,(wgt_mag1(i,nt),nt=2,ntime,2)
        enddo
c        write(90,903) ' ',i,(tot(nt),nt=1,ntime)
        write(90,*)
803   format(1a,30i5)
903   format(1a,i2,30f5.3)

C-----------by pattern---------------------------
       write(90,550)
       write(95,550)
550   format(41h    member from best (1) to worst (26) by,
     &37h "pattern measure" over forecast time)
        write(90,900) 'rnk ',(interv*(nt-1),nt=1,ntime)
        write(95,900) 'rnk ',(interv*(nt-1),nt=2,ntime,2)
        do i=mem,1,-1
      write(90,900) ' ',mem-i+1,(lorder_pat1(meth,i,nt),nt=1,ntime)
      write(95,900) ' ',i,(lorder_pat1(meth,i,nt),nt=2,ntime,2)
        enddo
        write(90,*)

       write(90,560)
560   format(46h    "corr coeff" for best (1) to worst (26) by,
     &37h "pattern measure" over forecast time)
        write(90,802) 'rnk',(interv*(nt-1),nt=1,ntime)
        do i=mem,1,-1
      write(90,902) ' ',mem-i+1,(value_pat1(meth,i,nt),nt=1,ntime)
        enddo
        write(90,*)
902   format(1a,i2,30f5.1)

       write(90,570)
       write(95,570)
570   format(37h    wgt for best (1) to worst (26) by,
     &37h "pattern measure" over forecast time)
        write(90,803) 'rnk',(interv*(nt-1),nt=1,ntime)
        write(95,803) 'rnk',(interv*(nt-1),nt=2,ntime,2)
        tot_pat1=0.
        tot=0.
        do nt=1,ntime
        do i=1,mem
         tot_pat1(nt)=tot_pat1(nt)+value_pat1(meth,i,nt)-
     &value_pat1(meth,1,nt)
        enddo
        enddo
        do i=mem,1,-1
         do nt=1,ntime
      wgt_pat1(i,nt)=(value_pat1(meth,i,nt)-
     &value_pat1(meth,1,nt))/tot_pat1(nt)
         if(wgt_pat1(i,nt).le.0.001) wgt_pat1(i,nt)=0.001
         tot(nt)=tot(nt)+wgt_pat1(i,nt)
         enddo
      write(90,903) ' ',mem-i+1,(wgt_pat1(i,nt),nt=1,ntime)
      write(95,903) ' ',mem-i+1,(wgt_pat1(i,nt),nt=2,ntime,2)
        enddo
c        write(90,903) ' ',mem-i+1,(tot(nt),nt=1,ntime)
        write(90,*)
        enddo   !meth in output

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C Member clustering:
c       write(90,580)
c580   format(51hfarest member pair by "magnitude measure" over time)
c        write(90,900) (interv*(nt-1),nt=1,ntime)
c        write(90,900) 
c     &(index_mag(lorder_mag((mem*mem-mem)/2,nt),1,nt),nt=1,ntime)
c        write(90,900) 
c     &(index_mag(lorder_mag((mem*mem-mem)/2,nt),2,nt),nt=1,ntime)
c        write(90,*)

c       write(90,590)
c590   format(49hfarest member pair by "pattern measure" over time)
c        write(90,900) (interv*(nt-1),nt=1,ntime)
c        write(90,900) 
c     &(index_pat(lorder_pat(1,nt),1,nt),nt=1,ntime)
c        write(90,900) 
c     &(index_pat(lorder_pat(1,nt),2,nt),nt=1,ntime)
c        write(90,*)

c       write(90,600)
c600   format(52hnearest member pair by "magnitude measure" over time)
c        write(90,900) (interv*(nt-1),nt=1,ntime)
c        write(90,900) 
c     &(index_mag(lorder_mag(1,nt),1,nt),nt=1,ntime)
c        write(90,900) 
c     &(index_mag(lorder_mag(1,nt),2,nt),nt=1,ntime)
c        write(90,*)
c       write(90,610)

c610   format(50hnearest member pair by "pattern measure" over time)
c        write(90,900) (interv*(nt-1),nt=1,ntime)
c        write(90,900) 
c     &(index_pat(lorder_pat((mem*mem-mem)/2,nt),1,nt),nt=1,ntime)
c        write(90,900) 
c     &(index_pat(lorder_pat((mem*mem-mem)/2,nt),2,nt),nt=1,ntime)
c        write(90,*)
c      write(90,620)

620   format(51hpairs from closest to farest by "magnitude measure")
c      write(90,801) (interv*(nt-1),nt=2,ntime,2)
c     do i=1,(mem*mem-mem)/2
c     write(90,910) (index_mag(lorder_mag(i,nt),1,nt),
c    &index_mag(lorder_mag(i,nt),2,nt),nt=2,ntime,2)
c     enddo
c       write(90,*)
801   format(30i7)
910   format(15(1h(,i2,1h,,i2,1h)))

c      write(90,630)
630   format(49hpairs from closest to farest by "pattern measure")
c      write(90,801) (interv*(nt-1),nt=2,ntime,2)
c     do i=(mem*mem-mem)/2,1,-1
c     write(90,910) (index_pat(lorder_pat(i,nt),1,nt),
c    &index_pat(lorder_pat(i,nt),2,nt),nt=2,ntime,2)
c     enddo
c       write(90,*)

c      write(90,640)
640   format(40hneighboring pairs by "magnitude measure")
c     do nt=ntime,ntime
c      do id=1,mem-1
c       do i=1,10 !just take first 10 closest pairs
c      if(                          index_mag(lorder_mag(i,nt),1,nt).
c    &eq.id) write(90,911) 3*(nt-1),index_mag(lorder_mag(i,nt),1,nt),
c    &  index_mag(lorder_mag(i,nt),2,nt)
c       do i=1,(mem*mem-mem)/2  !check all for close pairs with a threshold
c      if(value_mag(i,nt).lt..8.and.index_mag(lorder_mag(i,nt),1,nt).
c    &eq.id) write(90,911) 3*(nt-1),index_mag(lorder_mag(i,nt),1,nt),
c    &  index_mag(lorder_mag(i,nt),2,nt)
c       enddo
c      enddo
c     enddo
c     write(90,*)
911   format(i2,1hh,1h(,i2,1h,,i2,1h))

c      write(90,650)
650   format(38hneighboring pairs by "pattern measure")
c     do nt=ntime,ntime
c      do id=1,mem-1
c       do i=(mem*mem-mem)/2,(mem*mem-mem)/2-10,-1
c       if(                           index_pat(lorder_pat(i,nt),1,nt).
c    &eq.id) write(90,911) interv*(nt-1),
c    &  index_pat(lorder_pat(i,nt),1,nt),
c    &  index_pat(lorder_pat(i,nt),2,nt)
c       do i=(mem*mem-mem)/2,1,-1
c       if(value_pat(i,nt).gt.85..and.index_pat(lorder_pat(i,nt),1,nt). 
c    &eq.id) write(90,911) interv*(nt-1),
c    &  index_pat(lorder_pat(i,nt),1,nt),
c    &  index_pat(lorder_pat(i,nt),2,nt)
c       enddo
c      enddo
c     enddo


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

      return
      end

C......................................................
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
        if(r1.eq.0.0.and.r2.eq.0.0) r=1.0
        if(r1.eq.0.0.and.r2.ne.0.0) r1=0.01
        if(r2.eq.0.0.and.r1.ne.0.0) r2=0.01
        if(r1.ne.0.0.and.r2.ne.0.0)
     &  r=r/sqrt(r1)/sqrt(r2)*100.

        return
9999    end

c**********************************
      subroutine mean_median(n,fcsts,ave,xmed)

c...........................................................................
c  Purpose: computing ensemble mean and median for an n-member ensemble
c           at a fixed point
c
c  Input:
c	 n=the number of ensemble members
c	 fcsts=n forecasts in the ensemble
c  Output:
c	 ave=ensemble mean
c	 xmed=ensemble median
c Programer: Jun Du
c...........................................................................

      dimension fcsts(n),ens(n,3)
c     real fcsts(*),ens(n,3)

c 1. Initialization and ensemble mean
      ave=0.
      xmed=0.

      do 100 i=1,n
        ens(i,1)=i
        ens(i,2)=fcsts(i)
        ave=ave+fcsts(i)/float(n)
100   continue

c 2. Ordering ensemble data in an ascending manner
      call sortm(ens,n,3,2)

c 3. Finding median
      im=n/2
      if(mod(n,2).eq.1) then
       xmed=ens(im+1,2)
      else
       xmed=(ens(im,2)+ens(im+1,2))/2.
      endif

      return
      end

C**************************************
      subroutine sortm(a,n,nc,k)
      dimension a(n,nc),b(n,nc),js(n)
c
      do i1 = 1, n
      iless = 0
      imore = 0
      ieq   = 0
      aa=a(i1,k)
      do i2 = 1, n
       bb=a(i2,k)
       if ( aa.lt.bb ) iless = iless + 1
       if ( aa.gt.bb ) imore = imore + 1
       if ( aa.eq.bb ) then
          ieq   = ieq   + 1
          js(ieq) = i2
       endif
      enddo
       if ( ieq.eq.1) then
          b(imore+1,2)=aa
          b(imore+1,1)=i1
       else
        do i3 = 1, ieq
          b(imore+i3,2)=aa
          b(imore+i3,1)=js(i3)
        enddo
       endif
      enddo
      do jj= 1, n
        a(jj,3) = b(jj,1)
        a(jj,2) = b(jj,2)
      enddo
      return
      end

