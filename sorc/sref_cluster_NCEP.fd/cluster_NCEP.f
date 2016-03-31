       program cluster
c
c This is Steve Tracton's clustering algorithm rewitten by Jun Du to apply
c it to SREF forecasts. This program extracts the clusters based on either 
c z500 or mslp on AWIPS grid 212 (40km).
c
c The code first finds the farthest pair of ensemble members and clusters
c members close to them (AC> a threshold).
c
c cls(nm,6) = 6 clusters possible, nm forecasts assigned 0 or 1
c             (0=not in cluster; 1=included in cluster)
c nx(nm)=0 or 1; 0 if not used in particular cluster, 1 if used
c ny(nm)=0 or 1; 0 if not used in any cluster, 1 if used
c ncl(6)=number of members in each of the 6 possible clusters
c tcl=total number of forecasts used in all clusters (maybe < nm due to the 
c     limited numbers of allowed clusters)
c
c nm=size of ensemble (number of members)
c
c Change Log:
c xxx xx, 1996: Steve Tracton, orginal clustering algrithm
c May 25, 2011: Jun Du, modified for SREF forecast
c Mar 16, 2015: Jun Du, membership increased to 26 (SREF.v7.0.0)
c

      parameter(nm=26,nfile=nm+1,ntime=30)

C nfile - 1-26:ens mems; 27:climo
      parameter(jf=185*129,igx=185,igy=129,ngrid=igx*igy)  !for Grid 212
c     parameter(jf=512*256,igx=185,igy=129,ngrid=igx*igy) 

      dimension f(jf),data(ngrid,nfile),fcst(nfile,igx,igy)
     &,climo(igx,igy),pr(6,igx,igy),f1(igx,igy),f2(igx,igy)
     &,cross(nm,nm),prc(ngrid)
      integer cls(nm,6),ny(nm),nx(nm),ncl(6),nz(nm),kz(6)
      dimension kftime(nfile,ntime),jpds(25),jgds(22),kpds(25),kgds(22)
      integer tcl,tc,iclus,ifile
      integer yy1,mm1,dd1,hh1,iar,meth
      integer yy3(ntime),mm3(ntime),dd3(ntime),hh3(ntime)
      logical lb(jf)
      logical li(ngrid)
      character*6  mem(nm)
      character*2  cyc
      character*1  clus
      character*2  hr
      character*13  fname1
      character*15  fname2
      character*30  fname3
      common /j3/ix1,ix2,iy1,iy2
      common /cl/climo,clstd
      namelist/namin/yy1,mm1,dd1,hh1,iar,meth

      open(50,file='sref_cluster_info',status='unknown')
      open(51,file='sref_clustering_process',status='unknown')

       ix1=1
       ix2=igx  ! specifies region
       iy1=1
       iy2=igy 

c Passing over date information
      read(5,namin,end=1000)
1000  continue
      print*,yy1,mm1,dd1,hh1,iar,meth
      write(50,100) hh1,mm1,dd1,yy1
c     write(50,*)
      if (iar.eq.1.and.meth.eq.1) write(51,*) 
     &'NCEP method: clustering computed using corr coeff with Z500'
      if (iar.eq.1.and.meth.eq.2) write(51,*) 
     &'NCEP method: clustering computed using RMS diff with Z500'
      if (iar.eq.2.and.meth.eq.1) write(51,*) 
     &'NCEP method: clustering computed using corr coeff with MSLP'
      if (iar.eq.2.and.meth.eq.2) write(51,*) 
     &'NCEP method: clustering computed using RMS diff with MSLP'
100   format(22hForecast clusters for ,i2.2,2hz ,i2.2,1h/,i2.2,1h/,i4,
     *62h of 26-member SREF using NCEP method (min: 1; max: 6 clusters))

      read(7,*) (yy3(i),i=1,ntime),(mm3(i),i=1,ntime),
     &(dd3(i),i=1,ntime),(hh3(i),i=1,ntime)
      print*,(yy3(i),i=1,ntime),(mm3(i),i=1,ntime),
     &(dd3(i),i=1,ntime),(hh3(i),i=1,ntime)
      read(7,*) ((kftime(i,j),j=1,ntime),i=1,nfile)
      read(7,*) (mem(i),i=1,nm)
      do i=1,nfile
       print *,(kftime(i,j),j=1,ntime)
      enddo
       print *,(mem(i),i=1,nm)
c     write(50,*) (i,mem(i),i=1,nm)
      write(50,*) (i,mem(i),i=1,1)

      do 10 nt=1,ntime      !all fcst hours
c     do 10 nt=2,ntime,2    !selected fcst hours
      write(50,*) '************************************************'
      write(50,200) kftime(1,nt)
      write(50,*) '************************************************'
      write(51,*)
      write(51,*) '************************************************'
      write(51,200) kftime(1,nt)
      write(51,*) '************************************************'
200   format(1x,i2.2,14h Forecast Hour)

        print *,' '
        print *,'----------------------------------------------- '
        print *,'NEW TIME: nt = ',nt
        print *,'Corresponding forecast hour: nt = ',kftime(1,nt)
        print *,' '

        j=0
        data=0.0
        lugb=9       
        lugi=39      

CCCCCCCCCCCCCCCCCCCCCCCCCCCC
c part I: Geting data ready:
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
        do 5 ifile=1,nfile
c be careful of year and month when in the turn of time
        jpds=-1
        jgds=-1
C specify date information here:
       if(ifile.le.nm) then
c       jpds(8)=yy1
c       jpds(9)=mm1
        jpds(10)=dd1
        jpds(11)=hh1
       endif
       if(ifile.eq.nm+1) then
c       jpds(8)=yy3(nt)
c       jpds(9)=mm3(nt)
        jpds(10)=dd3(nt)
        jpds(11)=hh3(nt)
       endif
        jpds(14)=kftime(ifile,nt)
        jpds(15)=-1

        print *,' '
        print *,'----------'
        print *,'New File: nt = ',nt,' ifile= ',ifile
        print *,'Fcst Hour = ',kftime(ifile,nt)
        print *,' '

c define variables 
        if(iar.eq.1) then !500z
          jpds(7)=500
          jpds(6)=100
          jpds(5)=7
        else
          jpds(7)=0          
          jpds(6)=102        
          jpds(5)=2        !mslp
        endif

      print*, jpds(5)
      print*, jpds(6)
      print*, jpds(7)
      print*, jpds(8)
      print*, jpds(9)
      print*, jpds(10)
      print*, jpds(11)
      print*, jpds(12)
      print*, jpds(13)
      print*, jpds(14)
      print*, jpds(15)

      lugb=lugb+1
      lugi=lugi+1
      print *,'lugb= ',lugb,' lugi= ',lugi
      call getname(ifile,fname1,fname2)
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

      if (iret.eq.0) then
       do ipt=1,ngrid
        data(ipt,ifile) = f(ipt)
       enddo
      else
       print *, "Something wrong here!!!!"
      endif

c       print*, 'original data:',ifile, data(1,ifile)
5       continue    !files (mems, climo)
        print*, 'ok till here 1'

c Convert forecast from data(ngrid,nm) to fcst(nm,igx,igy)
       do ifile=1,nfile
        i=0
        j=1
        do 40 ipt=1,ngrid
         i=i+1
         fcst(ifile,i,j)=data(ipt,ifile)
c        print*,'i,j,ipt=',i,j,ipt
         if(i.eq.igx) then
          i=0
          j=j+1
         endif
40      continue 
c      print*, '2nd original data:',ifile, fcst(ifile,1,1)
       enddo

       do i=1,igx
        do j=1,igy
         climo(i,j)=fcst(nfile,i,j)
        enddo
       enddo
       print*,'3rd original climo=',climo(1,1)
c Now calculate the domain-averaged standard deviation (for meth=2 only)
       if(meth.eq.2) then
       clstd=0.
       do 847 j=1,igy
        sumc=0.
        do 708 i=1,igx
         sumc=sumc+climo(i,j)/igx
708     continue
        clstd=clstd+sumc/igy
847    continue
       print*,'The domain averaged sd is ',clstd
       endif

CCCCCCCCCCCCCCCCCCCCCCCCCCCC
c part II: Finding Clusters:
CCCCCCCCCCCCCCCCCCCCCCCCCCCC
       if(meth.eq.1) then
        thresh=0.93  ! ac threshold
c       thresh=0.99  ! for testing to have more clusters
       else
        thresh=0.75  ! rmse threshold
       endif
       print*,'thresh=',thresh

c initialization at every forecast hour
        do 11 l=1,6
        kz(l)=0
        ncl(l)=0
        do 12 k=1,nm
        cls(k,l)=0
        do 13 j=1,nm
        cross(j,k)=9999.0
13      continue
12      continue
11      continue
        tcl=0
        tc=0
        nf=0
        do 14 j=1,nm
        nx(j)=0
        ny(j)=0
        nz(j)=0
14      continue

c
c Get correlations amongst all possible pairs of ensemble members
c
       do 2000 m=1,nm
       do 2000 k=1,nm
        do 2012 i=1,igx
        do 2012 j=1,igy
        f1(i,j)=fcst(m,i,j)
        f2(i,j)=fcst(k,i,j)
2012    continue
        ans=9999.
        if(m.ge.k) goto 2050
c       print*,'m,k,Calling anglim ',m,k
        call anglim(ans,f1,f2,meth)
        cross(m,k)=ans
        goto 2051
2050    if(m.eq.k) cross(m,k)=1.0
        if(m.gt.k) cross(m,k)=cross(k,m)
2051    write(51,2001) m,k,cross(m,k),ans
2001    format(/2x,'  m,k,cross(m,k),ans,= ',2i4,2f10.3)
2000   continue
c
c Find pair with largest difference (smallest acc between pair members)
c
       amax=1.0
       do 2005 m=1,nm-1
       do 2005 k=m+1,nm
        if (cross(m,k).gt.9990.) goto 2005
        if (cross(m,k).lt.amax) goto 2006
        goto 2005
2006    mmax=m
        kmax=k
        amax=cross(m,k)
2005   continue
       write(51,2007) mmax,kmax,amax
2007   format(/2x,'Largest diff for forecasts m,k,amax ',2i4,f10.3)
       n1=mmax
       n2=kmax
c
c Find pair most similar (largest acc)
c
       amin=-1.0
       do 2008 m=1,nm-1
       do 2008 k=m+1,nm
        if(cross(m,k).gt.9990.) goto 2008
        if(ny(m).eq.1.or.ny(k).eq.1) goto 2008
        if(cross(m,k).gt.amin) goto 2009
        goto 2008
2009    mmin=m
        kmin=k
        amin=cross(m,k)
2008   continue
       write(51,2010) mmin,kmin,amin
2010   format(/2x,'Smallest diff for forecasts m,k,amin ',2i4,f10.3)
c
c Cluster ensemble members:  Start with 2 forecasts farthest apart.  First
c 2 clusters consist of forecasts closest to them. No clustering with 
c extremes of they are close (e.g., amax>thresh). Search for closest of 
c remaining forecasts and cluster. A given forecast can be in only one cluster.
c
       tcl=0
       do 8000 kl=1,6
       write(51,9876) kl
9876   format(//10x,' Cluster',i3,' coming up'//)
       do 6010 k=1,nm
       nx(k)=0
6010   continue
       if(kl.gt.2) goto 4009
       if(amax.le.thresh) goto 7011
c
c Note: Extremes allowed (if close) to be included in later clustering
c
       write(51,7012)
7012   format(//'*** Extremes are close; no clustering with them ***')
       write(50,300) kl,ncl(kl)
       goto 8000
7011   continue
       ncl(kl)=1
       if(kl.eq.2) goto 7015
       nx(n1)=1
       tcl=tcl+1
       ny(n1)=1
       cls(n1,kl)=1
       m=n1
       goto 7013
7015   nx(n2)=1
       tcl=tcl+1
       ny(n2)=1
       cls(n2,kl)=1
       m=n2
c
c Find all members close to each extreme; if close to both don't include
c them in either.
c
7013    do 4008 k=1,nm
         if(cross(m,k).gt.9990.) goto 4008
         if(k.eq.m) goto 4008
         if(ny(k).eq.1) goto 4008
         write(51,8010) kl,m,k,cross(m,k)
8010     format(5x,'kl,m,k,cross(m,k)= ',3i3,f10.3)
         if (cross(m,k).le.thresh) goto 4008
         if(kl.eq.1) nn=n2
         if(kl.eq.2) nn=n1
         write(51,8008) kl,nn,k,cross(nn,k)
         if(cross(nn,k).gt.thresh) goto 4008
8008     format(10x,'kl,nn,k,cross(nn,k) =',3i3,f10.3)
         cls(k,kl)=1
         nx(k)=1
         ny(k)=1
         ncl(kl)=ncl(kl)+1
         tcl=tcl+1
         write(51,8002) kl,tcl,m,k,cross(m,k),cls(k,kl)
8002     format(15x,'kl,tcl,m,k,cross(m,k),cls(k,kl)= ',4i3,f10.3,i3)
4008    continue
        goto 7001
4009    continue
c
c Cluster the remaining ensemble members. Start with pair closest together
c if amin>thresh. If amin<thresh, no clustering. Search with members
c with acc>thresh
c with respect to both members of this first pair.  Take only the one 
c closest in mean to the first pair. Scan others of this set for those
c which are close to these 3. Add only that which is closest on average.
c Repeat process with any remaining members for the next cluster.
c
c Find pair most similar of remaining forecasts (must be at least 2)
c
        if(tcl.ge.(nm-1)) goto 8005
        amin=-1.0
        do 3008 m=1,nm-1
        do 3008 k=m+1,nm
        if(cross(m,k).gt.9990.) goto 3008
        if(ny(m).eq.1.or.ny(k).eq.1) goto 3008
        if(cross(m,k).gt.amin) goto 3009
        goto 3008
3009    mmin=m
        kmin=k
        amin=cross(m,k)
3008    continue
        write(51,3010) kl,mmin,kmin,amin
3010    format(/2x,'Closest of remaining forecasts; kl,m,k,amin= ',
     *    3i4,f10.3)
        if(amin.le.thresh) goto 4000
        n1=mmin
        n2=kmin
        do 7010 k=1,nm
        nx(k)=0
        nz(k)=0
7010    continue
        nx(n1)=1
        nx(n2)=1
        ny(n1)=1
        ny(n2)=1
        ncl(kl)=2
        tcl=tcl+2
        cls(n1,kl)=1
        cls(n2,kl)=1
        tc=0
c
c Find forecasts close to both: tc=the number that are
c
        do 6000 m=1,nm
        if(cross(m,n1).gt.9990.) goto 6000
        if(cross(m,n2).gt.9990.) goto 6000
        if(m.eq.n1.or.m.eq.n2) goto 6000
        if(ny(m).eq.1) goto 6000
c       write(51,8003) kl,m,cls(m,kl)
8003    format(/2x,'kl,m,cls(m,kl)= ',3i3/)
        write(51,8007) kl,tcl,m,n1,cross(m,n1)
8007    format(10x,'kl,tcl,m,n1,cross(m,n1)= ',4i3,f10.3)
        write(51,9007) kl,tcl,m,n2,cross(m,n2)
9007    format(10x,'kl,tcl,m,n2,cross(m,n2)= ',4i3,f10.3)
        if(cross(m,n1).le.thresh) goto 6000
        if(cross(m,n2).le.thresh) goto 6000
        nz(m)=1
        tc=tc+1
6000    continue
6005    continue
        if(tc.eq.0) goto 7001
        if(tc.gt.1) goto 6009
        do 6008 m=1,nm
c       print*,'Here is ',m
        if(nz(m).eq.0) goto 6008
        nx(m)=1
        ny(m)=1
        cls(m,kl)=1
        ncl(kl)=ncl(kl)+1
        tcl=tcl+1
        goto 7001
6008    continue
6009    continue
c
c If more than one member close to forecasts already in cluster, get that
c one that's closes on average
c
        write(51,7016) kl,tc
7016    format(//8x,'Starting search for closest on avg for kl ',
     *    i3,' with ',i3,' members close')
        amx=0.0
        do 6001 m=1,nm
        if(m.eq.n1) goto 6001
        if(m.eq.n2) goto 6001
        if(nx(m).eq.1) goto 6001
        if(nz(m).eq.0) goto 6001
        t=0.
        do 6002 k=1,nm
        if(cross(m,k).gt.9990.) goto 6002
        if(k.eq.m) goto 6002
        if(nx(k).eq.0) goto 6002
        t=t+cross(m,k)/ncl(kl)
        write(51,7008) kl,m,k,ncl(kl),t,amx,cross(m,k)
7008    format(//5x,'In search for closest on avg '/8x,
     *    'kl,m,k,ncl(kl),t,amx,cross(m,k)= ',4i3,3f10.3/)
6002    continue
        if(t.le.amx) goto 6001
        amx=t
        n3=m
6001    continue
        tc=tc-1
        nx(n3)=1
        ny(n3)=1
        cls(n3,kl)=1
        ncl(kl)=ncl(kl)+1
        tcl=tcl+1
c
c Eliminated from consideration any remaining possiblilties not close to
c newest member of cluster.
c
         do 7005 m=1,nm
         if(cross(m,n3).gt.9990.) goto 7005
         if(m.eq.n3) goto 7005
         if(nz(m).eq.0) goto 7005
         if(cross(m,n3).gt.thresh) goto 7005
         tc=tc-1
         nz(m)=0
7005     continue
         if(tc.eq.0) goto 7001
         if(tc.gt.1) goto 6009
c
c One more left in first go around; add to cluster
c
         do 7004 m=1,nm
         if(nz(m).eq.0) goto 7004
         if(nx(m).eq.1) goto 7004
         nx(m)=1
         ny(m)=1
         cls(m,kl)=1
         ncl(kl)=ncl(kl)+1
         tcl=tcl+1
         write(51,7017) m,kl
7017     format(//5x,'In one more left, add forecast ',i3,'to cluster',
     *    i3)
         goto 7001
7004     continue
7001     do 5002 m=1,nm
         write(51,5003) kl,m,cls(m,kl)
5003     format(//2x,'Cluster ',i2,'; fcst # ',i2,'; cls(m,kl) = ',i2)
5002     continue
         write(51,5004) kl,ncl(kl),tcl,(nx(m),m=1,nm)
5004     format(//2x,'Cluster #',i2,', has ',i2,' forecasts of total',
     *    ' forecasts in clusters =',i3,/
     *    2x,'Forecasts in cluster are ',21i3/)
         write(50,300) kl,ncl(kl)
300      format(1x,9hCluster #,i2,5h has ,i2
     *,44h forecasts containing the following members:)
         do m=1,nm
          if(nx(m).eq.1) write(50,*) ' ',mem(m),' (',m,')'
         enddo
c
c Calculate cluster mean
c
         do 405 i=1,igx
         do 405 j=1,igy
         pr(kl,i,j)=0.0
405      continue
         nf=0
         do 6100 k=1,nm
         if(nx(k).eq.0) goto 6100
         nf=nf+1
         do 6101 i=1,igx
         do 6101 j=1,igy
         pr(kl,i,j)=pr(kl,i,j)+fcst(k,i,j)
6101     continue
c        print*,'nf,kl,pr(kl,20,20)=',nf,kl,pr(kl,20,20)
6100     continue
c        print*,'nf,kl,pr(kl,20,20)=',nf,kl,pr(kl,20,20)
         do 6102 i=1,igx
         do 6102 j=1,igy
          pr(kl,i,j)=pr(kl,i,j)/float(nf)
6102     continue
         goto 8000
4000     write(51,5000) kl,thresh
5000     format(12x,'For kl=',i2,' closest forecasts have ACC<',
     *     f4.2,' no clusters') 
8005     write(51,8004) kl,tcl
8004     format(//10x,'kl,tcl=',2i3,' only 1 or 0 forecasts left, so
     *     no more clusters possible')
         write(50,300) kl,0
8000     continue
c
c Write out cluster mean in grib format (array pr). Each cluster in its own separate file
c         
         do 888 k=1,6
          if(ncl(k).eq.0) goto 888
          iclus=50+k
c         print*,'iclus=',iclus
c
c convert 2-d array to 1-d array
c
           i=0
           j=1
           do 97 igrid=1,ngrid
            prc(igrid)=0.0
97         continue
           do 98 igrid=1,ngrid
            i=i+1
c           if(k.eq.3) then
            prc(igrid)=pr(k,i,j)
c           print*,'k,i,j,igrid,prc(igrid),pr(k,i,j)='
c    *      ,k,i,j,igrid,prc(igrid),pr(k,i,j)
c           endif
            if(i.eq.igx) then
             i=0
             j=j+1
            endif
98        continue

c Output in grib format
      write(clus,'(i1)') k
      write(hr,'(i2.2)') kftime(1,nt)
      write(cyc,'(i2.2)') hh1
      fname3='ref.t'//cyc//'z.pgrb212.clus'//clus//'mean.f'//hr
      call baopen(iclus,fname3,ierr)
      kpds(1)=7
      kpds(2)=130
      kpds(3)=212
      kpds(4)=0
c     kpds(5)=variable dependent
c     kpds(6)=variable dependent
c     kpds(7)=variable dependent
      kpds(8)=yy1-(yy1/100)*100
      kpds(9)=mm1
      kpds(10)=dd1
      kpds(11)=hh1
      kpds(12)=0
      kpds(13)=1
c     kpds(14)=variable dependent
c     kpds(15)=variable dependent
c     kpds(16)=variable dependent
      kpds(17)=-1
      kpds(18)=1
      kpds(19)=2
      kpds(20)=-1
      kpds(21)=21
c     kpds(22)=variable dependent (units decimal scale factor/precision)
c     kpds(23)=2  !if ensemble data
      kpds(23)=0  !if not ensemble data
      kpds(24)=128
      kpds(25)=-1

      kgds(20)=255

C-----------------z500 or mslp---------------------
      if(iar.eq.1) then !500z
       jpds(5)=7
       jpds(6)=100
       jpds(7)=500
      else
       jpds(5)=2        !mslp
       jpds(6)=102        
       jpds(7)=0          
      endif
      kpds(14)=kftime(1,nt)
      kpds(15)=0
      kpds(16)=0
      kpds(22)=1
      call putgb(iclus,jf,kpds,kgds,lb,prc,iret)
888      continue
      call baclose(iclus,ierr)
c     write(50,*)
c
c Now finished with forecast time period, move onto the next forecast time period.
10       continue
c        call baclose(iclus,ierr)
         stop
         end


         subroutine anglim(ans,f1,f2,meth)
c
c This is Steve Tracton's ANGLIM subroutine
c Gridpoint computations of the anomaly correlation coefficient and 
c the RMS error.
c
c Set for AWIPS grid 212 (1=211, 2=212)
c meth=1 is ACC, meth=2 is RMSE
c
         parameter(nm=26,igx=185,igy=129)

         dimension climo(igx,igy),f1(igx,igy),f2(igx,igy)
         common /j3/ix1,ix2,iy1,iy2
         common /cl/climo,clstd
    
         ipts=ix2-ix1+1
         jpts=iy2-iy1+1
         if(meth.eq.2) goto 205
         sumcov=0.
         sumvrx=0.
         sumvry=0.
         sumxb=0.
         sumyb=0.
         smxbyb=0.
         sumxbs=0.
         sumybs=0.
         cntla=0.
         sumx=0.
         sumy=0.
c        do ll=40,50
c        print*,'ll= ',ll
c        print*,'f1(ll,ll),climo(ll,ll)=',f1(ll,ll),climo(ll,ll)
c        print*,'f2(ll,ll),climo(ll,ll)=',f2(ll,ll),climo(ll,ll)
c        enddo
         do 37 la=iy1,iy2
          cntla=0.
          sumx=0.
          sumy=0.
          sumxy=0.
          sumxx=0.
          sumyy=0.
          do 73 lon=ix1,ix2
           lo=lon
           xa=f1(lo,la)-climo(lo,la)
           ya=f2(lo,la)-climo(lo,la)
c          if(xa.gt.1000.0) then
c            print*,'xa=',xa
c            goto 73
c          endif
c          print*,'xa,ya=',xa,ya
c          print*,'f1(lo,la),f2(lo,la),climo(lo,la)=',
c    *      f1(lo,la),f2(lo,la),climo(lo,la)
           sumx=sumx+xa/ipts
           sumy=sumy+ya/ipts
           sumxy=sumxy+xa*ya/ipts
           sumxx=sumxx+xa*xa/ipts
           sumyy=sumyy+ya*ya/ipts
73        continue
c         print*,'sumx,sumy,sumxy,sumxx,sumyy=',
c    *      sumx,sumy,sumxy,sumxx,sumyy
          xbar=sumx
          ybar=sumy
          xybar=sumxy
          xxbar=sumxx
          yybar=sumyy
c
c Get averages over current la belt
c
          covla=xybar-xbar*ybar
          varxla=xxbar-xbar*xbar
          varyla=yybar-ybar*ybar
c
c Increment sums over all la
c
          sumcov=sumcov+covla
          sumvrx=sumvrx+varxla
          sumvry=sumvry+varyla
          sumxb=sumxb+xbar
          sumyb=sumyb+ybar
          smxbyb=smxbyb+xbar*ybar
          sumxbs=sumxbs+xbar*xbar
          sumybs=sumybs+ybar*ybar
37       continue
c        stop
c        print*,'sumcov,smxbyb,sumxb,jpts=',sumcov,smxbyb,sumxb,jpts
c        print*,'sumvrx,sumxbs=',sumvrx,sumxbs
c        print*,'sumvry,sumybs,sumyb=',sumvry,sumybs,sumyb 
         cov=(sumcov+smxbyb-(sumxb*sumyb)/jpts)/jpts
         varx=(sumvrx+sumxbs-(sumxb*sumxb)/jpts)/jpts
         vary=(sumvry+sumybs-(sumyb*sumyb)/jpts)/jpts 
c        print*,'cov,varx,vary,varx*vary,sqrt(varx*vary)=',
c    *    cov,varx,vary,varx*vary,sqrt(varx*vary)
         if((varx.gt.0).or.(vary.gt.0)) goto 12
         write(6,101)
         ans=9999.
         return
12       ans=cov/sqrt(varx*vary)
         print*,'ans=',ans
101      format('%%%% Trouble in Ancoree %%%%')
         return
205      continue
c
c Calculate RMS difference between forecasts
c
         sumsq=0.
         do 87 lat=iy1,iy2
          sumlat=0.
          do 93 lon=ix1,ix2
            lo=lon
            diff=f2(lo,lat)-f1(lo,lat)
            sumlat=sumlat+diff**2/ipts
93         continue
           sumsq=sumsq+sumlat/jpts
87       continue
         ans=sqrt(sumsq)/clstd
         ans=1.0-ans
         return
         end

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

      return
      end

