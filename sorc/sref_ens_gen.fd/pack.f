	subroutine pack(itime,eps,vrbl_mn,vrbl_sp,vrbl_pr,
     +      derv_mn,derv_sp,derv_pr, ptype_mn,ptype_pr,
     +      derv_dtra,prcpmax,prcpmin,snowmax,snowmin,
     +      frznmax,frznmin,jf,iens,iyr,imon,idy,ihr,
     +      gribid,lb,kgds,interval) 

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
C KPDS(9)  - MONTH
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
C REMARKS:
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM
C
C$$$
C

       include 'parm.inc'

       dimension kpds(25)
       dimension kens(5),kprob(2),xprob(2),kclust(16),kmembr(80)
       character gds(42)
       character*20 mnout
       character*20 spout
       character*20 prout
       character*30 apcpmaxout
       character*30 apcpminout
       character*10 date
       character*4  eps     !ensemble prediction system, either sref or gens 

C for variable table:
        Integer numvar, nderiv
        Character*4 vname(maxvar)
        Integer k5(maxvar), k6(maxvar)
        Character*1 Msignal(maxvar), Psignal(maxvar)
        Integer Mlvl(maxvar), MeanLevel(maxvar,maxmlvl)
        Integer Plvl(maxvar), ProbLevel(maxvar,maxplvl)
        Integer Tlvl(maxvar)
        Character*1 op(maxvar)
        Real    Thrs(maxvar,maxtlvl)
                                                                                                                                                                
c    for derived variables
        Character*4 dvname(maxvar)
        Integer dk5(maxvar), dk6(maxvar)
        Character*1 dMsignal(maxvar), dPsignal(maxvar)
        Integer dMlvl(maxvar), dMeanLevel(maxvar,maxmlvl)
        Integer dPlvl(maxvar), dProbLevel(maxvar,maxplvl)
        Character*1 dop(maxvar)
        Integer dTlvl(maxvar)
        Real    dThrs(maxvar,maxtlvl)
        Integer MPairLevel(maxvar,maxmlvl,2)
        Integer PPairLevel(maxvar,maxplvl,2)
        
                                                                                                                                                        
        common /tbl/numvar,
     +              vname,k5,k6,Mlvl,Plvl,Tlvl,
     +              MeanLevel,ProbLevel,Thrs,
     +              Msignal,Psignal,op
                                                                                                                                                                
        common /dtbl/nderiv,
     +              dvname,dk5,dk6,dMlvl,dPlvl,dTlvl,
     +              dMeanLevel,dProbLevel,dThrs,
     +              dMsignal,dPsignal,MPairLevel,PPairLevel,dop
      
 
        INTEGER,intent(IN) :: itime,jf, iens, iyr,imon,idy,ihr, 
     +                        gribid, kgds(25)
        REAL,dimension(jf,nderiv,maxmlvl),intent(IN) :: derv_mn
        REAL,dimension(jf,nderiv,maxmlvl),intent(IN) :: derv_sp
        REAL,dimension(jf,nderiv,maxplvl,maxtlvl),intent(IN) ::
     +                                                  derv_pr
        REAL,dimension(jf,numvar,maxmlvl),intent(INOUT) :: vrbl_mn
        REAL,dimension(jf,numvar,maxmlvl),intent(IN) :: vrbl_sp
        REAL,dimension(jf,numvar,maxplvl,maxtlvl),intent(IN) ::
     +                                                  vrbl_pr
        REAL,dimension(jf,maxmlvl),intent(IN)        :: prcpmax
        REAL,dimension(jf,maxmlvl),intent(IN)        :: prcpmin
        REAL,dimension(jf,maxmlvl),intent(IN)        :: snowmax
        REAL,dimension(jf,maxmlvl),intent(IN)        :: snowmin
        REAL,dimension(jf,maxmlvl),intent(IN)        :: frznmax
        REAL,dimension(jf,maxmlvl),intent(IN)        :: frznmin
        logical, dimension(jf),intent(IN)            :: lb 

        REAL,dimension(jf,maxmlvl,4),intent(IN)      :: ptype_mn
        REAL,dimension(jf,maxmlvl,4),intent(IN)      :: ptype_pr

        REAL,dimension(jf,maxmlvl,10),intent(IN)      :: derv_dtra

        REAL var(jf)      


            
      if(itime.lt.100) then                                                                                                                           
       write(date,'(i2.2)') itime
      else
       write(date,'(i3.3)') itime
      end if

      write(*,*) 'packing ....'
      write(*,*) itime, jf,iens,numvar,maxplvl,maxtlvl    

      mnout='mean.'//eps//'.f' // trim(date)
      spout='spread.'//eps//'.f'// trim(date)
      prout='prob.'//eps//'.f' // trim(date)

      apcpmaxout='apcpmax.'//eps//'.f' // trim(date)      
      apcpminout='apcpmin.'//eps//'.f' // trim(date)      

      write(*,*) mnout, spout, prout, 
     +   apcpmaxout,apcpminout 

      write(*,*) 'Before baopen'
      call baopen(30,mnout,ierr)
      call baopen(40,spout,ierr)
      call baopen(50,prout,ierr)

      call baopen(60,apcpmaxout,ierr)
      call baopen(70,apcpminout,ierr)
      write(*,*) 'After baopen'
 
      kpds(1)=7
      kpds(2)=113
      kpds(3)=gribid
      kpds(4)=0
c     kpds(5)=variable dependent
c     kpds(6)=variable dependent
c     kpds(7)=variable dependent
      kpds(8)=iyr-(iyr/100)*100
      kpds(9)=imon
      kpds(10)=idy
      kpds(11)=ihr
      kpds(12)=0
      kpds(13)=1
c     kpds(14)=variable dependent
c     kpds(15)=variable dependent
c     kpds(16)=variable dependent
      kpds(17)=-1
      kpds(18)=1
      kpds(19)=2
      kpds(20)=255
      kpds(21)=21
c     kpds(22)=variable dependent (units decimal scale factor/precision)
      kpds(23)=2
      kpds(24)=128
      kpds(25)=-1

      kens(1)=1           !: OCT 41
      kens(2)=5           !: OCT 42
      kens(3)=0           !: OCT 43
c     kens(4)=variable dependent       !: OCT 44
      kens(5)=255         !: OCT 45    !: Original resolution retained
      kclust(1)=iens      !: OCT 61

c       kprob(1)=variable dependent         !: OCT 46
c       kprob(2)=variable dependent         !: OCT 47
c       xprob(1)=variable dependent        !: OCT 48-51
c       xprob(2)=variable dependent        !: OCT 52-55

cccccc 1: write mean and spread for direct variables into mean grib file
                                                                                                                            
        kprob(1) = -1         !: OCT 46
        kprob(2) = -1         !: OCT 47
        xprob(1) = -1.        !: OCT 48-51
        xprob(2) = -1.        !: OCT 52-55
        kpds(14) = itime
        kpds(15) = 0
        kpds(16) = 0
      write(*,*) 'All parameter are set 1'

       do nv = 1, numvar
        if (k5(nv).eq.7.and.k6(nv).eq.215) k1=nv
        if (k5(nv).eq.7.and.k6(nv).eq.3) k2=nv
       end do

       do kk=1,jf
        if(vrbl_mn(kk,k1,1).gt.vrbl_mn(kk,k2,1)) then           !Dec. 30, 2008: force "ceiling > cloudtop" case to be
         vrbl_mn(kk,k1,1)=vrbl_mn(kk,k2,1)-1000.0               ! ceiling=cloudtop-1000.0
         if(vrbl_mn(kk,k1,1).le.0.0) vrbl_mn(kk,k1,1)=0.
        end if
       end do

      write(*,*) 'All parameter are set 2'
       do 201 nv = 1, numvar
        if(trim(Msignal(nv)).eq.'M') then
           if(k5(nv).eq.24.or.k5(nv).eq.11.or.k5(nv).eq.33.or.
     +      k5(nv).eq.34.or.k5(nv).eq.7 .or.k5(nv).eq.32) then
            kpds(22) = 1
           else if (k5(nv).eq.41.or.k5(nv).eq.153) then
            kpds(22) = 6
           else if (k5(nv).eq.51) then     !JunDu, 07/19/2008
            kpds(22) = 4
           else
            kpds(22) = 2
           end if
                                                                                                                            
           kpds(5) = k5(nv)
           kpds(6) = k6(nv)

      write(*,*) 'All parameter are set 3'
           do lv = 1, Mlvl(nv)
 
             kpds(7)  = MeanLevel(nv,lv)
             
             if(k5(nv).eq.131.and.k6(nv).eq.101) then  !special case: kpds(5)(6)(7) in model
                kpds(5) = 24                           !       grib files are different from 
                kpds(6) = 116                          !       those of ensemble posted grib
                kpds(7) = 7680                         !       files set by Jun Du
             end if
                              
           !mean:
             kens(4) = 2
             var = vrbl_mn(:,nv,lv)

             if(k5(nv).eq.61.and.k6(nv).eq.1) then    !for direct/non-derived 1 hour precip
c              kpds(14)=itime-1
              kpds(14)=itime-interval                      !Jun Du: accumulation interval for cluster
              kpds(15)=itime
              kpds(16)= 4 
              kpds(22)= 4 
             end if

      write(*,*) 'Before putgbex 30'
             call putgbex(30,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
c             call print_kpds(nv,1,'_mn',kpds,kens,kprob,xprob)
           !spread:
             kens(4) = 11
             var = vrbl_sp(:,nv,lv)


      write(*,*) 'Before putgbex 40'
             call putgbex(40,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
c             call print_kpds(nv,1,'_sp',kpds,kens,kprob,xprob)
           end do
         end if
201      continue      
 
cccccc 2: write prob  for direct variables into prob grib file

         kpds(5)  = 191
         kpds(14) = itime
         kpds(15) = 0
         kpds(16) = 0
                                                                                                                            
         kens(4)  = 0
                                                                                                                            
         do 202 nv = 1, numvar

          if (k5(nv).eq.212.and.k6(nv).eq.200) then
            kpds(19)=129
          else
            kpds(19)=2
          end if

           if(k5(nv).eq.24.or.k5(nv).eq.11.or.k5(nv).eq.33.or.
     +      k5(nv).eq.34.or.k5(nv).eq.7 .or.k5(nv).eq.32) then
            kpds(22) = 1
           else if (k5(nv).eq.41) then
            kpds(22) = 6
           else
            kpds(22) = 2
           end if

           if(trim(Psignal(nv)).eq.'P') then


             kpds(6) = k6(nv)         
             kprob(1) = k5(nv)
                                                                                                                            
             if(trim(op(nv)).eq.'>') then
               kprob(2) = 2
             else if(trim(op(nv)).eq.'<') then
               kprob(2) = 1
             else if(trim(op(nv)).eq.'-') then
               kprob(2) =3
c             else if(trim(op(nv)).eq.'=') then
c               kprob(2) =4
              else if(trim(op(nv)).eq.'=') then
               kprob(2) = 3
             else
               kprob(2) = -1
             end if
                                                                                                                            
             do lv = 1, Plvl(nv)
                                                                                                                            
               kpds(7) = ProbLevel(nv,lv)
                
             if(k5(nv).eq.131.and.k6(nv).eq.101) then  !special case: kpds(5)(6)(7) in model
                kpds(6) = 116                          !       grib files are different from
                kprob(1) = 24                          !       those of ensemble posted grib
                kpds(7) = 7680                         !       files set by Jun Du
             end if

                                                                                                            
               do 101 lt = 1, Tlvl(nv)
                                                                                                                            
                 if(trim(op(nv)).eq.'>') then
                   xprob(1) = 0.0
                   xprob(2) = Thrs(nv,lt)
                 else if(trim(op(nv)).eq.'<') then
                   xprob(2) = 0.0
                   xprob(1) = Thrs(nv,lt)
                 else if(trim(op(nv)).eq.'-') then
                   if ( lt.LT.Tlvl(nv) ) then
                     xprob(1) = Thrs(nv,lt)
                     xprob(2) = Thrs(nv,lt+1)
                   end if 
cBinbin Ask:     else if(trim(op(nv)).eq.'=') then
c                     xprob(1) = Thrs(nv,lt)
c                     xprob(2) = 0.0

                 else if(trim(op(nv)).eq.'=') then
                     xprob(1) = Thrs(nv,lt)
                     xprob(2) = Thrs(nv,lt)

                 else
                   xprob(1) = -1.0
                   xprob(2) = -1.0
                 end if
                  
                 if(trim(op(nv)).eq.'-'.and.lt.eq.Tlvl(nv)) 
     +                goto 101          
         
                 var = vrbl_pr(:, nv, lv, lt)
                  do i=1,jf
                   if(var(i).lt.0.0) var(i)=0.0
                  end do 

                if(k5(nv).eq.61.and.k6(nv).eq.1) then  !for direct/non-derived 1 hour precip
                  kpds(14)=itime-1
                  kpds(15)=itime
                  kpds(22)= 2
                end if


      write(*,*) 'Before putgbex 50'
                 call putgbex(50,jf,kpds,kgds,kens,kprob,xprob,
     &           kclust,kmembr,lb,var,iret)

c                call print_kpds(nv,1,'_pr',kpds,kens,kprob,xprob)

101           continue
           end do

         end if
202     continue    

cccccc 3: write mean/spread for derived variables into grib files ccccccccccc

        kprob(1) = -1         !: OCT 46
        kprob(2) = -1         !: OCT 47
        xprob(1) = -1.        !: OCT 48-51
        xprob(2) = -1.        !: OCT 52-55

       do 203 nv = 1, nderiv
       if(dk5(nv).ne.501) then                    !DTRA case is put elsewhere

        !non-accumulated derived var mean/spread
        if(dk5(nv).ne.61.and.dk5(nv).ne.79.and.
     &     dk5(nv).ne.78) then

          kpds(14) = itime
          kpds(15) = 0
          kpds(16) = 0

           if(dk5(nv).eq.32.or.dk5(nv).eq.7) then
            kpds(22) = 1
           else
            kpds(22) = 2
           end if
           
           if(dk5(nv).eq.41) then
            kpds(22) = 5
           end if
                                                                                                                                    
           kpds(5) = dk5(nv)
           kpds(6) = dk6(nv)
                                                                                                                                               
           do lv = 1, dMlvl(nv)
                                                                                                                                               
             kpds(7)  = dMeanLevel(nv,lv)
                                                                                                                                               
           !mean:

             if(dk5(nv).eq.140) then   ! dominant precip type

               kens(4) = -1

               !pack rain mean
               kpds(5)=140
               var = ptype_mn(:,lv,4)


               call putgbex(30,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)

               !pack freezing rain mean
               kpds(5)=141
               var = ptype_mn(:,lv,1)


               call putgbex(30,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
 
               !pack sleet mean
               kpds(5)=142



               var = ptype_mn(:,lv,3)
               call putgbex(30,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)

               !pack snow mean
               kpds(5)=143


               var = ptype_mn(:,lv,2)
               call putgbex(30,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)             
              
             else                     !non-precip type variables

               kens(4) = 2
               var = derv_mn(:,nv,lv)



               call putgbex(30,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
c              call print_kpds(nv,2,'_mn',kpds,kens,kprob,xprob)
             end if

           !spread:
             if(dk5(nv).ne.140) then         !dominant precip has no spread computation
               kens(4) = 11
               var = derv_sp(:,nv,lv)


               call putgbex(40,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
c             call print_kpds(nv,2,'_sp',kpds,kens,kprob,xprob)

              end if
           end do

        end if

        !accumulated precip/snow/freezing rain mean/spread and mam/min
        if(dk5(nv).eq.61.or.dk5(nv).eq.79.or.
     &                     dk5(nv).eq.78) then

           kpds(22) = 1
           kpds(5)  = dk5(nv)
           kpds(6)  = dk6(nv)
           kpds(7)  = 0
           kpds(16) = 4

           do 2031 lv = 1, dMlvl(nv)
                                                                                                                                               
            kpds(14) = itime - dMeanLevel(nv,lv)
            if(kpds(14).lt.0) goto 2031                    !skip 
            kpds(15) = itime
           
            !mean  
            kens(4) = 2
            var = derv_mn(:,nv,lv)


            call putgbex(30,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
c            call print_kpds(nv,2,'_mn',kpds,kens,kprob,xprob)
 
            !spread:
             kens(4) = 11
             var = derv_sp(:,nv,lv) 


             call putgbex(40,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
c            call print_kpds(nv,2,'_sp',kpds,kens,kprob,xprob)      

            !maximum
             kens(4) = 21
             if(dk5(nv).eq.61) var=prcpmax(:,lv)
             if(dk5(nv).eq.79) var=snowmax(:,lv)
             call putgbex(60,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
             call print_kpds(nv,2,'_mn',kpds,kens,kprob,xprob)

            !minimum
             kens(4) = 22
             if(dk5(nv).eq.61) var=prcpmin(:,lv)
             if(dk5(nv).eq.79) var=snowmin(:,lv)
             call putgbex(70,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
             call print_kpds(nv,2,'_mn',kpds,kens,kprob,xprob)
2031       continue 
      
         end if

       end if
203    continue

cccccc 4: write prob for derived variables into grib files ccccccccccc
     
       write(*,*) 'In pack derived prob ... ','nderiv=',nderiv                                                                                                                                           
       do 204 nv = 1, nderiv
        write(*,*)'nv=',nv, 'dk5=', dk5(nv)

         !non-accumulated derived var prob
         if(dk5(nv).ne.61.and.dk5(nv).ne.79.and.
     &      dk5(nv).ne.78) then 

            kpds(14) = itime
            kpds(15) = 0
            kpds(16) = 0

          do lv = 1, dPlvl(nv)

          
          write(*,*)'lv=',lv, 'of ', dPlvl(nv)
            
           if(dk5(nv).eq.140) then  !precip type prob
             kens(4)  = -1
             kprob(1) = -1         !: OCT 46
             kprob(2) = 3         !: OCT 47
 
             kpds(6) = dk6(nv)
             kpds(7) = dProblevel(nv,lv)
             kpds(22) = 1
             xprob(1) = 1.0
             xprob(2) = 1.0

c            kpds(5)=140            
             kpds(5)=191 
             kprob(1) = 140          
             var = ptype_pr(:,lv,4)
              do i=1,jf
               if(var(i).lt.0.0) var(i)=0.0
              end do 


            write(*,*) "PrecpType at 15000:",ptype_pr(15000,lv,:)

             call putgbex(50,jf,kpds,kgds,kens,kprob,xprob,
     &                     kclust,kmembr,lb,var,iret)

c           kpds(5)=141
            kpds(5)=191        
            kprob(1) = 141    
            var = ptype_pr(:,lv,1)
              do i=1,jf
               if(var(i).lt.0.0) var(i)=0.0
              end do 

             call putgbex(50,jf,kpds,kgds,kens,kprob,xprob,
     &                     kclust,kmembr,lb,var,iret)

c           kpds(5)=142
            kpds(5)=191
            kprob(1) = 142            
            var = ptype_pr(:,lv,3) 
            
              do i=1,jf
               if(var(i).lt.0.0) var(i)=0.0
              end do 

             call putgbex(50,jf,kpds,kgds,kens,kprob,xprob,
     &                     kclust,kmembr,lb,var,iret)

c           kpds(5)=143
            kpds(5)=191
            kprob(1) = 143
 

            var = ptype_pr(:,lv,2)
              do i=1,jf
               if(var(i).lt.0.0) var(i)=0.0
              end do 

             call putgbex(50,jf,kpds,kgds,kens,kprob,xprob,
     &                     kclust,kmembr,lb,var,iret)

           else                     !other prob than precip type prob
             kens(4)  = 0
             kprob(1) = dk5(nv)
                                                                                                                                               
             if(trim(dop(nv)).eq.'>') then
               kprob(2) = 2
             else if(trim(dop(nv)).eq.'<') then
               kprob(2) = 1
             else if(trim(dop(nv)).eq.'-') then
               kprob(2) = 3
c             else if(trim(dop(nv)).eq.'=') then
c               kprob(2) = 4
             else if(trim(dop(nv)).eq.'=') then
                kprob(2) = 3
             else
               kprob(2) = -1
             end if

             kpds(5) = 191
             kpds(6) = dk6(nv)
             kpds(7) = dProblevel(nv,lv)
             kpds(22) = 1

                    
            do 102 lt = 1, dTlvl(nv)
              write(*,*) 'lt=',lt,' of', dTlvl(nv)

             if(trim(dop(nv)).eq.'>') then
                 xprob(1) = 0.0
                 xprob(2) = dThrs(nv,lt)
             else if(trim(dop(nv)).eq.'<') then
c     +               trim(dop(nv)).eq.'=') then
                 xprob(2) = 0.0
                 xprob(1) = dThrs(nv,lt)
             else if(trim(dop(nv)).eq.'-') then
                 if ( lt.LT.dTlvl(nv) ) then
                   xprob(1) = dThrs(nv,lt)
                   xprob(2) = dThrs(nv,lt+1)
                 end if
             else if(trim(dop(nv)).eq.'=') then
                 xprob(1) = dThrs(nv,lt)
                 xprob(2) = dThrs(nv,lt)
             else
                xprob(1) = -1.0
                xprob(2) = -1.0
             end if

                if(trim(dop(nv)).eq.'-'.and.lt.eq.dTlvl(nv))
     +             goto 102

              write(*,*) 'putgbex once'
              var = derv_pr(:, nv, lv, lt)
              do i=1,jf
               if(var(i).lt.0.0) var(i)=0.0
              end do

              call putgbex(50,jf,kpds,kgds,kens,kprob,xprob,
     &        kclust,kmembr,lb,var,iret)

c              call print_kpds(nv,2,'_pr',kpds,kens,kprob,xprob)
102          continue

           end if                 !end of prob than precip type prob

          end do     !end of Plevel

        end if       ! end of non-accumu precip/snow

        !accumulated precip/snow/freezing rain prob
        if(dk5(nv).eq.61.or.dk5(nv).eq.79.or.
     &                      dk5(nv).eq.78) then

          do 2041 lv = 1, dPlvl(nv)
                                                                                                                                               
             kens(4)  = 0
             kprob(1) = dk5(nv)
                                                                                                                                               
             if(trim(dop(nv)).eq.'>') then
               kprob(2) = 2
             else if(trim(dop(nv)).eq.'<') then
               kprob(2) = 1
             else if(trim(dop(nv)).eq.'-') then
               kprob(2) = 3
             else
               kprob(2) = -1
             end if
                                                                                                                                               
             kpds(5)  = 191
             kpds(6)  = dk6(nv)
             kpds(7)  = 0
             kpds(14) = itime - dProblevel(nv,lv)
             if(kpds(14).lt.0) goto 2041                    !skip
             kpds(15) = itime
             kpds(16) = 4
             kpds(22) = 1
                                                                                                                                               
            do 103 lt = 1, dTlvl(nv)
             if(trim(dop(nv)).eq.'>') then
                 xprob(1) = 0.0
                 xprob(2) = dThrs(nv,lt)
             else if(trim(dop(nv)).eq.'<') then
                 xprob(2) = 0.0
                 xprob(1) = dThrs(nv,lt)
             else if(trim(dop(nv)).eq.'-') then
                if ( lt.LT.dTlvl(nv) ) then
                   xprob(1) = dThrs(nv,lt)
                   xprob(2) = dThrs(nv,lt+1)
                 end if
             else
                xprob(1) = -1.0
                xprob(2) = -1.0
             end if
                  
             if(trim(dop(nv)).eq.'-'.and.lt.eq.dTlvl(nv))
     +          goto 103
                 var = derv_pr(:, nv, lv, lt)
                 do i=1,jf
                  if(var(i).lt.0.0) var(i)=0.0
                 end do

                 call putgbex(50,jf,kpds,kgds,kens,kprob,xprob,
     &           kclust,kmembr,lb,var,iret)
c                 call print_kpds(nv,2,'_pr',kpds,kens,kprob,xprob)
103        continue        

2041       continue     

        end if               !end of accumu precip/snow

204    continue


cccccc 5: write DTRA variances into spread files ccccccccccc
       
        kprob(1) = -1         !: OCT 46
        kprob(2) = -1         !: OCT 47
        xprob(1) = -1.        !: OCT 48-51
        xprob(2) = -1.        !: OCT 52-55
                                                                                                                              
      do 205 nv = 1, nderiv
       if(dk5(nv).eq.501) then
         kpds(14) = itime
         kpds(15) = 0
         kpds(16) = 0

         kpds(6) = dk6(nv)

           do lv = 1, dMlvl(nv)
                                                                                                                              
            kpds(7)  = dMeanLevel(nv,lv)
            kens(4) = 2
                           
            !pack (1) UU
            kpds(5)=203                 
            kpds(22) = 2
            var = derv_dtra(:,lv,1)
            call putgbex(40,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
          
            !pack (2) VV
            kpds(5)=204                 
            kpds(22) = 2
            var = derv_dtra(:,lv,2)
            call putgbex(40,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)


            !pack (3) UV-corr (i.e. Cov(U,V))
            kpds(5)=205
            kpds(22) = 2
            var(:) = derv_dtra(:,lv,3)*                       !UV-corr-coeff = Cov(U,V)/(sigU*sigV)
     &             derv_dtra(:,lv,4)*derv_dtra(:,lv,5)        !so, Cov(U,V) = UV-corr-coeff*sigU*sigV
            call putgbex(40,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
 
            !pack (4) SigmaU
            kpds(5)=33
            kpds(22) = 2
            kens(4) = 11
            var = derv_dtra(:,lv,4)
            call putgbex(40,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)
        
            !pack (5) SigmaV
            kpds(5)=34
            kpds(22) = 2
            var = derv_dtra(:,lv,5)
            call putgbex(40,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)

            !pack (6) UV corr-coef
            kpds(5)=205                 
            kpds(22) = 2
            var = derv_dtra(:,lv,3)
            call putgbex(30,jf,kpds,kgds,kens,kprob,xprob,
     &             kclust,kmembr,lb,var,iret)

 
          end do

       end if
205   continue

      call baclose(30,ierr)
      call baclose(40,ierr)
      call baclose(50,ierr)

      return
      end

