cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     
c     subroutine firewx: compute fire weather danger probability 
c     Fireweather prob: 1hr apcp <= 0.01(in) and RH2m <= 15% and 10m Wind speed >= 20MPH and T2m >= 60.0 F    
c
c     Author: Binbin Zhou, May 17, 2011
c     Modification: 
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
   	subroutine get_firewx (nv,rawdata,jf,iens,derv_pr)

            include 'parm.inc'

C for variable table:
        Integer numvar
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
                                                                                                                                                                

        INTEGER, intent(IN) :: nv, jf, iens
        REAL,dimension(jf,iens,numvar,maxmlvl),intent(IN) :: rawdata
        REAL,dimension(jf,nderiv,maxplvl,maxtlvl),intent(INOUT) :: 
     +               derv_pr

        real FWapnt(iens),T2mapnt(iens),RH2mapnt(iens),
     +     APCPapnt(iens),U10apnt(iens),V10apnt(iens),W10apnt(iens)
        real firewx(jf,iens),aprob

        ID_T2m  = index_table(k5,k6,11,105,maxvar)   !search index of total cloud in the table
        ID_RH2m = index_table(k5,k6,52,105,maxvar)      !search index of cloud base heightin the table
        ID_APCP = index_table(k5,k6,61,1,maxvar)      !search index of cloud cloud heightin the table
        ID_U10  = index_table(k5,k6,33,105,maxvar)      !search index of surface height in the table
        ID_V10  = index_table(k5,k6,34,105,maxvar)  !search index of cloud LWC in the table

        write(*,*) 'In firewx -->'


      if (ID_T2m .gt.0 .and. ID_RH2m .gt.0 .and.
     +    ID_APCP .gt.0 .and. ID_U10 .gt.0 .and.
     +    ID_V10  .gt.0 ) then

        do igrid = 1,jf

          T2mapnt=rawdata(igrid,:,ID_T2m,1)
          RH2mapnt=rawdata(igrid,:,ID_RH2m,1)
          APCPapnt=rawdata(igrid,:,ID_APCP,1)
          U10apnt=rawdata(igrid,:,ID_U10,1)
          V10apnt=rawdata(igrid,:,ID_V10,1)

          FWapnt = 0.

          do i = 1, iens

           W10apnt(i)=sqrt(U10apnt(i)*U10apnt(i) + 
     +      V10apnt(i)*V10apnt(i))
      
           TF=(T2mapnt(i)-273.15)*1.8 + 32.0        !from K to F 

            if(APCPapnt(i).le.0.25.and.RH2mapnt(i).le.15.0.and.
     +        TF.ge.60.0.and.W10apnt(i).ge.8.9) then
              FWapnt(i) = 1.0
            else
              FWapnt(i) = 0.0  
            end if

           end do

            firewx(igrid,:)=FWapnt(:)

            if(igrid.eq.2051) then
              write(*,*) 'firewx',(firewx(igrid,k),k=1,iens)
            end if
                     
        end do


           do 30 lv=1,dPlvl(nv)

             derv_pr(:,nv,lv,:)=0.0

             do lt = 1, dTlvl(nv)
 
              thr1 = dThrs(nv,lt)
                             
              do igrid = 1,jf

               FWapnt=firewx(igrid,:)

               call getprob(FWapnt,iens,
     +                thr1,thr2,dop(nv),aprob)
          
               derv_pr(igrid,nv,lv,lt)=aprob 
        
               if(igrid.eq.2051) then
                write(*,*) 'Thre', dop(nv),dThrs(nv,lt),':',
     +                      derv_pr(igrid,nv,lv,lt)
               end if

              end do
             end do

30          continue  

        end if

        return
        end
