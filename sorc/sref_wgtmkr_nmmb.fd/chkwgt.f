      PROGRAM WGTMKR
C    
      PARAMETER(IMAXOT=100,JJMAXOT=100,JMAXOT=IMAXOT*JJMAXOT)

      INTEGER KGDSIN(22),KGDSOUT2(22),KGDSOUT(22)

      INTEGER N11(JMAXOT),N21(JMAXOT),
     &              N12(JMAXOT),N22(JMAXOT),
     &              NPP(JMAXOT,25)
      INTEGER NV11(JMAXOT),NV21(JMAXOT),
     &              NV12(JMAXOT),NV22(JMAXOT)

      INTEGER N112(JMAXOT),N212(JMAXOT),
     &              N122(JMAXOT),N222(JMAXOT),
     &              NPP2(JMAXOT,25)
      INTEGER NV112(JMAXOT),NV212(JMAXOT),
     &              NV122(JMAXOT),NV222(JMAXOT)

      REAL RLAT(JMAXOT),RLON(JMAXOT)
      REAL CROT(JMAXOT),SROT(JMAXOT)
      REAL W11(JMAXOT),W21(JMAXOT),
     &     W12(JMAXOT),W22(JMAXOT)
      REAL WV11(JMAXOT),WV21(JMAXOT),
     &     WV12(JMAXOT),WV22(JMAXOT)
      REAL C11(JMAXOT),C21(JMAXOT),
     &     C12(JMAXOT),C22(JMAXOT)
      REAL S11(JMAXOT),S21(JMAXOT),
     &     S12(JMAXOT),S22(JMAXOT)
      REAL RLAT2(JMAXOT),RLON2(JMAXOT)
      REAL CROT2(JMAXOT),SROT2(JMAXOT)
      REAL W112(JMAXOT),W212(JMAXOT),
     &     W122(JMAXOT),W222(JMAXOT)
      REAL WV112(JMAXOT),WV212(JMAXOT),
     &     WV122(JMAXOT),WV222(JMAXOT)
      REAL C112(JMAXOT),C212(JMAXOT),
     &     C122(JMAXOT),C222(JMAXOT)
      REAL S112(JMAXOT),S212(JMAXOT),
     &     S122(JMAXOT),S222(JMAXOT)

      CHARACTER GDSO(400)

       LUNOUT=51

C
         READ(LUNOUT) KGRIDOT,NOUT
         READ(LUNOUT) (KGDSOUT(I),I=1,22)
         READ(LUNOUT) (N11(I),I=1,NOUT)
         READ(LUNOUT) (N12(I),I=1,NOUT)
         READ(LUNOUT) (N21(I),I=1,NOUT)
         READ(LUNOUT) (N22(I),I=1,NOUT)
         READ(LUNOUT) (NV11(I),I=1,NOUT)
         READ(LUNOUT) (NV12(I),I=1,NOUT)
         READ(LUNOUT) (NV21(I),I=1,NOUT)
         READ(LUNOUT) (NV22(I),I=1,NOUT)
         READ(LUNOUT) (C11(I),I=1,NOUT)
         READ(LUNOUT) (C12(I),I=1,NOUT)
         READ(LUNOUT) (C21(I),I=1,NOUT)
         READ(LUNOUT) (C22(I),I=1,NOUT)
         READ(LUNOUT) (S11(I),I=1,NOUT)
         READ(LUNOUT) (S12(I),I=1,NOUT)
         READ(LUNOUT) (S21(I),I=1,NOUT)
         READ(LUNOUT) (S22(I),I=1,NOUT)
         READ(LUNOUT) (W11(I),I=1,NOUT)
         READ(LUNOUT) (W12(I),I=1,NOUT)
         READ(LUNOUT) (W21(I),I=1,NOUT)
         READ(LUNOUT) (W22(I),I=1,NOUT)
         READ(LUNOUT) (WV11(I),I=1,NOUT)
         READ(LUNOUT) (WV12(I),I=1,NOUT)
         READ(LUNOUT) (WV21(I),I=1,NOUT)
         READ(LUNOUT) (WV22(I),I=1,NOUT)
         READ(LUNOUT) (RLAT(I),I=1,NOUT)
         READ(LUNOUT) (RLON(I),I=1,NOUT)
         READ(LUNOUT) (SROT(I),I=1,NOUT)
         READ(LUNOUT) (CROT(I),I=1,NOUT)
         READ(LUNOUT) ((NPP(I,J),I=1,NOUT),J=1,25)
C
       LUNOUT=52

C
         READ(LUNOUT) KGRIDOT2,NOUT2
         write(6,*) ' kgrid ',kgridot,kgridot2
         write(6,*) ' nout ',nout,nout2
         READ(LUNOUT) (KGDSOUT2(I),I=1,22)
         write(6,*) ' kgds1 ',kgdsout
         write(6,*) ' kgds2 ',kgdsout2
         READ(LUNOUT) (N112(I),I=1,NOUT)
        call diffi('n11',n11,n112,nout)
         READ(LUNOUT) (N122(I),I=1,NOUT)
        call diffi('n12',n12,n122,nout)
         READ(LUNOUT) (N212(I),I=1,NOUT)
        call diffi('n21',n21,n212,nout)
         READ(LUNOUT) (N222(I),I=1,NOUT)
        call diffi('n22',n22,n222,nout)
         READ(LUNOUT) (NV112(I),I=1,NOUT)
        call diffi('nv11',nv11,nv112,nout)
         READ(LUNOUT) (NV122(I),I=1,NOUT)
        call diffi('nv12',nv12,nv122,nout)
         READ(LUNOUT) (NV212(I),I=1,NOUT)
        call diffi('nv21',nv21,nv212,nout)
         READ(LUNOUT) (NV222(I),I=1,NOUT)
        call diffi('nv22',nv22,nv222,nout)
         READ(LUNOUT) (C112(I),I=1,NOUT)
        call diffr('c11',c11,c112,nout)
         READ(LUNOUT) (C122(I),I=1,NOUT)
        call diffr('c12',c12,c122,nout)
         READ(LUNOUT) (C212(I),I=1,NOUT)
        call diffr('c21',c21,c212,nout)
         READ(LUNOUT) (C222(I),I=1,NOUT)
        call diffr('c22',c22,c222,nout)
         READ(LUNOUT) (S112(I),I=1,NOUT)
        call diffr('s11',s11,s112,nout)
         READ(LUNOUT) (S122(I),I=1,NOUT)
        call diffr('s12',s12,s122,nout)
         READ(LUNOUT) (S212(I),I=1,NOUT)
        call diffr('s21',s21,s212,nout)
         READ(LUNOUT) (S222(I),I=1,NOUT)
        call diffr('s22',s22,s222,nout)
         READ(LUNOUT) (W112(I),I=1,NOUT)
        call diffr('w11',w11,w112,nout)
         READ(LUNOUT) (W122(I),I=1,NOUT)
        call diffr('w12',w12,w122,nout)
         READ(LUNOUT) (W212(I),I=1,NOUT)
        call diffr('w21',w21,w212,nout)
         READ(LUNOUT) (W222(I),I=1,NOUT)
        call diffr('w22',w22,w222,nout)
         READ(LUNOUT) (WV112(I),I=1,NOUT)
        call diffr('wv11',wv11,wv112,nout)
         READ(LUNOUT) (WV122(I),I=1,NOUT)
        call diffr('wv12',wv12,wv122,nout)
         READ(LUNOUT) (WV212(I),I=1,NOUT)
        call diffr('wv21',wv21,wv212,nout)
         READ(LUNOUT) (WV222(I),I=1,NOUT)
        call diffr('wv22',wv22,wv222,nout)
         READ(LUNOUT) (RLAT2(I),I=1,NOUT)
        call diffr('rlat',rlat,rlat2,nout)
         READ(LUNOUT) (RLON2(I),I=1,NOUT)
        call diffr('rlon',rlon,rlon2,nout)
         READ(LUNOUT) (SROT2(I),I=1,NOUT)
        call diffr('srot',srot,srot2,nout)
         READ(LUNOUT) (CROT2(I),I=1,NOUT)
        call diffr('crot',crot,crot2,nout)
         READ(LUNOUT) ((NPP2(I,J),I=1,NOUT),J=1,25)
        do j=1,25
        write(6,*) ' j= ',j
        call diffi('npp',npp(1,J),npp2(1,J),nout)
        enddo
C
      STOP1
      END
        subroutine diffi(fld,n11,n112,nout)
        dimension n11(nout),n112(nout)
        character fld*4
          do I=1,nout
           if(n11(i).ne.n112(i))  then
           write(6,*) fld,' = ',n11(i),fld,'2= ',n112(i),' at i= ',i
           endif
          enddo
          return
           end
        subroutine diffr(fld,w11,w112,nout)
        dimension w11(nout),w112(nout)
        character fld*4
          rmax=0.
          do I=1,nout
           diff=abs(w11(i)-w112(i))
           if(diff.gt.rmax)  then
            rmax=diff
            imax=i
           endif
          enddo
           write(6,*) fld,' max diff = ',rmax,' at i= ',imax
           write(6,*) fld,' = ',w11(imax),fld,'2= ',w112(imax)
          return
           end
