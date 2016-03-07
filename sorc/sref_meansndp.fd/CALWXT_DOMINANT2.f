       SUBROUTINE CALWXT_DOMINANT2(PREC,TOTSN,TOTIP,TOTZR,
     &              TOTR,SNOW,SLEET,FZR,RAIN) 
C
C     WRITTEN: 5 APRIL 2006, G MANIKIN 
C      
C     THIS ROUTINE TAKES THE PRECIP TYPE SOLUTIONS FROM DIFFERENT
C       ALGORITHMS AND DETERMINES THE MAJORITY TO GIVE A DOMINANT TYPE
C
      PARAMETER (PTHRESH=0.02)
      REAL PREC,TOTSN,TOTIP,TOTZR,TOTR
      REAL SNOW,SLEET,FZR,RAIN
C
        RAIN = 0.
        SLEET = 0.
        FZR = 0.
        SNOW = 0.
C
C   SKIP THIS POINT IF NO PRECIP THIS TIME STEP
       IF (PREC .LE. PTHRESH) RETURN 


C   TIES ARE BROKEN TO FAVOR THE MOST DANGEROUS FORM OF PRECIP
C     FREEZING RAIN > SNOW > SLEET > RAIN 
        IF (TOTSN .GT. TOTIP) THEN
         IF (TOTSN .GT. TOTZR) THEN
          IF (TOTSN .GE. TOTR) THEN
           SNOW = 1.
           GOTO 800 
          ELSE
           RAIN = 1. 
           GOTO 800 
          ENDIF
         ELSE IF (TOTZR .GE. TOTR) THEN
          FZR = 1.
          GOTO 800 
         ELSE
          RAIN = 1.
          GOTO 800 
         ENDIF 
        ELSE IF (TOTIP .GT. TOTZR) THEN
         IF (TOTIP .GE. TOTR) THEN
          SLEET = 1.
          GOTO 800 
         ELSE
           RAIN= 1.
          GOTO 800 
         ENDIF
        ELSE IF (TOTZR .GE. TOTR) THEN
         FZR = 1.
         GOTO 800 
         ELSE
           RAIN= 1.
          GOTO 800 
         ENDIF
 800   RETURN
      END
