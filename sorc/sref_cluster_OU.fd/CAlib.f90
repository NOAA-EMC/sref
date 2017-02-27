! ******************************************************************
! File Name: CAlib.f
! Author: Ahmad Alhamed 
!         University of Oklahoma
!         School of Computer Science
! Date: 2000
! Modified by: Nusrat Yussouf
!              CIMMS/OU/NSSL
! Date: September, 2002
!
! Purpose: This files contains several subroutines and functions 
! thar are used by other programs.
! ******************************************************************
!
!
!
! ==================================================================
! Subroutine Name: ReadDataMatrixFromFile
! Author: Ahmad Alhamad
! Date: 2000
! ==================================================================
       SUBROUTINE ReadDataMatrixFromFile(Row,Col,DataMatrix,InputFile,jt) 


        IMPLICIT NONE

          integer ntime ,nm, nfile,jf,igx,igy,ngrid
            parameter (nm=26,nfile=nm+1,ntime=30)
           parameter (jf=185*129,igx=185,igy=129,ngrid=igx*igy)

        INTEGER, INTENT(IN) :: Row, Col,jt
	CHARACTER(LEN=*), INTENT(IN) :: InputFile
!        REAL, DIMENSION(ntime,Row,Col), INTENT(OUT) :: DataMatrix 
         REAL, DIMENSION(Row,Col), INTENT(OUT) :: DataMatrix 
	INTEGER I, J, IERR, JERR
         integer kf,k,lugb,lugi
          character*13 fname1
          character*15 fname2
          

         integer nt

        real f(jf)
!     integer kftime(nfile,ntime),jpds(25),jgds(22),kpds(25),kgds(22)
      integer kftime(nfile,ntime),jpds(200),jgds(200),kpds(200),kgds(200)
      integer ifile
      integer yy1,mm1,dd1,hh1,iar,meth
      integer yy3(ntime),mm3(ntime),dd3(ntime),hh3(ntime)
         character*6 mem (nm)
      logical lb(jf)
       integer ipt,iret
        real dmin,dmax
         common /c1/yy1,mm1,dd1,hh1,iar,meth
         common /c2/yy3,mm3,dd3,hh3,kftime,mem
         common /c3/kpds,jpds,kgds,jgds,ierr,iret
         common /c4/lb
          DataMatrix=0.0 
!          do 10 nt =jt,jt
        nt=jt
        print *,' '
        print *,'----------------------------------------------- '
        print *,'NEW TIME: nt = ',nt
        print *,'Corresponding forecast hour: nt = ',kftime(1,nt)
        print *,' '
           
           j=0
        lugb=9
        lugi=39

! part I: Geting data ready:
!CCCCCCCCCCCCCCCCCCCCCCCCCCC
        do 5 ifile=1,nfile-1
! be careful of year and month when in the turn of time
        jpds=-1
        jgds=-1
!       jpds(1)=7      !center
!       jpds(2)=111    !model
!       jpds(3)=212    !grid
        jpds(4)=128    !table
 
! specify date information here:
       if(ifile.le.nm) then
!       jpds(8)=yy1
!       jpds(9)=mm1
        jpds(10)=dd1
        jpds(11)=hh1
       endif
!       if(ifile.eq.nm+1) then
!c       jpds(8)=yy3(nt)
!c       jpds(9)=mm3(nt)
!        jpds(10)=dd3(nt)
!        jpds(11)=hh3(nt)
!       endif


        jpds(14)=kftime(ifile,nt)
        jpds(15)=-1
!       jpds(15)=0

        print *,' '
        print *,'----------'
        print *,'New File: nt = ',nt,' ifile= ',ifile
        print *,'Fcst Hour = ',kftime(ifile,nt)
        print *,' '

! define variables
        if(iar.eq.1) then !500z
          jpds(7)=500
          jpds(6)=100
          jpds(5)=7
        else
          jpds(7)=0
          jpds(6)=102
          jpds(5)=2        !mslp
        endif

      print*, jpds(1)
      print*, jpds(2)
      print*, jpds(3)
      print*, jpds(4)
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
      call getgb(lugb,lugi,jf,j,jpds,jgds, &
     &             kf,k,kpds,kgds,lb,f,iret)
      call baclose(lugb,ierr)
      call baclose(lugi,ierr)
        call grange(kf,lb,f,dmin,dmax)
        print *,'immediately after call to getgb, iret= ',iret  &
     &,' j= ',j,' k= ',k,' ifile= ',ifile,' nt= ',nt  &
     &,' j1= ',jpds(1),' j2= ',jpds(2),' j3= ',jpds(3),' j4= ',jpds(4)  &
     &,' j5= ',jpds(5),' j6= ',jpds(6),' j7= ',jpds(7),' j8= ',jpds(8)  &
     &,' j9= ',jpds(9),' j10= ',jpds(10),' j11= ',jpds(11),' j12= ',jpds(12)   &
     &,' j13= ',jpds(13),' j14= ',jpds(14),' j15= ',jpds(15)    &
     &,' j16= ',jpds(16),' j17= ',jpds(17),' j18= ',jpds(18)    &
     &,' j19= ',jpds(19),' j20= ',jpds(20),' j21= ',jpds(21)    &
     &,' j22= ',jpds(22),' j23= ',jpds(23),' j24= ',jpds(24)    &
     &,' j25= ',jpds(25)
      print '(i4,2x,25i5,i8,2g12.4)',k,(kpds(i),i=1,25),kf,dmin,dmax

      if (iret.eq.0) then
       do ipt=1,ngrid
        DataMatrix(ipt,ifile) = f(ipt)
       enddo
      else
       print *, "Something wrong here!!!!"
      endif

!       print*, 'original data:',ifile, data(1,ifile)
5       continue    !files (mems, climo)
        print*, 'ok till here 1'
!10     continue

END SUBROUTINE ReadDataMatrixFromFile

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

END SUBROUTINE grange


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
!      return
END SUBROUTINE getname







! =======================================================================
! Subroutine Name: ReadResemblanceMatrixFromFile
! Author: Ahmad Alhamad
! Date: 2000
! =======================================================================
SUBROUTINE ReadResemblanceMatrixFromFile(N, ResemblanceMatrix, InputFile)
        IMPLICIT NONE            

        INTEGER, INTENT(IN) :: N
	CHARACTER(LEN=*), INTENT(IN) :: InputFile
        REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: ResemblanceMatrix
	INTEGER I, IERR, JERR

	OPEN (UNIT=9,FILE=InputFile,STATUS='OLD',IOSTAT=IERR)
        IF (IERR /= 0) THEN
                PRINT 8 , InputFile, IERR
8               FORMAT ( ' ERROR OPENING UNIT=9, FILE NAME = ',A15,&
                &  ', IOSTAT = ',I8)
	        STOP 8
        ENDIF
	READ (9, *) (ResemblanceMatrix(I), I=1, N*(N-1)/2 )
        CLOSE (UNIT=9) 
END SUBROUTINE ReadResemblanceMatrixFromFile
! ======================================================================
! Subroutine Name: PrintDataMatrix
! Author: Ahmad Alhamad
! Date: 2000
! ======================================================================
SUBROUTINE PrintDataMatrix(Row, Col, DataMatrix)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: Row, Col
        REAL, DIMENSION(Row,Col), INTENT(IN) :: DataMatrix
	INTEGER :: I, J

	WRITE (7, 10) Row, Col
	WRITE (7, *) '------------------------------'
10	FORMAT ('    DATA MATRIX (',I3,' x ',I3,')' )
	! print col and row number
	WRITE (7, '($,A)' ) '     '
	WRITE (7, 20) (I, I=1, Col)
20	FORMAT ('     [',I2,'] ',$)
	WRITE (7, *)		
        DO I=1, Row
		WRITE (7, 25) I
		WRITE (7,FMT='(50F10.2)') (DataMatrix(I,J),J=1,Col)
        END DO
25	FORMAT ('[',I2,'] ',$)
	WRITE (7, *)
END SUBROUTINE PrintDataMatrix
! =======================================================================
! Subroutine Name: PrintResemblanceMatrix
! Author: Ahmad Alhamad
! Date: 2000
! =======================================================================
SUBROUTINE PrintResemblanceMatrix(N, ResemblanceMatrix)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N
        REAL, DIMENSION(N*(N-1)/2), INTENT(IN) :: ResemblanceMatrix
	INTEGER :: Offset    ! Func. returns index in Resemblance Matrix
	INTEGER :: I, J

        WRITE (7, '($,A)' ) '     '
        WRITE (7,30) (I, I=1, N-1)
30      FORMAT ('       [',I2,'] ',$)
        WRITE (7, *)
	DO I=2, N
		WRITE (7, 35) I
		WRITE (7,FMT='(53F12.3)') &
		& ( ResemblanceMatrix(Offset(N,I,J)), J=1, I-1)
	END DO
35      FORMAT ('[',I2,'] ',$)
	WRITE (7, *)
END SUBROUTINE PrintResemblanceMatrix
! =================================================================================
! Subroutine Name: NormalizeMenu
! Author: Ahmad Alhamad
! Date: 2000
! Modified By: Nusrat Yussouf
! Date: September, 2002
! =================================================================================
INTEGER FUNCTION NormalizeMenu()
        IMPLICIT NONE

        INTEGER :: Norm

        PRINT *
        PRINT *,'---------------------------------------------------------'
        PRINT *,'Normalization:'
        PRINT *,'    (0) Do Not Normalize'
        PRINT *,'    (1) Centering (by the Mean of each Column)'
        PRINT *,'    (2) Normalize:Range zero-one(by the Min. and Max. of each Column)'
        PRINT *,'    (3) Standard Normalize (zero mean and Unit variance)'
        PRINT *,'----------------------------------------------------------'
        PRINT '($,A)','    Enter Choice ===>   '
        DO
                READ *, Norm
                IF (Norm >= 0 .AND. Norm <= 3)  EXIT
        END DO
        NormalizeMenu = Norm
END FUNCTION NormalizeMenu 
! ========================================================================
! Subroutine Name: CenterByColMean
! Author: Ahmad Alhamad
! Date: 2000 
! ========================================================================
SUBROUTINE CenterByColMean(M, N, X)
        IMPLICIT NONE            
	! ********************************************
	! Center matrix X by the mean of each col
	! Zij = Xij - MEANj   , 1<=I<=M  , 1<=J<=N
	! ********************************************

	INTEGER, INTENT(IN) :: M, N
	REAL, DIMENSION(M, N), INTENT(INOUT) :: X
	INTEGER :: I, J, K, AllocatStatus
	REAL, ALLOCATABLE, DIMENSION(:)   :: Mean

	ALLOCATE( Mean(N) , STAT = AllocatStatus )
	IF (AllocatStatus/= 0 )  STOP  "( Mean ) ** No Enough Memory **"
	! comput Mean of each col.
	DO J=1, N
		Mean(J) = SUM( X(1:M,J) ) / REAL(M)
	END DO
	! subtract Mean of each col from every element in this col.
	DO J=1, N
		X(1:M,J) = X(1:M,J) - Mean(J)
	END DO
	DEALLOCATE( Mean, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "( Mean )** Deallocation Error**"
END SUBROUTINE CenterByColMean
! ==========================================================================
! Subroutine Name: StdNormalize
! Author: Nusrat Yussouf
! Date: September, 2002
! ==========================================================================
SUBROUTINE StdNormalize(M, N, X)
        IMPLICIT NONE            
	! ********************************************
	! Standard Normalization
	! Zij = (Xij - MEANj)/std.deviation, 1<=I<=M  , 1<=J<=N
	! ********************************************

	INTEGER, INTENT(IN) :: M, N
	REAL, DIMENSION(M, N), INTENT(INOUT) :: X
	INTEGER :: I, J, K, AllocatStatus
	REAL, ALLOCATABLE, DIMENSION(:)   :: Mean, Var, Std
        REAL :: Unbias
	ALLOCATE( Mean(N) , STAT = AllocatStatus )
	IF (AllocatStatus/= 0 )  STOP  "( Mean ) ** Not Enough Memory **"
	ALLOCATE( Var(N) , STAT = AllocatStatus )
	IF (AllocatStatus/= 0 )  STOP  "( Var ) ** Not Enough Memory **"
	ALLOCATE( Std(N) , STAT = AllocatStatus )
	IF (AllocatStatus/= 0 )  STOP  "( Std ) ** Not Enough Memory **"
	! compute Mean of each col.
	DO J=1, N
		Mean(J) = SUM( X(1:M,J) ) / REAL(M)
	END DO
	! subtract Mean of each col from every element in this col.
	DO J=1, N
		X(1:M,J) = X(1:M,J) - Mean(J)
	END DO
	! compute unbiased variance of each col. 
	DO J=1, N
		Var(J) = DOT_PRODUCT( X(1:M,J), X(1:M,J) )		
	END DO	
	Unbias = REAL(M - 1)
	Unbias = REAL(1/Unbias)
	DO J=1, N
		Var(J) =  Unbias * Var(J)
	END DO
	! Compute the standard deviation
	DO J=1, N	         
		X(1:M,J) = X(1:M,J) / SQRT( Var(J) )
	END DO
	DEALLOCATE( Mean, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "( Mean )** Deallocation Error**"
	DEALLOCATE( Var, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "( Var )** Deallocation Error**"
	DEALLOCATE( Std, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "( Std )** Deallocation Error**"
END SUBROUTINE StdNormalize
! ===========================================================================
! Subroutine Name: NormalizeByColMax
! Author: Ahmad Alhamad
! Date: 2000
! ===========================================================================
SUBROUTINE NormalizeByColMax(M, N, X)
        IMPLICIT NONE            
	! ********************************************
	! Normalize matrix X by the max of each col
	! Zij = Xij/MAXj  , 1<=I<=M , 1<=J<=N
	! ********************************************

	INTEGER, INTENT(IN) :: M, N
	REAL, DIMENSION(M, N), INTENT(INOUT) :: X
	INTEGER :: I, J, K, AllocatStatus
	REAL, ALLOCATABLE, DIMENSION(:)   :: Mx

	ALLOCATE( Mx(N), STAT = AllocatStatus )
	IF (AllocatStatus/= 0 )  STOP  "( Mx ) ** No Enough Memory **"
	! comput Max of each col.
	DO J=1, N
		Mx(J) = MAXVAL( X(1:M,J) )
	END DO
	! divede every element of each col by its Max.
	DO J=1, N
		X(1:M,J) = X(1:M,J) / Mx(J)
	END DO
	DEALLOCATE( Mx, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "( Mx )** Deallocation Error**"
END SUBROUTINE NormalizeByColMax
! ============================================================================
! Subroutine Name: NormalizeByColMinMax
! Author: Ahmad Alhamad
! Date: 2000
! ============================================================================
SUBROUTINE NormalizeByColMinMax(M, N, X)
        IMPLICIT NONE            
	! ********************************************
	! Normalize matrix X by the min and max of 
	! each col Zij=(Xij-MINj)/(MAXj-MINj) ,1<=i<=M
	! ********************************************

	INTEGER, INTENT(IN) :: M, N
	REAL, DIMENSION(M, N), INTENT(INOUT) :: X
	INTEGER :: I, J, K, AllocatStatus
	REAL, ALLOCATABLE, DIMENSION(:)   :: Mn, Mx, Av

	ALLOCATE( Mn(N), STAT = AllocatStatus )
	IF (AllocatStatus/= 0 )  STOP  "( Mn ) ** No Enough Memory **"
	ALLOCATE( Mx(N), STAT = AllocatStatus )
	IF (AllocatStatus/= 0 )  STOP  "( Mx ) ** No Enough Memory **"
	! comput Min, Max of each col.
	DO J=1, N
		Mn(J) = MINVAL( X(1:M,J) )
		Mx(J) = MAXVAL( X(1:M,J) )
	END DO
	! normalize by Min, Max of each col.
	DO J=1, N
		X(1:M,J) = ( X(1:M,J)-Mn(J) ) / ( Mx(J)-Mn(J) )
	END DO
	DEALLOCATE( Mn, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "( Mn )** Deallocation Error**"
	DEALLOCATE( Mx, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "( Mx )** Deallocation Error**"
END SUBROUTINE NormalizeByColMinMax
! ===========================================================================
! Subroutine Name: NormalizeByMax
! Author: Ahmad Alhamad
! Date: 2000
! ===========================================================================
SUBROUTINE NormalizeByMax(M, N, X)
        IMPLICIT NONE            
	! ********************************************
	! Normalize matrix X by the max of X
	!  Zij = Xij/MAX , 1<=I<=M, 1<=J<=N
	! ********************************************

	INTEGER, INTENT(IN) :: M, N
	REAL, DIMENSION(M, N), INTENT(INOUT) :: X
	INTEGER :: I, J, K, AllocatStatus
	REAL :: Mx

	Mx = MAXVAL(X)
	X = X / Mx
END SUBROUTINE NormalizeByMax
! =============================================================================
! Subroutine Name:CalculateCentroid
! Author: Ahmad Alhamad
! Date: 2000
! =============================================================================
SUBROUTINE CalculateCentroid(M,N,NoOfVectors,Indecis,X, Centroid)
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: M, N, NoOfVectors
	INTEGER, DIMENSION(NoOfVectors), INTENT(IN) :: Indecis
	REAL, DIMENSION(N,M), INTENT(IN) :: X
	REAL, DIMENSION(M), INTENT(OUT) :: Centroid
	INTEGER :: I

	Centroid = 0.0
	DO I=1,NoOfVectors
		Centroid = Centroid + X( Indecis(I), 1:M)
	END DO
	Centroid = Centroid / REAL(NoOfVectors)
END SUBROUTINE CalculateCentroid
! ===============================================================================
! Subroutine Name: IsIn
! Author: Ahmad Alhamad
! Date: 2000
! ===============================================================================
LOGICAL FUNCTION IsIn(Num,Size,Array)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: Num, Size
        INTEGER, DIMENSION(Size), INTENT(IN) :: Array
        INTEGER :: I
        LOGICAL :: Found

        Found = .FALSE.
        DO I=1, Size
                IF(Array(I) == Num) THEN
                        Found = .TRUE.
                        EXIT
                END IF
        END DO
        IsIn = Found
END FUNCTION IsIn 
! ==============================================================================
! Subroutine Name: CalculateCenters
! Author: Ahmad Alhamad
! Date: 2000
! ==============================================================================
SUBROUTINE CalculateCenters(M, N, X, K, Membership, NC, Centers)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N, K
        REAL, DIMENSION(M,N), INTENT(IN) :: X
        INTEGER, DIMENSION(N), INTENT(IN) :: Membership
        INTEGER, DIMENSION(K), INTENT(IN) :: NC
        REAL, DIMENSION(M,K), INTENT(OUT) :: Centers
	INTEGER :: I, J

	Centers = 0.0
	DO I=1, N
	   Centers(1:M, Membership(I)) = Centers(1:M, Membership(I)) + &
		&  X(1:M,I)
	END DO
	DO I=1, K
	   Centers(1:M, I) = Centers(1:M, I) / REAL( NC(I) )
	END DO
END SUBROUTINE CalculateCenters
! =============================================================================
! Subroutine Name: EuclideanSimilarity
! Author: Ahmad Alhamad
! Date: 2000
! =============================================================================
SUBROUTINE EuclideanSimilarity(M, N, Data, ES)
	IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N
        REAL, DIMENSION(M,N), INTENT(IN) :: Data
        REAL, DIMENSION(N,N), INTENT(OUT) :: ES
        REAL, ALLOCATABLE, DIMENSION(:) :: Vector
        INTEGER:: J, K, AllocatStatus
        REAL :: Mx

        ALLOCATE( Vector(M), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** No Enough Memory **"

        DO J=1, N-1
                DO K=J+1, N
                        Vector = Data(1:M, J) - Data(1:M, K)
                        ES(K,J) = SQRT(DOT_PRODUCT(Vector,Vector))
			ES(J,K) = ES(K,J)
                END DO
        END DO
	DO J=1, N
		ES(J,J) = 0.0
	END DO

        Mx = MAXVAL(ES)
        ES = ES / Mx

       DO J=1, N
               DO K=1, N
                       ES(J,K) = 1.0 - ES(J,K)
               END DO
       END DO

        DEALLOCATE( Vector, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** Deallocation Error **"

END SUBROUTINE EuclideanSimilarity
! =======================================================================

