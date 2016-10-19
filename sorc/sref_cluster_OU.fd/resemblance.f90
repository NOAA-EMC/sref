! ***********************************************************************
! File Name: resemblance.f
! Author: Ahmad Alhamed 
!         University of Oklahoma
!         School of Computer Science
! Date: 2000
! Modified by: Nusrat Yussouf
!              CIMMS/OU/NSSL
! Date: September, 2002
! Purpose:
! Calculates resemblance matrix using different coefficients
! Input:
!	Coef  :  Coefficient of resemblance
!	     (1)  Euclidean Distance
!	     (2)  SquareEuclidean Distance
!	     (3)  Manhattan Distance
!	     (4)  Sup Distance
!	     (5)  Covariance Coefficient	
!	     (6)  Correlation Coefficient
!	     (7)  Euclidean Similarity
!	     (8)  Cosine Coefficient	
!	M     :  No. of attributes (rows) in the data matrix
!	N     :  No. of objects (columns) in the data matrix
!	Data  :  MxN data matrix (2D array)
! Output
!	Resemblance : 1D array contains the lower traingle of 
!		      resemblance matrix in a col-major order
! ******************************************************************
SUBROUTINE CalculateResemblanceMatrix(Coef, M, N, Data, Resemblance)
	USE Global
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: Coef, M, N
	REAL, DIMENSION(M,N), INTENT(IN) :: Data
	REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: Resemblance
	CHARACTER(LEN=50) :: CoefName

	SELECT CASE(Coef)
	CASE (1)
		CALL Euclidean(M, N, Data, Resemblance)
		CoefName = 'Euclidean Distance'
	CASE (2)
		CALL SquareEuclidean(M, N, Data, Resemblance)
		CoefName = 'SquareEuclidean Distance'
	CASE (3)
		CALL Manhattan(M, N, Data, Resemblance)
		CoefName = 'Manhattan Distance'
	CASE (4)
		CALL SupDistance(M, N, Data, Resemblance)
		CoefName = 'Sup Distance'
	CASE (5)
		CALL CovarianceCoef(M, N, Data, Resemblance)
		CoefName = 'Covariance Coefficient'	
	CASE (6)
		CALL CorrelationCoef(M, N, Data, Resemblance)
		CoefName = 'Correlation Coefficient'
	CASE (7)
		CALL EuclideanSimiCoef(M, N, Data, Resemblance)
		CoefName = 'Euclidean Similarity'
	CASE (8)
		CALL CosineCoef(M, N, Data, Resemblance)
		CoefName = 'Cosine Coefficient'	
	END SELECT
	IF (.NOT.ResWritten) THEN
		WRITE (7, 71) CoefName
		WRITE (7, *) '-------------------------------------&
		&------------------------------'
	END IF
71	FORMAT ('   Resemblance Matrix - computed using ', A )
END SUBROUTINE CalculateResemblanceMatrix
! ===========================================================
! Subroutine Name: Euclidean
! Author: Ahmad Alhamad
! Date: 2000
! ===========================================================
SUBROUTINE Euclidean(M, N, Data, Resemblance)
	IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N
        REAL, DIMENSION(M,N), INTENT(IN) :: Data
        REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: Resemblance
        REAL, ALLOCATABLE, DIMENSION(:) :: Vector
        INTEGER:: I, J, K, AllocatStatus

        ALLOCATE( Vector(M), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** No Enough Memory **"
        I = 1
        DO J=1, N-1
                DO K=J+1, N
                        Vector = Data(1:M, J) - Data(1:M, K)
                        Resemblance(I) = SQRT(DOT_PRODUCT(Vector,Vector))
                        I = I + 1
                END DO
        END DO  
        DEALLOCATE( Vector, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** Deallocation Error **"
END SUBROUTINE Euclidean
! ===========================================================
! Subroutine Name: SquareEuclidean
! Author: Ahmad Alhamad
! Date: 2000
! ===========================================================
SUBROUTINE SquareEuclidean(M, N, Data, Resemblance)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N
        REAL, DIMENSION(M,N), INTENT(IN) :: Data
        REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: Resemblance
        REAL, ALLOCATABLE, DIMENSION(:) :: Vector
        INTEGER:: I, J, K, AllocatStatus

        ALLOCATE( Vector(M), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** No Enough Memory **"
        I = 1
        DO J=1, N-1
                DO K=J+1, N  
                        Vector = Data(1:M, J) - Data(1:M, K)
                        Resemblance(I) = DOT_PRODUCT(Vector,Vector)
                        I = I + 1
                END DO
        END DO
        DEALLOCATE( Vector, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** Deallocation Error **"
END SUBROUTINE SquareEuclidean
! ===========================================================
! Subroutine Name: CorrelationCoef
! Author: Nusrat Yussouf
! Date: September, 2002 
! ===========================================================
SUBROUTINE CorrelationCoef(M, N, Data, Resemblance)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N
        REAL, DIMENSION(M,N), INTENT(IN) :: Data
        REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: Resemblance
	REAL, ALLOCATABLE, DIMENSION(:,:) :: C
        REAL, ALLOCATABLE, DIMENSION(:,:) :: X,Y
        REAL, ALLOCATABLE, DIMENSION(:)   :: Mean, Var
        INTEGER:: I, J, K, AllocatStatus
        REAL :: Unbias

        ALLOCATE( X(M,N), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "( X ) ** No Enough Memory **"
        ALLOCATE( Y(M,N), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "( Y ) ** No Enough Memory **"
        ALLOCATE( Var(N), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "( Var ) ** No Enough Memory **"     
	ALLOCATE( Mean(N) , STAT = AllocatStatus )
	IF (AllocatStatus/= 0 )  STOP  "( Mean ) ** No Enough Memory **"
	
	Unbias = REAL(M - 1)
	Unbias = REAL(1/Unbias)
	! compute Mean of each col.
	DO J=1, N
		Mean(J) = SUM( Data(1:M,J) ) / REAL(M)
	END DO
	! subtract Mean of each col from every element in this col.
	DO J=1, N
		X(1:M,J) = Data(1:M,J) - Mean(J)	
	END DO
	! calculating the variance
	DO J=1, N
	     Var(J) = Unbias*DOT_PRODUCT( X(1:M,J), X(1:M,J) )
	END DO   
	! calculating the correlation matrix
	DO J=1, N
		Y(1:M,J) = X(1:M,J) / SQRT( Var(J) )
	END DO
	I = 1
	DO J=1, N-1
           DO K=J+1, N
	     Resemblance(I) = Unbias*DOT_PRODUCT( Y(1:M,J),Y(1:M,K) )
             I = I + 1 
	   END DO
	END DO	

        DEALLOCATE( X, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "( X )** Deallocation Error**"
        DEALLOCATE( Y, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "( Y )** Deallocation Error**"
        DEALLOCATE( Var, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "( Var )** Deallocation Error**"
        DEALLOCATE( Mean, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "( Mean )** Deallocation Error**"
END SUBROUTINE CorrelationCoef
! ===========================================================
! Subroutine Name: CovarianceCoef
! Author: Nusrat Yussouf
! Date: September, 2002
! ===========================================================
SUBROUTINE CovarianceCoef(M, N, Data, Resemblance)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N
        REAL, DIMENSION(M,N), INTENT(IN) :: Data
        REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: Resemblance
        REAL, ALLOCATABLE, DIMENSION(:,:) :: X
        REAL, ALLOCATABLE, DIMENSION(:)   :: Mean
        INTEGER:: I, J, K, AllocatStatus
        REAL :: Unbias

        ALLOCATE( X(M,N), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "( X ) ** Not Enough Memory **"
	ALLOCATE( Mean(N) , STAT = AllocatStatus )
	IF (AllocatStatus/= 0 )  STOP  "( Mean ) ** Not Enough Memory **"
	
	! compute Mean of each col.
	DO J=1, N
		Mean(J) = SUM( Data(1:M,J) ) / REAL(M)
	END DO
	! subtract Mean of each col from every element in this col.
	DO J=1, N
		X(1:M,J) = Data(1:M,J) - Mean(J)	
	END DO
	Unbias = REAL(M - 1)
	Unbias = REAL(1/Unbias)
	! calculating the covariance matrix
        I = 1
        DO J=1, N-1
           DO K=J+1, N  
               Resemblance(I) = Unbias* DOT_PRODUCT(X(1:M,J), X(1:M,K))
               I = I + 1
           END DO
        END DO

        DEALLOCATE( X, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "( X )** Deallocation Error**"
        DEALLOCATE( Mean, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "( Mean )** Deallocation Error**"
END SUBROUTINE CovarianceCoef
! ===========================================================
! Subroutine Name: CosineCoef
! Author: Nusrat Yussouf
! Date: September, 2002
! ===========================================================
SUBROUTINE CosineCoef(M, N, Data, Resemblance)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N
        REAL, DIMENSION(M,N), INTENT(IN) :: Data
        REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: Resemblance
        REAL, ALLOCATABLE, DIMENSION(:)   :: Norm
        INTEGER:: I, J, K, AllocatStatus

        ALLOCATE( Norm(N), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "( Norm ) ** No Enough Memory **"     
	
	! calculating the norm
	DO J=1, N
	     Norm(J) = SQRT( DOT_PRODUCT( Data(1:M,J), Data(1:M,J) ) )
	END DO   
	I = 1
	DO J=1, N-1
        DO K=J+1, N
	Resemblance(I)=(DOT_PRODUCT(Data(1:M,J),Data(1:M,K)))/(Norm(J)*Norm(K))
        I = I + 1 
	END DO
	END DO	

        DEALLOCATE( Norm, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "( Norm )** Deallocation Error**"
END SUBROUTINE CosineCoef
! ===========================================================
! Subroutine Name: EuclideanSimiCoef
! Author: Ahmad Alhamad
! Date: 2000
! ===========================================================
SUBROUTINE EuclideanSimiCoef(M, N, Data, Resemblance)
	IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N
        REAL, DIMENSION(M,N), INTENT(IN) :: Data
        REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: Resemblance
        REAL, ALLOCATABLE, DIMENSION(:) :: Vector
        INTEGER:: I, J, K, AllocatStatus
        REAL :: Mx

        ALLOCATE( Vector(M), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** No Enough Memory **"
        I = 1
        DO J=1, N-1
                DO K=J+1, N
                        Vector = Data(1:M, J) - Data(1:M, K)
                        Resemblance(I) = SQRT(DOT_PRODUCT(Vector,Vector))
                        I = I + 1
                END DO
        END DO  
        Mx = MAXVAL(Resemblance)
        Resemblance = Resemblance / Mx
	DO I=1,N*(N-1)/2
		Resemblance(I) = 1.0 - Resemblance(I)
	END DO
        DEALLOCATE( Vector, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** Deallocation Error **"
END SUBROUTINE EuclideanSimiCoef
! ===========================================================
! Subroutine Name: Manhattan
! Author: Ahmad Alhamad
! Date: 2000 
! ===========================================================
SUBROUTINE Manhattan(M, N, Data, Resemblance)
	IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N
        REAL, DIMENSION(M,N), INTENT(IN) :: Data
        REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: Resemblance
        REAL, ALLOCATABLE, DIMENSION(:) :: Vector
        INTEGER:: I, J, K, AllocatStatus

        ALLOCATE( Vector(M), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** No Enough Memory **"
        I = 1
        DO J=1, N-1
                DO K=J+1, N
                        Vector = Data(1:M, J) - Data(1:M, K)
                        Resemblance(I) = SUM(ABS(Vector))
                        I = I + 1
                END DO
        END DO  
        DEALLOCATE( Vector, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** Deallocation Error **"
END SUBROUTINE Manhattan
! ===========================================================
! Subroutine Name: SupDistance
! Author: Ahmad Alhamad
! Date: 2000
! ===========================================================
SUBROUTINE SupDistance(M, N, Data, Resemblance)
	IMPLICIT NONE

        INTEGER, INTENT(IN) :: M, N
        REAL, DIMENSION(M,N), INTENT(IN) :: Data
        REAL, DIMENSION(N*(N-1)/2), INTENT(OUT) :: Resemblance
        REAL, ALLOCATABLE, DIMENSION(:) :: Vector
        INTEGER:: I, J, K, AllocatStatus

        ALLOCATE( Vector(M), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** No Enough Memory **"
        I = 1
        DO J=1, N-1
                DO K=J+1, N
                        Vector = Data(1:M, J) - Data(1:M, K)
                        Resemblance(I) = MAXVAL(ABS(Vector))
                        I = I + 1
                END DO
        END DO  
        DEALLOCATE( Vector, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Vector)** Deallocation Error **"
END SUBROUTINE SupDistance
! =======================================================================

