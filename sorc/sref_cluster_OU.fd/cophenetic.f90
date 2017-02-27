! ***********************************************************************
! File Name: cophenetic.f90
! Author: Ahmad Alhamed 
!         University of Oklahoma
!         School of Computer Science
! Date: 2000
! Modified by: Nusrat Yussouf
!              CIMMS/OU/NSSL
! Date: September, 2002
! Purpose:  Compute the cophenetic correlation coefficient
! Ref: Romesburg, H., "Cluster Analysis for Researchers",
!      Lifetime Learning Pub., 1984, PP 24-27
! Input:
!	N         : no. of clusters (col. in data matrix)
!	Hierarchy : hierarch of clusters as produced by agglomerative
!       Res    : Resemblance matrix
!    	Method : string containing the method names
! Output:
!       CopheneticCoeff: Cophenetic Correlation Coeff.
! ***********************************************************************
!
SUBROUTINE Cophenetic(N, Hierarchy, Res, CopheneticCoeff)
        USE Global
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        TYPE(Merge), DIMENSION(N-1), INTENT(IN) :: Hierarchy
	REAL, DIMENSION(N*(N-1)/2), INTENT(IN) :: Res
	REAL, INTENT(OUT) :: CopheneticCoeff
	INTEGER :: Offset    ! Func. returns index in CopheneticMatrix
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Clusters
        INTEGER, ALLOCATABLE, DIMENSION(:) :: Members      
	REAL, ALLOCATABLE, DIMENSION(:) :: CopheneticMatrix
	LOGICAL, ALLOCATABLE, DIMENSION(:) :: Exist
	REAL :: SUMxy, SUMx2, SUMy2, SUMx, SUMy
        INTEGER :: AllocatStatus, I, J, K, P, Q, NDX, L

	L = N*(N-1)/2
        ALLOCATE( Clusters(N,N), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "(Clusters) ** Not Enough Memory **"
        ALLOCATE( Members(N), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "(Members) ** Not Enough Memory **"
        ALLOCATE( CopheneticMatrix(L), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "(CopheneticMatrix) ** Not Enough Memory **"
        ALLOCATE( Exist(L), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Exist)** Not Enough Memory **"
        DO I=1, N
                Clusters(I, 1) = I
                Members(I) = 1
        END DO
	Exist = .FALSE.
        DO I=1, N-1
                J = Hierarchy(I)%MergedCluster1
                K = Hierarchy(I)%MergedCluster2
		! STEP 1: merge clusters
                DO Q=1,    Members(K)
                        Members(J) = Members(J) + 1
                        NDX = Members(J)
                        Clusters(J, NDX) = Clusters(K, Q)
                END DO
		! STEP 2: compute cophenetic matrix
                DO P=1,  Members(J)
		   DO Q=1, Members(K) 
			IF ( Clusters(J,P) /= Clusters(K,Q) ) THEN
                            NDX = Offset(N, Clusters(J,P), Clusters(K,Q) )
			    IF ( EXIST(NDX) .EQV. .FALSE. )  THEN
			       CopheneticMatrix(NDX) = Hierarchy(I)%ResDegree
			       EXIST(NDX) = .TRUE.
			    END IF
			END IF
		   END DO
                END DO
        END DO
	! compute Cophenetic correlation Coeff
	SUMxy = 0.0
	SUMx2 = 0.0
	SUMy2 = 0.0
	SUMx  = 0.0
	SUMy  = 0.0
	DO I=1, L
		SUMxy = SUMxy + ( Res(I) * CopheneticMatrix(I) )
		SUMx2 = SUMx2 + Res(I)**2.0
		SUMy2 = SUMy2 + CopheneticMatrix(I)**2.0
		SUMx  = SUMx + Res(I)
		SUMy  = SUMy + CopheneticMatrix(I)
	END DO
	CopheneticCoeff = ( SUMxy - (1/REAL(L))*SUMx*SUMy ) / &
			& ( SQRT( (SUMx2-(1/REAL(L))*SUMx**2.0 ) * &
			& (SUMy2-(1/REAL(L))*SUMy**2.0 ) ))
END SUBROUTINE Cophenetic
!=========================================================================
