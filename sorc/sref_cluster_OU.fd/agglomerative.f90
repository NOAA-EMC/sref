! **********************************************************************
! File Name: agglomerative.f
! Author: Ahmad Alhamed
!         University of Oklahoma
!         School of Computer Science
! Date: 2000
! Modified by: Nusrat Yussouf
!              CIMMS/OU/NSSL
! Date: September, 2002
!
! Purpose: Execute clustering method
! Input:
!	N      : no. of initial clusters (col. in data matrix)
!       Method : used method, it can be one of the following
!          (1) : Single-link
!          (2) : Complete-link
!          (3) : UPGMA (group average)
!          (4) : WPGMA (weighted average)
!          (5) : UPGMC (unweighted centroid)
!          (6) : WPGMC (weighted centroid)
!          (7) : Ward's method (minimum variance)
!	Coeff  : resemblance coefficient
!	Dissim : Dissimilarity coefficient or similarity coefficient
!	Res    : Resemblance matrix
! Ooutput:
! 	Hierarchy : result of clustering to be used to draw tree
! **********************************************************************
!
SUBROUTINE Agglomerative(N, Method, Coeff, Dissim, Res, Hierarchy,cluster)
	USE Global
	IMPLICIT NONE

        INTEGER, INTENT(IN) :: N, Method, Coeff
	LOGICAL, INTENT(IN) :: Dissim
        REAL, DIMENSION(N*(N-1)/2), INTENT(IN) :: Res
	TYPE(Merge), DIMENSION(N-1), INTENT(OUT) :: Hierarchy

! Local Variables
	REAL :: LargeDissimilarity, SmallSimilarity
	INTEGER :: Offset    ! Func. returns index in Resemblance Matrix
	REAL    :: UpdateResDeg ! Func. returns new degree of resemblance
	INTEGER :: Index, I, J, K, L, AllocatStatus, Len, iclus
	REAL    :: Degree
	REAL, ALLOCATABLE, DIMENSION(:) :: Temp, Resemblance
	INTEGER, ALLOCATABLE, DIMENSION(:) :: Members
	LOGICAL, ALLOCATABLE, DIMENSION(:) :: Exist
          integer ncluster
            integer nm,ntime,nfile
            parameter(nm=26,nfile=nm+1,ntime=30)
            parameter (ncluster = 6)

         integer, dimension(ncluster,nm), INTENT(out) :: cluster 
         integer index1,ntemp
         integer icount
         integer nidx(ncluster)
         integer yy3(ntime),mm3(ntime),dd3(ntime),hh3(ntime)
         integer kftime(nfile,ntime)
        character*6 mem (nm)
        common /c2/yy3,mm3,dd3,hh3,kftime,mem  


	Len = N*(N-1) / 2
	LargeDissimilarity = 32.0E+23
	SmallSimilarity    = -32.0E+23
	! Allocate memory space for arrays
        ALLOCATE( Members(N), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Members)** No Enough Memory **"
        ALLOCATE( Exist(N), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Exist)** No Enough Memory **"
        ALLOCATE( Resemblance(N*(N-1)/2), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Temp)** No Enough Memory **"
        ALLOCATE( Temp(N*(N-1)/2), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Resemblance)** No Enough Memory **"
	Resemblance = Res
	! Initialize Members and Exist
	Members = 1
	Exist = .TRUE.
	! Execute agglomerative clustering
	DO I = 1, N-1
		Temp = Resemblance
		! Find the most similar pair
		IF (Dissim)  THEN
			Degree = MINVAL(Resemblance)
			DO L=1, Len
				IF (Resemblance(L) == Degree) THEN
					Index = L
					EXIT
				END IF
			END DO			
		ELSE
                        Degree = MAXVAL(Resemblance)
			DO L=1, Len
                                IF (Resemblance(L) == Degree) THEN
                                        Index = L
                                        EXIT
                                END IF
			END DO
		END IF
		! Find corresponding cluster for this Index
		CALL ClustersOfIndex(N, Index, J, K)
		! Update Hierarchy
		Hierarchy(I)%MergedCluster1 = J
		Hierarchy(I)%MergedCluster2 = K
		Hierarchy(I)%NewLable = J
		IF (Method /= 7 ) THEN  ! for methods other than Ward's
		    Hierarchy(I)%ResDegree = Degree
		ELSE    ! for Ward's
		    IF (I>1) THEN
			IF (Coeff == 2) THEN  ! for square Euc
			   Hierarchy(I)%ResDegree = SQRT(Hierarchy(I-1)%ResDegree**2.0 &
						+ Degree/2.0)
			ELSE
			   Hierarchy(I)%ResDegree = Hierarchy(I-1)%ResDegree &
						+ Degree
			END IF
		    ELSE
			IF (Coeff == 2) THEN  ! for square Euc
			   Hierarchy(I)%ResDegree = SQRT(Degree/2.0)
			ELSE
			   Hierarchy(I)%ResDegree = Degree
			END IF
		    END IF
		END IF
		! Set the corresponding index of K in Exist array to False
		Exist(K) = .FALSE.
		! Update Resemblance matrix
		DO L = 1, N
			IF ( L/=J .AND. L/=K .AND. Exist(L) )  THEN
				Index = Offset(N, L, J)
				Resemblance(Index) = UpdateResDeg(N, L, J, K, &
					       &   Method, Members, Temp)
				Index = Offset(N, L, K)
				IF (Dissim)  THEN
					Resemblance(Index) = LargeDissimilarity
				ELSE
					Resemblance(Index) = SmallSimilarity
				END IF
			END IF 
		END DO
		! Cancel the value of the larger index of the newly merged
		! cluster in the Resemblance matrix
		Index = Offset(N, J, K)
		IF (Dissim)  THEN
			Resemblance(Index) = LargeDissimilarity
		ELSE
			Resemblance(Index) = SmallSimilarity
		END IF
		! Change no. of members with respect to the newly 
		! merged clusters
		Members(J) = Members(J) + Members(K)
		Members(K) = 0
	        WRITE (7, *) '---------------------------------------------------&
	             &----------------'
		WRITE (7, *) 'Resemblance at step:', I
		CALL PrintResemblanceMatrix(N, Resemblance)
		WRITE (7, 10) Hierarchy(I)%MergedCluster1, &
		&             Hierarchy(I)%MergedCluster2, &
		& Hierarchy(I)%NewLable, Hierarchy(I)%ResDegree
10 FORMAT ('Merge ',I2,' and ',I2,' to form ',I2,' at ==> ',F10.2)
	END DO

! six clusters  members of each cluster
          do I =1, ncluster
          do J =1, N
          cluster (I,j) = 0
          end do
          end do
          cluster (1,1) = Hierarchy(N-1)%MergedCluster1
          cluster (2,1) = Hierarchy(N-1)%MergedCluster2
          do I=2,ncluster-1
             index1 = 0
             do J=1,I   
          if (Hierarchy(N-I)%MergedCluster1 == cluster(J,1)) then
             index1 = 1
             end if
             enddo
          if (index1 == 1) then
             cluster (I+1,1) = Hierarchy(N-I)%MergedCluster2
          else 
             cluster (I+1,1) = Hierarchy(N-I)%MergedCluster1
          end if
          enddo
            do J = 1, 6
             K = 1
             icount = 1
1211         continue 
             do I=N-ncluster,1, -1
             if (Hierarchy(I)%NewLable == cluster(j,k) ) then
                icount = icount + 1
                if(Hierarchy(I)%Mergedcluster1 == cluster(j,k)) then
                cluster(j,icount) = Hierarchy(I)%MergedCluster2
                else
                cluster(j,icount) = Hierarchy(I)%MergedCluster1
                end if
              endif
              enddo
              if (icount > K) then            
              K=K+1
               go to 1211
              endif

            enddo
 
               do I=1,ncluster
                    do k=1,n-1
                     do j=k+1,n
                      if (cluster(I,k) > cluster(i,j)) then
                        ntemp = cluster (i,k)
                         cluster(I,K)=cluster(I,j)
                         cluster(I,J)=ntemp
                      endif
                     enddo
                    enddo
               enddo
                   do i =1, ncluster
                     nidx(i) = 0
                   enddo
 
              do I=1,ncluster
                       j=1
1202                  continue
                      if(cluster(i,j) == 0 ) then
                        j=j+1
                         goto 1202
                      else
                        nidx(i)=cluster(i,j)
                      end if
              enddo

                 do i = 1, ncluster-1
                 do j = i,ncluster
                      if(nidx(i) > nidx(j)) then
                         ntemp = nidx(i)
                         nidx(i) = nidx(j)
                         nidx(j)= ntemp
                      
                        do k =1, n
                         ntemp = cluster(i,k)
                         cluster(i,k) = cluster(j,k)
                         cluster(j,k) = ntemp
                       enddo
                      endif
                 enddo
                 enddo


      
                        

             write(7,*) '*********************************'
!            write(777,*) '*********************************'
            do I=1,ncluster 

             write(7,*) '*****************cluster ',I
            iclus=0
              do K=1,N
             if (cluster(I,K) /= 0) then
            iclus=iclus+1
            write(7,*) cluster(I,k),'----',mem(cluster(i,k))
             end if 
              enddo

!            write(777,*) ' Cluster #',I,'has',iclus,'forecasts containing the following members:'
             write(777,779) ' Cluster # ',I,' has ',iclus,' forecasts containing the following members:'
779          format(1a,i1,1a,i2,1a)
              do K=1,N
             if (cluster(I,K) /= 0) then
            write(777,*) '  ',mem(cluster(i,k)),' (',cluster(I,k),')'
             end if 
             enddo

            enddo
              
	! Deallocate arrays
        DEALLOCATE( Members, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Members)** Deallocation Error **"
        DEALLOCATE( Exist, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Exist)** Deallocation Error **"
        DEALLOCATE( Resemblance, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Resemblance)** Deallocation Error **"
        DEALLOCATE( Temp, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Temp)** Deallocation Error **"
END SUBROUTINE Agglomerative
! ======================================================================
! **********************************************************************
! ClustersOfIndex:
! Finds the corresponding clusters of an index in the Resemblance matrix
! INPUTS:
!	N  : No. of clusters in the Resemblance matrix
!	I  : Index in the Resemblance matrix
! OUTPUTS
!	J, K : Clusters correspond to I
! **********************************************************************
SUBROUTINE ClustersOfIndex(N, I, J, K)
	IMPLICIT NONE

	INTEGER, INTENT(IN)  :: N, I
	INTEGER, INTENT(OUT) :: J, K
	INTEGER :: Prev, Loc, R

	Prev = 0
	Loc = 0
	DO R =1, N-1
		Prev = Loc
		Loc = Loc + (N-R)
		IF (I <= Loc)  THEN
			J = R
			EXIT
		END IF
	END DO
	K = J + (I-Prev)
END SUBROUTINE ClustersOfIndex
! ====================================================================
! ********************************************************************
! UpdateResDeg: 
! Returns new resemblnace using Lance and Williams' formula which 
! given in :
!   A. Jain, R. Dubes, "Algorithms for Clustering Data", Prentice Hall, 
!   1988, PP. 79-80
! INPUTS:
!	N : Original no. of clusters = no. of col in Data matrix
!	K : Existing cluster
!	R, S : newly formed cluster
!	Method : used method, it can be one of the following
!	   (1) : Single-link
!	   (2) : Complete-link
!	   (3) : UPGMA (group average)
!	   (4) : WPGMA (weighted average)
!	   (5) : UPGMC (unweighted centroid)
!	   (6) : WPGMC (weighted centroid)
!	   (7) : Ward's method (minimum variance)
!	Members: no. of members in each cluster
!	Res    : Resemblance matrix
! OUTPUT:
!	returns the new resemblance degree
! ********************************************************************
REAL FUNCTION UpdateResDeg(N, K, R, S, Method, Members, Res)
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: N, K, R, S, Method
	INTEGER, DIMENSION(N), INTENT(IN) :: Members
	REAL, DIMENSION(N*(N-1)/2), INTENT(IN) :: Res
	REAL :: AlphaR, AlphaS, Beta, Gamma
        INTEGER :: Offset    ! Func. returns index in Resemblance Matrix

        SELECT CASE(Method)
        CASE (1)
		AlphaR = 0.50
		AlphaS = 0.50
		Beta   = 0.0
		Gamma  = -0.50
        CASE (2)
                AlphaR = 0.50
                AlphaS = 0.50
                Beta   = 0.0
                Gamma  = 0.50
	CASE (3)
                AlphaR = REAL(Members(R)) / ( REAL(Members(R)) &
			& + REAL(Members(S)) )
                AlphaS = REAL(Members(S)) / ( REAL(Members(R)) &
			& + REAL(Members(S)) )
                Beta   = 0.0
                Gamma  = 0.0	
	CASE (4)
                AlphaR = 0.50
                AlphaS = 0.50
                Beta   = 0.0
                Gamma  = 0.0
	CASE (5)
                AlphaR = REAL(Members(R)) / ( REAL(Members(R)) &
			& + REAL(Members(S)) )
                AlphaS = REAL(Members(S)) / ( REAL(Members(R)) &
			& + REAL(Members(S)) )
                Beta   = -1.0*REAL(Members(R))*REAL(Members(S)) &
			& / ( REAL(Members(R)) + REAL(Members(S)) )**2.0
                Gamma  = 0.0
	CASE (6)
                AlphaR = 0.50
                AlphaS = 0.50
                Beta   = -0.250
                Gamma  = 0.0
	CASE (7)
                AlphaR =( REAL(Members(R))+REAL(Members(K)) ) &
			& / ( REAL(Members(R))+REAL(Members(S)) &
			&     +REAL(Members(K)) )
                AlphaS =( REAL(Members(S))+REAL(Members(K)) ) &
			& / ( REAL(Members(R))+REAL(Members(S)) &
			&     +REAL(Members(K)) )
                Beta   = ( -1.0*REAL(Members(K)) )/( REAL(Members(R)) & 
			& +REAL(Members(S))+REAL(Members(K)) )
                Gamma  = 0.0
        END SELECT
	UpdateResDeg = AlphaR*Res(Offset(N,K,R)) + AlphaS*Res(Offset(N,K,S)) &
		  & + Beta*Res(Offset(N,R,S)) + Gamma*ABS(Res(Offset(N,K,R)) &
		  & - Res(Offset(N,K,S)) )
END FUNCTION UpdateResDeg
! ==================================================================
! ******************************************************************
! Return index of clusters J, and K  in the Resemblance matrix
! INPUTS:
!       N : No. of clusters in the Resemblance matrix
!       J, K : Clusters to find their index in the Resemblance matrix
! OUTPUT:
!       Returns index of J, and K in the Resemblance matrix
! *******************************************************************
INTEGER FUNCTION Offset(N, J, K)
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N, J, K
        INTEGER :: P, Q

        P = MIN(J, K)
        Q = MAX(J, K)
        Offset = Q + (P-1) * N - (P*(P+1)) / 2
END FUNCTION Offset
! ==================================================================
