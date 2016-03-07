! *****************************************************************
! File Name: tree.f
! Author: Ahmad Alhamed
!         University of Oklahoma
!         School of Computer Science
! Date: 2000
! Modified by: Nusrat Yussouf
!              CIMMS/OU/NSSL
! Date: September, 2002
!
! Purpose: Draw clustering tree or the Dendrogram
! *****************************************************************
SUBROUTINE Tree(N, Hierarchy, Title, Dissim,&
		& Scale, OriginX, OriginY, CophCorrCoef)
      USE Global
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: Dissim
      TYPE(Merge), DIMENSION(N-1), INTENT(IN) :: Hierarchy
      CHARACTER(LEN=*), INTENT(IN) :: Title
      REAL, INTENT(IN) :: Scale, OriginX, OriginY, CophCorrCoef
      
      IF (N >= 28) THEN
         CALL TreeLabelInTwoLines(N, Hierarchy, Title, Dissim,&
		& Scale, OriginX, OriginY, CophCorrCoef)
      ELSE IF( N <= 27) THEN
	 CALL TreeModelName(N, Hierarchy, Title, Dissim,&
		& Scale, OriginX, OriginY, CophCorrCoef)      
!         CALL TreeModelNumber(N, Hierarchy, Title, Dissim,&
!		& Scale, OriginX, OriginY, CophCorrCoef)
      END IF				     
END SUBROUTINE Tree   
!====================================================================
! Subroutine name: TreeModelName
! Author: Ahmad Alhamed
! Modified by: Nusrat Yussouf
! Date: 2002
! ===================================================================
SUBROUTINE TreeModelName(N, Hierarchy, Title, Dissim,&
		& Scale, OriginX, OriginY, CophCorrCoef)
      USE Global
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: Dissim
      TYPE(Merge), DIMENSION(N-1), INTENT(IN) :: Hierarchy
      CHARACTER(LEN=*), INTENT(IN) :: Title
      REAL, INTENT(IN) :: Scale, OriginX, OriginY, CophCorrCoef
      INTEGER :: appid, wid, gkswid, ierr
      REAL, DIMENSION(N) :: XLoc, YLoc
      REAL :: HUnit, VUnit, Height, Width, Offset, X, Y, Dot
      REAL :: Range1, Range2, FontScale
      INTEGER :: I, J, K
      REAL :: MaxDegree, DegreeScale, prevY, prevDegree
      INTEGER, DIMENSION(N) :: OrderedList
      CHARACTER*3, DIMENSION(N) :: XLabel1, XLabel2, & 
           & XLabel3, XLabel4
      CHARACTER(LEN=15) :: Degree
      CHARACTER(LEN=80) :: Buff

! order cluster on X axis
      CALL OrderClusters(N, Hierarchy, OrderedList)

! INITIALIZE
      IF (OutputFormat==ATOMOP .OR. OutputFormat==OTAMOP) THEN 
        IF (TPP == 4) THEN
	    Height = 0.45
	    Width = 0.48
	    Offset = 0.03
	    Dot = 0.0016
	    FontScale = 1.0
	END IF
        IF (TPP == 9) THEN
	    Height = 0.30
	    Width = 0.32
	    Offset = 0.03
	    Dot = 0.0016
	    FontScale = 1.0
	END IF
      ELSE
	! OutputFormat==OTOMOP .OR. OutputFormat==OTAMSP 
	Height = 1.0
	Width = 1.0
	Offset = 0.1
	Dot = 0.005
	FontScale = 2.0
     END IF

      IF (Dissim)  THEN
	SELECT CASE(OutputFormat)
	CASE(OTOMOP, OTAMSP, OTAMOP)
        	MaxDegree = Hierarchy(N-1)%ResDegree
	CASE(ATOMOP)
		MaxDegree = Scale
	END SELECT
      ELSE
         MaxDegree = 1.0
      END IF
      prevDegree = MaxDegree
      HUnit = (Width-(3.0*Offset)) / REAL(N)
      IF (HUnit == 0.0)   HUnit = Dot
      VUnit = (Height-(3.5*Offset)) / MaxDegree
	
! initialize X locations for N clusters
      X = OriginX + (2.0*Offset)
      DO I=1, N
	X = X + HUnit
      	XLoc(OrderedList(I)) = X
      END DO
      YLoc = OriginY + (2.0*Offset)    ! for all clusters

      ! display title
!      CALL PLCHLQ(OriginX+(Width/2.0), OriginY+(Height-0.3*Offset),&
!      & TRIM(Title), 9.0*FontScale,0.,0.)
      ! Set the line width
!      CALL GSLWSC(2.)
      ! draw scale to the left of tree
!      CALL LINE(OriginX+(2.0*Offset), OriginY+2.0*Offset, &
!      &  OriginX+(2.0*Offset), OriginY+(Height-0.7*Offset) )
!      CALL LINE(OriginX+(2.0*Offset-(3.0*Dot)),OriginY+2.0*Offset, &
!      &  OriginX+(2.0*Offset+(3.0*Dot)),OriginY+2.0*Offset)
!      CALL LINE(OriginX+(2.0*Offset-(3.0*Dot)),OriginY+(Height-0.7*Offset),&
!      &  OriginX+(2.0*Offset+(3.0*Dot)), OriginY+(Height-0.7*Offset))
      ! lable X axis
!      CALL LabelModelRun( N, OrderedList, XLabel1, XLabel2, XLabel3, XLabel4)
!      Do I = 1, N
!          XLabel1(I) = TRIM(XLabel1(I))
!      END do
!      DO I=1, N
!	CALL PLCHLQ(XLoc(OrderedList(I))+(2.0*Dot), OriginY+(1.7*Offset),&
!	& XLabel1(I) ,5*FontScale,0.,1.)
!	CALL PLCHLQ(XLoc(OrderedList(I))+(2.0*Dot), OriginY+(1.55*Offset),&
!	& XLabel2(I) ,5*FontScale,0.,1.)
!	CALL PLCHLQ(XLoc(OrderedList(I))+(2.0*Dot), OriginY+(1.4*Offset),&
!	& XLabel3(I) ,5*FontScale,0.,1.)
!	CALL PLCHLQ(XLoc(OrderedList(I))+(2.0*Dot), OriginY+(1.25*Offset),&
!	& XLabel4(I) ,5*FontScale,0.,1.)
!      END DO
      WRITE(Buff, 31)  Hierarchy(1)%ResDegree, Hierarchy(N-1)%ResDegree
31    FORMAT ('Range: ',F15.3,' , ',F15.3)
!      CALL PLCHLQ(OriginX+1.0*Offset, OriginY+(0.75*Offset),&
!	&  TRIM(Buff),9.0*FontScale,0.,-1.)
      WRITE(Buff, 32) CophCorrCoef
32    FORMAT ('Cophenetic Corr. Coeff. = ',F4.2)
!      CALL PLCHLQ(OriginX+1.0*Offset, OriginY+(.1*Offset),&
!	&  TRIM(Buff),9.*FontScale,0.,-1.)
      ! set the line width to smaller
!      CALL GSLWSC(1.)
      ! draw the tree
      prevY = OriginY
      Y = OriginY + (2.0*Offset)
      DO I=1, N-1
	! in case of tie, do not update Y
        IF(Hierarchy(I)%ResDegree /= prevDegree)  THEN
	   IF (Dissim)  THEN
	   	Y = OriginY + (2.0*Offset) + Hierarchy(I)%ResDegree*VUnit
           ELSE
                Y = OriginY+(2.0*Offset)+(1.0 - Hierarchy(I)%ResDegree)*VUnit
	   END IF
        END IF
	J = Hierarchy(I)%MergedCluster1
	K = Hierarchy(I)%MergedCluster2
	! NCARG code to draw line
!	CALL LINE(XLoc(J), YLoc(J), XLoc(J), Y)
!	CALL LINE(XLoc(K), YLoc(K), XLoc(K), Y)
!	CALL LINE(XLoc(J), Y, XLoc(K), Y)
        IF(Hierarchy(I)%ResDegree /= prevDegree)  THEN
		! draw tick on the scale, and display degree
!		CALL LINE(OriginX+(2.0*Offset-(3.0*Dot)), Y, &
!		&         OriginX+(2.0*Offset+(3.0*Dot)), Y)
		! if degree close to each other, they will be
		! overlaped on the plot. Skip if overlap
		IF (Dissim)  THEN
			WRITE(Degree, '(F15.3)' ) Hierarchy(I)%ResDegree
		ELSE
			WRITE(Degree, '(F15.5)' ) Hierarchy(I)%ResDegree
		END IF
                IF ( Y-prevY > (5.0*Dot) ) THEN
!                        CALL PLCHLQ(OriginX+(2.0*Offset-(6.0*Dot)), &
!                        & Y, Degree,6.0*FontScale,0.,1.)
                        prevY = Y
                END IF
	END IF
	! update X, Y locations
	XLoc(J) = XLoc(J) + (XLoc(K)-XLoc(J))/2.0
	YLoc(J) = Y
	prevDegree = Hierarchy(I)%ResDegree
      END DO
      ! draw line from last merge to the top
!      CALL LINE(XLoc(J), Y, XLoc(J), Y+0.4*Offset )	
END SUBROUTINE TreeModelName
! ===========================================================================
! File Name: tree.f
! Author: Ahmad Alhamed
! Modified by: Nusrat Yussouf
! Draw clustering tree
! *****************************************************************
SUBROUTINE TreeModelNumber(N, Hierarchy, Title, Dissim,&
		& Scale, OriginX, OriginY, CophCorrCoef)
      USE Global
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: Dissim
      TYPE(Merge), DIMENSION(N-1), INTENT(IN) :: Hierarchy
      CHARACTER(LEN=*), INTENT(IN) :: Title
      REAL, INTENT(IN) :: Scale, OriginX, OriginY, CophCorrCoef
      INTEGER :: appid, wid, gkswid, ierr
      REAL, DIMENSION(N) :: XLoc, YLoc
      REAL :: HUnit, VUnit, Height, Width, Offset, X, Y, Dot
      REAL :: Range1, Range2, FontScale
      INTEGER :: I, J, K
      REAL :: MaxDegree, DegreeScale, prevY, prevDegree
      INTEGER, DIMENSION(N) :: OrderedList
      CHARACTER(LEN=3) :: Label
      CHARACTER(LEN=15) :: Degree
      CHARACTER(LEN=80) :: Buff

! order cluster on X axis
      CALL OrderClusters(N, Hierarchy, OrderedList)

! INITIALIZE
      IF (OutputFormat==ATOMOP .OR. OutputFormat==OTAMOP) THEN 
        IF (TPP == 4) THEN
	    Height = 0.45
	    Width = 0.48
	    Offset = 0.03
	    Dot = 0.0016
	    FontScale = 1.0
	END IF
        IF (TPP == 9) THEN
	    Height = 0.30
	    Width = 0.32
	    Offset = 0.03
	    Dot = 0.0016
	    FontScale = 1.0
	END IF
      ELSE
	! OutputFormat==OTOMOP .OR. OutputFormat==OTAMSP 
	Height = 1.0
	Width = 1.0
	Offset = 0.1
	Dot = 0.005
	FontScale = 2.0
     END IF

      IF (Dissim)  THEN
	SELECT CASE(OutputFormat)
	CASE(OTOMOP, OTAMSP, OTAMOP)
        	MaxDegree = Hierarchy(N-1)%ResDegree
	CASE(ATOMOP)
		MaxDegree = Scale
	END SELECT
      ELSE
         MaxDegree = 1.0
      END IF
      prevDegree = MaxDegree
      HUnit = (Width-(3.0*Offset)) / REAL(N)
      IF (HUnit == 0.0)   HUnit = Dot
      VUnit = (Height-(3.5*Offset)) / MaxDegree
	
! initialize X locations for N clusters
      X = OriginX + (2.0*Offset)
      DO I=1, N
	X = X + HUnit
      	XLoc(OrderedList(I)) = X
      END DO
      YLoc = OriginY + (2.0*Offset)    ! for all clusters

      ! display title
!      CALL PLCHLQ(OriginX+(Width/2.0), OriginY+(Height-0.3*Offset),&
!      & TRIM(Title), 9.0*FontScale,0.,0.)
      ! Set the line width
!      CALL GSLWSC(2.)
      ! draw scale to the left of tree
!      CALL LINE(OriginX+(2.0*Offset), OriginY+2.0*Offset, &
!      &  OriginX+(2.0*Offset), OriginY+(Height-0.7*Offset) )
!      CALL LINE(OriginX+(2.0*Offset-(3.0*Dot)),OriginY+2.0*Offset, &
!      &  OriginX+(2.0*Offset+(3.0*Dot)),OriginY+2.0*Offset)
!      CALL LINE(OriginX+(2.0*Offset-(3.0*Dot)),OriginY+(Height-0.7*Offset),&
!      &  OriginX+(2.0*Offset+(3.0*Dot)), OriginY+(Height-0.7*Offset))
      ! lable X axis
      DO I=1, N
	WRITE(Label, '(A3)' ) TRIM(LABLES(OrderedList(I)))
!	CALL PLCHLQ(XLoc(OrderedList(I))+(2.0*Dot), OriginY+(1.7*Offset),&
!	& Label ,5*FontScale,0.,1.)
      END DO
      WRITE(Buff, 31)  Hierarchy(1)%ResDegree, Hierarchy(N-1)%ResDegree
31    FORMAT ('Range: ',F15.3,' , ',F15.3)
!      CALL PLCHLQ(OriginX+Offset, OriginY+(1.0*Offset),&
!	&  TRIM(Buff),9.0*FontScale,0.,-1.)
      WRITE(Buff, 32) CophCorrCoef
32    FORMAT ('Cophenetic Corr. Coeff. = ',F4.2)
!      CALL PLCHLQ(OriginX+Offset, OriginY+(.5*Offset),&
!	&  TRIM(Buff),9.*FontScale,0.,-1.)
      ! set the line width to smaller
!      CALL GSLWSC(1.)
      ! draw the tree
      prevY = OriginY
      Y = OriginY + (2.0*Offset)
      DO I=1, N-1
	! in case of tie, do not update Y
        IF(Hierarchy(I)%ResDegree /= prevDegree)  THEN
	   IF (Dissim)  THEN
	   	Y = OriginY + (2.0*Offset) + Hierarchy(I)%ResDegree*VUnit
           ELSE
                Y = OriginY+(2.0*Offset)+(1.0 - Hierarchy(I)%ResDegree)*VUnit
	   END IF
        END IF
	J = Hierarchy(I)%MergedCluster1
	K = Hierarchy(I)%MergedCluster2
	! NCARG code to draw line
!	CALL LINE(XLoc(J), YLoc(J), XLoc(J), Y)
!	CALL LINE(XLoc(K), YLoc(K), XLoc(K), Y)
!	CALL LINE(XLoc(J), Y, XLoc(K), Y)
        IF(Hierarchy(I)%ResDegree /= prevDegree)  THEN
		! draw tick on the scale, and display degree
!		CALL LINE(OriginX+(2.0*Offset-(3.0*Dot)), Y, &
!		&         OriginX+(2.0*Offset+(3.0*Dot)), Y)
		! if degree close to each other, they will be
		! overlaped on the plot. Skip if overlap
		IF (Dissim)  THEN
			WRITE(Degree, '(F15.3)' ) Hierarchy(I)%ResDegree
		ELSE
			WRITE(Degree, '(F15.5)' ) Hierarchy(I)%ResDegree
		END IF
                IF ( Y-prevY > (5.0*Dot) ) THEN
!                        CALL PLCHLQ(OriginX+(2.0*Offset-(6.0*Dot)), &
!                        & Y, Degree,6.0*FontScale,0.,1.)
                        prevY = Y
                END IF
	END IF
	! update X, Y locations
	XLoc(J) = XLoc(J) + (XLoc(K)-XLoc(J))/2.0
	YLoc(J) = Y
	prevDegree = Hierarchy(I)%ResDegree
      END DO
      ! draw line from last merge to the top
!      CALL LINE(XLoc(J), Y, XLoc(J), Y+0.4*Offset )	
END SUBROUTINE TreeModelNumber
! *****************************************************************
! File Name: tree.f
! Author: Ahmad Alhamed
! Modified by: Nusrat Yussouf
! Draw clustering tree
! *****************************************************************
SUBROUTINE TreeLabelInTwoLines(N, Hierarchy, Title, Dissim,&
		& Scale, OriginX, OriginY, CophCorrCoef)
      USE Global
      IMPLICIT NONE

      INTEGER, INTENT(IN) :: N
      LOGICAL, INTENT(IN) :: Dissim
      TYPE(Merge), DIMENSION(N-1), INTENT(IN) :: Hierarchy
      CHARACTER(LEN=*), INTENT(IN) :: Title
      REAL, INTENT(IN) :: Scale, OriginX, OriginY, CophCorrCoef
      INTEGER :: appid, wid, gkswid, ierr
      REAL, DIMENSION(N) :: XLoc, YLoc
      REAL :: HUnit, VUnit, Height, Width, Offset, X, Y, Dot
      REAL :: Range1, Range2, FontScale
      INTEGER :: I, J, K
      REAL :: MaxDegree, DegreeScale, prevY, prevDegree
      INTEGER, DIMENSION(N) :: OrderedList
      CHARACTER(LEN=3) :: Label
      CHARACTER(LEN=15) :: Degree
      CHARACTER(LEN=80) :: Buff

! order cluster on X axis
      CALL OrderClusters(N, Hierarchy, OrderedList)

! INITIALIZE
      IF (OutputFormat==ATOMOP .OR. OutputFormat==OTAMOP) THEN 
        IF (TPP == 4) THEN
	    Height = 0.45
	    Width = 0.48
	    Offset = 0.03
	    Dot = 0.0016
	    FontScale = 1.0
	END IF
        IF (TPP == 9) THEN
	    Height = 0.30
	    Width = 0.32
	    Offset = 0.03
	    Dot = 0.0016
	    FontScale = 1.0
	END IF
      ELSE
	! OutputFormat==OTOMOP .OR. OutputFormat==OTAMSP 
	Height = 1.0
	Width = 1.0
	Offset = 0.1
	Dot = 0.005
	FontScale = 2.5
     END IF

      IF (Dissim)  THEN
	SELECT CASE(OutputFormat)
	CASE(OTOMOP, OTAMSP, OTAMOP)
        	MaxDegree = Hierarchy(N-1)%ResDegree
	CASE(ATOMOP)
		MaxDegree = Scale
	END SELECT
      ELSE
         MaxDegree = 1.0
      END IF
      prevDegree = MaxDegree
      HUnit = (Width-(3.0*Offset)) / REAL(N)
      IF (HUnit == 0.0)   HUnit = Dot
      VUnit = (Height-(3.5*Offset)) / MaxDegree
	
! initialize X locations for N clusters
      X = OriginX + (2.0*Offset)
      DO I=1, N
	X = X + HUnit
      	XLoc(OrderedList(I)) = X
      END DO
      YLoc = OriginY + (2.0*Offset)    ! for all clusters

      ! display title
!      CALL PLCHLQ(OriginX+(Width/2.0), OriginY+(Height-0.3*Offset),&
!      & TRIM(Title), 9.0*FontScale,0.,0.)
      ! Set the line width
!      CALL GSLWSC(2.)
      ! draw scale to the left of tree
!      CALL LINE(OriginX+(2.0*Offset), OriginY+2.0*Offset, &
!      &  OriginX+(2.0*Offset), OriginY+(Height-0.7*Offset) )
!      CALL LINE(OriginX+(2.0*Offset-(3.0*Dot)),OriginY+2.0*Offset, &
!      &  OriginX+(2.0*Offset+(3.0*Dot)),OriginY+2.0*Offset)
!      CALL LINE(OriginX+(2.0*Offset-(3.0*Dot)),OriginY+(Height-0.7*Offset),&
!      &  OriginX+(2.0*Offset+(3.0*Dot)), OriginY+(Height-0.7*Offset))
      
      ! lable X axis in two line, so that it is readable
      DO I=1, N, 2
	WRITE(Label, '(A3)' ) TRIM(LABLES(OrderedList(I)))
!	CALL PLCHLQ(XLoc(OrderedList(I))+(2.0*Dot), OriginY+(1.7*Offset),&
!	& Label ,5*FontScale,0.,1.)
      END DO
      DO I=2, N, 2
	WRITE(Label, '(A3)' ) TRIM(LABLES(OrderedList(I)))
!	CALL PLCHLQ(XLoc(OrderedList(I))+(2.0*Dot), OriginY+(1.55*Offset),&
!	& Label ,5*FontScale,0.,1.)
      END DO
      
      WRITE(Buff, 31)  Hierarchy(1)%ResDegree, Hierarchy(N-1)%ResDegree
31    FORMAT ('Range: ',F15.3,' , ',F15.3)
!      CALL PLCHLQ(OriginX+Offset, OriginY+(1.0*Offset),&
!	&  TRIM(Buff),9.0*FontScale,0.,-1.)
      WRITE(Buff, 32) CophCorrCoef
32    FORMAT ('Cophenetic Corr. Coeff. = ',F4.2)
!      CALL PLCHLQ(OriginX+Offset, OriginY+(.5*Offset),&
!	&  TRIM(Buff),9.*FontScale,0.,-1.)
      ! set the line width to smaller
!      CALL GSLWSC(1.)
      ! draw the tree
      prevY = OriginY
      Y = OriginY + (2.0*Offset)
      DO I=1, N-1
	! in case of tie, do not update Y
        IF(Hierarchy(I)%ResDegree /= prevDegree)  THEN
	   IF (Dissim)  THEN
	   	Y = OriginY + (2.0*Offset) + Hierarchy(I)%ResDegree*VUnit
           ELSE
                Y = OriginY+(2.0*Offset)+(1.0 - Hierarchy(I)%ResDegree)*VUnit
	   END IF
        END IF
	J = Hierarchy(I)%MergedCluster1
	K = Hierarchy(I)%MergedCluster2
	! NCARG code to draw line
!	CALL LINE(XLoc(J), YLoc(J), XLoc(J), Y)
!	CALL LINE(XLoc(K), YLoc(K), XLoc(K), Y)
!	CALL LINE(XLoc(J), Y, XLoc(K), Y)
        IF(Hierarchy(I)%ResDegree /= prevDegree)  THEN
		! draw tick on the scale, and display degree
!		CALL LINE(OriginX+(2.0*Offset-(3.0*Dot)), Y, &
!		&         OriginX+(2.0*Offset+(3.0*Dot)), Y)
		! if degree close to each other, they will be
		! overlaped on the plot. Skip if overlap
		IF (Dissim)  THEN
			WRITE(Degree, '(F15.3)' ) Hierarchy(I)%ResDegree
		ELSE
			WRITE(Degree, '(F15.5)' ) Hierarchy(I)%ResDegree
		END IF
                IF ( Y-prevY > (5.0*Dot) ) THEN
!                        CALL PLCHLQ(OriginX+(2.0*Offset-(6.0*Dot)), &
!                        & Y, Degree,6.0*FontScale,0.,1.)
                        prevY = Y
                END IF
	END IF
	! update X, Y locations
	XLoc(J) = XLoc(J) + (XLoc(K)-XLoc(J))/2.0
	YLoc(J) = Y
	prevDegree = Hierarchy(I)%ResDegree
      END DO
      ! draw line from last merge to the top
!      CALL LINE(XLoc(J), Y, XLoc(J), Y+0.4*Offset )	
END SUBROUTINE TreeLabelInTwoLines

! ===========================================================
! Subroutine Name: OrderClusters
! Author: Ahmad Alhamad
! Order the clusters to appear nice on the tree
! ===========================================================
SUBROUTINE OrderClusters(N, Hierarchy, OrderedList)
        USE Global
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N
        TYPE(Merge), DIMENSION(N-1), INTENT(IN) :: Hierarchy
        INTEGER, DIMENSION(N), INTENT(OUT) :: OrderedList
        INTEGER, ALLOCATABLE, DIMENSION(:,:) :: Clusters
        INTEGER, ALLOCATABLE, DIMENSION(:) :: Members      
        INTEGER :: AllocatStatus, I, J, K, P, Q, NDX

        ALLOCATE( Clusters(N,N), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "(Clusters) ** No Enough Memory **"
        ALLOCATE( Members(N), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "(Members) ** No Enough Memory **"
        DO I=1, N
                Clusters(I, 1) = I
                Members(I) = 1
        END DO
        DO P=1, N-1
                J = Hierarchy(P)%MergedCluster1
                K = Hierarchy(P)%MergedCluster2
                DO Q=1,    Members(K)
                        Members(J) = Members(J) + 1
                        NDX = Members(J)
                        Clusters(J, NDX) = Clusters(K, Q)
                END DO
        END DO
        DO I=1, N
                OrderedList(I) = Clusters(1, I)
        END DO
END SUBROUTINE OrderClusters
! ====================================================================
! Subroutine Name: LabelModelRun
! Author: Nusrat Yussouf
! Date: September, 2002
! This subroutine returns the label of X axis in text format. The labels
! are named according to the 23 model runs of the Current project.
! ====================================================================
SUBROUTINE LabelModelRun(N, OrderedList, XLabel1, XLabel2, &
           & XLabel3, XLabel4)
        USE Global
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: N
        INTEGER, DIMENSION(N), INTENT(IN) :: OrderedList
        CHARACTER*3, DIMENSION(N), INTENT(OUT) :: XLabel1, & 
         & XLabel2,XLabel3,XLabel4
        INTEGER :: I

        DO I=1, N
      	    IF (OrderedList(I) == 1) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 1) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 1) XLabel3(I) = 'c'
      	    IF (OrderedList(I) == 1) XLabel4(I) = '0'
      	    IF (OrderedList(I) == 2) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 2) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 2) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 2) XLabel4(I) = '1'
      	    IF (OrderedList(I) == 3) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 3) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 3) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 3) XLabel4(I) = '1'
      	    IF (OrderedList(I) == 4) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 4) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 4) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 4) XLabel4(I) = '2'
      	    IF (OrderedList(I) == 5) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 5) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 5) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 5) XLabel4(I) = '2'
      	    IF (OrderedList(I) == 6) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 6) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 6) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 6) XLabel4(I) = '3'
      	    IF (OrderedList(I) == 7) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 7) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 7) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 7) XLabel4(I) = '3'
      	    IF (OrderedList(I) == 8) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 8) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 8) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 8) XLabel4(I) = '4'
      	    IF (OrderedList(I) == 9) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 9) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 9) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 9) XLabel4(I) = '4'
      	    IF (OrderedList(I) == 10) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 10) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 10) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 10) XLabel4(I) = '5'
      	    IF (OrderedList(I) == 11) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 11) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 11) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 11) XLabel4(I) = '4'
      	    IF (OrderedList(I) == 12) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 12) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 12) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 12) XLabel4(I) = '6'
      	    IF (OrderedList(I) == 13) XLabel1(I) = 'N'
      	    IF (OrderedList(I) == 13) XLabel2(I) = 'B'
      	    IF (OrderedList(I) == 13) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 13) XLabel4(I) = '6'
      	    IF (OrderedList(I) == 14) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 14) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 14) XLabel3(I) = 'c'
      	    IF (OrderedList(I) == 14) XLabel4(I) = '0'
      	    IF (OrderedList(I) == 15) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 15) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 15) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 15) XLabel4(I) = '1'
      	    IF (OrderedList(I) == 16) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 16) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 16) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 16) XLabel4(I) = '1'
      	    IF (OrderedList(I) == 17) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 17) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 17) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 17) XLabel4(I) = '2'
      	    IF (OrderedList(I) == 18) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 18) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 18) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 18) XLabel4(I) = '2'
      	    IF (OrderedList(I) == 19) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 19) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 19) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 19) XLabel4(I) = '3'
      	    IF (OrderedList(I) == 20) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 20) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 20) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 20) XLabel4(I) = '3'
      	    IF (OrderedList(I) == 21) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 21) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 21) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 21) XLabel4(I) = '4'
      	    IF (OrderedList(I) == 22) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 22) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 22) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 22) XLabel4(I) = '4'
      	    IF (OrderedList(I) == 23) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 23) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 23) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 23) XLabel4(I) = '5'
      	    IF (OrderedList(I) == 24) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 24) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 24) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 24) XLabel4(I) = '5'
      	    IF (OrderedList(I) == 25) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 25) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 25) XLabel3(I) = 'n'
      	    IF (OrderedList(I) == 25) XLabel4(I) = '6'
      	    IF (OrderedList(I) == 26) XLabel1(I) = 'A'
      	    IF (OrderedList(I) == 26) XLabel2(I) = 'W'
      	    IF (OrderedList(I) == 26) XLabel3(I) = 'p'
      	    IF (OrderedList(I) == 26) XLabel4(I) = '6'
        END DO
END SUBROUTINE LabelModelRun
! ===========================================================
