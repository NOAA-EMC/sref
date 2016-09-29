!************************************************************************
! Filename: Global.f90
! Author: Ahmad Alhamad 
!         University of Oklahoma
!         School of Computer Science
! Date: 2000
! Modified by: Nusrat Yussouf
!              CIMMS/OU/NSSL
! Date: September, 2002
!
! Purpose: This file defines the parameters used by other subroutines. 
!************************************************************************
!
MODULE Global
        IMPLICIT NONE

        TYPE Merge
                INTEGER :: MergedCluster1, MergedCluster2, NewLable
                REAL :: ResDegree
        END TYPE Merge

	LOGICAL :: DataWritten, ResWritten	
	INTEGER :: OutputFormat, TPP 
!				 TPP : Trees Per Page
	INTEGER, PARAMETER :: OTOMOP=1, OTAMSP=2, OTAMOP=3, ATOMOP=4
		! OTOMOP : One Time One Method One Page
		! OTAMSP : One Time All Methods Seperate Pages
		! OTAMOP : One Time All Methods One Page
		! ATOMOP : All Time One Method One Page

	INTEGER, PARAMETER :: NoOfMethods=7
!	INTEGER, PARAMETER :: IX = 145, JX = 115
!	INTEGER, PARAMETER :: Col = 54, Row = 1558, NumFields = 34
!	INTEGER, PARAMETER :: Col = 21, Row = i185*129 , NumFields = 34
	INTEGER :: StartingTime, NoOfTimes, Interval
	TYPE(Merge), ALLOCATABLE, DIMENSION(:) :: AllHierarchy
	REAL, ALLOCATABLE, DIMENSION(:)   :: AllCophCoef
	INTEGER :: AllHierarchyCntr, AllCophCoefCntr
	INTEGER :: ReadLable
	CHARACTER(LEN=3), ALLOCATABLE, DIMENSION(:)   :: LABLES

END MODULE Global
!************************************************************************
