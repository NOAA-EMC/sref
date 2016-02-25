!kftime,mem Hierarchical Cluster Analysis
! *********************************************************************
! File Name: main.f
! Author: Ahmad Alhamed
!         University of Oklahoma
!         School of Computer Science 
! Date: 2000 
! Modified by: Nusrat Yussouf
!              CIMMS/OU/NSSL
! Date: September, 2002
! Purpose: This program performs Hierarchical Cluster Analysis. It calls  
! other subroutines and functions from other files, resemblance.f,
! agglomerative.f, cophenetic.f, tree.f, and CAlib.f. When compiling 
! this program, it should be linked with the NCAR graphics library 
! which is used to draw the clustering trees. 
! *********************************************************************
!=========================================================================
! This is the main function. This function calls UserInterface subroutine
! to obtain input from the user.
!=========================================================================

PROGRAM HCA
        IMPLICIT NONE
	CALL UserInterface
END PROGRAM HCA

! =======================================================================
! This subroutine displays the user interface, and obtain all required 
! input.
! =======================================================================
SUBROUTINE UserInterface
        USE Global
	IMPLICIT NONE
	
	!Variable declarations
	INTEGER :: MainMenu, NormalizeMenu, ResCoeffMenu, &
		& ClusteringMethodMenu, SimCoeffMenu, SimDisimMenu  ! functions
	INTEGER :: M, N, Choice, Coef, Method, Normalize, IERR, &
	        &  SimDisim, AllocatStatus, I, J
	CHARACTER(LEN=50) :: OutputFile, InputFile
	CHARACTER(LEN=4) :: Prefix
	CHARACTER(LEN=70) :: Title
         integer nm,nfile,ntime,igx,igy,ngrid,jf,njt,ncluster,nt,k
        parameter (nm=26,nfile=nm+1,ntime=30,igx=185,igy=129,ngrid=185*129)
        parameter (jf=185*129)
        parameter (ncluster = 6)
          logical lb(jf)
         character*6 mem(nm)
         character*1 clus
         character*2 cyc
         character*2 hr 
         character*30 fname3
         integer cluster( ncluster, nm ) 
          integer kftime(nfile,ntime)
          integer iret
!         integer jpds(25),jgds(22),kpds(25),kgds(22)
          integer jpds(200),jgds(200),kpds(200),kgds(200)
            real dmax,dmin 
        integer yy1,mm1,dd1,hh1,iar,meth
        integer yy3(ntime),mm3(ntime),dd3(ntime),hh3(ntime)
         real pr(ncluster, igx,igy),prc(ngrid)
         REAL fcst(nm,igx,igy)
         integer kl,nf,iclus,igrid,kf

        namelist /namin/yy1,mm1,dd1,hh1,iar,meth
         common /c1/yy1,mm1,dd1,hh1,iar,meth
         common /c2/yy3,mm3,dd3,hh3,kftime,mem
         common /c3/kpds,jpds,kgds,jgds,ierr,iret
         common /c4/lb
        Title = 'ensemble'
	! main menu
!	Choice = MainMenu()
        Choice = 1
!	PRINT *
!	PRINT '($,A)','Enter No. of column in Data Matrix:  '
!	READ *, N
!	PRINT '($,A)','Enter No. of attributes OR No. of rows in Data Matrix:  '
!        READ *, M
        N = nm
        M = 185*129
        InputFile = 'dummy'
	! normalization menu
!	Normalize = NormalizeMenu()
        Normalize = 3
	! Res. coeff. menu
!	SimDisim = SimDisimMenu()
        SimDisim = 1
        Coef = 6
	Method = 3

!	IF (SimDisim == 1) THEN
!	   Coef = SimCoeffMenu()
!           Coef = 6
!	ELSE
!           Coef = ResCoeffMenu()
!	END IF	
!	! clustering method menu
!	Method = ClusteringMethodMenu()
!	Method = 3
!	PRINT '($,A)','Enter prefix for File Name To Save Tree on PS format:'
!	READ *, Prefix

        Prefix = 'sref'
	ALLOCATE( LABLES(N), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(LABLES) ** No EnoughMemory **"
	DO I=1,N
		WRITE(LABLES(I), '(I3)' ) I
	END DO

!!!!!  read namin & 
      read(5,namin,end=1000)
1000  continue
      print*,yy1,mm1,dd1,hh1,iar,meth
         read(7,*) (yy3(i),i=1,ntime),(mm3(i),i=1,ntime),  &
      &   (dd3(i),i=1,ntime),(hh3(i),i=1,ntime)
      print*,(yy3(i),i=1,ntime),(mm3(i),i=1,ntime),       &
      &   (dd3(i),i=1,ntime),(hh3(i),i=1,ntime)
      read(7,*) ((kftime(i,j),j=1,ntime),i=1,nfile)
      do i=1,nfile
       print *,(kftime(i,j),j=1,ntime)
      enddo
      read(7,*) (mem(i),i=1,nm)
       print *,(mem(i),i=1,nm)

!!!!!!!!!!!!!!!!!!

	OutputFormat = OTOMOP

	OutputFile = TRIM(Prefix)//'_clustering_process'
        
	OPEN (UNIT=7,FILE=OutputFile,STATUS='NEW',IOSTAT=IERR)
	OPEN (UNIT=777,FILE=TRIM(Prefix)//'_cluster_info',STATUS='NEW', &
        & IOSTAT=IERR)


       	IF (IERR /= 0) THEN
        	PRINT 7 , OutputFile, IERR
7               FORMAT (' ERROR OPENING UNIT=7, FILE NAME ='&
		& ,A13,', IOSTAT = ',I8)
                STOP 7
        END IF
          write(777,100) hh1,mm1,dd1,yy1
!         write(777,100) 'Forecast clusters for',hh1,'z',mm1,'/',dd1,'/',yy1,'of 26-member SREF using OU method (fixed to 6 clusters)'
100   format(22hForecast clusters for ,i2.2,2hz ,i2.2,1h/,i2.2,1h/,i4, &
      56h of 26-member SREF using OU method (fixed to 6 clusters))
!         write(777,*) (i,mem(i),i=1,nm)
          write(777,*) (i,mem(i),i=1,1)
            do njt = 1, ntime
          write(7,*)'************************************************'
          write(7,*) 'fcst hour=',(njt-1)*3,'Begin:::'
          write(7,*) '************************************************'
!         write(777,*)'************************************************'
!         write(777,*) 'fcst hour=',(njt-1)*3,'Begin:::'
          write(777,*) '************************************************'
          write(777,778) (njt-1)*3
          write(777,*) '************************************************'
778     format(1x,i2.2,14h Forecast Hour)

	! initialize counters to 1
	AllHierarchyCntr = 1
	AllCophCoefCntr = 1

	SELECT CASE(OutputFormat)
	CASE (OTOMOP)
	        ALLOCATE( AllHierarchy(N-1), STAT = AllocatStatus )
        	IF ( AllocatStatus /= 0 )  STOP  "(AllHierarchy) &
		& ** No EnoughMemory **"
	        ALLOCATE( AllCophCoef(1), STAT = AllocatStatus )
        	IF ( AllocatStatus /= 0 )  STOP  "(AllCophCoef) &
		& ** No EnoughMemory **"
		Call Start(Title, Choice, N, M, Normalize, Coef, &
		& Method, InputFile, njt, cluster, fcst )      
	CASE (OTAMOP)
	        ALLOCATE( AllHierarchy( (N-1)*NoOfMethods ), &
		& STAT = AllocatStatus )
        	IF ( AllocatStatus /= 0 )  STOP  "(AllHierarchy) &
		& ** No EnoughMemory **"
	        ALLOCATE( AllCophCoef(NoOfMethods), STAT = AllocatStatus )
        	IF ( AllocatStatus /= 0 )  STOP  "(AllCophCoef) &
		& ** No EnoughMemory **"
		Call Start(Title, Choice, N, M, Normalize, Coef, &
		& Method, InputFile, njt, cluster,fcst)
	END SELECT

!!!!!!!! get average field
         do 8000 kl =1, ncluster

         do 405 i=1,igx
         do 405 j=1,igy
         pr(kl,i,j)=0.0
405      continue

         nf=0
         do 6100 k=1,nm
         if(cluster(kl,k)== 0) goto 6100
         nf=nf+1
         do 6101 i=1,igx
         do 6101 j=1,igy
         pr(kl,i,j)=pr(kl,i,j)+fcst(cluster(kl,k),i,j)
6101     continue
!        print*,'nf,kl,pr(kl,20,20)=',nf,kl,pr(kl,20,20)
6100     continue

!        print*,'nf,kl,pr(kl,20,20)=',nf,kl,pr(kl,20,20)
         do 6102 i=1,igx
         do 6102 j=1,igy
          pr(kl,i,j)=pr(kl,i,j)/float(nf)
            if (pr(kl,i,j) > 150000.0) then
                  print*,'sth wrong!',pr(kl,i,j),kl,i,j
             endif

             if (pr(kl,i,j) < 20000.0) then
                  print*,'sth wrong too small!',pr(kl,i,j),kl,i,j
             endif

6102     continue
              print*,'pr=',pr(kl,20,20)

8000     continue
! Write out cluster mean in grib format (array pr). Each cluster in its own separate file
         do 888 k=1,6
          iclus=50+k
!         iclus=70+k

           i=0
           j=1
           do 97 igrid=1,ngrid
            prc(igrid)=0.0
97         continue
           do 98 igrid=1,ngrid
            i=i+1
            prc(igrid)=pr(k,i,j)
            if(i == igx) then
             i=0
             j=j+1
            endif
98        continue

           do igrid =1,ngrid
            if (prc(igrid) > 150000.0 ) then
                print*, prc(igrid),igrid,'too big'
             endif 

            if (prc(igrid) < 20000.0 ) then
                print*, prc(igrid),igrid,'too small'
             endif 

            end do 
               print*,'prc=', prc(1000)

! Output in grib format
      write(clus,'(i1)') k
      write(hr,'(i2.2)') kftime(1,njt)
      write(cyc,'(i2.2)') hh1
      fname3='ref.t'//cyc//'z.pgrb212.clus'//clus//'mean.f'//hr
       print*,'before='
      call baopen(iclus,fname3,ierr)
       print*,'after='

      kpds(1)=7
      kpds(2)=130
      kpds(3)=212
!     kpds(4)=0
      kpds(4)=128
!     kpds(5)=variable dependent
!     kpds(6)=variable dependent
!     kpds(7)=variable dependent
      kpds(8)=yy1-(yy1/100)*100
      kpds(9)=mm1
      kpds(10)=dd1
      kpds(11)=hh1
      kpds(12)=0
      kpds(13)=1
!     kpds(14)=variable dependent
!     kpds(15)=variable dependent
!     kpds(16)=variable dependent
      kpds(17)=-1
      kpds(18)=1
      kpds(19)=2
      kpds(20)=-1
      kpds(21)=21
!     kpds(22)=variable dependent (units decimal scale factor/precision)
!     kpds(23)=2  !if ensemble data
      kpds(23)=0  !if not ensemble data
      kpds(24)=128
      kpds(25)=-1


      kgds(20)=255
!-----------------z500 or mslp---------------------
      if(iar.eq.1) then !500z
       jpds(5)=7
       jpds(6)=100
       jpds(7)=500
      else
       jpds(5)=2        !mslp
       jpds(6)=102
       jpds(7)=0
      endif
      kpds(14)=kftime(1,njt)

         print*,'kpds(14)=',kpds(14)

      kpds(15)=0
      kpds(16)=0
      kpds(22)=1

      print*,'almost there!'

      print*, 'before putgb'
      print*, iclus,jf,iret,prc(100),prc(200),lb(100),lb(200)
 
      call putgb(iclus,jf,kpds,kgds,lb,prc,iret)
888      continue
      call baclose(iclus,ierr)

!!!!    

	 
	DEALLOCATE( AllHierarchy, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 ) STOP "(AllHierarchy) **DeallocationError**"
	DEALLOCATE( AllCophCoef, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 ) STOP "(AllCophCoef) **DeallocationError**"
          write(7,*)  'fcst hour=',(njt-1)*3,':::::end!!!' 
!         write(777,*)  'fcst hour=',(njt-1)*3,':::::end!!!' 
!         write(777,*)  
          enddo



	CLOSE (UNIT=7)
	CLOSE (UNIT=777)

END SUBROUTINE UserInterface
! =======================================================================
! =======================================================================
SUBROUTINE Start(Title, Choice, N, M, Normalize, Coef, Method, &
		& InputFileName,njt,cluster,fcst)  
	USE Global
	IMPLICIT NONE
           integer nm,nfile,ntime,igx,igy,ncluster
        parameter (nm=26,nfile=nm+1,ntime=30,igx=185,igy=129,ncluster=6)
	CHARACTER(LEN=*), INTENT(IN) :: InputFileName
	CHARACTER(LEN=*), INTENT(IN) :: Title
	INTEGER, INTENT(IN) :: N, M, Choice, Normalize, Coef, Method
	INTEGER :: AllocatStatus, I, J, NUM,NJT
	REAL, ALLOCATABLE, DIMENSION(:,:) :: DataMatrix
	REAL, ALLOCATABLE, DIMENSION(:)   :: ResemblanceMatrix
         REAL, DIMENSION(nm,igx,igy), INTENT(OUT) ::  fcst
         integer, dimension(ncluster,nm), INTENT(OUT) :: cluster

	CHARACTER(LEN=10) :: CoefAbbr
        integer yy1,mm1,dd1,hh1,iar,meth
        integer yy3(ntime),mm3(ntime),dd3(ntime),hh3(ntime)
        integer kftime(nfile,ntime)
          integer ifile, ipt  
        character*6 mem
           common /c1/yy1,mm1,dd1,hh1,iar,meth
           common /c2/yy3,mm3,dd3,hh3,kftime,mem
 
	DataWritten = .FALSE.
	ResWritten = .FALSE.	
	WRITE (7, *) '     ', Title
	WRITE (7, *) '------------------------------------------------------------&
	             &-------'
	WRITE (7, *)
	IF (Choice == 1 .OR. Choice == 2)  THEN
		ALLOCATE( DataMatrix(M,N), STAT = AllocatStatus )
		IF (AllocatStatus/= 0 )  STOP  "(DataMatrix) ** No Enough Memory **"
	END IF
        ALLOCATE( ResemblanceMatrix(N*(N-1)/2), STAT = AllocatStatus )
        IF (AllocatStatus/= 0 )  STOP  "(ResemblanceMatrix) ** No Enough Memory **" 
	SELECT CASE(Choice)
	CASE (1)
		!READ DATA MATRIX FROM INPUT FILE
	 CALL ReadDataMatrixFromFile(M, N, DataMatrix, InputFileName,njt)
	CASE (2)
		!READ DATA MATRIX FROM KEYBOARD
		PRINT *,'Enter Data Matrix elements in rowwise order:'
		READ *, ((DataMatrix(I,J), J=1,N),I=1,M)
	CASE (3)
		!READ RESEMBLANCE MATRIX FROM INPUT FILE
		CALL ReadResemblanceMatrixFromFile(N, ResemblanceMatrix, &
		& InputFileName)
	CASE (4)
		!READ RESEMBLANCE MATRIX FROM KEYBOARD
		PRINT *,'Enter Lower Triangular (Off diagnal) of Resemblance'
		PRINT *,'Matrix elements in columnwise order:'
		READ *, (ResemblanceMatrix(I),I=1, N*(N-1)/2 )
	CASE (5)
		STOP
	END SELECT
         print*, 'OK till now!(1)'
!!! convert forecast from dataMatrix(ngrid,nm) to fcst(nm,igx,igy)
       do ifile=1,nm
        i=0
        j=1
        do 40 ipt=1,M
         i=i+1
         fcst(ifile,i,j)=DataMatrix(ipt,ifile)
         if(i == igx) then
          i=0
          j=j+1
         endif
40      continue
       enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

           print *, 'ok till now (2)'



	IF (Choice==1.OR.Choice==2) THEN
		! print data matrix, and calculate resemblance matrix
		IF(M<=20.AND.N<=10.AND.(.NOT.DataWritten))   THEN
			CALL PrintDataMatrix(M, N, DataMatrix)
			DataWritten = .TRUE.
		END IF
		! normalize if requested
		SELECT CASE(Normalize)
		CASE (1)
			CALL CenterByColMean(M, N, DataMatrix)
		CASE (2)
			CALL NormalizeByColMinMax(M, N, DataMatrix)
		CASE (3)
			CALL StdNormalize(M, N, DataMatrix)
		END SELECT
		! calculating res matrix
		CALL CalculateResemblanceMatrix(Coef, M, N, DataMatrix, &
						& ResemblanceMatrix)
	END IF
        ! print resemblance matrix
	IF(.NOT.ResWritten) THEN
		CALL PrintResemblanceMatrix(N, ResemblanceMatrix)
		ResWritten = .TRUE.
	END IF
! ***************************************************
! if coefficient is correlation, we need the absolute
! values only. That's because the sign does not matter
! with respect to the significance of similarity.
	IF (Coef == 4)  THEN
		ResemblanceMatrix = ABS(ResemblanceMatrix)
	END IF
! ***************************************************
	SELECT CASE (OutputFormat)
	CASE (OTOMOP)
                 print *,'OK for now!!' 
		CALL DoClustering(N, Coef, Method, ResemblanceMatrix, cluster)
	CASE (OTAMSP, OTAMOP)
		DO I=1, 7
         	CALL DoClustering(N, Coef, I, ResemblanceMatrix, cluster)
		END DO
	CASE (ATOMOP)
		CALL DoClustering(N, Coef, Method, ResemblanceMatrix, cluster)
	END SELECT
	IF (Choice == 1 .OR. Choice == 2)  THEN
        	DEALLOCATE( DataMatrix, STAT = AllocatStatus )
        	IF ( AllocatStatus /= 0 )  STOP  "(DataMatrix)** Deallocation Error**"
	END IF
        DEALLOCATE( ResemblanceMatrix, STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(ResemblanceMatrix)** Deallocation Error**"
END SUBROUTINE Start
! =======================================================================
INTEGER FUNCTION MainMenu()
        IMPLICIT NONE

	INTEGER :: Choice

	PRINT *,'----------------------------------------------'
	PRINT *,'                Main Menu'
	PRINT *,'(1) Read Data Matrix from Input file'
	PRINT *,'(2) Read Data Matrix from Keyboard'
	PRINT *,'(3) Read Resemblance Matrix from Input file'
	PRINT *,'(4) Read Resemblance Matrix from Keyboard'
	PRINT *,'(5) Exit'
	PRINT *,'----------------------------------------------'
	PRINT '($,A)','    Enter Choice ===>  '
	DO
		READ *, Choice
		IF (Choice >= 1 .AND. Choice <= 5)  EXIT
	END DO
	IF (Choice == 5) STOP
	MainMenu = choice
END FUNCTION MainMenu
! =======================================================================
INTEGER FUNCTION SimDisimMenu()
        IMPLICIT NONE

	INTEGER :: Choice

	PRINT *,'----------------------------------------------'
	PRINT *,'  Similarity or Dissimilarity Measure'
	PRINT *,'(1) Similarity Measure'
	PRINT *,'(2) Dissimilarity Measure'
	PRINT *,'----------------------------------------------'
	PRINT '($,A)','    Enter Choice ===>  '
	DO
		READ *, Choice
		IF (Choice >= 1 .AND. Choice <= 2)  EXIT
	END DO
	SimDisimMenu = Choice
END FUNCTION SimDisimMenu
! =======================================================================
INTEGER FUNCTION SimCoeffMenu()
        IMPLICIT NONE

	INTEGER :: Coef

	PRINT *
        PRINT *,'-------------------------------------'
	PRINT *,' Similarity Measure: Resemblance Coefficient:'
	PRINT *,'    (1) Covariance Coefficient'
        PRINT *,'    (2) Correlation Coefficient'	
	PRINT *,'    (3) Euclidean Similarity'
	PRINT *,'    (4) Cosine Coefficient'
        PRINT *,'-------------------------------------'
	PRINT '($,A)','    Enter Choice ===>   '
        DO
                READ *, Coef
                IF (Coef >= 1 .AND. Coef <= 4)  EXIT
        END DO
        Coef = 4 + Coef
        SimCoeffMenu = Coef
END FUNCTION SimCoeffMenu
! =======================================================================
INTEGER FUNCTION ResCoeffMenu()
        IMPLICIT NONE

	INTEGER :: Coef

	PRINT *
        PRINT *,'-------------------------------------'
	PRINT *,'Resemblance Coefficient:'
	PRINT *,'    (1) Euclidean Distance'		
	PRINT *,'    (2) Square Euclidean Distance'	
	PRINT *,'    (3) Manhattan Distance'
	PRINT *,'    (4) Sup Distance'
        PRINT *,'-------------------------------------'
	PRINT '($,A)','    Enter Choice ===>   '
        DO
                READ *, Coef
                IF (Coef >= 1 .AND. Coef <= 4)  EXIT
        END DO
        ResCoeffMenu = Coef
END FUNCTION ResCoeffMenu
! =======================================================================
INTEGER FUNCTION ClusteringMethodMenu()
        IMPLICIT NONE

	INTEGER :: Method

        PRINT *,'-------------------------------------'
        PRINT *,'Clustering Method:'
        PRINT *,'    (1) SLINK (Single Linkage)'
        PRINT *,'    (2) CLINK (Complete Linkage)'
        PRINT *,'    (3) UPGMA (group average)'
	PRINT *,'    (4) WPGMA (weighted average)'
	PRINT *,'    (5) UPGMC (unweighted centroid)'
	PRINT *,'    (6) WPGMC (weighted centroid)'
	PRINT *,"    (7) Ward's method (minimum variance)"    
        PRINT *,'-------------------------------------'
	PRINT '($,A)','    Enter Choice ===>   '
        DO
        	READ *, Method
                IF (Method >= 1 .AND. Method <= 7)  EXIT
        END DO                                          
        ClusteringMethodMenu = Method
END FUNCTION ClusteringMethodMenu
! =======================================================================
SUBROUTINE DoClustering(N, Coef, Method, ResemblanceMatrix,cluster)
        USE Global
        IMPLICIT NONE

        INTEGER, INTENT(IN) :: N, Coef, Method
        REAL, DIMENSION(N*(N-1)/2), INTENT(IN) :: ResemblanceMatrix
	TYPE(Merge), ALLOCATABLE, DIMENSION(:) :: Hierarchy
	INTEGER :: I, AllocatStatus
	LOGICAL :: Dissim
	REAL :: CutAt, CophCorrCoef
              integer ncluster, nm
            parameter (ncluster=6,nm=26)
          integer, dimension(ncluster,nm), INTENT(OUT) :: cluster 

	CHARACTER(LEN=10) :: CoefAbbr, MethodAbbr
        ALLOCATE( Hierarchy(N-1), STAT = AllocatStatus )
        IF ( AllocatStatus /= 0 )  STOP  "(Hierarchy) ** No EnoughMemory **"
	CALL WhatCoefficient(Coef,CoefAbbr,Dissim)
	CALL WhatClusteringMethod( Method, MethodAbbr)
	WRITE (7, *) 'Clustering Method:  ' , MethodAbbr
	CALL Agglomerative(N, Method, Coef, Dissim, ResemblanceMatrix, &
			& Hierarchy,cluster)
	! add this Hierarchy to AllHierarchy
	DO I = 1, N-1
		AllHierarchy(AllHierarchyCntr) = Hierarchy(I)
		AllHierarchyCntr = AllHierarchyCntr + 1
	END DO
	! compute cophenetic correlation
	CALL Cophenetic(N, Hierarchy, ResemblanceMatrix, CophCorrCoef)
	! add this CophCorrCoef to AllCophCoef
	AllCophCoef(AllCophCoefCntr) = CophCorrCoef
	AllCophCoefCntr = AllCophCoefCntr + 1

	DEALLOCATE( Hierarchy, STAT = AllocatStatus )
	IF ( AllocatStatus /= 0 )  STOP  "(Hierarchy) ** DeallocationError**"
END SUBROUTINE DoClustering
! =======================================================================
! =======================================================================
! ========================================================
! =======================================================================
SUBROUTINE WhatCoefficient(CoefIndex,CoefAbb,Dis)
        IMPLICIT NONE

	INTEGER, INTENT(IN) :: CoefIndex
        CHARACTER(LEN=*), INTENT(OUT) :: CoefAbb
        LOGICAL, INTENT(OUT) :: Dis

        SELECT CASE(CoefIndex)
        CASE (1)
		CoefAbb = 'Euc'
		Dis = .TRUE.
        CASE (2)
		CoefAbb = 'sqEuc'
		Dis = .TRUE.
        CASE (3)
		CoefAbb = 'Manhattan'
		Dis = .TRUE.
        CASE (4)
		CoefAbb = 'Sup Distance'
		Dis = .TRUE.
	CASE (5)
		CoefAbb = 'Cov'
		Dis = .FALSE.	
        CASE (6)
		CoefAbb = 'Corr'
		Dis = .FALSE.
        CASE (7)
		CoefAbb = 'EucSimi'
		Dis = .FALSE.
	CASE (8)
		CoefAbb = 'Cosine'
		Dis = .FALSE.
        END SELECT 
END SUBROUTINE WhatCoefficient
! =======================================================================
SUBROUTINE WhatClusteringMethod( Method, ClusteringMethod)
	IMPLICIT NONE

	INTEGER, INTENT(IN) :: Method
	CHARACTER(LEN=*), INTENT(OUT) :: ClusteringMethod

        SELECT CASE(Method)
        CASE (1)
                ClusteringMethod = "SLINK"
        CASE (2)
                ClusteringMethod = "CLINK"
        CASE (3)
                ClusteringMethod = "UPGMA"
        CASE (4)
                ClusteringMethod = "WPGMA"
        CASE (5)
                ClusteringMethod = "UPGMC"
        CASE (6)
                ClusteringMethod = "WPGMC"
        CASE (7)
                ClusteringMethod = "Ward's"
        END SELECT 
END SUBROUTINE WhatClusteringMethod
! =======================================================================
! ========================================================================
REAL FUNCTION FindMaxOfAllHierarch(Size)
        USE Global
        IMPLICIT NONE

	INTEGER, INTENT(IN) :: Size
	REAL :: MaxDegree
	INTEGER :: I

	MaxDegree = AllHierarchy(1)%ResDegree
	DO I=2, Size
		IF( AllHierarchy(I)%ResDegree > MaxDegree ) THEN
			MaxDegree = AllHierarchy(I)%ResDegree
		END IF
	END DO
        FindMaxOfAllHierarch = MaxDegree
END FUNCTION FindMaxOfAllHierarch
! =======================================================================

