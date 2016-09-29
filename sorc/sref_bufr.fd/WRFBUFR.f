      PROGRAM WRFBUFR
C$$$  MAIN PROGRAM DOCUMENTATION BLOCK
C                .      .    .     
C MAIN PROGRAM: WRFBUFR
C   PRGMMR: PYLE             ORG: EMC/MMB    DATE: 2004-11-25
C     
C ABSTRACT:  
C     THIS PROGRAM DRIVES THE EXTERNAL WRF BUFR POST PROCESSOR.
C     
C PROGRAM HISTORY LOG (FOR THE PROF CODE):
C   99-04-22  T BLACK - ORIGINATOR
C   02-07-01  G MANIKIN - FIXED PROBLEM WITH DHCNVC AND DHRAIN
C                          COMPUTATIONS - SEE COMMENTS BELOW
C   03-04-01  M PYLE - BEGAN CONVERTING FOR WRF
C   04-05-26  M PYLE - MADE CHANGES FOR WRF-NMM
C   04-11-24  M PYLE - ELIMINATED USE OF PARMETA FILE, DIMENSIONS
C                      NOW READ IN FROM WRF OUTPUT FILE, WITH WORKING
C                      ARRAYS ALLOCATED TO NEEDED DIMENSIONS.  UNIFIED
C                      WRF-EM AND WRF-NMM VERSIONS INTO A SINGLE CODE
C                      THAT READS EITHER NETCDF OR BINARY OUTPUT FROM
C                      THE WRF MODEL.
C   05-08-29  M PYLE - ELIMINATE THE NEED TO RETAIN ALL WRF HISTORY FILES.
C   07-08-06  J Du & B Zhou - A new prefilename was defined to correctly
C                      calculate precip rate during INCHOUR interval                   
C   10-04-26  ROGERS - CHANGE ITAG and INCR to MINUTES FOR NMMB
C     
C USAGE:    WRFPOST
C   INPUT ARGUMENT LIST:
C     NONE     
C
C   OUTPUT ARGUMENT LIST: 
C     NONE
C     
C   SUBPROGRAMS CALLED:
C     UTILITIES:
C       NONE
C     LIBRARY:
C       COMMON - CTLBLK
C                RQSTFLD
C     
C   ATTRIBUTES:
C     LANGUAGE: FORTRAN 90
C     MACHINE : IBM RS/6000 SP
C$$$  
C
C
C
C     INCLUDE ARRAY DIMENSIONS.
C      INCLUDE "parmeta"
C      INCLUDE "mpif.h"
C
C     DECLARE VARIABLES.
C     
C     
C     INCLUDE COMMON BLOCKS.
!tst      INCLUDE "CTLBLK.comm"
C     
C     SET HEADER WRITER FLAGS TO TRUE.
c
C
      real rinc(5)
      integer jdate(8),idate(8)
          integer iii,IHR,IMINTOT
      character(len=6) :: IOFORM,model
      character(len=98) :: newname
      character(len=98) :: fileName
      character(len=98) :: prefileName
      character(len=19) :: DateStr
      character(len=3) :: ITAGLAB
      character(len=2) :: IMINLAB
      character(len=2):: hrpiece
      character(len=8):: minpiece
      integer :: DataHandle
!      integer, parameter:: INCR=3

c     print*,'ME,NUM_PROCS,NUM_SERVERS=',ME,NUM_PROCS,NUM_SERVERS

      DATA hrpiece /'h_'/
      DATA minpiece /'m_00.00s'/



C
C**************************************************************************
C
C     START PROGRAM WRFBUFR.
C
	write(6,*) 'to read statements'
       read(5,111) fileName
	write(6,*) 'initial filename= ', filename
       read(5,113) model
	write(6,*) 'model type= ', model
       read(5,113) IOFORM
	write(6,*) 'ioform= ', ioform
       read(5,112) DateStr
	write(6,*) 'datestr= ', datestr
       read(5,*) NFILES
       read(5,*) INCR

       if (model .ne. 'NMMB') then
         read(5,*) IHR
         read(5,111) prefileName
         write(6,*) 'previous filename= ', prefilename
       else
         read(5,*) IMINTOT
       endif



!!!! CHANGE THIS ASSUMPTION???

! assume for now that the first date in the stdin file is the start date

       read(DateStr,300) iyear,imn,iday,ihrst

C      write(*,*) 'in WRFPOST iyear,imn,iday,ihrst',iyear,imn,iday,ihrst
 300  format(i4,1x,i2,1x,i2,1x,i2)
	
	IDATE=0

         IDATE(2)=imn
         IDATE(3)=iday
         IDATE(1)=iyear
         IDATE(5)=ihrst

 111  format(a98)
 112  format(a19)
 113  format(a6)
C

	do N=1,NFILES

	len=index(filename,' ')-1

!	IHR=(N-1)*INCR
        IMIN=0

!	add forecast hour/minutes to start time

	if (model .eq. 'NMMB') then
        if (MOD(IMINTOT,60).EQ.0) THEN
          IHR=IMINTOT/60
          IMIN=0
        else
          IHR=IMINTOT/60
          IMIN=MOD(IMINTOT,60)
        endif
        endif


	RINC(1)=0.
	RINC(2)=float(IHR)
        RINC(3)=float(IMIN)
	RINC(4)=0.
	RINC(5)=0.

	call w3movdat(rinc,idate,jdate)

	if (model(1:4) .eq. 'NCEP') then
	write(DateStr,302) JDATE(1),JDATE(2),JDATE(3),JDATE(5)
	elseif (model(1:4) .eq. 'NCAR') then
	write(DateStr,302) JDATE(1),JDATE(2),JDATE(3),JDATE(5)
	endif

c20080707	filename=filename(1:len-19)//DateStr

!	if (model(1:4) .eq. 'NCAR') then
!	filename=filename(1:len-19)//DateStr
!	filename(len-2:len-2)='_'	
!	filename(len-5:len-5)='_'	
!	endif

 301  format(i4,'-',i2.2,'-',i2.2,'T',i2.2,':00:00')
 302  format(i4,'-',i2.2,'-',i2.2,'_',i2.2,':00:00')

	write(6,*) 'calling prof.... '
	write(6,*) 'datestr: ', datestr
	write(6,*) 'fileName ', fileName
	write(6,*) 'IHR: ', IHR
	write(6,*) '--------------------------------'

	
	if (ioform(1:6) .eq. 'binary') then 

	if (model(1:4) .eq. 'NCEP') then
        CALL PROF_NMM(fileName,prefileName,DateStr,IHR,INCR)
	elseif (model(1:4) .eq. 'NCAR') then
	CALL PROF_EM(fileName,prefileName,DateStr,IHR,INCR)
	endif

	endif

	if (ioform(1:6) .eq. 'netcdf') then 

	if (model(1:4) .eq. 'NCEP') then
        CALL PROF_NMM_NET(fileName,DateStr,IHR,INCR)
	elseif (model(1:4) .eq. 'NCAR') then
	CALL PROF_EM_NET(fileName,DateStr,IHR,INCR)
	endif

	endif

        if (ioform(1:6) .eq. 'nemsio') then
        write(0,*) 'call PROF_NMMB'
        CALL PROF_NMMB(fileName,IMINTOT,INCR)
        endif

	write(6,*) 'back from prof'
	END DO


      STOP0
      END

