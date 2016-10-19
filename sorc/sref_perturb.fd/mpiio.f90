MODULE module_mpiio

   IMPLICIT NONE

!! FROM MODULE_KINDS

!   The numerical data types defined in this module are:
!      i_byte    - specification kind for byte (1-byte) integer variable
!      i_short   - specification kind for short (2-byte) integer variable
!      i_long    - specification kind for long (4-byte) integer variable
!      i_llong   - specification kind for double long (8-byte) integer variable
!      r_single  - specification kind for single precision (4-byte) real variable
!      r_double  - specification kind for double precision (8-byte) real variable
!      r_quad    - specification kind for quad precision (16-byte) real variable
!
!      i_kind    - generic specification kind for default integer
!      r_kind    - generic specification kind for default floating point
!
!
! Integer type definitions below

! Integer types
  integer, parameter, public  :: i_byte  = selected_int_kind(1)      ! byte  integer
  integer, parameter, public  :: i_short = selected_int_kind(4)      ! short integer
  integer, parameter, public  :: i_long  = selected_int_kind(8)      ! long  integer
  integer, parameter, private :: llong_t = selected_int_kind(16)     ! llong integer
  integer, parameter, public  :: i_llong = max( llong_t, i_long )

! Expected 8-bit byte sizes of the integer kinds
  integer, parameter, public :: num_bytes_for_i_byte  = 1
  integer, parameter, public :: num_bytes_for_i_short = 2
  integer, parameter, public :: num_bytes_for_i_long  = 4
  integer, parameter, public :: num_bytes_for_i_llong = 8

! Define arrays for default definition
  integer, parameter, private :: num_i_kinds = 4
  integer, parameter, dimension( num_i_kinds ), private :: integer_types = (/ &
       i_byte, i_short, i_long,  i_llong  /) 
  integer, parameter, dimension( num_i_kinds ), private :: integer_byte_sizes = (/ &
       num_bytes_for_i_byte, num_bytes_for_i_short, &
       num_bytes_for_i_long, num_bytes_for_i_llong  /)

! Default values
! **** CHANGE THE FOLLOWING TO CHANGE THE DEFAULT INTEGER TYPE KIND ***
  integer, parameter, private :: default_integer = 2  ! 1=byte, 
                                                      ! 2=short, 
                                                      ! 3=long, 
                                                      ! 4=llong
  integer, parameter, public  :: i_kind = integer_types( default_integer )
  integer, parameter, public  :: num_bytes_for_i_kind = &
       integer_byte_sizes( default_integer )


! Real definitions below

! Real types
  integer, parameter, public  :: r_single = selected_real_kind(6)  ! single precision
  integer, parameter, public  :: r_double = selected_real_kind(15) ! double precision
  integer, parameter, private :: quad_t   = selected_real_kind(20) ! quad precision
  integer, parameter, public  :: r_quad   = max( quad_t, r_double )

! Expected 8-bit byte sizes of the real kinds
  integer, parameter, public :: num_bytes_for_r_single = 4
  integer, parameter, public :: num_bytes_for_r_double = 8
  integer, parameter, public :: num_bytes_for_r_quad   = 16

! Define arrays for default definition
  integer, parameter, private :: num_r_kinds = 3
  integer, parameter, dimension( num_r_kinds ), private :: real_kinds = (/ &
       r_single, r_double, r_quad    /) 
  integer, parameter, dimension( num_r_kinds ), private :: real_byte_sizes = (/ &
       num_bytes_for_r_single, num_bytes_for_r_double, &
       num_bytes_for_r_quad    /)

! Default values
! **** CHANGE THE FOLLOWING TO CHANGE THE DEFAULT REAL TYPE KIND ***
  integer, parameter, private :: default_real = 2  ! 1=single, 
                                                   ! 2=double, 
!! END FROM MODULE_KINDS

CONTAINS

subroutine retrieve_index(index,string,varname_all,nrecs,iret)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    retrieve_index  get record number of desired variable
!   prgmmr: parrish          org: np22                date: 2004-11-29
!
! abstract: by examining previously generated inventory of wrf binary restart file,
!             find record number that contains the header record for variable
!             identified by input character variable "string".
!
! program history log:
!   2004-11-29  parrish
!
!   input argument list:
!     string           - mnemonic for variable desired
!     varname_all      - list of all mnemonics obtained from inventory of file
!     nrecs            - total number of sequential records counted in wrf
!                        binary restart file
!
!   output argument list:
!     index            - desired record number
!     iret             - return status, set to 0 if variable was found,
!                        non-zero if not.
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

  implicit none

  integer,intent(out)::iret
  integer,intent(in)::nrecs
  integer,intent(out):: index
  character(*), intent(in):: string
  character(132),intent(in)::varname_all(nrecs)

  integer i

  iret=0

  do i=1,nrecs
   if(trim(string) == trim(varname_all(i))) then
      index=i
      return
   end if
  end do

  write(6,*)' problem reading wrf nmm binary file, rec id "',trim(string),'" not found'

  iret=-1

end subroutine retrieve_index
subroutine next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    next_buf    bring in next direct access block
!   prgmmr: parrish          org: np22                date: 2004-11-29
!
! abstract: bring in next direct access block when needed, as the file is scanned
!             from beginning to end during counting and inventory of records.
!             (subroutines count_recs_wrf_binary_file and inventory_wrf_binary_file)
!
! program history log:
!   2004-11-29  parrish
!
!   input argument list:
!     in_unit          - fortran unit number where input file is opened through.
!     nextbyte         - byte number from beginning of file that is desired
!     locbyte          - byte number from beginning of last block read for desired byt
!     lrecl            - direct access block length
!     nreads           - number of blocks read before now (for diagnostic information
!     lastbuf          - logical, if true, then no more blocks, so return
!
!   output argument list:
!     buf              - output array containing contents of next block
!     locbyte          - byte number from beginning of new block read for desired byte
!     thisblock        - number of new block being read by this routine
!     nreads           - number of blocks read now (for diagnostic information only)
!     lastbuf          - logical, if true, then at end of file.
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
!  use kinds, only: i_byte,i_llong
  implicit none

  integer(i_llong) lrecl
  integer in_unit,nreads
  integer(i_byte) buf(lrecl)
  integer(i_llong) nextbyte,locbyte,thisblock
  logical lastbuf

  integer ierr

  if(lastbuf) return

  ierr=0
  nreads=nreads+1

!  compute thisblock:

  thisblock = 1_i_llong + (nextbyte-1_i_llong)/lrecl

  locbyte = 1_i_llong+mod(locbyte-1_i_llong,lrecl)

  read(in_unit,rec=thisblock,iostat=ierr)buf
  lastbuf = ierr /= 0

end subroutine next_buf

subroutine inventory_wrf_binary_file(in_unit,wrf_ges_filename,nrecs, &
                                     datestr_all,varname_all,domainend_all, &
                                     start_block,end_block,start_byte,end_byte,file_offset)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    inventory_wrf_binary_file  get contents of wrf binary file
!   prgmmr: parrish          org: np22                date: 2004-11-29
!
! abstract: generate list of contents and map of wrf binary file which can be
!             used for reading and writing with mpi io routines.
!             same basic routine as count_recs_wrf_binary_file, except
!             now wrf unpacking routines are used to decode wrf header
!             records, and send back lists of variable mnemonics, dates,
!             grid dimensions, and byte addresses relative to start of
!             file for each field (this is used by mpi io routines).
!
! program history log:
!   2004-11-29  parrish
!
!
!   input argument list:
!     in_unit          - fortran unit number where input file is opened through.
!     wrf_ges_filename - filename of input wrf binary restart file
!     nrecs            - number of sequential records found on input wrf binary restart file.
!                          (obtained by a previous call to count_recs_wrf_binary_file)
!
!   output argument list:  (all following dimensioned nrecs)
!     datestr_all      - date character string for each field, where applicable (or else blanks)
!     varname_all      - wrf mnemonic for each variable, where applicable (or blank)
!     domainend_all    - dimensions of each field, where applicable (up to 3 dimensions)
!     start_block      - direct access block number containing 1st byte of record
!                            (after 4 byte record mark)
!     end_block        - direct access block number containing last byte of record
!                            (before 4 byte record mark)
!     start_byte       - relative byte address in direct access block of 1st byte of record
!     end_byte         - relative byte address in direct access block of last byte of record
!     file_offset      - absolute address of byte before 1st byte of record (used by mpi io)
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$
!   use kinds, only: r_single,i_byte,i_long,i_llong
  use module_internal_header_util
  implicit none

  integer,intent(in)::in_unit,nrecs
  character(*),intent(in)::wrf_ges_filename
  character(132),intent(out)::datestr_all(nrecs),varname_all(nrecs)
  integer,intent(out)::domainend_all(3,nrecs)
  integer,intent(out)::start_block(nrecs),end_block(nrecs)
  integer,intent(out)::start_byte(nrecs),end_byte(nrecs)
  integer(i_llong),intent(out)::file_offset(nrecs)

  integer irecs
  integer(i_llong) nextbyte,locbyte,thisblock
  integer(i_byte) lenrec4(4)
  integer(i_long) lenrec,lensave
  equivalence (lenrec4(1),lenrec)
  integer(i_byte) missing4(4)
  integer(i_long) missing
  equivalence (missing,missing4(1))
  integer(i_llong),parameter:: lrecl=2**20
  integer(i_byte) buf(lrecl)
  integer i,loc_count,nreads
  logical lastbuf
  integer(i_byte) hdrbuf4(2048)
  integer(i_long) hdrbuf(512)
  equivalence(hdrbuf(1),hdrbuf4(1))
  integer,parameter:: int_field       =       530
  integer,parameter:: int_dom_ti_char =       220
  integer,parameter:: int_dom_ti_real =       140
  integer,parameter:: int_dom_ti_integer =       180
  integer hdrbufsize
  integer inttypesize
  integer datahandle,count
  character(128) element,dumstr,strdata
  integer loccode
  character(132) blanks
  integer typesize
  integer fieldtype,comm,iocomm
  integer domaindesc
  character(132) memoryorder,stagger,dimnames(3)
  integer domainstart(3),domainend(3)
  integer memorystart(3),memoryend(3)
  integer patchstart(3),patchend(3)
  character(132) datestr,varname
  real(r_single) dummy_field(1)
!  integer dummy_field
!  real dummy_field
  integer itypesize
  integer idata(1)
  real rdata(100)

  call wrf_sizeof_integer(itypesize)
  inttypesize=itypesize

  blanks=trim(' ')

  open(in_unit,file=trim(wrf_ges_filename),access='direct',recl=lrecl)
  irecs=0
  missing=-9999
  nextbyte=0_i_llong
  locbyte=lrecl
  nreads=0
  lastbuf=.false.
  do

!   get length of next record

    do i=1,4
     nextbyte=nextbyte+1_i_llong
     locbyte=locbyte+1_i_llong
     if(locbyte > lrecl .and. lastbuf) go to 900
     if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
     lenrec4(i)=buf(locbyte)
    end do
    if(lenrec <= 0 .and. lastbuf) go to 900
    if(lenrec <= 0 .and. .not. lastbuf) go to 885
    nextbyte=nextbyte+1_i_llong
    locbyte=locbyte+1_i_llong
    if(locbyte > lrecl .and. lastbuf) go to 900
    if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)

    irecs=irecs+1
    start_block(irecs)=thisblock
    start_byte(irecs)=locbyte
    file_offset(irecs)=nextbyte-1_i_llong
    hdrbuf4(1)=buf(locbyte)
    hdrbuf4(2:4)=missing4(2:4)
    hdrbuf4(5:8)=missing4(1:4)
    datestr_all(irecs)=blanks
    varname_all(irecs)=blanks
    domainend_all(1:3,irecs)=0

    loc_count=1
    do i=2,8
       if(loc_count.ge.lenrec) exit
       loc_count=loc_count+1
       nextbyte=nextbyte+1_i_llong
       locbyte=locbyte+1_i_llong
       if(locbyte > lrecl .and. lastbuf) go to 900
       if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
       hdrbuf4(i)=buf(locbyte)
    end do

         if(lenrec==2048) write(6,*)' irecs,hdrbuf(2),int_dom_ti_char,int_field=', &
                                      irecs,hdrbuf(2),int_dom_ti_char,int_field
    if(lenrec==2048.and.(hdrbuf(2) == int_dom_ti_char .or. hdrbuf(2) == int_field &
    .or. hdrbuf(2) == int_dom_ti_real .or. hdrbuf(2) == int_dom_ti_integer)) then

!    bring in next full record, so we can unpack datestr, varname, and domainend
       do i=9,lenrec
          loc_count=loc_count+1
          nextbyte=nextbyte+1_i_llong
          locbyte=locbyte+1_i_llong
          if(locbyte > lrecl .and. lastbuf) go to 900
          if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
          hdrbuf4(i)=buf(locbyte)
       end do

       if(hdrbuf(2) == int_dom_ti_char) then

          call int_get_ti_header_char(hdrbuf,hdrbufsize,inttypesize, &
                   datahandle,element,dumstr,strdata,loccode)
          varname_all(irecs)=trim(element)
          datestr_all(irecs)=trim(strdata)
              write(6,*)' irecs,varname,datestr = ',irecs,trim(varname_all(irecs)),trim(datestr_all(irecs))

       else if(hdrbuf(2) == int_dom_ti_real) then

          call int_get_ti_header_real(hdrbuf,hdrbufsize,inttypesize,typesize, &
                   datahandle,element,rdata,count,loccode)
          varname_all(irecs)=trim(element)
!          datestr_all(irecs)=trim(strdata)
              write(6,*)' irecs,varname,datestr = ',irecs,trim(varname_all(irecs)),rdata(1:count)
         
       else if(hdrbuf(2) == int_dom_ti_integer) then

          call int_get_ti_header_integer(hdrbuf,hdrbufsize,inttypesize,typesize, &
                   datahandle,element,idata,count,loccode)
          varname_all(irecs)=trim(element)
!          datestr_all(irecs)=trim(strdata)
              write(6,*)' irecs,varname,datestr = ',irecs,trim(varname_all(irecs)),idata(1:count)

       else

          call int_get_write_field_header(hdrbuf,hdrbufsize,inttypesize,typesize, &
             datahandle,datestr,varname,dummy_field,fieldtype,comm,iocomm, &
             domaindesc,memoryorder,stagger,dimnames, &
             domainstart,domainend,memorystart,memoryend,patchstart,patchend)
          varname_all(irecs)=trim(varname)
          datestr_all(irecs)=trim(datestr)
          domainend_all(1:3,irecs)=domainend(1:3)
              write(6,*)' irecs,datestr,domend,varname = ', &
                  irecs,trim(datestr_all(irecs)),domainend_all(1:3,irecs),trim(varname_all(irecs))

       end if
    end if

    nextbyte=nextbyte-loc_count+lenrec
    locbyte=locbyte-loc_count+lenrec
    if(locbyte > lrecl .and. lastbuf) go to 900
    if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
    end_block(irecs)=thisblock
    end_byte(irecs)=locbyte
    lensave=lenrec
    do i=1,4
     nextbyte=nextbyte+1_i_llong
     locbyte=locbyte+1_i_llong
     if(locbyte > lrecl .and. lastbuf) go to 900
     if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
     lenrec4(i)=buf(locbyte)
    end do
    if(lenrec /= lensave) go to 890

  end do

880  continue
     write(6,*)' reached impossible place in inventory_wrf_binary_file'
     close(in_unit)
     return

885  continue
     write(6,*)' problem in inventory_wrf_binary_file, lenrec has bad value before end of file'
     write(6,*)'     lenrec =',lenrec
     close(in_unit)
     return

890  continue
     write(6,*)' problem in inventory_wrf_binary_file, beginning and ending rec len words unequal'
     write(6,*)'     begining reclen =',lensave
     write(6,*)'       ending reclen =',lenrec
     write(6,*)'               irecs =',irecs
     write(6,*)'               nrecs =',nrecs
     close(in_unit)
     return

900  continue
     write(6,*)' normal end of file reached in inventory_wrf_binary_file'
     write(6,*)'        nblocks=',thisblock
     write(6,*)'          irecs,nrecs=',irecs,nrecs
     write(6,*)'         nreads=',nreads
     close(in_unit)

end subroutine inventory_wrf_binary_file

SUBROUTINE wrf_sizeof_integer( retval )
  IMPLICIT NONE
  INTEGER retval
! 4 is defined by CPP
  retval = 4
  RETURN
END SUBROUTINE wrf_sizeof_integer

SUBROUTINE wrf_sizeof_real( retval )
  IMPLICIT NONE
  INTEGER retval
! 4 is defined by CPP
  retval = 4
  RETURN
END SUBROUTINE wrf_sizeof_real
subroutine count_recs_wrf_binary_file(in_unit,wrf_ges_filename,nrecs)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    count_recs_binary_file  count # recs on wrf binary file
!   prgmmr: parrish          org: np22                date: 2004-11-29
!
! abstract: count number of sequential records contained in wrf binary
!             file.  this is done by opening the file in direct access
!             mode with block length of 2**20, the size of the physical
!             blocks on ibm "blue" and "white" machines.  for optimal
!             performance, change block length to correspond to the
!             physical block length of host machine disk space.
!             records are counted by looking for the 4 byte starting
!             and ending sequential record markers, which contain the
!             record size in bytes.  only blocks are read which are known
!             by simple calculation to contain these record markers.
!             even though this is done on one processor, it is still
!             very fast, and the time will always scale by the number of
!             sequential records, not their size.  this step and the
!             following inventory step consistently take less than 0.1 seconds
!             to complete.
!
! program history log:
!   2004-11-29  parrish
!
!   input argument list:
!     in_unit          - fortran unit number where input file is opened through.
!     wrf_ges_filename - filename of input wrf binary restart file
!
!   output argument list:
!     nrecs            - number of sequential records found on input wrf binary restart fil
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

!   do an initial read through of a wrf binary file, and get total number of sequential fil

!   use kinds, only: r_single,i_byte,i_long,i_llong
  implicit none

  integer,intent(in)::in_unit
  character(*),intent(in)::wrf_ges_filename
  integer,intent(out)::nrecs

  integer(i_llong) nextbyte,locbyte,thisblock
  integer(i_byte) lenrec4(4)
  integer(i_long) lenrec,lensave
  equivalence (lenrec4(1),lenrec)
  integer(i_byte) missing4(4)
  integer(i_long) missing
  equivalence (missing,missing4(1))
  integer(i_llong),parameter:: lrecl=2**20
  integer(i_byte) buf(lrecl)
  integer i,loc_count,nreads
  logical lastbuf

  open(in_unit,file=trim(wrf_ges_filename),access='direct',recl=lrecl)
  nrecs=0
  missing=-9999
  nextbyte=0_i_llong
  locbyte=lrecl
  nreads=0
  lastbuf=.false.
  do

!   get length of next record

    do i=1,4
     nextbyte=nextbyte+1_i_llong
     locbyte=locbyte+1_i_llong
     if(locbyte > lrecl .and. lastbuf) go to 900
     if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
     lenrec4(i)=buf(locbyte)
    end do
    if(lenrec <= 0 .and. lastbuf) go to 900
    if(lenrec <= 0 .and. .not.lastbuf) go to 885
    nextbyte=nextbyte+1_i_llong
    locbyte=locbyte+1_i_llong
    if(locbyte > lrecl .and. lastbuf) go to 900
    if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)

    nrecs=nrecs+1

    loc_count=1
    do i=2,4
       if(loc_count.ge.lenrec) exit
       loc_count=loc_count+1
       nextbyte=nextbyte+1_i_llong
       locbyte=locbyte+1_i_llong
       if(locbyte > lrecl .and. lastbuf) go to 900
       if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
    end do
    do i=1,4
       if(loc_count.ge.lenrec) exit
       loc_count=loc_count+1
       nextbyte=nextbyte+1_i_llong
       locbyte=locbyte+1_i_llong
       if(locbyte > lrecl .and. lastbuf) go to 900
       if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
    end do
    nextbyte=nextbyte-loc_count+lenrec
    locbyte=locbyte-loc_count+lenrec
    if(locbyte > lrecl .and. lastbuf) go to 900
    if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
    lensave=lenrec
    do i=1,4
     nextbyte=nextbyte+1_i_llong
     locbyte=locbyte+1_i_llong
     if(locbyte > lrecl .and. lastbuf) go to 900
     if(locbyte > lrecl) call next_buf(in_unit,buf,nextbyte,locbyte,thisblock,lrecl,nreads,lastbuf)
     lenrec4(i)=buf(locbyte)
    end do
    if(lenrec /= lensave) go to 890

  end do

880  continue
     write(6,*)' reached impossible place in count_recs_wrf_binary_file'
     close(in_unit)
     return

885  continue
     write(6,*)' problem in count_recs_wrf_binary_file, lenrec has bad value before end of file'
     write(6,*)'     lenrec =',lenrec
     close(in_unit)
     return

890  continue
     write(6,*)' problem in count_recs_wrf_binary_file, beginning and ending rec len words unequal'
     write(6,*)'     begining reclen =',lensave
     write(6,*)'       ending reclen =',lenrec
     close(in_unit)
     return

900  continue
     write(6,*)' normal end of file reached in count_recs_wrf_binary_file'
     write(6,*)'        nblocks=',thisblock
     write(6,*)'          nrecs=',nrecs
     write(6,*)'         nreads=',nreads
     close(in_unit)

end subroutine count_recs_wrf_binary_file

subroutine retrieve_field(in_unit,wrfges,out,start_block,end_block,start_byte,end_byte)
!$$$  subprogram documentation block
!                .      .    .                                       .
! subprogram:    retrieve_field  retrieve field from wrf binary file
!   prgmmr: parrish          org: np22                date: 2004-11-29
!
! abstract: still using direct access, retrieve a field from the wrf binary restart file.
!
! program history log:
!   2004-11-29  parrish
!
!   input argument list:
!     in_unit          - fortran unit number where input file is opened through.
!     wrfges - filename of input wrf binary restart file
!     start_block      - direct access block number containing 1st byte of record
!                            (after 4 byte record mark)
!     end_block        - direct access block number containing last byte of record
!                            (before 4 byte record mark)
!     start_byte       - relative byte address in direct access block of 1st byte of record
!     end_byte         - relative byte address in direct access block of last byte of record
!
!   output argument list:
!     out              - output buffer where desired field is deposited
!
! attributes:
!   language: f90
!   machine:  ibm RS/6000 SP
!
!$$$

 ! use kinds, only: i_byte,i_kind
  implicit none

  integer(i_kind),intent(in)::in_unit
  character(50),intent(in)::wrfges
  integer(i_kind),intent(in)::start_block,end_block,start_byte,end_byte
  integer(i_byte),intent(out)::out(*)

  integer(i_kind),parameter:: lrecl=2**20
  integer(i_byte) buf(lrecl)
  integer(i_kind) i,ii,k,ibegin,iend,ierr

  open(in_unit,file=trim(wrfges),access='direct',recl=lrecl)

     write(6,*)' in retrieve_field, start_block,end_block=',start_block,end_block
     write(6,*)' in retrieve_field, start_byte,end_byte=',start_byte,end_byte
  ii=0
  do k=start_block,end_block
     read(in_unit,rec=k,iostat=ierr)buf
     ibegin=1 ; iend=lrecl
     if(k == start_block) ibegin=start_byte
     if(k == end_block) iend=end_byte
     do i=ibegin,iend
        ii=ii+1
        out(ii)=buf(i)
     end do
  end do
  close(in_unit)

end subroutine retrieve_field

END MODULE module_mpiio
