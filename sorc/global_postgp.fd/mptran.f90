!-------------------------------------------------------------------------------
subroutine mpstart(mpicomm,mpirank,mpisize,iret)
!$$$  Subprogram documentation block
!
! Subprogram:    mpstart     Initialize a distributed program
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram initializes a distributed program.
!   It returns a communicator and the rank and size of the distribution.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:    call mpstart(mpicomm,mpirank,mpisize,iret)
!   Output argument list:
!     mpicomm  integer(kint_mpi) MPI communicator
!     mpirank  integer(kint_mpi) rank of the process (from mpi_comm_rank)
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     iret     integer(kint_mpi) return code
!
! Subprograms called:
!   mpi_init      MPI initialization
!   mpi_comm_rank get MPI rank
!   mpi_comm_size get MPI size
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(out):: mpicomm,mpirank,mpisize,iret
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  call mpi_init(iret)
  if(iret.ne.0) return
  mpicomm=MPI_COMM_WORLD
  call mpi_comm_rank(mpicomm,mpirank,iret)
  call mpi_comm_size(mpicomm,mpisize,iret)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mpabort(iret)
!$$$  Subprogram documentation block
!
! Subprogram:    mpabort     Abort a distributed program
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram aborts a distributed program with a return code.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:    call mpabort(iret)
!   Input argument list:
!     iret     integer(kint_mpi) return code
!
! Subprograms called:
!   mpi_abort     MPI abort
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: iret
  integer(kint_mpi) ierr
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  call mpi_abort(MPI_COMM_WORLD,iret,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mpstop(iret)
!$$$  Subprogram documentation block
!
! Subprogram:    mpstop      Finalize a distributed program
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram finalizes a distributed program.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:    call mpstop(iret)
!   Output argument list:
!     iret     integer(kint_mpi) return code
!
! Subprograms called:
!   mpi_finalize  MPI finalization
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(out):: iret
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  call mpi_finalize(iret)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mpnexts(mpicomm,mpirank,mpisize,iflag,na,a,iret)
!$$$  Subprogram documentation block
!
! Subprogram:    mpnexts     Send integer data to the next process
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram sends integer data to the next process.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:    call mpnexts(mpicomm,mpirank,mpisize,iflag,na,a,iret)
!   Input argument list:
!     mpicomm  integer(kint_mpi) MPI communicator
!     mpirank  integer(kint_mpi) rank of the process (from mpi_comm_rank)
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     iflag    integer(kint_mpi) cyclic flag (=1 if last sends to the first)
!     na       integer(kint_mpi) number of words to send
!     a        integer(kint_mpi) (na) array of data to send
!   Output argument list:
!     iret     integer(kint_mpi) return code
!
! Subprograms called:
!   mpi_send      MPI send
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpirank,mpisize,iflag,na
  integer(kint_mpi),intent(in):: a(na)
  integer(kint_mpi),intent(out):: iret
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(iflag.eq.1.or.mpirank.lt.mpisize-1)&
  call mpi_send(a,na,MPI_INTEGER,mod(mpirank+1,mpisize),0,mpicomm,iret)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mpnextr(mpicomm,mpirank,mpisize,iflag,na,a,iret)
!$$$  Subprogram documentation block
!
! Subprogram:    mpnextr     Receive integer data from the previous process
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram receives integer data from the previous process.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:    call mpnextr(mpicomm,mpirank,mpisize,iflag,na,a,iret)
!   Input argument list:
!     mpicomm  integer(kint_mpi) MPI communicator
!     mpirank  integer(kint_mpi) rank of the process (from mpi_comm_rank)
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     iflag    integer(kint_mpi) cyclic flag (=1 if first receives from last)
!     na       integer(kint_mpi) number of words to receive
!   Output argument list:
!     a        integer(kint_mpi) (na) array to receive data
!     iret     integer(kint_mpi) return code
!
! Subprograms called:
!   mpi_recv      MPI receive
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpirank,mpisize,iflag,na
  integer(kint_mpi),intent(out):: a(na),iret
  integer(kint_mpi) status(MPI_STATUS_SIZE)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  if(iflag.eq.1.or.mpirank.gt.0)&
  call mpi_recv(a,na,MPI_INTEGER,MPI_ANY_SOURCE,MPI_ANY_TAG,mpicomm,status,iret)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mptgen(mpirank,mpisize,nd,jt1,jt2,j1,j2,jx,jm,jn)
!$$$  Subprogram documentation block
!
! Subprogram:    mptgen      Generate grid decomposition dimensions
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram decomposes total dimensions of a problem
!   into smaller domains to be managed on a distributed memory system.
!   The last dimension given is decomposed first.  If more decompositions
!   are possible, the next to last dimension is decomposed next, and so on.
!   The transpositions between decompositions should be done by mptran*.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:    call mptgen(mpirank,mpisize,nd,jt1,jt2,j1,j2,jx,jm,jn)
!   Input argument list:
!     mpirank  integer(kint_mpi) rank of the process (from mpi_comm_rank)
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     nd       integer(kint_mpi) number of dimensions to decompose
!     jt1      integer(kint_mpi) (nd) lower bounds of total dimensions
!     jt2      integer(kint_mpi) (nd) upper bounds of total dimensions
!   Output argument list:
!     j1       integer(kint_mpi) (nd) lower bounds of local decompositions
!     j2       integer(kint_mpi) (nd) upper bounds of local decompositions
!     jx       integer(kint_mpi) (nd) local size of decompositions
!     jm       integer(kint_mpi) (nd) maximum size of decompositions
!     jn       integer(kint_mpi) (nd) number of decompositions
!
! Attributes:
!   Language: Fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  integer(kint_mpi),intent(in):: mpirank,mpisize,nd,jt1(nd),jt2(nd)
  integer(kint_mpi),intent(out):: j1(nd),j2(nd),jx(nd),jm(nd),jn(nd)
  integer msize,mrank,msn,mrn,n
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  msize=mpisize
  mrank=mpirank
  do n=nd,1,-1
    if(jt2(n).ge.jt1(n)) then
      jm(n)=(jt2(n)-jt1(n))/msize+1
      msn=max(msize/(jt2(n)-jt1(n)+1),1)
      if(n.eq.1) msn=1
      jn(n)=msize/msn
      mrn=mrank/msn
      j1(n)=min(jt1(n)+jm(n)*mrn,jt2(n)+1)
      j2(n)=min(jt1(n)+jm(n)*mrn+jm(n)-1,jt2(n))
      jx(n)=j2(n)-j1(n)+1
      msize=msn
      mrank=mod(mrank,msn)
    else
      jm(n)=0
      jn(n)=1
      j1(n)=jt1(n)
      j2(n)=jt2(n)
      jx(n)=0
    endif
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mpfillc1(mpicomm,mpisize,im,ida,idb,km,kma,kmb,a,b)
!$$$  Subprogram documentation block
!
! Subprogram:    mpfillc1    Gather grid decomposition
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram gathers an array of data from a
!   grid decomposition to its full extent by using message passing.
!   The grid decomposition should be generated by mptgen.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:  call mpfillc1(mpicomm,mpisize,im,ida,idb,km,kma,kmb,a,b)
!   Input argument list:
!     mpicomm  integer(kint_mpi) mpi communicator
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     im       integer(kint_mpi) undecomposed range
!     ida      integer(kint_mpi) undecomposed input dimension
!     idb      integer(kint_mpi) undecomposed output dimension
!     km       integer(kint_mpi) input grid decomposition size
!     kma      integer(kint_mpi) input grid decomposed range
!     kmb      integer(kint_mpi) output grid undecomposed range
!     a        character(1) (ida,kma) input array
!   Output argument list:
!     b        character(1) (idb,kmb) output array
!
! Subprograms called:
!   mpi_allgather mpi exchange of data between every process pair
!
! Remarks:
!   While this routine serves a wide variety of scalable gather functions
!   for multidimensional grids,
!     (a) it does not work with nonrectanguloid grids;
!     (b) it does not do any load balancing;
!     (c) it does not do any communication hiding.
!
!   This subprogram must be used rather than mpi_allgather
!   in any of the following cases
!     (a) The undecomposed range is less than the respective dimension
!         (either im.lt.ida or im.lt.idb)
!     (b) The decomposed range is ever less than the decomposition size
!         (kma.lt.km for any process)
!     (c) The output grid range is not the full extent
!         (kmb.lt.mpisize*km)
!   If none of these conditions apply, mpi_allgather could be used directly
!   rather than this subprogram and would be more efficient.
!
! Attributes:
!   language: fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpisize
  integer(kint_mpi),intent(in):: im,ida,idb
  integer(kint_mpi),intent(in):: km,kma,kmb
  character(1),dimension(ida,kma),intent(in):: a
  character(1),dimension(idb,kmb),intent(out):: b
  character(1),dimension(im,km):: ta
  character(1),dimension(im,km,mpisize):: tb
  integer(kint_mpi)::i,k,l,ierr,kb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally gather input array
!$omp parallel do
  do k=1,kma
    do i=1,im
      ta(i,k)=a(i,k)
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  externally gather data
  call mpi_allgather(ta,im*km,MPI_CHARACTER,tb,im*km,MPI_CHARACTER,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally gather output array
!$omp parallel do private(kb)
  do l=1,mpisize
    do k=1,km
      kb=k+km*(l-1)
      if(kb.le.kmb) then
        do i=1,im
          b(i,kb)=tb(i,k,l)
        enddo
      endif
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mpfilli4(mpicomm,mpisize,im,ida,idb,km,kma,kmb,a,b)
!$$$  Subprogram documentation block
!
! Subprogram:    mpfilli4    Gather grid decomposition
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram gathers an array of data from a
!   grid decomposition to its full extent by using message passing.
!   The grid decomposition should be generated by mptgen.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:  call mpfilli4(mpicomm,mpisize,im,ida,idb,km,kma,kmb,a,b)
!   Input argument list:
!     mpicomm  integer(kint_mpi) mpi communicator
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     im       integer(kint_mpi) undecomposed range
!     ida      integer(kint_mpi) undecomposed input dimension
!     idb      integer(kint_mpi) undecomposed output dimension
!     km       integer(kint_mpi) input grid decomposition size
!     kma      integer(kint_mpi) input grid decomposed range
!     kmb      integer(kint_mpi) output grid undecomposed range
!     a        integer(kint_mpi) (ida,kma) input array
!   Output argument list:
!     b        integer(kint_mpi) (idb,kmb) output array
!
! Subprograms called:
!   mpi_allgather mpi exchange of data between every process pair
!
! Remarks:
!   While this routine serves a wide variety of scalable gather functions
!   for multidimensional grids,
!     (a) it does not work with nonrectanguloid grids;
!     (b) it does not do any load balancing;
!     (c) it does not do any communication hiding.
!
!   This subprogram must be used rather than mpi_allgather
!   in any of the following cases
!     (a) The undecomposed range is less than the respective dimension
!         (either im.lt.ida or im.lt.idb)
!     (b) The decomposed range is ever less than the decomposition size
!         (kma.lt.km for any process)
!     (c) The output grid range is not the full extent
!         (kmb.lt.mpisize*km)
!   If none of these conditions apply, mpi_allgather could be used directly
!   rather than this subprogram and would be more efficient.
!
! Attributes:
!   language: fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpisize
  integer(kint_mpi),intent(in):: im,ida,idb
  integer(kint_mpi),intent(in):: km,kma,kmb
  integer(kint_mpi),dimension(ida,kma),intent(in):: a
  integer(kint_mpi),dimension(idb,kmb),intent(out):: b
  integer(kint_mpi),dimension(im,km):: ta
  integer(kint_mpi),dimension(im,km,mpisize):: tb
  integer(kint_mpi)::i,k,l,ierr,kb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally gather input array
!$omp parallel do
  do k=1,kma
    do i=1,im
      ta(i,k)=a(i,k)
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  externally gather data
  call mpi_allgather(ta,im*km,MPI_INTEGER4,tb,im*km,MPI_INTEGER4,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally gather output array
!$omp parallel do private(kb)
  do l=1,mpisize
    do k=1,km
      kb=k+km*(l-1)
      if(kb.le.kmb) then
        do i=1,im
          b(i,kb)=tb(i,k,l)
        enddo
      endif
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mpfillr4(mpicomm,mpisize,im,ida,idb,km,kma,kmb,a,b)
!$$$  Subprogram documentation block
!
! Subprogram:    mpfillr4    Gather grid decomposition
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram gathers an array of data from a
!   grid decomposition to its full extent by using message passing.
!   The grid decomposition should be generated by mptgen.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:  call mpfillr4(mpicomm,mpisize,im,ida,idb,km,kma,kmb,a,b)
!   Input argument list:
!     mpicomm  integer(kint_mpi) mpi communicator
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     im       integer(kint_mpi) undecomposed range
!     ida      integer(kint_mpi) undecomposed input dimension
!     idb      integer(kint_mpi) undecomposed output dimension
!     km       integer(kint_mpi) input grid decomposition size
!     kma      integer(kint_mpi) input grid decomposed range
!     kmb      integer(kint_mpi) output grid undecomposed range
!     a        real(4) (ida,kma) input array
!   Output argument list:
!     b        real(4) (idb,kmb) output array
!
! Subprograms called:
!   mpi_allgather mpi exchange of data between every process pair
!
! Remarks:
!   While this routine serves a wide variety of scalable gather functions
!   for multidimensional grids,
!     (a) it does not work with nonrectanguloid grids;
!     (b) it does not do any load balancing;
!     (c) it does not do any communication hiding.
!
!   This subprogram must be used rather than mpi_allgather
!   in any of the following cases
!     (a) The undecomposed range is less than the respective dimension
!         (either im.lt.ida or im.lt.idb)
!     (b) The decomposed range is ever less than the decomposition size
!         (kma.lt.km for any process)
!     (c) The output grid range is not the full extent
!         (kmb.lt.mpisize*km)
!   If none of these conditions apply, mpi_allgather could be used directly
!   rather than this subprogram and would be more efficient.
!
! Attributes:
!   language: fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpisize
  integer(kint_mpi),intent(in):: im,ida,idb
  integer(kint_mpi),intent(in):: km,kma,kmb
  real(4),dimension(ida,kma),intent(in):: a
  real(4),dimension(idb,kmb),intent(out):: b
  real(4),dimension(im,km):: ta
  real(4),dimension(im,km,mpisize):: tb
  integer(kint_mpi)::i,k,l,ierr,kb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally gather input array
!$omp parallel do
  do k=1,kma
    do i=1,im
      ta(i,k)=a(i,k)
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  externally gather data
  call mpi_allgather(ta,im*km,MPI_REAL4,tb,im*km,MPI_REAL4,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally gather output array
!$omp parallel do private(kb)
  do l=1,mpisize
    do k=1,km
      kb=k+km*(l-1)
      if(kb.le.kmb) then
        do i=1,im
          b(i,kb)=tb(i,k,l)
        enddo
      endif
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mptran1l1(mpicomm,mpisize,jm,jma,jmb,jda,km,kma,kmb,kdb,a,b)
!$$$  Subprogram documentation block
!
! Subprogram:    mptran1l1   Transpose grid decompositions
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram transposes an array of data from one
!   grid decomposition to another by using message passing.
!   The grid decompositions should be generated by mptgen.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:  call mptran1l1(mpicomm,mpisize,jm,jma,jmb,jda,km,kma,kmb,kdb,a,b)
!   Input argument list:
!     mpicomm  integer(kint_mpi) mpi communicator
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     jm       integer(kint_mpi) output grid decomposition size
!     jma      integer(kint_mpi) input grid undecomposed range
!     jmb      integer(kint_mpi) output grid decomposed range
!     jda      integer(kint_mpi) input grid undecomposed dimension
!     km       integer(kint_mpi) input grid decomposition size
!     kma      integer(kint_mpi) input grid decomposed range
!     kmb      integer(kint_mpi) output grid undecomposed range
!     kdb      integer(kint_mpi) output grid undecomposed dimension
!     a        logical*1 (jda,kma) input array
!   Output argument list:
!     b        logical*1 (kdb,jmb) output array
!
! Subprograms called:
!   mpi_alltoall  mpi exchange of data between every process pair
!   sgetmo        essl transpose
!
! Remarks:
!   While this routine serves a wide variety of scalable transpose functions
!   for multidimensional grids,
!     (a) it does not work with nonrectanguloid grids;
!     (b) it does not do any load balancing;
!     (c) it does not do any communication hiding.
!
!   This subprogram must be used rather than mpi_alltoall
!   in any of the following cases:
!     (a) The decomposition size is greater than one
!         (either km.gt.1 or jm.gt.1)
!     (b) The decomposed range is ever zero
!         (either kma.eq.0 or jmb.eq.0 for any process)
!     (c) The output grid range is not the full extent
!         (either kmb.lt.mpisize or kmb.lt.kda or jma.lt.mpisize or jma.lt.jda)
!   If none of these conditions apply, mpi_alltoall could be used directly
!   rather than this subprogram and would be more efficient.
!
! Attributes:
!   language: fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpisize
  integer(kint_mpi),intent(in):: jm,jma,jmb,jda
  integer(kint_mpi),intent(in):: km,kma,kmb,kdb
  logical*1,dimension(jda,kma),intent(in):: a
  logical*1,dimension(kdb,jmb),intent(out):: b
  logical*1,dimension(jm,km,mpisize):: ta,tb
  integer(kint_mpi)::j,k,l,ierr
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  check single node case
  if(mpisize.eq.1) then
    do k=1,min(km,kma,kmb)
      do j=1,min(jm,jma,jmb)
        b(k,j)=a(j,k)
      enddo
    enddo
    return
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose input array
!$omp parallel do
  do l=1,mpisize
    do k=1,kma
      do j=1,min(jm,jma-jm*(l-1))
        ta(j,k,l)=a(j+jm*(l-1),k)
      enddo
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  externally transpose data
  call mpi_alltoall(ta,jm*km,MPI_LOGICAL1,tb,jm*km,MPI_LOGICAL1,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose output array
!$omp parallel do
  do l=1,mpisize
    do k=1,min(km,kmb-km*(l-1))
      do j=1,jmb
        b(k+km*(l-1),j)=tb(j,k,l)
      enddo
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mptran1r4(mpicomm,mpisize,jm,jma,jmb,jda,km,kma,kmb,kdb,a,b)
!$$$  Subprogram documentation block
!
! Subprogram:    mptran1r4   Transpose grid decompositions
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram transposes an array of data from one
!   grid decomposition to another by using message passing.
!   The grid decompositions should be generated by mptgen.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:  call mptran1r4(mpicomm,mpisize,jm,jma,jmb,jda,km,kma,kmb,kdb,a,b)
!   Input argument list:
!     mpicomm  integer(kint_mpi) mpi communicator
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     jm       integer(kint_mpi) output grid decomposition size
!     jma      integer(kint_mpi) input grid undecomposed range
!     jmb      integer(kint_mpi) output grid decomposed range
!     jda      integer(kint_mpi) input grid undecomposed dimension
!     km       integer(kint_mpi) input grid decomposition size
!     kma      integer(kint_mpi) input grid decomposed range
!     kmb      integer(kint_mpi) output grid undecomposed range
!     kdb      integer(kint_mpi) output grid undecomposed dimension
!     a        real(4) (jda,kma) input array
!   Output argument list:
!     b        real(4) (kdb,jmb) output array
!
! Subprograms called:
!   mpi_alltoall  mpi exchange of data between every process pair
!   sgetmo        essl transpose
!
! Remarks:
!   While this routine serves a wide variety of scalable transpose functions
!   for multidimensional grids,
!     (a) it does not work with nonrectanguloid grids;
!     (b) it does not do any load balancing;
!     (c) it does not do any communication hiding.
!
!   This subprogram must be used rather than mpi_alltoall
!   in any of the following cases:
!     (a) The decomposition size is greater than one
!         (either km.gt.1 or jm.gt.1)
!     (b) The decomposed range is ever zero
!         (either kma.eq.0 or jmb.eq.0 for any process)
!     (c) The output grid range is not the full extent
!         (either kmb.lt.mpisize or kmb.lt.kda or jma.lt.mpisize or jma.lt.jda)
!   If none of these conditions apply, mpi_alltoall could be used directly
!   rather than this subprogram and would be more efficient.
!
!   Example 1.  Transpose a 1000 x 10000 matrix.
!
!!!   include 'mpif.h'                                     ! use mpi
!!!   parameter(jt=1000,kt=10000)                          ! set problem size
!!!   real,allocatable:: a(:,:),b(:,:)                     ! declare arrays
!!!   call mpi_init(ierr)                                  ! initialize mpi
!!!   call mpi_comm_rank(MPI_COMM_WORLD,mpirank,ierr)      ! get mpi rank
!!!   call mpi_comm_size(MPI_COMM_WORLD,mpisize,ierr)      ! get mpi size
!!!   call mptgen(mpirank,mpisize,1,1,jt,j1,j2,jx,jm,jn)   ! decompose output
!!!   call mptgen(mpirank,mpisize,1,1,kt,k1,k2,kx,km,kn)   ! decompose input
!!!   allocate(a(jt,k1:k2),b(kt,j1:j2))                    ! allocate arrays
!!!   a=reshape((/((j+k,j=1,jt),k=k1,k2)/),(/jt,k2-k1+1/)) ! initialize input
!!!   call mptran1r4(MPI_COMM_WORLD,mpisize,               ! transpose arrays
!!!  &               jm,jt,j2-j1+1,jt,km,k2-k1+1,kt,kt,a,b)
!!!   print '(2i8,f16.1)',((k,j,b(k,j),k=2000,kt,2000),    ! print some values
!!!  &                    j=((j1-1)/200+1)*200,j2,200)
!!!   call mpi_finalize(ierr)                              ! finalize mpi
!!!   end
!
!   This 1000x10000 transpose took 0.57 seconds on 2 winterhawk processors;
!   it took 0.31 seconds on 4 processors and 0.18 seconds on 8 processors.
!   A 10000x10000 transpose took 0.84 seconds on 16 winterhawk processors;
!   it took 0.77 seconds on 20 processors and 0.40 seconds on 40 processors.
!
! Attributes:
!   language: fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpisize
  integer(kint_mpi),intent(in):: jm,jma,jmb,jda
  integer(kint_mpi),intent(in):: km,kma,kmb,kdb
  real(4),dimension(jda,kma),intent(in):: a
  real(4),dimension(kdb,jmb),intent(out):: b
  real(4),dimension(jm,km,mpisize):: ta,tb
  integer(kint_mpi)::j,k,l,ierr
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  check single node case
  if(mpisize.eq.1) then
    if(min(km,kma,kmb,jm,jma,jmb).gt.0)&
    call sgetmo(a,jda,min(jm,jma,jmb),min(km,kma,kmb),b,kdb)
    return
  endif
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose input array
!$omp parallel do
  do l=1,mpisize
    do k=1,kma
      do j=1,min(jm,jma-jm*(l-1))
        ta(j,k,l)=a(j+jm*(l-1),k)
      enddo
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  externally transpose data
  call mpi_alltoall(ta,jm*km,MPI_REAL4,tb,jm*km,MPI_REAL4,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose output array
  if(min(km,kmb,jmb).gt.0)&
  call sgetmo(tb,jm,jmb,kmb,b,kdb)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mptranl1(mpicomm,mpisize,im,ida,idb,&
                    jm,jma,jmb,jda,km,kma,kmb,kdb,a,b)
!$$$  Subprogram documentation block
!
! Subprogram:    mptranl1    Transpose grid decompositions
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram transposes an array of data from one
!   grid decomposition to another by using message passing.
!   The grid decompositions should be generated by mptgen.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:  call mptranl1(mpicomm,mpisize,im,ida,idb,&
!                       jm,jma,jmb,jda,km,kma,kmb,kdb,a,b)
!   Input argument list:
!     mpicomm  integer(kint_mpi) mpi communicator
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     im       integer(kint_mpi) undecomposed range
!     ida      integer(kint_mpi) undecomposed input dimension
!     idb      integer(kint_mpi) undecomposed output dimension
!     jm       integer(kint_mpi) output grid decomposition size
!     jma      integer(kint_mpi) input grid undecomposed range
!     jmb      integer(kint_mpi) output grid decomposed range
!     jda      integer(kint_mpi) input grid undecomposed dimension
!     km       integer(kint_mpi) input grid decomposition size
!     kma      integer(kint_mpi) input grid decomposed range
!     kmb      integer(kint_mpi) output grid undecomposed range
!     kdb      integer(kint_mpi) output grid undecomposed dimension
!     a        logical*1 (ida,jda,kma) input array
!   Output argument list:
!     b        logical*1 (idb,kdb,jmb) output array
!
! Subprograms called:
!   mpi_alltoall  mpi exchange of data between every process pair
!
! Remarks:
!   While this routine serves a wide variety of scalable transpose functions
!   for multidimensional grids,
!     (a) it does not work with nonrectanguloid grids;
!     (b) it does not do any load balancing;
!     (c) it does not do any communication hiding.
!
!   This subprogram must be used rather than mpi_alltoall
!   in any of the following cases:
!     (a) The undecomposed range is less than the respective dimension
!         (either im.lt.ida or im.lt.idb)
!     (b) The decomposition size is greater than one
!         (either km.gt.1 or jm.gt.1)
!     (c) The decomposed range is ever zero
!         (either kma.eq.0 or jmb.eq.0 for any process)
!     (d) The output grid range is not the full extent
!         (either kmb.lt.mpisize or kmb.lt.kda or jma.lt.mpisize or jma.lt.jda)
!   If none of these conditions apply, mpi_alltoall could be used directly
!   rather than this subprogram and would be more efficient.
!
! Attributes:
!   language: fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpisize
  integer(kint_mpi),intent(in):: im,ida,idb
  integer(kint_mpi),intent(in):: jm,jma,jmb,jda
  integer(kint_mpi),intent(in):: km,kma,kmb,kdb
  logical*1,dimension(ida,jda,kma),intent(in):: a
  logical*1,dimension(idb,kdb,jmb),intent(out):: b
  logical*1,dimension(im,jm,km,mpisize):: ta,tb
  integer(kint_mpi)::i,j,k,l,ierr,ja,kb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose input array
!$omp parallel do private(ja)
  do l=1,mpisize
    do k=1,kma
      do j=1,jm
        ja=j+jm*(l-1)
        if(ja.le.jma) then
          do i=1,im
            ta(i,j,k,l)=a(i,ja,k)
          enddo
        endif
      enddo
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  externally transpose data
  call mpi_alltoall(ta,im*jm*km,MPI_LOGICAL1,tb,im*jm*km,MPI_LOGICAL1,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose output array
!$omp parallel do private(kb)
  do l=1,mpisize
    do k=1,km
      kb=k+km*(l-1)
      if(kb.le.kmb) then
        do j=1,jmb
          do i=1,im
            b(i,kb,j)=tb(i,j,k,l)
          enddo
        enddo
      endif
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mptranr4(mpicomm,mpisize,im,ida,idb,&
                    jm,jma,jmb,jda,km,kma,kmb,kdb,a,b)
!$$$  Subprogram documentation block
!
! Subprogram:    mptranr4    Transpose grid decompositions
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram transposes an array of data from one
!   grid decomposition to another by using message passing.
!   The grid decompositions should be generated by mptgen.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:  call mptranr4(mpicomm,mpisize,im,ida,idb,&
!                       jm,jma,jmb,jda,km,kma,kmb,kdb,a,b)
!   Input argument list:
!     mpicomm  integer(kint_mpi) mpi communicator
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     im       integer(kint_mpi) undecomposed range
!     ida      integer(kint_mpi) undecomposed input dimension
!     idb      integer(kint_mpi) undecomposed output dimension
!     jm       integer(kint_mpi) output grid decomposition size
!     jma      integer(kint_mpi) input grid undecomposed range
!     jmb      integer(kint_mpi) output grid decomposed range
!     jda      integer(kint_mpi) input grid undecomposed dimension
!     km       integer(kint_mpi) input grid decomposition size
!     kma      integer(kint_mpi) input grid decomposed range
!     kmb      integer(kint_mpi) output grid undecomposed range
!     kdb      integer(kint_mpi) output grid undecomposed dimension
!     a        real(4) (ida,jda,kma) input array
!   Output argument list:
!     b        real(4) (idb,kdb,jmb) output array
!
! Subprograms called:
!   mpi_alltoall  mpi exchange of data between every process pair
!
! Remarks:
!   While this routine serves a wide variety of scalable transpose functions
!   for multidimensional grids,
!     (a) it does not work with nonrectanguloid grids;
!     (b) it does not do any load balancing;
!     (c) it does not do any communication hiding.
!
!   This subprogram must be used rather than mpi_alltoall
!   in any of the following cases:
!     (a) The undecomposed range is less than the respective dimension
!         (either im.lt.ida or im.lt.idb)
!     (b) The decomposition size is greater than one
!         (either km.gt.1 or jm.gt.1)
!     (c) The decomposed range is ever zero
!         (either kma.eq.0 or jmb.eq.0 for any process)
!     (d) The output grid range is not the full extent
!         (either kmb.lt.mpisize or kmb.lt.kda or jma.lt.mpisize or jma.lt.jda)
!   If none of these conditions apply, mpi_alltoall could be used directly
!   rather than this subprogram and would be more efficient.
!
!   Example 1.  Transpose a 1000 x 10000 matrix.
!
!!!   include 'mpif.h'                                     ! use mpi
!!!   parameter(jt=1000,kt=10000)                          ! set problem size
!!!   real,allocatable:: a(:,:),b(:,:)                     ! declare arrays
!!!   call mpi_init(ierr)                                  ! initialize mpi
!!!   call mpi_comm_rank(MPI_COMM_WORLD,mpirank,ierr)      ! get mpi rank
!!!   call mpi_comm_size(MPI_COMM_WORLD,mpisize,ierr)      ! get mpi size
!!!   call mptgen(mpirank,mpisize,1,1,jt,j1,j2,jx,jm,jn)   ! decompose output
!!!   call mptgen(mpirank,mpisize,1,1,kt,k1,k2,kx,km,kn)   ! decompose input
!!!   allocate(a(jt,k1:k2),b(kt,j1:j2))                    ! allocate arrays
!!!   a=reshape((/((j+k,j=1,jt),k=k1,k2)/),(/jt,k2-k1+1/)) ! initialize input
!!!   call mptranr4(MPI_COMM_WORLD,mpisize,1,1,1,          ! transpose arrays
!!!  &              jm,jt,j2-j1+1,jt,km,k2-k1+1,kt,kt,a,b)
!!!   print '(2i8,f16.1)',((k,j,b(k,j),k=2000,kt,2000),    ! print some values
!!!  &                    j=((j1-1)/200+1)*200,j2,200)
!!!   call mpi_finalize(ierr)                              ! finalize mpi
!!!   end
!
!   This transpose took 0.6 seconds on 4 2-way winterhawk nodes.
!   A 20000x10000 transpose took 3.4 seconds on 16 2-way winterhawk nodes.
!   Thus a transpose may take about 1 second for every 16 Mb per node.
!
! Attributes:
!   language: fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpisize
  integer(kint_mpi),intent(in):: im,ida,idb
  integer(kint_mpi),intent(in):: jm,jma,jmb,jda
  integer(kint_mpi),intent(in):: km,kma,kmb,kdb
  real(4),dimension(ida,jda,kma),intent(in):: a
  real(4),dimension(idb,kdb,jmb),intent(out):: b
  real(4),dimension(im,jm,km,mpisize):: ta,tb
  integer(kint_mpi)::i,j,k,l,ierr,ja,kb
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose input array
!$omp parallel do private(ja)
  do l=1,mpisize
    do k=1,kma
      do j=1,jm
        ja=j+jm*(l-1)
        if(ja.le.jma) then
          do i=1,im
            ta(i,j,k,l)=a(i,ja,k)
          enddo
        endif
      enddo
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  externally transpose data
  call mpi_alltoall(ta,im*jm*km,MPI_REAL4,tb,im*jm*km,MPI_REAL4,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose output array
!$omp parallel do private(kb)
  do l=1,mpisize
    do k=1,km
      kb=k+km*(l-1)
      if(kb.le.kmb) then
        do j=1,jmb
          do i=1,im
            b(i,kb,j)=tb(i,j,k,l)
          enddo
        enddo
      endif
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
!-------------------------------------------------------------------------------
subroutine mptrangr4(mpicomm,mpisize,ndi,ndj,ndk,im,isa,isb,&
                     jma,jsa,jmb,jsb,jm,jn,kma,ksa,kmb,ksb,km,kn,a,b)
!$$$  Subprogram documentation block
!
! Subprogram:    mptrangr4   Transpose grid decompositions
!   Prgmmr: Iredell    Org: W/NP23      Date: 1999-02-12
!
! Abstract: This subprogram transposes an array of data from one
!   grid decomposition to another by using message passing.
!   More than one dimension in the input or output array may be decomposed.
!   The respective grid decompositions should be generated by mptgen.
!
! Program history log:
!   1999-02-12  Mark Iredell
!
! Usage:  call mptrangr4(mpicomm,mpisize,ndi,ndj,ndk,im,isa,isb,&
!                        jma,jsa,jmb,jsb,jm,jn,kma,ksa,kmb,ksb,km,kn,a,b)
!   Input argument list:
!     mpicomm  integer(kint_mpi) mpi communicator
!     mpisize  integer(kint_mpi) size of the process (from mpi_comm_size)
!     ndi      integer(kint_mpi) number of undecomposed dimensions
!     ndj      integer(kint_mpi) number of dimensions decomposed in the output grid
!     ndk      integer(kint_mpi) number of dimensions decomposed in the input grid
!     im       integer(kint_mpi) (ndi) undecomposed dimension ranges
!     isa      integer(kint_mpi) (ndi) undecomposed dimension input grid skip values
!     isb      integer(kint_mpi) (ndi) undecomposed dimension output grid skip values
!     jma      integer(kint_mpi) (ndj) input grid undecomposed dimension ranges
!     jsa      integer(kint_mpi) (ndj) input grid undecomposed dimension skip values
!     jmb      integer(kint_mpi) (ndj) output grid decomposed dimension ranges
!     jsb      integer(kint_mpi) (ndj) output grid decomposed dimension skip values
!     jm       integer(kint_mpi) (ndj) output grid decomposition sizes
!     jn       integer(kint_mpi) (ndj) output grid decomposition divisions
!     kma      integer(kint_mpi) (ndk) input grid decomposed dimension ranges
!     ksa      integer(kint_mpi) (ndk) input grid decomposed dimension skip values
!     kmb      integer(kint_mpi) (ndk) output grid undecomposed dimension ranges
!     ksb      integer(kint_mpi) (ndk) output grid undecomposed dimension skip values
!     km       integer(kint_mpi) (ndk) input grid decomposition sizes
!     kn       integer(kint_mpi) (ndk) input grid decomposition divisions
!     a        real(4) (1+dot_product(im-1,isa)+&
!                         dot_product(jma-1,jsa)+&
!                         dot_product(kma-1,ksa)) input array
!   Output argument list:
!     b        real(4) (1+dot_product(im-1,isb)+&
!                         dot_product(jmb-1,jsb)+&
!                         dot_product(kmb-1,ksb)) output array
!
! Subprograms called:
!   mpi_alltoall  mpi exchange of data between every process pair
!
! Remarks:
!   While this routine serves a wide variety of scalable transpose functions
!   for multidimensional grids,
!     (a) it does not work with nonrectanguloid grids;
!     (b) it does not do any load balancing;
!     (c) it does not do any communication hiding.
!
! Remarks:
!   Example 1.  Transpose a 1000 x 10000 matrix.
!
!!!   include 'mpif.h'                                     ! use mpi
!!!   parameter(jt=1000,kt=10000)                          ! set problem size
!!!   real,allocatable:: a(:,:),b(:,:)                     ! declare arrays
!!!   call mpi_init(ierr)                                  ! initialize mpi
!!!   call mpi_comm_rank(MPI_COMM_WORLD,mpirank,ierr)      ! get mpi rank
!!!   call mpi_comm_size(MPI_COMM_WORLD,mpisize,ierr)      ! get mpi size
!!!   call mptgen(mpirank,mpisize,1,1,jt,j1,j2,jx,jm,jn)   ! decompose output
!!!   call mptgen(mpirank,mpisize,1,1,kt,k1,k2,kx,km,kn)   ! decompose input
!!!   allocate(a(jt,k1:k2),b(kt,j1:j2))                    ! allocate arrays
!!!   a=reshape((/((j+k,j=1,jt),k=k1,k2)/),(/jt,k2-k1+1/)) ! initialize input
!!!   call mptranr4g(MPI_COMM_WORLD,mpisize,1,1,1,1,1,1,   ! transpose arrays
!!!  &               jt,1,j2-j1+1,kt,jm,jn,
!!!  &               k2-k1+1,jt,kt,1,km,kn,a,b)
!!!   print '(2i8,f16.1)',((k,j,b(k,j),k=2000,kt,2000),    ! print some values
!!!  &                    j=((j1-1)/200+1)*200,j2,200)
!!!   call mpi_finalize(ierr)                              ! finalize mpi
!!!   end
!
!   This transpose took 2.4 seconds on 4 2-way winterhawk nodes.
!   A 20000x10000 transpose took 12.5 seconds on 16 2-way winterhawk nodes.
!   Thus a transpose may take about 1 second for every 4 Mb per node.
!
! Attributes:
!   language: fortran 90
!
!$$$
  use machine,only:kint_mpi
  implicit none
  include 'mpif.h'
  integer(kint_mpi),intent(in):: mpicomm,mpisize,ndi,ndj,ndk
  integer(kint_mpi),dimension(ndi),intent(in):: im,isa,isb
  integer(kint_mpi),dimension(ndj),intent(in):: jma,jsa,jmb,jsb,jm,jn
  integer(kint_mpi),dimension(ndk),intent(in):: kma,ksa,kmb,ksb,km,kn
  real(4),dimension(1+dot_product(im-1,isa)+&
                      dot_product(jma-1,jsa)+&
                      dot_product(kma-1,ksa)),intent(in):: a
  real(4),dimension(1+dot_product(im-1,isb)+&
                      dot_product(jmb-1,jsb)+&
                      dot_product(kmb-1,ksb)),intent(out):: b
  real(4),dimension(product(im),product(jm),product(km),mpisize):: ta,tb
  integer(kint_mpi)::iq,jq,kq
  integer(kint_mpi)::ip(ndi),jpa(ndj),jpb(ndj),kpa(ndk),kpb(ndk)
  integer(kint_mpi)::ix(ndi),jx(ndj),kx(ndk)
  integer(kint_mpi)::i,j,k,l,ierr,ks,js,is
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  compute dimension products
  iq=product(im)
  jq=product(jm)
  kq=product(km)
  ip=(/(product(im(1:i-1)),i=1,ndi)/)
  jpa=(/(product(jm(1:i-1)*jn(1:i-1)),i=1,ndj)/)
  jpb=(/(product(jm(1:i-1)),i=1,ndj)/)
  kpa=(/(product(km(1:i-1)),i=1,ndk)/)
  kpb=(/(product(km(1:i-1)*kn(1:i-1)),i=1,ndk)/)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose input array
!$omp parallel do private(kx,ks,jx,js,ix,is)
  do l=1,mpisize
    do k=1,kq
      kx=mod((k-1)/kpa(1:ndk),km)
      if(all(kx.lt.kma)) then
        ks=1+dot_product(kx,ksa)
        do j=1,jq
          jx=mod((j+jq*(l-1)-1)/jpa(1:ndj),jm*jn)
          if(all(jx.lt.jma)) then
            js=ks+dot_product(jx,jsa)
            do i=1,iq
              ix=mod((i-1)/ip(1:ndi),im)
              is=js+dot_product(ix,isa)
              ta(i,j,k,l)=a(is)
            enddo
          endif
        enddo
      endif
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  externally transpose data
  call mpi_alltoall(ta,iq*jq*kq,MPI_REAL4,tb,iq*jq*kq,MPI_REAL4,mpicomm,ierr)
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
!  internally transpose output array
!$omp parallel do private(kx,ks,jx,js,ix,is)
  do l=1,mpisize
    do k=1,kq
      kx=mod((k+kq*(l-1)-1)/kpb(1:ndk),km*kn)
      if(all(kx.lt.kmb)) then
        ks=1+dot_product(kx,ksb)
        do j=1,jq
          jx=mod((j-1)/jpb(1:ndj),jm)
          if(all(jx.lt.jmb)) then
            js=ks+dot_product(jx,jsb)
            do i=1,iq
              ix=mod((i-1)/ip(1:ndi),im)
              is=js+dot_product(ix,isb)
              b(is)=tb(i,j,k,l)
            enddo
          endif
        enddo
      endif
    enddo
  enddo
! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
end subroutine
