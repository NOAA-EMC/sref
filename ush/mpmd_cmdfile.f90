!
! ifort -o mpmd_cmdfile mpmd_cmdfile.f90 -lmpi
! mpiexec_mpt -prefix "[%g]" -np $PBS_NP ./mpmd_cmdfile
!
!  - Dusan Jovic
!
   program mpmd_cmdfile
     use mpi
     implicit none
     integer :: ierr, cmd_status
     integer :: nrank, nsize, ntask
     character(len=1024) :: cmd
     integer :: system

     call mpi_init(ierr)
     call mpi_comm_rank(mpi_comm_world,nrank,ierr)
     call mpi_comm_size(mpi_comm_world,nsize,ierr)

     open(10,file='cmdfile',form='formatted')

     ntask=0
     do while (.true.)
        read(10,"(a1024)",end=99) cmd
        if ( mod(ntask,nsize) == nrank ) then
!           call system(cmd, status=cmd_status)
           cmd_status = system(cmd)
           if ( cmd_status /= 0 ) then
             write(0,*) "command: ",trim(cmd)," failed. error code ",cmd_status/256
             call mpi_abort(mpi_comm_world,cmd_status/256,ierr)
           end if
        end if
        ntask=ntask+1
     end do

99   continue

     close(10)

     call mpi_finalize(ierr)

     stop
   end program mpmd_cmdfile
