!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains subroutines related to the Shot Gather
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine SG_function()

USE mod_parfile
USE mod_input_arrays
USE mod_SG_arrays
USE mod_PG_arrays

implicit none
include 'mpif.h'

integer :: numtasks,rank,ierr,status(MPI_STATUS_SIZE)

integer :: iDC
integer :: k1,k2,i1,i2,itimes
integer :: nPG,ntimes,nsamples
integer :: i,j,k,is,jrec,ishot
integer :: icount,ishotcount
integer :: perc_25,perc_50,perc_75

real :: SG(nt,NumRec)
real, allocatable :: trace(:)

DOUBLE PRECISION :: start, end

call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
ierr=0;

if(rank.eq.floor(numtasks/2.))start=MPI_Wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     MAIN LOOP DC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(trace(nt))

nPG=NumMax_PG

nsamples=nPG
ntimes=nsamples

if(numtasks.gt.1)	then
        if(nsamples.gt.numtasks)        then
                ntimes=ceiling(1.*nsamples/numtasks)
        endif
	if(nsamples.le.numtasks)        then
                ntimes=1
        endif
endif

perc_25=nsamples/4
perc_50=nsamples/2
perc_75=3*nsamples/4

do itimes=1,ntimes

        icount=(itimes-1)*numtasks+rank+1

        i1=(itimes-1)*numtasks+1
        i2=itimes*numtasks
        if(i2.gt.nsamples)i2=nsamples

        if(icount.eq.perc_25) write(*,*)
        if(icount.eq.perc_25) write(*,*)'25% DONE'
	if(icount.eq.perc_50) write(*,*)'50% DONE'
	if(icount.eq.perc_75) write(*,*)'75% DONE'

	if(icount.le.nsamples)     then

                do is=1,NumShots_PG(icount)

		        ishot=shot_grid_PG(is,icount)
			jrec=trace_grid_PG(is,icount)

			call WRITE_DATA_trace(rank,icount,is,ishot,jrec)

		enddo

        endif !icount

enddo !itimes

deallocate(trace)

call MPI_barrier(MPI_COMM_WORLD,ierr)

if(rank.eq.floor(numtasks/2.))   then
        end   = MPI_Wtime()
        write(*,*)
        write(*,*) 'The conversion from point gathers to shot gathers took',nint((end-start)/60),'minutes'
endif


end subroutine SG_function
