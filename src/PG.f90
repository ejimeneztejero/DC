!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains subroutines related to the Point Gather
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine PG_function()

USE mod_parfile
USE mod_SG_arrays
USE mod_input_arrays
USE mod_PG_arrays

implicit none

integer :: nPG

	call allocate_PG_arrays()

	call PG_parameters(nPG)

	NumMax_PG=nPG
	NumMaxShots_PG=maxval(NumShots_PG(1:nPG))

	call PG_write()

end subroutine PG_function

subroutine PG_parameters(num)

USE mod_parfile
USE mod_SG_arrays
USE mod_input_arrays
USE mod_PG_arrays

implicit none

integer :: nsamples,num
integer :: i,ix,j,k,k1,iPG,jj,ndim
integer :: ntimes,itimes,i1,i2,nh
integer :: ishot,icount,ixcount,ishotcount
integer :: file_PG,pos_byte


nh = size_su_header ! Size su header = 60*4 bytes

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     MAIN LOOP DC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

k1=pos_shot_grid(1)-far_offset_grid-ceiling(far_offset_grid/10.)  !!añado 10% por seguridad
if(k1.lt.1)k1=1

pos_byte=1
file_PG=1
ixcount=1
num=1

do icount=1,NumMax_PG

	ix=icount+k1-1
	ishotcount=0

        do ishot=1,NumShots

        	if(shotID_(ishot).ne.0)     then
	
		do j=1,NumRec

			if(pos_trace_grid(j,ishot).eq.ix)	then

				if(ishotcount.le.NumMaxShots_PG)	then

					ishotcount=ishotcount+1 
		
					shot_grid_PG(ishotcount,ixcount)=ishot
					trace_grid_PG(ishotcount,ixcount)=j
					pos_shot_grid_PG(ishotcount,ixcount)=pos_shot_grid(ishot)
					pos_trace_grid_PG(ishotcount,ixcount)=pos_trace_grid(j,ishot)

				endif

			endif

		enddo!j

		endif!ishot

	enddo!ishot

	if(ishotcount.ne.0)	then

		if(pos_byte.gt.maxbytes) then
			pos_byte=1
			file_PG=file_PG+1
		endif

		NumShots_PG(ixcount)=ishotcount
		PG_location(ixcount)=ix

		pos_byte_PG(ixcount)=pos_byte
		ifile_PG(ixcount)=file_PG

		ixcount=ixcount+1
		pos_byte=pos_byte+ishotcount*(nh+nt)*4

	endif	

enddo !icount

if(ixcount.gt.1)num=ixcount-1

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                        END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine PG_parameters

subroutine PG_write()

USE mod_parfile
USE mod_SG_arrays
USE mod_input_arrays
USE mod_PG_arrays

implicit none
include 'mpif.h'

integer :: rank,numtasks
integer :: ierr,status(MPI_STATUS_SIZE)
integer :: nPG,nsamples
integer :: i,ix,j,k,k1,k2,iPG,jj,iDC
integer :: ntimes,itimes,i1,i2
integer :: perc_25,perc_50,perc_75
integer :: ishot,jrec,icount,ixcount,ishotcount

real(4), allocatable :: source(:,:),trace(:)

DOUBLE PRECISION :: start,end

call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
ierr=0;


if(rank.eq.0)start=MPI_Wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

allocate(trace(nt))
allocate(source(nt,NumMaxShots_PG))
trace=0.;source=0.;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     MAIN LOOP DC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nsamples=NumMax_PG
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

		do ishotcount=1,NumShots_PG(icount)

			ishot=shot_grid_PG(ishotcount,icount)
			jrec=trace_grid_PG(ishotcount,icount)

			call READ_DATA_trace(rank,ishot,jrec,unit_DC1,ifile_su(ishot),pos_byte_su(ishot),trace)

			do k=1,nt
        	        	source(k,ishotcount)=trace(k)
                	enddo

		enddo

               call WRITE_PG(1,icount,NumShots_PG(icount),unit_PG1,source)

	endif !icount
 
enddo !itimes

deallocate(trace)
deallocate(source)

call MPI_barrier(MPI_COMM_WORLD,ierr)

if(rank.eq.0)   then
        end   = MPI_Wtime()
        write(*,*)
	write(*,*) 'To build the point gather took',nint((end-start)/60),'minutes'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                        END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine PG_write
