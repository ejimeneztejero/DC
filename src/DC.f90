!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains subroutines related to the main features of Downward Continuation
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DC_function(iDC,&
nSources,nSS,&
datum_i,datum_f)

USE mod_parfile
USE mod_input_arrays
USE mod_SG_arrays
USE mod_PG_arrays

implicit none
include 'mpif.h'
integer :: numtasks,rank,ierr,status(MPI_STATUS_SIZE)

integer :: i,j,k,iDC
integer :: nSources,nSS,nsamples
integer :: ntimes,itimes,icount,i1,i2
integer :: datum_i(nx),datum_f(nx)
integer :: unit_file1,unit_file2,pos_byte,ifile,nSSS

integer, allocatable :: pos_grid(:)
real, allocatable :: input(:,:),output(:,:)

DOUBLE PRECISION :: start, end

call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
ierr=0;

if(rank.eq.0)start=MPI_Wtime()

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     MAIN LOOP DC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nsamples=nSources

ntimes=nsamples
if(numtasks.gt.1)	then

        if(nsamples.gt.numtasks)        then
                ntimes=ceiling(1.*nsamples/numtasks)
        endif

	if(nsamples.le.numtasks)        then
                ntimes=1
        endif

endif

if(numtasks.eq.1.and.nsamples.gt.1) write(*,*)'NO PARALELIZATION, BAD IDEA'

do itimes=1,ntimes

icount=(itimes-1)*numtasks+rank+1        !!aquí sucede la paralelización

i1=(itimes-1)*numtasks+1
i2=itimes*numtasks
if(i2.gt.nsamples)i2=nsamples

if(rank.eq.0)   then

write(*,*)
write(*,*)"`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.             "
write(*,*)"`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.             "
write(*,*)"`'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.             "
write(*,*)
write(*,*)'ROUND ',itimes, 'OUT OF', ntimes
if(iDC.eq.1)write(*,*)'DC (STEP 1) OF SHOT GATHERS FROM: ',shotID_su(i1),' TO ',shotID_su(i2)
if(iDC.eq.2)write(*,*)'DC (STEP 2) OF POINT GATHERS FROM: ',i1,' TO ',i2

endif

if(icount.le.nsamples)     then

if(iDC.eq.1) then

	nSSS=nSS
	ifile=ifile_su(icount)
	pos_byte=pos_byte_su(icount)
	unit_file1=unit_DC0
	unit_file2=unit_DC1

	allocate(pos_grid(nSSS))
	pos_grid=0
	pos_grid(1:nSSS)=pos_trace_grid(1:nSSS,icount)

        if(rank.eq.0)write(*,*)'READING SHOT GATHERS ...'

endif

if(iDC.eq.2) then 

	nSSS=NumShots_PG(icount)
	ifile=ifile_PG(icount)
	pos_byte=pos_byte_PG(icount)
	unit_file1=unit_PG1
	unit_file2=unit_PG2

	allocate(pos_grid(nSSS))
	pos_grid=0
	pos_grid(1:nSSS)=pos_shot_grid_PG(1:nSSS,icount)

        if(rank.eq.0)write(*,*)'READING POINT GATHERS ...'

endif

allocate(input(nt,nSSS),output(nt,nSSS))
input=0.;output=0.;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     READING DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call READ_DATA(icount,unit_file1,pos_byte,ifile,nSSS,input)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!    PROPAGATION DC	                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0) write(*,*)'PROPAGATING ...'

call DC_propagation(iDC,nSSS,datum_i,datum_f,pos_grid,input,output)

deallocate(pos_grid)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!    FILTERING OUTPUT	                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!if(filtering.ne.0)	then
!!	if(rank.eq.0) write(*,*)'FILTERING RESULT BETWEEN ...', f1,' Hz AND ',f2,' Hz'
!!	call BP_filter(nSSS,nt,dt,f1,f2,output)
!!endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!    SAVING OUTPUT	                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(iDC.eq.1) then
        if(rank.eq.0)write(*,*)'SAVING SHOT GATHERS ...'
	call WRITE_SG(iDC,icount,0,nSSS,unit_file2,output)
endif

if(iDC.eq.2) then
        if(rank.eq.0)write(*,*)'SAVING POINT GATHERS FINALES ...'
	call WRITE_PG(iDC,icount,nSSS,unit_file2,output)
endif

deallocate(input,output)

endif   !icount

enddo   !itimes

call MPI_barrier(MPI_COMM_WORLD,ierr)

if(rank.eq.0)   then
        end   = MPI_Wtime()
        write(*,*)
        write(*,*) 'DC:',iDC,', took',nint((end-start)/60),'minutes'
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                        END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine DC_function


subroutine DC_propagation(iDC,nSS,&
datum_i,datum_f,pos_grid,input,output)

USE mod_parfile
USE mod_input_arrays
USE mod_SG_arrays
USE mod_PG_arrays

implicit none

integer :: i,j
integer :: iDC,nSS
integer :: datum_i(nx),datum_f(nx)
integer :: x_i,x_f,y_i,y_f
integer :: nx_DC,ny_DC
integer :: ddx,ddy
integer :: pos_grid(nSS)

real :: input(nt,nSS)
real :: output(nt,nSS)

integer, allocatable :: nxSou(:),nySou(:),nxRec(:),nyRec(:)
integer, allocatable :: nxx(:),nys(:),nyr(:)
real, allocatable :: source(:,:),model(:,:)

allocate(nxSou(nSS),nySou(nSS),nxRec(nSS),nyRec(nSS))
allocate(nxx(nSS),nys(nSS),nyr(nSS))
allocate(source(nt,nSS))

ddx=1+floor(added_space_model_X/dmodel) !!added model space in x-axis
ddy=1+floor(added_space_model_Y/dmodel) !!added model space in y-axis (points)
    
!!------ Geometry acquisicion, fuentes, receptores

        do i=1,nSS

                nxSou(i)=pos_grid(i)
                nxRec(i)=nxSou(i)

                nySou(i)=datum_i(nxSou(i)) !! fuentes datum_i
                nyRec(i)=datum_f(nxRec(i)) !! receivers datum_f

	        do j=1,nt
			source(j,i)=input(nt-j+1,i)
        	enddo

        enddo

!!---------- PROPAGATOR FOR A SUPERSHOT

        x_i=minval(nxRec)-ddx+1
	if(x_i.le.0)x_i=1

	x_f=maxval(nxRec)+ddx
	if(x_f.gt.nx)x_f=nx;	

        y_i=1;y_f=maxval(nyRec)+ddy;
        
	nx_DC=x_f-x_i+1;ny_DC=y_f-y_i+1

	allocate(model(ny_DC,nx_DC))
	model(1:ny_DC,1:nx_DC)=DC_model(y_i:y_f,x_i:x_f)

	nxx=nxRec-(x_i-1);
	nys=nySou-(y_i-1);nyr=nyRec-(y_i-1);

	call solver(iDC,nSS,nxx,nys,nyr,source,ny_DC,nx_DC,model,output)

deallocate(nxSou,nySou,nxRec,nyRec)
deallocate(nxx,nys,nyr)
deallocate(model)
deallocate(source)

end subroutine DC_propagation
