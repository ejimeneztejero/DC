!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains subroutines related to the main features of Downward Continuation
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine DC_function(iOBS,rank,datum_i,datum_f)

USE mod_parfile
USE mod_data_arrays
USE mod_model_arrays

implicit none

integer :: i,j,k,nSS,iOBS,rank
integer :: datum_i(nx),datum_f(nx)

real, allocatable :: input(:,:),output(:,:)

nSS=NumShots

allocate(input(nt,nSS),output(nt,nSS))
input=0.;output=0.;

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     READING DATA
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)write(*,*)'READING OBS FIELD DATA ...'
call READ_DATA(iOBS,nSS,input)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!    PROPAGATION DC	                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)write(*,*)'PROPAGATING ...'
call DC_propagation(nSS,datum_i,datum_f,pos_shot_grid,input,output)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!    SAVING OUTPUT	                !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)write(*,*)'SAVING RESULTS ...'
call WRITE_OUTPUT(iOBS,nSS,output)

deallocate(input,output)


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!                        END
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

end subroutine DC_function

subroutine DC_propagation(nSS,&
datum_i,datum_f,pos_grid,input,output)

USE mod_parfile
USE mod_data_arrays
USE mod_model_arrays

implicit none

integer :: i,j
integer :: nSS
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

                nySou(i)=datum_i(nxSou(i)) !! "acting like a fuente"
                nyRec(i)=datum_f(nxRec(i)) !! "acting like a receiver"

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

	call solver(nSS,nxx,nys,nyr,source,ny_DC,nx_DC,model,dt,nt,dmodel,output)

deallocate(nxSou,nySou,nxRec,nyRec)
deallocate(nxx,nys,nyr)
deallocate(model)
deallocate(source)

end subroutine DC_propagation
