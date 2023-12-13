!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This is the main file for the calculation of Downward Continuation
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

PROGRAM main

USE mod_parfile
USE mod_input_arrays
USE mod_SG_arrays
USE mod_PG_arrays

implicit none
include 'mpif.h'

integer :: numtasks,rank,ierr,status(MPI_STATUS_SIZE)
integer :: iDC,i,output

character(len=500) :: file_name,command
DOUBLE PRECISION :: start,end

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
ierr=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	ASCII ART WELLCOME MESSAGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(rank.eq.0)	then

		start=MPI_Wtime()
		call ascii_art(1)

        	write(*,*)
		write(*,*)'SHOT_GATHERS SHOULD BE:'
		write(*,*)'- BUTTER-WORTH FILTERED'
		write(*,*)'- MUTED TILL THE SEA-BOTTOM'
		!write(*,*)'- WITH SPHERICAL DIVERGENCE APPLIED (TO KEEP AMPLITUDES ALONG OFFSET)'
        	write(*,*)

	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	INPUTS AND DATA PREPARATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	call inputs()

	call MPI_barrier(MPI_COMM_WORLD,ierr)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	WRITE STREAMER DATA TO DC IN OUTPUT FOLDER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(DC.eq.1)	then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! START DC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'*******************************************'
	if(rank.eq.0)write(*,*)'DOWNWARD CONTINUATION'
	if(rank.eq.0)write(*,*)'*******************************************'

	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'FIRST STEP IN DC: redatuming receivers from datum 1 to datum 2'
	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'

	iDC=1
	call DC_function(iDC,NumShots,NumRec,datum1_rec,datum2)

        if(save_txt.eq.1) then
                if(rank.eq.0)write(*,*)'SAVING TXT DATA'
                call SAVE_SHOTS_TXT(iDC)
        endif

endif!DC=1


if(DC.eq.2)        then

	iDC=2

	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{'
	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'CONVERSION OF SHOT GATHERS TO POINT GATHERS ...'
	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}'

	call PG_function()

	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>'
	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'SECOND STEP IN DC: redatuming sources from datum 1 to datum 2'
	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<'
	if(rank.eq.0)write(*,*)

	call DC_function(iDC,NumMax_PG,NumMaxShots_PG,datum1_shot,datum2)

	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{'
	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'CONVERSION OF POINT GATHER TO SHOT GATHERS ...'
	if(rank.eq.0)write(*,*)
	if(rank.eq.0)write(*,*)'}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}'

	call SG_function()

	call MPI_barrier(MPI_COMM_WORLD,ierr)

        if(save_txt.eq.1) then
                if(rank.eq.0)write(*,*)'SAVING TXT DATA'
                call SAVE_SHOTS_TXT(iDC)
        endif

endif	!DC=2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! 			END DOWNWARD CONTINUATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)   then

        end   = MPI_Wtime()
        write(*,*)
	write(*,*) 'The whole process took',nint((end-start)/60),'minutes'

	call ascii_art(3)

endif

call close_su_files()

call deallocate_SG_arrays
call deallocate_input_arrays
if(DC.eq.2)call deallocate_PG_arrays

call MPI_FINALIZE(ierr)
		
end program main
