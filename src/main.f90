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

integer :: iDC,i,output,corr
integer :: numtasks,rank,ierr,status(MPI_STATUS_SIZE)
!!INTEGER ::  REQUEST, IERROR

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
		write(*,*)'- WITH SPHERICAL DIVERGENCE APPLIED (TO KEEP AMPLITUDES ALONG OFFSET)'
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

if(DC.eq.0.or.DC.eq.1)	then

if(rank.eq.0)write(*,*)
if(rank.eq.0)write(*,*)'*******************************************'
if(rank.eq.0)write(*,*)'PREPARING AND WRITING INPUT DATA, STEP 0'
if(rank.eq.0)write(*,*)'*******************************************'

	iDC=0
	call WRITE_SG0()
	call MPI_barrier(MPI_COMM_WORLD,ierr)

	if(save_txt.ne.0) then
		corr=0
		if(rank.eq.0)write(*,*)'SAVING TXT DATA'
		call SAVE_SHOTS_TXT(iDC,corr)

	endif

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
	call DC_function(iDC,&
	NumShots,NumRec,&
	datum1_rec,datum2)

        if(save_txt.eq.1) then
		corr=0
                if(rank.eq.0)write(*,*)'SAVING TXT DATA'
                call SAVE_SHOTS_TXT(iDC,corr)
        endif

	endif!DC=1

endif!DC=0,DC=1

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

	call DC_function(iDC,&
	NumMax_PG,NumMaxShots_PG,&
	datum1_shot,datum2)

if(rank.eq.0)write(*,*)
if(rank.eq.0)write(*,*)'{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{'
if(rank.eq.0)write(*,*)
if(rank.eq.0)write(*,*)'CONVERSION OF POINT GATHER TO SHOT GATHERS ...'
if(rank.eq.0)write(*,*)
if(rank.eq.0)write(*,*)'}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}'

	call SG_function()

	call MPI_barrier(MPI_COMM_WORLD,ierr)

        if(save_txt.eq.1) then
		corr=0
                if(rank.eq.0)write(*,*)'SAVING TXT DATA'
                call SAVE_SHOTS_TXT(iDC,corr)
        endif

endif	!DC=2

if(phase_correction.eq.1)	then

	if(DC.eq.1.or.DC.eq.2)	then

		iDC=DC
	        if(rank.eq.0)write(*,*)
	        if(rank.eq.0)write(*,*)'PHASE CORRECTION: ',iDC
		call WRITE_CORR(iDC)

		if(save_txt.eq.1) then
			corr=1
        		if(rank.eq.0)write(*,*)'SAVING TXT DATA'
        		call SAVE_SHOTS_TXT(iDC,corr)
		endif

	endif

	if(DC.lt.0)	then
	
		do iDC=2,2

	        if(rank.eq.0)write(*,*)
	        if(rank.eq.0)write(*,*)'PHASE CORRECTION: ',iDC
		call WRITE_CORR(iDC)

		if(save_txt.ne.0) then
			corr=1
        		if(rank.eq.0)write(*,*)'SAVING TXT DATA'
        		call SAVE_SHOTS_TXT(iDC,corr)
		endif

		enddo

	endif

endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! 			END DOWNWARD CONTINUATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)   then

        end   = MPI_Wtime()
        write(*,*)
	write(*,*) 'The whole process took',nint((end-start)/60),'minutes'

!	file_name = "cp" // " " // trim(adjustl(par_file))
!	command = trim(adjustl(file_name)) // " " // trim(adjustl(folder_output))
!	call system(command)

!	file_name= "cp job*"
!	command = trim(adjustl(file_name)) // " " // trim(adjustl(folder_output))
!	call system(command)

	call ascii_art(3)

endif

call close_su_files()

call deallocate_SG_arrays
call deallocate_input_arrays
if(DC.eq.2)call deallocate_PG_arrays

call MPI_FINALIZE(ierr)
		
end program main
