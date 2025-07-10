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

character(len=100) :: DC_in
character(len=500) :: file_name,command

DOUBLE PRECISION :: start1,end1
DOUBLE PRECISION :: start2,end2
DOUBLE PRECISION :: start0,end0

logical	:: exists

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call MPI_INIT(ierr)
call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
ierr=0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	ASCII ART WELLCOME MESSAGE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(rank.eq.0)	then

		start0=MPI_Wtime()
		call ascii_art(1)

	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Lectura de parametros Par_file
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        if(rank.eq.0)write(*,*)'*********************'
        if(rank.eq.0)write(*,*)'READ INPUT PARAMETERS'
        if(rank.eq.0)write(*,*)'*********************'

        call read_parfile(rank)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	WRITE STREAMER DATA TO DC IN OUTPUT FOLDER
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Transformacion datos SGY en SU, y particiones de menos de 2 GB

	DC_in=trim(adjustl(folder_input)) // "DC_in.dat"

	if(DC.eq.0.or.DC.eq.3)	then
	
		if(rank.eq.0)call Data_in()
	
	endif

	call MPI_barrier(MPI_COMM_WORLD,ierr)

	open(unit=10, file=trim(adjustl(DC_in)), status='old', action='read')
	read(10,'(A)') su_file0
	read(10,*) split_parts
	close(10)

	if(rank.eq.0)	then
		print *, 'SU files common name: ', trim(su_file0)
		print *, 'Number of files: ', split_parts
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	INPUTS AND DATA PREPARATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	call inputs()

	if(DC.eq.0.or.DC.eq.1.or.DC.eq.3)	then

		iDC=0
		call WRITE_SG0()
	
		if(save_txt.ne.0) then
			if(rank.eq.0) write(*,*)'SAVING TXT DATA'
			call SAVE_SHOTS_TXT(iDC)
		endif
	
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!! START DC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(DC.eq.1.or.DC.eq.3)	then

		call MPI_barrier(MPI_COMM_WORLD,ierr)

		if(rank.eq.0)start1=MPI_Wtime()

		iDC=1

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

		call DC_function(iDC,NumShots,NumRec,datum1_rec,datum2)

		if(rank.eq.0)   then
        		end1   = MPI_Wtime()
		endif

	        if(save_txt.eq.1) then
                	if(rank.eq.0)write(*,*)'SAVING TXT DATA'
                	call SAVE_SHOTS_TXT(iDC)
        	endif

	endif!DC=1

	if(DC.eq.2.or.DC.eq.3)        then

		call MPI_barrier(MPI_COMM_WORLD,ierr)

		if(rank.eq.0)start2=MPI_Wtime()

		iDC=2

		if(rank.eq.0)write(*,*)
		if(rank.eq.0)write(*,*)'{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{{'
		if(rank.eq.0)write(*,*)
		if(rank.eq.0)write(*,*)'CONVERSION OF SHOT GATHERS TO POINT GATHERS ...'
		if(rank.eq.0)write(*,*)
		if(rank.eq.0)write(*,*)'}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}}'

		call PG_function()

		call MPI_barrier(MPI_COMM_WORLD,ierr)

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

		if(rank.eq.0)   then
			end2   = MPI_Wtime()
		endif

	        if(save_txt.eq.1) then
	                if(rank.eq.0)write(*,*)'SAVING TXT DATA'
	                call SAVE_SHOTS_TXT(iDC)
	        endif

	endif	!DC=2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!			END DOWNWARD CONTINUATION
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call MPI_barrier(MPI_COMM_WORLD,ierr)

if(rank.eq.0)   then

        end0   = MPI_Wtime()

	write(*,*)
	if(DC.eq.1.or.DC.eq.3)write(*,*) 'DC1 took ',nint((end1-start1)/60),'minutes'
	if(DC.eq.2.or.DC.eq.3)write(*,*) 'DC2 took ',nint((end2-start2)/60),'minutes'
	write(*,*) 'The whole process took',nint((end0-start0)/60),'minutes'

	call ascii_art(3)

endif

call close_su_files()

call deallocate_SG_arrays
call deallocate_input_arrays
if(DC.eq.2)call deallocate_PG_arrays

if(rank.eq.0)	then
	if(DC.eq.0.or.DC.eq.1.or.DC.eq.3)call Data_out(0,su_file_DC0)
	if(DC.eq.1.or.DC.eq.3)call Data_out(1,su_file_DC1)
	if(DC.eq.2.or.DC.eq.3)call Data_out(2,su_file_DC2)
endif

call MPI_FINALIZE(ierr)
		
end program main
