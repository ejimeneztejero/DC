!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains axiliary subroutines, mostly for reading/writting data
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine fun_num_split(split_parts,ifile,num_split)

implicit none

	INTEGER :: ifile,split_parts
        CHARACTER(len=50) :: Str
        CHARACTER(len=50) :: num_split

        write(Str,*)ifile-1

        if(split_parts-1.lt.10) then
                num_split = trim(adjustl(Str))
	endif

        if(split_parts-1.ge.10.and.split_parts-1.lt.100)        then
                if(ifile-1.lt.10)                       num_split='0' // trim(adjustl(Str))
                if(ifile-1.ge.10)                       num_split= trim(adjustl(Str))
        endif

        if(split_parts-1.ge.100.and.split_parts-1.lt.1000)	then
                if(ifile-1.lt.10)                       num_split = '00' // trim(adjustl(Str))
                if(ifile-1.ge.10.and.ifile-1.lt.100)    num_split = '0' // trim(adjustl(Str))
                if(ifile-1.ge.100)                      num_split = trim(adjustl(Str))
        endif

        if(split_parts-1.ge.1000)	then
                if(ifile-1.lt.10)                       num_split = '000' // trim(adjustl(Str))
                if(ifile-1.ge.10.and.ifile-1.lt.100)    num_split = '00' // trim(adjustl(Str))
                if(ifile-1.ge.100.and.ifile-1.lt.1000)  num_split = '0' // trim(adjustl(Str))
                if(ifile-1.ge.1000)                     num_split = trim(adjustl(Str))
        endif

end subroutine fun_num_split

SUBROUTINE READ_DATA(icount,unit_new,pos_byte,ifile,nSS,input)

USE mod_parfile

IMPLICIT NONE

	INTEGER :: nh,unit_new
	INTEGER :: icount,nSS,i,j,k,kk
	INTEGER :: ifile
	INTEGER(4) :: pos_byte,pos_r

	REAL(4) :: input(nt,nSS) 	
	REAL(4), ALLOCATABLE :: sudata(:,:) 		

	CHARACTER(len=50) :: access,form
	CHARACTER(len=500) :: file_name

	nh = 60

	ALLOCATE(sudata(nt+nh,nSS))
	sudata=0.;

!	Loop trace by trace
	do j=1,nSS

		pos_r = pos_byte + (j-1)*(nh+nt)*4
		READ(unit_new+ifile, pos=pos_r)  sudata(1:nt+nh,j)	

	enddo

	do j=1,nSS

		do k=1,nt
			input(k,j)=sudata(k+nh,j)
		enddo

	enddo

	deallocate(sudata)

END SUBROUTINE READ_DATA

SUBROUTINE READ_DATA_trace(rank,icount,j,unit_new,ifile,pos_r,trace)

USE mod_parfile
USE mod_SG_arrays
USE mod_PG_arrays

IMPLICIT NONE

INTEGER :: nh,rank,unit_new
INTEGER :: icount,i,j,k,ifile
INTEGER(4) :: pos_byte,pos_r

REAL(4) :: trace(nt) 	
REAL(4), ALLOCATABLE :: sudata(:) 		

nh = 60 ! Size su header = 60*4 bytes

allocate(sudata(nt+nh))
sudata=0.;

pos_byte=pos_r+(j-1)*(nh+nt)*4

READ(unit_new+ifile,pos=pos_byte) sudata(1:nt+nh)	

do k=1,nt
	trace(k)=sudata(k+nh)
enddo

deallocate(sudata)

END SUBROUTINE READ_DATA_trace


SUBROUTINE WRITE_DATA_trace(rank,icount,is,ishot,jrec)

USE mod_parfile
USE mod_SG_arrays
USE mod_PG_arrays

IMPLICIT NONE

INTEGER :: nh,rank,unit_new
INTEGER :: icount,i,is,k,ifile,ishot,jrec
INTEGER(4) :: pos_byte,pos_r,pos_w,ifile_r,ifile_w

REAL(4) :: trace(nt) 	
REAL(4), ALLOCATABLE :: sudata(:) 		

nh = 60 ! Size su header = 60*4 bytes

allocate(sudata(nt+nh))
sudata=0.;

ifile_r=ifile_PG(icount)
pos_r = pos_byte_PG(icount) + (is-1)*(nt+nh)*4
READ(unit_PG2+ifile_r, pos=pos_r)  sudata(1:nh+nt)         !! read one trace from PG2

ifile_w=ifile_su(ishot)
pos_w = pos_byte_su(ishot) + (jrec-1)*(nt+nh)*4
WRITE(unit_DC2+ifile_w, pos=pos_w)  sudata(1:nh+nt)         !!write one trace from DC2


deallocate(sudata)

END SUBROUTINE WRITE_DATA_trace

subroutine WRITE_PG(iDC,icount,nSS,unit_PG,output)

USE mod_parfile
USE mod_input_arrays
USE mod_SG_arrays
USE mod_PG_arrays

implicit none

integer :: unit_PG
integer :: i,j,jrec,k,ishot,iDC,icount,rank
integer :: nSS,nh
integer :: ifile_r,ifile_w

INTEGER(4) :: pos_r,pos_w

real(4) :: output(nt,nSS)
REAL(4), ALLOCATABLE :: sudata(:,:),shot(:,:)

nh = 60 ! Size su header = 60*4 bytes

allocate(sudata(nh,nSS))
sudata=0.;
allocate(shot(nt,nSS))
shot=0.;

do j=1,nSS	!!specific shot

	ishot=shot_grid_PG(j,icount)	
	jrec=trace_grid_PG(j,icount)	

	ifile_r=ifile_su(ishot)
	pos_r = pos_byte_su(ishot) + (jrec-1)*(nt+nh)*4
	
	READ(unit_DC0+ifile_r, pos=pos_r)  sudata(1:nh,j)		!! read one trace from DC0

enddo

ifile_w=ifile_PG(icount)

if(iDC.eq.1)	then
do j=1,nSS

	do k=1,nt
                shot(k,j)=output(k,j) 
        enddo

enddo
endif

if(iDC.eq.2)	then
do j=1,nSS

	do k=1,nt
                shot(k,j)=output(nt-k+1,j) 
        enddo

enddo
endif

do j=1,nSS	!!specific shot

	pos_w = pos_byte_PG(icount) + (j-1)*(nt+nh)*4

	WRITE(unit_PG+ifile_w, pos=pos_w)  sudata(1:nh,j)	!! write header from the trace
	WRITE(unit_PG+ifile_w, pos=pos_w+nh*4) shot(1:nt,j)	!! write trace from point gather

enddo

deallocate(sudata)
deallocate(shot)

END SUBROUTINE WRITE_PG

subroutine write_SG0()

USE mod_parfile
USE mod_SG_arrays

implicit none
include 'mpif.h'

integer :: numtasks,rank,ierr,status(MPI_STATUS_SIZE)

integer :: nsamples,ntimes,itimes
integer :: i,j,k,nh,icount
integer :: ifile,jrec,fldr
INTEGER(4) :: pos_byte,pos_read

REAL(4), ALLOCATABLE :: sudata(:),sudata2(:,:)

call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
ierr=0;

pos_read=9	
nh = 60	! Size su header = 60*4 bytes

! Memory allocation
ALLOCATE(sudata(nh+nt))
sudata=0.;
ALLOCATE(sudata2(nh+nt,NumRec))
sudata2=0.;

nsamples=NumShots

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

if(shotID_(icount).ne.0.and.icount.le.nsamples)     then

	ifile=ifile_su(icount)

	if(reverse_streamer.eq.0)	then

		do j=1,NumRec

			pos_byte = pos_byte_su(icount) + (j-1)*(nh+nt)*4

			READ(unit0+ifile, pos=pos_byte)  sudata(1:nh+nt)
			WRITE(unit_DC0+ifile, pos=pos_byte)  sudata(1:nh+nt)

		enddo


	endif

	if(reverse_streamer.ne.0)	then

		do j=1,NumRec

			pos_byte = pos_byte_su(icount) + (j-1)*(nh+nt)*4
			READ(unit0+ifile, pos=pos_byte)  sudata2(1:nh+nt,j)

		enddo

		do j=1,NumRec

			pos_byte = pos_byte_su(icount) + (j-1)*(nh+nt)*4
			jrec=NumRec-j+1
			WRITE(unit_DC0+ifile, pos=pos_byte)  sudata2(1:nh+nt,jrec)	!! write header

		enddo

	endif

endif

enddo


deallocate(sudata)
deallocate(sudata2)


END SUBROUTINE write_SG0


subroutine WRITE_SG(iDC,icount,correction,nSS,unit_DC,SG)

USE mod_parfile
USE mod_SG_arrays

implicit none

integer :: i,j,k,nh,it,iDC,icount
integer :: ifile,unit_DC
INTEGER :: jrec,kdc,nSS
INTEGER :: correction
INTEGER(4) :: pos_byte	

real :: SG(nt,nSS)
REAL(4), ALLOCATABLE :: sudata(:,:),shot(:,:)

CHARACTER(len=500) :: file_name,Str,Str_DC

ifile=ifile_su(icount)

nh = 60	! Size su header = 60*4 bytes

ALLOCATE(sudata(nh,nSS))
ALLOCATE(shot(nt,nSS))
sudata=0.;shot=0.;

! Loop trace by trace

do j=1,nSS

	pos_byte = pos_byte_su(icount) + (j-1)*(nh+nt)*4
	READ(unit_DC0+ifile, pos=pos_byte)  sudata(1:nh,j)	!! read data

enddo

do j=1,nSS

	do k=1,nt
		shot(k,j)=SG(nt-k+1,j)
	enddo

enddo

do j=1,nSS

	pos_byte = pos_byte_su(icount) + (j-1)*(nh+nt)*4

	WRITE(unit_DC+ifile, pos=pos_byte)  sudata(1:nh,j)	!! write header
        WRITE(unit_DC+ifile, pos=pos_byte+nh*4) shot(1:nt,j)	!! real traces

enddo

deallocate(sudata)
deallocate(shot)

END SUBROUTINE WRITE_SG


subroutine SAVE_SHOTS_TXT(iDC,corr)

USE mod_parfile
USE mod_input_arrays
USE mod_SG_arrays
USE mod_PG_arrays

implicit none
include 'mpif.h'

integer :: numtasks,rank,ierr,status(MPI_STATUS_SIZE)
integer :: iDC,i,j,k,nsamples,ntimes
integer :: ifile_r,itimes,nh,icount
integer :: unit_DC,step,corr
integer(4) :: pos_r,pos_read,fldr

character(len=50) :: Str0,Str_,Str_txt,Str_DC,Str
character(len=50) :: Str_mat,Str_gnu
character(len=500) :: file_name,file_name2,file_name3
character(len=500) :: file_name4,file_name5,file_name6

integer, allocatable :: period(:)

REAL(4), ALLOCATABLE :: sudata(:,:),SG(:,:)

call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
ierr=0;

nh=60
pos_read=9

allocate(sudata(nt+nh,NumRec))
allocate(SG(nt,NumRec))
SG=0.;sudata=0.;

if(corr.eq.0)	then

	if(iDC.eq.0)unit_DC=unit_DC0
	if(iDC.eq.1)unit_DC=unit_DC1
	if(iDC.eq.2)unit_DC=unit_DC2
	Str_gnu= 'gnuplot_shot_DC'
	Str_mat= 'matlab_shot_DC'

endif

if(corr.ne.0)	then

	if(iDC.eq.1)unit_DC=unit_DC1_corr
	if(iDC.eq.2)unit_DC=unit_DC2_corr
	Str_gnu= 'gnuplot_corrected_shot_DC'
	Str_mat= 'matlab_corrected_shot_DC'

endif

write(Str_DC,*) iDC
Str_= '_'
Str_txt = '.txt'

nsamples=NumShots

ntimes=nsamples
if(numtasks.gt.1)	then

        if(nsamples.gt.numtasks)        then
                ntimes=ceiling(1.*nsamples/numtasks)
        endif

        if(nsamples.le.numtasks)        then
                ntimes=1
        endif

endif

step=floor(1.*NumShots/step_txt)
allocate(period(step))

do i=1,step

	period(i)=(i-1)*step_txt

enddo

do itimes=1,ntimes

	icount=(itimes-1)*numtasks+rank+1        !!aquí sucede la paralelización

	if(shotID_(icount).ne.0.and.icount.le.nsamples)	then

	do i=1,step

		if(icount.eq.period(i))	then

		ifile_r=ifile_su(icount)
	
		do j=1,NumRec

			pos_r=pos_byte_su(icount)+(j-1)*(nh+nt)*4
			READ(unit=unit_DC+ifile_r,pos=pos_r) sudata(1:nt+nh,j)

		enddo

	        do j=1,NumRec

		        do k=1,nt
				SG(k,j)=sudata(k+nh,j)
        		enddo

        	enddo

		write(Str,*) shotID_su(icount)

		if(save_matlab_txt.ne.0)	then

		file_name = trim(adjustl(folder_output)) // trim(adjustl(Str_mat))  
		file_name2 = trim(adjustl(file_name)) // trim(adjustl(Str_DC))
		file_name3 = trim(adjustl(file_name2)) // trim(adjustl(Str_))
		file_name4 = trim(adjustl(file_name3)) // trim(adjustl(Str)) 
		file_name5 = trim(adjustl(file_name4)) // trim(adjustl(Str_txt))


		open(12,FILE=file_name5,STATUS='unknown')
	        do k=1,nt
			write(12,'(20000(e12.5,2x))') (SG(k,j),j=1,NumRec)
        	enddo
		close(12)

		endif	!matlab

		if(save_gnuplot_txt.ne.0)	then

		file_name = trim(adjustl(folder_output)) // trim(adjustl(Str_gnu))  
		file_name2 = trim(adjustl(file_name)) // trim(adjustl(Str_DC))
		file_name3 = trim(adjustl(file_name2)) // trim(adjustl(Str_))
		file_name4 = trim(adjustl(file_name3)) // trim(adjustl(Str)) 
		file_name5 = trim(adjustl(file_name4)) // trim(adjustl(Str_txt))

		open(12,FILE=file_name5,STATUS='unknown')

		do j=1,NumRec
			do k=1,nt
				write(12,*)j,k,SG(k,j)
			enddo
		enddo
		close(12)

		endif	!gnuplot

		endif

	enddo

endif	!icount

enddo


deallocate(SG)
deallocate(sudata)
deallocate(period)

end subroutine SAVE_SHOTS_TXT


SUBROUTINE ascii_art(i)

USE mod_parfile

implicit none


integer :: i

if(i.eq.1)	then

write(*,*)"	______                                         _ 	"
write(*,*)"	|  _  \                                       | |	"
write(*,*)"	| | | |_____      ___ ____      ____ _ _ __ __| |	"
write(*,*)"	| | | / _ \ \ /\ / / '_ \ \ /\ / / _` | '__/ _` |	"
write(*,*)"	| |/ / (_) \ V  V /| | | \ V  V / (_| | | | (_| |	"
write(*,*)"	|___/ \___/ \_/\_/ |_| |_|\_/\_/ \__,_|_|  \__,_|	"                                                             
write(*,*)                                                             
write(*,*)"	 _____             _   _                   _   _		"
write(*,*)"	/  __ \           | | (_)                 | | (_)            	"
write(*,*)"	| /  \/ ___  _ __ | |_ _ _ __  _   _  __ _| |_ _  ___  _ __	"
write(*,*)"	| |    / _ \| '_ \| __| | '_ \| | | |/ _` | __| |/ _ \| '_ \	"
write(*,*)"	| \__/\ (_) | | | | |_| | | | | |_| | (_| | |_| | (_) | | | |	"
write(*,*)"	 \____/\___/|_| |_|\__|_|_| |_|\__,_|\__,_|\__|_|\___/|_| |_|	"
write(*,*)
write(*,*)                                                                      
write(*,*)"		First release (2021)			"
write(*,*)"		Author: Clara Estela Jiménez Tejero	"
write(*,*)"		email: ejimenez@icm.csic.es 		"
write(*,*)"		Barcelona Center for Subsurface Imaging "
write(*,*)"		Instituto de Ciencias Marinas (ICM-CSIC)"
write(*,*)
write(*,*)
write(*,*)"	                     |				"
write(*,*)"	                     |				"
write(*,*)"	            |        |				"
write(*,*)"	          |-|-|      |				"
write(*,*)"	            |        |				"
write(*,*)"	            | {O}    |				"
write(*,*)"	            '--|     |				"
write(*,*)"	              .|]_   |				"
write(*,*)"	        _.-=.' |     |				"	
write(*,*)"	       |    |  |]_   |				"
write(*,*)"	       |_.-='  |   __|__				"
write(*,*)"	        _.-='  |\   /|\				"
write(*,*)"	       |    |  |-'-'-'-'-.				"
write(*,*)"	       |_.-='  '========='				"
write(*,*)"	            `   |     |				"
write(*,*)"	             `. |    / \				"
write(*,*)"	               ||   /   \____.--=''''==--.._		"
write(*,*)"	               ||_.'--=='    |__  __  __  _.'	"
write(*,*)"	               ||  |    |    |\ ||  ||  || |                        ___	"	
write(*,*)"	  ____         ||__|____|____| \||__||__||_/    __________________/|   |	"
write(*,*)"	 |    |______  |===.---. .---.========''''=-._ |     |     |     / |   |	"
write(*,*)"	 |    ||     |\| |||   | |   |      '===' ||  \|_____|_____|____/__|___|	"
write(*,*)"	 |-.._||_____|_\___'---' '---'______....---===''======//=//////========|	"
write(*,*)"	 |--------------\------------------/-----------------//-//////---------/	"
write(*,*)"	 |               \                /                 // //////         /	"
write(*,*)"	 |                \______________/                 // //////         /	"
write(*,*)"	 |                                        _____===//=//////=========/	"
write(*,*)"	 |=================================================================/		"
write(*,*)"	  -----------------------------------------------------------------		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)
write(*,*)	
write(*,*)"				         ______		"
write(*,*)"				        /      \ 	"	
write(*,*)"				       /        \ 	"
write(*,*)"				       |        |	"
write(*,*)"				    )  o        o   (	"
write(*,*)"				   (    \      /    )	"
write(*,*)"				  _ \___/||||||\___/ _	"
write(*,*)"				   \____/ |||| \____/ `	"
write(*,*)"				   ,-.___/ || \__,-._	"
write(*,*)"				  /    ___/  \__	"
write(*,*)"				     _/         `--	"
write(*,*)
write(*,*)
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)"	``'-.,_,.-'``'-.,_,.='``'-.,_,.-'``'-.,_,.='````'-.,_,.-'``'-.,_,.='``'-.		"
write(*,*)
write(*,*)


endif

if(i.eq.2)	then

write(*,*)
write(*,*)
write(*,*)'Fight your bug'
write(*,*)'                                |     |'
write(*,*)'                                \\_V_//'
write(*,*)'                                \/=|=\/'
write(*,*)'                                 [=v=]'
write(*,*)'                               __\___/_____'
write(*,*)'                              /..[  _____  ]'
write(*,*)'                             /_  [ [  M /] ]'
write(*,*)'                            /../.[ [ M /@] ]'
write(*,*)'                           <-->[_[ [M /@/] ]'
write(*,*)'                          /../ [.[ [ /@/ ] ]'
write(*,*)'     _________________]\ /__/  [_[ [/@/ C] ]'
write(*,*)'    <_________________>>0---]  [=\ \@/ C / /'
write(*,*)'       ___      ___   ]/000o   /__\ \ C / /'
write(*,*)'          \    /              /....\ \_/ /'
write(*,*)'       ....\||/....           [___/=\___/'
write(*,*)'      .    .  .    .          [...] [...]'
write(*,*)'     .      ..      .         [___/ \___]'
write(*,*)'     .    0 .. 0    .         <---> <--->'
write(*,*)'  /\/\.    .  .    ./\/\      [..]   [..]'
write(*,*)' / / / .../|  |\... \ \ \    _[__]   [__]_'
write(*,*)'/ / /       \/       \ \ \  [____>   <____]'
write(*,*)
write(*,*)

endif

if(i.eq.3)	then

	write(*,*)
	write(*,*)
	write(*,*)
	write(*,*)
	write(*,*)"	 _______  ___   __    _  ___   _______  __   __  _______  ______  	"
	write(*,*)"	|       ||   | |  |  | ||   | |       ||  | |  ||       ||      | 	"
	write(*,*)"	|    ___||   | |   |_| ||   | |  _____||  |_|  ||    ___||  _    |	"
	write(*,*)"	|   |___ |   | |       ||   | | |_____ |       ||   |___ | | |   |	"
	write(*,*)"	|    ___||   | |  _    ||   | |_____  ||       ||    ___|| |_|   |	"
	write(*,*)"	|   |    |   | | | |   ||   |  _____| ||   _   ||   |___ |       |	"
	write(*,*)"	|___|    |___| |_|  |__||___| |_______||__| |__||_______||______|	"
	write(*,*)
	write(*,*)

endif

END SUBROUTINE ascii_art

SUBROUTINE open_su_files()

USE mod_parfile
USE mod_SG_arrays

IMPLICIT NONE

INTEGER :: ifile

CHARACTER(len=50) :: access,form,num_split
CHARACTER(len=500) :: file_name

access = 'stream'
form = 'unformatted'

do ifile=1,split_parts

       	call fun_num_split(split_parts,ifile,num_split)

	file_name = trim(adjustl(su_file0)) // trim(adjustl(num_split))
	file_name = trim(adjustl(folder_input)) // trim(adjustl(file_name))

	if(endianness_machine.eq.0)	then

	if(endianness_data.eq.0)open(unit0+ifile,FILE=file_name,ACCESS=access,FORM=form,STATUS='old')
	if(endianness_data.eq.1)open(unit0+ifile,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='old')

	endif

	if(endianness_machine.eq.1)	then

	if(endianness_data.eq.0)open(unit0+ifile,FILE=file_name,ACCESS=access,FORM=form,CONVERT='LITTLE_ENDIAN',STATUS='old')
	if(endianness_data.eq.1)open(unit0+ifile,FILE=file_name,ACCESS=access,FORM=form,STATUS='old')

	endif

	INQUIRE(FILE=file_name, SIZE=sizeof(ifile) )

	file_name = trim(adjustl(su_file_DC0)) // trim(adjustl(num_split))
        file_name = trim(adjustl(folder_output)) // trim(adjustl(file_name))
        open(unit_DC0+ifile,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='unknown')

enddo

if(DC.lt.0.or.DC.ge.1)	then

do ifile=1,split_parts

       	call fun_num_split(split_parts,ifile,num_split)

        file_name = trim(adjustl(su_file_DC1)) // trim(adjustl(num_split))
        file_name = trim(adjustl(folder_output)) // trim(adjustl(file_name))
        open(unit_DC1+ifile,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='unknown')

enddo

endif

if(DC.lt.0.or.DC.ge.2)	then

do ifile=1,split_parts

       	call fun_num_split(split_parts,ifile,num_split)

        file_name = trim(adjustl(su_file_PG1)) // trim(adjustl(num_split))
        file_name = trim(adjustl(folder_output)) // trim(adjustl(file_name))
        open(unit_PG1+ifile,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='unknown')

        file_name = trim(adjustl(su_file_PG2)) // trim(adjustl(num_split))
        file_name = trim(adjustl(folder_output)) // trim(adjustl(file_name))
        open(unit_PG2+ifile,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='unknown')

        file_name = trim(adjustl(su_file_DC2)) // trim(adjustl(num_split))
        file_name = trim(adjustl(folder_output)) // trim(adjustl(file_name))
        open(unit_DC2+ifile,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='unknown')

enddo

endif

if(phase_correction.eq.1)	then

do ifile=1,split_parts

       	call fun_num_split(split_parts,ifile,num_split)

	if(DC.eq.1.or.DC.lt.0)	then

        file_name = trim(adjustl(su_file_DC1_corr)) // trim(adjustl(num_split))
        file_name = trim(adjustl(folder_output)) // trim(adjustl(file_name))
        open(unit_DC1_corr+ifile,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='unknown')

	endif

	if(DC.eq.2.or.DC.lt.0)	then

        file_name = trim(adjustl(su_file_DC2_corr)) // trim(adjustl(num_split))
        file_name = trim(adjustl(folder_output)) // trim(adjustl(file_name))
        open(unit_DC2_corr+ifile,FILE=file_name,ACCESS=access,FORM=form,CONVERT='BIG_ENDIAN',STATUS='unknown')

	endif

enddo

endif

!maxbytes=maxval(sizeof(:))

END SUBROUTINE open_su_files

subroutine close_su_files()
USE mod_parfile

implicit none

INTEGER :: k

do k=1,split_parts
        close(unit_DC0+k)
        close(unit_DC1+k)
        close(unit_DC2+k)
        close(unit_PG1+k)
        close(unit_PG2+k)
enddo

if(phase_correction.eq.1)	then

	do k=1,split_parts
	        if(DC.lt.0.or.DC.eq.1)close(unit_DC1_corr+k)
	        if(DC.lt.0.or.DC.eq.2)close(unit_DC2_corr+k)
	enddo

endif

end subroutine close_su_files

subroutine WRITE_CORR(iDC)

USE mod_parfile
USE mod_SG_arrays

implicit none
include 'mpif.h'

integer :: numtasks,rank,ierr,status(MPI_STATUS_SIZE)

integer :: i,j,k,nh,it,iDC,ishot,k1,i1,i2
integer :: ifile,FA_error
integer :: itimes,ntimes,nsamples
integer :: unit_read,unit_write
INTEGER(4) :: pos_r, pos_w, pos_byte

real    :: FA_st,FA_obs,FA_dc
real    :: t_1,t_2,depth_st,depth_shot
REAL(4), ALLOCATABLE :: sudata(:,:),input(:,:),shot(:,:)

CHARACTER(len=500) :: file_new

nh = 60 ! Size su header = 60*4 bytes

ALLOCATE(sudata(nh+nt,NumRec))
ALLOCATE(input(nt,NumRec))
ALLOCATE(shot(nt,NumRec))

call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)
ierr=0;

if(iDC.eq.1)	then
	unit_read=unit_DC1
	unit_write=unit_DC1_corr
endif

if(iDC.eq.2)	then
	unit_read=unit_DC2
	unit_write=unit_DC2_corr
endif


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     MAIN LOOP DC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

nsamples=NumShots

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

ishot=(itimes-1)*numtasks+rank+1        !!aquí sucede la paralelización

i1=(itimes-1)*numtasks+1
i2=itimes*numtasks
if(i2.gt.nsamples)i2=nsamples

if(ishot.le.nsamples)     then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Calculo error de fase
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	FA_error=1
	ifile=ifile_su(ishot)
        pos_byte=pos_byte_su(ishot)

	sudata=0.;shot=0.;

        depth_st=abs(pos_bat(ishot)-streamer_depth)
	depth_shot=abs(pos_bat(ishot)-shot_depth)

        FA_st=sqrt((near_offset/2)**2.+depth_st**2.)+ &
        sqrt((near_offset/2)**2.+depth_shot**2.)
        FA_st=FA_st*(1./water_velocity) !!real

        FA_obs = sqrt(near_offset**2.+depth_shot**2.)*(1./water_velocity) !!real
        FA_dc= near_offset/land_velocity

        if(iDC.eq.1)     then
	        t_1 =  FA_obs
        	t_2 =  FA_st-depth_st/water_velocity            !!real
         endif

         if(iDC.eq.2)     then
         	t_1 = FA_dc                             !!real
         	!!t_2 = FA_obs-depth_shot/water_velocity        !!real
         	t_2 = FA_st-(depth_shot+depth_st)/water_velocity        !!real
         endif

	FA_error=1+ceiling(abs(t_1-t_2)/dt)  !!integer

!	if(rank.eq.0)   then
!		write(123,*)"ishot,FA_error, FA_error(s): ",i,FA_error(ishot),abs(t_1-t_2)
!		write(456,*)"ishot,DC,depth,depth,FA_st,t_1,t_2: ",i,DC,depth_st,depth_shot,FA_st,t_1,t_2
!	endif


! Loop trace by trace

        do j=1,NumRec
        	pos_r = pos_byte + (j-1)*(nh+nt)*4
                READ(unit_read+ifile, pos=pos_r)  sudata(1:nt+nh,j)
        enddo

	do j=1,NumRec
                do k=1,nt
                        input(k,j)=sudata(k+nh,j)
                enddo
        enddo

	k1=FA_error
	do j=1,NumRec
		do k=1,k1
        	        shot(k,j)=input(1,j)
        	enddo
        	do k=k1,nt
        	        shot(k,j)=input(k-k1+1,j)
        	enddo
	enddo

	do j=1,NumRec
        	pos_w = pos_byte + (j-1)*(nh+nt)*4
        	WRITE(unit_write+ifile, pos=pos_w)  sudata(1:nh,j)	!! write header
		WRITE(unit_write+ifile, pos=pos_w+nh*4) shot(1:nt,j)    !! real traces
	enddo


endif	!ishot
enddo	!itimes

deallocate(sudata)
deallocate(input)
deallocate(shot)

END SUBROUTINE WRITE_CORR
