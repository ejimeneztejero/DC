!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains the modules:
!!	(mod_parfile, mod_PG_arrays, mod_SG_arrays, mod_input_arrays)
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_parfile


implicit none


  ! Control file variables

  integer, parameter :: unit0=100,unit_DC0=200
  integer, parameter :: unit_DC1=300,unit_DC2=400
  integer, parameter :: unit_PG1=500,unit_PG2=600

  integer, parameter :: LITTLE_ENDIAN = 0, BIG_ENDIAN = 1
 
  INTEGER, parameter :: size_su_header = 60
  INTEGER(4), parameter :: byte_fldr = 9
  INTEGER(4), parameter :: byte_offset = 37
  INTEGER(4), parameter :: byte_scalco = 71
  INTEGER(4), parameter :: byte_sx = 73
  INTEGER(4), parameter :: byte_sy = 77
  
  INTEGER(4) :: byte_shotnumber
  integer :: data_type
  integer :: DC
  integer :: reverse_streamer,reg_grid
  integer :: nt,NumRec,NumShots,print_bat,print_geom
  integer :: shot_init,shot_fin
  integer :: far_offset_grid

  integer :: endianness_data,endianness_machine
  integer :: sx_sy_header,offset_header
  integer :: TWT_option,split_parts
  integer :: save_txt,step_txt
  integer :: save_gmt,save_matlab

  integer :: NumMax_PG,NumMaxShots_PG
  integer :: maxbytes
  
  real :: dmodel,dt,near_offset,offset_unit,far_offset
  real :: shot_depth,streamer_depth
  real :: added_space_model_X,added_space_model_Y,water_velocity
  real :: dshots,drec,time

  real :: f1,f2

  character(len=500) :: folder_input,folder_output,folder_data,folder_data_parts
  character(len=500) :: folder_gmt,folder_matlab

  character(len=500) :: data_file
  character(len=500) :: sgy_file,su_file
  character(len=500) :: su_file0
  character(len=500) :: su_file_DC0
  character(len=500) :: su_file_DC1,su_file_DC2
  character(len=500) :: su_file_PG1,su_file_PG2
  character(len=500) :: temp_DC1,temp_DC2
  character(len=500) :: temp_PG1,temp_PG2
  character(len=500) :: nav_file,vp_file
  character(len=500) :: par_file,dc_str

  contains

  subroutine read_parfile(rank)


  implicit none
!  use ISO_FORTRAN_ENV, only: IACHAR

  include 'mpif.h'

  integer :: numtasks,rank,ierr,errcode,status(MPI_STATUS_SIZE)

  ! Input related variables
  character(len=200) :: buffer,label
  integer :: pos,icount,i,ifile,interval,DC_b
  integer, parameter :: fh = 10
  integer :: ios = 0
  integer :: line = 0
  integer :: NumMax_PG_,NumMaxShots_PG_
  character(len=500) :: file_name,command0,command
  character(len=50) :: Str,access,form,num_split
  character(len=100) :: arg1, arg2

  logical :: su_exist,nav_exist,vp_exist
  logical :: input_exist, output_exist

  call MPI_COMM_SIZE(MPI_COMM_WORLD,numtasks,ierr)
  call MPI_COMM_RANK(MPI_COMM_WORLD,rank,ierr)

  ierr=0;
  DC=-10;

  access = 'STREAM'
  form = 'UNFORMATTED'

  icount = iargc()  ! The number does not includes the executable name, so if user passed one argument, 1 is returned.

!  if (icount < 2) then
  if (icount < 1) then
	  if(rank.eq.0) then
	    write(*,*)
	    write(*,*)
	    write(*,*) 'ERROR: parfile argument is required in execution line'
	    write(*,*)
	    write(*,*)
	  endif

	  if(rank.eq.0)call ascii_art(2)
	  call MPI_barrier(MPI_COMM_WORLD,ierr)
	  call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	  stop

  endif

  if (icount.ge.1)  call getarg(1, par_file)
  if (icount.eq.2) then
	call getarg(2, dc_str)
        read(dc_str, *) DC
  endif

  if(rank.eq.0)write(*,*)'name par_file: ',trim(adjustl(par_file))
  file_name = trim(adjustl(par_file))
  INQUIRE(file=file_name,EXIST=su_exist)
  if(.NOT. su_exist)     then
	if(rank.eq.0)write(*,*)'ERROR: Par_file named: ', trim(adjustl(par_file)),' does not exist'
       	if(rank.eq.0)call ascii_art(2)
	call MPI_barrier(MPI_COMM_WORLD,ierr)
	call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
  endif

  data_type = 0

  data_file = 'null'
  sgy_file = 'null'
  su_file = 'null'
  
  su_file0 = 'null'

  su_file_DC0 = 'su_DC0'
  su_file_DC1 = 'su_DC1'
  su_file_DC2 = 'su_DC2'

  su_file_PG1 = 'su_PG1'
  su_file_PG2 = 'su_PG2'

  temp_DC1 = 'temp_DC1_'
  temp_DC2 = 'temp_DC2_'
  temp_PG1 = 'temp_PG1_'
  temp_PG2 = 'temp_PG2_'

  nav_file = 'null'
  vp_file = 'null'

  print_bat=1
  print_geom=1
  reg_grid=0
  byte_shotnumber= byte_fldr
  sx_sy_header=0;offset_header=0;offset_unit=1
  TWT_option=0;
  
  endianness_data=1; 
  endianness_machine=0; 

  near_offset=0;reverse_streamer=0;
  added_space_model_X=0.
  added_space_model_Y=0.
  far_offset=0.;far_offset_grid=0;

  water_velocity=1500.
  drec=0;dmodel=0;dt=0;dshots=0;
  shot_depth=0;streamer_depth=0;
  NumRec=0;nt=0;split_parts=1;
  step_txt=1;
  save_gmt=0;save_matlab=0;
  DC_b=-10;

  open(fh, file=par_file)

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected. It is positive if an error was
  ! detected.  ios is zero otherwise.

  do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

	 buffer = adjustl(buffer)  ! Remove leading spaces
	 if (len_trim(buffer) == 0) cycle  ! <-- Skip empty lines


        ! Verificar si la línea es un comentario
        if (buffer(1:1) .ne. '#') then

        ! Find the first instance of whitespace.  Split label and data.
        pos = scan(buffer,' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)

        case ('DC:')
           read(buffer, *, iostat=ios) DC_b
        case ('reg_grid:')
           read(buffer, *, iostat=ios) reg_grid
        case ('byte_shotnumber:')
           read(buffer, *, iostat=ios) byte_shotnumber
        case ('sx_sy_header:')
           read(buffer, *, iostat=ios) sx_sy_header
        case ('offset_header:')
           read(buffer, *, iostat=ios) offset_header
        case ('offset_unit:')
           read(buffer, *, iostat=ios) offset_unit
        case ('TWT_option:')
           read(buffer, *, iostat=ios) TWT_option
        case ('endianness_data:')
           read(buffer, *, iostat=ios) endianness_data
        case ('endianness_machine:')
           read(buffer, *, iostat=ios) endianness_machine
        case ('data_format:')
           read(buffer, *, iostat=ios) data_type
        case ('data_file:')
           read(buffer, *, iostat=ios) data_file
        case ('nav_file:')
           read(buffer, *, iostat=ios) nav_file
        case ('vp_file:')
           read(buffer, *, iostat=ios) vp_file
        case ('input_folder:')
           read(buffer, *, iostat=ios) folder_input
        case ('output_folder:')
           read(buffer, *, iostat=ios) folder_output
        case ('folder:')
           read(buffer, *, iostat=ios) folder_output
        case ('dt:')
           read(buffer, *, iostat=ios) dt
        case ('drec:')
           read(buffer, *, iostat=ios) drec
        case ('dmodel:')
           read(buffer, *, iostat=ios) dmodel
        case ('dshots:')
           read(buffer, *, iostat=ios) dshots
        case ('NumRec:')
           read(buffer, *, iostat=ios) NumRec
        case ('shot_init:')
           read(buffer, *, iostat=ios) shot_init
        case ('shot_fin:')
           read(buffer, *, iostat=ios) shot_fin
        case ('ns:')
           read(buffer, *, iostat=ios) nt
        case ('nt:')
           read(buffer, *, iostat=ios) nt
        case ('near_offset:')
           read(buffer, *, iostat=ios) near_offset
        case ('shot_depth:')
           read(buffer, *, iostat=ios) shot_depth
        case ('streamer_depth:')
           read(buffer, *, iostat=ios) streamer_depth
        case ('reverse_streamer:')
           read(buffer, *, iostat=ios) reverse_streamer
        case ('water_velocity:')
           read(buffer, *, iostat=ios) water_velocity
        case ('save_gmt:')
           read(buffer, *, iostat=ios) save_gmt
        case ('save_matlab:')
           read(buffer, *, iostat=ios) save_matlab
        case ('shot_step_txt:')
           read(buffer, *, iostat=ios) step_txt
        case default
           if(rank.eq.0)print *, 'WARNING in file ',trim(adjustl(par_file)),': skipping invalid label at line', line

        end select

        end if  ! Fin del if para verificar comentarios

     end if
  end do

 close(fh)

file_name = trim(adjustl(folder_input)) 
INQUIRE(FILE=file_name, EXIST=input_exist)
if(.NOT. input_exist)	then
	if(rank.eq.0)write(*,*)'ERROR: input folder not found'
        if(rank.eq.0)call ascii_art(2)
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

file_name = trim(adjustl(folder_output)) 
INQUIRE(FILE=file_name, EXIST=output_exist)
if(.NOT. output_exist)	then
	if(rank.eq.0)write(*,*)'WARNING: output folder not found, creating output folder named as: output'
	folder_output="output"
endif

folder_input = trim(adjustl(folder_input)) // '/'
folder_output = trim(adjustl(folder_output)) // '/'
folder_data = trim(adjustl(folder_output)) // 'DATA/'
folder_data_parts = trim(adjustl(folder_data)) // 'PARTS/'

if(rank.eq.0)	then
	write(*,*)
	write(*,*)"CREATING DATA OUTPUT FOLDERS"
endif
command="mkdir " // trim(adjustl(folder_output))
if(rank.eq.0)call system(command)
command="mkdir " // trim(adjustl(folder_data))
if(rank.eq.0)call system(command)
command="mkdir " // trim(adjustl(folder_data_parts))
if(rank.eq.0)call system(command)

time=(nt-1)*dt

if(dmodel.eq.0)dmodel=drec

if(dmodel.gt.drec)	then
	if(rank.eq.0)write(*,*)'dmodel cannot be greather than drec, instead, it is set to drec'
	dmodel=drec
endif

NumShots=shot_fin-shot_init+1

far_offset=(NumRec-1)*drec+near_offset
far_offset_grid=1+ceiling(far_offset/dmodel)

!NumMaxShots_PG_=1+ceiling(far_offset/dshots)	!!numero maximo de shots por PG
!NumMaxShots_PG=NumMaxShots_PG_+ceiling(NumMaxShots_PG_/5.)

!NumMax_PG_=1+ceiling((NumShots-1)*dshots+far_offset)/dmodel
!NumMax_PG=NumMax_PG_+ceiling(NumMax_PG_/5.)

if(added_space_model_X.eq.0)added_space_model_X=20.*dmodel
if(added_space_model_Y.eq.0)added_space_model_Y=20.*dmodel

if(DC.lt.0) then
	if(DC_b.ge.0)DC=DC_b
endif

if(DC.ne.0.and.DC.ne.1.and.DC.ne.2.and.DC.ne.3) then
if(rank.eq.0)write(*,*)'DC value must be 0, 1, 2 or 3: '
if(rank.eq.0)write(*,*)'Introduce flag 0, 1, 2 or 3 at the command line, or at the parfile'
        if(rank.eq.0)call ascii_art(2)
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(DC.eq.0)	then
if(data_type.ne.1.and.data_type.ne.2)	then
	if(rank.eq.0)write(*,*)'WARNING: data_type should be set to 1 (SGY file) or 2 (SU file)'
!        if(rank.eq.0)call ascii_art(2)
!        call MPI_barrier(MPI_COMM_WORLD,ierr)
!        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
!	stop
endif
endif

if(sx_sy_header.ne.0.and.sx_sy_header.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: sx_sy_header should be set to 0 or 1 in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(offset_header.ne.0.and.offset_header.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: offset_header should be set to 0 or 1 in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(TWT_option.ne.0.and.TWT_option.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: TWT_option should be set to 0 or 1 in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(endianness_data.ne.0.and.endianness_data.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: endianness_data should be set to 0 (little endian)&
	 or 1 (big endian) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(endianness_machine.ne.0.and.endianness_machine.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: endianness_machine should be set to 0 (little endian)&
	 or 1 (big endian) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(reverse_streamer.ne.0.and.reverse_streamer.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: reverse_streamer should be set to 0 (first channel the closest to shot) &
	or 1 (last channel the closes to shot) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(offset_header.eq.0.and.near_offset.eq.0)	then
	if(rank.eq.0)write(*,*)'ERROR: If offset_header: 0, near_offset should be set &
	to a certain value (in meters) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

file_name = trim(adjustl(folder_input)) // trim(adjustl(nav_file))
INQUIRE(FILE=file_name, EXIST=nav_exist)
if(.NOT. nav_exist)	then
	if(rank.eq.0)write(*,*)'ERROR: nav_file not found'
        if(rank.eq.0)call ascii_art(2)
        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

file_name = trim(adjustl(folder_input)) // trim(adjustl(vp_file))
INQUIRE(FILE=file_name, EXIST=vp_exist)
if(.NOT. vp_exist)	then
if(water_velocity.le.1400)	then
	if(rank.eq.0)write(*,*) 'ERROR: vp_file not found in your your input folder &
	and water_velocity is set to a wrong value in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif
endif

if(dshots.eq.0) then
!if(dshots.eq.0.and.reg_grid.eq.1) then
	if(rank.eq.0)write(*,*) 'ERROR: Please give an average value to dshots (meters) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(dt.eq.0) then
	if(rank.eq.0)write(*,*) 'ERROR: Please give a value to dt (seconds) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(nt.eq.0) then
	if(rank.eq.0)write(*,*)'ERROR: Please give a value to nt in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif
if(shot_depth.eq.0) then
	if(rank.eq.0)write(*,*)'ERROR: Please give a value to shot_depth (meters) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif
if(streamer_depth.eq.0) then
	if(rank.eq.0)write(*,*)'ERROR: Please give a value to streamer_depth (meters) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif
if(NumRec.eq.0) then
	if(rank.eq.0)write(*,*) 'ERROR: Please give a value to NumRec in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif
if(drec.eq.0) then
	if(rank.eq.0)write(*,*) 'ERROR: Please give a value to drec in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

if(shot_init.le.0) then
	if(rank.eq.0)write(*,*) 'ERROR: Please give a correct value to shot_init in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

        call MPI_barrier(MPI_COMM_WORLD,ierr)
        call MPI_Abort(MPI_COMM_WORLD, errcode, ierr)
	stop
endif

maxbytes=1900000000

if(save_matlab.ne.0.or.save_gmt.ne.0)save_txt=1


if(rank.eq.0)	then

write(*,*)
write(*,*)'*******************************************'
write(*,*)'RUNNING PROGRAM FOR NEXT SET OF PARAMETERS:'
write(*,*)'*******************************************'

if(DC.ge.0)	then
	write(*,*)'Case DC:', DC
	write(*,*)'-------------------------------------'
endif
if(sx_sy_header.eq.1)	then	
	write(*,*)'sx_sy_header: ',sx_sy_header
	write(*,*)'The .sgy file should contain shotgather geometry (UTM) in headers (variables: sx, sy)'
	write(*,*)'TWT_option: ',TWT_option
	if(TWT_option.eq.1)write(*,*)'nav_file should contain: ShotID	Z(m)'
	if(TWT_option.eq.1)write(*,*)'nav_file should contain: ShotID	TWT(s)'
endif

if(sx_sy_header.eq.0)	then
	write(*,*)'sx_sy_header: ',sx_sy_header
	write(*,*)'No shot gathers geometry in headers'
	write(*,*)'TWT_option: ',TWT_option
	if(reg_grid.eq.1.and.TWT_option.eq.0)write(*,*)'nav_file minimum should contain: ShotID Z(m)'
	if(reg_grid.eq.1.and.TWT_option.eq.1)write(*,*)'nav_file minimum should contain: ShotID TWT(s)'
	if(reg_grid.eq.0)	then
		if(TWT_option.eq.0)write(*,*)'nav_file should contain: ShotID	Xshot(UTM)	Yshot(UTM)	Z(m)'
		if(TWT_option.eq.1)write(*,*)'nav_file should contain: ShotID	Xshot(UTM)	Yshot(UTM)	TWT(s)'
	endif
endif

write(*,*)'offset_header: ',offset_header

if(offset_header.eq.1) then
	write(*,*)
	write(*,*)'The .sgy file should contain offset information in headers (variable offset)'
	write(*,*)'offset_unit: ',offset_unit
else if(offset_header.eq.0) then
	write(*,*)
	write(*,*)'Streamer is considered rigid -fix distance between channels (dmodel)- &
	and fix value for near offset (near_offset)'
	write(*,*)'near_offset: ',near_offset
endif

write(*,*)'Input folder: ',adjustl(trim(folder_input))
write(*,*)'Output folder: ',adjustl(trim(folder_output))
write(*,*)'name of data file (SU): ',adjustl(trim(su_file))
write(*,*)'split_parts: ',split_parts
write(*,*)'name of navigation file: ',adjustl(trim(nav_file))
if(reg_grid.eq.1)write(*,*)'REGULAR GRID FOR SHOTS AND RECEIVERS'
if(reg_grid.eq.1)write(*,*)'dshots: ',dshots
write(*,*)'drec: ',drec
write(*,*)'NumRec: ',NumRec
write(*,*)'NumShots: ',NumShots
!write(*,*)'Num of Point Gathers: approx ',NumMax_PG_
write(*,*)'shot_init: ',shot_init
write(*,*)'shot_fin: ',shot_fin
write(*,*)'dt: ',dt
write(*,*)'nt: ',nt
write(*,*)'dmodel: ',dmodel
write(*,*)'shot_depth: ',shot_depth
write(*,*)'streamer_depth: ',streamer_depth
if(reverse_streamer.ne.0)write(*,*)'reverse_streamer: ',reverse_streamer
if(vp_exist)write(*,*)'vp_file: ',adjustl(trim(vp_file))
if(.NOT. vp_exist)write(*,*)'water_velocity constant: ',water_velocity
!write(*,*)'NumMaxShots_PG,NumMax_PG: ',NumMaxShots_PG,NumMax_PG

endif

if(save_matlab.eq.1)    then
                if(rank.eq.0)write(*,*)
                if(rank.eq.0)write(*,*)"CREATING FOLDER: MATLAB"
                folder_matlab=trim(adjustl(folder_output))// 'MATLAB/'
                command="mkdir "  // trim(adjustl(folder_matlab))
                if(rank.eq.0)write(*,*)trim(adjustl(command))
                if(rank.eq.0)call system(command)
endif
if(save_gmt.eq.1)       then
                if(rank.eq.0)write(*,*)
                if(rank.eq.0)write(*,*)"CREATING FOLDER: GMT"
                folder_gmt=trim(adjustl(folder_output))// 'GMT/'
                command="mkdir "  // trim(adjustl(folder_gmt))
                if(rank.eq.0)write(*,*)trim(adjustl(command))
                if(rank.eq.0)call system(command)
endif

end subroutine read_parfile

end module mod_parfile

module mod_SG_arrays
implicit none

        integer, allocatable :: shotID_1(:),shotID_n(:),sizeof(:)
	integer, allocatable :: shotID_nav(:),shotID_su(:),shotID_(:) 
	real, allocatable :: sx_su(:),sy_su(:)
	real, allocatable :: offset_su(:,:)
	real, allocatable :: pos_bat(:)

	integer, allocatable :: ifile_su(:)
	integer(4), allocatable :: pos_byte_su(:)

	integer, allocatable :: pos_shot_grid(:) !!position shot, no dimensions (points)
	integer, allocatable :: pos_trace_grid(:,:) !!position trace, no dimensions (points)

	real, allocatable :: pos_shot(:) !!position shot
	real, allocatable :: pos_trace(:,:) !!position trace

contains

subroutine allocate_SG_arrays()

USE mod_parfile
implicit none

        allocate(shotID_1(split_parts),shotID_n(split_parts),sizeof(split_parts))
        allocate(shotID_nav(NumShots))
        allocate(sx_su(NumShots),sy_su(NumShots))
        allocate(offset_su(NumRec,NumShots))
        allocate(pos_bat(NumShots))
	allocate(shotID_su(NumShots),shotID_(NumShots))
        allocate(pos_byte_su(NumShots),ifile_su(NumShots))
        allocate(pos_shot(NumShots),pos_shot_grid(NumShots))
        allocate(pos_trace(NumRec,NumShots),pos_trace_grid(NumRec,NumShots))

        shotID_1=0;shotID_n=0;sizeof=0;
        shotID_nav=0;
        sx_su=0;sy_su=0;
	offset_su=0;
	pos_bat=0;
	shotID_su=0;shotID_=0;
        pos_byte_su=0;ifile_su=0;
        pos_shot=0;pos_shot_grid=0;
        pos_trace=0;pos_trace_grid=0;

end subroutine	allocate_SG_arrays


subroutine deallocate_SG_arrays

        deallocate(shotID_1,shotID_n,sizeof)
        deallocate(shotID_nav)
        deallocate(sx_su,sy_su)
	deallocate(offset_su)
        deallocate(pos_bat)
	deallocate(shotID_su,shotID_)
        deallocate(pos_byte_su,ifile_su)
        deallocate(pos_shot,pos_shot_grid)
        deallocate(pos_trace,pos_trace_grid)


end subroutine deallocate_SG_arrays

end module mod_SG_arrays


module mod_input_arrays

implicit none

	integer	:: nx,ny
        integer, allocatable :: datum1_shot(:)
        integer, allocatable :: datum1_rec(:)
        integer, allocatable :: datum2(:)  
        real, allocatable :: DC_model(:,:)

contains

subroutine allocate_input_arrays(i)

implicit none
	
	integer	:: i

	if(i.eq.1)	then

	        allocate(datum1_shot(nx))
        	allocate(datum1_rec(nx))
	        allocate(datum2(nx))
	
		datum1_shot=0;
		datum1_rec=0;
		datum2=0;
	endif

	if(i.eq.2)	then

        	allocate(DC_model(ny,nx))
		DC_model=0;
	endif

end subroutine  allocate_input_arrays


subroutine deallocate_input_arrays()

        deallocate(datum1_shot)
        deallocate(datum1_rec)
        deallocate(datum2)
	deallocate(DC_model)

end subroutine deallocate_input_arrays

end module mod_input_arrays


module mod_PG_arrays
implicit none

	integer(4), allocatable :: pos_byte_PG(:)
	integer, allocatable :: NumShots_PG(:)
	integer, allocatable :: shot_grid_PG(:,:)
	integer, allocatable :: trace_grid_PG(:,:)
	integer, allocatable :: pos_trace_grid_PG(:,:)
	integer, allocatable :: PG_location(:)
	integer, allocatable :: pos_shot_grid_PG(:,:)
	integer, allocatable :: ifile_PG(:)


contains

subroutine allocate_PG_arrays()

USE mod_parfile

implicit none
integer :: NumPG,NumS

	NumPG=NumMax_PG
	NumS=NumMaxShots_PG

	allocate(pos_byte_PG(NumPG))
	allocate(NumShots_PG(NumPG))
	allocate(PG_location(NumPG))
	allocate(ifile_PG(NumPG))

	allocate(shot_grid_PG(NumS,NumPG))
	allocate(trace_grid_PG(NumS,NumPG))
	allocate(pos_trace_grid_PG(NumS,NumPG))
	allocate(pos_shot_grid_PG(NumS,NumPG))

	pos_byte_PG=0
	NumShots_PG=0
	PG_location=0
	ifile_PG=0

	shot_grid_PG=0
	trace_grid_PG=0
	pos_trace_grid_PG=0
	pos_shot_grid_PG=0

end subroutine	allocate_PG_arrays


subroutine deallocate_PG_arrays

	deallocate(pos_byte_PG)
	deallocate(NumShots_PG)
	deallocate(PG_location)
	deallocate(ifile_PG)

	deallocate(shot_grid_PG)
	deallocate(trace_grid_PG)
	deallocate(pos_trace_grid_PG)
	deallocate(pos_shot_grid_PG)
	
end subroutine deallocate_PG_arrays

end module mod_PG_arrays
