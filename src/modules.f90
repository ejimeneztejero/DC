!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains the modules:
!!	(mod_parfile, mod_data_arrays, mod_model_arrays)
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module mod_parfile

implicit none

  ! Control file variables
  integer, parameter :: unit0=100
  integer, parameter :: unit_DC=300

  INTEGER(4), parameter :: size_su_header = 60
  INTEGER(4), parameter :: byte_fldr = 9
  INTEGER(4), parameter :: byte_tracr = 5
  INTEGER(4), parameter :: byte_scalco = 71
  INTEGER(4), parameter :: byte_sx = 73
  INTEGER(4), parameter :: byte_sy = 77

  INTEGER(4) :: byte_shotnumber
  
  integer :: nt,NumShots,NumOBS
  integer :: shot_init,shot_fin

  integer :: endianness_data,endianness_machine
  integer :: sx_sy_header
  integer :: TWT_option
  integer :: save_txt
  integer :: save_gmt,save_matlab

  integer :: maxbytes
  
  real :: dmodel,dt
  real :: shot_depth
  real :: added_space_model_X,added_space_model_Y,water_velocity
  real :: dshots,time

  real :: f1,f2

  character(len=500) :: par_file
  character(len=500) :: folder_input,folder_output
  character(len=500) :: su_file0
  character(len=500) :: su_file_DC
  character(len=500) :: temp_DC
  character(len=500) :: nav_file,vp_file, obs_data

  character(len=500), allocatable :: original_file(:)

  contains

subroutine read_parfile(rank)

  implicit none

  ! Input related variables
  character(len=200) :: buffer,label
  integer :: pos,icount,i,interval,rank,nlines
  integer, parameter :: fh = 10
  integer :: ios = 0
  integer :: line = 0

  character(len=500) :: command0,command
  character(len=500) :: file_name,file_name2,file_name3
  character(len=50) :: Str,Str_su,access,form,OBS_num

  logical :: su_exist,nav_exist,vp_exist,obs_data_exist
  logical :: input_exist, output_exist

  access = 'STREAM'
  form = 'UNFORMATTED'

  icount = iargc()  ! The number does not includes the executable name, so if user passed one argument, 1 is returned.
  if(rank.eq.0)write(*,*)'Number of arguments passed after executable name: ',icount
  if(rank.eq.0)write(*,*)'Your Parfile containing input parameteres must be passed the last one in the command execution line'

!  if ( icount.eq.3 ) then
	call getarg(icount, par_file)	! The file name of the executable.
	if(rank.eq.0)write(*,*)'name par_file: ',trim(adjustl(par_file))
	file_name = trim(adjustl(par_file))
	INQUIRE(file=file_name,EXIST=su_exist)
        if(.NOT. su_exist)     then
                if(rank.eq.0)write(*,*)'ERROR: Par_file named: ', trim(adjustl(par_file)),' does not exist'
                if(rank.eq.0)call ascii_art(2)
                stop
        endif
!  endif

  write(*,*)'par_file: ',trim(adjustl(par_file))

  su_file0 = 'null'
  nav_file = 'null'
  vp_file = 'null'
  su_file_DC = 'null'
  temp_DC = 'temp_DC'

  NumOBS=1
  byte_shotnumber=byte_tracr
  sx_sy_header=0;
  TWT_option=0;
  endianness_data=1;endianness_machine=0;
  added_space_model_X=0.
  added_space_model_Y=0.
  water_velocity=1500.
  dmodel=25;
  dt=0;
  shot_depth=0;
  nt=0;
  save_gmt=0;save_matlab=0;

  open(fh, file=par_file)

  ! ios is negative if an end of record condition is encountered or if
  ! an endfile condition was detected. It is positive if an error was
  ! detected.  ios is zero otherwise.

  do while (ios == 0)
     read(fh, '(A)', iostat=ios) buffer
     if (ios == 0) then
        line = line + 1

        ! Find the first instance of whitespace.
        pos = scan(buffer,' ')
        label = buffer(1:pos)
        buffer = buffer(pos+1:)

        select case (label)

        case ('byte_shotnumber:')
           read(buffer, *, iostat=ios) byte_shotnumber
        case ('sx_sy_header:')
           read(buffer, *, iostat=ios) sx_sy_header
        case ('TWT_option:')
           read(buffer, *, iostat=ios) TWT_option
        case ('endianness_data:')
           read(buffer, *, iostat=ios) endianness_data
        case ('endianness_machine:')
           read(buffer, *, iostat=ios) endianness_machine
        case ('obs_file_list:')
           read(buffer, *, iostat=ios) obs_data
        case ('nav_file:')
           read(buffer, *, iostat=ios) nav_file
        case ('su_file:')
           read(buffer, *, iostat=ios) su_file0
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
        case ('dmodel:')
           read(buffer, *, iostat=ios) dmodel
        case ('dshots:')
           read(buffer, *, iostat=ios) dshots
        case ('shot_init:')
           read(buffer, *, iostat=ios) shot_init
        case ('shot_fin:')
           read(buffer, *, iostat=ios) shot_fin
        case ('NumOBS:')
           read(buffer, *, iostat=ios) NumOBS
        case ('ns:')
           read(buffer, *, iostat=ios) nt
        case ('nt:')
           read(buffer, *, iostat=ios) nt
        case ('shot_depth:')
           read(buffer, *, iostat=ios) shot_depth
        case ('water_velocity:')
           read(buffer, *, iostat=ios) water_velocity
        case ('save_gmt:')
           read(buffer, *, iostat=ios) save_gmt
        case ('save_matlab:')
           read(buffer, *, iostat=ios) save_matlab
!        case default
!           if(rank.eq.0)print *, 'WARNING in file ',trim(adjustl(par_file)),': skipping invalid label at line', line
        end select
     end if
  end do

 close(fh)

folder_input = trim(adjustl(folder_input)) // '/'
folder_output = trim(adjustl(folder_output)) // '/'

time=(nt-1)*dt

if(dmodel.eq.0)dmodel=dshots

if(dmodel.gt.dshots)	then
	if(rank.eq.0)write(*,*)'dmodel cannot be greather than dshots, instead, it is set to dshots'
	dmodel=dshots
endif

NumShots=shot_fin-shot_init+1

if(added_space_model_X.eq.0)added_space_model_X=20.*dmodel
if(added_space_model_Y.eq.0)added_space_model_Y=20.*dmodel

if(sx_sy_header.ne.0.and.sx_sy_header.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: sx_sy_header should be set to 0 or 1 in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)
	stop
endif

if(TWT_option.ne.0.and.TWT_option.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: TWT_option should be set to 0 or 1 in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)
	stop
endif

if(endianness_data.ne.0.and.endianness_data.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: endianness_data should be set to 0 (little endian)&
	 or 1 (big endian) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

	stop
endif

if(endianness_machine.ne.0.and.endianness_machine.ne.1)	then
	if(rank.eq.0)write(*,*)'ERROR: endianness_machine should be set to 0 (little endian)&
	 or 1 (big endian) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

	stop
endif

file_name = trim(adjustl(folder_input)) // trim(adjustl(nav_file))
INQUIRE(FILE=file_name, EXIST=nav_exist)
if(.NOT. nav_exist)	then
	if(rank.eq.0)write(*,*)'ERROR: nav_file not found'
        if(rank.eq.0)call ascii_art(2)

	stop
endif

file_name = trim(adjustl(folder_input)) // trim(adjustl(vp_file))
INQUIRE(FILE=file_name, EXIST=vp_exist)
if(.NOT. vp_exist)	then
if(water_velocity.le.1400)	then
	if(rank.eq.0)write(*,*) 'ERROR: vp_file not found in your your input folder &
	and water_velocity is set to a wrong value in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

	stop
endif
endif


file_name = trim(adjustl(folder_input)) // trim(adjustl(obs_data))
INQUIRE(FILE=file_name, EXIST=obs_data_exist)
if(.NOT. obs_data_exist)     then
        write(*,*)'ERROR: "file_data_list: " file not found'
        !error=1
else
    	!!      OBS DATA FILE
        file_name= trim(adjustl(folder_input)) // trim(adjustl(obs_data))
        open(unit=10,file=file_name,status='old')
        nlines=0
        do
          	read(10,*, END=10)
                nlines = nlines + 1
        enddo
	10 close (10)

        if(nlines.lt.NumOBS)    then
                write(*,*)'ERROR: ', abs(nlines-NumOBS),' OBS data files are missing in obs_file_list'
!                error=1
        	if(rank.eq.0)call ascii_art(2)
        	stop
        endif

endif

if(dt.eq.0) then
	if(rank.eq.0)write(*,*) 'ERROR: Please give a value to dt (seconds) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)
	stop
endif

if(nt.eq.0) then
	if(rank.eq.0)write(*,*)'ERROR: Please give a value to nt in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

	stop
endif

if(shot_depth.eq.0) then
	if(rank.eq.0)write(*,*)'ERROR: Please give a value to shot_depth (meters) in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

	stop
endif

if(shot_init.le.0) then
	if(rank.eq.0)write(*,*) 'ERROR: Please give a correct value to shot_init in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

	stop
endif

if(shot_fin.lt.shot_init) then
	if(rank.eq.0)write(*,*) 'ERROR: Please give a correct value to shot_fin in ',trim(adjustl(par_file))
        if(rank.eq.0)call ascii_art(2)

	stop
endif

allocate(original_file(NumOBS))
file_name= trim(adjustl(folder_input)) // trim(adjustl(obs_data))
open(unit=10,file=file_name,status='old')
do i=1,NumOBS
        read(10,*)file_name2
        original_file(i) = trim(adjustl(file_name2))
enddo
close(10)

maxbytes=1900000000

if(save_matlab.ne.0.or.save_gmt.ne.0)save_txt=1

if(rank.eq.0)	then

write(*,*)
write(*,*)'*******************************************'
write(*,*)'RUNNING PROGRAM FOR NEXT SET OF PARAMETERS:'
write(*,*)'*******************************************'

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
	if(TWT_option.eq.0)write(*,*)'nav_file should contain: ShotID	Xshot(UTM)	Yshot(UTM)	Z(m)'
	if(TWT_option.eq.1)write(*,*)'nav_file should contain: ShotID	Xshot(UTM)	Yshot(UTM)	TWT(s)'
endif

write(*,*)'-------------------------------------'
write(*,*)'Input folder: ',adjustl(trim(folder_input))
write(*,*)'Output folder: ',adjustl(trim(folder_output))
write(*,*)'name of file with the name of the obs data files: ',adjustl(trim(obs_data))
write(*,*)'name of navigation file for shotgathers: ',adjustl(trim(nav_file))
write(*,*)'byte_shotnumber: ',byte_shotnumber
write(*,*)'NumOBS: ',NumOBS
write(*,*)'NumShots: ',NumShots
write(*,*)'shot_init: ',shot_init
write(*,*)'shot_fin: ',shot_fin
write(*,*)'nt: ',nt
write(*,*)'dt (s): ',dt
write(*,*)'dshots (m): ',dshots
write(*,*)'dmodel (m): ',dmodel
write(*,*)'shot_depth (m): ',shot_depth
if(.NOT. vp_exist)write(*,*)'water_velocity constant: ',water_velocity

endif

end subroutine read_parfile

subroutine read_OBSfile(iOBS,rank)

  implicit none

  integer :: rank,iOBS
  character(len=500) :: file_name

  logical :: su_exist

  file_name = trim(adjustl(folder_input)) // original_file(iOBS)

  INQUIRE(file=file_name,EXIST=su_exist)
  if(.NOT. su_exist)     then
       write(*,*)'WARNING: SU Field data: ',trim(adjustl(file_name)),' not found'
        if(rank.eq.0)call ascii_art(2)
	stop
  endif

end subroutine read_OBSfile

end module mod_parfile

module mod_data_arrays
implicit none

	integer, allocatable :: shotID_nav(:),shotID_su(:),shotID_(:) 

	real, allocatable :: sx_su(:),sy_su(:)
	real, allocatable :: pos_bat(:)

	integer(4), allocatable :: pos_byte_su(:)

	integer, allocatable :: pos_shot_grid(:) !!position shot, no dimensions (points)

	real, allocatable :: pos_shot(:) !!position shot

	INTEGER :: shotID_1, shotID_n, sizeof

contains

subroutine allocate_data_arrays()

USE mod_parfile
implicit none

        allocate(shotID_nav(NumShots))
        allocate(sx_su(NumShots),sy_su(NumShots))
        allocate(pos_bat(NumShots))
	allocate(shotID_su(NumShots),shotID_(NumShots))
        allocate(pos_byte_su(NumShots))
        allocate(pos_shot(NumShots),pos_shot_grid(NumShots))

        shotID_1=0;shotID_n=0;sizeof=0;
        shotID_nav=0;
        sx_su=0;sy_su=0;
	pos_bat=0;
	shotID_su=0;shotID_=0;
        pos_byte_su=0;
        pos_shot=0;pos_shot_grid=0;

end subroutine	allocate_data_arrays

subroutine deallocate_data_arrays

        deallocate(shotID_nav)
        deallocate(sx_su,sy_su)
        deallocate(pos_bat)
	deallocate(shotID_su,shotID_)
        deallocate(pos_byte_su)
        deallocate(pos_shot,pos_shot_grid)

end subroutine deallocate_data_arrays

end module mod_data_arrays

module mod_model_arrays

implicit none

	integer	:: nx,ny
        integer, allocatable :: datum1(:)
        integer, allocatable :: datum2(:)  
        real, allocatable :: DC_model(:,:)

contains

subroutine allocate_model_arrays(i)

implicit none
	
	integer	:: i

	if(i.eq.1)	then

	        allocate(datum1(nx))
	        allocate(datum2(nx))
	
		datum1=0;
		datum2=0;
	endif

	if(i.eq.2)	then

        	allocate(DC_model(ny,nx))
		DC_model=0;
	endif

end subroutine  allocate_model_arrays


subroutine deallocate_model_arrays()

        deallocate(datum1)
        deallocate(datum2)
	deallocate(DC_model)

end subroutine deallocate_model_arrays

end module mod_model_arrays
