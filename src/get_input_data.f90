!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains subroutines for extracting input data features
!!	Author:	Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine inputs(iOBS,rank,numtasks)

USE mod_parfile
USE mod_model_arrays
USE mod_data_arrays

implicit none

        integer :: i,iOBS,rank,numtasks
        integer :: NumBat,nmodel
        integer :: shotID_first_su,shotID_last_su	
        integer :: shotID_first_nav,shotID_last_nav     

        real    :: length_model,add1

	integer, allocatable :: bat_model_grid(:)

        character(len=500) :: file_name_water,file_name

        logical :: file_exists


call read_OBSfile(iOBS,rank)

call open_su_files(iOBS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!	SHOTGATHERS NAVIGATION

if(rank.eq.0)	then
        write(*,*)
        write(*,*)'*********************'
        write(*,*)'READING NAVIGATION'
        write(*,*)'*********************'

endif

call get_navigation_file(iOBS,shotID_first_nav,shotID_last_nav)	!fuera del loop

if(rank.eq.0)	then
        write(*,*)
        write(*,*)'*********************'
        write(*,*)'READING FIELD DATA'
        write(*,*)'*********************'
endif

call get_su_files1(iOBS,shotID_first_su,shotID_last_su)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Warnings
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(shotID_first_su.ne.shotID_first_nav.or.shotID_last_su.ne.shotID_last_nav)    then
	write(*,*)'ERROR: first/last shotID in .su and navigation file must coincide'
        call ascii_art(2)
        stop
endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

call get_su_files2(iOBS)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Lectura de geometria_ shots, traces, dimensiones modelo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

if(rank.eq.0)	then
        write(*,*)
        write(*,*)'********************************************'
        write(*,*)'GETTING GEOMETRY: POSITION of SHOTS'
        write(*,*)'********************************************'
endif

call get_geometry(NumBat,length_model,nmodel,add1)	!!fuera del loop

do i=1,NumShots

	if(shotID_su(i).ne.shotID_nav(i))   then
                        write(*,*)'ERROR: shotID must be coincident in su and nav file'
                        call ascii_art(2)
                        stop
                endif
        enddo

        shotID_=shotID_nav*shotID_su    !!multiplication, to just work with shotID_ different than 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Lectura e interpolacion de batimetria
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	if(rank.eq.0)	then

        write(*,*)
        write(*,*)'****************************'
        write(*,*)'CALCULATING BATHYMETRY MODEL'
        write(*,*)'****************************'

	endif

        nx=nmodel

        allocate(bat_model_grid(nx))
        bat_model_grid=0;

	call allocate_model_arrays(1)

        call get_bathymetry_model(NumBat,length_model,nmodel,bat_model_grid)	!!fuera del loop

        datum1=1+ceiling(shot_depth/dmodel)        	!! superficie 1 shots
        datum2=bat_model_grid                         	!! superficie 2

        ny=maxval(datum2)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!     Superficies datum
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	call allocate_model_arrays(2)

        DC_model=water_velocity

	if(rank.eq.0)	then
	        write(*,*)'Full line model dimensions (points): ',&
	        nx,' length, ',ny,' depth'
	        write(*,*)'In km: ',&
	        (nx-1)*dmodel/1000.,' length, ',&
       		(ny-1)*dmodel/1000.,' depth'
	endif

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!------Lectura de modelo Vp agua, constante o interpolando XBT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        file_name_water= trim(adjustl(folder_input)) // trim(adjustl(vp_file))
        INQUIRE(FILE=file_name_water, EXIST=file_exists)

        if(file_exists) then

		if(rank.eq.0)	then

	                write(*,*)
        	        write(*,*)'Building Vp model for water column ...'
                	write(*,*)

		endif

                call get_Vp_model(add1)

        endif

        deallocate(bat_model_grid)

end subroutine inputs

subroutine get_bathymetry_model(NumBat,length_model,nmodel,bat_model_grid)

USE mod_parfile
USE mod_data_arrays

implicit none

	integer :: i,j,ERR
	integer :: nmodel,NumBat
	real :: length_model

	integer :: bat_model_grid(nmodel)

	real, allocatable :: bat_model(:),xmodel(:),xmodel_(:),xbat_(:)
	real, allocatable :: bat_pos_shot(:),bat_model_shot(:)
	character (len=500) :: file_name

	allocate(bat_model(nmodel),xmodel(nmodel))
	allocate(xmodel_(NumBat+2),xbat_(NumBat+2))
	allocate(bat_model_shot(NumBat),bat_pos_shot(NumBat))

	i=0
	do j=1,NumShots

		if(pos_bat(j).gt.0)	then
			i=i+1			
			bat_model_shot(i)=pos_bat(j)
			bat_pos_shot(i)=pos_shot(j)	
		endif

	enddo

	do j=2,NumBat+1
		xbat_(j)=bat_model_shot(j-1)
		xmodel_(j)=bat_pos_shot(j-1)
	enddo

	xmodel_(1)=0;xmodel_(NumBat+2)=length_model;!!x
	xbat_(1)=xbat_(2);xbat_(NumBat+2)=xbat_(NumBat+1);!!y

	do j=1,nmodel
		bat_model(j)=0
		xmodel(j)=(j-1)*dmodel
	enddo

        ERR=0
	call INTRPL(NumBat+2,xmodel_,xbat_,nmodel,xmodel,bat_model,ERR)

	do j=1,nmodel
		bat_model_grid(j)=1+ceiling(bat_model(j)/dmodel)
	enddo

	file_name=trim(adjustl(folder_output))//'bathymetry_meters.txt'
	open(unit=12,file=file_name,status='unknown')
	do j=1,nmodel
		write(12,*)xmodel(j),-bat_model(j)
	enddo
	close(12)

	deallocate(bat_model,xmodel)
	deallocate(xmodel_,xbat_)
	deallocate(bat_model_shot,bat_pos_shot)

end subroutine get_bathymetry_model

subroutine get_geometry(NumBat,length_model,nmodel,add1)

USE mod_parfile
USE mod_data_arrays

implicit none

        integer :: j,k,ii,ishot
	integer :: iline,itr,NumBat,ERR,nlines
        integer :: shotID,nmodel

        real :: xs,ys,bat,x1,y1,twt,add1
        real :: length_model

        character (len=500) :: file_name

!!	LECTURA DATOS SHOTGATHERS

        file_name= trim(adjustl(folder_input)) // trim(adjustl(nav_file))
        open(unit=10,file=file_name,status='old')
	nlines=0
	do
		read(10,*, END=20)
		nlines = nlines+1
	enddo
	20 close(10)

        open(unit=10,file=file_name,status='old')

        NumBat=0

        do ishot=1,nlines !!number of lines nav_file

		shotID=0;xs=0;ys=0;twt=0;bat=0;

!!!		POSITION (UTM) OF SHOT GATHERS, AND BAT 

		if(sx_sy_header.eq.0)	then

			if(TWT_option.eq.0)	then
				read(10,*)shotID,xs,ys,bat
			endif

			if(TWT_option.eq.1)	then
				read(10,*)shotID,xs,ys,twt
				bat=water_velocity*twt/2.
			endif

	                ii=shotID-shot_init+1
	       	        shotID_nav(ii)=shotID

		endif!sx_sy_header

		if(sx_sy_header.eq.1)	then

			if(TWT_option.eq.0)	then
				read(10,*)shotID,bat
			endif

			if(TWT_option.eq.1)	then
				read(10,*)shotID,twt
				bat=water_velocity*twt/2.
			endif

	                ii=shotID-shot_init+1	!!ii, "icount"
	       	        shotID_nav(ii)=shotID
			xs=sx_su(ii)
			ys=sy_su(ii)

		endif

		add1=added_space_model_X

                if(ii.eq.1)  then
                        x1=xs;y1=ys
                endif

                pos_shot(ii) = add1+sqrt((xs-x1)**2.+(ys-y1)**2.) !! distancia entre cada shot y el primero
                pos_shot_grid(ii)=1+ceiling(pos_shot(ii)/dmodel)

                if(bat.ne.bat.or.bat.eq.0)      then    !localiza NaN y CEROS
                        bat=0
                        pos_bat(ii)=0;
                else
                        pos_bat(ii)=abs(bat)
                        NumBat=NumBat+1
                endif

        enddo

        file_name=trim(adjustl(folder_output))//'shots_position.txt'
        open(unit=12,file=file_name,status='unknown')
        do j=1,NumShots
                write(12,*)j,pos_shot(j),-shot_depth
        enddo
        close(12)

	length_model=pos_shot(NumShots)+add1

        nmodel=1+ceiling(length_model/dmodel)

	CLOSE(10)

end subroutine get_geometry

subroutine get_navigation_file(iOBS,ID_first_nav,ID_last_nav)

USE mod_parfile
USE mod_data_arrays

implicit none

	integer :: i,nlines,num_id,iOBS,rank
	integer :: ID_first_nav,ID_last_nav

	CHARACTER(len=50) :: access,form,Str
	CHARACTER(len=500) :: file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	NAVIGATION FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	file_name= trim(adjustl(folder_input)) // trim(adjustl(nav_file))

	open(unit=10,file=file_name,status='old')
	nlines=0
	do
	    read(10,*, END=10) 
	    nlines = nlines + 1 
	enddo 
	10 close (10) 

	open(unit=10,file=file_name,status='old')
	do i=1,nlines
		read(10,*)num_id
		if(i.eq.1)ID_first_nav=num_id
		if(i.eq.nlines)ID_last_nav=num_id
	enddo
	close(10)

	write(*,*)'Navigation file contains shots from ',ID_first_nav,' to ',ID_last_nav	
	write(*,*)'NumShots total:', ID_last_nav-ID_first_nav+1,', Num Lines Nav_file:',nlines!

	if(nlines.ne.ID_last_nav-ID_first_nav+1)	then
		write(*,*)'WARNING: navigation for ',abs(nlines-(ID_last_nav-ID_first_nav+1)),' shots are missing in nav_file'
	endif

end subroutine get_navigation_file

subroutine get_su_files1(iOBS,shotID_first_su,shotID_last_su)

USE mod_parfile
USE mod_data_arrays

implicit none

	INTEGER :: i,j,k,iOBS
	INTEGER :: shotID_first_su,shotID_last_su
	INTEGER :: nh

	INTEGER(4) :: pos_read

	CHARACTER(len=50) :: Str
	CHARACTER(len=500) :: file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	.SU FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! TO READ ShotID:

pos_read=byte_shotnumber
nh = size_su_header	!Size su header = 60*4 bytes

! READ SU FILES
!file_name = original_file
file_name = trim(adjustl(folder_input)) // original_file(iOBS)

READ(unit0+iOBS,pos=pos_read) shotID_1
READ(unit0+iOBS,pos=sizeof-(nh+nt)*4+pos_read) shotID_n

write(*,*)trim(adjustl(file_name)),&
', contains shots from ',shotID_1,' to ',shotID_n

READ(unit0+iOBS,pos=pos_read) shotID_first_su
READ(unit0+iOBS,pos=sizeof-(nh+nt)*4+pos_read) shotID_last_su

end subroutine get_su_files1

subroutine get_su_files2(iOBS)

USE mod_parfile
USE mod_data_arrays

implicit none

	integer :: i,j,k,icount,iOBS
	INTEGER :: nh

	INTEGER*2 :: pos_scalco,scalco

	INTEGER*4 :: shotID,pos_read
	INTEGER*4 :: pos_byte,pos_byte2 

	INTEGER*4 :: pos_sx,pos_sy
	INTEGER*4 :: sx,sy

	CHARACTER(len=500) :: file_name

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!	.SU FILE
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! TO READ ShotID:

nh = size_su_header     !Size su header = 60*4 bytes
pos_read=byte_shotnumber
pos_scalco= byte_scalco
pos_sx = byte_sx
pos_sy = byte_sy

!------------------------------------------------------------------------
! Loop shot-by-shot

pos_byte=1
icount=1

READ(unit0+iOBS,pos=pos_read) shotID

if(sx_sy_header.eq.1)	then

	READ(unit0+iOBS,pos=pos_byte+pos_sx-1) sx
	READ(unit0+iOBS,pos=pos_byte+pos_sy-1) sy
	READ(unit0+iOBS,pos=pos_byte+pos_scalco-1) scalco

endif

do while(shotID.ge.shot_init.and.shotID.le.shot_fin.and.icount.le.NumShots)

	shotID_su(icount)=shotID	
	pos_byte_su(icount)=pos_byte

	if(sx_sy_header.eq.1)	then

        	if(sx.eq.0)   then

                	write(*,*)'ERROR: sx in header is 0. &
                        set sx_sy_header: 0 in file: ',trim(adjustl(par_file)),&
                        ' or include geometry in .sgy file'
                	call ascii_art(2)

			stop

        	endif

        	if(sy.eq.0)   then

                	write(*,*)'ERROR: sy in header is 0. &
                        set sx_sy_header: 0 in file: ',trim(adjustl(par_file)),&
                        ' or include geometry in .sgy file'
                	call ascii_art(2)

			stop

        	endif

		sx_su(icount)=abs(sx)
		sy_su(icount)=abs(sy)

		if(scalco.gt.0)	then
			sx_su(icount)=abs(sx)*scalco
			sy_su(icount)=abs(sy)*scalco
		endif

		if(scalco.lt.0)	then
			sx_su(icount)=abs(sx)/abs(scalco)
			sy_su(icount)=abs(sy)/abs(scalco)
		endif

	endif

	pos_byte=pos_byte+(nh+nt)*4	

        if(pos_byte.ge.sizeof) then

		exit

	endif

	if(pos_byte.gt.0.and.pos_byte.lt.sizeof)	then

		READ(unit0+iOBS,pos=pos_byte+pos_read-1) shotID
		if(sx_sy_header.eq.1)	then
			READ(unit0+iOBS,pos=pos_byte+pos_sx-1) sx
			READ(unit0+iOBS,pos=pos_byte+pos_sy-1) sy
			READ(unit0+iOBS,pos=pos_byte+pos_scalco-1) scalco
		endif

		icount=shotID-shot_init+1

	endif

enddo

do k=1,NumShots
	if(shotID_su(k).eq.0)	then
		write(*,*)'WARNING, not recorded in .su file, shotID: ',k+shot_init-1
	endif
enddo

end subroutine get_su_files2

subroutine get_Vp_model(add1)

USE mod_parfile
USE mod_model_arrays
USE mod_data_arrays

implicit none

	integer :: ddx,i,j,topo(nx),nlines,ERR,imax
	integer :: maxlines,shotline,maxl
	integer :: nxx,nyy,addx
	integer :: npos
	real ::    add1
        character(len=500) :: file_name

	integer, allocatable :: nlines_(:)
	real, allocatable :: x(:),z(:,:),x2(:),y2(:),zb(:,:)
	real, allocatable :: vpm(:,:),vpmb(:,:),vpm2(:,:),vpm3(:,:)
	real, allocatable :: xpos(:),vpm4(:,:)

	character(len=500), allocatable :: vpf(:)

	ERR=0

!	addx=nint(add1/dmodel)
!	ddx=1+floor(added_space_model_X/dmodel)
	addx=ceiling(add1/dmodel)
	ddx=1+ceiling(added_space_model_X/dmodel)

	nxx=NumShots
	nyy=ny

	allocate(x2(nxx),y2(nyy))

	do i=1,nxx
		x2(i)=i
	enddo
	do j=1,nyy
		y2(j)=(j-1)*dmodel
	enddo

	DC_model=water_velocity

	file_name=trim(adjustl(folder_input)) // trim(adjustl(vp_file))

	open(10,file=file_name,status='old')
        nlines=0
	do
		read(10,*, END=10)
		nlines=nlines+1
	enddo
	10 close(10)

	allocate(vpf(nlines),nlines_(nlines),x(nlines))

	open(10,file=file_name,status='old')
	do i=1,nlines

		read(10,*)shotline,vpf(i)
	        x(i)=shotline-shot_init+1
	enddo
	close(10)

	do i=1,nlines

		file_name=trim(adjustl(folder_input))//trim(adjustl(vpf(i)))

		open(10,file=file_name,status='old')
	     	nlines_(i)=0
	     	do
	    		read(10,*, END=11)
			nlines_(i) = nlines_(i) + 1
		enddo
		11 close(10)
		
	enddo

	maxlines=maxval(nlines_)

	allocate(z(maxlines,nlines),vpm(maxlines,nlines))

	do i=1,nlines

		maxl=nlines_(i)

		if(maxl.eq.maxlines)imax=i

		file_name=trim(adjustl(folder_input))//trim(adjustl(vpf(i)))
		open(10,file=file_name,status='old')
		do j=1,nlines_(i)

			read(10,*)z(j,i),vpm(j,i)

		enddo	
		close(10)	
	enddo

	do i=1,nlines
		maxl=nlines_(i)
		if(maxl.lt.maxlines) z(maxl+1:maxlines,i)=z(maxl+1:maxlines,imax)
		if(maxl.lt.maxlines) vpm(maxl+1:maxlines,i)=vpm(maxl,i)
	enddo

	if(z(maxlines,imax).gt.y2(nyy))	then

!!		write(*,*)'z>ny2'

	!!	Interpolacion de las lineas de (z,vp) a todo el espacio (nx,ny) con resolucion dmodel

		allocate(vpm2(nyy,nlines))

	!!	Interpolacion en ny 
		do i=1,nlines
			call  INTRPL(maxlines,z(:,i),vpm(:,i),nyy,y2,vpm2(:,i),ERR)
		enddo

	endif

	!!	Si no llegan las lineas hasta el fondo, el programa las prolonga con velocidades constantes

	if(z(maxlines,imax).lt.y2(nyy))	then

!!		write(*,*)'z<ny2'

		allocate(zb(maxlines+1,nlines),vpmb(maxlines+1,nlines))
		
		zb(1:maxlines,:)=z(:,:)
		vpmb(1:maxlines,:)=vpm(:,:)

		zb(maxlines+1,:)=y2(nyy)
		vpmb(maxlines+1,:)=vpm(maxlines,:)

	!!	Interpolacion de las lineas de (z,vp) a todo el espacio (nx,ny) con resolucion dmodel

		allocate(vpm2(nyy,nlines))

	!!	Interpolacion en ny
		do i=1,nlines
			call  INTRPL(maxlines+1,zb(:,i),vpmb(:,i),nyy,y2,vpm2(:,i),ERR)
		enddo

		deallocate(zb,vpmb)

	endif

	allocate(vpm3(nyy,nxx))

!!	Interpolacion en shotID
	do j=1,nyy
		call INTRPL(nlines,x,vpm2(j,:),nxx,x2,vpm3(j,:),ERR)
	enddo

	do i=1,nxx

		x2(i)=pos_shot(i)

		if(x2(i).eq.0)	then
			x2(i)=pos_shot(i-1)+dshots
			write(*,*)'x2=0,line nav_file:',i
		endif

	enddo

!	npos= 1+nint( (pos_shot(NumShots)-pos_shot(1))/dmodel )
	npos= 1+ceiling( (pos_shot(NumShots)-pos_shot(1))/dmodel )
	allocate(xpos(npos),vpm4(ny,npos))

	do i=1,npos
		xpos(i)=pos_shot(1)+(i-1)*dmodel
	enddo


	if(xpos(npos).gt.x2(nxx))	then
		xpos(npos)=x2(nxx)
	endif

	!!	Interpolacion en grid
	do j=1,ny
		call  INTRPL(nxx,x2,vpm3(j,:),npos,xpos,vpm4(j,:),ERR)
	enddo

	do i=1,addx
		do j=1,ny
			DC_model(j,i)=vpm4(j,1)
		enddo
	enddo

	do i=1,npos
		do j=1,ny
			DC_model(j,addx+i)=vpm4(j,i)
		enddo
	enddo
	
	do i=addx+npos+1,nx
		do j=1,ny
			DC_model(j,i)=vpm4(j,npos)
		enddo
	enddo

		do i=1,nx
			do j=1,ny
				if(DC_model(j,i).eq.0)write(*,*)'WARNING Vp model is wrong'
			enddo
		enddo

	        file_name=trim(adjustl(folder_output))//'Vp_water_model.txt'
		open(unit=12,file=file_name,status='unknown') 	
		do j=1,ny
			write(12,'(20000(e12.5,2x))') (DC_model(j,i),i=1,nx)
		enddo
		close(12)

	deallocate(vpm,vpm2,vpm3,vpm4)
	deallocate(z,xpos,x2,y2)
        deallocate(vpf,nlines_,x)


end subroutine get_Vp_model


SUBROUTINE  INTRPL(L,X,Y,N,U,V,ERR)

! INTERPOLATION OF A SINGLE-VALUED FUNCTION
! THIS SUBROUTINE INTERPOLATES, FROM VALUES OF THE FUNCTION
! GIVEN AS ORDINATES OF INPUT DATA POINTS IN AN X-Y PLANE
! AND FOR A GIVEN SET OF X VALUES (ABSCISSAS), THE VALUES OF
! A SINGLE-VALUED FUNCTION Y = Y(X).
! THE INPUT PARAMETERS ARE
!     L  = NUMBER OF INPUT DATA POINTS
!          (MUST BE 2 OR GREATER)
!     X  = ARRAY OF DIMENSION L STORING THE X VALUES
!          (ABSCISSAS) OF INPUT DATA POINTS
!          (IN ASCENDING ORDER)
!     Y  = ARRAY OF DIMENSION L STORING THE Y VALUES
!          (ORDINATES) OF INPUT DATA POINTS
!     N  = NUMBER OF POINTS AT WHICH INTERPOLATION OF THE
!          Y VALUE (ORDINATE) IS DESIRED
!          (MUST BE 1 OR GREATER)
!     U  = ARRAY OF DIMENSION N STORING THE X VALUES
!          (ABSCISSAS) OF DESIRED POINTS
! THE OUTPUT PARAMETER IS
!     V  = ARRAY OF DIMENSION N WHERE THE INTERPOLATED Y
!          VALUES (ORDINATES) ARE TO BE DISPLAYED
!   ERR  = ERROR CODE (added by ThPe)
! DECLARATION STATEMENTS
      INTEGER :: L, N, ERR
      REAL  ::  X(L),Y(L),U(N),V(N)
      EQUIVALENCE  (P0,X3),(Q0,Y3),(Q1,T3)
      REAL ::    M1,M2,M3,M4,M5
      REAL ::    A2,A4
      EQUIVALENCE  (UK,DX),(IMN,X2,A1,M1),(IMX,X5,A5,M5),(J,SW,SA),(Y2,W2,W4,Q2),(Y5,W3,Q3)
! PRELIMINARY PROCESSING
    10 L0=L
      LM1=L0-1
      LM2=LM1-1
      LP1=L0+1
      N0=N
      M2=0
      A2=0
      A4=0
      IF(LM2.LT.0)        GO TO 90
      IF(N0.LE.0)         GO TO 91
      DO 11  I=2,L0
!        IF(X(I-1)-X(I))   11,95,96
        IF(X(I-1)-X(I).LT.0)   GO TO 11
        IF(X(I-1)-X(I).EQ.0)   GO TO 95
        IF(X(I-1)-X(I).GT.0)   GO TO 96
   11   CONTINUE
      IPV=0
! MAIN DO-LOOP
      DO 80  K=1,N0
        UK=U(K)
! ROUTINE TO LOCATE THE DESIRED POINT
   20   IF(LM2.EQ.0)      GO TO 27
        IF(UK.GE.X(L0))   GO TO 26
        IF(UK.LT.X(1))    GO TO 25
        IMN=2
        IMX=L0
   21   I=(IMN+IMX)/2
        IF(UK.GE.X(I))    GO TO 23
   22   IMX=I
        GO TO 24
   23   IMN=I+1
   24   IF(IMX.GT.IMN)    GO TO 21
        I=IMX
        GO TO 30
   25   I=1
        GO TO 30
   26   I=LP1
        GO TO 30
   27   I=2
! CHECK IF I = IPV
   30   IF(I.EQ.IPV)      GO TO 70
        IPV=I
! ROUTINES TO PICK UP NECESSARY X AND Y VALUES AND
!          TO ESTIMATE THEM IF NECESSARY
   40   J=I
        IF(J.EQ.1)        J=2
        IF(J.EQ.LP1)      J=L0
        X3=X(J-1)
        Y3=Y(J-1)
        X4=X(J)
        Y4=Y(J)
        A3=X4-X3
        M3=(Y4-Y3)/A3
        IF(LM2.EQ.0)      GO TO 43
        IF(J.EQ.2)        GO TO 41
        X2=X(J-2)
        Y2=Y(J-2)
        A2=X3-X2
        M2=(Y3-Y2)/A2
        IF(J.EQ.L0)       GO TO 42
   41   X5=X(J+1)
        Y5=Y(J+1)
        A4=X5-X4
        M4=(Y5-Y4)/A4
        IF(J.EQ.2)        M2=M3+M3-M4
        GO TO 45
   42   M4=M3+M3-M2
        GO TO 45
   43   M2=M3
        M4=M3
   45   IF(J.LE.3)        GO TO 46
        A1=X2-X(J-3)
        M1=(Y2-Y(J-3))/A1
        GO TO 47
   46   M1=M2+M2-M3
   47   IF(J.GE.LM1)      GO TO 48
        A5=X(J+2)-X5
        M5=(Y(J+2)-Y5)/A5
        GO TO 50
   48   M5=M4+M4-M3
! NUMERICAL DIFFERENTIATION
   50   IF(I.EQ.LP1)      GO TO 52
        W2=ABS(M4-M3)
        W3=ABS(M2-M1)
        SW=W2+W3
        IF(SW.NE.0.0)     GO TO 51
        W2=0.5
        W3=0.5
        SW=1.0
   51   T3=(W2*M2+W3*M3)/SW
        IF(I.EQ.1)        GO TO 54
   52   W3=ABS(M5-M4)
        W4=ABS(M3-M2)
        SW=W3+W4
        IF(SW.NE.0.0)     GO TO 53
        W3=0.5
        W4=0.5
        SW=1.0
   53   T4=(W3*M3+W4*M4)/SW
        IF(I.NE.LP1)      GO TO 60
        T3=T4
        SA=A2+A3
        T4=0.5*(M4+M5-A2*(A2-A3)*(M2-M3)/(SA*SA))
        X3=X4
        Y3=Y4
        A3=A2
        M3=M4
        GO TO 60
   54   T4=T3
        SA=A3+A4
        T3=0.5*(M1+M2-A4*(A3-A4)*(M3-M4)/(SA*SA))
        X3=X3-A4
        Y3=Y3-M2*A4
        A3=A4
        M3=M2
! DETERMINATION OF THE COEFFICIENTS
   60   Q2=(2.0*(M3-T3)+M3-T4)/A3
        Q3=(-M3-M3+T3+T4)/(A3*A3)
! COMPUTATION OF THE POLYNOMIAL
   70   DX=UK-P0
   80   V(K)=Q0+DX*(Q1+DX*(Q2+DX*Q3))
      RETURN
! ERROR EXIT
   90 ERR = 1
      RETURN
   91 ERR = 2
      RETURN
   95 ERR = 6
      RETURN
   96 ERR = 7
      RETURN
   99 ERR = 10
      RETURN
      END
