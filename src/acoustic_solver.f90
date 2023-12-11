!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!	This file contains subroutines related to the acoustic wave propagation
!!	Author: Clara Estela Jimenez Tejero. 
!!	License: GNU General Public License v3.0
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine solver(nSS,nxx,nySou,nyRec,S,ny,nx,ctmp,dt,nt,dmodel,shot_synth)

!USE mod_parfile, only: dt,nt,dmodel

implicit none

	!! Forward propagation

	!!!! Simulation Geometry
	integer i,ny,nx,j,k,dn,ix,nt
	integer k1,k2,nSS
	integer l,m,acc
	integer dPML,dny,dnx,dny2,nxt,nyt
	integer nxx(nSS),nySou(nSS),nyRec(nSS)
	integer nxreco(nSS),nyreco(nSS)

	!!!! Simulation parameters
	real ctmp(ny,nx),dmodel,dt
	real S(nt,nSS),shot_synth(nt,nSS)
	integer, allocatable :: nxs(:),nys(:)
	real, allocatable :: coef1(:),coef2(:)
		
	!!!! Variables ecuacion acustica
	real, allocatable :: source(:,:)
	real, allocatable :: p(:,:),p_1(:,:),p_3(:,:)
	real, allocatable :: c(:,:),c2(:,:)
	real, allocatable :: deriv2_p_x(:,:),deriv2_p_y(:,:)
	real, allocatable :: deriv_p_x(:,:),deriv_p_y(:,:)

	!!!! Variables PML
	integer, allocatable :: vyt(:),vxt(:),my(:),mx(:),yu(:),yd(:),xl(:),xr(:) 
	real, allocatable :: factorA(:),factorB(:)
	real, allocatable :: ld_factorA(:,:),ld_factorB(:,:)
	real, allocatable :: rd_factorA(:,:),rd_factorB(:,:)
	real, allocatable :: lu_factorA(:,:),lu_factorB(:,:)
	real, allocatable :: ru_factorA(:,:),ru_factorB(:,:)

	real, allocatable :: x_prop_pml(:,:),y_prop_pml(:,:),corner_prop_pml(:,:)

	!!!! Variables PML y
	real, allocatable :: yu_pml_p_1(:,:),yu_pml_p(:,:),yu_pml_p_3(:,:)
	real, allocatable :: yu_phi_p_y_1(:,:),yu_phi_p_y(:,:),yu_chi_p_y_1(:,:),yu_chi_p_y(:,:)

	real, allocatable :: yd_pml_p_1(:,:),yd_pml_p(:,:),yd_pml_p_3(:,:)
	real, allocatable :: yd_phi_p_y_1(:,:),yd_phi_p_y(:,:),yd_chi_p_y_1(:,:),yd_chi_p_y(:,:)

	!!!! Variables PML x direction
	real, allocatable :: xl_pml_p_1(:,:),xl_pml_p(:,:),xl_pml_p_3(:,:)
	real, allocatable :: xl_phi_p_x_1(:,:),xl_phi_p_x(:,:),xl_chi_p_x_1(:,:),xl_chi_p_x(:,:)

	real, allocatable :: xr_pml_p_1(:,:),xr_pml_p(:,:),xr_pml_p_3(:,:)
	real, allocatable :: xr_phi_p_x_1(:,:),xr_phi_p_x(:,:),xr_chi_p_x_1(:,:),xr_chi_p_x(:,:)
	
	!!!! Variables PML corner up direction
	real, allocatable :: lu_pml_p_1(:,:),lu_pml_p(:,:),lu_pml_p_3(:,:)
	real, allocatable :: lu_phi_p_x_1(:,:),lu_phi_p_x(:,:),lu_chi_p_x_1(:,:),lu_chi_p_x(:,:)
	real, allocatable :: lu_phi_p_y_1(:,:),lu_phi_p_y(:,:),lu_chi_p_y_1(:,:),lu_chi_p_y(:,:)

	real, allocatable :: ru_pml_p_1(:,:),ru_pml_p(:,:),ru_pml_p_3(:,:)
	real, allocatable :: ru_phi_p_x_1(:,:),ru_phi_p_x(:,:),ru_chi_p_x_1(:,:),ru_chi_p_x(:,:)
	real, allocatable :: ru_phi_p_y_1(:,:),ru_phi_p_y(:,:),ru_chi_p_y_1(:,:),ru_chi_p_y(:,:)
	
	!!!! Variables PML corner down direction
	real, allocatable :: ld_pml_p_1(:,:),ld_pml_p(:,:),ld_pml_p_3(:,:)
	real, allocatable :: ld_phi_p_x_1(:,:),ld_phi_p_x(:,:),ld_chi_p_x_1(:,:),ld_chi_p_x(:,:)
	real, allocatable :: ld_phi_p_y_1(:,:),ld_phi_p_y(:,:),ld_chi_p_y_1(:,:),ld_chi_p_y(:,:)

	real, allocatable :: rd_pml_p_1(:,:),rd_pml_p(:,:),rd_pml_p_3(:,:)
	real, allocatable :: rd_phi_p_x_1(:,:),rd_phi_p_x(:,:),rd_chi_p_x_1(:,:),rd_chi_p_x(:,:)
	real, allocatable :: rd_phi_p_y_1(:,:),rd_phi_p_y(:,:),rd_chi_p_y_1(:,:),rd_chi_p_y(:,:)

        character(len=500) :: file_name,file_name2,file_name3
        character(len=1000) :: Strtmp
        character(len=1000) :: Strtmp2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	acc=3
	dPML=20		!! PML layers and points discretization

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!! PML layers and points discretization
	
	dn=dPML+2*acc
	dnx=dn
	dny=dn	
        dny2=3+2*acc

	allocate(nxs(nSS),nys(nSS))
	nxs=dnx+nxx
	nys=dny2+nySou
	nxreco=dnx+nxx
	nyreco=dny2+nyRec	

	nxt=nx+2*dnx
	nyt=dny2+ny+dny

	allocate(vxt(nxt),vyt(nyt),mx(nx),my(ny),yu(dny2),yd(dny),xl(dnx),xr(dnx))
	vxt=0;vyt=0;mx=0;my=0;yu=0;yd=0;xl=0;xr=0;
	allocate(coef1(acc+1),coef2(acc+1))
	coef1=0.;coef2=0.;	
	allocate(source(nyt,nxt))
	allocate(p_1(nyt,nxt),p(nyt,nxt),p_3(nyt,nxt))
	source=0.;p_1=0.;p=0.;p_3=0.;

	allocate(c(nyt,nxt),c2(nyt,nxt))
	c=0.;c2=0.; 

	allocate(deriv2_p_y(nyt,nxt),deriv2_p_x(nyt,nxt))
	allocate(deriv_p_y(nyt,nxt),deriv_p_x(nyt,nxt))		
	deriv2_p_y=0.;deriv2_p_x=0.;
	deriv_p_y=0.;deriv_p_x=0.;	
			
	!!! Variables PML
	allocate(x_prop_pml(ny,dn),y_prop_pml(dn,nx),corner_prop_pml(dn,dn))
	x_prop_pml=0.;y_prop_pml=0.;corner_prop_pml=0.;
	
	allocate(factorA(dn),factorB(dn))
	allocate(lu_factorA(dn,dn),lu_factorB(dn,dn),ru_factorA(dn,dn),ru_factorB(dn,dn))
	allocate(ld_factorA(dn,dn),ld_factorB(dn,dn),rd_factorA(dn,dn),rd_factorB(dn,dn))
	factorA=0.;factorB=0.;
	lu_factorA=0;lu_factorB=0;ru_factorA=0;ru_factorB=0;
	ld_factorA=0;ld_factorB=0;rd_factorA=0;rd_factorB=0;
	
	!!! Variables PML y	
	allocate(yu_pml_p_1(dn,nx),yu_pml_p(dn,nx),yu_pml_p_3(dn,nx))
	allocate(yu_phi_p_y_1(dn,nx),yu_phi_p_y(dn,nx),yu_chi_p_y_1(dn,nx),yu_chi_p_y(dn,nx))
	yu_pml_p_1=0.;yu_pml_p=0.;yu_pml_p_3=0.;
	yu_phi_p_y_1=0.;yu_phi_p_y=0.;yu_chi_p_y_1=0.;yu_chi_p_y=0.;
	
	allocate(yd_pml_p_1(dn,nx),yd_pml_p(dn,nx),yd_pml_p_3(dn,nx))
	allocate(yd_phi_p_y_1(dn,nx),yd_phi_p_y(dn,nx),yd_chi_p_y_1(dn,nx),yd_chi_p_y(dn,nx))
	yd_pml_p_1=0.;yd_pml_p=0.;yd_pml_p_3=0.;
	yd_phi_p_y_1=0.;yd_phi_p_y=0.;yd_chi_p_y_1=0.;yd_chi_p_y=0.;

	!!! Variables PML x 
	allocate(xl_pml_p_1(ny,dn),xl_pml_p(ny,dn),xl_pml_p_3(ny,dn))
	allocate(xl_phi_p_x_1(ny,dn),xl_phi_p_x(ny,dn),xl_chi_p_x_1(ny,dn),xl_chi_p_x(ny,dn))
	xl_pml_p_1=0.;xl_pml_p=0.;xl_pml_p_3=0.;
	xl_phi_p_x_1=0.;xl_phi_p_x=0.;xl_chi_p_x_1=0.;xl_chi_p_x=0.;
	
	allocate(xr_pml_p_1(ny,dn),xr_pml_p(ny,dn),xr_pml_p_3(ny,dn))
	allocate(xr_phi_p_x_1(ny,dn),xr_phi_p_x(ny,dn),xr_chi_p_x_1(ny,dn),xr_chi_p_x(ny,dn))
	xr_pml_p_1=0.;xr_pml_p=0.;xr_pml_p_3=0.;
	xr_phi_p_x_1=0.;xr_phi_p_x=0.;xr_chi_p_x_1=0.;xr_chi_p_x=0.;
	
	!!! Variables PML corner up
	allocate(lu_pml_p_1(dn,dn),lu_pml_p(dn,dn),lu_pml_p_3(dn,dn))
	allocate(lu_phi_p_x_1(dn,dn),lu_phi_p_x(dn,dn),lu_chi_p_x_1(dn,dn),lu_chi_p_x(dn,dn))
	allocate(lu_phi_p_y_1(dn,dn),lu_phi_p_y(dn,dn),lu_chi_p_y_1(dn,dn),lu_chi_p_y(dn,dn))
	lu_pml_p_1=0.;lu_pml_p=0.;lu_pml_p_3=0.;
	lu_phi_p_x_1=0.;lu_phi_p_x=0.;lu_chi_p_x_1=0.;lu_chi_p_x=0.;
	lu_phi_p_y_1=0.;lu_phi_p_y=0.;lu_chi_p_y_1=0.;lu_chi_p_y=0.;
	
	allocate(ru_pml_p_1(dn,dn),ru_pml_p(dn,dn),ru_pml_p_3(dn,dn))
	allocate(ru_phi_p_x_1(dn,dn),ru_phi_p_x(dn,dn),ru_chi_p_x_1(dn,dn),ru_chi_p_x(dn,dn))
	allocate(ru_phi_p_y_1(dn,dn),ru_phi_p_y(dn,dn),ru_chi_p_y_1(dn,dn),ru_chi_p_y(dn,dn))
	ru_pml_p_1=0.;ru_pml_p=0.;ru_pml_p_3=0.;
	ru_phi_p_x_1=0.;ru_phi_p_x=0.;ru_chi_p_x_1=0.;ru_chi_p_x=0.;
	ru_phi_p_y_1=0.;ru_phi_p_y=0.;ru_chi_p_y_1=0.;ru_chi_p_y=0.;

	!!! Variables PML corner down
	allocate(ld_pml_p_1(dn,dn),ld_pml_p(dn,dn),ld_pml_p_3(dn,dn))
	allocate(ld_phi_p_x_1(dn,dn),ld_phi_p_x(dn,dn),ld_chi_p_x_1(dn,dn),ld_chi_p_x(dn,dn))
	allocate(ld_phi_p_y_1(dn,dn),ld_phi_p_y(dn,dn),ld_chi_p_y_1(dn,dn),ld_chi_p_y(dn,dn))
	ld_pml_p_1=0.;ld_pml_p=0.;ld_pml_p_3=0.;
	ld_phi_p_x_1=0.;ld_phi_p_x=0.;ld_chi_p_x_1=0.;ld_chi_p_x=0.;
	ld_phi_p_y_1=0.;ld_phi_p_y=0.;ld_chi_p_y_1=0.;ld_chi_p_y=0.;
	
	allocate(rd_pml_p_1(dn,dn),rd_pml_p(dn,dn),rd_pml_p_3(dn,dn))
	allocate(rd_phi_p_x_1(dn,dn),rd_phi_p_x(dn,dn),rd_chi_p_x_1(dn,dn),rd_chi_p_x(dn,dn))
	allocate(rd_phi_p_y_1(dn,dn),rd_phi_p_y(dn,dn),rd_chi_p_y_1(dn,dn),rd_chi_p_y(dn,dn))
	rd_pml_p_1=0.;rd_pml_p=0.;rd_pml_p_3=0.;
	rd_phi_p_x_1=0.;rd_phi_p_x=0.;rd_chi_p_x_1=0.;rd_chi_p_x=0.;
	rd_phi_p_y_1=0.;rd_phi_p_y=0.;rd_chi_p_y_1=0.;rd_chi_p_y=0.;

	!!!!!!
	!!!!!! Define matrix pointers
	!!!!!!

	m=0
	do i=1,nxt
		m=m+1
		vxt(m)=i
	enddo
	m=0
	do i=1,nyt
		m=m+1
		vyt(m)=i
	enddo
	
	!!!! Pointers bottom PML 

	m=0
	do i=dnx+1,dnx+nx
		m=m+1
		mx(m)=i
	enddo

	m=0
	do i=dny2+1,dny2+ny
		m=m+1
		my(m)=i
	enddo
	
	!!!! Pointers PML right-y corner 
	m=0
	do i=1,dny2
		m=m+1
		yu(m)=i
	enddo
	
	m=0
	do i=dny2+ny+1,dny2+ny+dny
		m=m+1
		yd(m)=i
	enddo

	!!!! Pointers PML x-y corner

	m=0
	do i=1,dnx
		m=m+1
		xl(m)=i
	enddo

	m=0
	do i=dnx+nx+1,dnx+nx+dnx
		m=m+1
		xr(m)=i
	enddo
	
	!!!!!	Define Acquisition Geometry
	!!!!!

	c(dny2+1:dny2+ny,dnx+1:dnx+nx)=ctmp;
	do i=1,dny2
	c(i,dnx+1:dnx+nx)=ctmp(1,:)		!! y up
	do j=1,dnx
		c(i,j)=ctmp(1,1)		!! corner left up
		c(i,dnx+nx+j)=ctmp(1,nx) 	!! corner right up
	enddo
	enddo

	do i=dny2+ny+1,dny2+ny+dny
	c(i,dn+1:dnx+nx)=ctmp(ny,:);		
	do j=1,dn
		c(i,j)=ctmp(ny,1);		
		c(i,dnx+nx+j)=ctmp(ny,nx);	
	enddo
	enddo
	
	do i=1,dnx
		c(dny2+1:dny2+ny,i)=ctmp(:,1);		!!layer xl
		c(dny2+1:dny2+ny,dnx+nx+i)=ctmp(:,nx);	!!layer xr
	enddo
	
	c2=c**2*dt**2/dmodel**2;
	
!!!!!!!!!!!!!!   PML'S ABSORTION COEFFICIENTS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	call PML_coeffs(acc,dPML,dt,dmodel,dn,factorA,factorB,&
	ld_factorA,rd_factorA,lu_factorA,ru_factorA,ld_factorB,rd_factorB,lu_factorB,ru_factorB)
	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	!!!!! Start Computation !!!!!

	l=1

	call finite_difference_coef(1,acc,coef1)			
	call finite_difference_coef(2,acc,coef2)			

	do k=2,nt

	l=l+1

	do j=1,nSS	
		source(nys(j),nxs(j))=S(k-1,j)
	enddo

        call space_deriv2x(acc,coef2,nyt,nxt,p,deriv2_p_x)
        call space_deriv2y(acc,coef2,nyt,nxt,p,deriv2_p_y)  
   
!!!		START PML'S CALCULATION
        
	call space_derivx(acc,coef1,nyt,nxt,p,deriv_p_x)
	call space_derivy(acc,coef1,nyt,nxt,p,deriv_p_y)
		
!!!		START PML'S CALCULATION

!!!		PML x- CALCULATION

	call pml_matrices_x(1,acc,coef1,ny,dn,nyt,nxt,my,xr,factorA,factorB,&
	deriv_p_x,deriv2_p_x,xr_phi_p_x,xr_phi_p_x_1,xr_chi_p_x,xr_chi_p_x_1,x_prop_pml)		
	xr_pml_p_3=-xr_pml_p_1+2*xr_pml_p+c2(my,xr)*&
	(deriv2_p_x(my,xr)+deriv2_p_y(my,xr)+x_prop_pml)
		
	call pml_matrices_x(2,acc,coef1,ny,dn,nyt,nxt,my,xl,factorA,factorB,&
	deriv_p_x,deriv2_p_x,xl_phi_p_x,xl_phi_p_x_1,xl_chi_p_x,xl_chi_p_x_1,x_prop_pml)		
	xl_pml_p_3=-xl_pml_p_1+2*xl_pml_p+c2(my,xl)*&
	(deriv2_p_x(my,xl)+deriv2_p_y(my,xl)+x_prop_pml)		

!!!		PML y- CALCULATION
	
	call pml_matrices_y(1,acc,coef1,dn,nx,nyt,nxt,yd,mx,factorA,factorB,&
	deriv_p_y,deriv2_p_y,yd_phi_p_y,yd_phi_p_y_1,yd_chi_p_y,yd_chi_p_y_1,y_prop_pml)		
	yd_pml_p_3=-yd_pml_p_1+2*yd_pml_p+(deriv2_p_x(yd,mx)+deriv2_p_y(yd,mx)+y_prop_pml)*c2(yd,mx)
		
!!!	PML corner down CALCULATION
			
	call pml_matrices_xycorner(1,acc,coef1,dn,dn,nyt,nxt,yd,xl,ld_factorA,ld_factorB,&
	deriv_p_x,deriv2_p_x,ld_phi_p_x,ld_phi_p_x_1,ld_chi_p_x,ld_chi_p_x_1,corner_prop_pml)
	ld_pml_p_3=-ld_pml_p_1+2*ld_pml_p+(deriv2_p_x(yd,xl)+deriv2_p_y(yd,xl)+corner_prop_pml)*c2(yd,xl)
		
	call pml_matrices_xycorner(2,acc,coef1,dn,dn,nyt,nxt,yd,xl,ld_factorA,ld_factorB,&
	deriv_p_y,deriv2_p_y,ld_phi_p_y,ld_phi_p_y_1,ld_chi_p_y,ld_chi_p_y_1,corner_prop_pml)
	ld_pml_p_3=ld_pml_p_3+c2(yd,xl)*corner_prop_pml				

	call pml_matrices_xycorner(1,acc,coef1,dn,dn,nyt,nxt,yd,xr,rd_factorA,rd_factorB,&
	deriv_p_x,deriv2_p_x,rd_phi_p_x,rd_phi_p_x_1,rd_chi_p_x,rd_chi_p_x_1,corner_prop_pml)
	rd_pml_p_3=-rd_pml_p_1+2*rd_pml_p+(deriv2_p_x(yd,xr)+deriv2_p_y(yd,xr)+corner_prop_pml)*c2(yd,xr)
		
	call pml_matrices_xycorner(2,acc,coef1,dn,dn,nyt,nxt,yd,xr,rd_factorA,rd_factorB,&
	deriv_p_y,deriv2_p_y,rd_phi_p_y,rd_phi_p_y_1,rd_chi_p_y,rd_chi_p_y_1,corner_prop_pml)
	rd_pml_p_3=rd_pml_p_3+c2(yd,xr)*corner_prop_pml

	
        p_3(1:dny2,:)=0.

	yu_pml_p_3(1:dny2,:)=0.
	lu_pml_p_3(1:dny2,:)=0.
	ru_pml_p_3(1:dny2,:)=0.	


	p_3(my,mx)=-p_1(my,mx)+p(my,mx)*2.+&
	(deriv2_p_x(my,mx)+deriv2_p_y(my,mx)+source(my,mx))*c2(my,mx)

	!!!!! Add PML layers in propagator model

	p_3(yu,mx)=yu_pml_p_3                   !!capa "y" arriba
	p_3(yu,xl)=lu_pml_p_3            	!!esquina "x-y" izda arriba
	p_3(yu,xr)=ru_pml_p_3            	!!esquina "x-y" derecha arriba

	p_3(yd,mx)=yd_pml_p_3                  !!capa "y" abajo
	p_3(yd,xl)=ld_pml_p_3            	!!esquina "x-y" izda abajo
	p_3(yd,xr)=rd_pml_p_3            	!!esquina "x-y" derecha abajo

	p_3(my,xl)=xl_pml_p_3                	!!capa "x" izda
	p_3(my,xr)=xr_pml_p_3                	!!capa "x" derecha

	!!!!! Time actualization
	
	p_1=p
	p=p_3

	yu_pml_p_1=yu_pml_p
	yu_pml_p=yu_pml_p_3
		
	yd_pml_p_1=yd_pml_p
	yd_pml_p=yd_pml_p_3

	xl_pml_p_1=xl_pml_p
	xl_pml_p=xl_pml_p_3
		
	xr_pml_p_1=xr_pml_p
	xr_pml_p=xr_pml_p_3

	lu_pml_p_1=lu_pml_p
	lu_pml_p=lu_pml_p_3
	
	ru_pml_p_1=ru_pml_p
	ru_pml_p=ru_pml_p_3

	ld_pml_p_1=ld_pml_p
	ld_pml_p=ld_pml_p_3
						
	rd_pml_p_1=rd_pml_p
	rd_pml_p=rd_pml_p_3
			
	yu_phi_p_y_1=yu_phi_p_y
	yu_chi_p_y_1=yu_chi_p_y
		
	yd_phi_p_y_1=yd_phi_p_y
	yd_chi_p_y_1=yd_chi_p_y

	xl_phi_p_x_1=xl_phi_p_x		
	xl_chi_p_x_1=xl_chi_p_x
	xr_phi_p_x_1=xr_phi_p_x		
	xr_chi_p_x_1=xr_chi_p_x

	lu_phi_p_x_1=lu_phi_p_x		
	lu_phi_p_y_1=lu_phi_p_y
	lu_chi_p_x_1=lu_chi_p_x
	lu_chi_p_y_1=lu_chi_p_y
		
	ru_phi_p_x_1=ru_phi_p_x		
	ru_phi_p_y_1=ru_phi_p_y
	ru_chi_p_x_1=ru_chi_p_x
	ru_chi_p_y_1=ru_chi_p_y
						
	ld_phi_p_x_1=ld_phi_p_x		
	ld_phi_p_y_1=ld_phi_p_y
	ld_chi_p_x_1=ld_chi_p_x
	ld_chi_p_y_1=ld_chi_p_y
		
	rd_phi_p_x_1=rd_phi_p_x		
	rd_phi_p_y_1=rd_phi_p_y
	rd_chi_p_x_1=rd_chi_p_x
	rd_chi_p_y_1=rd_chi_p_y
		
	do i=1,nSS
		shot_synth(k-1,i)=p_3(nyreco(i),nxreco(i))
	enddo

	enddo

	shot_synth(nt,:)=shot_synth(nt-1,:)
	
	!!!! Variables propagador
	deallocate(coef1,coef2)
	deallocate(vyt,vxt,my,mx,yd,xl,xr)
	deallocate(c,c2)
	deallocate(deriv2_p_x,deriv2_p_y)
	deallocate(deriv_p_x,deriv_p_y)
  	deallocate(source)
	deallocate(p,p_1,p_3)

	!!!! Variables PML 
	deallocate(factorA,factorB)
	deallocate(lu_factorA,lu_factorB,ru_factorA,ru_factorB)
	deallocate(ld_factorA,ld_factorB,rd_factorA,rd_factorB)
	deallocate(x_prop_pml,y_prop_pml,corner_prop_pml)
	
	!!!! Variables PML y	
	deallocate(yu_pml_p_1,yu_pml_p,yu_pml_p_3)
	deallocate(yu_phi_p_y_1,yu_phi_p_y,yu_chi_p_y_1,yu_chi_p_y)
	
	deallocate(yd_pml_p_1,yd_pml_p,yd_pml_p_3)
	deallocate(yd_phi_p_y_1,yd_phi_p_y,yd_chi_p_y_1,yd_chi_p_y)

	!!!! Variables PML x direction
	deallocate(xl_pml_p_1,xl_pml_p,xl_pml_p_3)
	deallocate(xl_phi_p_x_1,xl_phi_p_x,xl_chi_p_x_1,xl_chi_p_x)

	deallocate(xr_pml_p_1,xr_pml_p,xr_pml_p_3)
	deallocate(xr_phi_p_x_1,xr_phi_p_x,xr_chi_p_x_1,xr_chi_p_x)
	
	!!!! Variables PML corner up direction
	deallocate(lu_pml_p_1,lu_pml_p,lu_pml_p_3)
	deallocate(lu_phi_p_x_1,lu_phi_p_x,lu_chi_p_x_1,lu_chi_p_x)
	deallocate(lu_phi_p_y_1,lu_phi_p_y,lu_chi_p_y_1,lu_chi_p_y)

	deallocate(ru_pml_p_1,ru_pml_p,ru_pml_p_3)
	deallocate(ru_phi_p_x_1,ru_phi_p_x,ru_chi_p_x_1,ru_chi_p_x)
	deallocate(ru_phi_p_y_1,ru_phi_p_y,ru_chi_p_y_1,ru_chi_p_y)
	
	!!!! Variables PML corner down direction
	deallocate(ld_pml_p_1,ld_pml_p,ld_pml_p_3)
	deallocate(ld_phi_p_x_1,ld_phi_p_x,ld_chi_p_x_1,ld_chi_p_x)
	deallocate(ld_phi_p_y_1,ld_phi_p_y,ld_chi_p_y_1,ld_chi_p_y)

	deallocate(rd_pml_p_1,rd_pml_p,rd_pml_p_3)
	deallocate(rd_phi_p_x_1,rd_phi_p_x,rd_chi_p_x_1,rd_chi_p_x)
	deallocate(rd_phi_p_y_1,rd_phi_p_y,rd_chi_p_y_1,rd_chi_p_y)

	return
	end

subroutine finite_difference_coef(num,accuracy,coef)
implicit none

integer	num,accuracy
real coef(accuracy+1)

	if(num.eq.1)	then

        if(accuracy.eq.3)  then
                coef(1)=0.
                coef(2)=3./4.
                coef(3)=-3./20.
                coef(4)=1./60.  
        endif

	if(accuracy.eq.4)  then
                coef(1)=0.
                coef(2)=4./5.
                coef(3)=-1./5.
                coef(4)=4./105.
                coef(5)=-1./280.
        endif

	!! coef(1)=-3.+3./10.-1./40. 	!!modificacion 
	!! coef(4)=1./80.          	!!modificacion 

	endif


	if(num.eq.2)	then

        if(accuracy.eq.3)  then

                coef(1)=-49./18.
                coef(2)=3./2.
                coef(3)=-3./20.	
                coef(4)=1./90.                  

        endif

	if(accuracy.eq.4)  then
                coef(1)=-205./72.     
                coef(2)=8./5.
                coef(3)=-1./5.
                coef(4)=8./315.     
                coef(5)=-1./560.
        endif


	coef(1)=coef(1)/2.

	endif


return
end
	
	subroutine PML_coeffs(acc,dPML,dt,dmodel,dn,factorA,factorB,&
	ld_factorA,rd_factorA,lu_factorA,ru_factorA,ld_factorB,rd_factorB,lu_factorB,ru_factorB)
	implicit none

	integer dn,dPML,acc,i,j,ix
	real pi,dt,dmodel,Refl,wPML
	real factorA(dn),factorB(dn)
	real lu_factorA(dn,dn),lu_factorB(dn,dn),ru_factorA(dn,dn),ru_factorB(dn,dn)
	real ld_factorA(dn,dn),ld_factorB(dn,dn),rd_factorA(dn,dn),rd_factorB(dn,dn)
	
	real, allocatable :: yPML(:),alpha(:),sigma(:)
	allocate(yPML(dn),alpha(dn),sigma(dn))

	pi=3.14159265			
	Refl=0.0000001
	wPML=(dPML-1)*dmodel

	do i=1+acc,dPML+acc	!!from interior to exterior
		yPML(i)=(i-(acc+1))*dmodel
		alpha(i)=30.*pi*(1.-(yPML(i)/wPML)**1)
		sigma(i)=(-(2+1)*2000*log(Refl)/(2.*wPML))*(yPML(i)/wPML)**2
	enddo	

	do i=1+acc,dPML+acc	!!from interior to exterior
		factorA(i)=exp(-(sigma(i)+alpha(i))*dt)
		factorB(i)=(factorA(i)-1)*sigma(i)/(sigma(i)+alpha(i))
	enddo

	do i=dPML+acc+1,dn	!!from interior to exterior
		factorA(i)=factorA(dPML+acc)
		factorB(i)=factorB(dPML+acc)	
	enddo
		
	do i=1,dn
		do j=1,i
			rd_factorA(j,i)=factorA(i)
			rd_factorB(j,i)=factorB(i)
		enddo
	enddo
	do i=1,dn
		if(i+1.le.dn)then
			do j=i+1,dn
				rd_factorA(j,i)=rd_factorA(i,j)	
				rd_factorB(j,i)=rd_factorB(i,j)	
			enddo
		endif
	enddo
		
	do i=1,dn
		do j=1,dn-i+1
			ld_factorA(j,i)=factorA(dn-i+1)
			ld_factorB(j,i)=factorB(dn-i+1)
		enddo
	enddo
	do ix=1,dn-1
		j=2+(ix-1)
		do i=dn-(ix-1),dn
			ld_factorA(j,i)=factorA(j)
			ld_factorB(j,i)=factorB(j)
		enddo
	enddo

	do j=1,dn
		lu_factorA(j,:)=ld_factorA(dn-j+1,:)		
		lu_factorB(j,:)=ld_factorB(dn-j+1,:)	
		ru_factorA(j,:)=rd_factorA(dn-j+1,:)		
		ru_factorB(j,:)=rd_factorB(dn-j+1,:)	
	enddo
	
	deallocate(yPML,alpha,sigma)
	
	return
	end

	subroutine pml_matrices_x(idir,acc,coef1,nvy,nvx,nyt,nxt,vy,vx,factorA,factorB,&
	deriv,deriv2,phi,phi_1,chi,chi_1,prop_pml)

	implicit none

	integer :: idir,acc,nvy,nvx,nxt,nyt,i,j,vx(nvx),vy(nvy)
	real :: deriv(nyt,nxt),deriv2(nyt,nxt),factorA(nvx),factorB(nvx)
	real :: prop_pml(nvy,nvx),phi_1(nvy,nvx),phi(nvy,nvx),chi_1(nvy,nvx),chi(nvy,nvx)
        real :: coef1(acc+1)
	
	real, allocatable :: deriv_phi(:,:)
	allocate(deriv_phi(nvy,nvx))	

	prop_pml=0.;deriv_phi=0.;

	call space_derivx(acc,coef1,nvy,nvx,phi,deriv_phi)
           
        if(idir.eq.1)	then 
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(i)*phi_1(j,i)+factorB(i)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(i)*chi_1(j,i)+(deriv2(vy(j),vx(i))+ deriv_phi(j,i) )*factorB(i)	
	enddo
	enddo
	endif
	
        if(idir.eq.2)	then 
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(nvx-i+1)*phi_1(j,i)+factorB(nvx-i+1)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(nvx-i+1)*chi_1(j,i)+(deriv2(vy(j),vx(i))+deriv_phi(j,i))*factorB(nvx-i+1)	
	enddo
	enddo
	endif

	prop_pml=deriv_phi+chi
		
	deallocate(deriv_phi)	
				
	return
	end
	
	subroutine pml_matrices_y(idir,acc,coef1,nvy,nvx,nyt,nxt,vy,vx,factorA,factorB,&
	deriv,deriv2,phi,phi_1,chi,chi_1,prop_pml)

	implicit none

	integer :: idir,acc,nvy,nvx,nxt,nyt,i,j,vx(nvx),vy(nvy)
	real :: deriv(nyt,nxt),deriv2(nyt,nxt),factorA(nvy),factorB(nvy)
	real :: prop_pml(nvy,nvx),phi_1(nvy,nvx),phi(nvy,nvx),chi_1(nvy,nvx),chi(nvy,nvx)
        real :: coef1(acc+1)
	
	real, allocatable :: deriv_phi(:,:)
	allocate(deriv_phi(nvy,nvx))	

	prop_pml=0.;deriv_phi=0.;
	
	call space_derivy(acc,coef1,nvy,nvx,phi,deriv_phi)
                       
        if(idir.eq.1)	then 
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(j)*phi_1(j,i)+factorB(j)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(j)*chi_1(j,i)+(deriv2(vy(j),vx(i))+deriv_phi(j,i))*factorB(j)	
	enddo
	enddo
	endif
	
        if(idir.eq.2)	then 
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(nvy-j+1)*phi_1(j,i)+factorB(nvy-j+1)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(nvy-j+1)*chi_1(j,i)+(deriv2(vy(j),vx(i))+deriv_phi(j,i))*factorB(nvy-j+1)	
	enddo
	enddo
	endif
	
	prop_pml=deriv_phi+chi
		
	deallocate(deriv_phi)	
				
	return
	end
	
	subroutine pml_matrices_xycorner(id,acc,coef1,nvy,nvx,nyt,nxt,vy,vx,factorA,factorB,&
	deriv,deriv2,phi,phi_1,chi,chi_1,prop_pml)

	implicit none

	integer :: id,acc,nvy,nvx,nxt,nyt,i,j,vx(nvx),vy(nvy)
	real :: deriv(nyt,nxt),deriv2(nyt,nxt),factorA(nvy,nvx),factorB(nvy,nvx)
	real :: prop_pml(nvy,nvx),phi_1(nvy,nvx),phi(nvy,nvx),chi_1(nvy,nvx),chi(nvy,nvx)
        real :: coef1(acc+1)
	
	real, allocatable :: deriv_phi(:,:)
	allocate(deriv_phi(nvy,nvx))	

	prop_pml=0.;deriv_phi=0.;
	
	if(id.eq.1)call space_derivx(acc,coef1,nvy,nvx,phi,deriv_phi)
	if(id.eq.2)call space_derivy(acc,coef1,nvy,nvx,phi,deriv_phi)
	               
	do i=1,nvx
	do j=1,nvy
		phi(j,i)=factorA(j,i)*phi_1(j,i)+factorB(j,i)*deriv(vy(j),vx(i))
		chi(j,i)=factorA(j,i)*chi_1(j,i)+(deriv2(vy(j),vx(i))+deriv_phi(j,i))*factorB(j,i)	
	enddo
	enddo
	
	prop_pml=deriv_phi+chi
		
	deallocate(deriv_phi)	
				
	return
	end
	
	subroutine space_deriv2x(acc,coef,nyt,nxt,p,deriv2_px)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2,k
	real p(nyt,nxt),deriv2_px(nyt,nxt)
	real coef(acc+1)

	j1=1+acc
	j2=nxt-acc

	do j=j1,j2
		do k=1,nyt
			deriv2_px(k,j)=0.
		enddo

		do i=1,acc+1
		do k=1,nyt
			deriv2_px(k,j)=deriv2_px(k,j)+(p(k,j+(i-1))+p(k,j-(i-1)))*coef(i)
		enddo
		enddo

	enddo

	do j=1,j1-1
	do k=1,nyt
		deriv2_px(k,j)=deriv2_px(k,acc+1)
	enddo
	enddo

	do j=j2+1,nxt
	do k=1,nyt
		deriv2_px(k,j)=deriv2_px(k,nxt-acc)
	enddo
	enddo				

	return
	end
	
	subroutine space_deriv2y(acc,coef,nyt,nxt,p,deriv2_py)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2,k
	real p(nyt,nxt),deriv2_py(nyt,nxt)
	real coef(acc+1)

	j1=1+acc
	j2=nyt-acc

	do k=1,nxt
		do j=j1,j2
			deriv2_py(j,k)=0.
			do i=1,acc+1
				deriv2_py(j,k)=deriv2_py(j,k)+(p(j+(i-1),k)+p(j-(i-1),k))*coef(i)
			enddo
		enddo
	enddo

	do k=1,nxt
		do j=1,j1-1
			deriv2_py(j,k)=deriv2_py(acc+1,k)
		enddo

		do j=j2+1,nyt
			deriv2_py(j,k)=deriv2_py(nyt-acc,k)
		enddo
	enddo

	return
	end

	subroutine space_derivx(acc,coef,nyt,nxt,p,deriv_px)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2
	real p(nyt,nxt),deriv_px(nyt,nxt)
	real coef(acc+1)

	j1=1+acc
	j2=nxt-acc

	do j=j1,j2

		do k=1,nyt
			deriv_px(k,j)=0.
		enddo

		do i=1,acc+1
		do k=1,nyt
			deriv_px(k,j)=deriv_px(k,j)+(p(k,j+(i-1))-p(k,j-(i-1)))*coef(i)
		enddo
		enddo

	enddo

	do j=1,j1-1
	do k=1,nyt
		deriv_px(k,j)=deriv_px(k,acc+1)
	enddo
	enddo

	do j=j2+1,nxt
	do k=1,nyt
		deriv_px(k,j)=deriv_px(k,nxt-acc)
	enddo
	enddo				    

	return
	end
		
	subroutine space_derivy(acc,coef,nyt,nxt,p,deriv_py)

	integer nyt,nxt,acc
	integer i,i1,i2,j,j1,j2
	real p(nyt,nxt),deriv_py(nyt,nxt)
	real coef(acc+1)
	
	j1=1+acc
	j2=nyt-acc

	do k=1,nxt
		do j=j1,j2
			deriv_py(j,k)=0.
			do i=1,acc+1
				deriv_py(j,k)=deriv_py(j,k)+(p(j+(i-1),k)-p(j-(i-1),k))*coef(i)
			enddo
		enddo
	enddo

	do k=1,nxt
		do j=1,j1-1
			deriv_py(j,k)=deriv_py(acc+1,k)
		enddo
	enddo

	do k=1,nxt
		do j=j2+1,nyt
			deriv_py(j,k)=deriv_py(nyt-acc,k)
		enddo
	enddo
	
	return
	end

