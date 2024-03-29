module var

use decomp_2d
USE variables
USE param

! define all major arrays here

real(mytype), save, allocatable, dimension(:,:,:) :: ux1, ux2, ux3, po3, pp3, dv3, ppx
real(mytype), save, allocatable, dimension(:,:,:) :: uy1, uy2, uy3
real(mytype), save, allocatable, dimension(:,:,:) :: uz1, uz2, uz3
real(mytype), save, allocatable, dimension(:,:,:) :: phi1, phi2, phi3
real(mytype), save, allocatable, dimension(:,:,:) :: gx1, gy1, gz1, hx1, hy1, hz1, phis1,phiss1
real(mytype), save, allocatable, dimension(:,:,:) :: px1, py1, pz1
real(mytype), save, allocatable, dimension(:,:,:) :: ep1
!real(mytype), save, allocatable, dimension(:,:,:) :: epsfull,epsthi,epsaxi

!arrays for statistic collection
real(mytype), save, allocatable, dimension(:,:,:) :: umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
real(mytype), save, allocatable, dimension(:,:,:) :: phimean, phiphimean,pmean
!real(mytype), save, allocatable, dimension(:,:,:) :: epsfullm,epsthim,epsaxim

!real(mytype), save, allocatable, dimension(:,:,:) :: dudxm,dudym,dudzm
!real(mytype), save, allocatable, dimension(:,:,:) :: dvdxm,dvdym,dvdzm
!real(mytype), save, allocatable, dimension(:,:,:) :: dwdxm,dwdym,dwdzm

!real(mytype), save, allocatable, dimension(:,:,:) :: dudx2m,dudy2m,dudz2m
!real(mytype), save, allocatable, dimension(:,:,:) :: dvdx2m,dvdy2m,dvdz2m
!real(mytype), save, allocatable, dimension(:,:,:) :: dwdx2m,dwdy2m,dwdz2m

!real(mytype), save, allocatable, dimension(:,:,:) :: dudydvdxm,dudzdwdxm,dvdzdwdym

!arrays for visualization
real(mytype), save, allocatable, dimension(:,:,:) :: uvisu

!arrays for probe collection
real(mytype), save, allocatable, dimension(:,:,:) :: uprobe
real(mytype), save, allocatable, dimension(:,:,:,:) :: uprobe_write,vprobe_write,wprobe_write,phiprobe_write

! define all work arrays here
real(mytype), save, allocatable, dimension(:,:,:) :: ta1,tb1,tc1,td1,&
     te1,tf1,tg1,th1,ti1,di1,tempgrad
real(mytype), save, allocatable, dimension(:,:,:) :: ta2,tb2,tc2,td2,&
     te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype), save, allocatable, dimension(:,:,:) :: ta3,tb3,tc3,td3,&
     te3,tf3,tg3,th3,ti3,di3

! 
integer, save :: nxmsize, nymsize, nzmsize 


contains

  
  subroutine init_variables

    TYPE(DECOMP_INFO), save :: ph  ! decomposition object

    if (nclx==0) then
       nxmsize = xsize(1)
    else
       nxmsize = xsize(1) -1
    endif
    if (ncly==0) then
       nymsize = ysize(2)
    else
       nymsize = ysize(2) -1
    endif
    if (nclz==0) then
       nzmsize = zsize(3)
    else
       nzmsize = zsize(3) -1
    endif
    call decomp_info_init(nxmsize, nymsize, nzmsize, ph)
    

!X PENCILS
!    call alloc_x0(ux1);call alloc_x0(uy1);call alloc_x0(uz1)
    allocate (ux1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    allocate (uy1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    allocate (uz1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3))) 
    allocate (px1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    allocate (py1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    allocate (pz1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    allocate (phi1(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
    allocate (ppx(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))    
!!$    allocate (epsfull(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
!!$    allocate (epsthi(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
!!$    allocate (epsaxi(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))

    call alloc_x0(gx1);call alloc_x0(gy1);call alloc_x0(gz1);call alloc_x0(phis1) 
    call alloc_x0(hx1);call alloc_x0(hy1);call alloc_x0(hz1);call alloc_x0(phiss1)
    call alloc_x0(ta1);call alloc_x0(tb1);call alloc_x0(tc1)
    call alloc_x0(td1);call alloc_x0(te1);call alloc_x0(tf1)
    call alloc_x0(tg1);call alloc_x0(th1);call alloc_x0(ti1)
    call alloc_x0(di1);call alloc_x0(ep1);call alloc_x0(tempgrad)
    allocate(sx(xsize(2),xsize(3)),vx(xsize(2),xsize(3)))
    !inflow/ouflow 2d arrays
    allocate(bxx1(xsize(2),xsize(3)),bxy1(xsize(2),xsize(3)))
    allocate(bxz1(xsize(2),xsize(3)),bxxn(xsize(2),xsize(3)))
    allocate(bxyn(xsize(2),xsize(3)),bxzn(xsize(2),xsize(3)))
    allocate(bxo(xsize(2),xsize(3)),byo(xsize(2),xsize(3)))
    allocate(bzo(xsize(2),xsize(3)))
    allocate(byx1(xsize(1),xsize(3)),byy1(xsize(1),xsize(3)))
    allocate(byz1(xsize(1),xsize(3)),byxn(xsize(1),xsize(3)))
    allocate(byyn(xsize(1),xsize(3)),byzn(xsize(1),xsize(3)))   
    allocate(bzx1(xsize(1),xsize(2)),bzy1(xsize(1),xsize(2)))
    allocate(bzz1(xsize(1),xsize(2)),bzxn(xsize(1),xsize(2)))
    allocate(bzyn(xsize(1),xsize(2)),bzzn(xsize(1),xsize(2)))
    !pre_correc 2d array
    allocate(dpdyx1(xsize(2),xsize(3)),dpdyxn(xsize(2),xsize(3)))
    allocate(dpdzx1(xsize(2),xsize(3)),dpdzxn(xsize(2),xsize(3)))
    allocate(dpdxy1(xsize(1),xsize(3)),dpdxyn(xsize(1),xsize(3)))
    allocate(dpdzy1(xsize(1),xsize(3)),dpdzyn(xsize(1),xsize(3)))
    allocate(dpdxz1(xsize(1),xsize(2)),dpdxzn(xsize(1),xsize(2)))
    allocate(dpdyz1(xsize(1),xsize(2)),dpdyzn(xsize(1),xsize(2)))

    !arrays for statistic collection!pay attention to the size!
    !allocate (umean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    !allocate (vmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    !allocate (wmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    !allocate (uumean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    !allocate (vvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    !allocate (wwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    !allocate (uvmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    !allocate (uwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    !allocate (vwmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    !allocate (tmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))    
    !allocate (pmean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    if (iscalar==1) then
       allocate (phimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
       allocate (phiphimean(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))
    else
       allocate (phimean(1,1,1))
       allocate (phiphimean(1,1,1))
    endif
!!$    allocate (epsfullm(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
!!$    allocate (epsthim(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
!!$    allocate (epsaxim(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 

    !allocate (dudxm(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dudym(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dudzm(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 

    !allocate (dvdxm(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dvdym(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dvdzm(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 

    !allocate (dwdxm(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dwdym(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dwdzm(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 

    !allocate (dudx2m(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dudy2m(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dudz2m(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 

    !allocate (dvdx2m(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dvdy2m(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dvdz2m(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 

    !allocate (dwdx2m(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dwdy2m(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dwdz2m(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))

    !allocate (dudydvdxm(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dudzdwdxm(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3))) 
    !allocate (dvdzdwdym(xstS(1):xenS(1),xstS(2):xenS(2),xstS(3):xenS(3)))

    !arrays for visualization!pay attention to the size!
    allocate (uvisu(xstV(1):xenV(1),xstV(2):xenV(2),xstV(3):xenV(3)))

!arrays for probe collection!pay attention to the size!
    allocate (uprobe(xstP(1):xenP(1),xstP(2):xenP(2),xstP(3):xenP(3)))
    allocate (uprobe_write(xstP(1):xenP(1),xstP(2):xenP(2),xstP(3):xenP(3),nlength))
    allocate (vprobe_write(xstP(1):xenP(1),xstP(2):xenP(2),xstP(3):xenP(3),nlength))
    allocate (wprobe_write(xstP(1):xenP(1),xstP(2):xenP(2),xstP(3):xenP(3),nlength))
    allocate (phiprobe_write(xstP(1):xenP(1),xstP(2):xenP(2),xstP(3):xenP(3),nlength))
!Y PENCILS
    call alloc_y0(ux2);call alloc_y0(uy2);call alloc_y0(uz2)
    call alloc_y0(ta2);call alloc_y0(tb2);call alloc_y0(tc2)
    call alloc_y0(td2);call alloc_y0(te2);call alloc_y0(tf2)
    call alloc_y0(tg2);call alloc_y0(th2);call alloc_y0(ti2)
    call alloc_y0(tj2)
    call alloc_y0(di2);call alloc_y0(phi2)
    allocate(sy(ysize(1),ysize(3)),vy(ysize(1),ysize(3)))
!Z PENCILS
    call alloc_z0(ux3);call alloc_z0(uy3);call alloc_z0(uz3)
    call alloc_z0(ta3);call alloc_z0(tb3);call alloc_z0(tc3)
    call alloc_z0(td3);call alloc_z0(te3);call alloc_z0(tf3)
    call alloc_z0(tg3);call alloc_z0(th3);call alloc_z0(ti3)
    call alloc_z0(di3);call alloc_z0(phi3)
    allocate(sz(zsize(1),zsize(2)),vz(zsize(1),zsize(2)))

 ! if all periodic
!    allocate (pp3(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
!    allocate (dv3(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
!    allocate (po3(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
    allocate (pp3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
    allocate (dv3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))
    allocate (po3(ph%zst(1):ph%zen(1),ph%zst(2):ph%zen(2),ph%zst(3):ph%zen(3)))


    return
  end subroutine init_variables

  
  ! a set of allocation routines to allocate arrays defined above
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Memory allocation routines to define variables using 2D decomposition
  !    alloc_x? -- define x-pensil variable
  !    alloc_y? -- define y-pensil variable
  !    alloc_z? -- define z-pensil variable
  ! The dimension all local in memory may have additional boundary points.
  !  for example:
  !    alloc_x0 -- size of x dimension is nx
  !    alloc_x1 -- size of x dimension is nx+1
  !    alloc_x2 -- size of x dimension is nx+2
  !    alloc_xm -- size of x dimension is nx-1
  ! The code uses 10 combinations - x0,x1,x2,y0,y1,z0,z1,xm,ym,zm
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  subroutine alloc_x0(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(xsize(1), xsize(2), xsize(3)))
    return
  end subroutine alloc_x0

  subroutine alloc_x1(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(xsize(1)+1, xsize(2), xsize(3)))
    return
  end subroutine alloc_x1

  subroutine alloc_x2(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(xsize(1)+2, xsize(2), xsize(3)))
    return
  end subroutine alloc_x2

  subroutine alloc_xm(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(xsize(1)-1, xsize(2), xsize(3)))
    return
  end subroutine alloc_xm

  subroutine alloc_y0(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(ysize(1), ysize(2), ysize(3)))
    return
  end subroutine alloc_y0

  subroutine alloc_y1(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(ysize(1), ysize(2)+1, ysize(3)))
    return
  end subroutine alloc_y1

  subroutine alloc_ym(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(ysize(1), ysize(2)-1, ysize(3)))
    return
  end subroutine alloc_ym

  subroutine alloc_z0(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(zsize(1), zsize(2), zsize(3)))
    return
  end subroutine alloc_z0

  subroutine alloc_z1(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(zsize(1), zsize(2), zsize(3)+1))
    return
  end subroutine alloc_z1

  subroutine alloc_zm(var)
    implicit none
    real(mytype), allocatable, dimension(:,:,:) :: var
    allocate(var(zsize(1), zsize(2), zsize(3)-1))
    return
  end subroutine alloc_zm

end module var

