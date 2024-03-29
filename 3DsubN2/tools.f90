!*******************************************************************
!
subroutine test_scalar_min_max(phi)
!
!*******************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI


implicit none

integer :: i,j,k
real(mytype) :: phimax,phimin,cfl
real(mytype) :: phimax1,phimin1

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: phi
integer :: code

phimax=-1609.
phimin=1609.
do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)
   if (phi(i,j,k).gt.phimax) phimax=phi(i,j,k)
   if (phi(i,j,k).lt.phimin) phimin=phi(i,j,k)
enddo
enddo
enddo


call MPI_REDUCE(phimax,phimax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(phimin,phimin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
if (nrank==0) then
   print *,'PHI max=',phimax1
   print *,'PHI min=',phimin1
endif


return
end subroutine test_scalar_min_max

!*******************************************************************
!
subroutine test_speed_min_max(ux,uy,uz)
!
!*******************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI


implicit none

integer :: i,j,k
real(mytype) :: uxmax,uymax,uzmax,uxmin,uymin,uzmin,cfl
real(mytype) :: uxmax1,uymax1,uzmax1,uxmin1,uymin1,uzmin1

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux,uy,uz
integer :: code,ierror

uxmax=-1609.
uymax=-1609.
uzmax=-1609.
uxmin=1609.
uymin=1609.
uzmin=1609.
do k=1,xsize(3)
do j=1,xsize(2)
do i=1,xsize(1)
   if (ux(i,j,k).gt.uxmax) uxmax=ux(i,j,k)
   if (uy(i,j,k).gt.uymax) uymax=uy(i,j,k)
   if (uz(i,j,k).gt.uzmax) uzmax=uz(i,j,k)
   if (ux(i,j,k).lt.uxmin) uxmin=ux(i,j,k)
   if (uy(i,j,k).lt.uymin) uymin=uy(i,j,k)
   if (uz(i,j,k).lt.uzmin) uzmin=uz(i,j,k)
enddo
enddo
enddo


call MPI_REDUCE(uxmax,uxmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uymax,uymax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uzmax,uzmax1,1,real_type,MPI_MAX,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uxmin,uxmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uymin,uymin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
call MPI_REDUCE(uzmin,uzmin1,1,real_type,MPI_MIN,0,MPI_COMM_WORLD,code)
if (nrank==0) then
   print *,'U,V,W max=',uxmax1,uymax1,uzmax1
   print *,'U,V,W min=',uxmin1,uymin1,uzmin1

!VERIF ERREUR CODE DYNAMIQUE
if( (abs(uxmin1).gt.500.).or.(uxmax1 == -1609.d0).or.(uymax1 == -1609.d0).or.(uzmax1 == -1609.d0) ) then
   print *, 'PLANTAGE DU CODE MAN (VITESSE)'
   call MPI_ABORT( MPI_COMM_WORLD, code ,ierror) 
endif
endif

return
end subroutine test_speed_min_max

!*******************************************************************
!
subroutine restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,px1,py1,pz1,phis1,&
                   hx1,hy1,hz1,phiss1,phG,irestart)
!
!*******************************************************************

USE decomp_2d
USE decomp_2d_io
USE variables
USE param
USE MPI

implicit none

TYPE(DECOMP_INFO) :: phG
integer :: i,j,k,irestart,nzmsize,fh,ierror
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,ep1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: gx1,gy1,gz1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: hx1,hy1,hz1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: px1,py1,pz1
real(mytype), dimension(xsize(1),xsize(2),xsize(3)) :: phi1, phis1, phiss1
real(mytype), dimension(phG%zst(1):phG%zen(1),phG%zst(2):phG%zen(2),phG%zst(3):phG%zen(3)) :: pp3
integer (kind=MPI_OFFSET_KIND) :: filesize, disp
real(mytype) :: xdt



!!$if(nrank==0) then
!!$   print *, 'iscalar= ', iscalar
!!$   print *, 'nscheme= ', nscheme
!!$   print *, 'irestart= ', irestart
!!$endif
!!$stop

if (iscalar==0) then
   if (nscheme.ne.4) then
      if (irestart==1) then !write
!         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
!              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
!              fh, ierror)
!         filesize = 0_MPI_OFFSET_KIND
!         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
!         disp = 0_MPI_OFFSET_KIND
!         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)
!         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
!         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
!         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
!         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
!         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
!         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
!         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
!         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
!         call decomp_2d_write_var(fh,disp,phG,1,pp3)
!         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve01.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve02.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve03.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve04.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve05.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve06.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve07.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve08.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve09.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
         call MPI_FILE_CLOSE(fh,ierror)
!

         !!DEBUG
!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve10.dat', &
!!$              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         filesize = 0_MPI_OFFSET_KIND
!!$         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
!!$         disp = 0_MPI_OFFSET_KIND
!!$         call decomp_2d_write_var(fh,disp,phG,1,pp3)
!!$         call MPI_FILE_CLOSE(fh,ierror)
         !!

      else
         if (nrank==0) print *,'RESTART'
!         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
!              MPI_MODE_RDONLY, MPI_INFO_NULL, &
!              fh, ierror)
!         disp = 0_MPI_OFFSET_KIND
!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)   
!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
!         call decomp_2d_read_var(fh,disp,phG,1,pp3)
!         call MPI_FILE_CLOSE(fh,ierror)
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve01.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)   
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve02.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve03.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve04.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve05.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve06.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve07.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve08.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve09.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         !!DEBUG
!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve10.dat', &
!!$              MPI_MODE_RDONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         disp = 0_MPI_OFFSET_KIND
!!$         call decomp_2d_read_var(fh,disp,phG,1,pp3)
!!$         call MPI_FILE_CLOSE(fh,ierror)
         !!

      endif
   else !AB3
     if (irestart==1) then !write
!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
!!$              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         filesize = 0_MPI_OFFSET_KIND
!!$         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
!!$         disp = 0_MPI_OFFSET_KIND
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hx1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hy1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pp3)
!!$         call MPI_FILE_CLOSE(fh,ierror)

         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve01.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve02.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve03.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve04.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve05.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve06.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve07.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve08.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve09.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve10.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve11.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve12.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
         call MPI_FILE_CLOSE(fh,ierror)
      else
         if (nrank==0) print *,'RESTART'
!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
!!$              MPI_MODE_RDONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         disp = 0_MPI_OFFSET_KIND
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)   
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hx1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hy1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pp3)
!!$         call MPI_FILE_CLOSE(fh,ierror)

         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve01.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)   
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve02.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve03.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve04.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve05.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve06.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve07.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve08.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve09.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hz1)
         call MPI_FILE_CLOSE(fh,ierror)
!      
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve10.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve11.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve12.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
         call MPI_FILE_CLOSE(fh,ierror)
      endif
   endif
else !SCALAR 
if (nscheme.ne.4) then
      if (irestart==1) then !write

!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
!!$              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         filesize = 0_MPI_OFFSET_KIND
!!$         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
!!$         disp = 0_MPI_OFFSET_KIND
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pp3)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phi1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phis1)
!!$         call MPI_FILE_CLOSE(fh,ierror)


         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve01.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve02.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve03.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve04.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve05.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve06.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve07.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve08.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve09.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
         call MPI_FILE_CLOSE(fh,ierror)
         !

         !!DEBUG
!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve10.dat', &
!!$              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         filesize = 0_MPI_OFFSET_KIND
!!$         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
!!$         disp = 0_MPI_OFFSET_KIND
!!$         !DEBUG
!!$         call decomp_2d_write_var(fh,disp,phG,1,pp3)
!!$
!!$         !call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pp3)
!!$         call MPI_FILE_CLOSE(fh,ierror)
         !call decomp_2d_write_one(nx_global,ny_global,nz_global,3,pp3,'sauve10.dat')
         !
         !!
      
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve11.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phi1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve12.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phis1)         
         call MPI_FILE_CLOSE(fh,ierror)
!
       
      else
         if (nrank==0) print *,'RESTART'
!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
!!$              MPI_MODE_RDONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         disp = 0_MPI_OFFSET_KIND
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)   
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pp3)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phi1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phis1)
!!$         call MPI_FILE_CLOSE(fh,ierror)
         
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve01.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)   
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve02.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve03.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve04.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve05.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve06.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve07.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve08.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve09.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
         call MPI_FILE_CLOSE(fh,ierror)
!


         !!DEBUG
!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve10.dat', &
!!$              MPI_MODE_RDONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         disp = 0_MPI_OFFSET_KIND
!!$         call decomp_2d_read_var(fh,disp,phG,1,pp3)
!!$
!!$!         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pp3)
!!$         call MPI_FILE_CLOSE(fh,ierror)
         !call decomp_2d_read_one(nx_global,ny_global,nz_global,3,pp3,'sauve10.dat')
         !!
   
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve11.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phi1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve12.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phis1)
         call MPI_FILE_CLOSE(fh,ierror)      
 
      endif
   else !SCALAR + AB3
     if (irestart==1) then !write
!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
!!$              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         filesize = 0_MPI_OFFSET_KIND
!!$         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
!!$         disp = 0_MPI_OFFSET_KIND
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hx1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hy1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pp3)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phi1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phis1)
!!$         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phiss1)
!!$         call MPI_FILE_CLOSE(fh,ierror)

         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve01.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve02.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve03.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve04.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve05.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve06.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve07.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve08.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve09.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,hz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve10.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve11.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve12.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve13.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phi1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve14.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phis1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve15.dat', &
              MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
              fh, ierror)
         filesize = 0_MPI_OFFSET_KIND
         call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_write_var(fh,disp,nx_global,ny_global,nz_global,1,phiss1)
         call MPI_FILE_CLOSE(fh,ierror)

      else
         if (nrank==0) print *,'RESTART'
!!$         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve.dat', &
!!$              MPI_MODE_RDONLY, MPI_INFO_NULL, &
!!$              fh, ierror)
!!$         disp = 0_MPI_OFFSET_KIND
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)   
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hx1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hy1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pp3)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phi1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phis1)
!!$         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phiss1)
!!$         call MPI_FILE_CLOSE(fh,ierror)

         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve01.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,ux1)   
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve02.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve03.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,uz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve04.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve05.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve06.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,gz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve07.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hx1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve08.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hy1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve09.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,hz1)
         call MPI_FILE_CLOSE(fh,ierror)
!      
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve10.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,px1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve11.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,py1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve12.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,pz1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve13.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phi1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve14.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phis1)
         call MPI_FILE_CLOSE(fh,ierror)
!
         call MPI_FILE_OPEN(MPI_COMM_WORLD, 'sauve15.dat', &
              MPI_MODE_RDONLY, MPI_INFO_NULL, &
              fh, ierror)
         disp = 0_MPI_OFFSET_KIND
         call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,phiss1)
         call MPI_FILE_CLOSE(fh,ierror)

      endif 
   endif
endif

if (irestart==0) then
! reconstruction of the dp/dx, dp/dy and dp/dz from px1,py1 and pz1
! Temporal scheme (1:AB2, 2: RK3, 3:RK4, 4:AB3)
   if (nscheme==1) xdt=dt
   if (nscheme==2) xdt=(3./4.)*dt + (-5./12.)*dt
   if (nscheme==3) xdt=0.041717869325*dt
   if (nscheme==4) xdt=dt
   
   do k=1,xsize(3)
   do j=1,xsize(2)
      dpdyx1(j,k)=py1(1,j,k)/xdt
      dpdzx1(j,k)=pz1(1,j,k)/xdt
      dpdyxn(j,k)=py1(nx,j,k)/xdt
      dpdzxn(j,k)=pz1(nx,j,k)/xdt
   enddo
   enddo
   
   if (xsize(3)==1) then
      do j=1,xsize(2)
      do i=1,xsize(1)
         dpdxz1(i,j)=px1(i,j,1)/xdt
         dpdyz1(i,j)=py1(i,j,1)/xdt
      enddo
      enddo
   endif
   if (xsize(3)==nz) then
      do j=1,xsize(2)
      do i=1,xsize(1)
         dpdxzn(i,j)=px1(i,j,nz)/xdt
         dpdyzn(i,j)=py1(i,j,nz)/xdt
      enddo
      enddo
   endif

   if (xsize(2)==1) then
      do k=1,xsize(3)
      do i=1,xsize(1)
         dpdxy1(i,k)=px1(i,1,k)/xdt
         dpdzy1(i,k)=pz1(i,1,k)/xdt
      enddo
      enddo
   endif
   if (xsize(2)==ny) then
      do k=1,xsize(3)
      do i=1,xsize(1)
         dpdxyn(i,k)=px1(i,ny,k)/xdt
         dpdzyn(i,k)=pz1(i,ny,k)/xdt
      enddo
      enddo
   endif
   if (nrank==0) print *,'reconstruction pressure gradients done!'
endif

if (irestart==0) then
   if (ivirt==2) then
      call MPI_FILE_OPEN(MPI_COMM_WORLD, 'epsilon.dat', &
           MPI_MODE_RDONLY, MPI_INFO_NULL, &
           fh, ierror)
      disp = 0_MPI_OFFSET_KIND
      call decomp_2d_read_var(fh,disp,nx_global,ny_global,nz_global,1,ep1) 
      call MPI_FILE_CLOSE(fh,ierror)
      if (nrank==0) print *,'read epsilon file done from restart'
   endif
endif

end subroutine restart
 
!
!*******************************************************************
subroutine stretching()
!
!*******************************************************************
!
USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

real(mytype) :: yinf,den,xnum,xcx,den1,den2,den3,den4,xnum1,cst
integer :: j

yinf=-yly/2.
den=2.*beta*yinf
xnum=-yinf-sqrt(pi*pi*beta*beta+yinf*yinf)
alpha=abs(xnum/den)
xcx=1./beta/alpha
if (alpha.ne.0.) then
   if (istret.eq.1) yp(1)=0.
   if (istret.eq.2) yp(1)=0.
   if (istret.eq.1) yeta(1)=0.
   if (istret.eq.2) yeta(1)=-1./2.
   if (istret.eq.3) yp(1)=0.
   if (istret.eq.3) yeta(1)=-1./2.
   do j=2,ny
      if (istret==1) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yeta(j)=(j-1.)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=(j-1.)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yeta(j)=((j-1.)*(1./2./ny)-0.5)
         if ((ncly.eq.1).or.(ncly.eq.2)) yeta(j)=((j-1.)*(1./2./(ny-1.))-0.5)
      endif
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yeta(j)))*(sin(pi*yeta(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yeta(j))+1.
      xnum1=(atan(xnum*tan(pi*yeta(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst-yinf
         if (yeta(j).eq.0.5) yp(j)=0.-yinf
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yeta(j).lt.0.5) yp(j)=xnum1-cst+yly
         if (yeta(j).eq.0.5) yp(j)=0.+yly
         if (yeta(j).gt.0.5) yp(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yeta(j).lt.0.5) yp(j)=(xnum1-cst+yly)*2.
         if (yeta(j).eq.0.5) yp(j)=(0.+yly)*2.
         if (yeta(j).gt.0.5) yp(j)=(xnum1+cst+yly)*2.
      endif
   enddo
endif
if (alpha.eq.0.) then
   yp(1)=-1.e10
   do j=2,ny
      yeta(j)=(j-1.)*(1./ny)
      yp(j)=-beta*cos(pi*yeta(j))/sin(yeta(j)*pi)
   enddo
endif
if (alpha.ne.0.) then
   do j=1,ny
      if (istret==1) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))
      endif
      if (istret==2) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./(ny-1.))-0.5
      endif
      if (istret==3) then
         if (ncly.eq.0) yetai(j)=(j-0.5)*(1./2./ny)-0.5
         if ((ncly.eq.1).or.(ncly.eq.2)) yetai(j)=(j-0.5)*(1./2./(ny-1.))-0.5
      endif
      
      den1=sqrt(alpha*beta+1.)
      xnum=den1/sqrt(alpha/pi)/sqrt(beta)/sqrt(pi)
      den=2.*sqrt(alpha/pi)*sqrt(beta)*pi*sqrt(pi)
      den3=((sin(pi*yetai(j)))*(sin(pi*yetai(j)))/beta/pi)+alpha/pi
      den4=2.*alpha*beta-cos(2.*pi*yetai(j))+1.
      xnum1=(atan(xnum*tan(pi*yetai(j))))*den4/den1/den3/den
      cst=sqrt(beta)*pi/(2.*sqrt(alpha)*sqrt(alpha*beta+1.))
      if (istret==1) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst-yinf
         if (yetai(j).eq.0.5) ypi(j)=0.-yinf
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst-yinf
      endif
      if (istret==2) then
         if (yetai(j).lt.0.5) ypi(j)=xnum1-cst+yly
         if (yetai(j).eq.0.5) ypi(j)=0.+yly
         if (yetai(j).gt.0.5) ypi(j)=xnum1+cst+yly
      endif
      if (istret==3) then
         if (yetai(j).lt.0.5) ypi(j)=(xnum1-cst+yly)*2.
         if (yetai(j).eq.0.5) ypi(j)=(0.+yly)*2.
         if (yetai(j).gt.0.5) ypi(j)=(xnum1+cst+yly)*2.
      endif
   enddo
endif
if (alpha.eq.0.) then
   ypi(1)=-1.e10
   do j=2,ny
      yetai(j)=(j-1.)*(1./ny)
      ypi(j)=-beta*cos(pi*yetai(j))/sin(yetai(j)*pi)
   enddo
endif

do j=1,ny
   ppy(j)=yly*(alpha/pi+(1./pi/beta)*sin(pi*yeta(j))* &
        sin(pi*yeta(j)))
   pp2y(j)=ppy(j)*ppy(j)
   pp4y(j)=(-2./beta*cos(pi*yeta(j))*sin(pi*yeta(j)))
enddo
do j=1,ny
   ppyi(j)=yly*(alpha/pi+(1./pi/beta)*sin(pi*yetai(j))* &
        sin(pi*yetai(j)))
   pp2yi(j)=ppyi(j)*ppyi(j)
   pp4yi(j)=(-2./beta*cos(pi*yetai(j))*sin(pi*yetai(j)))
enddo

end subroutine stretching

!*****************************************************************
!
subroutine inversion5_v1(aaa,eee,spI)
!
!*****************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

! decomposition object for spectral space
TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

complex(mytype),dimension(spI%yst(1):spI%yen(1),ny/2,spI%yst(3):spI%yen(3),5) :: aaa
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(2):spI%yen(2),spI%yst(3):spI%yen(3)) :: eee
integer :: i,j,k,m,mi,jc
integer,dimension(2) :: ja,jb
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

real(mytype) :: tmp1,tmp2,tmp3,tmp4

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,ny/2-2
   do i=1,2
      mi=m+i
      do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1) 
         if (real(aaa(j,m,k,3), kind=mytype).ne.0.) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
         if (aimag(aaa(j,m,k,3)).ne.0.)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
         sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
         eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
              aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
      enddo
      enddo
      do jc=ja(i),jb(i)
      do k=spI%yst(3),spI%yen(3)
      do j=spI%yst(1),spI%yen(1) 
         aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
              aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
      enddo
      enddo
      enddo
   enddo
enddo


do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1) 
   if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=real(aaa(j,ny/2,k,2), kind=mytype)/real(aaa(j,ny/2-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
      tmp2=aimag(aaa(j,ny/2,k,2))/aimag(aaa(j,ny/2-1,k,3))
   else
      tmp2=0.
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,ny/2,k,3), kind=mytype)-tmp1*real(aaa(j,ny/2-1,k,4), kind=mytype),&
        aimag(aaa(j,ny/2,k,3))-tmp2*aimag(aaa(j,ny/2-1,k,4)), kind=mytype)

   if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
      tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
      tmp3=real(eee(j,ny/2,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,ny/2-1,k), kind=mytype)
   else
      tmp1=0.
      tmp3=0.
   endif
   if (abs(aimag(b1(j,k))).gt.epsilon) then
      tmp2=aimag(sr(j,k))/aimag(b1(j,k))
      tmp4=aimag(eee(j,ny/2,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,ny/2-1,k))
   else
      tmp2=0.
      tmp4=0.
   endif
   a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   eee(j,ny/2,k)=cmplx(tmp3,tmp4, kind=mytype)
   
   if (abs(real(aaa(j,ny/2-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,ny/2-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,ny/2-1,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,ny/2-1,k,3))
   else
      tmp2=0.
   endif
   b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,ny/2-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
        aimag(aaa(j,ny/2-1,k,4))*aimag(b1(j,k)), kind=mytype)
   eee(j,ny/2-1,k)=cmplx(real(eee(j,ny/2-1,k))*real(b1(j,k))-real(a1(j,k))*real(eee(j,ny/2,k)),&
        aimag(eee(j,ny/2-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,ny/2,k)), kind=mytype)
enddo
enddo

do i=ny/2-2,1,-1
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1) 
   if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,i,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,i,k,3))
   else
      tmp2=0.
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
   eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
        real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
        aimag(eee(j,i,k))*aimag(sr(j,k))-&
        aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
enddo
enddo
enddo

return

end subroutine inversion5_v1 

!*****************************************************************
!
subroutine inversion5_v2(aaa,eee,spI)
!
!*****************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

! decomposition object for spectral space
TYPE(DECOMP_INFO) :: spI

#ifdef DOUBLE_PREC
  real(mytype), parameter :: epsilon = 1.e-16
#else
  real(mytype), parameter :: epsilon = 1.e-8
#endif

complex(mytype),dimension(spI%yst(1):spI%yen(1),ny,spI%yst(3):spI%yen(3),5) :: aaa
complex(mytype),dimension(spI%yst(1):spI%yen(1),ny,spI%yst(3):spI%yen(3)) :: eee
integer :: i,j,k,m,mi,jc
integer,dimension(2) :: ja,jb
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: sr
complex(mytype),dimension(spI%yst(1):spI%yen(1),spI%yst(3):spI%yen(3)) :: a1,b1

real(mytype) :: tmp1,tmp2,tmp3,tmp4

do i=1,2
   ja(i)=4-i
   jb(i)=5-i
enddo
do m=1,ny-2
   do i=1,2
   mi=m+i
   do k=spI%yst(3),spI%yen(3)
   do j=spI%yst(1),spI%yen(1) 
      if (real(aaa(j,m,k,3), kind=mytype).ne.0.) tmp1=real(aaa(j,mi,k,3-i), kind=mytype)/real(aaa(j,m,k,3), kind=mytype)
      if (aimag(aaa(j,m,k,3)).ne.0.)tmp2=aimag(aaa(j,mi,k,3-i))/aimag(aaa(j,m,k,3))
      eee(j,mi,k)=cmplx(real(eee(j,mi,k), kind=mytype)-tmp1*real(eee(j,m,k), kind=mytype),&
           aimag(eee(j,mi,k))-tmp2*aimag(eee(j,m,k)), kind=mytype)
   enddo
   enddo
   do jc=ja(i),jb(i)
   do k=spI%yst(3),spI%yen(3)
   do j=spI%yst(1),spI%yen(1) 
      aaa(j,mi,k,jc)=cmplx(real(aaa(j,mi,k,jc), kind=mytype)-real(sr(j,k), kind=mytype)*real(aaa(j,m,k,jc+i), kind=mytype),&
           aimag(aaa(j,mi,k,jc))-aimag(sr(j,k))*aimag(aaa(j,m,k,jc+i)), kind=mytype)
   enddo
   enddo
   enddo
enddo
enddo
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1) 
   if (abs(real(aaa(j,ny-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=real(aaa(j,ny,k,2), kind=mytype)/real(aaa(j,ny-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,ny-1,k,3))).gt.epsilon) then
      tmp2=aimag(aaa(j,ny,k,2))/aimag(aaa(j,ny-1,k,3))
   else
      tmp2=0.
   endif
   b1(j,k)=cmplx(real(aaa(j,ny,k,3), kind=mytype)-tmp1*real(aaa(j,ny-1,k,4), kind=mytype),&
        aimag(aaa(j,ny,k,3))-tmp2*aimag(aaa(j,ny-1,k,4)), kind=mytype)
   if (abs(real(b1(j,k), kind=mytype)).gt.epsilon) then
      tmp1=real(sr(j,k), kind=mytype)/real(b1(j,k), kind=mytype)
      tmp3=real(eee(j,ny,k), kind=mytype)/real(b1(j,k), kind=mytype)-tmp1*real(eee(j,ny-1,k), kind=mytype)
   else
      tmp1=0.
      tmp3=0.
   endif
   if (abs(aimag(b1(j,k))).gt.epsilon) then
      tmp2=aimag(sr(j,k))/aimag(b1(j,k))
      tmp4=aimag(eee(j,ny,k))/aimag(b1(j,k))-tmp2*aimag(eee(j,ny-1,k))
   else
      tmp2=0.
      tmp4=0.
   endif
   a1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   eee(j,ny,k)=cmplx(tmp3,tmp4, kind=mytype)
   
   if (abs(real(aaa(j,ny-1,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,ny-1,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,ny-1,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,ny-1,k,3))
   else
      tmp2=0.
   endif
   b1(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,ny-1,k,4), kind=mytype)*real(b1(j,k), kind=mytype),&
        aimag(aaa(j,ny-1,k,4))*aimag(b1(j,k)), kind=mytype)
   eee(j,ny-1,k)=cmplx(real(eee(j,ny-1,k), kind=mytype)*real(b1(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,ny,k), kind=mytype),&
        aimag(eee(j,ny-1,k))*aimag(b1(j,k))-aimag(a1(j,k))*aimag(eee(j,ny,k)), kind=mytype)
enddo
enddo

do i=ny-2,1,-1
do k=spI%yst(3),spI%yen(3)
do j=spI%yst(1),spI%yen(1) 
   if (abs(real(aaa(j,i,k,3), kind=mytype)).gt.epsilon) then
      tmp1=1./real(aaa(j,i,k,3), kind=mytype)
   else
      tmp1=0.
   endif
   if (abs(aimag(aaa(j,i,k,3))).gt.epsilon) then
      tmp2=1./aimag(aaa(j,i,k,3))
   else
      tmp2=0.
   endif
   sr(j,k)=cmplx(tmp1,tmp2, kind=mytype)
   a1(j,k)=cmplx(real(aaa(j,i,k,4), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,4))*aimag(sr(j,k)), kind=mytype)
   b1(j,k)=cmplx(real(aaa(j,i,k,5), kind=mytype)*real(sr(j,k), kind=mytype),&
        aimag(aaa(j,i,k,5))*aimag(sr(j,k)), kind=mytype)
   eee(j,i,k)=cmplx(real(eee(j,i,k), kind=mytype)*real(sr(j,k), kind=mytype)-&
        real(a1(j,k), kind=mytype)*real(eee(j,i+1,k), kind=mytype)-&
        real(b1(j,k), kind=mytype)*real(eee(j,i+2,k), kind=mytype),&
        aimag(eee(j,i,k))*aimag(sr(j,k))-&
        aimag(a1(j,k))*aimag(eee(j,i+1,k))-aimag(b1(j,k))*aimag(eee(j,i+2,k)), kind=mytype)
enddo
enddo
enddo

return

end subroutine inversion5_v2

!********************************************************************
!
subroutine channel (ux)
!
!********************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ux

integer :: j,i,k,code
real(mytype) :: can,ut3,ut,ut4

ut3=0.
do k=1,ysize(3)
do i=1,ysize(1)
   ut=0.
   do j=1,ny-1
      if (istret.ne.0) ut=ut+(yp(j+1)-yp(j))*(ux(i,j+1,k)-0.5*(ux(i,j+1,k)-ux(i,j,k)))
      if (istret.eq.0) ut=ut+(yly/(ny-1))*(ux(i,j+1,k)-0.5*(ux(i,j+1,k)-ux(i,j,k)))
   enddo
   ut=ut/yly
   ut3=ut3+ut
enddo
enddo
ut3=ut3/ysize(1)/ysize(3)

call MPI_ALLREDUCE(ut3,ut4,1,real_type,MPI_SUM,MPI_COMM_WORLD,code)
ut4=ut4/nproc

!can=-2.*xnu*gdt(itr) ! Poisseuille    
can=-(2./3.-ut4) ! constant flow rate

if (nrank==0) print *,nrank,'UT',ut4,can

do k=1,ysize(3)
do i=1,ysize(1)
do j=2,ny-1
   ux(i,j,k)=-can+ux(i,j,k)
enddo
enddo
enddo

return
end subroutine channel



!*****************************************************************
!
subroutine collect_data()
!
!*****************************************************************

USE decomp_2d
USE decomp_2d_poisson
USE variables
USE param
USE var
USE MPI

implicit none

integer :: i,j,imin,ii,code
real(mytype), dimension(2000) :: ta
integer,dimension(2000) :: tai,tbi

if (nrank==0) then
   call random_number(ta)
   do i=1,2000
      tai(i)=int(ta(i)*4700./dt)
   enddo
   do j=1,2000
      imin=999998
      ii=0
      do i=1,2000
         if (tai(i).le.imin) then
            imin=tai(i)
            ii=i
         endif
      enddo
      tbi(j)=imin
      tai(ii)=999999
   enddo
   idata=tbi+120001
   do i=1,2000
      print *,i,idata(i),'VALUE COLLECTION DATA'
   enddo
endif

call MPI_BCAST(idata,2000,MPI_INTEGER,0,MPI_COMM_WORLD,code)

end subroutine collect_data
