module decomp_2d_io

  use decomp_2d
  use MPI

  implicit none

  private        ! Make everything private unless declared public

  public :: decomp_2d_write_one, &
       decomp_2d_write_var, decomp_2d_read_var, &
       decomp_2d_write_plane, decomp_2d_write_every, sub_domain

  ! Generic interface to handle multiple data types
  interface decomp_2d_write_one
     module procedure mpiio_write_real
     module procedure mpiio_write_complex
     module procedure mpiio_write_real_decomp
     module procedure mpiio_write_real_coarse
     module procedure mpiio_write_real_probe
  end interface

  interface decomp_2d_write_var
     module procedure write_var ! write a variable of default global size
     module procedure write_var_decomp !         of specified global size
  end interface

  interface decomp_2d_read_var
     module procedure read_var  ! write a variable of default global size
     module procedure read_var_decomp  !         of specified global size
  end interface

  interface decomp_2d_write_plane
     module procedure write_plane_3d
     module procedure write_plane_2d
  end interface

  interface decomp_2d_write_every
     module procedure write_every
  end interface
  
contains
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Using MPI-IO library to write a single 3D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine mpiio_write_real(nx,ny,nz,ipencil,var,filename)
    
    implicit none
    
    integer, intent(IN) :: nx,ny,nz
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh
    
    sizes(1) = nx
    sizes(2) = ny
    sizes(3) = nz
    
    if (ipencil == 1) then
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
    else if (ipencil == 2) then
       subsizes(1) = ysize(1)
       subsizes(2) = ysize(2)
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = ystart(2)-1
       starts(3) = ystart(3)-1
    else if (ipencil == 3) then
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = zsize(3)
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = zstart(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine mpiio_write_real


  subroutine mpiio_write_real_decomp(decomp,ipencil,var,filename)
    
    implicit none
    
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh
    
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine mpiio_write_real_decomp


  subroutine mpiio_write_real_coarse(ipencil,var,filename,icoarse)
    
    USE param
    USE variables

    implicit none
    
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    integer, intent(IN) :: icoarse !(nstat=1; nvisu=2)
    real(mytype), dimension(:,:,:), intent(IN) :: var
  
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh

    if (icoarse==1) then
       sizes(1) = xszS(1)
       sizes(2) = yszS(2)
       sizes(3) = zszS(3)
       
       if (ipencil == 1) then
          subsizes(1) = xszS(1)
          subsizes(2) = xszS(2)
          subsizes(3) = xszS(3)
          starts(1) = xstS(1)-1  ! 0-based index
          starts(2) = xstS(2)-1
          starts(3) = xstS(3)-1
       else if (ipencil == 2) then
          subsizes(1) = yszS(1)
          subsizes(2) = yszS(2)
          subsizes(3) = yszS(3)
          starts(1) = ystS(1)-1
          starts(2) = ystS(2)-1
          starts(3) = ystS(3)-1
       else if (ipencil == 3) then
          subsizes(1) = zszS(1)
          subsizes(2) = zszS(2)
          subsizes(3) = zszS(3)
          starts(1) = zstS(1)-1
          starts(2) = zstS(2)-1
          starts(3) = zstS(3)-1
       endif
    endif

    if (icoarse==2) then
       sizes(1) = xszV(1)
       sizes(2) = yszV(2)
       sizes(3) = zszV(3)
       
       if (ipencil == 1) then
          subsizes(1) = xszV(1)
          subsizes(2) = xszV(2)
          subsizes(3) = xszV(3)
          starts(1) = xstV(1)-1  ! 0-based index
          starts(2) = xstV(2)-1
          starts(3) = xstV(3)-1
       else if (ipencil == 2) then
          subsizes(1) = yszV(1)
          subsizes(2) = yszV(2)
          subsizes(3) = yszV(3)
          starts(1) = ystV(1)-1
          starts(2) = ystV(2)-1
          starts(3) = ystV(3)-1
       else if (ipencil == 3) then
          subsizes(1) = zszV(1)
          subsizes(2) = zszV(2)
          subsizes(3) = zszV(3)
          starts(1) = zstV(1)-1
          starts(2) = zstV(2)-1
          starts(3) = zstV(3)-1
       endif
    endif

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    
    return
  end subroutine mpiio_write_real_coarse

  subroutine mpiio_write_real_probe(ipencil,var,filename)
    
    USE param
    USE variables

    implicit none
    
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:,:), intent(IN) :: var
  
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(4) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh

    
    sizes(1) = xszP(1)
    sizes(2) = yszP(2)
    sizes(3) = zszP(3)
    sizes(4) = nlength
    if (ipencil == 1) then
       subsizes(1) = xszP(1)
       subsizes(2) = xszP(2)
       subsizes(3) = xszP(3)
       subsizes(4) = nlength
       starts(1) = xstP(1)-1  ! 0-based index
       starts(2) = xstP(2)-1
       starts(3) = xstP(3)-1
       starts(4) = 0
    else if (ipencil == 2) then
       subsizes(1) = yszP(1)
       subsizes(2) = yszP(2)
       subsizes(3) = yszP(3)
       starts(1) = ystP(1)-1
       starts(2) = ystP(2)-1
       starts(3) = ystP(3)-1
    else if (ipencil == 3) then
       subsizes(1) = zszP(1)
       subsizes(2) = zszP(2)
       subsizes(3) = zszP(3)
       starts(1) = zstP(1)-1
       starts(2) = zstP(2)-1
       starts(3) = zstP(3)-1
    endif
!   print *,nrank,starts(1),starts(2),starts(3),starts(4)
    call MPI_TYPE_CREATE_SUBARRAY(4, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3)*subsizes(4), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    
    return
  end subroutine mpiio_write_real_probe


  subroutine mpiio_write_complex(nx,ny,nz,ipencil,var,filename)
    
    implicit none
    
    integer, intent(IN) :: nx,ny,nz
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh
    
    sizes(1) = nx
    sizes(2) = ny
    sizes(3) = nz
    
    if (ipencil == 1) then
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
    else if (ipencil == 2) then
       subsizes(1) = ysize(1)
       subsizes(2) = ysize(2)
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = ystart(2)-1
       starts(3) = ystart(3)-1
    else if (ipencil == 3) then
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = zsize(3)
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = zstart(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, complex_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         complex_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine mpiio_write_complex


  subroutine mpiio_write_complex_decomp(decomp,ipencil,var,filename)
    
    implicit none
    
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    complex(mytype), dimension(:,:,:), intent(IN) :: var
    character(len=*) :: filename
    
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh
    
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, complex_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         complex_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    return
  end subroutine mpiio_write_complex_decomp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 3D array to a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next writing.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_var(fh,disp,nx,ny,nz,ipencil,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: nx,ny,nz       ! global size
    integer, intent(IN) :: ipencil        ! pencil orientation
    real(mytype), dimension(:,:,:), &
         intent(IN) :: var                ! distributed 3D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype

    ! Create file type and set file view
    sizes(1) = nx
    sizes(2) = ny
    sizes(3) = nz
    if (ipencil == 1) then
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
    else if (ipencil == 2) then
       subsizes(1) = ysize(1)
       subsizes(2) = ysize(2)
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = ystart(2)-1
       starts(3) = ystart(3)-1
    else if (ipencil == 3) then
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = zsize(3)
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = zstart(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + nx*ny*nz*mytype

    return
  end subroutine write_var


  subroutine write_var_decomp(fh,disp,decomp,ipencil,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: ipencil        ! pencil orientation
    real(mytype), dimension(:,:,:), &
         intent(IN) :: var                ! distributed 3D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype

    ! Create file type and set file view
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype

    return
  end subroutine write_var_decomp


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Read a 3D array from a big MPI-IO file, starting from 
  !  displacement 'disp'; 'disp' will be updated at end to prepare
  !  next reading.
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine read_var(fh,disp,nx,ny,nz,ipencil,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    integer, intent(IN) :: nx,ny,nz       ! global size
    integer, intent(IN) :: ipencil        ! pencil orientation
    real(mytype), dimension(:,:,:), &
         intent(INOUT) :: var                ! distributed 3D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype

    ! Create file type and set file view
    sizes(1) = nx
    sizes(2) = ny
    sizes(3) = nz
    if (ipencil == 1) then
       subsizes(1) = xsize(1)
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = xstart(1)-1  ! 0-based index
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
    else if (ipencil == 2) then
       subsizes(1) = ysize(1)
       subsizes(2) = ysize(2)
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = ystart(2)-1
       starts(3) = ystart(3)-1
    else if (ipencil == 3) then
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = zsize(3)
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = zstart(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)

    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + nx*ny*nz*mytype

    return
  end subroutine read_var


  subroutine read_var_decomp(fh,disp,decomp,ipencil,var)

    implicit none

    integer, intent(IN) :: fh             ! file handle
    TYPE(DECOMP_INFO), intent(IN) :: decomp
    integer, intent(IN) :: ipencil        ! pencil orientation
    real(mytype), dimension(:,:,:), &
         intent(INOUT) :: var                ! distributed 3D data array
    integer(KIND=MPI_OFFSET_KIND), &
         intent(INOUT) :: disp            ! displacement

    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype

    ! Create file type and set file view
    sizes(1) = decomp%xsz(1)
    sizes(2) = decomp%ysz(2)
    sizes(3) = decomp%zsz(3)
    if (ipencil == 1) then
       subsizes(1) = decomp%xsz(1)
       subsizes(2) = decomp%xsz(2)
       subsizes(3) = decomp%xsz(3)
       starts(1) = decomp%xst(1)-1  ! 0-based index
       starts(2) = decomp%xst(2)-1
       starts(3) = decomp%xst(3)-1
    else if (ipencil == 2) then
       subsizes(1) = decomp%ysz(1)
       subsizes(2) = decomp%ysz(2)
       subsizes(3) = decomp%ysz(3)
       starts(1) = decomp%yst(1)-1
       starts(2) = decomp%yst(2)-1
       starts(3) = decomp%yst(3)-1
    else if (ipencil == 3) then
       subsizes(1) = decomp%zsz(1)
       subsizes(2) = decomp%zsz(2)
       subsizes(3) = decomp%zsz(3)
       starts(1) = decomp%zst(1)-1
       starts(2) = decomp%zst(2)-1
       starts(3) = decomp%zst(3)-1
    endif
    
    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_READ_ALL(fh, var, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)

    call MPI_TYPE_FREE(newtype,ierror)

    ! update displacement for next write
    disp = disp + sizes(1)*sizes(2)*sizes(3)*mytype

    return
  end subroutine read_var_decomp

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 2D slice of the 3D data to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_plane_3d(ipencil,var,iplane,n,filename)

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iplane !(x-plane=1; y-plane=2; z-plane=3)
    integer, intent(IN) :: n ! which plane to write (global coordinate)
    character(len=*) :: filename

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    real(mytype), allocatable, dimension(:,:,:) :: wk2d

    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh

    ! if needed transpose 3D data so that all mpi rank participate I/O
    if (iplane==1) then
       allocate(wk(xsize(1),xsize(2),xsize(3)))
       if (ipencil==1) then
          wk = var
       else if (ipencil==2) then
          call transpose_y_to_x(var,wk)
       else if (ipencil==3) then
          allocate(wk2(ysize(1),ysize(2),ysize(3)))
          call transpose_z_to_y(var,wk2)
          call transpose_y_to_x(wk2,wk)
          deallocate(wk2)
       end if
       allocate(wk2d(1,xsize(2),xsize(3)))
       do k=1,xsize(3)
          do j=1,xsize(2)
             wk2d(1,j,k)=wk(n,j,k)
          end do
       end do
       sizes(1) = 1
       sizes(2) = ysize(2)
       sizes(3) = zsize(3)
       subsizes(1) = 1
       subsizes(2) = xsize(2)
       subsizes(3) = xsize(3)
       starts(1) = 0
       starts(2) = xstart(2)-1
       starts(3) = xstart(3)-1
    else if (iplane==2) then
       allocate(wk(ysize(1),ysize(2),ysize(3)))
       if (ipencil==1) then
          call transpose_x_to_y(var,wk)
       else if (ipencil==2) then
          wk = var
       else if (ipencil==3) then
          call transpose_z_to_y(var,wk)
       end if
       allocate(wk2d(ysize(1),1,ysize(3)))
       do k=1,ysize(3)
          do i=1,ysize(1)
             wk2d(i,1,k)=wk(i,n,k)
          end do
       end do
       sizes(1) = xsize(1)
       sizes(2) = 1
       sizes(3) = zsize(3)
       subsizes(1) = ysize(1)
       subsizes(2) = 1
       subsizes(3) = ysize(3)
       starts(1) = ystart(1)-1
       starts(2) = 0
       starts(3) = ystart(3)-1
    else if (iplane==3) then
       allocate(wk(zsize(1),zsize(2),zsize(3)))
       if (ipencil==1) then
          allocate(wk2(ysize(1),ysize(2),ysize(3)))
          call transpose_x_to_y(var,wk2)
          call transpose_y_to_z(wk2,wk)
          deallocate(wk2)
       else if (ipencil==2) then
          call transpose_y_to_z(var,wk)
       else if (ipencil==3) then
          wk = var
       end if
       allocate(wk2d(zsize(1),zsize(2),1))
       do j=1,zsize(2)
          do i=1,zsize(1) 
             wk2d(i,j,1)=wk(i,j,n)
          end do
       end do
       sizes(1) = xsize(1)
       sizes(2) = ysize(2)
       sizes(3) = 1
       subsizes(1) = zsize(1)
       subsizes(2) = zsize(2)
       subsizes(3) = 1
       starts(1) = zstart(1)-1
       starts(2) = zstart(2)-1
       starts(3) = 0
    end if

    call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
         MPI_ORDER_FORTRAN, real_type, newtype, ierror)
    call MPI_TYPE_COMMIT(newtype,ierror)
    call MPI_FILE_OPEN(MPI_COMM_WORLD, filename, &
         MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
         fh, ierror)
    filesize = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
    disp = 0_MPI_OFFSET_KIND
    call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
         newtype,'native',MPI_INFO_NULL,ierror)
    call MPI_FILE_WRITE_ALL(fh, wk2d, &
         subsizes(1)*subsizes(2)*subsizes(3), &
         real_type, MPI_STATUS_IGNORE, ierror)
    call MPI_FILE_CLOSE(fh,ierror)
    call MPI_TYPE_FREE(newtype,ierror)
    
    deallocate(wk,wk2d);

    return
  end subroutine write_plane_3d

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write a 2D array to a file
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !************** TO DO ***************
  subroutine write_plane_2d(ipencil,var,filename)
    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:), intent(IN) :: var ! 2D array
    character(len=*) :: filename

    if (ipencil==1) then
       ! var should be defined as var(xsize(2)

    else if (ipencil==2) then

    else if (ipencil==3) then

    end if

    return
  end subroutine write_plane_2d


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write 3D array data for every specified mesh point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine write_every(ipencil,var,iskip,jskip,kskip,filename, from1)

    integer, intent(IN) :: ipencil !(x-pencil=1; y-pencil=2; z-pencil=3)
    real(mytype), dimension(:,:,:), intent(IN) :: var
    integer, intent(IN) :: iskip,jskip,kskip 
    character(len=*) :: filename
    logical, intent(IN) :: from1  ! .true.  - save 1,n+1,2n+1...
                                  ! .false. - save n,2n,3n...

    real(mytype), allocatable, dimension(:,:,:) :: wk, wk2
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, key, color, newcomm

    integer, dimension(3) :: xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip

    ! work out the distribution parameters, which may be different from 
    ! the default distribution used by the decomposition library
    !  For exmample if nx=17 and p_row=4
    !    distribution is: 4 4 4 5

    ! If writing from the 1st element
    !  If saving every 3 points, then 5 points to be saved (17/3)
    !    default distribution would be 1 1 1 2
    !    However, 1st block (1-4) contains the 3rd point
    !             2nd block (5-8) contains the 6th point
    !             3rd block (9-12) contains the 9th and 12th point
    !             4th block (13-17) contains then 15th point
    !    giving a 1 1 2 1 distribution
    !    So cannot use the base decomposition library for such IO

    ! If writing from the n-th element (n=?skip)
    !  If saving every 3 points, then 6 points to be saved
    !    However, 1st block (1-4) contains the 1st & 4th point
    !             2nd block (5-8) contains the 7th point
    !             3rd block (9-12) contains the 10th point
    !             4th block (13-17) contains then 12th & 15th point
    !    giving a 1 2 2 1 distribution

    skip(1)=iskip
    skip(2)=jskip
    skip(3)=kskip

    do i=1,3
       if (from1) then
          xst(i) = (xstart(i)+skip(i)-1)/skip(i)
          if (mod(xstart(i)+skip(i)-1,skip(i))/=0) xst(i)=xst(i)+1
          xen(i) = (xend(i)+skip(i)-1)/skip(i)
       else
          xst(i) = xstart(i)/skip(i)
          if (mod(xstart(i),skip(i))/=0) xst(i)=xst(i)+1
          xen(i) = xend(i)/skip(i)
       end if
       xsz(i) = xen(i)-xst(i)+1
    end do
       
    do i=1,3
       if (from1) then
          yst(i) = (ystart(i)+skip(i)-1)/skip(i)
          if (mod(ystart(i)+skip(i)-1,skip(i))/=0) yst(i)=yst(i)+1
          yen(i) = (yend(i)+skip(i)-1)/skip(i)
       else
          yst(i) = ystart(i)/skip(i)
          if (mod(ystart(i),skip(i))/=0) yst(i)=yst(i)+1
          yen(i) = yend(i)/skip(i)
       end if
       ysz(i) = yen(i)-yst(i)+1
    end do

    do i=1,3
       if (from1) then
          zst(i) = (zstart(i)+skip(i)-1)/skip(i)
          if (mod(zstart(i)+skip(i)-1,skip(i))/=0) zst(i)=zst(i)+1
          zen(i) = (zend(i)+skip(i)-1)/skip(i)
       else
          zst(i) = zstart(i)/skip(i)
          if (mod(zstart(i),skip(i))/=0) zst(i)=zst(i)+1
          zen(i) = zend(i)/skip(i)
       end if
       zsz(i) = zen(i)-zst(i)+1
    end do

    ! if 'skip' value is large it is possible that some ranks do not 
    ! contain any points to be written. Subarray constructor requires 
    ! nonzero size so it is not possible to use MPI_COMM_WORLD for IO.
    ! Create a sub communicator for this...
    color = 1
    key = 0  ! rank order doesn't matter
    if (ipencil==1) then
       if (xsz(1)==0 .or. xsz(2)==0 .or. xsz(3)==0) then
          color = 2
       end if
    else if (ipencil==2) then
       if (ysz(1)==0 .or. ysz(2)==0 .or. ysz(3)==0) then
          color = 2
       end if
    else if (ipencil==3) then
       if (zsz(1)==0 .or. zsz(2)==0 .or. zsz(3)==0) then
          color = 2
       end if
    end if
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,newcomm,ierror)

    if (color==1) then ! only ranks in this group do IO collectively
       
       ! generate subarray information
       sizes(1) = xsz(1)
       sizes(2) = ysz(2)
       sizes(3) = zsz(3)
       if (ipencil==1) then
          subsizes(1) = xsz(1)
          subsizes(2) = xsz(2)
          subsizes(3) = xsz(3)
          starts(1) = xst(1)-1
          starts(2) = xst(2)-1
          starts(3) = xst(3)-1
       else if (ipencil==2) then
          subsizes(1) = ysz(1)
          subsizes(2) = ysz(2)
          subsizes(3) = ysz(3)
          starts(1) = yst(1)-1
          starts(2) = yst(2)-1
          starts(3) = yst(3)-1
       else if (ipencil==3) then
          subsizes(1) = zsz(1)
          subsizes(2) = zsz(2)
          subsizes(3) = zsz(3)
          starts(1) = zst(1)-1
          starts(2) = zst(2)-1
          starts(3) = zst(3)-1
       end if
       
       ! copy data from original array
       ! needs a copy of original array in global coordinate 
       if (ipencil==1) then
          allocate(wk(xst(1):xen(1),xst(2):xen(2),xst(3):xen(3)))
          allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
          wk2=var
          if (from1) then
             do k=xst(3),xen(3)
                do j=xst(2),xen(2)
                   do i=xst(1),xen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=xst(3),xen(3)
                do j=xst(2),xen(2)
                   do i=xst(1),xen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if   
       else if (ipencil==2) then
          allocate(wk(yst(1):yen(1),yst(2):yen(2),yst(3):yen(3)))
          allocate(wk2(ystart(1):yend(1),ystart(2):yend(2),ystart(3):yend(3)))
          wk2=var
          if (from1) then
             do k=yst(3),yen(3)
                do j=yst(2),yen(2)
                   do i=yst(1),yen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=yst(3),yen(3)
                do j=yst(2),yen(2)
                   do i=yst(1),yen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if
       else if (ipencil==3) then
          allocate(wk(zst(1):zen(1),zst(2):zen(2),zst(3):zen(3)))
          allocate(wk2(zstart(1):zend(1),zstart(2):zend(2),zstart(3):zend(3)))
          wk2=var
          if (from1) then
             do k=zst(3),zen(3)
                do j=zst(2),zen(2)
                   do i=zst(1),zen(1)
                      wk(i,j,k) = wk2((i-1)*iskip+1,(j-1)*jskip+1,(k-1)*kskip+1)
                   end do
                end do
             end do
          else
             do k=zst(3),zen(3)
                do j=zst(2),zen(2)
                   do i=zst(1),zen(1)
                      wk(i,j,k) = wk2(i*iskip,j*jskip,k*kskip)
                   end do
                end do
             end do
          end if
       end if
       deallocate(wk2)

       ! MPI-IO
       call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
            MPI_ORDER_FORTRAN, real_type, newtype, ierror)
       call MPI_TYPE_COMMIT(newtype,ierror)
       call MPI_FILE_OPEN(newcomm, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
            newtype,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_WRITE_ALL(fh, wk, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            real_type, MPI_STATUS_IGNORE, ierror)
       call MPI_FILE_CLOSE(fh,ierror)
       call MPI_TYPE_FREE(newtype,ierror)

       deallocate(wk)

    end if ! color==1

    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    
    return
  end subroutine write_every
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Write 3D array data for every specified mesh point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  subroutine sub_domain(var0,var1,var2,var3,filename)

    real(mytype), dimension(:,:,:), intent(IN) :: var0, var1, var2, var3
    character(len=*) :: filename

    real(4), allocatable, dimension(:,:,:) :: wk2
    real(4), allocatable, dimension(:,:,:) :: wk
    integer (kind=MPI_OFFSET_KIND) :: filesize, disp
    integer, dimension(3) :: sizes, subsizes, starts
    integer :: i,j,k, ierror, newtype, fh, key, color, newcomm

    integer, dimension(3) :: xsz,ysz,zsz,xst,yst,zst,xen,yen,zen,skip
    integer, dimension(3) :: ssz,sst,sen ! size of the subdomain
    integer, dimension(1) :: ss2,se2     ! The second collection location$B!!(B
    color = 0
    key = 0  ! rank order doesn't matter    
    ! The current subroutine only work for the X Pencil

  ! the process to get it is 3841/120=32 then 3841/2=1921 =50Lb from this -20Lb (32*2) and same for other 
  ! interval

    ! global locations   
    sst(1) = 1281; sen(1) = 1601 ! (30Lb,40Lb)
    sst(2) = 241; sen(2) = 337 ! (0Lb,3Lb)
    sst(3) = 145; sen(3) = 337 ! (-3Lb,3Lb)

    ss2(1) = 2561; se2(1) = 2881 ! (70Lb,80Lb)

    ! local locations     
    do i=1,3
      ssz(i) = sen(i)-sst(i)+1
    end do

    ! local locations     
    do i=1,1
      xst(i) = 1
      xen(i) = sen(i)-sst(i)+1
      xsz(i) = sen(i)-sst(i)+1
    end do

    ! local locations     
    do i=2,2
       If(xstart(i).LE.sst(i) .And. xend(i).GE.sst(i))        then
          xst(i) = 1
          xen(i) = xend(i)-sst(i)+1 
          xsz(i) = xen (i)-xst(i)+1
          Color = 1
       Else if (xstart(i).GE.sst(i) .And. xend(i).lE.sen(i))  then  
          xst(i) = xstart(i)-sst(i)+1
          xen(i) = xend  (i)-sst(i)+1
          xsz(i) = xen   (i)-xst(i)+1  
          Color = 1
       Else if (xstart(i).LE.sen(i) .And. xend(i).GE. sen(i))  then
          xst(i) = xstart(i)-sst(i)+1
          xen(i) = sen   (i)-sst(i)+1
          xsz(i) = xen   (i)-xst(i)+1
          Color = 1
       Else 
          xst(i) = 0
          xen(i) = 0 
          xsz(i) = 0 
          Color  = 0
       End If  
    end do
    
    do i=3,3
      If (color==1) then
       If(xstart(i).LE.sst(i) .And. xend(i).GE.sst(i))        then
          xst(i) = 1
          xen(i) = xend(i)-sst(i)+1 
          xsz(i) = xen (i)-xst(i)+1
          Color = 2
       Else if (xstart(i).GE.sst(i) .And. xend(i).lE.sen(i))  then  
          xst(i) = xstart(i)-sst(i)+1
          xen(i) = xend  (i)-sst(i)+1
          xsz(i) = xen   (i)-xst(i)+1  
          Color = 2
       Else if (xstart(i).LE.sen(i) .And. xend(i).GE. sen(i))  then 
          xst(i) = xstart(i)-sst(i)+1
          xen(i) = sen   (i)-sst(i)+1
          xsz(i) = xen   (i)-xst(i)+1
          Color = 2
       Else 
          xst(i) = 0
          xen(i) = 0 
          xsz(i) = 0 
          Color  = 0
       End If 
     End If 
    end do    
    ! if 'skip' value is large it is possible that some ranks do not 
    ! contain any points to be written. Subarray constructor requires 
    ! nonzero size so it is not possible to use MPI_COMM_WORLD for IO.
    ! Create a sub communicator for this...
    call MPI_COMM_SPLIT(MPI_COMM_WORLD,color,key,newcomm,ierror)

    if (color==2) then ! only ranks in this group do IO collectively
       !write(*,*) xsz(1), xsz(2), xsz(3), nrank
       !write(*,*) xst(1), xst(2), xst(3), nrank
       !write(*,*) xen(1), xen(2), xen(3)
       !write(*,*)  color, nrank
       ! generate subarray information
       sizes(1) = 8*ssz(1)
       sizes(2) = ssz(2)
       sizes(3) = ssz(3)
          subsizes(1) = 8*xsz(1)
          subsizes(2) = xsz(2)
          subsizes(3) = xsz(3)
          starts(1) = xst(1)-1
          starts(2) = xst(2)-1
          starts(3) = xst(3)-1
       
       ! copy data from original array
       ! needs a copy of original array in global coordinate 
        allocate(wk(xst(1):xen(1)+7*xsz(1),xst(2):xen(2),xst(3):xen(3)))
        allocate(wk2(xstart(1):xend(1),xstart(2):xend(2),xstart(3):xend(3)))
        
        wk2=var0
             
        Do k=xst(3),xen(3)
        Do j=xst(2),xen(2)
        Do i=xst(1),xen(1)
          wk(i,j,k) = wk2(i+sst(1)-1-0*xsz(1),j+sst(2)-1,k+sst(3)-1)
        End do;End do
        End do

        wk2=var1
             
        Do k=xst(3),xen(3)
        Do j=xst(2),xen(2)
        Do i=xst(1)+1*xsz(1),xen(1)+1*xsz(1)
          wk(i,j,k) = wk2(i+sst(1)-1-1*xsz(1),j+sst(2)-1,k+sst(3)-1)
        End do;End do
        End do

        wk2=var2
             
        Do k=xst(3),xen(3)
        Do j=xst(2),xen(2)
        Do i=xst(1)+2*xsz(1),xen(1)+2*xsz(1)
          wk(i,j,k) = wk2(i+sst(1)-1-2*xsz(1),j+sst(2)-1,k+sst(3)-1)
        End do;End do
        End do

        wk2=var3
             
        Do k=xst(3),xen(3)
        Do j=xst(2),xen(2)
        Do i=xst(1)+3*xsz(1),xen(1)+3*xsz(1)
          wk(i,j,k) = wk2(i+sst(1)-1-3*xsz(1),j+sst(2)-1,k+sst(3)-1)
        End do;End do
        End do


        wk2=var0
             
        Do k=xst(3),xen(3)
        Do j=xst(2),xen(2)
        Do i=xst(1)+4*xsz(1),xen(1)+4*xsz(1)
          wk(i,j,k) = wk2(i+ss2(1)-1-4*xsz(1),j+sst(2)-1,k+sst(3)-1)
        End do;End do
        End do

        wk2=var1
        
        Do k=xst(3),xen(3)
        Do j=xst(2),xen(2)
        Do i=xst(1)+5*xsz(1),xen(1)+5*xsz(1)
          wk(i,j,k) = wk2(i+ss2(1)-1-5*xsz(1),j+sst(2)-1,k+sst(3)-1)
        End do;End do
        End do

        wk2=var2
             
        Do k=xst(3),xen(3)
        Do j=xst(2),xen(2)
        Do i=xst(1)+6*xsz(1),xen(1)+6*xsz(1)
          wk(i,j,k) = wk2(i+ss2(1)-1-6*xsz(1),j+sst(2)-1,k+sst(3)-1)
        End do;End do
        End do

        wk2=var3
             
        Do k=xst(3),xen(3)
        Do j=xst(2),xen(2)
        Do i=xst(1)+7*xsz(1),xen(1)+7*xsz(1)
          wk(i,j,k) = wk2(i+ss2(1)-1-7*xsz(1),j+sst(2)-1,k+sst(3)-1)
        End do;End do
        End do

       ! MPI-IO
       call MPI_TYPE_CREATE_SUBARRAY(3, sizes, subsizes, starts,  &
            MPI_ORDER_FORTRAN,  MPI_REAL, newtype, ierror)
       call MPI_TYPE_COMMIT(newtype,ierror)
       call MPI_FILE_OPEN(newcomm, filename, &
            MPI_MODE_CREATE+MPI_MODE_WRONLY, MPI_INFO_NULL, &
            fh, ierror)
       filesize = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_SIZE(fh,filesize,ierror)  ! guarantee overwriting
       disp = 0_MPI_OFFSET_KIND
       call MPI_FILE_SET_VIEW(fh,disp,MPI_INTEGER, &
            newtype,'native',MPI_INFO_NULL,ierror)
       call MPI_FILE_WRITE_ALL(fh, wk, &
            subsizes(1)*subsizes(2)*subsizes(3), &
            MPI_REAL, MPI_STATUS_IGNORE, ierror)
       call MPI_FILE_CLOSE(fh,ierror)
       call MPI_TYPE_FREE(newtype,ierror)

       deallocate(wk,wk2)

    end if ! color==1

    Call MPI_COMM_FREE(newcomm, IERROR)
    call MPI_BARRIER(MPI_COMM_WORLD, ierror)
    
    return
  end subroutine sub_domain
end module decomp_2d_io
