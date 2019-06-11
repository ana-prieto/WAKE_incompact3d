!############################################################################
!
subroutine VISU_INSTA (ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 


integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

nvect1=xsize(1)*xsize(2)*xsize(3)
!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)
!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1


if(itime==1) then
   print *,nrank,xstart(2),xend(2),xstart(3),xend(3)
endif


!!$!############################################################################
!!VORTICITY
!di1=0.
!do ijk=1,nvect1
!   di1(ijk,1,1)=sqrt((tf1(ijk,1,1)-th1(ijk,1,1))**2+&
!        (tg1(ijk,1,1)-tc1(ijk,1,1))**2+&
!        (tb1(ijk,1,1)-td1(ijk,1,1))**2)
!enddo
!uvisu=0.
!call fine_to_coarseV(1,di1,uvisu)
!990 format('vort',I3.3)
!990 format('vort', I4.4)
!write(filename, 990) itime/imodulo
!call decomp_2d_write_one(1,uvisu,filename,2)
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!     1,di1,filename)
!!$!############################################################################

!!$!############################################################################
!!$!R=-(1./3.)*(ta1**3+te1**3+ti1**3+2*td1*th1*tc1+2*tb1*tf1*tg1+
!!$!  tb1*tg1*tf1+tc1*td1*th1)-ta1*td1*tb1-ta1*tg1*tc1-
!!$!  td1*te1*tb1-tg1*ti1*tc1-te1*th1*tf1-th1*ti1*tf1
!!$di1=0.
!!$do ijk=1,nvect1
!!$   di1(ijk,1,1)=-(1./3.)*(ta1(ijk,1,1)**3+te1(ijk,1,1)**3+ti1(ijk,1,1)**3+&
!!$        2.*td1(ijk,1,1)*th1(ijk,1,1)*tc1(ijk,1,1)+&
!!$        2.*tb1(ijk,1,1)*tf1(ijk,1,1)*tg1(ijk,1,1)+&
!!$        tb1(ijk,1,1)*tg1(ijk,1,1)*tf1(ijk,1,1)+&
!!$        tc1(ijk,1,1)*td1(ijk,1,1)*th1(ijk,1,1))-&
!!$        ta1(ijk,1,1)*td1(ijk,1,1)*tb1(ijk,1,1)-&
!!$        ta1(ijk,1,1)*tg1(ijk,1,1)*tc1(ijk,1,1)-&
!!$        td1(ijk,1,1)*te1(ijk,1,1)*tb1(ijk,1,1)-&
!!$        tg1(ijk,1,1)*ti1(ijk,1,1)*tc1(ijk,1,1)-&
!!$        te1(ijk,1,1)*th1(ijk,1,1)*tf1(ijk,1,1)-&
!!$        th1(ijk,1,1)*ti1(ijk,1,1)*tf1(ijk,1,1)
!!$enddo
!!$uvisu=0.
!!$call fine_to_coarseV(1,di1,uvisu)
!!$991 format('invar_r',I3.3)
!!$write(filename, 991) itime/imodulo
!!$call decomp_2d_write_one(1,uvisu,filename,2)
!!$!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!!$!     1,di1,filename)
!!$!############################################################################

!!$!############################################################################
!!$!Q=-0.5*(ta1**2+te1**2+ti1**2)-td1*tb1-tg1*tc1-th1*tf1
!!$di1=0.
!!$do ijk=1,nvect1
!!$   di1(ijk,1,1)=-0.5*(ta1(ijk,1,1)**2+te1(ijk,1,1)**2+ti1(ijk,1,1)**2)-&
!!$        td1(ijk,1,1)*tb1(ijk,1,1)-&
!!$        tg1(ijk,1,1)*tc1(ijk,1,1)-&
!!$        th1(ijk,1,1)*tf1(ijk,1,1)
!!$enddo
!!$uvisu=0.
!!$call fine_to_coarseV(1,di1,uvisu)
!!$992 format('Qcrit',I4.4)
!!$write(filename, 992) itime/imodulo
!!$call decomp_2d_write_one(1,uvisu,filename,2)
!!$!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!!$!     1,di1,filename)
!!$!############################################################################

!!$!############################################################################
!!$!VELOCITY
!!$uvisu=0.
!!$call fine_to_coarseV(1,ux1,uvisu)
!!$993 format('ux',I4.4)
!!$      write(filename, 993) itime/imodulo
!!$call decomp_2d_write_one(1,uvisu,filename,2)
!!$
!!$uvisu=0.
!!$call fine_to_coarseV(1,uy1,uvisu)
!!$994 format('uy',I4.4)
!!$      write(filename, 994) itime/imodulo
!!$call decomp_2d_write_one(1,uvisu,filename,2)
!!$
!!$uvisu=0.
!!$call fine_to_coarseV(1,uz1,uvisu)
!!$995 format('uz',I4.4)
!!$      write(filename, 995) itime/imodulo
!!$call decomp_2d_write_one(1,uvisu,filename,2)

!############################################################################

!!$!############################################################################
!!$!PASSIVE SCALAR
!!$if (iscalar==1) then
!!$uvisu=0.
!!$call fine_to_coarseV(1,phi1,uvisu)
!!$996 format('phi',I4.4)
!!$   write(filename, 996) itime/imodulo
!!$   call decomp_2d_write_one(1,uvisu,filename,2)
!!$!   call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!!$!        1,phi1,filename)
!!$endif
!!$!############################################################################




!############################################################################
923 format('sondes_U_DU',I4.4)
922 format('Csondes_U_DU',I4.4)


if (nrank==5808) then
   !CENTER (241,241) 
   if (itime==ifirst) then
      write(filename, 922) nrank
      open(nrank,file=filename,form='unformatted')
      print *,nrank, ' cores to',xstart(2),xend(2),xstart(3),xend(3)
   endif
   write(nrank) (((ux1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uy1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uz1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ta1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tb1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tc1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((td1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((te1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tf1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tg1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((th1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ti1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1)
   if (itime==ilast) close(nrank)  
endif


!(241+8,241)
if (nrank==6000) then
   if (itime==ifirst) then
      write(filename, 923) nrank
      open(nrank,file=filename,form='unformatted')
      print *,nrank, ' cores to',xstart(2),xend(2),xstart(3),xend(3)
   endif
   write(nrank) (((ux1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uy1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uz1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ta1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tb1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tc1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((td1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((te1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tf1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tg1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((th1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ti1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1)
   if (itime==ilast) close(nrank)
endif
!(241,241+8)
if (nrank==5809) then
   if (itime==ifirst) then
      write(filename, 923) nrank
      open(nrank,file=filename,form='unformatted')
      print *,nrank, ' cores to',xstart(2),xend(2),xstart(3),xend(3)
   endif
   write(nrank) (((ux1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((uy1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((uz1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((ta1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((tb1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((tc1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((td1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((te1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((tf1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((tg1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((th1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4),&
        (((ti1(i,j,k),i=1,xsize(1),8),j=1,1),k=4,4)
   if (itime==ilast) close(nrank)
endif

!(241+16,241)
if (nrank==6192) then
   if (itime==ifirst) then
      write(filename, 923) nrank
      open(nrank,file=filename,form='unformatted')
      print *,nrank, ' cores to',xstart(2),xend(2),xstart(3),xend(3)
   endif
   write(nrank) (((ux1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uy1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uz1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ta1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tb1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tc1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((td1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((te1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tf1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tg1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((th1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ti1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1)
   if (itime==ilast) close(nrank)
endif
!(241,241+16)
if (nrank==5811) then
   if (itime==ifirst) then
      write(filename, 923) nrank
      open(nrank,file=filename,form='unformatted')
      print *,nrank, ' cores to',xstart(2),xend(2),xstart(3),xend(3)
   endif
   write(nrank) (((ux1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((uy1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((uz1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((ta1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((tb1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((tc1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((td1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((te1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((tf1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((tg1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((th1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2),&
        (((ti1(i,j,k),i=1,xsize(1),8),j=1,1),k=2,2)
   if (itime==ilast) close(nrank)
endif

!(241+24,241)
if (nrank==6384) then
   if (itime==ifirst) then
      write(filename, 923) nrank
      open(nrank,file=filename,form='unformatted')
      print *,nrank, ' cores to',xstart(2),xend(2),xstart(3),xend(3)
   endif
   write(nrank) (((ux1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uy1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uz1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ta1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tb1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tc1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((td1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((te1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tf1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tg1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((th1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ti1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1)
   if (itime==ilast) close(nrank)
endif
!(241,241+24)
if (nrank==5812) then
   if (itime==ifirst) then
      write(filename, 923) nrank
      open(nrank,file=filename,form='unformatted')
      print *,nrank, ' cores to',xstart(2),xend(2),xstart(3),xend(3)
   endif
   write(nrank) (((ux1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((uy1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((uz1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((ta1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((tb1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((tc1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((td1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((te1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((tf1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((tg1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((th1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5),&
        (((ti1(i,j,k),i=1,xsize(1),8),j=1,1),k=5,5)
   if (itime==ilast) close(nrank)
endif

!(241+32,241)
if (nrank==6576) then
   if (itime==ifirst) then
      write(filename, 923) nrank
      open(nrank,file=filename,form='unformatted')
      print *,nrank, ' cores to',xstart(2),xend(2),xstart(3),xend(3)
   endif
   write(nrank) (((ux1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uy1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((uz1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ta1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tb1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tc1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((td1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((te1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tf1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((tg1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((th1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1),&
        (((ti1(i,j,k),i=1,xsize(1),8),j=1,1),k=1,1)
   if (itime==ilast) close(nrank)
endif
!(241,241+32)
if (nrank==5814) then
   if (itime==ifirst) then
      write(filename, 923) nrank
      open(nrank,file=filename,form='unformatted')
      print *,nrank, ' cores to',xstart(2),xend(2),xstart(3),xend(3)
   endif
   write(nrank) (((ux1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((uy1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((uz1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((ta1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((tb1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((tc1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((td1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((te1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((tf1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((tg1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((th1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3),&
        (((ti1(i,j,k),i=1,xsize(1),8),j=1,1),k=3,3)
   if (itime==ilast) close(nrank)
endif



!############################################################################
!PRESSURE
!IT IS IN A SEPARATE SUBROUTINE
!############################################################################
end subroutine VISU_INSTA


!############################################################################
!
subroutine PROBE(ux1,uy1,uz1,phi1,uprobe,uprobe_write,vprobe_write,wprobe_write,phiprobe_write)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xszP(1),xszP(2),xszP(3)) :: uprobe
real(mytype),dimension(xszP(1),xszP(2),xszP(3),nlength) :: uprobe_write
real(mytype),dimension(xszP(1),xszP(2),xszP(3),nlength) :: vprobe_write
real(mytype),dimension(xszP(1),xszP(2),xszP(3),nlength) :: wprobe_write
real(mytype),dimension(xszP(1),xszP(2),xszP(3),nlength) :: phiprobe_write
character(len=20) :: filename
integer :: i,j,k

if (itime-ifirst+1==1) uprobe_write=123456.

call fine_to_coarseP(1,ux1,uprobe)
do k=1,xszP(3)
do j=1,xszP(2)
do i=1,xszP(1)
   uprobe_write(i,j,k,itime-ifirst+1)=uprobe(i,j,k)
enddo
enddo
enddo
call fine_to_coarseP(1,uy1,uprobe)
do k=1,xszP(3)
do j=1,xszP(2)
do i=1,xszP(1)
   vprobe_write(i,j,k,itime-ifirst+1)=uprobe(i,j,k)
enddo
enddo
enddo
call fine_to_coarseP(1,uz1,uprobe)
do k=1,xszP(3)
do j=1,xszP(2)
do i=1,xszP(1)
   wprobe_write(i,j,k,itime-ifirst+1)=uprobe(i,j,k)
enddo
enddo
enddo
if (iscalar==1) then
   call fine_to_coarseP(1,phi1,uprobe)
   do k=1,xszP(3)
   do j=1,xszP(2)
   do i=1,xszP(1)
      phiprobe_write(i,j,k,itime-ifirst+1)=uprobe(i,j,k)
   enddo
   enddo
   enddo
endif

if (mod(itime,isave)==0) then
   call decomp_2d_write_one(1,uprobe_write,'probe_u')       
   call decomp_2d_write_one(1,vprobe_write,'probe_v')       
   call decomp_2d_write_one(1,wprobe_write,'probe_w')       
   if (iscalar==1) then
      call decomp_2d_write_one(1,phiprobe_write,'probe_phi')     
   endif
endif

end subroutine PROBE


!############################################################################
!
subroutine STATISTIC(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
     uvmean,uwmean,vwmean,phiphimean,tmean)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: phimean, phiphimean
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1

!umean=ux1
call fine_to_coarseS(1,ux1,tmean)
umean(:,:,:)=umean(:,:,:)+tmean(:,:,:)

!vmean=uy1
call fine_to_coarseS(1,uy1,tmean)
vmean(:,:,:)=vmean(:,:,:)+tmean(:,:,:)

!wmean=uz1
call fine_to_coarseS(1,uz1,tmean)
wmean(:,:,:)=wmean(:,:,:)+tmean(:,:,:)

if (iscalar==1) then
   !phimean=phi1
   call fine_to_coarseS(1,phi1,tmean)
   phimean(:,:,:)=phimean(:,:,:)+tmean(:,:,:)
endif

!uumean=ux1*ux1
ta1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uumean(:,:,:)=uumean(:,:,:)+tmean(:,:,:)

!vvmean=uy1*uy1
ta1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vvmean(:,:,:)=vvmean(:,:,:)+tmean(:,:,:)

!wwmean=uz1*uz1
ta1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
wwmean(:,:,:)=wwmean(:,:,:)+tmean(:,:,:)

!uvmean=ux1*uy1
ta1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uvmean(:,:,:)=uvmean(:,:,:)+tmean(:,:,:)

!uwmean=ux1*uz1
ta1(:,:,:)=ux1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uwmean(:,:,:)=uwmean(:,:,:)+tmean(:,:,:)

!vwmean=uy1*uz1
ta1(:,:,:)=uy1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vwmean(:,:,:)=vwmean(:,:,:)+tmean(:,:,:)

if (iscalar==1) then
   !phiphimean=phi1*phi1
   ta1(:,:,:)=phi1(:,:,:)*phi1(:,:,:)
   call fine_to_coarseS(1,ta1,tmean)
   phiphimean(:,:,:)=phiphimean(:,:,:)+tmean(:,:,:)
endif

!for a verification
!call decomp_2d_write_one(nx_global,ny_global,nz_global,&
!           1,ta1,'compa.dat')

if (mod(itime,isave)==0) then

   call decomp_2d_write_one(1,umean,'umean3.dat',1)
   call decomp_2d_write_one(1,vmean,'vmean3.dat',1)
   call decomp_2d_write_one(1,wmean,'wmean3.dat',1)
   call decomp_2d_write_one(1,uumean,'uumean3.dat',1)
   call decomp_2d_write_one(1,vvmean,'vvmean3.dat',1)
   call decomp_2d_write_one(1,wwmean,'wwmean3.dat',1)
   call decomp_2d_write_one(1,uvmean,'uvmean3.dat',1)
   call decomp_2d_write_one(1,uwmean,'uwmean3.dat',1)
   call decomp_2d_write_one(1,vwmean,'vwmean3.dat',1)
   if (nrank==0) print *,'write stat arrays velocity done!'
   if (iscalar==1) then
      call decomp_2d_write_one(1,phimean,'phimean3.dat',1)
      call decomp_2d_write_one(1,phiphimean,'phiphimean3.dat',1)
      if (nrank==0) print *,'write stat arrays scalar done!'
   endif
!   call decomp_2d_write_one(nx_global,ny_global,nz_global,1,ux1,'compa.dat')

endif

end subroutine STATISTIC

!############################################################################
!
subroutine STAT_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
     ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,pmean,tmean)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG,ph2,ph3

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3 
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ta1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1 
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: pmean,tmean

integer :: code,icomplet,nxmsize,nymsize,nzmsize
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

!WORK Z-PENCILS
call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!WORK X-PENCILS
call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
call interi6(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
!The pressure field on the main mesh is in tb1
!pmean
call fine_to_coarseS(1,tb1,tmean)
pmean(:,:,:)=pmean(:,:,:)+tmean(:,:,:)

if (mod(itime,isave)==0) then
   call decomp_2d_write_one(1,pmean,'pmean3.dat',1)
   if (nrank==0) print *,'write stat arrays pressure done!'
endif

end subroutine STAT_PRE

!############################################################################
!
subroutine STAT_EPS (ux1,uy1,uz1,&
     dudxm,dudym,dudzm,&
     dvdxm,dvdym,dvdzm,&
     dwdxm,dwdym,dwdzm,&
     dudx2m,dudy2m,dudz2m,&
     dvdx2m,dvdy2m,dvdz2m,&
     dwdx2m,dwdy2m,dwdz2m,&
     dudydvdxm,dudzdwdxm,dvdzdwdym,tempgrad,&
     tmean,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,tempgrad
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: dudxm,dudym,dudzm,tmean
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: dvdxm,dvdym,dvdzm
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: dwdxm,dwdym,dwdzm
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: dudx2m,dudy2m,dudz2m
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: dvdx2m,dvdy2m,dvdz2m
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: dwdx2m,dwdy2m,dwdz2m
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: dudydvdxm,dudzdwdxm,dvdzdwdym

integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

nvect1=xsize(1)*xsize(2)*xsize(3)
!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)
!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1

!!$!ENERGY DISSIPATION (FULL, THI and AXI HYPOTHESIS)
!!$do k=1,xsize(3)
!!$   do j=1,xsize(2)
!!$      do i=1,xsize(1)
!!$         epsfull(i,j,k)=xnu*( 2.*(ta1(i,j,k)*ta1(i,j,k)+te1(i,j,k)*te1(i,j,k)&
!!$              +ti1(i,j,k)*ti1(i,j,k)) &
!!$              +(td1(i,j,k)*td1(i,j,k)+tb1(i,j,k)*tb1(i,j,k)+tg1(i,j,k)*tg1(i,j,k)&
!!$              +tc1(i,j,k)*tc1(i,j,k)+th1(i,j,k)*th1(i,j,k)+tf1(i,j,k)*tf1(i,j,k))&
!!$              +2.*(td1(i,j,k)*tb1(i,j,k)+tg1(i,j,k)*tc1(i,j,k)+th1(i,j,k)*tf1(i,j,k)))
!!$
!!$         epsthi(i,j,k)=15.*xnu*ta1(i,j,k)*ta1(i,j,k)
!!$         
!!$         epsaxi(i,j,k)=xnu*(-ta1(i,j,k)*ta1(i,j,k)+2.*td1(i,j,k)*td1(i,j,k)&
!!$              +2.*tb1(i,j,k)*tb1(i,j,k)+8.*te1(i,j,k)*te1(i,j,k))
!!$         !xnu*((5./3.)*ta1(i,j,k)*ta1(i,j,k)+2.*tg1(i,j,k)*tg1(i,j,k)&
!!$         !     +2.*tb1(i,j,k)*tb1(i,j,k)+(8./3.)*th1(i,j,k)*th1(i,j,k))
!!$      enddo
!!$   enddo
!!$enddo
!!$!epsmean=epsilon1 
!!$call fine_to_coarseS(1,epsfull,tmean)
!!$epsfullm(:,:,:)=epsfullm(:,:,:)+tmean(:,:,:)
!!$
!!$call fine_to_coarseS(1,epsthi,tmean)
!!$epsthim(:,:,:)=epsthim(:,:,:)+tmean(:,:,:)
!!$
!!$call fine_to_coarseS(1,epsaxi,tmean)
!!$epsaxim(:,:,:)=epsaxim(:,:,:)+tmean(:,:,:)
!!$
!!$call fine_to_coarseS(1,ta1,tmean)
!!$dudxm(:,:,:)=dudxm(:,:,:)+tmean(:,:,:)
!
!!$if (mod(itime,isave)==0) then
!!$   call decomp_2d_write_one(1,epsfullm,'epsmean_full1.dat',1)
!!$   call decomp_2d_write_one(1,epsthim,'epsmean_thi1.dat',1)
!!$   call decomp_2d_write_one(1,epsaxim,'epsmean_axi1.dat',1)
!!$   call decomp_2d_write_one(1,dudxm,'dudxm1.dat',1)
!!$   if (nrank==0) print *,'write stat array epsilon done!'  
!!$endif

!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
! 9 MEAN GRADIENTS
!dudx mean
call fine_to_coarseS(1,ta1,tmean)
dudxm(:,:,:)=dudxm(:,:,:)+tmean(:,:,:)
!dudy mean
call fine_to_coarseS(1,td1,tmean)
dudym(:,:,:)=dudym(:,:,:)+tmean(:,:,:)
!dudz mean
call fine_to_coarseS(1,tg1,tmean)
dudzm(:,:,:)=dudzm(:,:,:)+tmean(:,:,:)

!dvdx mean
call fine_to_coarseS(1,tb1,tmean)
dvdxm(:,:,:)=dvdxm(:,:,:)+tmean(:,:,:)
!dvdy mean
call fine_to_coarseS(1,te1,tmean)
dvdym(:,:,:)=dvdym(:,:,:)+tmean(:,:,:)
!dvdz mean
call fine_to_coarseS(1,th1,tmean)
dvdzm(:,:,:)=dvdzm(:,:,:)+tmean(:,:,:)

!dwdx mean
call fine_to_coarseS(1,tc1,tmean)
dwdxm(:,:,:)=dwdxm(:,:,:)+tmean(:,:,:)
!dwdy mean
call fine_to_coarseS(1,tf1,tmean)
dwdym(:,:,:)=dwdym(:,:,:)+tmean(:,:,:)
!dwdz mean
call fine_to_coarseS(1,ti1,tmean)
dwdzm(:,:,:)=dwdzm(:,:,:)+tmean(:,:,:)

!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
! 9 MEAN SQUARE GRADIENTS
!dudx2 mean
tempgrad(:,:,:)=ta1(:,:,:)*ta1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dudx2m(:,:,:)=dudx2m(:,:,:)+tmean(:,:,:)
!dudy2 mean
tempgrad(:,:,:)=td1(:,:,:)*td1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dudy2m(:,:,:)=dudy2m(:,:,:)+tmean(:,:,:)
!dudz2 mean
tempgrad(:,:,:)=tg1(:,:,:)*tg1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dudz2m(:,:,:)=dudz2m(:,:,:)+tmean(:,:,:)

!dvdx2 mean
tempgrad(:,:,:)=tb1(:,:,:)*tb1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dvdx2m(:,:,:)=dvdx2m(:,:,:)+tmean(:,:,:)
!dvdy2 mean
tempgrad(:,:,:)=te1(:,:,:)*te1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dvdy2m(:,:,:)=dvdy2m(:,:,:)+tmean(:,:,:)
!dvdz2 mean
tempgrad(:,:,:)=th1(:,:,:)*th1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dvdz2m(:,:,:)=dvdz2m(:,:,:)+tmean(:,:,:)

!dwdx2 mean
tempgrad(:,:,:)=tc1(:,:,:)*tc1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dwdx2m(:,:,:)=dwdx2m(:,:,:)+tmean(:,:,:)
!dwdy2 mean
tempgrad(:,:,:)=tf1(:,:,:)*tf1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dwdy2m(:,:,:)=dwdy2m(:,:,:)+tmean(:,:,:)
!dwdz2 mean
tempgrad(:,:,:)=ti1(:,:,:)*ti1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dwdz2m(:,:,:)=dwdz2m(:,:,:)+tmean(:,:,:)

!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
! 3 MEAN CROSS PRODUCTS GRADIENTS
!dudydvdx mean
tempgrad(:,:,:)=td1(:,:,:)*tb1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dudydvdxm(:,:,:)=dudydvdxm(:,:,:)+tmean(:,:,:)
!dudzdwdx mean
tempgrad(:,:,:)=tg1(:,:,:)*tc1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dudzdwdxm(:,:,:)=dudzdwdxm(:,:,:)+tmean(:,:,:)
!dvdzdwdy mean
tempgrad(:,:,:)=th1(:,:,:)*tf1(:,:,:)
call fine_to_coarseS(1,tempgrad,tmean)
dvdzdwdym(:,:,:)=dvdzdwdym(:,:,:)+tmean(:,:,:)

!WRITE FILES
if (mod(itime,isave)==0) then
   !Mean gradients
   call decomp_2d_write_one(1,dudxm,'dudxm3.dat',1)
   call decomp_2d_write_one(1,dudym,'dudym3.dat',1)
   call decomp_2d_write_one(1,dudzm,'dudzm3.dat',1)

   call decomp_2d_write_one(1,dvdxm,'dvdxm3.dat',1)
   call decomp_2d_write_one(1,dvdym,'dvdym3.dat',1)
   call decomp_2d_write_one(1,dvdzm,'dvdzm3.dat',1)

   call decomp_2d_write_one(1,dwdxm,'dwdxm3.dat',1)
   call decomp_2d_write_one(1,dwdym,'dwdym3.dat',1)
   call decomp_2d_write_one(1,dwdzm,'dwdzm3.dat',1)

   !Mean square gradients
   call decomp_2d_write_one(1,dudx2m,'dudx2m3.dat',1)
   call decomp_2d_write_one(1,dudy2m,'dudy2m3.dat',1)
   call decomp_2d_write_one(1,dudz2m,'dudz2m3.dat',1)

   call decomp_2d_write_one(1,dvdx2m,'dvdx2m3.dat',1)
   call decomp_2d_write_one(1,dvdy2m,'dvdy2m3.dat',1)
   call decomp_2d_write_one(1,dvdz2m,'dvdz2m3.dat',1)

   call decomp_2d_write_one(1,dwdx2m,'dwdx2m3.dat',1)
   call decomp_2d_write_one(1,dwdy2m,'dwdy2m3.dat',1)
   call decomp_2d_write_one(1,dwdz2m,'dwdz2m3.dat',1)

   !Mean cross products
   call decomp_2d_write_one(1,dudydvdxm,'dudydvdxm3.dat',1)
   call decomp_2d_write_one(1,dudzdwdxm,'dudzdwdxm3.dat',1)
   call decomp_2d_write_one(1,dvdzdwdym,'dvdzdwdym3.dat',1)

   if (nrank==0) print *,'write stat array gradient done!'  
endif

end subroutine STAT_EPS

!############################################################################
!
subroutine EXTRA_STAT (ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

integer :: code,icomplet
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

nvect1=xsize(1)*xsize(2)*xsize(3)
!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)
!du/dx=ta1, du/dy=td1 and du/dz=tg1
!dv/dx=tb1, dv/dy=te1 and dv/dz=th1
!dw/dx=tc1, dw/dy=tf1 and dw/dz=ti1

!FIRST DERIVATIVE
!call fine_to_coarseS(1,ta1,tmean)
!dudx(:,:,:)=dudx(:,:,:)+tmean(:,:,:)

!THIRD ORDER MOMENTS
!di1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)*ux1(:,:,:)
!call fine_to_coarseS(1,di1,tmean)
!uuu(:,:,:)=uuu(:,:,:)+tmean(:,:,:)

!di1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)*uy1(:,:,:)
!call fine_to_coarseS(1,di1,tmean)
!vvv(:,:,:)=vvv(:,:,:)+tmean(:,:,:)

!di1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)*uz1(:,:,:)
!call fine_to_coarseS(1,di1,tmean)
!www(:,:,:)=www(:,:,:)+tmean(:,:,:)
!

!FOURTH ORDER MOMENTS
!di1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)*ux1(:,:,:)*ux1(:,:,:)
!call fine_to_coarseS(1,di1,tmean)
!uuuu(:,:,:)=uuuu(:,:,:)+tmean(:,:,:)

!di1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)*uy1(:,:,:)*uy1(:,:,:)
!call fine_to_coarseS(1,di1,tmean)
!vvvv(:,:,:)=vvvv(:,:,:)+tmean(:,:,:)

!di1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)*uz1(:,:,:)*uz1(:,:,:)
!call fine_to_coarseS(1,di1,tmean)
!wwww(:,:,:)=wwww(:,:,:)+tmean(:,:,:)
!
end subroutine EXTRA_STAT

!############################################################################
!
subroutine VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
     ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG,ph2,ph3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3 
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ta1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1 

integer :: code,icomplet,nxmsize,nymsize,nzmsize
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filenameP

!WORK Z-PENCILS
call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!WORK X-PENCILS
call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
call interi6(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
!The pressure field on the main mesh is in tb1

!PRESSURE
1901 format('InstP',I7.7)
     write(filenameP,1901) itime
     call decomp_2d_write_plane(1,tb1,1,150,trim(filenameP)//".dat")

   if (nrank==0) then
        write(*,1391) itime
1391    format('Vorticity_Collection at Time step =',i7)
     endif

end subroutine VISU_PRE


!############################################################################
!
subroutine VISU_SLICE  (ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
!
!############################################################################

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

   implicit none
  
   TYPE(DECOMP_INFO) :: phG
   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: Oa1,Ob1,Oc1
   real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
   real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
   real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
   real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

   real(mytype) :: xp,xvis
   integer :: code,icomplet,ivis
   integer :: ijk,nvect1,nvect2,nvect3,i,j,k
   character(len=20) nfichier,nfichier1
   character(len=20) :: filename
   character(len=50) :: filenameU,  filenameV,  filenameW, filenameO
   character(len=50) :: filenameO1, filenameO2, filenameO3, filenameD

   nvect1=xsize(1)*xsize(2)*xsize(3)


!	1911 format('InstU',I7.7)   
!	1912 format('InstV',I7.7)
!	1913 format('InstW',I7.7)

!     write(filenameU,1911) itime
!     write(filenameV,1912) itime
!     write(filenameW,1913) itime
      
!     call decomp_2d_write_plane(1,ux1,1,17,trim(filenameU)//".dat")
!     call decomp_2d_write_plane(1,uy1,1,17,trim(filenameV)//".dat")
!     call decomp_2d_write_plane(1,uz1,1,17,trim(filenameW)//".dat")

  !x-derivatives
  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  !y-derivatives
  call transpose_x_to_y(ux1,td2)
  call transpose_x_to_y(uy1,te2)
  call transpose_x_to_y(uz1,tf2)
  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  !!z-derivatives
  call transpose_y_to_z(td2,td3)
  call transpose_y_to_z(te2,te3)
  call transpose_y_to_z(tf2,tf3)
  call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
  !!all back to x-pencils
  call transpose_z_to_y(ta3,td2)
  call transpose_z_to_y(tb3,te2)
  call transpose_z_to_y(tc3,tf2)
  call transpose_y_to_x(td2,tg1)
  call transpose_y_to_x(te2,th1)
  call transpose_y_to_x(tf2,ti1)
  call transpose_y_to_x(ta2,td1)
  call transpose_y_to_x(tb2,te1)
  call transpose_y_to_x(tc2,tf1)
  !du/dx=ta1 du/dy=td1 and du/dz=tg1
  !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
  !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
  !############################################################################
  !COMPUTE VORTICITY MAGNITUDE IN di1
  !di1=0.

  !Oa1, Omega1; Oa2, Omega2;  Oa3, Omega3
  !############################################################################
  !COMPUTE VORTICITY MAGNITUDE IN di1  
  ! do k=1,xsize(3)
  ! do j=1,xsize(2)
  ! do i=1,xsize(1)
  !   Oa1(i,j,k) = (tf1(i,j,k)-th1(i,j,k))
  !   Ob1(i,j,k) = (tg1(i,j,k)-tc1(i,j,k))
  !   Oc1(i,j,k) = (tb1(i,j,k)-td1(i,j,k))
  ! enddo; enddo
  ! enddo 

!1921 format('InstO1',I7.7)   
!1922 format('InstO2',I7.7)
!1923 format('InstO3',I7.7)

!     write(filenameO1,1921) itime
!     write(filenameO2,1922) itime
!     write(filenameO3,1923) itime

!     call decomp_2d_write_plane(1,Oa1,1,17,trim(filenameO1)//".dat")
!     call decomp_2d_write_plane(1,Ob1,1,17,trim(filenameO2)//".dat")
!    call decomp_2d_write_plane(1,Oc1,1,17,trim(filenameO3)//".dat")

!############################################################################
!COMPUTE VORTICITY MAGNITUDE IN di1
  di1=0.

  !  do k=1,xsize(3)
  !  do j=1,xsize(2)
  !  do i=1,xsize(1)
  !    di1(i,j,k)=sqrt((tf1(i,j,k)-th1(i,j,k))**2+&
  !                    (tg1(i,j,k)-tc1(i,j,k))**2+&
  !                    (tb1(i,j,k)-td1(i,j,k))**2)
  !  enddo; enddo
  !  enddo

  do ijk=1,nvect1
    di1(ijk,1,1)=sqrt((tf1(ijk,1,1)-th1(ijk,1,1))**2+&
         (tg1(ijk,1,1)-tc1(ijk,1,1))**2+&
         (tb1(ijk,1,1)-td1(ijk,1,1))**2)
  enddo

!1941 format('InstO',I7.7)
!     write(filenameO,1941) itime
!     call decomp_2d_write_plane(1,di1,1,17,trim(filenameO)//".dat")

  !Curl of Omega
  !############################################################################
  !COMPUTE VORTICITY MAGNITUDE IN di1  
! x-derivatives
! call derx (ta1,Oa1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
! call derx (tb1,Ob1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!  call derx (tc1,Oc1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  !y-derivatives
!  call transpose_x_to_y(Oa1,td2)
!  call transpose_x_to_y(Ob1,te2)
!  call transpose_x_to_y(Oc1,tf2)
!  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!  call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
!  call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  !!z-derivatives
!  call transpose_y_to_z(td2,td3)
!  call transpose_y_to_z(te2,te3)
!  call transpose_y_to_z(tf2,tf3)
! call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
!  call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
! call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
  !!all back to x-pencils
!  call transpose_z_to_y(ta3,td2)
!  call transpose_z_to_y(tb3,te2)
!  call transpose_z_to_y(tc3,tf2)
!  call transpose_y_to_x(td2,tg1)
!  call transpose_y_to_x(te2,th1)
!  call transpose_y_to_x(tf2,ti1)
! call transpose_y_to_x(ta2,td1)
!  call transpose_y_to_x(tb2,te1)
!  call transpose_y_to_x(tc2,tf1)
  !du/dx=ta1 du/dy=td1 and du/dz=tg1
  !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
  !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1
  !############################################################################
  !COMPUTE VORTICITY MAGNITUDE IN di1
!   di1=0.; Oa1=0.;Ob1=0.;Oc1=0.;

!   do k=1,xsize(3)
!   do j=1,xsize(2)
!   do i=1,xsize(1)
!     Oa1(i,j,k) = (tf1(i,j,k)-th1(i,j,k))
!     Ob1(i,j,k) = (tg1(i,j,k)-tc1(i,j,k))
!     Oc1(i,j,k) = (tb1(i,j,k)-td1(i,j,k))
!   enddo; enddo
!   enddo 

!   do k=1,xsize(3)
!   do j=1,xsize(2)
!   do i=1,xsize(1)
!     di1(i,j,k) = ux1(i,j,k)*Oa1(i,j,k) + uy1(i,j,k)*Ob1(i,j,k) + &
!                 uz1(i,j,k)*Oc1(i,j,k)
!   enddo; enddo
!   enddo
!
!1931 format('InstD',I7.7)
!    write(filenameD,1931) itime
!    call decomp_2d_write_plane(1,di1,1,17,trim(filenameD)//".dat")

     
!     if (nrank==0) then
!        write(*,1291) itime
!1291    format('Vorticity_Collection at Time step =',i7)
  !   endif   
  !!
  !############################################################################
  !PRESSURE
  !IT IS IN A SEPARATE SUBROUTINE
  !############################################################################


993 format('ux_X',I2.2,'_',I4.4)
994 format('uy_X',I2.2,'_',I4.4)
995 format('uz_X',I2.2,'_',I4.4)
931 format ('vort_X',I2.2,'_',I4.4)

    xp=10
    xvis=20.
    do i=1,3
       ivis = int((xvis+xp)/dx)+1
       write(filename,993) int(xvis),itime/imodulo
       call decomp_2d_write_plane(1,ux1,1,ivis,filename)

       write(filename,994) int(xvis), itime/imodulo
       call decomp_2d_write_plane(1,uy1,1,ivis,filename)

       write(filename,995) int(xvis), itime/imodulo
       call decomp_2d_write_plane(1,uz1,1,ivis,filename)

       write(filename,931) int(xvis), itime/imodulo
       call decomp_2d_write_plane(1,di1,1,ivis,filename)
       xvis=xvis+20.
    enddo


!! Zc =0
997 format('ux_Zc',I4.4)
  write(filename,997) itime/imodulo
  call decomp_2d_write_plane(1,ux1,3,nz/2,filename)

998 format('uy_Zc',I4.4)
  write(filename,998) itime/imodulo
  call decomp_2d_write_plane(1,uy1,3,nz/2,filename)

999 format('uz_Zc',I4.4)
  write (filename,999) itime/imodulo
  call decomp_2d_write_plane(1,uz1,3,nz/2,filename)

930 format('vort_Zc',I4.4)
  write(filename,930) itime/imodulo
  call decomp_2d_write_plane(1,di1,3,nz/2,filename)

do ijk=1,nvect1
    di1(ijk,1,1)=di1(ijk,1,1)*di1(ijk,1,1)
enddo

893 format('enstrophy_X',I2.2,'_',I4.4)
    xp=10
    xvis=20.
    do i=1,3
       ivis=int((xvis+xp)/dx)+1
       write(filename,893) int(xvis), itime/imodulo
       call decomp_2d_write_plane(1,di1,1,ivis,filename)
       xvis=xvis+20
    enddo



end subroutine VISU_SLICE

!############################################################################
!
subroutine PRE_SLICE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
     ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG,ph2,ph3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3 
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real(mytype),dimension(nxmsize,xsize(2),xsize(3)) :: ta1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1 

integer :: code,icomplet,nxmsize,nymsize,nzmsize
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename

!WORK Z-PENCILS
call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!WORK X-PENCILS
call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
call interi6(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
!The pressure field on the main mesh is in tb1

!!X=0 (xlx/2) PLANE 
!PRESSURE
990 format('pp_x0',I4.4)
write(filename,990) itime/imodulo
call decomp_2d_write_plane(1,tb1,1,(nx+1)/2,filename)

!!Z=0
991 format('pp_z0',I4.4)
write(filename,991) itime/imodulo
call decomp_2d_write_plane(1,tb1,3,(nz+1)/2,filename)

!!Y=0
992 format ('pp_y0',I4.4)
write(filename,992) itime/imodulo
call decomp_2d_write_plane(1,tb1,2,1,filename)



end subroutine PRE_SLICE

!#############################################################################




subroutine VISU_SUBDOMAIN  (ux1,uy1,uz1,ppx,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
!
!############################################################################

  USE param
  USE variables
  USE decomp_2d
  USE decomp_2d_io

  implicit none

  TYPE(DECOMP_INFO) :: phG
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1,ppx
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: Oa1,Ob1,Oc1
  real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
  real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
  real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
  real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

  real(mytype) :: xp,xvis
  integer :: code,icomplet,ivis
  integer :: ijk,nvect1,nvect2,nvect3,i,j,k
  character(len=20) nfichier,nfichier1
  character(len=20) :: filename
  character(len=50) :: filenameU,  filenameV,  filenameW, filenameO
  character(len=50) :: filenameO1, filenameO2, filenameO3, filenameD

  nvect1=xsize(1)*xsize(2)*xsize(3)

  !x-derivatives
  call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
  call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
  !y-derivatives
  call transpose_x_to_y(ux1,td2)
  call transpose_x_to_y(uy1,te2)
  call transpose_x_to_y(uz1,tf2)
  call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
  call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
  !!z-derivatives
  call transpose_y_to_z(td2,td3)
  call transpose_y_to_z(te2,te3)
  call transpose_y_to_z(tf2,tf3)
  call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
  call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
  !!all back to x-pencils
  call transpose_z_to_y(ta3,td2)
  call transpose_z_to_y(tb3,te2)
  call transpose_z_to_y(tc3,tf2)
  call transpose_y_to_x(td2,tg1)
  call transpose_y_to_x(te2,th1)
  call transpose_y_to_x(tf2,ti1)
  call transpose_y_to_x(ta2,td1)
  call transpose_y_to_x(tb2,te1)
  call transpose_y_to_x(tc2,tf1)
  !du/dx=ta1 du/dy=td1 and du/dz=tg1
  !dv/dx=tb1 dv/dy=te1 and dv/dz=th1
  !dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1



  !!VORTICITY COMPUTATION
  di1=0.

  do ijk=1,nvect1
    di1(ijk,1,1)=sqrt((tf1(ijk,1,1)-th1(ijk,1,1))**2+&
         (tg1(ijk,1,1)-tc1(ijk,1,1))**2+&
         (tb1(ijk,1,1)-td1(ijk,1,1))**2)
  enddo

1912 format('InstN23D_UVWP',I7.7)
      write(filename,1912) itime

     ! call sub_domain(ux1,trim(filenameU)//".dat")
     call sub_domain(ux1,uy1,uz1,di1,trim(filename)//".dat")


     if (nrank==0) then
        write(*,1291) itime
1291    format('Vorticity_Collection at Time step =',i7)
     endif   
  !!
  !############################################################################
  !PRESSURE
  !IT IS IN A SEPARATE SUBROUTINE
  !############################################################################
end subroutine VISU_SUBDOMAIN

!############################################################################
###########################################################################
!
subroutine VISU_STG (pp3,ppx,ta1,tb1,di1,ta2,tb2,di2,&
     ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG,ph2,ph3
real(mytype),dimension(xszV(1),xszV(2),xszV(3)) :: uvisu 

real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),nzmsize) :: pp3 
!Z PENCILS NXM NYM NZM-->NXM NYM NZ
real(mytype),dimension(ph3%zst(1):ph3%zen(1),ph3%zst(2):ph3%zen(2),zsize(3)) :: ta3,di3
!Y PENCILS NXM NYM NZ -->NXM NY NZ
real(mytype),dimension(ph3%yst(1):ph3%yen(1),nymsize,ysize(3)) :: ta2
real(mytype),dimension(ph3%yst(1):ph3%yen(1),ysize(2),ysize(3)) :: tb2,di2
!X PENCILS NXM NY NZ  -->NX NY NZ
real(mytype),dimension(nxmsize,xsize(2),xsize(3))  :: ta1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: tb1,di1, ppx 

integer :: code,icomplet,nxmsize,nymsize,nzmsize
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filenameP

!WORK Z-PENCILS
call interiz6(ta3,pp3,di3,sz,cifip6z,cisip6z,ciwip6z,cifz6,cisz6,ciwz6,&
     (ph3%zen(1)-ph3%zst(1)+1),(ph3%zen(2)-ph3%zst(2)+1),nzmsize,zsize(3),1)
!WORK Y-PENCILS
call transpose_z_to_y(ta3,ta2,ph3) !nxm nym nz
call interiy6(tb2,ta2,di2,sy,cifip6y,cisip6y,ciwip6y,cify6,cisy6,ciwy6,&
     (ph3%yen(1)-ph3%yst(1)+1),nymsize,ysize(2),ysize(3),1)
!WORK X-PENCILS
call transpose_y_to_x(tb2,ta1,ph2) !nxm ny nz
call interi6(tb1,ta1,di1,sx,cifip6,cisip6,ciwip6,cifx6,cisx6,ciwx6,&
     nxmsize,xsize(1),xsize(2),xsize(3),1)
!The pressure field on the main mesh is in tb1
   ppx(:,:,:) = tb1(:,:,:) 
     
!PRESSURE
!1901 format('InstP',I7.7)
!     write(filenameP,1901) itime
!     call decomp_2d_write_plane(1,tb1,1,17,trim(filenameP)//".dat")

!   if (nrank==0) then
!        write(*,1391) itime
!1391    format('Vorticity_Collection at Time step =',i7)
!     endif

end subroutine VISU_STG


!############################################################################
!
subroutine STAT_EPS_SLICE (ux1,uy1,uz1,epsfull,epsthi,epsaxi,epsfullm,epsthim,epsaxim,tmean,&
     ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG)
!
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: epsfull,epsthi,epsaxi
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: epsfullm,epsthim,epsaxim,tmean

integer :: code,icomplet,ivis
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename
real(mytype) :: xvis,xp

!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)
!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1


!ENERGY DISSIPATION (FULL, THI and AXI HYPOTHESIS)
do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         epsfull(i,j,k)= xnu*( 2.*(ta1(i,j,k)*ta1(i,j,k)+te1(i,j,k)*te1(i,j,k)&
              +ti1(i,j,k)*ti1(i,j,k)) &
              +td1(i,j,k)*td1(i,j,k)+tb1(i,j,k)*tb1(i,j,k)+tg1(i,j,k)*tg1(i,j,k)&
              +tc1(i,j,k)*tc1(i,j,k)+th1(i,j,k)*th1(i,j,k)+tf1(i,j,k)*tf1(i,j,k)&
              +2.*(td1(i,j,k)*tb1(i,j,k)+tg1(i,j,k)*tc1(i,j,k)+th1(i,j,k)*tf1(i,j,k)))

         epsthi(i,j,k)=15.*xnu*ta1(i,j,k)*ta1(i,j,k)
         
         epsaxi(i,j,k)= xnu*(-ta1(i,j,k)*ta1(i,j,k)+2.*td1(i,j,k)*td1(i,j,k)&
              +2.*tb1(i,j,k)*tb1(i,j,k)+8.*te1(i,j,k)*te1(i,j,k))
         !xnu*((5./3.)*ta1(i,j,k)*ta1(i,j,k)+2.*tg1(i,j,k)*tg1(i,j,k)&
         !     +2.*tb1(i,j,k)*tb1(i,j,k)+(8./3.)*th1(i,j,k)*th1(i,j,k))
      enddo
   enddo
enddo


call fine_to_coarseS(1,epsfull,tmean)
epsfullm(:,:,:)=epsfullm(:,:,:)+tmean(:,:,:)

call fine_to_coarseS(1,epsthi,tmean)
epsthim(:,:,:)=epsthim(:,:,:)+tmean(:,:,:)

call fine_to_coarseS(1,epsaxi,tmean)
epsaxim(:,:,:)=epsaxim(:,:,:)+tmean(:,:,:)

!
if (mod(itime,isave)==0) then
  
909 format('epsfullm_X',I2.2)
910 format('epsthim_X',I2.2)
911 format('epsaxim_X',I2.2)

   xp=10.
   !!X/l=5,10,15,...,45 = 9 PLANES (xp=reference)
   xvis=5.
   do i=1,3
      ivis=int((xvis+xp)/dx)+1

      write(filename, 909) int(xvis)
      call decomp_2d_write_plane(1,epsfullm,1,ivis,filename)
      
      write(filename, 910) int(xvis)
      call decomp_2d_write_plane(1,epsthim,1,ivis,filename)
      
      write(filename, 911) int(xvis)
      call decomp_2d_write_plane(1,epsaxim,1,ivis,filename)
      
      xvis=xvis+5.
   enddo
   if (nrank==0) print *,'X/l=5,10,15,... planes: write Epsilon done!'
   !!

   !!Z=0 (zlz/2) PLANE
   call decomp_2d_write_plane(1,epsfullm,3,nz/2,'epsfullm_Z0')
   call decomp_2d_write_plane(1,epsthim,3,nz/2,'epsthim_Z0')
   call decomp_2d_write_plane(1,epsaxim,3,nz/2,'epsaxim_Z0')
   if (nrank==0) print *,'z/l=0 plane: write Epsilon done!'
   !!
   !!Y=0 (yly/2) PLANE
   call decomp_2d_write_plane(1,epsfullm,2,ny/2,'epsfullm_Y0')
   call decomp_2d_write_plane(1,epsthim,2,ny/2,'epsthim_Y0')
   call decomp_2d_write_plane(1,epsaxim,2,ny/2,'epsaxim_Y0')
   if (nrank==0) print *,'y/l=0 plane: write Epsilon done!'
   !!

   
endif

end subroutine STAT_EPS_SLICE


!############################################################################
!
subroutine STATISTIC_SLICE(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
     uvmean,uwmean,vwmean,phiphimean,tmean)
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1,phi1
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: umean,vmean,wmean,uumean,vvmean,wwmean,uvmean,uwmean,vwmean,tmean
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: phimean, phiphimean
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1
real(mytype) :: xp,xvis
integer :: ivis,i
character(len=20) :: filename

!umean=ux1
call fine_to_coarseS(1,ux1,tmean)
umean(:,:,:)=umean(:,:,:)+tmean(:,:,:)

!vmean=uy1
call fine_to_coarseS(1,uy1,tmean)
vmean(:,:,:)=vmean(:,:,:)+tmean(:,:,:)

!wmean=uz1
call fine_to_coarseS(1,uz1,tmean)
wmean(:,:,:)=wmean(:,:,:)+tmean(:,:,:)

if (iscalar==1) then
   !phimean=phi1
   call fine_to_coarseS(1,phi1,tmean)
   phimean(:,:,:)=phimean(:,:,:)+tmean(:,:,:)
endif

!uumean=ux1*ux1
ta1(:,:,:)=ux1(:,:,:)*ux1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uumean(:,:,:)=uumean(:,:,:)+tmean(:,:,:)

!vvmean=uy1*uy1
ta1(:,:,:)=uy1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vvmean(:,:,:)=vvmean(:,:,:)+tmean(:,:,:)

!wwmean=uz1*uz1
ta1(:,:,:)=uz1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
wwmean(:,:,:)=wwmean(:,:,:)+tmean(:,:,:)

!uvmean=ux1*uy1
ta1(:,:,:)=ux1(:,:,:)*uy1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uvmean(:,:,:)=uvmean(:,:,:)+tmean(:,:,:)

!uwmean=ux1*uz1
ta1(:,:,:)=ux1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
uwmean(:,:,:)=uwmean(:,:,:)+tmean(:,:,:)

!vwmean=uy1*uz1
ta1(:,:,:)=uy1(:,:,:)*uz1(:,:,:)
call fine_to_coarseS(1,ta1,tmean)
vwmean(:,:,:)=vwmean(:,:,:)+tmean(:,:,:)

if (iscalar==1) then
   !phiphimean=phi1*phi1
   ta1(:,:,:)=phi1(:,:,:)*phi1(:,:,:)
   call fine_to_coarseS(1,ta1,tmean)
   phiphimean(:,:,:)=phiphimean(:,:,:)+tmean(:,:,:)
endif

if (mod(itime,isave)==0) then

900 format('umean_X',I2.2)
901 format('vmean_X',I2.2)
902 format('wmean_X',I2.2)
   
903 format('uumean_X',I2.2)
904 format('vvmean_X',I2.2)
905 format('wwmean_X',I2.2)
   
906 format('uvmean_X',I2.2)
907 format('uwmean_X',I2.2)
908 format('vwmean_X',I2.2)
   
   xp=10.
   !!X/l=5,10,15,...,45 = 9 PLANES (xp=reference)
   xvis=5.
   do i=1,3
      ivis=int((xvis+xp)/dx)+1

      write(filename, 900) int(xvis)
      call decomp_2d_write_plane(1,umean,1,ivis,filename)

      write(filename, 901) int(xvis)
      call decomp_2d_write_plane(1,vmean,1,ivis,filename)

      write(filename, 902) int(xvis)
      call decomp_2d_write_plane(1,wmean,1,ivis,filename)

      write(filename, 903) int(xvis)
      call decomp_2d_write_plane(1,uumean,1,ivis,filename)

      write(filename, 904) int(xvis)
      call decomp_2d_write_plane(1,vvmean,1,ivis,filename)

      write(filename, 905) int(xvis)
      call decomp_2d_write_plane(1,wwmean,1,ivis,filename)

      write(filename, 906) int(xvis)
      call decomp_2d_write_plane(1,uvmean,1,ivis,filename)

      write(filename, 907) int(xvis)
      call decomp_2d_write_plane(1,uwmean,1,ivis,filename)

      write(filename, 908) int(xvis)
      call decomp_2d_write_plane(1,vwmean,1,ivis,filename)

      xvis=xvis+5.
   enddo
   if (nrank==0) print *,'x/l=5,10,15,... planes: write stat arrays done!'
   !!

   !!Z=0 (zlz/2) PLANE
   call decomp_2d_write_plane(1,umean,3,nz/2,'umean_Z0')
   call decomp_2d_write_plane(1,vmean,3,nz/2,'vmean_Z0')
   call decomp_2d_write_plane(1,wmean,3,nz/2,'wmean_Z0')
   call decomp_2d_write_plane(1,uumean,3,nz/2,'uumean_Z0')
   call decomp_2d_write_plane(1,vvmean,3,nz/2,'vvmean_Z0')
   call decomp_2d_write_plane(1,wwmean,3,nz/2,'wwmean_Z0')
   call decomp_2d_write_plane(1,uvmean,3,nz/2,'uvmean_Z0')
   call decomp_2d_write_plane(1,uwmean,3,nz/2,'uwmean_Z0')
   call decomp_2d_write_plane(1,vwmean,3,nz/2,'vwmean_Z0')
   if (nrank==0) print *,'z/l=0 plane: write stat arrays done!'
   !!

   !!Y=0 (ylz/2) PLANE
   call decomp_2d_write_plane(1,umean,2,ny/2,'umean_Y0')
   call decomp_2d_write_plane(1,vmean,2,ny/2,'vmean_Y0')
   call decomp_2d_write_plane(1,wmean,2,ny/2,'wmean_Y0')
   call decomp_2d_write_plane(1,uumean,2,ny/2,'uumean_Y0')
   call decomp_2d_write_plane(1,vvmean,2,ny/2,'vvmean_Y0')
   call decomp_2d_write_plane(1,wwmean,2,ny/2,'wwmean_Y0')
   call decomp_2d_write_plane(1,uvmean,2,ny/2,'uvmean_Y0')
   call decomp_2d_write_plane(1,uwmean,2,ny/2,'uwmean_Y0')
   call decomp_2d_write_plane(1,vwmean,2,ny/2,'vwmean_Y0')
   if (nrank==0) print *,'y/l=0 plane: write stat arrays done!'
   !!
  
endif

end subroutine STATISTIC_SLICE


!############################################################################
!
subroutine STAT_EPS_SLICE (ux1,uy1,uz1,epsfull,epsthi,epsaxi,epsfullm,epsthim,epsaxim,tmean,&
     ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG)
!
!
!############################################################################

USE param
USE variables
USE decomp_2d
USE decomp_2d_io

implicit none

TYPE(DECOMP_INFO) :: phG
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ux1,uy1,uz1
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: epsfull,epsthi,epsaxi
real(mytype),dimension(xsize(1),xsize(2),xsize(3)) :: ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1
real(mytype),dimension(ysize(1),ysize(2),ysize(3)) :: ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2
real(mytype),dimension(zsize(1),zsize(2),zsize(3)) :: ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3
real(mytype),dimension(xszS(1),xszS(2),xszS(3)) :: epsfullm,epsthim,epsaxim,tmean

integer :: code,icomplet,ivis
integer :: ijk,nvect1,nvect2,nvect3,i,j,k
character(len=20) nfichier,nfichier1
character(len=20) :: filename
real(mytype) :: xvis,xp

!x-derivatives
call derx (ta1,ux1,di1,sx,ffx,fsx,fwx,xsize(1),xsize(2),xsize(3),0)
call derx (tb1,uy1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
call derx (tc1,uz1,di1,sx,ffxp,fsxp,fwxp,xsize(1),xsize(2),xsize(3),1)
!y-derivatives
call transpose_x_to_y(ux1,td2)
call transpose_x_to_y(uy1,te2)
call transpose_x_to_y(uz1,tf2)
call dery (ta2,td2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
call dery (tb2,te2,di2,sy,ffy,fsy,fwy,ppy,ysize(1),ysize(2),ysize(3),0)
call dery (tc2,tf2,di2,sy,ffyp,fsyp,fwyp,ppy,ysize(1),ysize(2),ysize(3),1)
!!z-derivatives
call transpose_y_to_z(td2,td3)
call transpose_y_to_z(te2,te3)
call transpose_y_to_z(tf2,tf3)
call derz (ta3,td3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tb3,te3,di3,sz,ffzp,fszp,fwzp,zsize(1),zsize(2),zsize(3),1)
call derz (tc3,tf3,di3,sz,ffz,fsz,fwz,zsize(1),zsize(2),zsize(3),0)
!!all back to x-pencils
call transpose_z_to_y(ta3,td2)
call transpose_z_to_y(tb3,te2)
call transpose_z_to_y(tc3,tf2)
call transpose_y_to_x(td2,tg1)
call transpose_y_to_x(te2,th1)
call transpose_y_to_x(tf2,ti1)
call transpose_y_to_x(ta2,td1)
call transpose_y_to_x(tb2,te1)
call transpose_y_to_x(tc2,tf1)
!du/dx=ta1 du/dy=td1 and du/dz=tg1
!dv/dx=tb1 dv/dy=te1 and dv/dz=th1
!dw/dx=tc1 dw/dy=tf1 and dw/dz=ti1


!ENERGY DISSIPATION (FULL, THI and AXI HYPOTHESIS)
do k=1,xsize(3)
   do j=1,xsize(2)
      do i=1,xsize(1)
         epsfull(i,j,k)= xnu*( 2.*(ta1(i,j,k)*ta1(i,j,k)+te1(i,j,k)*te1(i,j,k)&
              +ti1(i,j,k)*ti1(i,j,k)) &
              +td1(i,j,k)*td1(i,j,k)+tb1(i,j,k)*tb1(i,j,k)+tg1(i,j,k)*tg1(i,j,k)&
              +tc1(i,j,k)*tc1(i,j,k)+th1(i,j,k)*th1(i,j,k)+tf1(i,j,k)*tf1(i,j,k)&
              +2.*(td1(i,j,k)*tb1(i,j,k)+tg1(i,j,k)*tc1(i,j,k)+th1(i,j,k)*tf1(i,j,k)))

         epsthi(i,j,k)=15.*xnu*ta1(i,j,k)*ta1(i,j,k)
         
         epsaxi(i,j,k)= xnu*(-ta1(i,j,k)*ta1(i,j,k)+2.*td1(i,j,k)*td1(i,j,k)&
              +2.*tb1(i,j,k)*tb1(i,j,k)+8.*te1(i,j,k)*te1(i,j,k))
         !xnu*((5./3.)*ta1(i,j,k)*ta1(i,j,k)+2.*tg1(i,j,k)*tg1(i,j,k)&
         !     +2.*tb1(i,j,k)*tb1(i,j,k)+(8./3.)*th1(i,j,k)*th1(i,j,k))
      enddo
   enddo
enddo


call fine_to_coarseS(1,epsfull,tmean)
epsfullm(:,:,:)=epsfullm(:,:,:)+tmean(:,:,:)

call fine_to_coarseS(1,epsthi,tmean)
epsthim(:,:,:)=epsthim(:,:,:)+tmean(:,:,:)

call fine_to_coarseS(1,epsaxi,tmean)
epsaxim(:,:,:)=epsaxim(:,:,:)+tmean(:,:,:)

!
if (mod(itime,isave)==0) then
  
909 format('epsfullm_X',I2.2)
910 format('epsthim_X',I2.2)
911 format('epsaxim_X',I2.2)

   xp=10.
   !!X/l=5,10,15,...,45 = 9 PLANES (xp=reference)
   xvis=5.
   do i=1,3
      ivis=int((xvis+xp)/dx)+1

      write(filename, 909) int(xvis)
      call decomp_2d_write_plane(1,epsfullm,1,ivis,filename)
      
      write(filename, 910) int(xvis)
      call decomp_2d_write_plane(1,epsthim,1,ivis,filename)
      
      write(filename, 911) int(xvis)
      call decomp_2d_write_plane(1,epsaxim,1,ivis,filename)
      
      xvis=xvis+5.
   enddo
   if (nrank==0) print *,'X/l=5,10,15,... planes: write Epsilon done!'
   !!

   !!Z=0 (zlz/2) PLANE
   call decomp_2d_write_plane(1,epsfullm,3,nz/2,'epsfullm_Z0')
   call decomp_2d_write_plane(1,epsthim,3,nz/2,'epsthim_Z0')
   call decomp_2d_write_plane(1,epsaxim,3,nz/2,'epsaxim_Z0')
   if (nrank==0) print *,'z/l=0 plane: write Epsilon done!'
   !!
   !!Y=0 (yly/2) PLANE
   call decomp_2d_write_plane(1,epsfullm,2,ny/2,'epsfullm_Y0')
   call decomp_2d_write_plane(1,epsthim,2,ny/2,'epsthim_Y0')
   call decomp_2d_write_plane(1,epsaxim,2,ny/2,'epsaxim_Y0')
   if (nrank==0) print *,'y/l=0 plane: write Epsilon done!'
   !!

   
endif

end subroutine STAT_EPS_SLICE
