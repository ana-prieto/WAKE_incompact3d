
PROGRAM incompact3d

  USE decomp_2d
  USE decomp_2d_poisson
  use decomp_2d_io
  USE variables
  USE param
  USE var
  USE MPI
  USE IBM
  USE derivX
  USE derivZ

  implicit none

  integer :: code,nlock,i,j,k,ii,iii,bcx,bcy,bcz,fh,ierror, ncollect
  real(mytype) :: x,y,z,tmp1
  double precision :: t1,t2
  character(len=20) :: filename

  TYPE(DECOMP_INFO) :: phG,ph1,ph2,ph3,ph4

  CALL MPI_INIT(code)
  call decomp_2d_init(nx,ny,nz,p_row,p_col)
  call init_coarser_mesh_statS(nstat,nstat,nstat,.true.) !start from 1 == true
  call init_coarser_mesh_statV(nvisu,nvisu,nvisu,.true.) !start from 1 == true
  call init_coarser_mesh_statP(nprobe,nprobe,nprobe,.true.) !start from 1 == true
  call parameter()

  call init_variables

  call schemes()

  if (ifilter.eq.1) call filter()

  if (nclx==0) then
     bcx=0
  else
     bcx=1
  endif
  if (ncly==0) then
     bcy=0
  else
     bcy=1
  endif
  if (nclz==0) then
     bcz=0
  else
     bcz=1
  endif

  call decomp_2d_poisson_init(bcx,bcy,bcz)

  call decomp_info_init(nxm,nym,nzm,phG)


!!$  !if you want to collect N snapshots randomly on IT time steps
!!$  call collect_data() !it will generate N random time steps

  if (ilit==0) call init(ux1,uy1,uz1,ep1,phi1,gx1,gy1,gz1,phis1,hx1,hy1,hz1,phiss1)  
  if (ilit==1) call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
       px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,0)
  call test_speed_min_max(ux1,uy1,uz1)
  if (iscalar==1) call test_scalar_min_max(phi1)

  !phi1=0.;phis1=0.;phiss1=0.

  !array for stat to zero
  !umean=0.;vmean=0.;wmean=0.
  !uumean=0.;vvmean=0.;wwmean=0.
  !uvmean=0.;uwmean=0.;vwmean=0.
  !phimean=0.;phiphimean=0.
  !pmean=0.

  !dudxm=0.;dudym=0.;dudzm=0.
  !dvdxm=0.;dvdym=0.;dvdzm=0.
  !dwdxm=0.;dwdym=0.;dwdzm=0.

  !dudx2m=0.;dudy2m=0.;dudz2m=0.
  !dvdx2m=0.;dvdy2m=0.;dvdz2m=0.
  !dwdx2m=0.;dwdy2m=0.;dwdz2m=0.

  !dudydvdxm=0.;dudzdwdxm=0.;dvdzdwdym=0.

  t1 = MPI_WTIME()

  !div: nx ny nz --> nxm ny nz --> nxm nym nz --> nxm nym nzm
  call decomp_info_init(nxm, nym, nzm, ph1)
  call decomp_info_init(nxm, ny, nz, ph4)

  !gradp: nxm nym nzm -> nxm nym nz --> nxm ny nz --> nx ny nz
  call decomp_info_init(nxm, ny, nz, ph2)  
  call decomp_info_init(nxm, nym, nz, ph3) 

  !for the wake case
  if(itype==1) then
     call plate(ep1)
  endif
  
  ncollect = 2 ! saving frequency

  do itime=ifirst,ilast
     t=(itime-1)*dt
     if (nrank==0) then
        write(*,1001) itime,t
1001    format('Time step =',i7,', Time unit =',F9.3)
     endif

     do itr=1,iadvance_time

        if (nclx.eq.2) then
           call inflow (ux1,uy1,uz1,phi1) !X PENCILS
           call outflow(ux1,uy1,uz1,phi1) !X PENCILS 
        endif

        !X-->Y-->Z-->Y-->X
        call convdiff(ux1,uy1,uz1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
             ux2,uy2,uz2,ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
             ux3,uy3,uz3,ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3)

        if (iscalar==1) then
           call scalar(ux1,uy1,uz1,phi1,phis1,phiss1,di1,tg1,th1,ti1,td1,&
                uy2,uz2,phi2,di2,ta2,tb2,tc2,td2,uz3,phi3,di3,ta3,tb3,ep1) 
        endif

        !BUFFER ZONE
        !X PENCILS
        !call fringe(ux1,uy1,uz1,ta1,tb1,tc1)
        !   

        !X PENCILS
        call intt (ux1,uy1,uz1,gx1,gy1,gz1,hx1,hy1,hz1,ta1,tb1,tc1) 

        call pre_correc(ux1,uy1,uz1)

        if (ivirt==1) then !solid body old school
           !we are in X-pencil
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,1)
           call force_plate(ux1,uy1,uz1,ep1)
           call corgp_IBM(ux1,uy1,uz1,px1,py1,pz1,2)
        endif

        !X-->Y-->Z
        call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
             td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,pp3,&
             nxmsize,nymsize,nzmsize,ph1,ph3,ph4,1)       

        !POISSON Z-->Z 
        call decomp_2d_poisson_stg(pp3,bcx,bcy,bcz)

        !Z-->Y-->X
        call gradp(px1,py1,pz1,di1,td2,tf2,ta2,tb2,tc2,di2,&
             ta3,tc3,di3,pp3,nxmsize,nymsize,nzmsize,ph2,ph3)

        !X PENCILS
        call corgp(ux1,ux2,uy1,uz1,px1,py1,pz1) 

        !does not matter -->output=DIV U=0 (in dv3)
        call divergence (ux1,uy1,uz1,ep1,ta1,tb1,tc1,di1,td1,te1,tf1,&
             td2,te2,tf2,di2,ta2,tb2,tc2,ta3,tb3,tc3,di3,td3,te3,tf3,dv3,&
             nxmsize,nymsize,nzmsize,ph1,ph3,ph4,2)

        call test_speed_min_max(ux1,uy1,uz1)
        if (iscalar==1) call test_scalar_min_max(phi1)

     enddo

     if (mod(itime,ncollect)-1 .EQ. 0) then
       call VISU_STG  (pp3,ppx,ta1,tb1,di1,ta2,tb2,di2,&
            ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu) !pression in tb1 & ppx on cell surfaces

       call VISU_SUBDOMAIN  (ux1,uy1,uz1,ppx,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
                ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
                ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)
     endif


     !!3D STATS 
!     call STATISTIC(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
!          uvmean,uwmean,vwmean,phiphimean,tmean)
!
!     call STAT_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
!          ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,pmean,tmean)
!
!     call STAT_EPS (ux1,uy1,uz1,&
!          dudxm,dudym,dudzm,&
!          dvdxm,dvdym,dvdzm,&
!          dwdxm,dwdym,dwdzm,&
!          dudx2m,dudy2m,dudz2m,&
!          dvdx2m,dvdy2m,dvdz2m,&
!          dwdx2m,dwdy2m,dwdz2m,&
!          dudydvdxm,dudzdwdxm,dvdzdwdym,tempgrad,&
!          tmean,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
!          ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
!          ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG)
     !!




!!$   !!SLICES STATS
!!$   call  STATISTIC_SLICE(ux1,uy1,uz1,phi1,ta1,umean,vmean,wmean,phimean,uumean,vvmean,wwmean,&
!!$        uvmean,uwmean,vwmean,phiphimean,tmean)
!!$
!!$   call STAT_EPS_SLICE  (ux1,uy1,uz1,epsfull,epsthi,epsaxi,epsfullm,epsthim,epsaxim,tmean,&
!!$     ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
!!$     ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
!!$     ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG)
!!$   !!


1042 format(I6)
     if (mod(itime,isave)==0) then 
        call restart(ux1,uy1,uz1,ep1,pp3,phi1,gx1,gy1,gz1,&
             px1,py1,pz1,phis1,hx1,hy1,hz1,phiss1,phG,1)

        !Sauvegarde de l'iter de restart dans "iter.dat"
        if(nrank==0) then
           open(41,file='iter.dat',status="unknown",form='formatted')
           write(41,1042) itime+1
           close(41)
        endif

     endif

   if (mod(itime,imodulo)==0) then
      
    !  call VISU_INSTA(ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
    !       ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
    !       ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)

    !  call VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
    !       ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu)

      call VISU_SLICE  (ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
           ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
           ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)

   endif



     !PROBES 
!     call VISU_INSTA(ux1,uy1,uz1,phi1,ta1,tb1,tc1,td1,te1,tf1,tg1,th1,ti1,di1,&
!          ta2,tb2,tc2,td2,te2,tf2,tg2,th2,ti2,tj2,di2,&
!          ta3,tb3,tc3,td3,te3,tf3,tg3,th3,ti3,di3,phG,uvisu)


!!$     !collection of u,v,w,phi (and p eventually) using idata
!!$     !1-setup istart to 1, collection for 50000.
!!$     if (itime==588001) then
!!$        ii=962
!!$        iii=959
!!$     endif
!!$     if (itime==idata(ii)) then
!!$925     format('ux_global',I4.4)
!!$        write(filename, 925) iii
!!$        call decomp_2d_write_one(nx_global,ny_global,nz_global,1,ux1,filename)
!!$926     format('uy_global',I4.4)
!!$        write(filename, 926) iii
!!$        call decomp_2d_write_one(nx_global,ny_global,nz_global,1,uy1,filename)
!!$927     format('uz_global',I4.4)
!!$        write(filename, 927) iii
!!$        call decomp_2d_write_one(nx_global,ny_global,nz_global,1,uz1,filename)
!!$        if (iscalar==1) then
!!$997        format('phi_global',I4.4)
!!$           write(filename, 997) iii
!!$           call decomp_2d_write_one(nx_global,ny_global,nz_global,1,phi1,filename)
!!$        endif
!!$        call VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
!!$             ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu) !pression in tb1
!!$928     format('pp_global',I4.4)
!!$        write(filename, 928) iii
!!$        call decomp_2d_write_one(nx_global,ny_global,nz_global,1,tb1,filename)
!!$        if (nrank==0) print *,'COLLECTION RANDOM DATA FOR ITIME=',itime
!!$        iii=iii+1
!!$        ii=ii+1
!!$        if (itime==idata(ii)) ii=ii+1 !in case of the generation of 2 same numbers
!!$     endif


!!$     !Regularly spaced collection of u,v,w,phi (and p eventually) using imodulo
!!$     if(itime-ifirst+1==1) then
!!$        iii=1
!!$     endif
!!$     if (mod(itime,imodulo)==0) then
!!$925     format('ux_global',I4.4)
!!$        write(filename, 925) iii
!!$        call decomp_2d_write_one(nx_global,ny_global,nz_global,1,ux1,filename)
!!$926     format('uy_global',I4.4)
!!$        write(filename, 926) iii
!!$        call decomp_2d_write_one(nx_global,ny_global,nz_global,1,uy1,filename)
!!$927     format('uz_global',I4.4)
!!$        write(filename, 927) iii
!!$        call decomp_2d_write_one(nx_global,ny_global,nz_global,1,uz1,filename)
!!$        if (iscalar==1) then
!!$997        format('phi_global',I4.4)
!!$           write(filename, 997) iii
!!$           call decomp_2d_write_one(nx_global,ny_global,nz_global,1,phi1,filename)
!!$        endif
!!$        call VISU_PRE (pp3,ta1,tb1,di1,ta2,tb2,di2,&
!!$             ta3,di3,nxmsize,nymsize,nzmsize,phG,ph2,ph3,uvisu) !pression in tb1
!!$928     format('pp_global',I4.4)
!!$        write(filename, 928) iii
!!$        call decomp_2d_write_one(nx_global,ny_global,nz_global,1,tb1,filename)
!!$        if (nrank==0) print *,'COLLECTION RANDOM DATA FOR ITIME=',itime
!!$        iii=iii+1
!!$     endif

  enddo

  t2=MPI_WTIME()-t1
  call MPI_ALLREDUCE(t2,t1,1,MPI_REAL8,MPI_SUM, &
       MPI_COMM_WORLD,code)
  if (nrank==0) print *,'time per time_step: ', &
       t1/float(nproc)/(ilast-ifirst+1),' seconds'
  if (nrank==0) print *,'simulation with nx*ny*nz=',nx,ny,nz,'mesh nodes'
  if (nrank==0) print *,'Mapping p_row*p_col=',p_row,p_col

  !call decomp_2d_poisson_finalize
  call decomp_2d_finalize
  CALL MPI_FINALIZE(code)

end PROGRAM incompact3d
