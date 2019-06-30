Program Main
  implicit none


  !number of grid points with ghost cells
  integer ( kind = 4 ), parameter :: nx = 302
  integer ( kind = 4 ), parameter :: ny = 66
  integer ( kind = 4 ), parameter :: nz = 1282

  integer ( kind = 4 ) :: nxi1
  integer ( kind = 4 ) :: nyi1
  integer ( kind = 4 ) :: nzi1

  integer ( kind = 4 ) :: nxi2
  integer ( kind = 4 ) :: nyi2
  integer ( kind = 4 ) :: nzi2

  real    ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nz)
  real    ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nz)
  real    ( kind = 8 ) :: ru(nx),rp(nx)

  integer ( kind = 4 ) :: i,j,k,jp,nstep,imaxval,jmaxval,kmaxval,tmp,kstart,kend,s1,jend,expu
  integer ( kind = 4 ) :: tag
  real    ( kind = 8 ) :: time,dtm1,grav
  character(len=128)   :: u_filename,u_filename_out,v_filename,v_filename_out,w_filename,w_filename_out
  character(len=128)   :: p_filename,p_filename_out
  character(len=128)   :: d_filename,d_filename_out

  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: u,v,w,p,d,debug
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: ui1,vi1,wi1,pi1,di1
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: ui2,vi2,wi2,pi2,di2
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: ui3,vi3,wi3,pi3,di3
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: ui4,vi4,wi4,pi4,di4
  real    ( kind = 8 ), allocatable, dimension(:,:,:) :: ui5,vi5,wi5,pi5,di5
  real    ( kind = 8 ) :: r(nx),theta(ny),z(nz),dtheta

  integer ( kind = 4 ) :: iu,iv,iw,ip,id


  ! interpolate on == 1 / off==0

  iu = 1
  iv = 1
  iw = 1
  ip = 1
  id = 1

  ! points in the first interpolation

  !nxi1 = 2*(nx-2)+2

  nxi1 = nx

  nyi1 = 2*(ny-2)+2

  !nyi1 = ny
  nzi1 = 2*(nz-2)+2

  !nzi1=nz

  ! points in the second interpolation

  !nxi2 = 2*(nxi1-2)+2
  nxi2 = nxi1
  nzi2 = 2*(nzi1-2)+2
  !nzi2 =nzi1

  u_filename = 'u_00030001.res'
  u_filename_out = 'u_00030001_in_yzz.res'

  v_filename = 'v_00030001.res'
  v_filename_out = 'v_00030001_in_yzz.res'

  w_filename = 'w_00030001.res'
  w_filename_out = 'w_00030001_in_yzz.res'

  p_filename = 'p_00030001.res'
  p_filename_out = 'p_00030001_in_yzz.res'

  d_filename = 'densp_00030001.res'
  d_filename_out = 'densp_00030001_in_yzz.res'


   write(*,*) 'Original restartfile size: nx=',nx,' ny=',ny,' nz=',nz
   write(*,*) 'Interpolated restartfile size: nx=',nxi2,' ny=',nyi2,' nzi2=',nzi2


  ! ------------- U velocity  ----------------

  if (iu == 1) then


   write(*,*) 'Interpolating and imposing the BC of: U'

  allocate(u(nx,ny,nz))
  allocate(ui1(nxi1,ny,nz))
  allocate(ui2(nxi2,ny,nz))
  allocate(ui3(nxi2,nyi1,nz))
  allocate(ui4(nxi2,nyi1,nzi1))
  allocate(ui5(nxi2,nyi1,nzi2))

  u=0.0d0;ui1=0.0d0;ui2=0.0d0;ui3=0.0d0;ui4=0.0d0;ui5=0.0d0

  call read_restart(u_filename,nx,ny,nz,u,time,dtm1,grav,jp,nstep)


   ! -- 1st X interpolation --

   !do  i=1,nx-1
   !    ui1(2*i-1,:,:)=u(i,:,:)
   !enddo

   !do  i=2,nxi1-1,2
   !    ui1(i,:,:)=0.5*(ui1(i+1,:,:)+ui1(i-1,:,:))
   !enddo

   !ghost
   !ui1(nxi1,:,:)=0.5*(u(nx,:,:)+u(nx-1,:,:))


   ui1 = u

   deallocate(u)

   ! -- 2nd X interpolation --


     ui2=ui1

   ! do  i=1,nxi1-1
   !     ui2(2*i-1,:,:)=ui1(i,:,:)
   ! enddo
   !
   ! do  i=2,nxi2-1,2
   !     ui2(i,:,:)=0.5*(ui2(i+1,:,:)+ui2(i-1,:,:))
   ! enddo
   !
   ! !ghost
   ! ui2(nxi2,:,:)=0.5*(ui1(nxi1,:,:)+ui1(nxi1-1,:,:))
   !

   deallocate(ui1)

   ! -- Y interpolation --

   ui2(:,1,nz)=ui2(:,ny-1,nz)
   ui2(:,ny,nz)=ui2(:,2,nz)

   do  j=1,ny-1
      ui3(:,2*j-1,:)=0.75*ui2(:,j,:)+0.25*ui2(:,j+1,:)
   enddo

   do  j=2,ny
      ui3(:,2*j-2,:)=0.25*ui2(:,j-1,:)+0.75*ui2(:,j,:)
   enddo

   !ui3=ui2

    ui3(:,1,nz)=0.0
    ui3(:,nyi1,nz)=0.0

   deallocate(ui2)

   ! -- 1st Z interpolation --

   do  k=1,nz-1
      ui4(:,:,2*k-1)=0.75*ui3(:,:,k)+0.25*ui3(:,:,k+1)
   enddo

   do  k=2,nz
      ui4(:,:,2*k-2)=0.25*ui3(:,:,k-1)+0.75*ui3(:,:,k)
   enddo

!   ui4=ui3

   deallocate(ui3)

   ! -- 2nd Z interpolation --

   do  k=1,nzi1-1
      ui5(:,:,2*k-1)=0.75*ui4(:,:,k)+0.25*ui4(:,:,k+1)
   enddo

   do  k=2,nzi1
      ui5(:,:,2*k-2)=0.25*ui4(:,:,k-1)+0.75*ui4(:,:,k)
   enddo

!    ui5=ui4

   deallocate(ui4)


   ! boundary conditions

   !ui5(nxi2,:,:)=ui5(nxi2-1,:,:)
   !ui5(:,nyi1,:)=ui5(:,2,:)
   !ui5(:,1,:)=ui5(:,nyi1-1,:)
   ui5(nxi2,:,nzi2)=0.0

   call write_restart(u_filename_out,nxi2,nyi1,nzi2,ui5,time,dtm1,grav,jp,nstep)

   deallocate(ui5)

   endif




   ! ------------- V velocity  ----------------

   if (iv==1) then


   write(*,*) 'Interpolating and imposing the BC of: V'

   allocate(v(nx,ny,nz))
   allocate(vi1(nx,nyi1,nz))
   allocate(vi2(nxi1,nyi1,nz))
   allocate(vi3(nxi2,nyi1,nz))
   allocate(vi4(nxi2,nyi1,nzi1))
   allocate(vi5(nxi2,nyi1,nzi2))


   v=0.0d0;vi1=0.0d0;vi2=0.0d0;vi3=0.0d0;vi4=0.0d0;vi5=0.0d0

   call read_restart(v_filename,nx,ny,nz,v,time,dtm1,grav,jp,nstep)

   !  -- Y interpolation --

   !v(:,1,nz)=v(:,ny-1,nz)
   !v(:,ny,nz)=v(:,2,nz)

   do  j=1,ny-1
       vi1(:,2*j-1,:)=v(:,j,:)
   enddo

   do  j=2,nyi1-1,2
       vi1(:,j,:)=0.5*(vi1(:,j+1,:)+vi1(:,j-1,:))
   enddo

   !ghost
   vi1(:,nyi1,:)=0.5*(v(:,ny,:)+v(:,ny-1,:))

   !not updated cells
   vi1(:,nyi1,nz)=0.0


   !vi1=v

   deallocate(v)

   ! -- 1st X interpolation --


!   do  i=1,nx-1
!      vi2(2*i-1,:,:)=0.75*vi1(i,:,:)+0.25*vi1(i+1,:,:)
!   enddo
!
!   do  i=2,nx
!      vi2(2*i-2,:,:)=0.25*vi1(i-1,:,:)+0.75*vi1(i,:,:)
!   enddo

   vi2 =vi1

   deallocate(vi1)

   ! -- 2nd X interpolation --


    vi3=vi2
   ! do  i=1,nxi1-1
   !    vi3(2*i-1,:,:)=0.75*vi2(i,:,:)+0.25*vi2(i+1,:,:)
   ! enddo
   !
   ! do  i=2,nxi1
   !    vi3(2*i-2,:,:)=0.25*vi2(i-1,:,:)+0.75*vi2(i,:,:)
   ! enddo

   deallocate(vi2)


   ! -- 1st Z interpolation --

   do  k=1,nz-1
      vi4(:,:,2*k-1)=0.75*vi3(:,:,k)+0.25*vi3(:,:,k+1)
   enddo

   do  k=2,nz
      vi4(:,:,2*k-2)=0.25*vi3(:,:,k-1)+0.75*vi3(:,:,k)
   enddo


    !vi4=vi3

   deallocate(vi3)

   ! -- 2nd Z interpolation --


   do  k=1,nzi1-1
      vi5(:,:,2*k-1)=0.75*vi4(:,:,k)+0.25*vi4(:,:,k+1)
   enddo

   do  k=2,nzi1
      vi5(:,:,2*k-2)=0.25*vi4(:,:,k-1)+0.75*vi4(:,:,k)
   enddo


   ! vi5=vi4

   deallocate(vi4)


   ! boundary conditions

   !vi5(nxi2,:,:)=vi5(nxi2-1,:,:)
   !vi5(:,nyi1,:)=vi5(:,2,:)
   !vi5(:,1,:)=vi5(:,nyi1-1,:)
   vi5(nxi2,:,nzi2)=0.0

   call write_restart(v_filename_out,nxi2,nyi1,nzi2,vi5,time,dtm1,grav,jp,nstep)


   deallocate(vi5)

   endif



  ! ------------- W velocity  ----------------

  if (iw==1) then


   write(*,*) 'Interpolating and imposing the BC of: W'
   write(*,*) 'The bug in the update of the azymuthal nz-1 gc is fixed'

   allocate(w(nx,ny,nz))
   allocate(wi1(nx,ny,nzi1))
   allocate(wi2(nx,ny,nzi2))
   allocate(wi3(nxi1,ny,nzi2))
   allocate(wi4(nxi2,ny,nzi2))
   allocate(wi5(nxi2,nyi1,nzi2))

   w=0.0d0;wi1=0.0d0;wi2=0.0d0;wi3=0.0d0;wi4=0.0d0;wi5=0.0d0

   call read_restart(w_filename,nx,ny,nz,w,time,dtm1,grav,jp,nstep)

   write(*,*) 'Solving problem with w(nx,:,nz-1) not being updated'
   write(*,*) 'Solving problem with w(:,0,nz-1) not being updated'

   w(nx,:,nz-1)=w(nx,:,nz-2)
   w(:,1,nz-1)=w(:,1,nz-2)

   ! -- 1st Z interpolation --

   do  k=1,nz-1
       wi1(:,:,2*k-1)=w(:,:,k)
   enddo

   do  k=2,nzi1-1,2
       wi1(:,:,k)=0.5*(wi1(:,:,k+1)+wi1(:,:,k-1))
   enddo

   !ghost
   wi1(:,:,nzi1)=0.5*(w(:,:,nz)+w(:,:,nz-1))

    !wi1=w


   deallocate(w)

   ! -- 2nd Z interpolation --

   do  k=1,nzi1-1
       wi2(:,:,2*k-1)=wi1(:,:,k)
   enddo

   do  k=2,nzi2-1,2
       wi2(:,:,k)=0.5*(wi2(:,:,k+1)+wi2(:,:,k-1))
   enddo

   !ghost
   !wi2(:,:,nzi2)=0.5*(wi1(:,:,nzi1)+wi1(:,:,nzi1-1))
   wi2(:,:,nzi2)=1.0d0

   !wi2=wi1


   deallocate(wi1)


   ! -- 1st X interpolation --

!   do  i=1,nx-1
!      wi3(2*i-1,:,:)=0.75*wi2(i,:,:)+0.25*wi2(i+1,:,:)
!   enddo
!
!   do  i=2,nx
!      wi3(2*i-2,:,:)=0.25*wi2(i-1,:,:)+0.75*wi2(i,:,:)
!   enddo

    wi3=wi2

   deallocate(wi2)

   ! -- 2nd X interpolation --

     wi4=wi3
   ! do  i=1,nxi1-1
   !    wi4(2*i-1,:,:)=0.75*wi3(i,:,:)+0.25*wi3(i+1,:,:)
   ! enddo
   !
   ! do  i=2,nxi1
   !    wi4(2*i-2,:,:)=0.25*wi3(i-1,:,:)+0.75*wi3(i,:,:)
   ! enddo
   !


   deallocate(wi3)

   ! -- Y interpolation --

   wi4(:,1,:)=wi4(:,ny-1,:)
   wi4(:,ny,:)=wi4(:,2,:)

   do  j=1,ny-1
      wi5(:,2*j-1,:)=0.75*wi4(:,j,:)+0.25*wi4(:,j+1,:)
   enddo

   do  j=2,ny
      wi5(:,2*j-2,:)=0.25*wi4(:,j-1,:)+0.75*wi4(:,j,:)
   enddo

   wi5(:,1,nzi2)=1.0
   wi5(:,nyi2,nzi2)=1.0

   deallocate(wi4)


   !boundary conditions
   wi5(:,:,1)=1.0d0
   wi5(:,:,nzi2)=1.0d0
   !wi5(nxi2,:,:)=wi5(nxi2-1,:,:)
   !wi5(:,nyi1,:)=wi5(:,2,:)
   !wi5(:,1,:)=wi5(:,nyi1-1,:)

   call write_restart(w_filename_out,nxi2,nyi1,nzi2,wi5,time,dtm1,grav,jp,nstep)


   deallocate(wi5)

   endif


  ! ------------- P  ----------------

   if (ip==1) then


   write(*,*) 'Interpolating and imposing the BC of: Pressure'

   allocate(p(nx,ny,nz))
   allocate(pi1(nxi1,ny,nz))
   allocate(pi2(nxi2,ny,nz))
   allocate(pi3(nxi2,nyi1,nz))
   allocate(pi4(nxi2,nyi1,nzi1))
   allocate(pi5(nxi2,nyi1,nzi2))


   p=0.0d0;pi1=0.0d0;pi2=0.0d0;pi3=0.0d0;pi4=0.0d0;pi5=0.0d0

   call read_restart(p_filename,nx,ny,nz,p,time,dtm1,grav,jp,nstep)


   ! -- 1st X interpolation --

!   do  i=1,nx-1
!      pi1(2*i-1,:,:)=0.75*p(i,:,:)+0.25*p(i+1,:,:)
!   enddo
!
!   do  i=2,nx
!      pi1(2*i-2,:,:)=0.25*p(i-1,:,:)+0.75*p(i,:,:)
!   enddo

   pi1 = p

   deallocate(p)

   ! -- 2nd X interpolation --

     pi2=pi1
   ! do  i=1,nxi1-1
   !    pi2(2*i-1,:,:)=0.75*pi1(i,:,:)+0.25*pi1(i+1,:,:)
   ! enddo
   !
   ! do  i=2,nxi1
   !    pi2(2*i-2,:,:)=0.25*pi1(i-1,:,:)+0.75*pi1(i,:,:)
   ! enddo

   deallocate(pi1)

   ! -- Y interpolation --

   !pi2(:,1,nz)=pi2(:,ny-1,nz)
   !pi2(:,ny,nz)=pi2(:,2,nz)

   do  j=1,ny-1
      pi3(:,2*j-1,:)=0.75*pi2(:,j,:)+0.25*pi2(:,j+1,:)
   enddo

   do  j=2,ny
      pi3(:,2*j-2,:)=0.25*pi2(:,j-1,:)+0.75*pi2(:,j,:)
   enddo

   !pi3=pi2

   deallocate(pi2)

   ! -- 1st Z interpolation --

   do  k=1,nz-1
      pi4(:,:,2*k-1)=0.75*pi3(:,:,k)+0.25*pi3(:,:,k+1)
   enddo

   do  k=2,nz
      pi4(:,:,2*k-2)=0.25*pi3(:,:,k-1)+0.75*pi3(:,:,k)
   enddo

   !pi4=pi3

   deallocate(pi3)

   ! -- 2nd Z interpolation --


   do  k=1,nzi1-1
      pi5(:,:,2*k-1)=0.75*pi4(:,:,k)+0.25*pi4(:,:,k+1)
   enddo

   do  k=2,nzi1
      pi5(:,:,2*k-2)=0.25*pi4(:,:,k-1)+0.75*pi4(:,:,k)
   enddo

    !pi5=pi4

   deallocate(pi4)

   !boundary conditions
   pi5(:,:,1)=pi5(:,:,2)
   !pi5(nxi2,:,:)=pi5(nxi2-1,:,:)
   !pi5(:,nyi1,:)=pi5(:,2,:)
   !pi5(:,1,:)=pi5(:,nyi1-1,:)


   call write_restart(p_filename_out,nxi2,nyi1,nzi2,pi5,time,dtm1,grav,jp,nstep)


   deallocate(pi5)

   endif





   ! ------------- D  ----------------


   if (id==1) then

   write(*,*) 'Interpolating and imposing the BC of: Density'

   allocate(d(nx,ny,nz))
   allocate(di1(nxi1,ny,nz))
   allocate(di2(nxi2,ny,nz))
   allocate(di3(nxi2,nyi1,nz))
   allocate(di4(nxi2,nyi1,nzi1))
   allocate(di5(nxi2,nyi1,nzi2))


   d=0.0d0;di1=0.0d0;di2=0.0d0;di3=0.0d0;di4=0.0d0;di5=0.0d0


   call read_restart(d_filename,nx,ny,nz,d,time,dtm1,grav,jp,nstep)


   !d(:,:,1)=d(:,:,2)
   !d(:,:,nz)=d(:,:,nz-1)

   ! -- 1st X interpolation --

!   do  i=1,nx-1
!      di1(2*i-1,:,:)=0.75*d(i,:,:)+0.25*d(i+1,:,:)
!   enddo
!
!   do  i=2,nx
!      di1(2*i-2,:,:)=0.25*d(i-1,:,:)+0.75*d(i,:,:)
!   enddo

   di1=d

   deallocate(d)

   ! -- 2nd X interpolation --

   di2=di1

   ! do  i=1,nxi1-1
   !    di2(2*i-1,:,:)=0.75*di1(i,:,:)+0.25*di1(i+1,:,:)
   ! enddo
   !
   ! do  i=2,nxi1
   !    di2(2*i-2,:,:)=0.25*di1(i-1,:,:)+0.75*di1(i,:,:)
   ! enddo


   deallocate(di1)

   ! -- Y interpolation --

   !di2(:,1,:)=di2(:,ny-1,:)
   !di2(:,2,:)=di2(:,ny,:)

   do  j=1,ny-1
      di3(:,2*j-1,:)=0.75*di2(:,j,:)+0.25*di2(:,j+1,:)
   enddo

   do  j=2,ny
      di3(:,2*j-2,:)=0.25*di2(:,j-1,:)+0.75*di2(:,j,:)
   enddo

!   di3=di2

   deallocate(di2)

   ! -- 1st Z interpolation --

!   do  k=1,nz-1
!      di4(:,:,2*k-1)=0.75*di3(:,:,k)+0.25*di3(:,:,k+1)
!   enddo
!
!   do  k=2,nz
!      di4(:,:,2*k-2)=0.25*di3(:,:,k-1)+0.75*di3(:,:,k)
!   enddo


!new interp
   do  k=1,nz-2
      di4(:,:,2*k)=di3(:,:,k+1)
      di4(:,:,2*k+1)=di3(:,:,k+1)
   enddo

  di4(:,:,1)=di3(:,:,1)
  di4(:,:,nzi1)=di3(:,:,nz)

   !di4=di3
   
   deallocate(di3)

   ! -- 2nd Z interpolation --


!   do  k=1,nzi1-1
!      di5(:,:,2*k-1)=0.75*di4(:,:,k)+0.25*di4(:,:,k+1)
!   enddo
!
!   do  k=2,nzi1
!      di5(:,:,2*k-2)=0.25*di4(:,:,k-1)+0.75*di4(:,:,k)
!   enddo

   do  k=1,nzi1-2
      di5(:,:,2*k)=di4(:,:,k+1)
      di5(:,:,2*k+1)=di4(:,:,k+1)
   enddo

  di5(:,:,1)=di4(:,:,1)
  di5(:,:,nzi2)=di4(:,:,nzi1)


!   di5 = di4

   deallocate(di4)

   !boundary conditions
   !di5(:,nyi1,:)=di5(:,2,:)
   !di5(:,1,:)=di5(:,nyi1-1,:)
   !di5(1:nxi2-1,:,nzi2)=0.0
   !di5(1:nxi2-1,:,1)=0.0

   call write_restart(d_filename_out,nxi2,nyi1,nzi2,di5,time,dtm1,grav,jp,nstep)


   deallocate(di5)

   endif




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!      write(*,*)'maxval(u):',maxval(u),maxloc(u)
!      write(*,*)'maxval(w):',maxval(w),maxloc(w)
!      write(*,*)'maxval(v):',maxval(v),maxloc(v)
!      write(*,*)'maxval(p):',maxval(p),maxloc(p)
!      write(*,*)'maxval(d):',maxval(d),maxloc(d)

!      write(*,*)'maxval(ui1):',maxval(ui1),maxloc(ui1)
!      write(*,*)'maxval(vi1):',maxval(vi1),maxloc(vi1)
!      write(*,*)'maxval(wi1):',maxval(wi1),maxloc(wi1)




!      write(*,*)'maxval(ui2):',maxval(ui2),maxloc(ui2)
!      write(*,*)'maxval(vi2):',maxval(vi2),maxloc(vi2)
!      write(*,*)'maxval(wi2):',maxval(wi2),maxloc(wi2)


!      write(*,*)'maxval(ui4):',maxval(ui4),maxloc(ui4)
!      write(*,*)'maxval(vi4):',maxval(vi4),maxloc(vi4)
!      write(*,*)'maxval(wi4):',maxval(wi4),maxloc(wi4)



!      write(*,*)'maxval(ui5):',maxval(ui5),maxloc(ui5)
!      write(*,*)'maxval(vi5):',maxval(vi5),maxloc(vi5)
!      write(*,*)'maxval(wi5):',maxval(wi5),maxloc(wi5)
!      write(*,*)'maxval(pi5):',maxval(pi5),maxloc(pi5)
!      write(*,*)'maxval(di5):',maxval(di5),maxloc(di5)

!  call write_restart(v_filename_out,nxi2,nyi1,nzi2,vi5,time,dtm1,grav,jp,nstep)
!  call write_restart(w_filename_out,nxi2,nyi1,nzi2,wi5,time,dtm1,grav,jp,nstep)
!  call write_restart(p_filename_out,nxi2,nyi1,nzi2,pi5,time,dtm1,grav,jp,nstep)
!  call write_restart(d_filename_out,nxi2,nyi1,nzi2,di5,time,dtm1,grav,jp,nstep)


  stop
end Program Main

subroutine read_grid(xu,yv,zw,zwg,xc,yc,zc,zcg,nx,ny,nz,nzg,ru,rp,tag)
  implicit none

  INTEGER nx,ny,nz,nzg,tag
  REAL ( kind = 8 ) :: xu(nx),yv(ny),zw(nz),zwg(nzg)
  REAL ( kind = 8 ) :: xc(nx),yc(ny),zc(nz),zcg(nzg)
  REAL ( kind = 8 ) :: ru(nx),rp(nx)

  INTEGER i,j,k

  real rdelx,rdely,rdelz, dtheta
  real, allocatable, dimension(:) :: cug,cvg

  ALLOCATE(cug(nzg),cvg(nzg))

  ! ! READ GRID

  OPEN(UNIT=1,FILE='./x1_grid.in',STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nx-1
     read(1,*) j, xu(i)
     !write(6,*) "xu(",i,") = ", xu(i)
  enddo
  close(1)
  xc(2:nx-1) = .5*(xu(1:nx-2)+xu(2:nx-1))
  xc(1 ) = 2.*xu(1  )-xc(2  )
  xc(nx) = 2.*xu(nx-1)-xc(nx-1)

  do i= 2, nx-1
     write(6,*) "xc(",i,") = ", xc(i)
  enddo


  OPEN(UNIT=1,FILE='./x2_grid.in',STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, ny-1
     read(1,*) j, yv(i)
     !write(6,*) "yv(",i,") = ", yv(i)
  enddo
  close(1)
  yc(2:ny-1) = .5*(yv(1:ny-2)+yv(2:ny-1))
  yc(1 ) = 2.*yv(1  )-yc(2  )
  yc(ny) = 2.*yv(ny-1)-yc(ny-1)

  do i= 2, ny-1
     write(6,*) "yc(",i,") = ", yc(i)
  enddo


  OPEN(UNIT=1,FILE='./x3_grid.in',STATUS='OLD',FORM='FORMATTED')
  read(1,*) j
  do i= 1, nz-1
     read(1,*) j, zwg(i)
     !write(6,*) "zwg(",i,") = ", zwg(i)
  enddo
  close(1)
  zcg(2:nz-1) = .5*(zwg(1:nz-2)+zwg(2:nz-1))
  zcg(1  )= 2.*zwg(1  )-zcg(2  )
  zcg(nz) = 2.*zwg(nz-1)-zcg(nz-1)

  do i= 2, nz-1
     write(6,*) "zcg(",i,") = ", zcg(i)
  enddo


  ! ru(1:nx)=xu(1:nx)
  ! rp(1:nx)=xc(1:nx)

  close(1)

  write(6,*) "Read Grid: done"

  return
end subroutine read_grid

subroutine read_restart(filename,nx,ny,nz,var,time,dtm1,grav,jp,nstep)
  implicit none

  character(len=128)   :: filename
  integer ( kind = 4 ) :: i,j,k,jp,nx,ny,nz,nstep
  real    ( kind = 8 ) :: var(nx,ny,nz),time,DTM1,grav

  !write(6,*) "nx,ny,nz = ", nx,ny,nz

  ! READ RESTART FILE
  OPEN(19,FILE=filename,STATUS='UNKNOWN',FORM='UNFORMATTED')
  READ(19) I,J,K,JP
  !write(6,*) "I,J,K,JP = ", I,J,K,JP
  DO K=1,NZ
  !   write(6,*) " READ K = ", K
     READ(19) ((var(I,J,K),I=1,NX),J=1,NY)
  ENDDO
  READ(19) nstep
  READ(19) TIME
  !write(6,*) 'time=',time
  READ(19) DTM1,grav
  CLOSE(19)
  !write(6,*) "Read ",filename, " restartfile: done"

  return
end subroutine read_restart


subroutine write_restart(filename,nx,ny,nz,var,time,dtm1,grav,jp,nstep)
  implicit none

  character(len=128)   :: filename
  integer ( kind = 4 ) :: i,j,k,jp,nx,ny,nz,nstep
  real    ( kind = 8 ) :: var(nx,ny,nz),time,DTM1,grav

  write(6,*) "nx,ny,nz = ", nx,ny,nz

  ! WRITE RESTART FILE
  OPEN(20,FILE=trim(filename),STATUS='UNKNOWN',FORM='UNFORMATTED')
  i=nx; j=ny; k=nz
  WRITE(20) i,j,k,JP
  !write(6,*) "I,J,K,JP = ", I,J,K,JP
  DO K=1,NZ
  !   write(6,*) " WRITE K = ", K
     write(20) ((var(I,J,K),I=1,NX),J=1,NY)
  ENDDO
  WRITE(20) nstep
  WRITE(20) TIME
  !WRITE(6,*) 'time=',time
  WRITE(20) DTM1,grav
  CLOSE(20)
  write(6,*) "Write ", filename," restartfile: done"

  return
end subroutine write_restart
