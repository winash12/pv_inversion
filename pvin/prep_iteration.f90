PROGRAM prep_iteration

!     ************************************************************************
!     * Prepare the next step for the PV inversion                           *
!     * Michael Sprenger / Summer, Autumn 2006                               *
!     ************************************************************************

  use kind_parameters,ONLY:&
       sp
  use netcdflibrary
  use cdfreadwrite
  implicit none

!    ------------------------------------------------------------------------
!     Declaration of variables and parameters
!     ------------------------------------------------------------------------
      
!     Input and output file
  character*80   anomafile
  character*80   iterafile
  
  !    Grid parameters
  integer        nx,ny,nz
  real(kind=sp)           xmin,ymin,zmin,xmax,ymax,zmax
  real(kind=sp)           dx,dy,dz
  real(kind=sp)           mdv

!    Numerical epsilon and other variables
  real(kind=sp)           eps
  parameter      (eps=0.01)
  real(kind=sp)           alpha
  
!    3d arrays
  real(kind=sp),allocatable,dimension (:,:,:) :: v_iter,v_anom
  real(kind=sp),allocatable,dimension (:,:,:) :: u_iter,u_anom
  real(kind=sp),allocatable,dimension (:,:,:) :: t_iter,t_anom
  real(kind=sp),allocatable,dimension (:,:,:) :: p_iter,p_anom

!    Auciliary variables 
  integer      i,j,k
  integer      stat
  character*80 varname
  character*80 name

!    --------------------------------------------------------------------------------
!     Input
!     --------------------------------------------------------------------------------

  print*,'********************************************************'
  print*,'* PREP_ITERATION                                       *'
  print*,'********************************************************'

!    Read parameter file
  open(10,file='fort.10')
  read(10,*) iterafile
  read(10,*) anomafile
  read(10,*) name,alpha
  if ( trim(name).ne.'ALPHA' ) stop
  close(10) 
  print*
  print*,trim(anomafile)
  print*,trim(iterafile)
  print*
  
!    Get lat/lon gid parameters from input file
  call read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,anomafile)
  print*,'Read_Dim: nx,ny,nz         = ',nx,ny,nz
  print*,'          dx,dy,dz         = ',dx,dy,dz
  print*,'          xmin,ymin,zmin   = ',xmin,ymin,zmin
  print*,'          mdv              = ',mdv
  print*
  
!     Count from 0, not from 1: consistent with <inv_cart.f>.
  nx=nx-1
  ny=ny-1
  nz=nz-1
  
!     Allocate memory for 3d arrays 
  allocate(u_anom (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating u_anom'
  allocate(v_anom (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating v_anom'
  allocate(t_anom  (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating t_anom'
  allocate(p_anom  (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating p_anom'
  
  allocate(u_iter (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating u_iter'
  allocate(v_iter (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating v_iter'
  allocate(t_iter  (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating t_iter'
  allocate(p_iter  (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating p_iter'
  
!     Read anomaly and iteration fields
  varname='U'
  call read_inp_prep(u_anom,varname,anomafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='V'
  call read_inp_prep(v_anom,varname,anomafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='T'
  call read_inp_prep(t_anom,varname,anomafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='P'
  call read_inp_prep(p_anom,varname,anomafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='U'
  call read_inp_prep(u_iter,varname,iterafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='V'
  call read_inp_prep(v_iter,varname,iterafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='T'
  call read_inp_prep(t_iter,varname,iterafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='P'
  call read_inp_prep(p_iter,varname,iterafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)

!     --------------------------------------------------------------------------------
!     Modify the iteration fields
!     --------------------------------------------------------------------------------

  do i=0,nx
     do j=0,ny
        do k=0,nz
           
!              Update zonal velocity
           if ((abs(u_anom(i,j,k)-mdv).gt.eps).and. &
                (abs(u_iter(i,j,k)-mdv).gt.eps)) then
              u_iter(i,j,k)=u_iter(i,j,k)-alpha*u_anom(i,j,k)
           else
              u_iter(i,j,k)=mdv
           endif

           !              Update meridional velocity
           if ((abs(v_anom(i,j,k)-mdv).gt.eps).and.&
                (abs(v_iter(i,j,k)-mdv).gt.eps)) then
              v_iter(i,j,k)=v_iter(i,j,k)-alpha*v_anom(i,j,k)
           else
              v_iter(i,j,k)=mdv
           endif

!              Update temperature
           if ((abs(t_anom(i,j,k)-mdv).gt.eps).and. &
                (abs(t_iter(i,j,k)-mdv).gt.eps)) then
              t_iter(i,j,k)=t_iter(i,j,k)-alpha*t_anom(i,j,k)
           else
              t_iter(i,j,k)=mdv
           endif

!              Update pressure
           if ((abs(p_anom(i,j,k)-mdv).gt.eps).and. &
                (abs(p_iter(i,j,k)-mdv).gt.eps)) then
              p_iter(i,j,k)=p_iter(i,j,k)-alpha*p_anom(i,j,k)
           else
              p_iter(i,j,k)=mdv
           endif
           
        enddo
     enddo
  enddo


!     --------------------------------------------------------------------------------
!     Write output
!     --------------------------------------------------------------------------------

  varname='U'
  call write_inp (u_iter,varname,iterafile,nx,ny,nz)
  varname='V'
  call write_inp (v_iter,varname,iterafile,nx,ny,nz)
  varname='T'
  call write_inp (t_iter,varname,iterafile,nx,ny,nz)
  varname='P'
  call write_inp (p_iter,varname,iterafile,nx,ny,nz)
      
      
end PROGRAM prep_iteration




!     ********************************************************************************
!     * NETCDF INPUT AND OUTPUT                                                      *
!     ********************************************************************************


!     --------------------------------------------------------------------------------
!     Read input fields for reference profile
!     --------------------------------------------------------------------------------

!     --------------------------------------------------------------------------------
!     Check whether variable is found on netcdf file
!     --------------------------------------------------------------------------------

subroutine check_varok (isok,varname,varlist,nvars)

!     Check whether the variable <varname> is in the list <varlist(nvars)>. If this is 
!     the case, <isok> is incremented by 1. Otherwise <isok> keeps its value.

  implicit none
  
!    Declaraion of subroutine parameters
  integer      isok
  integer      nvars
  character*80 varname
  character*80 varlist(nvars)

!    Auxiliary variables
  integer      i
  
  !    Main
  do i=1,nvars
     if (trim(varname).eq.trim(varlist(i))) isok=isok+1
  enddo
  
end subroutine check_varok

!     --------------------------------------------------------------------------------
!     Get grid parameters
!     --------------------------------------------------------------------------------

subroutine read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,pvsrcfile)

  !     Get the grid parameters from the variable <THETA> on the input file <pvsrcfile>.
  !     The grid parameters are
  !        nx,ny,nz                : Number of grid points in x, y and z direction
  !        xmin,ymin,zmin          : Minimum domain coordinates in x, y and z direction
  !        xmax,ymax,zmax          : Maximal domain coordinates in x, y and z direction
  !dx,dy,dz                : Horizontal and vertical resolution
!     Additionally, it is checked whether the vertical grid is equally spaced. If ok,
!     the grid paramters are transformed from lon/lat to distance (in meters)
  use netcdflibrary
  implicit none
  
  !    Declaration of subroutine parameters
  character*80   pvsrcfile
  integer        nx,ny,nz
  real(kind=sp)           dx,dy,dz
  real(kind=sp)           xmin,ymin,zmin,xmax,ymax,zmax
  real(kind=sp)           mdv

!     Numerical epsilon and other physical/geoemtrical parameters
  real(kind=sp)           eps
  parameter      (eps=0.01)
  
!     Auxiliary variables
  integer        cdfid,cstid
  integer        ierr
  character*80   vnam(100),varname
  integer        nvars
  integer        isok
  integer        vardim(4)
  real(kind=sp)           misdat
!  real(kind=sp)           varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  real(kind=sp)           aklev(1000),bklev(1000),aklay(1000),bklay(1000)
  real(kind=sp)           dh
  character*80   csn
  integer        ndim
  integer        i

!    Get all grid parameters
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  call cdfopn(pvsrcfile,cdfid)
  call getvars(cdfid,nvars,vnam)
  isok=0
  varname='TH'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 998
  call getcfn(cdfid,csn)
  call cdfopn(csn,cstid)
  call getdef(cdfid,varname,ndim,misdat,vardim,varmin,varmax,stag)
  nx=vardim(1)
  ny=vardim(2)
  nz=vardim(3)
  xmin=varmin(1)
  ymin=varmin(2)
  zmin=varmin(3)
  call getlevs(cstid,nz,aklev,bklev,aklay,bklay)
  call getgrid(cstid,dx,dy)
  xmax=varmax(1)
  ymax=varmax(2)
  zmax=varmax(3)
  dz=(zmax-zmin)/real(nz-1)
  call clscdf(cstid)
  call clscdf(cdfid)
  
!    Check whether the grid is equallay spaced in the vertical
  do i=1,nz-1
     dh=aklev(i+1)-aklev(i)
     if (abs(dh-dz).gt.eps) then
        print*,'Aklev: Vertical grid must be equally spaced... Stop'
        print*,(aklev(1:nz))
        stop
     endif
     dh=aklay(i+1)-aklay(i)
     if (abs(dh-dz).gt.eps) then
        print*,'Aklay: Vertical grid must be equally spaced... Stop'
        print*,(aklay(1:nz))
        stop
     endif
  enddo
  
  !     Set missing data value
  mdv=misdat
  
  return
  
  !     Exception handling
998 print*,'Read_Dim: Problem with input netcdf file... Stop'
  stop
  
end subroutine read_dim


