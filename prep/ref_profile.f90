PROGRAM inv_prep
  
  !     ********************************************************************************
  !     * CALCULATE REFERENCE PROFILE, CORIOLIS PARAMETER AND GRID PARAMETERS          *
  !     * Rene Fehlmann 1994 / Code re-organization: Michael Sprenger, 2006            *
  !     * Modified by Ashwin Dinakar 2017 & 2018 conversion to f90, use of modues and *
  !     * Fortran netcdf
  !     ********************************************************************************

  !     --------------------------------------------------------------------------------
  !     Declaration of variables, parameters, externals and common blocks
  !     --------------------------------------------------------------------------------
  use kind_parameters,ONLY:&
       sp
  use netcdflibrary
  use cdfreadwrite
  implicit none
  
  !     Input and output file
  character*80   pvsrcfile,referfile
  integer        mode
  real(kind=sp)           radius
  
  !     Grid parameters
  integer        nx,ny,nz
  real(kind=sp)           xmin,ymin,zmin,xmax,ymax,zmax
  real(kind=sp)           dx,dy,dz
  real(kind=sp)           mdv
  
!     Reference state
  real(kind=sp),   allocatable,dimension (:) :: nsqref
  real(kind=sp),   allocatable,dimension (:) :: thetaref
  real(kind=sp),   allocatable,dimension (:) :: rhoref
  real(kind=sp),   allocatable,dimension (:) :: pressref
  real(kind=sp),   allocatable,dimension (:) :: zref
  real(kind=sp)    pressn,thetan

  !  3d fields for calculation of reference profile
  real(kind=sp),allocatable,dimension (:,:,:) :: th,rho,nsq,p,z
  
  !  2d weight for mean
  real(kind=sp),allocatable,dimension (:,:) :: weight
  
!  Auxiliary variables
  integer      i,j,k
  integer      stat
  integer      jj,kk
  character*80 varname
  integer      istep
  real(kind=sp)         sum,max,min
  integer      cnt,nmd
  integer      step,ntimes
  integer      i0,j0
  real(kind=sp)         lon0,lat0,lon1,lat1,dist
  integer      vardim(4),ndim,cdfid,ierr
!  real(kind=sp)         varmin(4),varmax(4),stag(4)
  real(kind=sp),   allocatable,dimension (:) :: varmin,varmax,stag
  real(kind=sp),dimension(:),allocatable::time
  !   real(kind=sp)         mdv
  character*80 name
  
!  Externals
  real(kind=sp)      sdis
  external  sdis
  
!  --------------------------------------------------------------------------------
!  Input
!  --------------------------------------------------------------------------------
  allocate(time(1))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  print*,'********************************************************'
  print*,'* ref_grid                                             *'
  print*,'********************************************************'
  
  !  Read parameter file
  open(10,file='fort.10')
  read(10,*) pvsrcfile
  read(10,*) referfile
  read(10,*) name,radius
  if ( trim(name).ne.'REF_R') stop
  close(10) 
  print*
  print*,trim(pvsrcfile)
  print*,trim(referfile)
  print*,radius
  print*


  call cdfopn(referfile,cdfid)
  call gettimes(cdfid,time,ntimes)
  call clscdf(cdfid)
!  Get lat/lon gid parameters from input file
  call read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,&
       pvsrcfile)
  print*,'Read_Dim: nx,ny,nz         = ',nx,ny,nz
  print*,'          dx,dy,dz         = ',dx,dy,dz
  print*,'          xmin,ymin,zmin   = ',xmin,ymin,zmin
  print*,'          mdv              = ',mdv
  print*
  xmax = xmin + real(nx-1) * dx
  ymax = ymin + real(ny-1) * dy
  
  !  Count from 0, not from 1: consistent with <inv_cart.f>.
  nx=nx-1
  ny=ny-1
  nz=nz-1

  !  Allocate memory for reference profile
  allocate(rhoref  (0:2*nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating rhoref'
  allocate(pressref(0:2*nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating pressref'
  allocate(thetaref(0:2*nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating thetaref'
  allocate(nsqref  (0:2*nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating nsqref'
  allocate(zref    (0:2*nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating zref'
  
  !  Allocate memory for calculatation of reference profile
  allocate(th (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating th'
  allocate(rho(0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating rho'
  allocate(p  (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating p'
  allocate(nsq(0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating nsq'
  allocate(z(0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating z'
  
  !  Allocate memory for weight
  allocate(weight(0:nx,0:ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating weight'
  
  !  --------------------------------------------------------------------------------
!  Calculate the reference profile and put it onto netcdf file
!  --------------------------------------------------------------------------------

!  Read data from file
  varname='TH'
  call read_inp (th,varname,pvsrcfile,&
       nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='RHO'
  call read_inp (rho,varname,pvsrcfile,&
       nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='NSQ'
  call read_inp (nsq,varname,pvsrcfile,&
       nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='P'
  call read_inp (p,varname,pvsrcfile,&
       nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  
  !  Init the height field (not real(kind=sp)ly necessary, but code becomes more symmetrical)
  do i=0,nx
     do j=0,ny
        do k=0,nz
           z(i,j,k)=zmin+real(k)*dz
        enddo
     enddo
  enddo
  
  !  Define the weight
  if ( radius.lt.0 ) then
     do i=0,nx
        do j=0,ny
           weight(i,j) = 1.
        enddo
     enddo
  else
     i0 = nx/2
     j0 = ny/2
     lon0 = xmin + real(i0-1) * dx
     lat0 = ymin + real(j0-1) * dy
     weight(i0,j0)=1.
     do i=0,nx
        do j=0,ny
           lon1 = xmin + real(i-1) * dx
           lat1 = ymin + real(j-1) * dy
           dist = sdis(lon0,lat0,lon1,lat1)
           if ( dist.lt.radius ) then
              weight(i,j) = 1.
           else
              weight(i,j) = 0.
           endif
        enddo
     enddo
  endif

  !  Determine the reference profile (mean over domain, split levels and layers)

  call mean(zref,    z,  nx,ny,nz,mdv,weight)
  call mean(nsqref,  nsq,nx,ny,nz,mdv,weight)
  call mean(rhoref,  rho,nx,ny,nz,mdv,weight)
  call mean(thetaref,th, nx,ny,nz,mdv,weight)
  call mean(pressref,p,  nx,ny,nz,mdv,weight)
  
  !  Write reference file to netcdf file
  call write_ref (nsqref,rhoref,thetaref,pressref,zref,&
       nz,referfile)      
  
  !  Write some info
  print*
  print*,'Ref:            z         p       rho       nsq     theta'
  step=2*nz/10
  if (step.lt.1) step=1
  do k=0,2*nz,step
     write(*,'(8x,f10.1,2f10.2,f10.6,f10.2)')&
          zref(k),pressref(k),rhoref(k),nsqref(k),thetaref(k) 
  enddo
  print*
  
  !  Write weighht to REF file
  call cdfwopn(trim(referfile),cdfid)
  varname='WEIGHT'
  vardim(1)=nx+1
  vardim(2)=ny+1
  vardim(3)=1
  vardim(4)=1
  varmin(1)=xmin
  varmin(2)=ymin
  varmin(3)=0.
  varmax(1)=xmax
  varmax(2)=ymax
  varmax(3)=0.
  ndim=3
  mdv=-999.

  call putdef(cdfid,varname,ndim,mdv,vardim,&
       varmin,varmax,stag)

  call putdatRank2(cdfid,varname,time,1,weight)
  
  !  --------------------------------------------------------------------------------
  !  Format specifications
  !  --------------------------------------------------------------------------------
  
111 format (5f20.9)
106 format (2f20.9)
  
end PROGRAM inv_prep

!  --------------------------------------------------------------------------
!  Spherical distance between lat/lon points
!  --------------------------------------------------------------------------

real(kind=sp) function sdis(xp,yp,xq,yq)
  use kind_parameters,ONLY:&
       sp
  
  !  calculates spherical distance (in km) between two points given
  !  by their spherical coordinates (xp,yp) and (xq,yq), respectively.
  
  real(kind=sp)      re
  parameter (re=6370.)
  real(kind=sp)      pi180
  parameter (pi180=3.14159/180.)
  real(kind=sp)      xp,yp,xq,yq,arg
  
  arg=sin(pi180*yp)*sin(pi180*yq)+&
       cos(pi180*yp)*cos(pi180*yq)*cos(pi180*(xp-xq))
  if (arg.lt.-1.) arg=-1.
  if (arg.gt.1.) arg=1.
  
  sdis=re*acos(arg)
  
end function sdis
      
!  ********************************************************************************
!  * NETCDF INPUT AND OUTPUT                                                      *
!  ********************************************************************************


!  --------------------------------------------------------------------------------
!  Write reference file to netcdf
!  --------------------------------------------------------------------------------



!  --------------------------------------------------------------------------------
!  Read input fields for reference profile
!  --------------------------------------------------------------------------------


!  --------------------------------------------------------------------------------
!  Check whether variable is found on netcdf file
!  --------------------------------------------------------------------------------

subroutine check_varok (isok,varname,varlist,nvars)
  
  !  Check whether the variable <varname> is in the list <varlist(nvars)>. If this is 
  !  the case, <isok> is incremented by 1. Otherwise <isok> keeps its value.
  
  implicit none
  
  !  Declaraion of subroutine parameters
  integer      isok
  integer      nvars
  character*80 varname
  character*80 varlist(nvars)
  
  !  Auxiliary variables
  integer      i
  
  !  Main
  do i=1,nvars
     if (trim(varname).eq.trim(varlist(i))) isok=isok+1
  enddo
  
end subroutine check_varok

!  --------------------------------------------------------------------------------
!  Get grid parameters
!  --------------------------------------------------------------------------------

subroutine read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,&
     pvsrcfile)
  
  !  Get the grid parameters from the variable <TH> on the input file <pvsrcfile>.
  !  The grid parameters are
  !     nx,ny,nz                : Number of grid points in x, y and z direction
  !     xmin,ymin,zmin          : Minimum domain coordinates in x, y and z direction
  !     xmax,ymax,zmax          : Maximal domain coordinates in x, y and z direction
  !     dx,dy,dz                : Horizontal and vertical resolution
  !  Additionally, it is checked whether the vertical grid is equally spaced. If ok,
  !  the grid paramters are transformed from lon/lat to distance (in meters)
  use kind_parameters,ONLY:&
       sp

  use netcdflibrary
  implicit none
  
  !  Declaration of subroutine parameters
  character*80   pvsrcfile
  integer        nx,ny,nz
  real(kind=sp)           dx,dy,dz
  real(kind=sp)           xmin,ymin,zmin,xmax,ymax,zmax
  real(kind=sp)           mdv
  
  !  Numerical epsilon and other physical/geoemtrical parameters
  real(kind=sp)           eps
  parameter      (eps=0.01)
  
  !  Auxiliary variables
  integer        cdfid,cstid
  integer        ierr
  character*80   vnam(100),varname
  integer        nvars
  integer        isok
  integer        vardim(4)
  real(kind=sp)           misdat
!  real(kind=sp)           varmin(4),varmax(4),stag(4)
  real(kind=sp),   allocatable,dimension (:) :: varmin,varmax,stag
  real(kind=sp)           aklev(1000),bklev(1000),aklay(1000),bklay(1000)
  real(kind=sp)           dh
  character*80   csn
  integer        ndim
  integer        i

  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  !  Get all grid parameters
  call cdfopn(pvsrcfile,cdfid)

  call getvars(cdfid,nvars,vnam)

  isok=0
  varname='T'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 998
  call getcfn(cdfid,csn)

  call cdfopn(csn,cstid)

  call getdef(cdfid,varname,ndim,misdat,vardim,varmin,varmax,&
       stag)

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

  
  !  Check whether the grid is equallay spaced in the vertical
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
  
  !  Set missing data value
  mdv=misdat
  
  return
  
  !  Exception handling
 998  print*,'Read_Dim: Problem with input netcdf file... Stop'
  stop
  
end subroutine read_dim


!  ********************************************************************************
!  * DEFINE REFERENCE PROFILE                                                     *
!  ********************************************************************************

!  --------------------------------------------------------------------------------
!  Calculate area mean
!  --------------------------------------------------------------------------------

SUBROUTINE mean(a,m,nx,ny,nz,mdv,weight)
  
  !  Calculate the area-mean of <m> and save it on <a>.
  use kind_parameters,ONLY:&
       sp

  implicit none
  
  !  Declaration of subroutine parameters
  real(kind=sp)       mdv
  integer    nx,ny,nz
  real(kind=sp)       m(0:nx,0:ny,0:nz),a(0:2*nz)
  real(kind=sp)       weight(0:nx,0:ny)
  
  !  Numerical epsilon
  real(kind=sp)       eps
  parameter  (eps=0.01)
  
  !  Auxiliary varaibles
  real(kind=sp)       mea(0:nz)
  integer    i,j,k,kk
  real(kind=sp)       counter
  
  !  Determine the mean over all levels (handle missing data)
  do k=0,nz
     mea(k)=0.
     counter=0.
     do i=0,nx
        do j=0,ny
           if (abs(m(i,j,k)-mdv).gt.eps) then
              mea(k)=mea(k)+m(i,j,k)*weight(i,j)
              counter=counter+weight(i,j)
           endif
           
        enddo
     enddo
     if (counter.gt.0) then
        mea(k)=mea(k)/counter
     else
        mea(k)=mdv
     endif
  enddo
  
  !  Prepare the output array: split layers and levels
  do k=0,nz-1
     kk=2*k
     a(kk)=mea(k)
     a(kk+1)=0.5*(mea(k)+mea(k+1))
  enddo
  a(2*nz)=mea(nz)
      
end SUBROUTINE mean
