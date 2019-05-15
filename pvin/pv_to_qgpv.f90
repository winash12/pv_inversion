PROGRAM pv_to_qgpv

  !     ********************************************************************************
  !     * TRANSFORM ERTEL'S PV TO QUASI-GEOSTROPHIC PV                                 *
  !     * Rene Fehlmann 1994 / Code re-organization: Michael Sprenger, 2006            *
  !     ********************************************************************************
  !     * Conversion to f90, netcdf fortran by Ashwin Dinakar 2017
  !     --------------------------------------------------------------------------------
  !     Declaration of variables, parameters, externals and common blocks
  !     --------------------------------------------------------------------------------
  use netcdflibrary
  use cdfreadwrite
  implicit none
  
!     Input and output file
  character*80   pvsrcfile
  character*80   referfile
  character*80   anomafile
  
!     Grid parameters
  integer        nx,ny,nz
  real           xmin,ymin,zmin,xmax,ymax,zmax
  real           dx,dy,dz
  real           mdv

!    Numerical and physical parameters
  real           pi180                           ! Pi/180
  parameter      (pi180=3.141592654/180.)
  real           rerd                            ! Earth's radius
  parameter      (rerd=6.371229e6)
  real           eps                             ! Numerical epsilon
  parameter      (eps=0.01)
  real           scale                           ! Scale for PV unit
  parameter      (scale=1e6)
  real           minagl                          ! No PV and qgPV below this height AGL
  parameter      (minagl=1000.)
      
!     Reference state and grid parameters
  real,   allocatable,dimension (:)   :: nsqref
  real,   allocatable,dimension (:)   :: thetaref
  real,   allocatable,dimension (:)   :: rhoref
  real,   allocatable,dimension (:)   :: pressref
  real,   allocatable,dimension (:)   :: zref
  real,   allocatable,dimension (:,:) :: coriol
  real,   allocatable,dimension (:,:) :: oro
  real    deltax,deltay,deltaz

!     3d fields for calculation of qgPV and Ertel's PV
  real,allocatable,dimension (:,:,:) :: qgpv,pv1,pv2,pv
  
  !   Auxiliary variables
  real         zpos
  integer      i,j,k
  integer      stat
  character*80 varname
  integer      istep
  real         mean,rmsq,min,max
  integer      step
  real,allocatable,dimension (:,:) :: tmp
  
!   --------------------------------------------------------------------------------
!   Input
!   --------------------------------------------------------------------------------

  print*,'********************************************************'
  print*,'* PV_TO_QGPV                                           *'
  print*,'********************************************************'
  
  
  !   Read parameter file
  open(10,file='fort.10')
  read(10,*) pvsrcfile
  read(10,*) referfile
  read(10,*) anomafile
  close(10) 
  print*
  print*,trim(pvsrcfile)
  print*,trim(referfile)
  print*,trim(anomafile)
  print*
  
  !   Get lat/lon gid parameters from input file
  call read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,pvsrcfile)
  print*,'Read_Dim: nx,ny,nz         = ',nx,ny,nz
  print*,'          dx,dy,dz         = ',dx,dy,dz
  print*,'          xmin,ymin,zmin   = ',xmin,ymin,zmin
  print*,'          mdv              = ',mdv
  print*
  
!   Count from 0, not from 1: consistent with <inv_cart.f>.
  nx=nx-1
  ny=ny-1
  nz=nz-1

!   Allocate memory for reference profile and grid parameters
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
  allocate(coriol  (0:nx,0:ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating coriol'
  allocate(oro     (0:nx,0:ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating oro'
  
  
  !   Allocate memory for calculation of qgPV and Ertel's PV
  allocate(pv1 (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating pv1'
  allocate(pv2 (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating pv2'
  allocate(pv  (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating pv'      
  allocate(qgpv(0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating qgpv'
  
  !   Allocate memory for temporary array
  allocate(tmp(0:nx,0:ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating tmp'
  
!   --------------------------------------------------------------------------------
!   Calculate the qgPV from Ertel's PV and put it onto file
!   --------------------------------------------------------------------------------

!   Read data from file
  varname='PV'
  call read_inp_pv_to_qgpv(pv1,varname,pvsrcfile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='PV_AIM'
  call read_inp_pv_to_qgpv(pv2,varname,pvsrcfile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
      
!   Read reference profile and grid parameters
  call read_ref_pv_to_qgpv(nsqref,rhoref,thetaref,pressref,zref,nx,ny,nz,deltax,deltay,deltaz,coriol,oro,referfile)

  !   If the PV is negative, set it to zero
  do i=0,nx
     do j=0,ny
        do k=0,nz
           if (pv1(i,j,k).lt.0.) pv1(i,j,k)=0.
           if (pv2(i,j,k).lt.0.) pv2(i,j,k)=0.
        enddo
     enddo
  enddo
  
!   Get the difference of Ertel's PV and set all missing values to 0
  do i=0,nx
     do j=0,ny
        do k=0,nz
           if ((abs(pv1(i,j,k)-mdv).gt.eps).and. &
                (abs(pv2(i,j,k)-mdv).gt.eps)) then
              pv(i,j,k)=pv1(i,j,k)-pv2(i,j,k)
           else
              pv(i,j,k)=0.
           endif
        enddo
     enddo
  enddo

!   Calculate qgPV
  call epv_to_qgpv (qgpv,pv,rhoref,pressref,nsqref,thetaref,nx,ny,nz,mdv)
  
!   Set values on the boundaries to zero
  do i=0,nx
     do j=0,ny
        qgpv(i,j, 0)=0.
        qgpv(i,j,nz)=0.
     enddo
  enddo
  do i=0,nx
     do k=0,nz
        qgpv(i, 0,k)=0.
        qgpv(i,ny,k)=0.
     enddo
  enddo
  do j=0,ny
     do k=0,nz
        qgpv( 0,j,k)=0.
        qgpv(nx,j,k)=0.
     enddo
  enddo

!   Set all values to zero which are too near to the surface
  do i=0,nx
     do j=0,ny
        do k=0,nz
           zpos=zmin+real(k)*dz
           if (zpos.lt.(oro(i,j)+minagl)) then
              pv(i,j,k)=0.
              qgpv(i,j,k)=0.
           endif
        enddo
     enddo
  enddo

!   Write result to netcdf file
  varname='QGPV'
  call write_inp_pv_to_qgpv(qgpv,varname,anomafile,nx,ny,nz)
  varname='PV'
  call write_inp_pv_to_qgpv(pv,varname,anomafile,nx,ny,nz)
  
  !   Write some info
  print*
  print*,'PV -> qgPV:     k    z     min     max    mean      rmsq'
  step=nz/10
  if (step.lt.1) step=1
  do k=0,nz,step
     do i=0,nx
        do j=0,ny
           tmp(i,j)=qgpv(i,j,k)
        enddo
     enddo
     call calc_error(min,max,mean,rmsq,tmp,nx,ny)
     write(*,'(8x,i3,f10.1,4f10.2)')&
          k,zmin+real(k)*dz,scale*min,scale*max,&
          scale*mean,scale*rmsq
  enddo
  print*         
  
  !   --------------------------------------------------------------------------------
  !   Format specifications
  !   --------------------------------------------------------------------------------
  
111 format (5f20.9)
106 format (2f20.9)
  
end PROGRAM pv_to_qgpv


!   ********************************************************************************
!   * NETCDF INPUT AND OUTPUT                                                      *
!   ********************************************************************************



!   --------------------------------------------------------------------------------
!   Read refernece profile from netcdf
!   --------------------------------------------------------------------------------

 
!   --------------------------------------------------------------------------------
!   Check whether variable is found on netcdf file
!   --------------------------------------------------------------------------------

subroutine check_varok (isok,varname,varlist,nvars)
  
  !   Check whether the variable <varname> is in the list <varlist(nvars)>. If this is 
  !   the case, <isok> is incremented by 1. Otherwise <isok> keeps its value.
  
  implicit none
  
  !   Declaraion of subroutine parameters
  integer      isok
  integer      nvars
  character*80 varname
  character*80 varlist(nvars)
  
  !   Auxiliary variables
  integer      i
  
  !   Main
  do i=1,nvars
     if (trim(varname).eq.trim(varlist(i))) isok=isok+1
  enddo
  
end subroutine check_varok

!   --------------------------------------------------------------------------------
!   Get grid parameters
!   --------------------------------------------------------------------------------

subroutine read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,pvsrcfile)

  !   Get the grid parameters from the variable <THETA> on the input file <pvsrcfile>.
  !   The grid parameters are
  !      nx,ny,nz                : Number of grid points in x, y and z direction
  !      xmin,ymin,zmin          : Minimum domain coordinates in x, y and z direction
  !      xmax,ymax,zmax          : Maximal domain coordinates in x, y and z direction
  !      dx,dy,dz                : Horizontal and vertical resolution
  !   Additionally, it is checked whether the vertical grid is equally spaced. If ok,
  !   the grid paramters are transformed from lon/lat to distance (in meters)
  use netcdflibrary
  implicit none
  
  !   Declaration of subroutine parameters
  character*80   pvsrcfile
  integer        nx,ny,nz
  real           dx,dy,dz
  real           xmin,ymin,zmin,xmax,ymax,zmax
  real           mdv

  !   Numerical epsilon and other physical/geoemtrical parameters
  real           eps
  parameter      (eps=0.01)
  
  !   Auxiliary variables
  integer        cdfid,cstid
  integer        ierr
  character*80   vnam(100),varname
  integer        nvars
  integer        isok
  integer        vardim(4)
  real           misdat
!  real           varmin(4),varmax(4),stag(4)
  real,dimension(:),allocatable::varmin,varmax,stag
  real           aklev(1000),bklev(1000),aklay(1000),bklay(1000)
  real           dh
  character*80   csn
  integer        ndim
  integer        i

  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  !   Get all grid parameters
  call cdfopn(pvsrcfile,cdfid)
  call getvars(cdfid,nvars,vnam)
  isok=0
  varname='PV'
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
  
  !   Check whether the grid is equallay spaced in the vertical
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
  
!   Set missing data value
  mdv=misdat
  
  return
  
  !   Exception handling
998 print*,'Read_Dim: Problem with input netcdf file... Stop'
  stop
  
end subroutine read_dim


!   ********************************************************************************
!   * CALCULATION                                                                  *
!   ********************************************************************************

!   --------------------------------------------------------------------------------
!   Calculate qgPV from Ertels's PV 
!   --------------------------------------------------------------------------------               

subroutine epv_to_qgpv (qgpv,pv,rhoref,pressref,nsqref,thetaref,nx,ny,nz,mdv)

!   Calculate the qgPV from Ertel's PV according to equation 2.11 p16, Thesis 
!   from Rene Fehlmann.

  implicit none
  
  !   Declaration of subroutine parameters
  integer   nx,ny,nz
  real      qgpv(0:nx,0:ny,0:nz),pv(0:nx,0:ny,0:nz)
  real      rhoref  (0:2*nz)
  real      nsqref  (0:2*nz)
  real      thetaref(0:2*nz)
  real      pressref(0:2*nz)
  real      mdv
  
!   Numerical epsilon
  real       g
  parameter  (g=9.81)
  real       eps
  parameter  (eps=0.01)
  real       scale
  parameter  (scale=1e6)
  
  !   Auxiliary variables
  integer i,j,k
  integer kk
  
  !   Calculation
  do i=0,nx
     do j=0,ny
        do k=0,nz
           
           kk=2*k
           
           if (( abs(rhoref(kk)  -mdv).gt.eps).and. &
                ( abs(thetaref(kk)-mdv).gt.eps).and. &
                ( abs(nsqref(kk)  -mdv).gt.eps).and. &
                ( abs(pv(i,j,k)   -mdv).gt.eps)) then
              
              qgpv(i,j,k)=rhoref(kk)*g*pv(i,j,k)/thetaref(kk)/ nsqref(kk)/scale
           else
              qgpv(i,j,k)=0.
           endif
           
        enddo
     enddo
  enddo
  
end subroutine epv_to_qgpv

!   --------------------------------------------------------------------------------
!   Calculate error statistics
!   --------------------------------------------------------------------------------               

subroutine calc_error (min,max,mean,rmsq,tmp,nx,ny)
  
  !   Calculate the error statistics for the two-dimensional error field <tmp>. The
  !   following error measures are calculated: the minimum <min>, the maximum <max>,
  !   the mean <mean>, the root-mean square <rmsq>
  
  implicit none
  
  !   Declaration of subroutine parameters
  integer nx,ny
  real    tmp(0:nx,0:ny)
  real    mean,rmsq
  real    min,max
  
  !   Auxiliary variables
  integer i,j
  real    sum
  integer cnt
  
  !   Calculate the minimum and maximum
  min=tmp(0,0)
  max=tmp(0,0)
  do i=0,nx
     do j=0,ny
        if (tmp(i,j).lt.min) min=tmp(i,j)
        if (tmp(i,j).gt.max) max=tmp(i,j)
     enddo
  enddo
  
  !   Calculate the mean
  sum=0.
  cnt=0
  do i=0,nx
     do j=0,ny
        cnt=cnt+1
        sum=sum+tmp(i,j)
     enddo
  enddo
  if (cnt.ge.1) then
     mean=sum/real(cnt)
  else
     mean=0.
  endif

!   Calculate rmsq
  rmsq=0.
  cnt=0
  do i=0,nx
     do j=0,ny
        cnt=cnt+1
        rmsq=rmsq+(tmp(i,j)-mean)**2
     enddo
  enddo
  if (cnt.ge.1) then
     rmsq=1./real(cnt)*sqrt(rmsq)
  else
     rmsq=0.
  endif
  
end subroutine calc_error
      


      
      
      
         
