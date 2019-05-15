program inv_cart

!     ********************************************************************************
!     * CALCULATES FOR A PV AND THETA DISTRIBUTION OTHER PROGNOSTIC FIELDS BY        *
!     * MEANS OF A PV INVERSION                                                      *
!     *                                                                              *
!     * Rene Fehlmann 1994 / Code re-organization: Michael Sprenger, 2006            *
!     ********************************************************************************

!     --------------------------------------------------------------------------------
!     Declaration of variables, parameters, externals and common blocks            
!     --------------------------------------------------------------------------------
  use kind_parameters,ONLY:&
       sp
  use netcdflibrary
  use cdfreadwrite
  implicit none

  !     Grid parameters
  integer        nx,ny,nz
  real(kind=sp)           xmin,ymin,zmin,xmax,ymax,zmax
  real(kind=sp)           dx,dy,dz
  real(kind=sp)           deltax,deltay,deltaz
  real(kind=sp)           mdv
  real(kind=sp),          allocatable,dimension (:,:) :: coriol

!     Reference state
  real(kind=sp),   allocatable,dimension (:) :: nsqref
  real(kind=sp),   allocatable,dimension (:) :: thetaref
  real(kind=sp),   allocatable,dimension (:) :: rhoref
  real(kind=sp),   allocatable,dimension (:) :: pressref
  real(kind=sp),   allocatable,dimension (:) :: zref


!     Boundary conditions
  real(kind=sp),   allocatable,dimension (:,:) :: thetatop
  real(kind=sp),   allocatable,dimension (:,:) :: thetabot
  real(kind=sp),   allocatable,dimension (:,:) :: ua
  real(kind=sp),   allocatable,dimension (:,:) :: ub
  real(kind=sp),   allocatable,dimension (:,:) :: va
  real(kind=sp),   allocatable,dimension (:,:) :: vb

!     Potentiual vorticity and stream function
  real(kind=sp),   allocatable,dimension (:,:,:) :: psi
  real(kind=sp),   allocatable,dimension (:,:,:) :: pv
  
  !     Auxiliary arrays for numerics
  real(kind=sp),   allocatable,dimension (:,:,:) :: a
  real(kind=sp),   allocatable,dimension (:,:,:) :: b
  real(kind=sp),   allocatable,dimension (:,:,:) :: c
  
  !     Input parameters
  character*160         pvsrcfile
  character*160         referfile
  
!     Auxiliary variables
  integer             i,j,k
  integer             stat

!     --------------------------------------------------------------------------------
!     Preparations
!     --------------------------------------------------------------------------------

  print*,'********************************************************'
  print*,'* INV_CART                                             *'
  print*,'********************************************************'

!     Read parameter file

  open(10,file='fort.10')
!  read(10,*,end =100) pvsrcfile
  read(10,*)pvsrcfile
  read(10,*) referfile
!100 close(10)
  print*
  print*,trim(pvsrcfile)
  print *,referfile
  !     Get lat/lon gid parameters from input file
  call read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,pvsrcfile)
  print*
  print*,'Read_Dim: nx,ny,nz         = ',nx,ny,nz
  print*,'          dx,dy,dz         = ',dx,dy,dz
  print*,'          xmin,ymin,zmin   = ',xmin,ymin,zmin
  print*,'          mdv              = ',mdv
  print*

      !     Count from 0, not from 1
  nx=nx-1
  ny=ny-1
  nz=nz-1
  
      !     Allocate memory for boundary conditions
  allocate(thetatop(0:nx,0:ny),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array thetatop ***'
  allocate(thetabot(0:nx,0:ny),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array thetabot ***'
  allocate(ua(0:nx,0:nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array ua ***'
  allocate(ub(0:nx,0:nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array ub ***'
  allocate(va(0:ny,0:nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array va ***'
  allocate(vb(0:ny,0:nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array vb ***'

  !     Allocate memory for 3d PV and stream function
  allocate(psi(0:nx,0:ny,0:nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array psi ***'
  allocate(pv(0:nx,0:ny,0:nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array pv ***'
  
  !     Alllocate memory for matrix elements for inversion operator
  allocate(a(0:nx,0:ny,0:nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array a ***'
  allocate(b(0:nx,0:ny,0:nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array b ***'
  allocate(c(0:nx,0:ny,0:nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array c ***'
  
  !     Allocate memory for reference profile
  allocate(nsqref(0:2*nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array nsqref ***'
  allocate(thetaref(0:2*nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array thetaref ***'
  allocate(rhoref(0:2*nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array rhoref ***'
  allocate(pressref(0:2*nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array pressref ***'
  allocate(zref(0:2*nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array zref ***'
  
      !     Allocate memory for Coriolis parameter
  allocate(coriol(0:nx,0:ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating coriol'

!     --------------------------------------------------------------------------------
!     Input
!     --------------------------------------------------------------------------------

!     Read reference profile and grid parameters
  call read_ref (nsqref,rhoref,thetaref,pressref,zref,nx,ny,nz,deltax,deltay,deltaz,coriol,referfile)
  deltax=1000.*deltax
  deltay=1000.*deltay
  print *,pvsrcfile
  !     Read input fields from netcdf 
  call read_inp (pv,thetabot,thetatop,ua,ub,va,vb,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,pvsrcfile)

!     --------------------------------------------------------------------------------
!     Perform the inversion
!     --------------------------------------------------------------------------------

!     Init matrix elements for inversion
  call matrixel(a,b,c,coriol,nx,ny,nz,nsqref,rhoref,deltax,deltay,deltaz)
  
  !     Inversion
  call sor(psi,nsqref,thetatop,thetabot,thetaref,rhoref,coriol,pv,ua,ub,va,vb,a,b,c,nx,ny,nz,deltax,deltay,deltaz)

!     --------------------------------------------------------------------------------
!     Output
!     --------------------------------------------------------------------------------

!     Write output to netcdf
  print*
  call write_out (psi,thetabot,thetatop,ua,ub,va,vb,nx,ny,nz,deltax,deltay,deltaz,coriol,thetaref,rhoref,pressref,pvsrcfile)
      
  
end program inv_cart


!     ********************************************************************************
!     * NETCDF AND ASCII INPUT AND OUTPUT                                            *
!     ********************************************************************************





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
  use kind_parameters,ONLY : &
       sp
  
  implicit none

!   Declaration of subroutine parameters
  character*80   pvsrcfile
  integer        nx,ny,nz
  real(kind=sp)  dx,dy,dz
  real(kind=sp)  xmin,ymin,zmin,xmax,ymax,zmax
  real(kind=sp)  mdv

!   Numerical epsilon and other physical/geoemtrical parameters
  real(kind=sp)           eps
  parameter      (eps=0.01)

!   Auxiliary variables
  integer        cdfid,cstid
  integer        ierr
  character*80   vnam(100),varname
  integer        nvars
  integer        isok
  integer        vardim(4)
  real(kind=sp)  misdat
!  real(kind=sp)(kind=sp)  varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  real(kind=sp)  aklev(1000),bklev(1000),aklay(1000),bklay(1000)
  real(kind=sp)  dh
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
  varname='QGPV'
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


!   -------------------------------------------------------------------------------
!   Read the input netcdf file
!   --------------------------------------------------------------------------------

SUBROUTINE read_inp (pv,thetabot,thetatop,ua,ub,va,vb,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,pvsrcfile)
  
  !   Read all needed field from netcdf file <pvsrcfile>. The input fields are:
  !      pv                : quasigeostrophic potential vorticity
  !      thetabot,thetatop : potential temperature at lower and upper boundary
  !      ua,ub             : Zonal wind at western and eastern boundary
  !      va,vb             : Meridional wind at southern and northern boundary
  !   The grid is specified by <nx,ny,nz,dx,dy,dz,xmin,ymin,zmin>. A check is performed
  !   whether the input files are consitent with this grid. The input netcdf file must
  !   contain the variables <QGPV,THETA,U,V>. If the netcdf file also contains the fields
  !   <DQGPV,DTHETA,DU,DV>, these increments are added to <QGPV,THETA,U,V>.
  
  use netcdflibrary
  use kind_parameters,ONLY : &
       sp
  implicit   none

  !   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp)        pv(0:nx,0:ny,0:nz)
  real(kind=sp)        thetatop(0:nx,0:ny)
  real(kind=sp)        thetabot(0:nx,0:ny)
  real(kind=sp)        ua(0:nx,0:nz)
  real(kind=sp)        ub(0:nx,0:nz)
  real(kind=sp)        va(0:ny,0:nz)
  real(kind=sp)        vb(0:ny,0:nz)
  character*(80)       pvsrcfile
  real(kind=sp)        dx,dy,dz,xmin,ymin,zmin

  !   Numerical and physical parameters
  real(kind=sp)                 eps
  parameter            (eps=0.01)
  
  !   Auxiliary variables
  integer              cdfid
  integer              vardim(4)
  real(kind=sp)        misdat
!  real(kind=sp)(kind=sp)        varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  integer              ndimin,outid,i,j,k
  real(kind=sp)                 max_th

  integer              ntimes
  real(kind=sp),allocatable,dimension(:)::time
  real(kind=sp),allocatable,dimension(:,:,:)::tmp
  integer              nvars
  character*80         vnam(100),varname
  integer              isok
  allocate(tmp(0:nx,0:ny,0:nz))
  allocate(time(1))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))

  !   Open the input netcdf file
  call cdfopn(pvsrcfile,cdfid)
  
  !   Check whether needed variables are on file
  call getvars(cdfid,nvars,vnam)
  isok=0
  varname='TH'
  call check_varok(isok,varname,vnam,nvars)
  varname='U'
  call check_varok(isok,varname,vnam,nvars)
  varname='V'
  call check_varok(isok,varname,vnam,nvars)
  varname='QGPV'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.ne.4) goto 998

!   Get the grid parameters 
  varname='QGPV'
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)    
  time(1)=0.
  call gettimes(cdfid,time,ntimes)  
  
  !   Check whether grid parameters are consitent
  if ( (vardim(1).ne.nx+1).or.(vardim(2).ne.ny+1).or.(vardim(3).ne.nz+1).or.(abs(varmin(1)-xmin).gt.eps).or. &
       (abs(varmin(2)-ymin).gt.eps).or.(abs(varmin(3)-zmin).gt.eps).or. (abs((varmax(1)-varmin(1))/real(nx)-dx).gt.eps).or.&
       (abs((varmax(2)-varmin(2))/real(ny)-dy).gt.eps).or.(abs((varmax(3)-varmin(3))/real(nz)-dz).gt.eps) ) then
     print*,'Input grid inconsitency...'
     print*,'xmin : ',xmin,varmin(1)
     print*,'ymin : ',ymin,varmin(2)
     print*,'zmin : ',zmin,varmin(3)
     print*,'dx   : ',dx,(varmax(1)-varmin(1))/real(nx)
     print*,'dy   : ',dy,(varmax(2)-varmin(2))/real(ny)
     print*,'dz   : ',dz,(varmax(3)-varmin(3))/real(nz)
     print*,'nx   : ',nx
     print*,'ny   : ',ny
     print*,'nz   : ',nz
     goto 998
  endif

!   THETA: Load upper and lower boundary values
  varname='TH'
  call getdat(cdfid,varname,time,0,tmp)  
  print*,'R TH       ',trim(pvsrcfile)
  do i=0,nx
     do j=0,ny
        if (abs(tmp(i,j,0)-misdat).lt.eps ) then
           thetabot(i,j)=0.
        else
           thetabot(i,j)=tmp(i,j,0)     
        endif
        if ( abs(tmp(i,j,nz)-misdat).lt.eps ) then
           thetatop(i,j)=0.
        else
           thetatop(i,j)=tmp(i,j,nz)
        endif
     enddo
  enddo

!   U: Load zonal velocity at southern and northern boundary
  varname='U'
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)
  call getdat(cdfid,varname,time,0,tmp)
  print*,'R U        ',trim(pvsrcfile)
  do i=0,nx
     do k=0,nz
        if ( abs(tmp(i,0,k)-misdat).lt.eps ) then
           ua(i,k)=0.
        else
           ua(i,k)=tmp(i,0,k)
        endif
        if ( abs(tmp(i,ny,k)-misdat).lt.eps ) then
           ub(i,k)=0.
        else
           ub(i,k)=tmp(i,ny,k)
        endif
     enddo
  enddo
      
!   Load meridional velocity at western and eastern boundary
  varname='V'
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)
  call getdat(cdfid,varname,time,0,tmp)
  print*,'R V        ',trim(pvsrcfile)
  do j=0,ny
     do k=0,nz
        if ( abs(tmp(0,j,k)-misdat).lt.eps ) then
           va(j,k)=0.
        else
           va(j,k)=tmp(0,j,k)
        endif
        if ( abs(tmp(nx,j,k)-misdat).lt.eps ) then
           vb(j,k)=0.
        else
           vb(j,k)=tmp(nx,j,k)
        endif
     enddo
  enddo
  
!   Load qgPV
  varname='QGPV'
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)
  call getdat(cdfid,varname,time,0,tmp)
  print*,'R QGPV     ',trim(pvsrcfile)
  do i=0,nx
     do j=0,ny
        do k=0,nz
           if ( abs(tmp(i,j,k)-misdat).lt.eps ) then
              pv(i,j,k)=0.
           else
              pv(i,j,k)=tmp(i,j,k)
           endif
        enddo
     enddo
  enddo

!   Close input netcdf file
  call clscdf(cdfid)
  
  return
  
  !   Exception handling
998 print*,'Problem with input netcdf file... Stop'
  stop
  
end subroutine read_inp

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


!   ********************************************************************************
!   * INVERSION ROUTINES                                                           *
!   ********************************************************************************

!   --------------------------------------------------------------------------------
!   SOR algorithm (successive over relaxation)
!   --------------------------------------------------------------------------------

SUBROUTINE sor(psi,nsq,thetatop,thetabot,thetaref,rhoref,coriol,pv,ua,ub,va,vb,a,b,c,nx,ny,nz,dx,dy,dz)
  
  !   Solve the qgPV equation by succesive over relaxation (SOR). The subroutine
  !   parameters are:
  !
  !      psi                 : Streamfunction, i.e. result of the PV inversion
  !      nsq,rhoref,thetaref : Reference profile
  !      thetatop,thetabot   : Upper and lower boundary condition
  !      pv                  : quasigeostrophic potential vorticity (qgPV)
  !      ua,ub,va,vb         : lateral boundary condition for wind
  !      a,b,!             : Matrices for the inversion operator
  !      nx,ny,nz,dx,dy,dz   : Grid specification
  !      coriol              : Coriolis parameter
  use kind_parameters,ONLY : &
       sp
  implicit none
  
  !   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp)                 dx,dy,dz
  real(kind=sp)                 psi     (0:nx,0:ny,0:nz)
  real(kind=sp)                 nsq     (0:2*nz)
  real(kind=sp)                 thetatop(0:nx,0:ny)
  real(kind=sp)                 thetabot(0:nx,0:ny)
  real(kind=sp)                 thetaref(0:2*nz)
  real(kind=sp)                 rhoref  (0:2*nz)
  real(kind=sp)                 pv      (0:nx,0:ny,0:nz)
  real(kind=sp)                 ua      (0:nx,0:nz)
  real(kind=sp)                 ub      (0:nx,0:nz)
  real(kind=sp)                 va      (0:ny,0:nz)
  real(kind=sp)                 vb      (0:ny,0:nz)
  real(kind=sp)                 a       (0:nx,0:ny,0:nz)
  real(kind=sp)                 b       (0:nx,0:ny,0:nz)
  real(kind=sp)                 c       (0:nx,0:ny,0:nz)
  real(kind=sp)                 coriol  (0:nx,0:ny)
      
  !   Numerical and physical parameters
  real(kind=sp)                 maxspec
  parameter            (maxspec=2.0)
  integer              nofiter
  parameter            (nofiter=500)
  real(kind=sp)                 omega
  parameter            (omega=1.81)
  
  !   Auxiliary variables
  integer              counter
  integer              i,j,k
  real(kind=sp)                 deltasq,psigauge
  real(kind=sp)                 specx,specy,specz
  real(kind=sp)                 helpx,helpy,helpz
      
  !   Init the output array
  do i=0,nx
     do j=0,ny
        do k=0,nz
           psi(i,j,k)=0.
        enddo
     enddo
  enddo
  
      !   Calculate the spectrum of the matrix 
  i=nx/2
  
  specx=4.*a(i,0,0)/(2.*a(i,0,0)+b(i,0,0)+b(i,1,0)+c(i,0,0)+c(i,0,1))
  specy=2.*(b(i,0,0)+b(i,1,0))/(2.*a(i,0,0)+b(i,0,0)+b(i,1,0)+c(i,0,0)+c(i,0,1))
  specz=2.*(c(i,0,0)+c(i,0,1))/(2.*a(i,0,0)+b(i,0,0)+b(i,1,0)+c(i,0,0)+c(i,0,1))
  
  do k=1,nz-2
     do j=1,ny-2
        
        helpx=4.*a(i,j,k)/ (2.*a(i,j,k)+b(i,j,k)+b(i,j-1,k)+c(i,j,k)+c(i,j,k-1))
        
        if (helpx.gt.specx) specx=helpx
        
        helpy=2.*(b(i,j,k)+b(i,j+1,k))/(2.*a(i,j,k)+b(i,j,k)+b(i,j-1,k)+c(i,j,k)+c(i,j,k-1))
        if (helpy.gt.specy) specy=helpy
        
        helpz=2.*(c(i,j,k)+c(i,j,k+1))/(2.*a(i,j,k)+b(i,j,k)+b(i,j-1,k)+c(i,j,k)+c(i,j,k-1))
        if (helpz.gt.specz) specz=helpz
        
     enddo
  enddo
      
  !   Check whether the dimensions of the grid are sufficient
  print *
  print *, 'Spectrum of the matrix in each direction '
  print *, 'Spectrum = ', specx, specy, specz
  print *
  if ((maxspec*specx.lt.specy).or.(maxspec*specx.lt.specz)) then
     print*,' Nx too small... Stop'
     stop
  endif
  if ((maxspec*specy.lt.specx).or.(maxspec*specy.lt.specz)) then
     print*,'Ny too small... Stop'
     stop
  endif
  if ((maxspec*specz.lt.specx).or.(maxspec*specz.lt.specy)) then
     print*,'Nz too small... Stop'
     stop
  endif
  
!   Calculate error: control variable for the iteration
  psigauge=0.
  deltasq=0.
  do k=1,nz-1
     do i=1,nx-1
        do j=1,ny-1
           deltasq=deltasq+(-pv(i,j,k)+(a(i,j,k)* &
                (psi(i+1,j,k)+psi(i-1,j,k)- &
                2.*psi(i,j,k)) + &
                b(i,j,k)*(psi(i,j+1,k)-psi(i,j,k))- &
                b(i,j-1,k)*(psi(i,j,k)-psi(i,j-1,k))+ &
                c(i,j,k)*(psi(i,j,k+1)-psi(i,j,k))- &
                c(i,j,k-1)*(psi(i,j,k)-psi(i,j,k-1)) &
                )/(dx*dy*dz*rhoref(2*k)))**2.
               
        enddo
     enddo
  enddo
  print 102, 'psigauge', psigauge, 'deltasq',deltasq/(real(nx)*real(ny)*real(nz))
  
  !   Iterations
  do counter=1,nofiter
     !      Perform one iteration step
     call psiappsor(omega,pv,psi,nsq,rhoref,thetatop,thetabot,thetaref,coriol,ua,ub,va,vb,a,b,c,nx,ny,nz,dx,dy,dz)
     
     
     !      Adjustment
     if (mod(counter,100).eq.0) then
        psigauge=0.
        do i=0,nx
           do j=0,ny
              if (psi(i,j,0).lt.psigauge) then 
                 psigauge=psi(i,j,0)
              endif
           enddo
        enddo
        do k=0,nz
           do i=0,nx
              do j=0,ny
                 psi(i,j,k)=psi(i,j,k)-psigauge
              enddo
           enddo
        enddo
     endif
     
     !      Calculate error: control variable for the iteration        
     if (mod(counter,nofiter/10).eq.0) then
        deltasq=0.
        do k=1,nz-1
           do i=1,nx-1
              do j=1,ny-1
                 deltasq=deltasq+(-pv(i,j,k)+( &
                      a(i,j,k)*(psi(i+1,j,k)+psi(i-1,j,k)- &
                      2.*psi(i,j,k)) + &
                      b(i,j,k)*(psi(i,j+1,k)-psi(i,j,k))- &
                      b(i,j-1,k)*(psi(i,j,k)-psi(i,j-1,k))+ &
                      c(i,j,k)*(psi(i,j,k+1)-psi(i,j,k))- &
                      c(i,j,k-1)*(psi(i,j,k)-psi(i,j,k-1)) &
                      )/(dx*dy*dz*rhoref(2*k)))**2. 
              enddo
           enddo
        enddo
        print 102, 'psigauge', psigauge, 'deltasq',deltasq/(real(nx)*real(ny)*real(nz))
     endif
     
  enddo
     
  return
         
         !   Format specifications
102  format (a11, ' = ',e10.3,a11, ' = ',e10.3)
     
end subroutine sor

!   --------------------------------------------------------------------------------
!   SOR algorithm (successive over relaxation)
!   --------------------------------------------------------------------------------
       
subroutine psiappsor(omega,pv,psi,nsq,rhoref,thetatop,thetabot,thetaref,coriol,ua,ub,va,vb,a,b,c,nx,ny,nz,dx,dy,dz)
  
  !   Perform one relaxation step
  !
  !      psi                 : Streamfunction, i.e. result of the PV inversion
  !      nsq,rhoref,thetaref : Reference profile
  !      thetatop,thetabot   : Upper and lower boundary condition
  !      pv                  : quasigeostrophic potential vorticity (qgPV)
  !      ua,ub,va,vb         : lateral boundary condition for wind
  !      a,b,!             : Matrices for the inversion operator
  !      nx,ny,nz,dx,dy,dz   : Grid specification
  !      nofiter             : Number of iterations
  !      omega               : Relaxation parameter
  !      coriol              : Coriolis parameter
  use ieee_arithmetic
  use ieee_features
  use kind_parameters,ONLY : &
       sp
  
  implicit none
  logical :: underflow_support, gradual, underflow

  !   Declaration of subroutine parameters
  integer             nx,ny,nz
  real(kind=sp)       pv(0:nx,0:ny,0:nz)
  real(kind=sp)       psi(0:nx,0:ny,0:nz)
  real(kind=sp)       nsq(0:2*nz)
  real(kind=sp)       rhoref(0:2*nz)
  real(kind=sp)       thetatop(0:nx,0:ny)
  real(kind=sp)       thetabot(0:nx,0:ny)
  real(kind=sp)       thetaref(0:2*nz)
  real(kind=sp)       ua(0:nx,0:nz)
  real(kind=sp)       ub(0:nx,0:nz)
  real(kind=sp)       va(0:ny,0:nz)
  real(kind=sp)       vb(0:ny,0:nz)
  real(kind=sp)       a(0:nx,0:ny,0:nz)
  real(kind=sp)       b(0:nx,0:ny,0:nz)
  real(kind=sp)       c(0:nx,0:ny,0:nz)
  real(kind=sp)       coriol(0:nx,0:ny)
  
  real(kind=sp)        dx,dy,dz
  real(kind=sp)        omega
  
  !   Numerical and physical parameters
  real(kind=sp)                 g
  real(kind=sp)        expr0,expr1,expr2,expr3
  real(kind=sp)                 expr4,expr5,expr6,expr7
  parameter            (g=9.81)
  
  !   Auxiliary variables
  integer              i,j,k
  real(kind=sp)        dxy,dxz,dyz,dxyz

  
  !   Set the area and volume infinitesimals for integration
  dxy=dx*dy
  dxz=dx*dz
  dyz=dy*dz
  dxyz=dx*dy*dz

  call ieee_set_underflow_mode(.false.)
  call ieee_get_underflow_mode(gradual)
  if (.not.gradual) then 
!     print *,'Able to set abrupt underflow mode'
  else
 !    stop 'error setting underflow mode'
  end if
  !   Inner

  do k=1,nz-1
     do i=1,nx-1
        do j=1,ny-1

           psi(i,j,k)=omega*(-dxyz*rhoref(2*k)*pv(i,j,k)+a(i,j,k)*&
                (psi(i+1,j,k)+psi(i-1,j,k))&
                +b(i,j,k)*psi(i,j+1,k)+b(i,j-1,k)&
                *psi(i,j-1,k)+c(i,j,k)*&
                psi(i,j,k+1)+c(i,j,k-1)&
                *psi(i,j,k-1))/&
                (2.*a(i,j,k)+b(i,j,k)+b(i,j-1,k)&
                +c(i,j,k-1)+&
                c(i,j,k))&
                +(1.-omega)*psi(i,j,k)
           call ieee_get_flag(ieee_underflow,underflow)                    
        enddo
     enddo
  enddo
!   ZY plane
  do k=1,nz-1
     do j=1,ny-1
        psi(0,j,k)=omega*(-dyz*&
             rhoref(2*k)*(dx*pv(0,j,k)+va(j,k))+&
             a(0,j,k)*psi(1,j,k)+&
             b(0,j,k)*psi(0,j+1,k)+b(0,j-1,k)*psi(0,j-1,k)+&
             c(0,j,k)*psi(0,j,k+1)+c(0,j,k-1)*psi(0,j,k-1))/&
             (a(0,j,k)+b(0,j,k)+b(0,j-1,k)+c(0,j,k-1)+c(0,j,k))&
             +(1.-omega)*psi(0,j,k)
        !
        expr4 = omega*(-dyz*&
             rhoref(2*k)*(dx*pv(nx,j,k)-vb(j,k))+&
             a(nx,j,k)*psi(nx-1,j,k)+&
             b(nx,j,k)*psi(nx,j+1,k)+b(nx,j-1,k)*psi(nx,j-1,k)+&
             c(nx,j,k)*psi(nx,j,k+1)+c(nx,j,k-1)*psi(nx,j,k-1))

        
        expr5 = (a(nx,j,k)+b(nx,j-1,k)+b(nx,j,k)+c(nx,j,k-1)+c(nx,j,k))
        expr6=   (1.-omega)*psi(nx,j,k)
        expr7 = expr4/expr5 + expr6
        psi(nx,j,k)= expr7
     enddo
  enddo

!   ZX plane
  do k=1,nz-1
     do i=1,nx-1
        psi(i,0,k)=omega*(-dxz*&
             rhoref(2*k)*(dy*pv(i,0,k)-ua(i,k))+&
                a(i,0,k)*(psi(i+1,0,k)+psi(i-1,0,k))+&
                b(i,0,k)*psi(i,1,k)+&
                c(i,0,k)*psi(i,0,k+1)+c(i,0,k-1)*psi(i,0,k-1))/&
                (2.*a(i,0,k)+b(i,0,k)+c(i,0,k-1)+c(i,0,k))&
                +(1.-omega)*psi(i,0,k)
        !   
        psi(i,ny,k)=omega*(-dxz*&
             rhoref(2*k)*(dy*pv(i,ny,k)+ub(i,k))+&
             a(i,ny-1,k)*(psi(i+1,ny,k)+psi(i-1,ny,k))+&
             b(i,ny-1,k)*psi(i,ny-1,k)+&
             c(i,ny-1,k)*psi(i,ny,k+1)+c(i,ny-1,k-1)*&
             psi(i,ny,k-1))/(2.*a(i,ny-1,k)+b(i,ny-1,k)+&
             c(i,ny-1,k-1)+c(i,ny-1,k))&
             +(1.-omega)*psi(i,ny,k)
     enddo
  enddo
!   XY plane
  do i=1,nx-1
     do j=1,ny-1
        psi(i,j,0)=omega*(-dxy*rhoref(0)*(&
             dz*pv(i,j,0)+g*coriol(i,j)*thetabot(i,j)/&
             (nsq(0)*thetaref(0)))+&
             a(i,j,0)*(psi(i+1,j,0)+psi(i-1,j,0))+&
             b(i,j,0)*psi(i,j+1,0)+b(i,j-1,0)&
             *psi(i,j-1,0)+&
             c(i,j,0)*psi(i,j,1))/&
             (2.*a(i,j,0)+b(i,j-1,0)+b(i,j,0)+&
             c(i,j,0))&
             +(1.-omega)*psi(i,j,0)
        !   
        psi(i,j,nz)=omega*(-dxy*rhoref(2*nz)*(&
             dz*pv(i,j,nz)-g*coriol(i,j)*thetatop(i,j)/&
             (nsq(2*nz)*thetaref(2*nz)))+&
             a(i,j,nz)*(psi(i+1,j,nz)+psi(i-1,j,nz))+&
             b(i,j,nz)*psi(i,j+1,nz)+b(i,j-1,nz)*psi(i,j-1,nz)+&
             c(i,j,nz-1)*psi(i,j,nz-1))/&
             (2.*a(i,j,nz)+b(i,j-1,nz)+b(i,j,nz)+c(i,j,nz-1))&
             +(1.-omega)*psi(i,j,nz)
     enddo
  enddo

!   Y edges
  do j=1,ny-1
     psi(0,j,0)=omega*(-dy*rhoref(0)*(dxz*pv(0,j,0)+&
          dz*va(j,0)+dx*g*coriol(0,j)*thetabot(0,j)/&
          (nsq(0)*thetaref(0)))+&
          a(0,j,0)*psi(1,j,0)+&
          b(0,j,0)*psi(0,j+1,0)+b(0,j-1,0)*psi(0,j-1,0)+&
          c(0,j,0)*psi(0,j,1))/&
          (a(0,j,0)+b(0,j-1,0)+b(0,j,0)+c(0,j,0))&
          +(1.-omega)*psi(0,j,0)
     !   
     psi(nx,j,0)=omega*(-dy*rhoref(0)*(dxz*pv(nx,j,0)-&
          dz*vb(j,0)+dx*g*coriol(nx,j)*thetabot(nx,j)/&
          (nsq(0)*thetaref(0)))+&
          a(nx,j,0)*psi(nx-1,j,0)+&
          b(nx,j,0)*psi(nx,j+1,0)+b(nx,j-1,0)*psi(nx,j-1,0)+&
          c(nx,j,0)*psi(nx,j,1))/&
          (a(nx,j,0)+b(nx,j-1,0)+b(nx,j,0)+c(nx,j,0))&
          +(1.-omega)*psi(nx,j,0)
!
     psi(0,j,nz)=omega*(-dy*rhoref(2*nz)*(dxz*pv(0,j,nz)+&
          dz*va(j,nz)-dx*g*coriol(0,j)*thetatop(0,j)/&
          (nsq(2*nz)*thetaref(2*nz)))+&
          a(0,j,nz)*psi(1,j,nz)+&
          b(0,j,nz)*psi(0,j+1,nz)+b(0,j-1,nz)*psi(0,j-1,nz)+&
          c(0,j,nz-1)*psi(0,j,nz-1))/&
          (a(0,j,nz)+b(0,j-1,nz)+b(0,j,nz)+c(0,j,nz-1))&
          +(1.-omega)*psi(0,j,nz)
     !   
     psi(nx,j,nz)=omega*(-dy*rhoref(2*nz)*(dxz*pv(nx,j,nz)-&
          dz*vb(j,nz)-dx*g*coriol(nx,j)*thetatop(nx,j)/&
          (nsq(2*nz)*thetaref(2*nz)))+&
          a(nx,j,nz)*psi(nx-1,j,nz)+&
          b(nx,j,nz)*psi(nx,j+1,nz)+b(nx,j-1,nz)*psi(nx,j-1,nz)+&
          c(nx,j,nz-1)*psi(nx,j,nz-1))/&
          (a(nx,j,nz)+b(nx,j-1,nz)+b(nx,j,nz)+c(nx,j,nz-1))&
          +(1.-omega)*psi(nx,j,nz)
  enddo
  
  !   X edges
  do i=1,nx-1
     psi(i,0,0)=omega*(-dx*rhoref(0)*(dyz*pv(i,0,0)-&
          dz*ua(i,0)+dy*g*coriol(i,0)*thetabot(i,0)/&
          (nsq(0)*thetaref(0)))+&
          a(i,0,0)*(psi(i+1,0,0)+psi(i-1,0,0))+&
          b(i,0,0)*psi(i,1,0)+&
          c(i,0,0)*psi(i,0,1))/&
          (2.*a(i,0,0)+b(i,0,0)+c(i,0,0))&
          +(1.-omega)*psi(i,0,0)
     !   
     psi(i,ny,0)=omega*(-dx*rhoref(0)*(dyz*pv(i,ny,0)+&
          dz*ub(i,0)+dy*g*coriol(i,ny)*thetabot(i,ny)/&
          (nsq(0)*thetaref(0)))+&
          a(i,ny,0)*(psi(i+1,ny,0)+psi(i-1,ny,0))+&
          b(i,ny-1,0)*psi(i,ny-1,0)+&
          c(i,ny,0)*psi(i,ny,1))/&
          (2.*a(i,ny,0)+b(i,ny-1,0)+c(i,ny,0))&
          +(1.-omega)*psi(i,ny,0)
     !   
     psi(i,0,nz)=omega*(-dx*rhoref(2*nz)*(dyz*pv(i,0,nz)-&
          dz*ua(i,nz)-dy*g*coriol(i,0)*thetatop(i,0)/&
          (nsq(2*nz)*thetaref(2*nz)))+&
          a(i,0,nz)*(psi(i+1,0,nz)+psi(i-1,0,nz))+&
          b(i,0,nz)*psi(i,1,nz)+&
          c(i,0,nz-1)*psi(i,0,nz-1))/&
          (2.*a(i,0,nz)+b(i,0,nz)+c(i,0,nz-1))&
          +(1.-omega)*psi(i,0,nz)
     !   
     psi(i,ny,nz)=omega*(-dx*rhoref(2*nz)*(dyz*pv(i,ny,nz)+&
          dz*ub(i,nz)-dy*g*coriol(i,ny)*thetatop(i,ny)/&
          (nsq(2*nz)*thetaref(2*nz)))+&
          a(i,ny,nz)*(psi(i+1,ny,nz)+psi(i-1,ny,nz))+&
          b(i,ny-1,nz)*psi(i,ny-1,nz)+&
          c(i,ny,nz-1)*psi(i,ny,nz-1))/&
          (2.*a(i,ny,nz)+b(i,ny-1,nz)+c(i,ny,nz-1))&
          +(1.-omega)*psi(i,ny,nz)
  enddo
  !   Z edges
  do k=1,nz-1
     psi(0,0,k)=omega*(-dz*rhoref(2*k)*(dxy*pv(0,0,k)+&
          dy*va(0,k)-dx*ua(0,k))+&
          a(0,0,k)*psi(1,0,k)+&
          b(0,0,k)*psi(0,1,k)+&
          c(0,0,k)*psi(0,0,k+1)+c(0,0,k-1)*psi(0,0,k-1))/&
          (a(0,0,k)+b(0,0,k)+c(0,0,k-1)+c(0,0,k))&
          +(1.-omega)*psi(0,0,k)
     !
     psi(nx,0,k)=omega*(-dz*rhoref(2*k)*(dxy*pv(nx,0,k)-&
          dy*vb(0,k)-dx*ua(nx,k))+&
          a(nx,0,k)*psi(nx-1,0,k)+&
          b(nx,0,k)*psi(nx,1,k)+&
          c(nx,0,k)*psi(nx,0,k+1)+c(nx,0,k-1)*psi(nx,0,k-1))/&
          (a(nx,0,k)+b(nx,0,k)+c(nx,0,k-1)+c(nx,0,k))&
          +(1.-omega)*psi(nx,0,k)
     !
     psi(0,ny,k)=omega*(-dz*rhoref(2*k)*(dxy*pv(0,ny,k)+&
          dy*va(ny,k)+dx*ub(0,k))+&
          a(0,ny,k)*psi(1,ny,k)+&
          b(0,ny-1,k)*psi(0,ny-1,k)+&
          c(0,ny,k)*psi(0,ny,k+1)+c(0,ny,k-1)*psi(0,ny,k-1))/&
          (a(0,ny,k)+b(0,ny-1,k)+c(0,ny,k-1)+c(0,ny,k))&
          +(1.-omega)*psi(0,ny,k)
     !
     psi(nx,ny,k)=omega*(-dz*rhoref(2*k)*(dxy*pv(nx,ny,k)-&
          dy*vb(ny,k)+dx*ub(nx,k))+&
          a(nx,ny,k)*psi(nx-1,ny,k)+&
          b(nx,ny-1,k)*psi(nx,ny-1,k)+&
          c(nx,ny,k)*psi(nx,ny,k+1)+c(nx,ny,k-1)*psi(nx,ny,k-1))/&
          (a(nx,ny,k)+b(nx,ny-1,k)+c(nx,ny,k-1)+c(nx,ny,k))&
          +(1.-omega)*psi(nx,ny,k)
  enddo
  !   Points
  psi(0,0,0)=omega*(-rhoref(0)*(dxyz*pv(0,0,0)+dyz*va(0,0)-&
       dxz*ua(0,0)+dxy*g*coriol(0,0)*thetabot(0,0)/&
       (nsq(0)*thetaref(0)))+&
       a(0,0,0)*psi(1,0,0)+&
       b(0,0,0)*psi(0,1,0)+&
       c(0,0,0)*psi(0,0,1))/&
       (a(0,0,0)+b(0,0,0)+c(0,0,0))+&
       (1.-omega)*psi(0,0,0)
  !
  psi(nx,0,0)=omega*(-rhoref(0)*(dxyz*pv(nx,0,0)-dyz*vb(0,0)-&
       dxz*ua(nx,0)+dxy*g*coriol(nx,0)*thetabot(nx,0)/&
             (nsq(0)*thetaref(0)))+&
             a(nx,0,0)*psi(nx-1,0,0)+&
             b(nx,0,0)*psi(nx,1,0)+&
             c(nx,0,0)*psi(nx,0,1))/&
             (a(nx,0,0)+b(nx,0,0)+c(nx,0,0))+&
             (1.-omega)*psi(nx,0,0)
  !
  psi(0,ny,0)=omega*(-rhoref(0)*(dxyz*pv(0,ny,0)+dyz*va(ny,0)+&
       dxz*ub(0,0)+dxy*g*coriol(0,ny)*thetabot(0,ny)/&
       (nsq(0)*thetaref(0)))+&
       a(0,ny,0)*psi(1,ny,0)+&
       b(0,ny-1,0)*psi(0,ny-1,0)+&
       c(0,ny,0)*psi(0,ny,1))/&
       (a(0,ny,0)+b(0,ny-1,0)+c(0,ny,0))+&
       (1.-omega)*psi(0,ny,0)
  !
  psi(nx,ny,0)=omega*(-rhoref(0)*(dxyz*pv(nx,ny,0)-dyz*vb(ny,0)+&
       dxz*ub(nx,0)+dxy*g*coriol(nx,ny)*thetabot(nx,ny)/&
       (nsq(0)*thetaref(0)))+&
       a(nx,ny,0)*psi(nx-1,ny,0)+&
       b(nx,ny-1,0)*psi(nx,ny-1,0)+&
       c(nx,ny,0)*psi(nx,ny,1))/&
       (a(nx,ny,0)+b(nx,ny-1,0)+c(nx,ny,0))+&
       (1.-omega)*psi(nx,ny,0)
  !
  psi(0,0,nz)=omega*(-rhoref(2*nz)*(dxyz*pv(0,0,nz)+dyz*va(0,nz)-&
       dxz*ua(0,nz)-dxy*g*coriol(0,0)*thetatop(0,0)/&
       (nsq(2*nz)*thetaref(2*nz)))+&
       a(0,0,nz)*psi(1,0,nz)+&
       b(0,0,nz)*psi(0,1,nz)+&
       c(0,0,nz-1)*psi(0,0,nz-1))/&
       (a(0,0,nz)+b(0,0,nz)+c(0,0,nz-1))+&
       (1.-omega)*psi(0,0,nz)
  !
  psi(nx,0,nz)=omega*(-rhoref(2*nz)*(dxyz*pv(nx,0,nz)-dyz*vb(0,nz)&
       -dxz*ua(nx,nz)-dxy*g*coriol(nx,0)*thetatop(nx,0)/&
       (nsq(2*nz)*thetaref(2*nz)))+&
       a(nx,0,nz)*psi(nx-1,0,nz)+&
       b(nx,0,nz)*psi(nx,1,nz)+&
       c(nx,0,nz-1)*psi(nx,0,nz-1))/&
       (a(nx,0,nz)+b(nx,0,nz)+c(nx,0,nz-1))+&
       (1.-omega)*psi(nx,0,nz)
  !
  psi(0,ny,nz)=omega*(-rhoref(2*nz)*(dxyz*pv(0,ny,nz)+&
       dyz*va(ny,nz)+dxz*ub(0,nz)-&
       dxy*g*coriol(0,ny)*thetatop(0,ny)/&
       (nsq(2*nz)*thetaref(2*nz)))+&
       a(0,ny,nz)*psi(1,ny,nz)+&
       b(0,ny-1,nz)*psi(0,ny-1,nz)+&
       c(0,ny,nz-1)*psi(0,ny,nz-1))/&
       (a(0,ny,nz)+b(0,ny-1,nz)+c(0,ny,nz-1))+&
       (1.-omega)*psi(0,ny,nz)
  !
  psi(nx,ny,nz)=omega*(-rhoref(2*nz)*(dxyz*pv(nx,ny,nz)-&
       dyz*vb(ny,nz)+dxz*ub(nx,nz)-&
       dxy*g*coriol(nx,ny)*thetatop(nx,ny)/&
       (nsq(2*nz)*thetaref(2*nz)))+&
       a(nx,ny,nz)*psi(nx-1,ny,nz)+&
       b(nx,ny-1,nz)*psi(nx,ny-1,nz)+&
       c(nx,ny,nz-1)*psi(nx,ny,nz-1))/&
       (a(nx,ny,nz)+b(nx,ny-1,nz)+c(nx,ny,nz-1))+&
       (1.-omega)*psi(nx,ny,nz)
  
end subroutine psiappsor

!   --------------------------------------------------------------------------------
!   Init matrix elements for the inversion
!   --------------------------------------------------------------------------------

subroutine matrixel(a,b,c,coriol,nx,ny,nz,nsq,rhoref,dx,dy,dz)
  
  !   Define the coefficients for the inversion problem (see page 119ff in Rene's 
  !   dissertation).
  

  use kind_parameters,ONLY : &
       sp
  implicit none  
  !   Declaration of subroutine parameters
  integer   nx,nz,ny
  real(kind=sp)      a     (0:nx,0:ny,0:nz)
  real(kind=sp)      b     (0:nx,0:ny,0:nz)
  real(kind=sp)      c     (0:nx,0:ny,0:nz)
  real(kind=sp)      coriol(0:nx,0:ny)
  real(kind=sp)      nsq   (0:2*nz)
  real(kind=sp)      rhoref(0:2*nz)
  real(kind=sp)      dx,dy,dz
  
  !   Auxiliary variables
  integer   i,j,k
  
  !   Calculate coefficients
  do i=0,nx
     do j=0,ny
        do k=0,nz
           
           a(i,j,k)=dy*dz*rhoref(2*k)/(dx*coriol(i,j))
           
           if (j.lt.ny) then
              b(i,j,k)=dx*dz*rhoref(2*k)/(dy*0.5*(coriol(i,j)+coriol(i,j+1)))
           else
              b(i,j,k)=0.
           endif
           if (k.lt.nz) then
              c(i,j,k)=dx*dy*rhoref(2*k+1)*coriol(i,j)/(dz*nsq(2*k+1))
           else
              c(i,j,k)=0.
           endif
        enddo
     enddo
  enddo
  
end subroutine matrixel




