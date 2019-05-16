  PROGRAM set_boundcon

!  ************************************************************************
!  * Set boundary conditions for inversion; lower and upper boundary      *
!  * conditions for potential temperature; lateral boundary conditions    *
!  * for zonal and meridional wind; in particular, missing data values    *
!  * are removed.                                                         *
!  *                                                                      *
!  * Michael Sprenger / Summer 2006                                       *
!  ************************************************************************

!  --------------------------------------------------------------------------------
!  Declaration of variables, parameters, externals and common blocks
!  --------------------------------------------------------------------------------
  use kind_parameters,ONLY:&
       sp
  implicit none
  
  !  Input and output file
  character*80   anomafile
  character*80   referfile
  
  !  Grid parameters
  integer        nx,ny,nz
  real(kind=sp)           xmin,ymin,zmin,xmax,ymax,zmax
  real(kind=sp)           dx,dy,dz
  real(kind=sp)           mdv
  real(kind=sp)           deltax,deltay,deltaz
  real(kind=sp),          allocatable,dimension (:,:) :: coriol
  
  !  Reference state
  real(kind=sp),   allocatable,dimension (:) :: nsqref
  real(kind=sp),   allocatable,dimension (:) :: thetaref
  real(kind=sp),   allocatable,dimension (:) :: rhoref
  real(kind=sp),   allocatable,dimension (:) :: pressref
  real(kind=sp),   allocatable,dimension (:) :: zref
  
  !  Boundary conditions
  real(kind=sp),   allocatable,dimension (:,:) :: thetatop
  real(kind=sp),   allocatable,dimension (:,:) :: thetabot
  real(kind=sp),   allocatable,dimension (:,:) :: ua
  real(kind=sp),   allocatable,dimension (:,:) :: ub
  real(kind=sp),   allocatable,dimension (:,:) :: va
  real(kind=sp),   allocatable,dimension (:,:) :: vb
  
  !  3d arrays
  real(kind=sp),allocatable,dimension (:,:,:) :: th_anom,pv_anom
  real(kind=sp),allocatable,dimension (:,:,:) :: uu_anom,vv_anom
  
  !  Auxiliary variables
  integer      i,j,k,kk
  integer      stat
  character*80 varname
  integer      n1,n2
  
!  --------------------------------------------------------------------------------
!  Input
!  --------------------------------------------------------------------------------

  print*,'********************************************************'
  print*,'* CHECK_BOUNDCON                                       *'
  print*,'********************************************************'
      
  !  Read parameter file
  open(10,file='fort.10')
  read(10,*) anomafile
  read(10,*) referfile
  close(10) 
  print*
  print*,trim(anomafile)
  print*,trim(referfile)
  print*

  !  Get lat/lon gid parameters from input file
  call read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,&
       anomafile)
  print*,'Read_Dim: nx,ny,nz         = ',nx,ny,nz
  print*,'          dx,dy,dz         = ',dx,dy,dz
  print*,'          xmin,ymin,zmin   = ',xmin,ymin,zmin
  print*,'          mdv              = ',mdv
  print*
  
!  Count from 0, not from 1: consistent with <inv_cart.f>.
  nx=nx-1
  ny=ny-1
  nz=nz-1
  
  !  Allocate memory for 3d arrays 
  allocate(pv_anom (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating pv_anom'
  allocate(th_anom (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating th_anom'
  allocate(uu_anom (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating uu_anom'
  allocate(vv_anom (0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating vv_anom'
  
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
  
!  Allocate memory for Coriolis parameter
  allocate(coriol(0:nx,0:ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating f'

!  Allocate memory for boundary conditions 
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
  
!  Read reference profile and ngrid parameters
  call read_ref (nsqref,rhoref,thetaref,pressref,zref,&
       nx,ny,nz,deltax,deltay,deltaz,coriol,&
       referfile)
  deltax=1000.*deltax
  deltay=1000.*deltay
  print*,'Deltax,deltay,deltaz =',deltax,deltay,deltaz
  
!  Read data from MOD file
  varname='QGPV'
  call read_inp (pv_anom,varname,anomafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='TH'
  call read_inp (th_anom,varname,anomafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='U'
  call read_inp (uu_anom,varname,anomafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='V'
  call read_inp (vv_anom,varname,anomafile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)


!  --------------------------------------------------------------------------------
!  Consistency check for boundary conditions and adaptions
!  --------------------------------------------------------------------------------

!  Copy 3d to boundary conditions
  do i=0,nx
     do j=0,ny
        thetatop(i,j)=th_anom(i,j,nz)
        thetabot(i,j)=th_anom(i,j,0)
     enddo
  enddo
  do i=0,nx
     do k=0,nz
        ua(i,k)=uu_anom(i, 0,k)
        ub(i,k)=uu_anom(i,ny,k)
     enddo
  enddo
  do j=0,ny
     do k=0,nz
        va(j,k)=vv_anom( 0,j,k)
        vb(j,k)=vv_anom(nx,j,k)
     enddo
  enddo
  
!  Check the lower and upper boundary condition for consistency check
  print*
  call combouncon(pv_anom,nsqref,rhoref,thetatop,&
       thetabot,thetaref,coriol,ua,ub,va,vb,&
       nx,ny,nz,deltax,deltay,deltaz)
  print*
      

end PROGRAM set_boundcon


!  ********************************************************************************
!  * NETCDF INPUT AND OUTPUT                                                      *
!  ********************************************************************************

!  --------------------------------------------------------------------------------
!  Read input fields for reference profile
!  --------------------------------------------------------------------------------

SUBROUTINE read_inp (field,fieldname,pvsrcfile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)

!  Read <fieldname> from netcdf file <pvsrcfile> into <field>. The grid is specified 
!  by <nx,ny,nz,dx,dy,dz,xmin,ymin,zmin>. A check is performed whether the input 
!  files are consitent with this grid. The missing data value is set to <mdv>.
  use kind_parameters,ONLY:&
       sp
  implicit   none
  
  !  Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp)                 field(0:nx,0:ny,0:nz)
  character*80         fieldname
  character*80         pvsrcfile
  real(kind=sp)                 dx,dy,dz,xmin,ymin,zmin
  real(kind=sp)                 mdv
  
  !  Numerical and physical parameters
  real(kind=sp)                 eps
  parameter            (eps=0.01)
      
      !  Auxiliary variables
  integer              cdfid,stat,cdfid99
  integer              vardim(4)
  real(kind=sp)                 misdat
  real(kind=sp)                 varmin(4),varmax(4),stag(4)
  integer              ndimin,outid,i,j,k
  real(kind=sp)                 max_th
  real(kind=sp)                 tmp(nx,ny,nz)
  integer              ntimes
  real(kind=sp)                 time(200)       
  integer              nvars
  character*80         vnam(100),varname
  integer              isok
  
  !  Open the input netcdf file
  call cdfopn(pvsrcfile,cdfid,stat)


!  Check whether needed variables are on file
  call getvars(cdfid,nvars,vnam,stat)

  isok=0
  varname=trim(fieldname)
  call check_varok(isok,varname,vnam,nvars)


!  Get the grid parameters from theta     
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag,stat)    
  
  time(1)=0.
  call gettimes(cdfid,time,ntimes,stat)  


!  Check whether grid parameters are consistent
      if ( (vardim(1).ne.(nx+1)).or.
     >     (vardim(2).ne.(ny+1)).or.
     >     (vardim(3).ne.(nz+1)).or.
     >     (abs(varmin(1)-xmin).gt.eps).or.
     >     (abs(varmin(2)-ymin).gt.eps).or.
     >     (abs(varmin(3)-zmin).gt.eps).or.
     >     (abs((varmax(1)-varmin(1))/real(vardim(1)-1)-dx).gt.eps).or.
     >     (abs((varmax(2)-varmin(2))/real(vardim(2)-1)-dy).gt.eps).or.
     >     (abs((varmax(3)-varmin(3))/real(vardim(3)-1)-dz).gt.eps) ) 
     >then
         print*,'Input grid inconsitency...'
         print*,'  Nx      = ',vardim(1),nx+1
         print*,'  Ny      = ',vardim(2),ny+1
         print*,'  Nz      = ',vardim(3),nz+1
         print*,'  Varminx = ',varmin(1),xmin
         print*,'  Varminy = ',varmin(2),ymin
         print*,'  Varminz = ',varmin(3),zmin
         print*,'  Dx      = ',(varmax(1)-varmin(1))/real(nx-1),dx
         print*,'  Dy      = ',(varmax(2)-varmin(2))/real(ny-1),dy
         print*,'  Dz      = ',(varmax(3)-varmin(3))/real(nz-1),dz
         goto 998
      endif

!  Load variables
      call getdef(cdfid,varname,ndimin,misdat,vardim,
     >            varmin,varmax,stag,stat)
      if (stat.ne.0) goto 998
      call getdat(cdfid,varname,time(1),0,field,stat)
      print*, 'R ',trim(varname),' ',trim(pvsrcfile)
      if (stat.ne.0) goto 998

!  Close input netcdf file
      call clscdf(cdfid,stat)
      if (stat.ne.0) goto 998

!  Set missing data value to <mdv>
      do i=1,nx
         do j=1,ny
            do k=1,nz
               if (abs(field(i,j,k)-misdat).lt.eps) then
                  field(i,j,k)=mdv
               endif
            enddo
         enddo
      enddo

      return

!  Exception handling
 998  print*,'Read_Inp: Problem with input netcdf file... Stop'
      stop

      end

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

      end

!  --------------------------------------------------------------------------------
!  Get grid parameters
!  --------------------------------------------------------------------------------

      subroutine read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,
     >                     pvsrcfile)

!  Get the grid parameters from the variable <THETA> on the input file <pvsrcfile>.
!  The grid parameters are
!     nx,ny,nz                : Number of grid points in x, y and z direction
!     xmin,ymin,zmin          : Minimum domain coordinates in x, y and z direction
!     xmax,ymax,zmax          : Maximal domain coordinates in x, y and z direction
!     dx,dy,dz                : Horizontal and vertical resolution
!  Additionally, it is checked whether the vertical grid is equally spaced. If ok,
!  the grid paramters are transformed from lon/lat to distance (in meters)

  use kind_parameters,ONLY:&
       sp

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
      real(kind=sp)           varmin(4),varmax(4),stag(4)
      real(kind=sp)           aklev(1000),bklev(1000),aklay(1000),bklay(1000)
      real(kind=sp)           dh
      character*80   csn
      integer        ndim
      integer        i

!  Get all grid parameters
      call cdfopn(pvsrcfile,cdfid,ierr)
      if (ierr.ne.0) goto 998
      call getvars(cdfid,nvars,vnam,ierr)
      if (ierr.ne.0) goto 998
      isok=0
      varname='QGPV'
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 998
      call getcfn(cdfid,csn,ierr)
      if (ierr.ne.0) goto 998
      call cdfopn(csn,cstid,ierr)
      if (ierr.ne.0) goto 998
      call getdef(cdfid,varname,ndim,misdat,vardim,varmin,varmax,
     >            stag,ierr)
      if (ierr.ne.0) goto 998
      nx=vardim(1)
      ny=vardim(2)
      nz=vardim(3)
      xmin=varmin(1)
      ymin=varmin(2)
      zmin=varmin(3)
      call getlevs(cstid,nz,aklev,bklev,aklay,bklay,ierr)
      if (ierr.ne.0) goto 998
      call getgrid(cstid,dx,dy,ierr)
      if (ierr.ne.0) goto 998
      xmax=varmax(1)
      ymax=varmax(2)
      zmax=varmax(3)
      dz=(zmax-zmin)/real(nz-1)
      call clscdf(cstid,ierr)
      if (ierr.ne.0) goto 998
      call clscdf(cdfid,ierr)
      if (ierr.ne.0) goto 998

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

      end


!  --------------------------------------------------------------------------------
!  Read refernece profile from netcdf
!  --------------------------------------------------------------------------------

      SUBROUTINE read_ref (nsqref,rhoref,thetaref,pressref,zref,
     >                     nx,ny,nz,deltax,deltay,deltaz,coriol,
     >                     pvsrcfile) 

!  Read the reference profile from file
c
!      thetaref             : Reference potential temperature (K)
!      pressref             : Reference pressure (Pa)
!      rhoref               : Reference density (kg/m^3)
!      nsqref               : Stratification (s^-1)
!      zref                 : Reference height (m)
!      nx,nny,nz            : Grid dimension in x,y,z direction
!      deltax,deltay,deltaz : Grid spacings used for calculations (m)
!      coriol               : Coriolis parameter (s^-1)
!      pvsrcfile            : Input file
  use kind_parameters,ONLY:&
       sp

  implicit   none

!  Declaration of subroutine parameters
      integer              nx,ny,nz
      real(kind=sp)                 nsqref  (0:2*nz)
      real(kind=sp)                 thetaref(0:2*nz)
      real(kind=sp)                 rhoref  (0:2*nz)
      real(kind=sp)                 pressref(0:2*nz)
      real(kind=sp)                 zref    (0:2*nz)
      real(kind=sp)                 deltax,deltay,deltaz
      real(kind=sp)                 coriol  (0:nx,0:ny)
      character*80         pvsrcfile

!  Numerical and physical parameters
      real(kind=sp)                 eps
      parameter            (eps=0.01)

!  Auxiliary variables
      integer              cdfid,stat
      integer              vardim(4)
      real(kind=sp)                 misdat
      integer              ndimin
      real(kind=sp)                 varmin(4),varmax(4),stag(4)
      integer              i,j,k,nf1
      integer              ntimes
      real(kind=sp)                 time(200)  
      character*80         vnam(100),varname
      integer              nvars
      integer              isok,ierr
      real(kind=sp)                 x(0:nx,0:ny),y(0:nx,0:ny)
      real(kind=sp)                 mean,count

!  Get grid description from topography
      call cdfopn(pvsrcfile,cdfid,stat)
      if (stat.ne.0) goto 997
      call getvars(cdfid,nvars,vnam,stat)
      if (stat.ne.0) goto 997
      isok=0
      varname='ORO'
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 997
      call getdef(cdfid,varname,ndimin,misdat,vardim,
     >            varmin,varmax,stag,stat)    
      if (stat.ne.0) goto 997
      time(1)=0.
      call gettimes(cdfid,time,ntimes,stat)  
      if (stat.ne.0) goto 997
      call clscdf(cdfid,stat)
      if (stat.ne.0) goto 997

!  Open output netcdf file
      call cdfopn(pvsrcfile,cdfid,stat)
      if (stat.ne.0) goto 997
      
!  Create the variable if necessary
      call getvars(cdfid,nvars,vnam,stat)
      if (stat.ne.0) goto 997

!  Read data from netcdf file
      isok=0
      varname='NSQREF'
      print*,'R ',trim(varname),' ',trim(pvsrcfile)
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 997
      call getdat(cdfid,varname,time(1),0,nsqref,stat)
      if (stat.ne.0) goto 997

      isok=0
      varname='RHOREF'
      print*,'R ',trim(varname),' ',trim(pvsrcfile)
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 997
      call getdat(cdfid,varname,time(1),0,rhoref,stat)
      if (stat.ne.0) goto 997
 
      isok=0
      varname='THETAREF'
      print*,'R ',trim(varname),' ',trim(pvsrcfile)
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 997
      call getdat(cdfid,varname,time(1),0,thetaref,stat)
      if (stat.ne.0) goto 997

      isok=0
      varname='PREREF'
      print*,'R ',trim(varname),' ',trim(pvsrcfile)
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 997
      call getdat(cdfid,varname,time(1),0,pressref,stat)
      if (stat.ne.0) goto 997

      isok=0
      varname='ZREF'
      print*,'R ',trim(varname),' ',trim(pvsrcfile)
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 997
      call getdat(cdfid,varname,time(1),0,zref,stat)
      if (stat.ne.0) goto 997

      isok=0
      varname='CORIOL'
      print*,'R ',trim(varname),' ',trim(pvsrcfile)
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 997
      call getdat(cdfid,varname,time(1),0,coriol,stat)
      if (stat.ne.0) goto 997

      isok=0
      varname='X'
      print*,'R ',trim(varname),' ',trim(pvsrcfile)
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 997
      call getdat(cdfid,varname,time(1),0,x,stat)
      if (stat.ne.0) goto 997

      isok=0
      varname='Y'
      print*,'R ',trim(varname),' ',trim(pvsrcfile)
      call check_varok(isok,varname,vnam,nvars)
      if (isok.eq.0) goto 997
      call getdat(cdfid,varname,time(1),0,y,stat)
      if (stat.ne.0) goto 997

!  Close  netcdf file
      call clscdf(cdfid,stat)
      if (stat.ne.0) goto 997

!  Determine the grid spacings <deltax, deltay, deltaz>
      mean=0.
      count=0.
      do i=1,nx
         do j=0,ny
            mean=mean+abs(x(i,j)-x(i-1,j))
            count=count+1.
         enddo
      enddo
      deltax=mean/count

      mean=0.
      count=0.
      do j=1,ny
         do i=0,nx
            mean=mean+abs(y(i,j)-y(i,j-1))
            count=count+1.
         enddo
      enddo
      deltay=mean/count

      mean=0.
      count=0.
      do k=1,nz-1
         mean=mean+abs(zref(k+1)-zref(k-1))
         count=count+1.
      enddo
      deltaz=mean/count

      return

!  Exception handling
 997  print*,'Read_Ref: Problem with input netcdf file... Stop'
      stop

      end


!  ********************************************************************************
!  * BOUNDARY CONDITIONS - CONSISTENCY CHECK AND ADAPTIONS                        *
!  ********************************************************************************

!  --------------------------------------------------------------------------------
!  Boundary condition 
!  --------------------------------------------------------------------------------

      subroutine combouncon(pv,nsq,rhoref,thetatop,
     >                      thetabot,thetaref,coriol,
     >                      ua,ub,va,vb,nx,ny,nz,dx,dy,dz)

!  Evaluate the consistency integrals A.7 from Rene's dissertation. This inegral
!  is a necessary condition that the von Neumann problem has a unique solution.   
!  Adjust the upper and lower boundary conditions on <thetabot> and <thetatop>, so
!  that the consitency check is ok.
  use kind_parameters,ONLY:&
       sp

  implicit none

!  Declaration of subroutine parameters
      integer              nx,ny,nz
      real(kind=sp)                 dx,dy,dz
      real(kind=sp)                 pv(0:nx,0:ny,0:nz)
      real(kind=sp)                 nsq(0:2*nz)
      real(kind=sp)                 rhoref(0:2*nz)
      real(kind=sp)                 thetatop(0:nx,0:ny)
      real(kind=sp)                 thetabot(0:nx,0:ny)
      real(kind=sp)                 thetaref(0:2*nz)
      real(kind=sp)                 coriol(0:nx,0:ny)
      real(kind=sp)                 ua(0:nx,0:nz)
      real(kind=sp)                 ub(0:nx,0:nz)
      real(kind=sp)                 va(0:ny,0:nz)
      real(kind=sp)                 vb(0:ny,0:nz)

!  Numerical and physical parameters
      real(kind=sp)                 g
      parameter            (g=9.81)

!  Auxiliary variables
      integer              i,j,k
      real(kind=sp)                 dxy,dxz,dyz,dxyz
      real(kind=sp)                 integr,denombot,denomtop,denom
      real(kind=sp)                 shifttop,shiftbot

!  Set the area and volume infinitesimals for integration
      dxy =dx*dy
      dxz =dx*dz
      dyz =dy*dz
      dxyz=dx*dy*dz

!  Init integration variables
      integr=0.

!  Inner
      do k=1,nz-1
         do i=1,nx-1
            do j=1,ny-1
               integr=integr+dxyz*rhoref(2*k)*pv(i,j,k)          
            enddo
         enddo
      enddo

!  ZY plane
      do k=1,nz-1
         do j=1,ny-1
            integr=integr+dyz*
     >             rhoref(2*k)*(dx*pv(0, j,k)+va(j,k))
!  
            integr=integr+dyz*
     >             rhoref(2*k)*(dx*pv(nx,j,k)-vb(j,k))
         enddo
      enddo

!  ZX plane
      do k=1,nz-1
         do i=1,nx-1
            integr=integr+dxz*
     >             rhoref(2*k)*(dy*pv(i,0,k)-ua(i,k))
!        
            integr=integr+dxz*
     >             rhoref(2*k)*(dy*pv(i,ny,k)+ub(i,k))           
         enddo         
      enddo

!  XY plane
      do i=1,nx-1
         do j=1,ny-1
            integr=integr+dxy*rhoref(0)*(
     >             dz*pv(i,j,0)+coriol(i,j)*g*thetabot(i,j)/
     >             (nsq(0)*thetaref(0)))
c
            integr=integr+dxy*rhoref(2*nz)*(
     >             dz*pv(i,j,nz)-coriol(i,j)*g*thetatop(i,j)/
     >             (nsq(2*nz)*thetaref(2*nz)))
c
         enddo
      enddo

!  X edges
      do i=1,nx-1
         integr=integr+dx*
     >          rhoref(0)*(dyz*pv(i,0,0)-
     >          dz*ua(i,0)+dy*coriol(i,0)*g*thetabot(i,0)/
     >          (nsq(0)*thetaref(0)))
c
         integr=integr+dx*
     >          rhoref(0)*(dyz*pv(i,ny,0)+
     >          dz*ub(i,0)+dy*coriol(i,ny)*g*thetabot(i,ny)/
     >          (nsq(0)*thetaref(0)))
c
         integr=integr+dx*
     >          rhoref(2*nz)*(dyz*pv(i,0,nz)-
     >          dz*ua(i,nz)-dy*coriol(i,0)*g*thetatop(i,0)/
     >          (nsq(2*nz)*thetaref(2*nz)))

         integr=integr+dx*
     >          rhoref(2*nz)*(dyz*pv(i,ny,nz)+
     >          dz*ub(i,nz)-dy*coriol(i,ny)*g*thetatop(i,ny)/
     >          (nsq(2*nz)*thetaref(2*nz)))
c
      enddo

!  Y edges
      do j=1,ny-1
         integr=integr+dy*
     >          rhoref(0)*(dxz*pv(0,j,0)+
     >          dz*va(j,0)+dx*coriol(0,j)*g*thetabot(0,j)/
     >          (nsq(0)*thetaref(0)))
!  
         integr=integr+dy*
     >        rhoref(0)*(dxz*pv(nx,j,0)-
     >        dz*vb(j,0)+dx*coriol(nx,j)*g*thetabot(nx,j)/
     >        (nsq(0)*thetaref(0)))
c
         integr=integr+dy*
     >        rhoref(2*nz)*(dxz*pv(0,j,nz)+
     >        dz*va(j,nz)-dx*coriol(0,j)*g*thetatop(0,j)/
     >        (nsq(2*nz)*thetaref(2*nz)))
c
         integr=integr+dy*
     >        rhoref(2*nz)*(dxz*pv(nx,j,nz)-
     >        dz*vb(j,nz)-dx*coriol(nx,j)*g*thetatop(nx,j)/
     >        (nsq(2*nz)*thetaref(2*nz)))
c
      enddo

!  Z edges
      do k=1,nz-1
         integr=integr+dz*
     >          rhoref(2*k)*(dxy*pv(0,0,k)+
     >          dy*va(0,k)-dx*ua(0,k))
!  
         integr=integr+dz*
     >          rhoref(2*k)*(dxy*pv(nx,0,k)-
     >          dy*vb(0,k)-dx*ua(nx,k))
c
         integr=integr+dz*
     >          rhoref(2*k)*(dxy*pv(0,ny,k)+
     >          dy*va(ny,k)+dx*ub(0,k))
c
         integr=integr+dz*
     >          rhoref(2*k)*(dxy*pv(nx,ny,k)-
     >          dy*vb(ny,k)+dx*ub(nx,k))
      enddo

!  Points
      integr=integr+rhoref(0)*(dxyz*pv(0,0,0)+
     >       dyz*va(0,0)-dxz*ua(0,0)+
     >       dxy*coriol(0,0)*g*thetabot(0,0)/
     >       (nsq(0)*thetaref(0)))
c
      integr=integr+rhoref(0)*(dxyz*pv(nx,0,0)-
     >       dyz*vb(0,0)-dxz*ua(nx,0)+
     >       dxy*coriol(nx,0)*g*thetabot(nx,0)/
     >       (nsq(0)*thetaref(0)))
c
      integr=integr+rhoref(0)*(dxyz*pv(0,ny,0)+
     >       dyz*va(ny,0)+dxz*ub(0,0)+
     >       dxy*coriol(0,ny)*g*thetabot(0,ny)/
     >       (nsq(0)*thetaref(0)))
!  
      integr=integr+rhoref(0)*(dxyz*pv(nx,ny,0)-
     >       dyz*vb(ny,0)+dxz*ub(nx,0)+
     >       dxy*coriol(nx,ny)*g*thetabot(nx,ny)/
     >       (nsq(0)*thetaref(0)))
c 
      integr=integr+rhoref(2*nz)*(dxyz*pv(0,0,nz)+
     >       dyz*va(0,nz)-dxz*ua(0,nz)-
     >       dxy*coriol(0,0)*g*thetatop(0,0)/
     >       (nsq(2*nz)*thetaref(2*nz)))
c
      integr=integr+rhoref(2*nz)*(dxyz*pv(nx,0,nz)-
     >       dyz*vb(0,nz)-dxz*ua(nx,nz)-
     >       dxy*coriol(nx,0)*g*thetatop(nx,0)/
     >       (nsq(2*nz)*thetaref(2*nz)))
c
      integr=integr+rhoref(2*nz)*(dxyz*pv(0,ny,nz)+
     >       dyz*va(ny,nz)+dxz*ub(0,nz)-
     >       dxy*coriol(0,ny)*g*thetatop(0,ny)/
     >       (nsq(2*nz)*thetaref(2*nz)))
c
      integr=integr+rhoref(2*nz)*(dxyz*pv(nx,ny,nz)-
     >       dyz*vb(ny,nz)+dxz*ub(nx,nz)-
     >       dxy*coriol(nx,ny)*g*thetatop(nx,ny)/
     >       (nsq(2*nz)*thetaref(2*nz)))
c

!  Get the integrals from the reference state at bottom and top 
      denombot=0.
      denomtop=0.
      do i=0,nx
         do j=0,ny
            denombot=denombot+dxy*
     >               rhoref(0)*coriol(i,j)*g/
     >               (nsq(0)*thetaref(0))
!  
            denomtop=denomtop+dxy*
     >               rhoref(2*nz)*coriol(i,j)*g/
     >               (nsq(2*nz)*thetaref(2*nz))
         enddo
      enddo
      denom=denomtop-denombot
      
!  Determine the deviation of potential temperature from reference profile
      shiftbot=0.
      shifttop=0.
      do i=0,nx
         do j=0,ny
            shifttop=shifttop+thetatop(i,j)
            shiftbot=shiftbot+thetabot(i,j)
         enddo
      enddo
      shifttop=shifttop/real((nx+1)*(ny+1))
      shiftbot=shiftbot/real((nx+1)*(ny+1))

!  Write some information about the consitency integrals
      print*,'Consistency Check for boundary' 
      print*,'       integ                      = ', integr 
      print*,'       denombot                   = ', denombot
      print*,'       denomtop                   = ', denomtop
      print*,'       denom                      = ', denom
      print*,'       theta adjustment           = ', integr/denom
      print*,'       theta shift @ top          = ', shifttop,
     >                                               thetaref(2*nz)
      print*,'       theta shift @ bot          = ', shiftbot,
     >                                               thetaref(0)

   end subroutine combouncon
      
