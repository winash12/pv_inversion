PROGRAM main_calc_qgpv

!   ********************************************************************************
!   * CALCULATE QUASI-GEOSTROPHIC PV ACCORDING TO ITS DEFINITION                   *
!   * Michael Sprenger / Summer, Autumn 2006                                       *
!   ********************************************************************************

!   --------------------------------------------------------------------------------
!   Declaration of variables, parameters, externals and common blocks
!   --------------------------------------------------------------------------------
  use kind_parameters,ONLY:&
       sp
  implicit none

!   Input and output file
  character*80   diagfile
  character*80   referfile
  
  !   Grid parameters
  integer        nx,ny,nz
  real(kind=sp)           xmin,ymin,zmin,xmax,ymax,zmax
  real(kind=sp)           dx,dy,dz
  real(kind=sp)           deltax,deltay,deltaz
  real(kind=sp)           mdv
  real(kind=sp),          allocatable,dimension (:,:) :: coriol
  
!   Reference state 
  real(kind=sp),   allocatable,dimension (:) :: nsqref
  real(kind=sp),   allocatable,dimension (:) :: thetaref
  real(kind=sp),   allocatable,dimension (:) :: rhoref
  real(kind=sp),   allocatable,dimension (:) :: pressref
  real(kind=sp),   allocatable,dimension (:) :: zref
  
  !   3d fields for calculation of qgPV 
  real(kind=sp),allocatable,dimension (:,:,:) :: qgpv,th,uu,vv
  
  !   Auxiliary variables
  integer      i,j,k
  integer      stat
  character*80 varname
  
  !   --------------------------------------------------------------------------------
  !   Input
  !   --------------------------------------------------------------------------------
  
  print*,'********************************************************'
  print*,'* CALC_QGPV                                            *'
  print*,'********************************************************'
  
  !   Read parameter file
  open(10,file='fort.10')
  read(10,*) diagfile
  read(10,*) referfile
  close(10) 
  print*
  print*,trim(diagfile)
  print*,trim(referfile)
  print*
  
  !   Get lat/lon gid parameters from input file
  call read_dim (nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv,diagfile)
  print*,'Read_Dim: nx,ny,nz         = ',nx,ny,nz
  print*,'          dx,dy,dz         = ',dx,dy,dz
  print*,'          xmin,ymin,zmin   = ',xmin,ymin,zmin
  print*,'          mdv              = ',mdv
  print*
  
  !   Count from 0, not from 1: consistent with <inv_cart.f>.
  nx=nx-1
  ny=ny-1
  nz=nz-1
  
  !   Allocate memory for Coriolis parameters
  allocate(coriol(0:nx,0:ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating coriol'
  
  !   Allocate memory for reference profile
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
         
  !   Allocate memory for 3d fields
  allocate(qgpv(0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating qgpv'
  allocate(uu(0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating uu'
  allocate(vv(0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating vv'
  allocate(th(0:nx,0:ny,0:nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating th'
  
  
  !   --------------------------------------------------------------------------------
!   Calculate the qgPV from definition and put it onto file
!   --------------------------------------------------------------------------------

!   Read data from file
  varname='U'
  call read_inp (uu,varname,diagfile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='V'
  call read_inp (vv,varname,diagfile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)
  varname='TH'
  call read_inp (th,varname,diagfile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)

!   Read reference profile and grid parameters
  call read_ref (nsqref,rhoref,thetaref,pressref,zref,nx,ny,nz,deltax,deltay,deltaz,coriol,referfile)
  deltax=1000.*deltax
  deltay=1000.*deltay
  print*,'Deltax,deltay,deltaz =',deltax,deltay,deltaz
  
  !   Calculate qgPV
  print*,'C qgPV'
  call calc_qgpv (qgpv,uu,vv,th,&
       rhoref,pressref,nsqref,thetaref,coriol,&
       nx,ny,nz,deltax,deltay,deltaz,mdv)
  
  !   Write result to netcdf file
  varname='QGPV_DIA'
  call write_inp (qgpv,varname,diagfile,nx,ny,nz)
  
end PROGRAM main_calc_qgpv

!   ********************************************************************************
!   * NETCDF OUTPUT                                                                *
!   ********************************************************************************

!   --------------------------------------------------------------------------------
!   Write input field to netcdf
!   --------------------------------------------------------------------------------

SUBROUTINE write_inp (field,fieldname,pvsrcfile,nx,ny,nz)
        
  !   Read <fieldname> from netcdf file <pvsrcfile> into <field>. The grid is specified 
  !   by <nx,ny,nz,dx,dy,dz,xmin,ymin,zmin>. A check is performed whether the input 
  !   files are consitent with this grid. 
  
  implicit   none
  
  !   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp)                 field (0:nx,0:ny,0:nz)
  character*80         fieldname
  character*80         pvsrcfile
  
  !   Auxiliary variables
  integer              cdfid,stat
  integer              vardim(4)
  real(kind=sp)                 misdat
  real(kind=sp)                 varmin(4),varmax(4),stag(4)
  integer              ndimin,outid,i,j,k
  real(kind=sp)                 max_th
  real(kind=sp)                 tmp(0:nx,0:ny,0:nz)
  integer              ntimes
  real(kind=sp)                 time(200)       
  integer              nvars
  character*80         vnam(100),varname
  integer              isok
  
  !   Get grid parameters from THETA
  call cdfopn(pvsrcfile,cdfid,stat)

  call getvars(cdfid,nvars,vnam,stat)

  isok=0
  varname='TH'
  call check_varok(isok,varname,vnam,nvars)

  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag,stat)    

  time(1)=0.
  call gettimes(cdfid,time,ntimes,stat)  

  call clscdf(cdfid,stat)


!   Save variables (write definition, if necessary)
  call cdfwopn(pvsrcfile,cdfid,stat)

  isok=0
  varname=fieldname
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndimin,misdat,vardim,&
          varmin,varmax,stag,stat)
  endif
  call putdat(cdfid,varname,time(1),0,field,stat)
  print*,'W ',trim(varname),' ',trim(pvsrcfile)

  
  !   Close input netcdf file
  call clscdf(cdfid,stat)

      
  return

!   Exception handling
 998  print*,'Write_Inp: Problem with input netcdf file... Stop'
  stop
  
end SUBROUTINE write_inp

!   --------------------------------------------------------------------------------
!   Read input fields 
!   --------------------------------------------------------------------------------

SUBROUTINE read_inp (field,fieldname,pvsrcfile,&
  nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)

!   Read <fieldname> from netcdf file <pvsrcfile> into <field>. The grid is specified 
!   by <nx,ny,nz,dx,dy,dz,xmin,ymin,zmin>. A check is performed whether the input 
!   files are consitent with this grid. The missing data value is set to <mdv>.
  
  implicit   none
  
  !   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp)                 field(0:nx,0:ny,0:nz)
  character*80         fieldname
  character*80         pvsrcfile
  real(kind=sp)                 dx,dy,dz,xmin,ymin,zmin
  real(kind=sp)                 mdv
  
  !   Numerical and physical parameters
  real(kind=sp)                 eps
  parameter            (eps=0.01)
  
  !   Auxiliary variables
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
  
  !   Open the input netcdf file
  call cdfopn(pvsrcfile,cdfid,stat)

  
  !   Check whether needed variables are on file
  call getvars(cdfid,nvars,vnam,stat)

  isok=0
  varname=trim(fieldname)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 998

!   Get the grid parameters from theta     
  call getdef(cdfid,varname,ndimin,misdat,vardim,&
       varmin,varmax,stag,stat)    

  time(1)=0.
  call gettimes(cdfid,time,ntimes,stat)  
!   Check whether grid parameters are consistent
  if ( (vardim(1).ne.(nx+1)).or.&
       (vardim(2).ne.(ny+1)).or.&
       (vardim(3).ne.(nz+1)).or.&
       (abs(varmin(1)-xmin).gt.eps).or.&
       (abs(varmin(2)-ymin).gt.eps).or.&
       (abs(varmin(3)-zmin).gt.eps).or.&
       (abs((varmax(1)-varmin(1))/real(kind=sp)(vardim(1)-1)-dx).gt.eps).or.&
       (abs((varmax(2)-varmin(2))/real(kind=sp)(vardim(2)-1)-dy).gt.eps).or.&
       (abs((varmax(3)-varmin(3))/real(kind=sp)(vardim(3)-1)-dz).gt.eps) ) &
       then
     print*,'Input grid inconsitency...'
     print*,'  Nx      = ',vardim(1),nx+1
     print*,'  Ny      = ',vardim(2),ny+1
     print*,'  Nz      = ',vardim(3),nz+1
     print*,'  Varminx = ',varmin(1),xmin
     print*,'  Varminy = ',varmin(2),ymin
     print*,'  Varminz = ',varmin(3),zmin
     print*,'  Dx      = ',(varmax(1)-varmin(1))/real(kind=sp)(nx-1),dx
     print*,'  Dy      = ',(varmax(2)-varmin(2))/real(kind=sp)(ny-1),dy
     print*,'  Dz      = ',(varmax(3)-varmin(3))/real(kind=sp)(nz-1),dz
     goto 998
  endif
  
  !   Load variables
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag,stat)
  call getdat(cdfid,varname,time(1),0,field,stat)
  print*, 'R ',trim(varname),' ',trim(pvsrcfile)
  
  !   Close input netcdf file
  call clscdf(cdfid,stat)
  
  
  !   Set missing data value to <mdv>
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
  
  !   Exception handling
998 print*,'Read_Inp: Problem with input netcdf file... Stop'
  stop
  
end SUBROUTINE read_inp

!   --------------------------------------------------------------------------------
!   Read refernece profile from netcdf
!   --------------------------------------------------------------------------------

SUBROUTINE read_ref (nsqref,rhoref,thetaref,pressref,zref,&
     nx,ny,nz,deltax,deltay,deltaz,coriol,&
     pvsrcfile) 
  
!   Read the reference profile from file

  !       thetaref             : Reference potential temperature (K)
!       pressref             : Reference pressure (Pa)
  !       rhoref               : Reference density (kg/m^3)
!       nsqref               : Stratification (s^-1)
!       zref                 : Reference height (m)
!       nx,nny,nz            : Grid dimension in x,y,z direction
!       deltax,deltay,deltaz : Grid spacings used for calculations (m)
!       coriol               : Coriolis parameter (s^-1)
!       pvsrcfile            : Input file

  implicit   none
      
  !   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp)                 nsqref  (0:2*nz)
  real(kind=sp)                 thetaref(0:2*nz)
  real(kind=sp)                 rhoref  (0:2*nz)
  real(kind=sp)                 pressref(0:2*nz)
  real(kind=sp)                 zref    (0:2*nz)
  real(kind=sp)                 deltax,deltay,deltaz
  real(kind=sp)                 coriol  (0:nx,0:ny)
  character*80         pvsrcfile
  
  !   Numerical and physical parameters
  real(kind=sp)                 eps
  parameter            (eps=0.01)
  
  !   Auxiliary variables
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
  
  !   Get grid description from topography
  call cdfopn(pvsrcfile,cdfid,stat)
  call getvars(cdfid,nvars,vnam,stat)
  isok=0
  varname='ORO'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdef(cdfid,varname,ndimin,misdat,vardim,&
                 varmin,varmax,stag)    
  time(1)=0.
  call gettimes(cdfid,time,ntimes,stat)  
  call clscdf(cdfid,stat)


!   Open output netcdf file
  call cdfopn(pvsrcfile,cdfid,stat)
  !   Create the variable if necessary
  call getvars(cdfid,nvars,vnam,stat)

  
!   Read data from netcdf file
  isok=0
  varname='NSQREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdat(cdfid,varname,time(1),0,nsqref,stat)
  isok=0
  varname='RHOREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdat(cdfid,varname,time(1),0,rhoref,stat)
  isok=0
  varname='THETAREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdat(cdfid,varname,time(1),0,thetaref,stat)
  
  
  isok=0
  varname='PREREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdat(cdfid,varname,time(1),0,pressref,stat)


  isok=0
  varname='ZREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdat(cdfid,varname,time(1),0,zref,stat)
  isok=0
  varname='CORIOL'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdat(cdfid,varname,time(1),0,coriol,stat)
  isok=0
  varname='X'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdat(cdfid,varname,time(1),0,x,stat)


  isok=0
  varname='Y'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdat(cdfid,varname,time(1),0,y,stat)

!   Close  netcdf file
  call clscdf(cdfid,stat)


!   Determine the grid spacings <deltax, deltay, deltaz>
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
      
!   Exception handling
997 print*,'Read_Ref: Problem with input netcdf file... Stop'
  stop
  
end SUBROUTINE read_ref
      

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

  implicit none
  
  !   Declaration of subroutine parameters
  character*80   pvsrcfile
  integer        nx,ny,nz
  real(kind=sp)           dx,dy,dz
  real(kind=sp)           xmin,ymin,zmin,xmax,ymax,zmax
  real(kind=sp)           mdv
  
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
  real(kind=sp)           misdat
  real(kind=sp)           varmin(4),varmax(4),stag(4)
  real(kind=sp)           aklev(1000),bklev(1000),aklay(1000),bklay(1000)
  real(kind=sp)           dh
  character*80   csn
  integer        ndim
  integer        i
  
  !   Get all grid parameters
  call cdfopn(pvsrcfile,cdfid,ierr)
  if (ierr.ne.0) goto 998
  call getvars(cdfid,nvars,vnam,ierr)
  if (ierr.ne.0) goto 998
  isok=0
  varname='TH'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 998
  call getcfn(cdfid,csn,ierr)
  if (ierr.ne.0) goto 998
  call cdfopn(csn,cstid,ierr)
  if (ierr.ne.0) goto 998
  call getdef(cdfid,varname,ndim,misdat,vardim,varmin,varmax,stag)
  
  nx=vardim(1)
  ny=vardim(2)
  nz=vardim(3)
  xmin=varmin(1)
  ymin=varmin(2)
  zmin=varmin(3)
  call getlevs(cstid,nz,aklev,bklev,aklay,bklay)

  call getgrid(cstid,dx,dy,ierr)
  
  xmax=varmax(1)
  ymax=varmax(2)
  zmax=varmax(3)
  dz=(zmax-zmin)/real(kind=sp)(nz-1)
  call clscdf(cstid,ierr)
  
  call clscdf(cdfid,ierr)

  
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



!   --------------------------------------------------------------------------------
!   Calculate qgPV  from wind and theta 
!   --------------------------------------------------------------------------------       

subroutine calc_qgpv (qgpv,uu,vv,th,&
     rhoref,pressref,nsqref,thetaref,coriol,&
     nx,ny,nz,deltax,deltay,deltaz,mdv)
  
  !   Calculate qgPV from  wind and potential temperature according to
  !   equation 2.9 p 16 Thesis Rene Fehlmann. Note a cartesian grid with
  !   equidistant grid spacings <deltax,deltay,deltaz> is assumend. No
  !   'correction' is made for spherical geoemtry.
  
  implicit none
 
  !   Declaration of subroutine parameters
  integer   nx,ny,nz
  real(kind=sp)      qgpv(0:nx,0:ny,0:nz)
  real(kind=sp)      uu(0:nx,0:ny,0:nz)
  real(kind=sp)      vv(0:nx,0:ny,0:nz)
  real(kind=sp)      th(0:nx,0:ny,0:nz)
  real(kind=sp)      rhoref(0:2*nz)
  real(kind=sp)      nsqref(0:2*nz)
  real(kind=sp)      thetaref(0:2*nz)
  real(kind=sp)      pressref(0:2*nz)
  real(kind=sp)      deltax,deltay,deltaz
  real(kind=sp)      mdv
  real(kind=sp)      coriol(0:nx,0:ny)
  
  !   Numerical epsilon and physical constants
  real(kind=sp)       g
  parameter  (g=9.81)
  real(kind=sp)       eps
  parameter  (eps=0.01)
  real(kind=sp)       scale
  parameter  (scale=1e6)
  
  !   Auxiliary variables
  integer  i,j,k
  integer  kk,jj
  real(kind=sp)     dvdx(0:nx,0:ny,0:nz)
  real(kind=sp)     dudy(0:nx,0:ny,0:nz)
  real(kind=sp)     dtdz(0:nx,0:ny,0:nz)
  real(kind=sp)     t1,t2
  
  !   Calculate horizontal derivatives dudy and dvdx of velocity
  do k=0,nz
     do j=0,ny
        do i=0,nx
           
           !            Calculate dudy
           if (j.eq.0) then
              if ( (abs(uu(i,1,k)-mdv).gt.eps).and.&
                   (abs(uu(i,0,k)-mdv).gt.eps)) then
                 dudy(i,0,k)=(uu(i,1,k)-uu(i,0,k))/deltay
              else
                 dudy(i,0,k)=mdv
              endif
           elseif (j.eq.ny) then
              if ( (abs(uu(i,ny,  k)-mdv).gt.eps).and.&
                   (abs(uu(i,ny-1,k)-mdv).gt.eps)) then
                 dudy(i,ny,k)=(uu(i,ny,k)-uu(i,ny-1,k))/deltay
              else
                 dudy(i,ny,k)=mdv
              endif
           else
              if ( (abs(uu(i,j+1,k)-mdv).gt.eps).and.&
                   (abs(uu(i,j-1,k)-mdv).gt.eps)) then
                 dudy(i,j,k)=(uu(i,j+1,k)-uu(i,j-1,k))/(2.*deltay)
              else
                 dudy(i,j,k)=mdv
              endif
           endif
           
           !            Calculate dvdx
           if (i.eq.0) then
              if ((abs(vv(1,j,k)-mdv).gt.eps).and.&
                   (abs(vv(0,j,k)-mdv).gt.eps)) then
                 dvdx(0,j,k)=(vv(1,j,k)-vv(0,j,k))/deltax
              else
                 dvdx(0,j,k)=mdv
              endif
           elseif (i.eq.nx) then
              if ((abs(vv(nx,  j,k)-mdv).gt.eps).and. &
                   (abs(vv(nx-1,j,k)-mdv).gt.eps)) then
                 dvdx(nx,j,k)=(vv(nx,j,k)-vv(nx-1,j,k))/deltax
              else
                 dvdx(nx,j,k)=mdv
              endif
           else
              if ((abs(vv(i+1,j,k)-mdv).gt.eps).and.&
                   (abs(vv(i-1,j,k)-mdv).gt.eps)) then
                 dvdx(i,j,k)=(vv(i+1,j,k)-vv(i-1,j,k))/(2.*deltax)
              else
                 dvdx(i,j,k)=mdv
              endif
           endif
           
        enddo
     enddo
  enddo

  !   Calculate vertical derivative of potential temperature
  do i=0,nx
     do j=0,ny
        do k=0,nz
           
           if (k.eq.0) then
              if ((abs(th(i,j,2)-mdv).gt.eps).and.&
                   (abs(th(i,j,1)-mdv).gt.eps)) then
                 t1=rhoref(2)*th(i,j,1)/thetaref(2)/nsqref(2)*g
                 t2=rhoref(0)*th(i,j,0)/thetaref(0)/nsqref(0)*g
                 dtdz(i,j,0)=(t1-t2)/deltaz
              else
                 dtdz(i,j,0)=mdv
              endif
           else if (k.eq.nz) then
              if ((abs(th(i,j,nz  )-mdv).gt.eps).and.
                 (abs(th(i,j,nz-1)-mdv).gt.eps)) then
                 kk=2*nz
                 t1=rhoref(kk  )*th(i,j,nz  )/thetaref(kk)/nsqref(kk)*g
                 t2=rhoref(kk-2)*th(i,j,nz-1)/thetaref(kk-2)/nsqref(kk-2)*9.8
                 dtdz(i,j,nz)=(t1-t2)/deltaz
              else
                 dtdz(i,j,nz)=mdv
              endif
           else
              if ((abs(th(i,j,k+1)-mdv).gt.eps).and. &
                   (abs(th(i,j,k  )-mdv).gt.eps).and.&
                   (abs(th(i,j,k-1)-mdv).gt.eps)) then
                 kk=2*k
                 t1=rhoref(kk+1)*(th(i,j,k)+th(i,j,k+1))/2./thetaref(kk+1)/nsqref(kk+1)*g
                 t2=rhoref(kk-1)*(th(i,j,k)+th(i,j,k-1))/2./thetaref(kk-1)/nsqref(kk-1)*g
                 dtdz(i,j,k)=(t1-t2)/deltaz
              else
                 dtdz(i,j,k)=mdv
              endif
           endif
           
        enddo
     enddo
  enddo
  
!   Calculate qgPV
  do i=0,nx
     do j=0,ny
        do k=0,nz
           
           kk=2*k
           
           if ((abs(dudy(i,j,k)-mdv).gt.eps).and. &
                (abs(dvdx(i,j,k)-mdv).gt.eps).and.&
                (abs(dtdz(i,j,k)-mdv).gt.eps)) then
              
              qgpv(i,j,k)=-dudy(i,j,k)+dvdx(i,j,k)+ coriol(i,j)*dtdz(i,j,k)/rhoref(kk)
           else
              qgpv(i,j,k)=mdv
           endif
           
        enddo
     enddo
  enddo
  
end subroutine calc_qgpv

