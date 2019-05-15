module cdfreadwrite
contains


!   --------------------------------------------------------------------------------
!   Read input fields for reference profile
!   --------------------------------------------------------------------------------

SUBROUTINE read_inp_pv_to_qgpv (field,fieldname,pvsrcfile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)

!   Read <fieldname> from netcdf file <pvsrcfile> into <field>. The grid is specified 
!   by <nx,ny,nz,dx,dy,dz,xmin,ymin,zmin>. A check is performed whether the input 
!   files are consitent with this grid. The missing data value is set to <mdv>.
  use kind_parameters,ONLY:&
       sp
  use netcdflibrary
  implicit   none
  
  !   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp),   allocatable,dimension (:,:,:)   :: field

  character*80         fieldname
  character*80         pvsrcfile
  real(kind=sp)                 dx,dy,dz,xmin,ymin,zmin
  real(kind=sp)                 mdv
  
  !   Numerical and physical parameters
  real(kind=sp)                 eps
  parameter            (eps=0.01)
  
!   Auxiliary variables
  integer              cdfid,cdfid99
  integer              vardim(4)
  real(kind=sp)                 misdat
!  real(kind=sp)                 varmin(4),varmax(4),stag(4)
  real(kind=sp),   allocatable,dimension (:) :: varmin,varmax,stag
  integer              ndimin,outid,i,j,k
  real(kind=sp)                 tmp(nx,ny,nz)
  integer              ntimes
  real(kind=sp),   allocatable,dimension (:)   :: time

  integer              nvars
  character*80         vnam(100),varname
  integer              isok
  
  !   Open the input netcdf file
  allocate(time(1))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
!  allocate(field(0:nx,0:ny,0:nz))
  call cdfopn(pvsrcfile,cdfid)

  !   Check whether needed variables are on file
  call getvars(cdfid,nvars,vnam)
  isok=0
  varname=trim(fieldname)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 998
  
!   Get the grid parameters from theta     
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)    
  time(1)=0.
  call gettimes(cdfid,time,ntimes)  
  
  !   Check whether grid parameters are consistent
  if ( (vardim(1).ne.(nx+1)).or.&
       (vardim(2).ne.(ny+1)).or. &
       (vardim(3).ne.(nz+1)).or.&
       (abs(varmin(1)-xmin).gt.eps).or.&
       (abs(varmin(2)-ymin).gt.eps).or.&
       (abs(varmin(3)-zmin).gt.eps).or.&
       (abs((varmax(1)-varmin(1))/real(vardim(1)-1)-dx).gt.eps).or.&
       (abs((varmax(2)-varmin(2))/real(vardim(2)-1)-dy).gt.eps).or.&
       (abs((varmax(3)-varmin(3))/real(vardim(3)-1)-dz).gt.eps) ) &
       then
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
  
  !   Load variables
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)
  call getdat(cdfid,varname,time,0,field)
  print*, 'R ',trim(varname),' ',trim(pvsrcfile)
  
  !   Close input netcdf file
  call clscdf(cdfid)
  
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
  
end SUBROUTINE read_inp_pv_to_qgpv


  
SUBROUTINE read_ref_pv_to_qgpv (nsqref,rhoref,thetaref,pressref,zref,nx,ny,nz,deltax,deltay,deltaz,coriol,oro,pvsrcfile) 

        !   Read the reference profile from file
        !
        !       thetaref             : Reference potential temperature (K)
        !       pressref             : Reference pressure (Pa)
        !       rhoref               : Reference density (kg/m^3)
        !       nsqref               : Stratification (s^-1)
        !       zref                 : Reference height (m)
        !       nx,nny,nz            : Grid dimension in x,y,z direction
        !       deltax,deltay,deltaz : Grid spacings used for calculations (m)
        !       coriol               : Coriolis parameter (s^-1)
        !       oro                  : Height of orography (m)
  !       pvsrcfile            : Input file
  use netcdflibrary
  implicit   none
  
  !   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp),allocatable,dimension(:)::nsqref,thetaref,rhoref,pressref,zref
  real(kind=sp)                 deltax,deltay,deltaz
  real(kind=sp),allocatable,dimension(:,:)::coriol,oro

  character*80         pvsrcfile
  
  !   Numerical and physical parameters
  real(kind=sp)                 eps
  parameter            (eps=0.01)
  
  !   Auxiliary variables
  integer              cdfid
  integer              vardim(4)
  real(kind=sp)                 misdat
  integer              ndimin
!  real(kind=sp)                 varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  integer              i,j,k,nf1
  integer              ntimes
  real(kind=sp),allocatable,dimension(:)::time

  character*80         vnam(100),varname
  integer              nvars
  integer              isok,ierr
  real(kind=sp),allocatable,dimension(:,:)::x,y
  real(kind=sp)                 mean,count

!   Get grid description from topography
  allocate(time(1))
  allocate(x(0:nx,0:ny))
  allocate(y(0:nx,0:ny))
!  allocate(coriol(0:nx,0:ny))
!  allocate(oro(0:nx,0:ny))
!  allocate(nsqref(0:2*nz))
!  allocate(thetaref(0:2*nz))
!  allocate(rhoref(0:2*nz))
!  allocate(pressref(0:2*nz))
!  allocate(zref(0:2*nz))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))

  call cdfopn(pvsrcfile,cdfid)
  

  call getvars(cdfid,nvars,vnam)
  

  isok=0
  varname='ORO'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)    
      

  time(1)=0.
  call gettimes(cdfid,time,ntimes)  
  

  call clscdf(cdfid)
  
 
  
  !   Open output netcdf file
  call cdfopn(pvsrcfile,cdfid)
  
 
  
  !   Create the variable if necessary
  call getvars(cdfid,nvars,vnam)
  
 
  
  !   Read data from netcdf file
  isok=0
  varname='NSQREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  
  if (isok.eq.0) goto 997
  call getdatRank1(cdfid,varname,time,0,nsqref)
  
 
  
  isok=0
  varname='RHOREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  
  if (isok.eq.0) goto 997
  call getdatRank1(cdfid,varname,time,0,rhoref)
  
 
  
  isok=0
  varname='THETAREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  
  if (isok.eq.0) goto 997
  call getdatRank1(cdfid,varname,time,0,thetaref)
  
 
  
  isok=0
  varname='PREREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  
 
  call getdatRank1(cdfid,varname,time,0,pressref)
  
 
  
  isok=0
  varname='ZREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  
 
  call getdatRank1(cdfid,varname,time,0,zref)
  
 
  
  isok=0
  varname='CORIOL'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  
  if (isok.eq.0) goto 997
  call getdatRank2(cdfid,varname,time,0,coriol)
  
 
  
  isok=0
  varname='ORO'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  
  if (isok.eq.0) goto 997
  call getdatRank2(cdfid,varname,time,0,oro)
  

  
  isok=0
  varname='X'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdatRank2(cdfid,varname,time,0,x)

  
  isok=0
  varname='Y'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  
  if (isok.eq.0) goto 997
  call getdatRank2(cdfid,varname,time,0,y)
  
 
  
  !   Close  netcdf file
  call clscdf(cdfid)


  
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
        mean=mean+abs(y(j,i)-y(j-1,i))
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
  
end SUBROUTINE read_ref_pv_to_qgpv


  
!   --------------------------------------------------------------------------------
!   Write input field to netcdf
!   --------------------------------------------------------------------------------

SUBROUTINE write_inp_pv_to_qgpv (field,fieldname,pvsrcfile,nx,ny,nz)
  
  !   Read <fieldname> from netcdf file <pvsrcfile> into <field>. The grid is specified 
  !   by <nx,ny,nz,dx,dy,dz,xmin,ymin,zmin>. A check is performed whether the input 
  !   files are consitent with this grid. 
  use netcdflibrary
  implicit   none
  
  !   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp),allocatable,dimension(:,:,:)::field

  character*80         fieldname
  character*80         pvsrcfile
  
  !   Auxiliary variables
  integer              cdfid,stat
  integer              vardim(4)
  real(kind=sp)                 misdat
!  real(kind=sp)                 varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  integer              ndimin,outid,i,j,k
  real(kind=sp)                 tmp(0:nx,0:ny,0:nz)
  integer              ntimes

  real(kind=sp),allocatable,dimension(:)::time
  integer              nvars
  character*80         vnam(100),varname
  integer              isok

  !   Get grid parameters from PV
  allocate(time(1))
  allocate(varmax(4))
  allocate(varmin(4))
  allocate(stag(4))
!  allocate(field(0:nx,0:ny,0:nz))
  call cdfopn(pvsrcfile,cdfid)
  call getvars(cdfid,nvars,vnam)
  isok=0
  varname='PV'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 998
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)    
  time(1)=0.
  call gettimes(cdfid,time,ntimes)  
  call clscdf(cdfid)
  
!   Save variables (write definition, if necessary)
  call cdfwopn(pvsrcfile,cdfid)
  isok=0
  varname=fieldname
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)
  endif
  call putdat(cdfid,varname,time,0,field)
  print*,'W ',trim(varname),' ',trim(pvsrcfile)
      
!   Close input netcdf file
  call clscdf(cdfid)
  
  return
  
  !   Exception handling
998 print*,'Write_Inp: Problem with input netcdf file... Stop'
  stop
  
end SUBROUTINE write_inp_pv_to_qgpv

  

  !     --------------------------------------------------------------------------------
!     Write input field to netcdf
!     --------------------------------------------------------------------------------

SUBROUTINE write_inp (field,fieldname,pvsrcfile,nx,ny,nz)

!     Read <fieldname> from netcdf file <pvsrcfile> into <field>. The grid is specified 
!     by <nx,ny,nz,dx,dy,dz,xmin,ymin,zmin>. A check is performed whether the input 
!     files are consitent with this grid. 
  use netcdflibrary
  implicit   none
  
!     Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp),allocatable,dimension(:,:,:)::field

  character*80         fieldname
  character*80         pvsrcfile

!    Auxiliary variables
  integer              cdfid,stat
  integer              vardim(4)
  real(kind=sp)                 misdat
!  real(kind=sp)                 varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  integer              ndimin,outid,i,j,k
  real(kind=sp)                 max_th
  real(kind=sp)                 tmp(0:nx,0:ny,0:nz)
  integer              ntimes
  real(kind=sp),allocatable,dimension(:)::time

  integer              nvars
  character*80         vnam(100),varname
  integer              isok

!    Get grid parameters
!  allocate(field(0:nx,0:ny,0:nz))
  allocate(time(1))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  call cdfopn(pvsrcfile,cdfid)

  call getvars(cdfid,nvars,vnam)

  isok=0
  varname='TH'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 998
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)    

  time(1)=0.
  call gettimes(cdfid,time,ntimes)  

  call clscdf(cdfid)

  
!     Save variables (write definition, if necessary)
  call cdfwopn(pvsrcfile,cdfid)

  isok=0
  varname=fieldname
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)
  endif
  call putdat(cdfid,varname,time,0,field)
  print*,'W ',trim(varname),' ',trim(pvsrcfile)

  
!    Close input netcdf file
  call clscdf(cdfid)
  
  return
  
!     Exception handling
998 print*,'Write_Inp: Problem with input netcdf file... Stop'
  stop
  
end SUBROUTINE write_inp



!     --------------------------------------------------------------------------------
!     Output to netcdf file
!     --------------------------------------------------------------------------------

SUBROUTINE write_out (psi,thetabot,thetatop,ua,ub,va,vb,nx,ny,nz,dx,dy,dz,coriol,thetaref,rhoref,pref,pvsrcfile)

  !   Write the result of the inversion to output netcdf
  !   
  !      psi               : streamm function as calculated from inversion
  !      thetabot,thetatop : potential temperature at lower and upper boundary
  !      ua,ub             : Zonal wind at western and eastern boundary
  !      va,vb             : Meridional wind at southern and northern boundary
  !      nx,ny,nz          : Grid dimension in x, y and z direction
  !      dx,dy,dz          : Grid resolution
  !      coriol            : Coriolis parameter
  !      thetaref          : Reference profile of potential temperature
  !      rhoref            : Reference profile of density
  !      pref              : Reference profile of pressure
  !      pvsrcfile         : Name of the output netcdf file
  use netcdflibrary
  implicit   none
  
  !   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp),allocatable,dimension(:,:,:)::psi

  real(kind=sp)                 thetatop(0:nx,0:ny)
  real(kind=sp)                 thetabot(0:nx,0:ny)
  real(kind=sp)                 ua(0:nx,0:nz)
  real(kind=sp)                 ub(0:nx,0:nz)
  real(kind=sp)                 va(0:ny,0:nz)
  real(kind=sp)                 vb(0:ny,0:nz)
  character*(160)       pvsrcfile
  real(kind=sp)                 dx,dy,dz
  real(kind=sp)                 thetaref(0:2*nz)
  real(kind=sp)                 rhoref(0:2*nz)
  real(kind=sp)                 pref(0:2*nz)
  real(kind=sp)                 coriol(0:nx,0:ny)
      
!   Numerical and physical parameters
  real(kind=sp)                 eps
  parameter            (eps=0.01)
  real(kind=sp)                 g
  parameter            (g=9.81)
  real(kind=sp)                 preunit
  parameter            (preunit=0.01)
  real(kind=sp)                 kappa
  parameter            (kappa=0.286)

!   Auxiliary variables
  integer              cdfid,stat
  integer              vardim(4)
  real(kind=sp)                 misdat
  integer              ndimin
  !  real(kind=sp)                 varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  integer              i,j,k,nf1,jj,kk
  real(kind=sp),allocatable,dimension(:,:,:)::tmp,pr,th



  integer              ntimes
  real(kind=sp),allocatable,dimension(:)::time

  character*80         vnam(100),varname
  integer              nvars
  integer              isok,ierr
  real(kind=sp)                 meanpsi,meancnt

!   Get grid description 

  allocate(time(1))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
!  allocate(psi(0:nx,0:ny,0:nz))
  allocate(tmp(0:nx,0:ny,0:nz))
  allocate(pr(0:nx,0:ny,0:nz))
  allocate(th(0:nx,0:ny,0:nz))
  call cdfopn(pvsrcfile,cdfid)

  call getvars(cdfid,nvars,vnam)

  isok=0
  varname='QGPV'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)

  time(1)=0.
  call gettimes(cdfid,time,ntimes)  

  call clscdf(cdfid)


!   Open output netcdf file
  call cdfwopn(pvsrcfile,cdfid)

  
  !   Create the variable if necessary
  call getvars(cdfid,nvars,vnam)

  isok=0
  varname='PSI'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)

  endif
  isok=0
  varname='U'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)

  endif
  isok=0
  varname='V'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)

  endif
  isok=0
  varname='TH'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)

  endif
  isok=0
  varname='T'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)

  endif
  isok=0
  varname='P'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)

  endif
  
!   Write stream function
  varname='PSI'
  call putdat(cdfid,varname,time,0,psi)

  print*,'W PSI      ',trim(pvsrcfile)

!   Calculate and write velocity U: keep southern and northern boundary
  do k=0,nz
     do i=0,nx
        do j=1,ny-1
           tmp(i,j,k)=(psi(i,j-1,k)-psi(i,j+1,k))/(2.*dy*coriol(i,j))
        enddo
        tmp(i,0,k) =ua(i,k)
        tmp(i,ny,k)=ub(i,k)
     enddo
  enddo
  varname='U'
  call putdat(cdfid,varname,time,0,tmp)

  print*,'W U        ',trim(pvsrcfile)

!   Calculate and write velocity V: keep western and eastern boundary
  do k=0,nz
     do j=0,ny
        do i=1,nx-1
           tmp(i,j,k)=(psi(i+1,j,k)-psi(i-1,j,k))/(2.*dx*coriol(i,j))
        enddo
        tmp(0,j,k)=va(j,k)
        tmp(nx,j,k)=vb(j,k)
     enddo
  enddo
  varname='V'
  call putdat(cdfid,varname,time,0,tmp)

  print*,'W V        ',trim(pvsrcfile)
  
!   Calculate and write potential temperature: keep lower and upper boundary
  !   Potential temperature is needed for calculation of temperature: keep it
  do i=0,nx
     do j=0,ny
        th(i,j,0)=thetabot(i,j)
        do k=1,nz-1      
           th(i,j,k)=thetaref(2*k)*(psi(i,j,k+1)-psi(i,j,k-1))/(2.*dz*g)
           
        enddo
        th(i,j,nz)=thetatop(i,j)
     enddo
  enddo
  varname='TH'
  call putdat(cdfid,varname,time,0,th)

  print*,'W TH       ',trim(pvsrcfile)

!   Calculate and write pressure: The pressure is directly proportional to the
!   streamfunction. But the streamfunction is determined only up to an additive
!   constant. Shift the streamfunction in such a way that it vanish in the mean
!   on the boundaries. Pressure is needed for calculation of temperature: keep it
  meanpsi=0.
  meancnt=0.
  do i=0,nx
     do j=0,ny
        meanpsi=meanpsi+psi(i,j,0)+psi(i,j,nz)
        meancnt=meancnt+2
     enddo
  enddo
  do i=0,nx
     do k=0,nz
        meanpsi=meanpsi+psi(i,0,k)+psi(i,ny,k)
        meancnt=meancnt+2 
     enddo
  enddo
  do j=0,ny
     do k=0,nz
        meanpsi=meanpsi+psi(0,j,k)+psi(nx,j,k)
        meancnt=meancnt+2 
     enddo
  enddo
  meanpsi=meanpsi/meancnt
  do i=0,nx
     do j=0,ny
        do k=0,nz
           kk=2*k
           pr(i,j,k)=preunit*rhoref(kk)*(psi(i,j,k)-meanpsi)
        enddo
     enddo
  enddo
  varname='P'
  call putdat(cdfid,varname,time,0,pr)

  print*,'W P        ',trim(pvsrcfile)

  !   Calculate and write temperature
  do i=0,nx
     do j=0,ny
        do k=0,nz      
           kk=2*k
           tmp(i,j,k)=((pref(kk)/1000.)**kappa) * (th(i,j,k)+kappa*thetaref(kk)*pr(i,j,k)/pref(kk))
        enddo
     enddo
  enddo
  varname='T'
  call putdat(cdfid,varname,time,0,tmp)

  print*,'W T        ',trim(pvsrcfile)

!   Close output netcdf file
  call clscdf(cdfid)

  
  return

!   Exception handling
997 print*,'Problem with output netcdf file... Stop'
  stop

end subroutine write_out

SUBROUTINE read_inp_prep(field,fieldname,pvsrcfile,nx,ny,nz,dx,dy,dz,xmin,ymin,zmin,mdv)

!     Read <fieldname> from netcdf file <pvsrcfile> into <field>. The grid is specified 
!     by <nx,ny,nz,dx,dy,dz,xmin,ymin,zmin>. A check is performed whether the input 
!     files are consitent with this grid. The missing data value is set to <mdv>.
  use netcdflibrary
  implicit   none
  
!     Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp),allocatable,dimension(:,:,:)::field

  character*80         fieldname
  character*80         pvsrcfile
  real(kind=sp)                 dx,dy,dz,xmin,ymin,zmin
  real(kind=sp)                 mdv
  
!     Numerical and physical parameters
  real(kind=sp)                 eps
  parameter            (eps=0.01)

!     Auxiliary variables
  integer              cdfid,cdfid99
  integer              vardim(4)
  real(kind=sp)                 misdat
!  real(kind=sp)                 varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  integer              ndimin,outid,i,j,k
  real(kind=sp)                 max_th
  real(kind=sp)                 tmp(nx,ny,nz)
  integer              ntimes
  real(kind=sp),allocatable,dimension(:)::time

  integer              nvars
  character*80         vnam(100),varname
  integer              isok

!    Open the input netcdf file
!  allocate(field(0:nx,0:ny,0:nz))
  allocate(time(200))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  call cdfopn(pvsrcfile,cdfid)

  
  !    Check whether needed variables are on file
  call getvars(cdfid,nvars,vnam)

  isok=0
  varname=trim(fieldname)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 998
  
!     Get the grid parameters from theta     
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)    

  time(1)=0.
  call gettimes(cdfid,time,ntimes)  

  
!     Check whether grid parameters are consistent
  if ( (vardim(1).ne.(nx+1)).or.&
       (vardim(2).ne.(ny+1)).or.&
       (vardim(3).ne.(nz+1)).or.&
       (abs(varmin(1)-xmin).gt.eps).or.&
       (abs(varmin(2)-ymin).gt.eps).or.&
       (abs(varmin(3)-zmin).gt.eps).or.&
          (abs((varmax(1)-varmin(1))/real(vardim(1)-1)-dx).gt.eps).or.&
          (abs((varmax(2)-varmin(2))/real(vardim(2)-1)-dy).gt.eps).or.&
          (abs((varmax(3)-varmin(3))/real(vardim(3)-1)-dz).gt.eps) ) &
          then
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
       
!     Load variables
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)

  call getdat(cdfid,varname,time,0,field)
  print*, 'R ',trim(varname),' ',trim(pvsrcfile)


!     Close input netcdf file
  call clscdf(cdfid)

  
  !     Set missing data value to <mdv>
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
      
  !     Exception handling
998 print*,'Read_Inp: Problem with input netcdf file... Stop'
  stop
  
end SUBROUTINE read_inp_prep


  
!   --------------------------------------------------------------------------------
!   Read the reference file
!   --------------------------------------------------------------------------------

SUBROUTINE read_ref (nsqref,rhoref,thetaref,pressref,zref,nx,ny,nz,deltax,deltay,deltaz,coriol,pvsrcfile) 
  
  !   Read the reference profile from file
  !
  !       thetaref             : Reference potential temperature (K)
  !       pressref             : Reference pressure (Pa)
  !       rhoref               : Reference density (kg/m^3)
  !       nsqref               : Stratification (s^-1)
  !       zref                 : Reference height (m)
  !       nx,nny,nz            : Grid dimension in x,y,z direction
  !       deltax,deltay,deltaz : Grid spacings used for calculations (m)
  !       coriol               : Coriolis parameter (s^-1)
  !       pvsrcfile            : Input file
  use netcdflibrary
  use kind_parameters,ONLY : &
       sp
  implicit   none

!   Declaration of subroutine parameters
  integer              nx,ny,nz
  real(kind=sp),dimension(:),allocatable::nsqref,thetaref,rhoref,pressref,zref
  real(kind=sp)        deltax,deltay,deltaz
  real(kind=sp),dimension(:,:),allocatable::coriol
  character*80         pvsrcfile

!   Numerical and physical parameters
  real(kind=sp)        eps
  parameter            (eps=0.01)
  
  !   Auxiliary variables
  integer              cdfid
  integer              vardim(4)
  real(kind=sp)        misdat
  integer              ndimin
!  real(kind=sp)(kind=sp)        varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  integer              i,j,k,nf1
  integer              ntimes
  real(kind=sp),dimension(:),allocatable::time
  character*80         vnam(100),varname
  integer              nvars
  integer              isok
  real(kind=sp),dimension(:,:),allocatable::x,y

  real(kind=sp)        mean,count

  allocate(time(1))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  !  allocate(nsqref(0:2*nz))
  ! allocate(thetaref(0:2*nz))
 ! allocate(rhoref(0:2*nz))
 ! allocate(pressref(0:2*nz))
 ! allocate(zref(0:2*nz))
 ! allocate(coriol(0:nx,0:ny))
  allocate(x(0:nx,0:ny))
  allocate(y(0:nx,0:ny))

  !   Get grid description from topography
  call cdfopn(pvsrcfile,cdfid)

  call getvars(cdfid,nvars,vnam)

  isok=0
  varname='ORO'
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdef(cdfid,varname,ndimin,misdat,vardim,varmin,varmax,stag)
  

  time(1)=0.
  call gettimes(cdfid,time,ntimes)  

  call clscdf(cdfid)


!   Open output netcdf file
  call cdfopn(pvsrcfile,cdfid)

  
  !   Create the variable if necessary
  call getvars(cdfid,nvars,vnam)

  
  !   Read data from netcdf file
  isok=0
  varname='NSQREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdatRank1(cdfid,varname,time,0,nsqref)


  isok=0
  varname='RHOREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdatRank1(cdfid,varname,time,0,rhoref)

 
  isok=0
  varname='THETAREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdatRank1(cdfid,varname,time,0,thetaref)


  isok=0
  varname='PREREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdatRank1(cdfid,varname,time,0,pressref)

  
  isok=0
  varname='ZREF'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdatRank1(cdfid,varname,time,0,zref)

  
  isok=0
  varname='CORIOL'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdatRank2(cdfid,varname,time,0,coriol)


  isok=0
  varname='X'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdatRank2(cdfid,varname,time,0,x)

  
  isok=0
  varname='Y'
  print*,'R ',trim(varname),' ',trim(pvsrcfile)
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) goto 997
  call getdatRank2(cdfid,varname,time,0,y)

  
!   Close  netcdf file
  call clscdf(cdfid)

  
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
  
end subroutine read_ref

  

end module cdfreadwrite
