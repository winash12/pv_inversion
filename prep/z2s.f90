Program z2s
  
  !  Calculate secondary fields on z levels
  !  Michael Sprenger / Summer 2006
  !  Modified by Ashwin Dinakar 2017,2018
  !  Converted code from Fortran 77 to fortran 99 and fortran netcdf

  use kind_parameters,ONLY:&
       sp
  use netcdflibrary
  implicit none
  
  !  ---------------------------------------------------------------
  !  Declaration of variables
  !  ---------------------------------------------------------------
  
  !  Variables for output Z file  : height level
  character*80 cfn
!  real(kind=sp)         varmin(4),varmax(4),stag(4)
  real(kind=sp),   allocatable,dimension (:) :: varmin,varmax,stag
  integer      vardim(4)
  real(kind=sp)         mdv
  integer      ndim
  integer      nx,ny,nz
  real(kind=sp)         xmin,xmax,ymin,ymax,dx,dy
  integer      ntimes
  real(kind=sp)         aklev(1000),bklev(1000)
  real(kind=sp)         aklay(1000),bklay(1000)
  real(kind=sp),dimension(:),allocatable::         time
  real(kind=sp)         pollon,pollat
  integer      idate(5)
  integer      nfields
  character*80 field_nam(100)

  real(kind=sp),allocatable, dimension (:,:,:,:) :: field_dat
  real(kind=sp),allocatable, dimension (:,:,:)   :: z3
  real(kind=sp),allocatable, dimension (:,:)     :: x2,y2,f2
  real(kind=sp),allocatable, dimension (:,:,:)   :: out
  real(kind=sp),allocatable, dimension (:,:,:)   :: inp
  integer      nvars
  character*80 vnam(100),varname
  integer      isok
  integer      cdfid,cstid
  
!  Parameter file
  character*80 fieldname
  integer      nfilt
  character*80 ofn,gri
  integer ios
  !  Auxiliary variables
  integer      ierr,stat
  integer      i,j,k,n
  real(kind=sp),allocatable, dimension (:,:)   :: out2
  character*80 vnam2(100)
  integer      nvars2
  
  !  ---------------------------------------------------------------
  !  Preparations
  !  ---------------------------------------------------------------
  allocate(time(1))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  print*,'*********************************************************'
  print*,'* z2s                                                   *'
  print*,'*********************************************************'
  
  !  Read parameter file
  open(10,file='fort.10')
  read(10,*,IOSTAT=ios) fieldname
  
  read(10,*,IOSTAT=ios) ofn

  
  read(10,*,IOSTAT=ios) gri

  
  read(10,*,IOSTAT=ios) nfilt

  close(10)
  
  !  Get grid description for Z file : height level
  call cdfopn(ofn,cdfid)
  call getcfn(cdfid,cfn)

  call cdfopn(cfn,cstid)

  call getdef(cdfid,'T',ndim,mdv,vardim,&
       varmin,varmax,stag)
  nx  =vardim(1)
  print *,'nx is  ' , nx
  ny  =vardim(2)
  print *,'ny is  ' , ny
  nz  =vardim(3)
  print *,'nz is  ' , nz
  xmin=varmin(1)
  ymin=varmin(2)
  call getlevs(cstid,nz,aklev,bklev,aklay,bklay)
  call getgrid(cstid,dx,dy)
  xmax=xmin+real(nx-1)*dx
  ymax=ymin+real(ny-1)*dy
  call gettimes(cdfid,time,ntimes)

  call getstart(cstid,idate)
  call getpole(cstid,pollon,pollat)
  call getvars(cdfid,nvars,vnam)
  call clscdf(cstid)
  call clscdf(cdfid)
  
!  Get a list of all variables on the GRID file
  if (trim(gri).ne.trim(ofn)) then
     call cdfopn(gri,cdfid)
     call getvars(cdfid,nvars2,vnam2)
     do i=1,nvars2
        vnam(nvars+i)=vnam2(i)
     enddo
     nvars=nvars+nvars2
  endif
  
  !  Check whether calculation of <fieldname> is supported
  if ( (fieldname.ne.'TH' ).and.&
       (fieldname.ne.'NSQ')  .and.&
       (fieldname.ne.'PV' )   .and.&
       (fieldname.ne.'RHO')) then
     print*,'Variable not supported ',trim(fieldname)
     stop
  endif
  
!  Set dependencies
  if (fieldname.eq.'TH') then
     nfields=2
     field_nam(1)='T'
     field_nam(2)='P'
  else if (fieldname.eq.'RHO') then
     nfields=2
     field_nam(1)='T'
     field_nam(2)='P'
  else if (fieldname.eq.'NSQ') then
     nfields=2
     field_nam(1)='T'
     field_nam(2)='Q'
  else if (fieldname.eq.'PV') then
     nfields=4
     field_nam(1)='T'
     field_nam(2)='P'
     field_nam(3)='U'
     field_nam(4)='V'
  endif

  !  Allocate memory
  allocate(field_dat(nfields,nx,ny,nz),stat=stat)
  if (stat.ne.0) print*,'error allocating field_dat'
  allocate(out(nx,ny,nz),stat=stat)
  if (stat.ne.0) print*,'error allocating out'
  allocate(inp(nx,ny,nz),stat=stat)
  if (stat.ne.0) print*,'error allocating inp'
  allocate(z3(nx,ny,nz),stat=stat)
  if (stat.ne.0) print*,'error allocating z3'
  allocate(x2(nx,ny),stat=stat)
  if (stat.ne.0) print*,'error allocating x2'
  allocate(y2(nx,ny),stat=stat)
  if (stat.ne.0) print*,'error allocating y2'
  allocate(f2(nx,ny),stat=stat)
  if (stat.ne.0) print*,'error allocating f2'
  
  !  Allocate auxiliary fields
  allocate(out2(nx,ny),stat=stat)
  if (stat.ne.0) print*,'error allocating out2'
  
  !  Read X grid
  varname='X'
  isok=0
  call check_varok (isok,varname,vnam,nvars)
  if (isok.eq.0) then
     print*,'Variable ',trim(varname),' missing... Stop'
     stop
  endif
  call cdfopn(gri,cdfid)
  
  call getdatRank2(cdfid,varname,time,0,x2)
  print*,'R ',trim(varname),' ',trim(gri)
  call clscdf(cdfid)
  
  !  Read Y grid
  varname='Y'
  isok=0
  call check_varok (isok,varname,vnam,nvars)
  if (isok.eq.0) then
     print*,'Variable ',trim(varname),' missing... Stop'
     stop
  endif
  call cdfopn(gri,cdfid)
  call getdatRank2(cdfid,varname,time,0,y2)
  print*,'R ',trim(varname),' ',trim(gri)
  call clscdf(cdfid)
  
  !  Read Coriolis parametzer
  varname='CORIOL'
  isok=0
  call check_varok (isok,varname,vnam,nvars)
  if (isok.eq.0) then
     print*,'Variable ',trim(varname),' missing... Stop'
     stop
  endif
  call cdfopn(gri,cdfid)
  call getdatRank2(cdfid,varname,time,0,f2)
  print*,'R ',trim(varname),' ',trim(gri)
  call clscdf(cdfid)
  
  !  Init height levels
  do i=1,nx
     do j=1,ny
        do k=1,nz
           z3(i,j,k)=aklay(k)
        enddo
     enddo
  enddo
  
  !  Load needed variables
  do n=1,nfields
     
     !     Check whether variable is available on file
     varname=field_nam(n)
     isok=0
     call check_varok (isok,varname,vnam,nvars)
     
     !     Variable is available on file
     if (isok.eq.1) then
        
        call cdfopn(ofn,cdfid)
        call getdat(cdfid,varname,time,0,inp)
        print*,'R ',trim(varname),' ',trim(ofn)
        call clscdf(cdfid)
        
        do i=1,nx
           do j=1,ny
              do k=1,nz
                 field_dat(n,i,j,k)=inp(i,j,k)
              enddo
           enddo
        enddo
        
     else            
        print*,'Variable missing : ',trim(varname)
        stop
     endif
     
  enddo
  
  !  Change unit of pressure from hPa to Pa
  n=0
  do i=1,nfields
     if (field_nam(i).eq.'P') n=i
  enddo
  if (n.ne.0) then
     do i=1,nx
        do j=1,ny
           do k=1,nz
              field_dat(n,i,j,k)=100.*field_dat(n,i,j,k)
           enddo
        enddo
     enddo
  endif
  
  !  ----------------------------------------------------------------
  !  Calculations
  !  ----------------------------------------------------------------
  
  !  Call to the defining routines
  if (fieldname.eq.'RHO') then
     
     print*,'C ',trim(fieldname)
     call calc_rho (out,  &                               ! RHO
          field_dat(1,:,:,:),&                  ! T
          field_dat(2,:,:,:), &                 ! P
          nx,ny,nz,mdv)
     
  else if (fieldname.eq.'TH') then
     
     print*,'C ',trim(fieldname)
     call calc_theta (out, &                              ! TH
          field_dat(1,:,:,:),&                ! T
          field_dat(2,:,:,:), &               ! P
          nx,ny,nz,mdv)
     
  else if (fieldname.eq.'NSQ') then
     
     print*,'C ',trim(fieldname)
     call calc_nsq (out,&                                 ! NSQ
          field_dat(1,:,:,:),&                  ! T
          field_dat(2,:,:,:), &                 ! Q
          z3,&                                  ! Z      
          nx,ny,nz,mdv)
     
  else if (fieldname.eq.'PV') then
     
     print*,'C ',trim(fieldname)
     
     call calc_pv (out,&                                  ! PV
          field_dat(1,:,:,:),&                   ! T
          field_dat(2,:,:,:),&                   ! P
          field_dat(3,:,:,:),&                   ! U
          field_dat(4,:,:,:),&                   ! V
          z3,x2,y2,f2,&                          ! Z,X,Y,CORIOL
          nx,ny,nz,mdv)
     
  endif
  
  !   Horizontal filtering of the resulting fields
  print*,'F ',trim(fieldname)
  do k=1,nz
     
     do i=1,nx
        do j=1,ny
           out2(i,j)=out(i,j,k)
        enddo
     enddo
     do n=1,nfilt     
        call filt2d (out2,out2,nx,ny,1.,mdv,0,0,0,0)
     enddo
     
     do i=1,nx
        do j=1,ny
           out(i,j,k)=out2(i,j)
        enddo
     enddo
     
  enddo
  
!  ----------------------------------------------------------------
!  Save result onto netcdf file
!  ----------------------------------------------------------------
       
  call cdfwopn(ofn,cdfid)
  isok=0
  varname=fieldname
  call check_varok(isok,varname,vnam,nvars)
  if (isok.eq.0) then
     call putdef(cdfid,varname,ndim,mdv,vardim,&
          varmin,varmax,stag)
  endif
  call putdat(cdfid,varname,time,0,out)
  print*,'W ',trim(varname),' ',trim(ofn)
  call clscdf(cdfid)
  
  
  !  ----------------------------------------------------------------
  !  Exception handling
  !  ----------------------------------------------------------------
  stop
  
end Program z2s

       

!  ****************************************************************
!  * SUBROUTINE SECTION: AUXILIARY ROUTINES                       *
!  ****************************************************************     

!  ----------------------------------------------------------------
!  Check whether variable is found on netcdf file
!  ----------------------------------------------------------------

subroutine check_varok (isok,varname,varlist,nvars)
  
  !  Check whether the variable <varname> is in the list <varlist(nvars)>. 
  !  If this is the case, <isok> is incremented by 1. Otherwise <isok> 
  !  keeps its value.
  
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


!  ****************************************************************
!  * SUBROUTINE SECTION: CALCULATE SECONDARY FIELDS               *
!  ****************************************************************

!  -----------------------------------------------------------------
!  Calculate density
!  -----------------------------------------------------------------

subroutine calc_rho (rho3,t3,p3,nx,ny,nz,mdv)

!  Calculate the density <rho3> (in kg/m^3) from temperature <t3> 
!  (in deg C) and pressure <p3> (in hPa).
  use kind_parameters,ONLY:&
       sp

  implicit none
      
  !  Declaration of subroutine parameters
  integer   nx,ny,nz
  real(kind=sp)      mdv
  real(kind=sp)      rho3(nx,ny,nz)
  real(kind=sp)      t3  (nx,ny,nz)
  real(kind=sp)      p3  (nx,ny,nz)
  
  !  Declaration of physical constants
  real(kind=sp)      eps 
  parameter (eps=0.01)
  real(kind=sp)      rd
  parameter (rd=287.)
  real(kind=sp)      tzero
  parameter (tzero=273.15)
  
  !  Auxiliary variables
  integer   i,j,k
  
  !  Calculation
  do i=1,nx
     do j=1,ny
        do k=1,nz
           
           if ((abs(t3(i,j,k)-mdv).gt.eps).and.&
                (abs(p3(i,j,k)-mdv).gt.eps)) then
              
              rho3(i,j,k)=1./rd*p3(i,j,k)/(t3(i,j,k)+tzero)
              
           else
              rho3(i,j,k)=mdv
           endif
           
        enddo
     enddo
  enddo
  
end subroutine calc_rho


!  -----------------------------------------------------------------
!  Calculate potential temperature
!  -----------------------------------------------------------------

subroutine calc_theta (th3,t3,p3,nx,ny,nz,mdv)

  !  Calculate potential temperature <th3> (in K) from temperature <t3> 
  !  (in deg C) and pressure <p3> (in hPa). 
  use kind_parameters,ONLY:&
       sp
    
  implicit none
  
  !  Declaration of subroutine parameters
  integer   nx,ny,nz
  real(kind=sp)      mdv
  real(kind=sp)      th3(nx,ny,nz)
  real(kind=sp)      t3 (nx,ny,nz)
  real(kind=sp)      p3 (nx,ny,nz)
  
  !  Declaration of physical constants
  real(kind=sp)      rdcp,tzero,p0
  parameter (rdcp=0.286)
  parameter (tzero=273.15)
  parameter (p0=100000.)
  real(kind=sp)      eps
  parameter (eps=0.01)
  
  !  Auxiliary variables
  integer  i,j,k
  
  !  Calculation
  
  do i=1,nx
     do j=1,ny
        do k=1,nz
           
           if ((abs(t3(i,j,k)-mdv).gt.eps).and.&
                  (abs(p3(i,j,k)-mdv).gt.eps)) then
              
              th3(i,j,k)=(t3(i,j,k)+tzero)*( (p0/p3(i,j,k))**rdcp )
              
           else
              th3(i,j,k)=mdv
              
              
           endif
           
        enddo
     enddo
  enddo
  
end subroutine calc_theta

!  -----------------------------------------------------------------
!  Calculate stratification
!  -----------------------------------------------------------------

subroutine calc_nsq (nsq3,t3,q3,z3,nx,ny,nz,mdv)
        
  !  Calculate stratification <nsq3> on the model grid. The grid is
  !  specified in the horizontal by <xmin,ymin,dx,dy,nx,ny>. The number
  !  of vertical levels is <nz>. The input field are: temperature <t3>, 
  !  specific humidity <q3>, height <z3> and horizontal grid <x2>, <y2>. 
  use kind_parameters,ONLY:&
        sp
  
  implicit none
  
  !  Declaration of subroutine parameters
  integer   nx,ny,nz
  real(kind=sp)      nsq3(nx,ny,nz)
  real(kind=sp)      t3  (nx,ny,nz)
  real(kind=sp)      q3  (nx,ny,nz)
  real(kind=sp)      z3  (nx,ny,nz)
  real(kind=sp)      x2  (nx,ny)
  real(kind=sp)      y2  (nx,ny)
  real(kind=sp)      mdv
  
  !  Physical and numerical parameters
  real(kind=sp)      scale
  parameter (scale=1.)
  real(kind=sp)      g
  parameter (g=9.80665)
  real(kind=sp)      eps
  parameter (eps=0.01)
  real(kind=sp)      tzero
  parameter (tzero=273.15)
  real(kind=sp)      kappa
  parameter (kappa=0.6078)
  real(kind=sp)      cp
  parameter (cp=1005.)
  real(kind=sp)      zerodiv
  parameter (zerodiv=0.0000000001)
  real(kind=sp)      R
  parameter (R=287.)
  
  !  Auxiliary variables
  real(kind=sp)      tv   (nx,ny,nz)
  real(kind=sp)      dtvdz(nx,ny,nz)
  integer   i,j,k
  real(kind=sp)      scaledz
  
  !  Calculate 3d virtual temperature
  do k=1,nz
     do i=1,nx
        do j=1,ny
           if ((abs(t3(i,j,k)-mdv).gt.eps).and.&
                (abs(q3(i,j,k)-mdv).gt.eps)) &
                then
              tv(i,j,k) = (t3(i,j,k)+tzero)*(1.+kappa*q3(i,j,k))
           else
              tv(i,j,k) = mdv
           endif
        enddo
     enddo
  enddo
  
!  Vertical derivative of virtual temperature
  call  deriv (dtvdz,tv,'z',z3,x2,y2,nx,ny,nz,mdv)
  
  !  Stratification
  do i=1,nx
     do j=1,ny
        do k=1,nz
           if ((abs(dtvdz(i,j,k)-mdv).gt.eps).and.&
                (abs(tv   (i,j,k)-mdv).gt.eps)) &
                then
              nsq3(i,j,k) = g/tv(i,j,k) * (dtvdz(i,j,k) + g/cp)  
           else
              nsq3(i,j,k) = mdv
           endif
        enddo
     enddo
  enddo
  
end subroutine calc_nsq
          
!  -----------------------------------------------------------------
!  Calculate potential vorticity
!  -----------------------------------------------------------------

subroutine calc_pv (pv3,t3,p3,u3,v3,z3,x2,y2,f2,nx,ny,nz,mdv)

!  Calculate potential vorticity <pv3> on the model grid. The grid is
!  specified in the horizontal by <xmin,ymin,dx,dy,nx,ny>. The number
!  of vertical levels is <nz>. The input field are: potential temperature
!  <th3>, horizontal wind <u3> and <v3>, density <rho3>, height <z3>. 
!  The terms including the vertical velocity are neglected in the calculation 
!  of the PV. 
!  use,intrinsic:: ieee_arithmetic,only:IEEE_SELECTED_REAL_KIND
!  use,intrinsic :: ieee_exceptions
!  use ieee_arithmetic
!  use ieee_features
  use kind_parameters,ONLY:&
       sp
  
  implicit none
!  integer,parameter::sp=IEEE_SELECTED_REAL_KIND(6,37)
!  logical should_halt,was_flagged
!  logical underflow_support,gradual,underflow
  !  Declaration of subroutine parameters
  integer   nx,ny,nz
  real(kind=sp)      pv3 (nx,ny,nz)
  real(kind=sp)      t3  (nx,ny,nz)
  real(kind=sp)      p3  (nx,ny,nz)
  real(kind=sp)      u3  (nx,ny,nz)
  real(kind=sp)      v3  (nx,ny,nz)
  real(kind=sp)      z3  (nx,ny,nz)
  real(kind=sp)      x2  (nx,ny)
  real(kind=sp)      y2  (nx,ny)
  real(kind=sp)      f2  (nx,ny)
  
  real(kind=sp)      mdv
  
  !  Physical and numerical parameters
  real(kind=sp)      scale
  parameter (scale=1.E6)
  real(kind=sp)      omega
  parameter (omega=7.292E-5)
  real(kind=sp)      pi180
  parameter (pi180=3.141592654/180.)
  real(kind=sp)      eps
  parameter (eps=1.12E-07)
  
  !  Auxiliary variables
  real(kind=sp)      dtdz(nx,ny,nz)
  real(kind=sp)      dtdx(nx,ny,nz)
  real(kind=sp)      dtdy(nx,ny,nz)
  real(kind=sp)      dudz(nx,ny,nz)
  real(kind=sp)      dvdz(nx,ny,nz)
  real(kind=sp)      dvdx(nx,ny,nz)
  real(kind=sp)      dudy(nx,ny,nz)
  real(kind=sp)      rho3(nx,ny,nz)
  real(kind=sp)      th3 (nx,ny,nz)
  
  integer   i,j,k
  logical :: underflow_support, gradual, underflow


!  call ieee_set_underflow_mode(.false.)
!  call ieee_get_underflow_mode(gradual)

!  call ieee_get_halting_mode(IEEE_UNDERFLOW,should_halt)
!  call ieee_get_flag(IEEE_UNDERFLOW,was_flagged)
!  call ieee_set_halting_mode(IEEE_UNDERFLOW,.FALSE.)
  !  Calculate density and potential temperature
  call calc_rho   (rho3,t3,p3,nx,ny,nz,mdv)
  call calc_theta (th3 ,t3,p3,nx,ny,nz,mdv)
  
  !  Calculate needed derivatives
  call deriv (dudz, u3,'z',z3,x2,y2,nx,ny,nz,mdv)
  
  call deriv (dvdz, v3,'z',z3,x2,y2,nx,ny,nz,mdv)
  
  call deriv (dtdz,th3,'z',z3,x2,y2,nx,ny,nz,mdv)
  
  call deriv (dtdx,th3,'x',z3,x2,y2,nx,ny,nz,mdv)
  
  call deriv (dtdy,th3,'y',z3,x2,y2,nx,ny,nz,mdv)     
  
  call deriv (dvdx, v3,'x',z3,x2,y2,nx,ny,nz,mdv)
  
  call deriv (dudy, u3,'y',z3,x2,y2,nx,ny,nz,mdv)     
  
!  Calculate potential vorticity
  
  do j=1,ny
     do i=1,nx
        do k=1,nz
           
           !           Evaluate PV formula with missing data check
           if ((abs(dtdz(i,j,k)-mdv).gt.eps).and.&
                (abs(dudz(i,j,k)-mdv).gt.eps).and.&
                (abs(dtdy(i,j,k)-mdv).gt.eps).and.&
                (abs(dvdz(i,j,k)-mdv).gt.eps).and.&
                (abs(rho3(i,j,k)-mdv).gt.eps).and. &             
                (abs(dtdx(i,j,k)-mdv).gt.eps).and.&
                (abs(dvdx(i,j,k)-mdv).gt.eps).and.&
                (abs(dudy(i,j,k)-mdv).gt.eps)) then


              pv3(i,j,k)=scale*1./rho3(i,j,k)*&
                   ((dvdx(i,j,k)-dudy(i,j,k)+f2(i,j))*dtdz(i,j,k)&
                   +dudz(i,j,k)*dtdy(i,j,k)*&
                   -dvdz(i,j,k)*dtdx(i,j,k))

!              call ieee_get_flag(ieee_underflow,underflow)

           else
              pv3(i,j,k)=mdv
              
           endif
           
        enddo
     enddo
  enddo
!  call ieee_set_halting_mode(IEEE_UNDERFLOW,should_halt)
!  call ieee_set_flag(IEEE_UNDERFLOW,was_flagged)
end subroutine calc_pv


!  ****************************************************************
!  * SUBROUTINE SECTION: GRID HANDLING                            *
!  ****************************************************************

!  -----------------------------------------------------------------
!  Horizontal and vertical derivatives for 3d fields
!  -----------------------------------------------------------------

subroutine deriv (df,f,direction,z3,x2,y2,nx,ny,nz,mdv)

!  Calculate horizontal and vertical derivatives of the 3d field <f>.
!  The direction of the derivative is specified in <direction> 
!      'x','y'          : Horizontal derivative in x and y direction
!      'p','z','t','m'  : Vertical derivative (pressure, height, theta, model)
!  The 3d field <z3> specifies the isosurfaces along which the horizontal 
!  derivatives are calculated or the levels for the vertical derivatives.
  use ieee_arithmetic
  use ieee_features
  use kind_parameters,ONLY:&
       sp
  implicit none

!  Input and output parameters
  integer    nx,ny,nz
  real(kind=sp)      df (nx,ny,nz)
  real(kind=sp)       f  (nx,ny,nz)
  real(kind=sp)       z3 (nx,ny,nz)
  real(kind=sp)       x2 (nx,ny)
  real(kind=sp)       y2 (nx,ny)
  character  direction
  real(kind=sp)       mdv
  
  !  Numerical and physical parameters
  real(kind=sp)       pi180
  parameter  (pi180=4.D0*DATAN(1.D0)/180.)
  
  real(kind=sp)       zerodiv
  parameter  (zerodiv=0.00000001)
  real(kind=sp)       eps
  parameter  (eps=0.01) 
  
  !  Auxiliary variables
  integer    i,j,k
  real(kind=sp)       vmin,vmax
  real(kind=sp)       scale,lat
  real(kind=sp)       vu,vl,vuvl,vlvu
  integer    o,u,w,e,n,s
  logical :: underflow_support, gradual, underflow
  call ieee_set_underflow_mode(.false.)
  call ieee_get_underflow_mode(gradual)

  !  Vertical derivative
  if ((direction.eq.'z').or.&
         (direction.eq.'th').or.&
         (direction.eq.'p').or.&
         (direction.eq.'m').and.&
         (nz.gt.1)) then
     
     do i=1,nx
        do j=1,ny
           do k=1,nz
              
              o=k+1
              if (o.gt.nz) o=nz
              u=k-1
              if (u.lt.1) u=1                  
              
              if ((abs(f(i,j,o)-mdv).gt.eps).and.&
                   (abs(f(i,j,u)-mdv).gt.eps).and.&
                   (abs(f(i,j,k)-mdv).gt.eps)) then
                 
                 vu = z3(i,j,k)-z3(i,j,o)
                 vl = z3(i,j,u)-z3(i,j,k)
                 vuvl = vu/(vl+zerodiv)
                 vlvu = 1./(vuvl+zerodiv)

                 df(i,j,k) = 1./(vu+vl)&
                      * (vuvl*(f(i,j,u)-f(i,j,k)) &
                      +  vlvu*(f(i,j,k)-f(i,j,o))) 
                 call ieee_get_flag(ieee_underflow,underflow) 
              else
                 df(i,j,k) = mdv
              endif
              
           enddo
        enddo
     enddo
     
!  Horizontal derivative in the y direction: 3d
  elseif (direction.eq.'y') then
     
     do i=1,nx
        do j=1,ny
           do k=1,nz
              
              s=j-1
              if (s.lt.1) s=1
              n=j+1
              if (n.gt.ny) n=ny
              
              if ((abs(f(i,n,k)-mdv).gt.eps).and.&
                   (abs(f(i,j,k)-mdv).gt.eps).and.&
                   (abs(f(i,s,k)-mdv).gt.eps)) then  
                 
                 vu = 1000.*(y2(i,j)-y2(i,n))
                 vl = 1000.*(y2(i,s)-y2(i,j))
                 vuvl = vu/(vl+zerodiv)
                 vlvu = 1./(vuvl+zerodiv)
 
                 df(i,j,k) =  1./(vu+vl)&
                      * (vuvl*(f(i,s,k)-f(i,j,k))&
                      +  vlvu*(f(i,j,k)-f(i,n,k)))
                 
              else
                 df(i,j,k) = mdv
              endif
              
           enddo
        enddo
     enddo
     
!  Horizontal derivative in the x direction: 3d
  elseif (direction.eq.'x') then
     do i=1,nx
        do j=1,ny
           do k=1,nz
              !  what is w ?
              w=i-1
              if (w.lt.1) w=1
              e=i+1
              if (e.gt.nx) e=nx
              
              if ((abs(f(w,j,k)-mdv).gt.eps).and.&
                   (abs(f(i,j,k)-mdv).gt.eps).and.&
                   (abs(f(e,j,k)-mdv).gt.eps)) then  
                 
                 vu = 1000.*(x2(i,j)-x2(e,j))
                 vl = 1000.*(x2(w,j)-x2(i,j))
                 
                 vuvl = vu/(vl+zerodiv)
                 vlvu = 1./(vuvl+zerodiv)
                 
                 
                 df(i,j,k) =  1./(vu+vl)&
                      * (vuvl*(f(w,j,k)-f(i,j,k))&
                      +  vlvu*(f(i,j,k)-f(e,j,k)))
                 
              else
                 df(i,j,k) = mdv
              endif
           enddo
        enddo
     enddo
     
!  Undefined direction for derivative
  else
     
     print*,'Invalid direction of derivative... Stop'
     stop
     
  endif
  
end subroutine deriv

!  -----------------------------------------------------------------
!  Horizontal filter
!  -----------------------------------------------------------------

subroutine filt2d (a,af,nx,ny,fil,misdat,iperx,ipery,ispol,inpol)

!  Apply a conservative diffusion operator onto the 2d field a,
!  with full missing data checking.
!  a     real(kind=sp)   inp  array to be filtered, dimensioned (nx,ny)
!  af    real(kind=sp)   out  filtered array, dimensioned (nx,ny), can be
!                    equivalenced with array a in the calling routine
!  f1    real(kind=sp)        workarray, dimensioned (nx+1,ny)
!  f2    real(kind=sp)        workarray, dimensioned (nx,ny+1)
!  fil   real(kind=sp)   inp  filter-coeff., 0<afil<=1. Maximum filtering with afil=1
!                    corresponds to one application of the linear filter.
!  misdat real(kind=sp)  inp  missing-data value, a(i,j)=misdat indicates that
!                    the corresponding value is not available. The
!                    misdat-checking can be switched off with with misdat=0.
!  iperx int    inp  periodic boundaries in the x-direction (1=yes,0=no)
!  ipery int    inp  periodic boundaries in the y-direction (1=yes,0=no)
!  inpol int    inp  northpole at j=ny  (1=yes,0=no)
!  ispol int    inp  southpole at j=1   (1=yes,0=no)


  use kind_parameters,ONLY:&
       sp

  implicit none

!  argument declaration
  integer     nx,ny
  real(kind=sp)        a(nx,ny),af(nx,ny),fil,misdat
  integer     iperx,ipery,inpol,ispol
  
  !  local variable declaration
  integer     i,j,is
  real(kind=sp)        fh
  real(kind=sp)        f1(nx+1,ny),f2(nx,ny+1)
  
  !  compute constant fh
  fh=0.125*fil
  
  !  compute fluxes in x-direction
  if (misdat.eq.0.) then
     do j=1,ny
        do i=2,nx
           f1(i,j)=a(i-1,j)-a(i,j)
        enddo
     enddo
  else
     do j=1,ny
        do i=2,nx
           if ((a(i,j).eq.misdat).or.(a(i-1,j).eq.misdat)) then
              f1(i,j)=0.
           else
              f1(i,j)=a(i-1,j)-a(i,j)
           endif
        enddo
     enddo
  endif
  if (iperx.eq.1) then
     !    do periodic boundaries in the x-direction
     do j=1,ny
        f1(1,j)=f1(nx,j)
        f1(nx+1,j)=f1(2,j)
     enddo
  else
     !    set boundary-fluxes to zero
     do j=1,ny
        f1(1,j)=0.
        f1(nx+1,j)=0.
     enddo
  endif
  
  !  compute fluxes in y-direction
  if (misdat.eq.0.) then
     do j=2,ny
        do i=1,nx
           f2(i,j)=a(i,j-1)-a(i,j)
        enddo
       enddo
    else
       do j=2,ny
          do i=1,nx
             if ((a(i,j).eq.misdat).or.(a(i,j-1).eq.misdat)) then
                f2(i,j)=0.
             else
                f2(i,j)=a(i,j-1)-a(i,j)
             endif
          enddo
       enddo
    endif
    !  set boundary-fluxes to zero
    do i=1,nx
       f2(i,1)=0.
       f2(i,ny+1)=0.
    enddo
    if (ipery.eq.1) then
       !    do periodic boundaries in the x-direction
       do i=1,nx
          f2(i,1)=f2(i,ny)
          f2(i,ny+1)=f2(i,2)
       enddo
    endif
    if (iperx.eq.1) then
       if (ispol.eq.1) then
          !      do south-pole
          is=(nx-1)/2
          do i=1,nx
             f2(i,1)=-f2(mod(i-1+is,nx)+1,2)
          enddo
       endif
       if (inpol.eq.1) then
          !      do north-pole
          is=(nx-1)/2
          do i=1,nx
             f2(i,ny+1)=-f2(mod(i-1+is,nx)+1,ny)
          enddo
       endif
    endif
    
    !  compute flux-convergence -> filter
    if (misdat.eq.0.) then
       do j=1,ny
          do i=1,nx
             af(i,j)=a(i,j)+fh*(f1(i,j)-f1(i+1,j)+f2(i,j)-f2(i,j+1))
          enddo
       enddo
    else
       do j=1,ny
          do i=1,nx
             if (a(i,j).eq.misdat) then
                af(i,j)=misdat
             else
                af(i,j)=a(i,j)+fh*(f1(i,j)-f1(i+1,j)+f2(i,j)-f2(i,j+1))
             endif
          enddo
       enddo
    endif
  end subroutine filt2d
