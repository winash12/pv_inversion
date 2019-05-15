module cdfwriter
  implicit none
contains

!   ********************************************************************************
!   * NETCDF INPUT/OUTPUT                                                          *
!   ********************************************************************************

!   ---------------------------------------------------------
!   Subroutines to write the netcdf output file
!   ---------------------------------------------------------

  
subroutine writecdf2D(cdfname,cstname,varname,arr,ctime,dx,dy,xmin,ymin,nx,ny,crefile,crevar)
    
    !   Create and write to the netcdf file <cdfname>. The variable
    !   with name <varname> and with time <time> is written. The data
    !   are in the two-dimensional array <arr>. The list <dx,dy,xmin,
    !   ymin,nx,ny> specifies the output grid. The flags <crefile> and
    !   <crevar> determine whether the file and/or the variable should
    !   be created. 1: create / 0: not created
  use kind_parameters,ONLY:&
       sp
  use netcdflibrary
  IMPLICIT NONE
  
  !   Declaration of input parameters
  character*80 cdfname,cstname,varname
  integer      nx,ny
  real(kind=sp),dimension(:,:),allocatable::arr

  real(kind=sp)         dx,dy,xmin,ymin
  real(kind=sp)      ctime
  integer      crefile,crevar
  
  !   Further variables
  !  real(kind=sp)         varmin(4),varmax(4),stag(4)
   real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  integer         ndim,vardim(4)
  integer      cdfid
  real(kind=sp)         mdv
  integer      datar(14),stdate(5)
  integer      i
  real(kind=sp),dimension(:),allocatable::time  
  
  !   Definitions 
  allocate(time(1))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  varmin(1)=xmin
  varmin(2)=ymin
  varmin(3)=1050.
  varmax(1)=xmin+real(nx-1)*dx
  varmax(2)=ymin+real(ny-1)*dy
  varmax(3)=1050.
  ndim=4
  vardim(1)=nx
  vardim(2)=ny
  vardim(3)=1
  stag(1)=0.
  stag(2)=0.
  stag(3)=0.
  mdv=-999.98999
  time(1)=real(ctime)
  
  !   Create the file
  if (crefile.eq.0) then
     call cdfwopn(cdfname,cdfid)
  else if (crefile.eq.1) then
     call crecdf(cdfname,cdfid,varmin,varmax,ndim,cstname)
     
     !      Write the constants file
     datar(1)=vardim(1)
     datar(2)=vardim(2)
     datar(3)=int(1000.*varmax(2))
     datar(4)=int(1000.*varmin(1))
     datar(5)=int(1000.*varmin(2))
     datar(6)=int(1000.*varmax(1))
     datar(7)=int(1000.*dx)
     datar(8)=int(1000.*dy)
     datar(9)=1
     datar(10)=1
     datar(11)=0            ! data version
     datar(12)=0            ! cstfile version
     datar(13)=0            ! longitude of pole
     datar(14)=90000        ! latitude of pole     
     do i=1,5
        stdate(i)=0
     enddo
     call wricstScalar(cstname,datar,0.,0.,0.,0.,stdate)
  endif
  
  !   Write the data to the netcdf file, and close the file
  if (crevar.eq.1) then
     call putdef(cdfid,varname,ndim,mdv,vardim,varmin,varmax,stag)
  endif
  call putdatRank2(cdfid,varname,time,0,arr)
  call clscdf(cdfid)
  
  return
  
  !   Error handling
903 print*,'*** Problem to create netcdf file ***'
  stop
904 print*,'*** Problem to write definition ***'
  stop
905 print*,'*** Problem to write data ***'
  stop
906 print*,'*** Problem to open netcdf file ***'
  stop
  
END subroutine writecdf2D

subroutine writecdf2DRank3(cdfname,cstname,varname,arr,ctime,dx,dy,xmin,ymin,nx,ny,crefile,crevar)
    
    !   Create and write to the netcdf file <cdfname>. The variable
    !   with name <varname> and with time <time> is written. The data
    !   are in the two-dimensional array <arr>. The list <dx,dy,xmin,
    !   ymin,nx,ny> specifies the output grid. The flags <crefile> and
    !   <crevar> determine whether the file and/or the variable should
    !   be created. 1: create / 0: not created
  use netcdflibrary
  IMPLICIT NONE
  
  !   Declaration of input parameters
  character*80 cdfname,cstname,varname
  integer      nx,ny
  real(kind=sp),dimension(:,:,:),allocatable::arr

  real(kind=sp)         dx,dy,xmin,ymin
  real(kind=sp)      ctime
  integer      crefile,crevar
  
  !   Further variables
!  real(kind=sp)         varmin(4),varmax(4),stag(4)
  real(kind=sp),dimension(:),allocatable::varmin,varmax,stag
  integer         ndim,vardim(4)
  integer      cdfid
  real(kind=sp)         mdv
  integer      datar(14),stdate(5)
  integer      i
  real(kind=sp),dimension(:),allocatable::         time


  allocate(time(1))
  allocate(varmin(4))  
  allocate(varmax(4))
  allocate(stag(4))  
  !   Definitions 
  
  varmin(1)=xmin
  varmin(2)=ymin
  varmin(3)=1050.
  varmax(1)=xmin+real(nx-1)*dx
  varmax(2)=ymin+real(ny-1)*dy
  varmax(3)=1050.
  ndim=4
  vardim(1)=nx
  vardim(2)=ny
  vardim(3)=1
  stag(1)=0.
  stag(2)=0.
  stag(3)=0.
  mdv=-999.98999
  time(1)=real(ctime)
  
  !   Create the file
  if (crefile.eq.0) then
     call cdfwopn(cdfname,cdfid)
  else if (crefile.eq.1) then
     call crecdf(cdfname,cdfid,varmin,varmax,ndim,cstname)
     
     !      Write the constants file
     datar(1)=vardim(1)
     datar(2)=vardim(2)
     datar(3)=int(1000.*varmax(2))
     datar(4)=int(1000.*varmin(1))
     datar(5)=int(1000.*varmin(2))
     datar(6)=int(1000.*varmax(1))
     datar(7)=int(1000.*dx)
     datar(8)=int(1000.*dy)
     datar(9)=1
     datar(10)=1
     datar(11)=0            ! data version
     datar(12)=0            ! cstfile version
     datar(13)=0            ! longitude of pole
     datar(14)=90000        ! latitude of pole     
     do i=1,5
        stdate(i)=0
     enddo
     call wricstScalar(cstname,datar,0.,0.,0.,0.,stdate)
  endif
  
  !   Write the data to the netcdf file, and close the file
  if (crevar.eq.1) then
     call putdef(cdfid,varname,ndim,mdv,vardim,varmin,varmax,stag)
  endif
  call putdatRank3(cdfid,varname,time,0,arr)
  call clscdf(cdfid)
  
  return
  
  !   Error handling
903 print*,'*** Problem to create netcdf file ***'
  stop
904 print*,'*** Problem to write definition ***'
  stop
905 print*,'*** Problem to write data ***'
  stop
906 print*,'*** Problem to open netcdf file ***'
  stop
  
END subroutine writecdf2DRank3


  
end module cdfwriter
