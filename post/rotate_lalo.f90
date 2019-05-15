PROGRAM rotate_lalo

!   *******************************************************************
!   * Rotate to a geographical latitude/longitude grid                *
!   * Michael Sprenger / Autumn 2006                                  *
!   *******************************************************************
  use kind_parameters,ONLY:&
       sp
  use netcdflibrary
  implicit none
  
  !   -------------------------------------------------------------------
  !   Declaration of parameters and variables
  !   -------------------------------------------------------------------
  
!   Specification of input parameters
  real(kind=sp)         clon,clat,rotation
  integer      nx,ny
  real(kind=sp)         dx,dy,xmin,ymin
  integer      nfield
  character*80 fieldname(100)
  
  !   Numerical and physical parameters
  real(kind=sp)                 degrad
  parameter            (degrad=0.0174532925)
  real(kind=sp)                 omegaerd
  parameter            (omegaerd=7.292e-5 )
  real(kind=sp)                 eps
  parameter            (eps=0.001)
  
  !   Variables for input Z file : height level
  character*80 in_cfn
  real(kind=sp)         in_varmin(3),in_varmax(3),in_stag(3)
  integer      in_vardim(3)
  real(kind=sp)         in_mdv
  integer      in_ndim
  integer      in_nx,in_ny,in_nz
  real(kind=sp)         in_xmin,in_xmax,in_ymin,in_ymax,in_dx,in_dy
  integer      in_ntimes
  real(kind=sp)         in_aklev(500),in_bklev(500)
  real(kind=sp)         in_aklay(500),in_bklay(500)
  real(kind=sp),dimension(:),allocatable::in_time,out_time
  real(kind=sp)         in_pollon,in_pollat
  integer      in_nvars
  character*80 in_vnam(100)
  integer      in_idate(5)
  real(kind=sp),allocatable, dimension (:,:,:) :: in_field3
  real(kind=sp),allocatable, dimension (:,:,:) :: in_vect1
  real(kind=sp),allocatable, dimension (:,:,:) :: in_vect2
  
  !   Variables for output Z file : height level
  character*80 out_cfn
  real(kind=sp)         out_varmin(3),out_varmax(3),out_stag(3)
  integer      out_vardim(3)
  real(kind=sp)         out_mdv
  integer      out_ndim
  integer      out_nx,out_ny,out_nz
  real(kind=sp)         out_xmin,out_xmax,out_ymin,out_ymax,out_dx,out_dy
  integer      out_ntimes
  real(kind=sp)         out_aklev(500),out_bklev(500)
  real(kind=sp)         out_aklay(500),out_bklay(500)

  real(kind=sp)         out_pollon,out_pollat
  integer      out_idate(5)
  real(kind=sp),allocatable, dimension (:,:,:) :: out_field3
  real(kind=sp),allocatable, dimension (:,:,:) :: out_vect1
  real(kind=sp),allocatable, dimension (:,:,:) :: out_vect2
  real(kind=sp),allocatable, dimension (:,:)   :: out_lat,out_lon
  real(kind=sp),allocatable, dimension (:,:)   :: out_x,out_y
  real(kind=sp),allocatable, dimension (:,:)   :: out_coriol
  
  !   Auxiliary variables
  integer      ifield
  integer      i,j,k
  integer      cdfid,cstid
  character*80 cfn
  integer      stat,ierr,isok
  real(kind=sp)         time
  character*80 varname,cdfname,varname1,varname2
  integer      idate(5),stdate(5),datar(14)
  integer      rotmode
  integer      i1,i2,i3,i4
  integer      isvector
  real(kind=sp)         lat
  character*80 name
  
  !   Externals
  real(kind=sp)         sdis
  external     sdis
  
!   -------------------------------------------------------------------
!   Preparations
!   -------------------------------------------------------------------
  allocate(in_time(1))
  allocate(out_time(1))
  print*,'*********************************************************'
  print*,'* rotate_lalo                                          *'
  print*,'*********************************************************'
  
  !   Read parameter file
  open(10,file='fort.10')
  
  read(10,*) in_cfn
  read(10,*) out_cfn

  read(10,*) name,nx
  if ( trim(name).ne.'GEO_NX'  ) stop
  read(10,*) name,ny
  if ( trim(name).ne.'GEO_NY'  ) stop
  read(10,*) name,dx
  if ( trim(name).ne.'GEO_DX'  ) stop
  read(10,*) name,dy
  if ( trim(name).ne.'GEO_DY'  ) stop
  read(10,*) name,xmin
  if ( trim(name).ne.'GEO_XMIN') stop
  read(10,*) name,ymin
  if ( trim(name).ne.'GEO_YMIN') stop
  read(10,*) name,clon
  if ( trim(name).ne.'CLON'    ) stop
  read(10,*) name,clat
  if ( trim(name).ne.'CLAT'    ) stop
  read(10,*) name,rotation
  if ( trim(name).ne.'CROT'    ) stop
  
  read(10,*) nfield
  do i=1,nfield
     read(10,*) fieldname(i)
  enddo
  close(10)
  
  print*,clon,clat,rotation
  print*,nx,ny,dx,dy,xmin,ymin
  print*,trim(in_cfn),' -> ',trim(out_cfn)
  
  !   Get grid description for input Z file : height level
  call cdfopn(in_cfn,cdfid)
  
  call getcfn(cdfid,cfn)
  
  call cdfopn(cfn,cstid)

  call getvars(cdfid,in_nvars,in_vnam)

  isok=0
  varname=fieldname(1)
  call check_varok (isok,varname,in_vnam,in_nvars)
  if (isok.eq.0) goto 998
  call getdef(cdfid,varname,in_ndim,in_mdv,in_vardim,in_varmin,in_varmax,in_stag)

  in_nx  =in_vardim(1)
  in_ny  =in_vardim(2)
  in_xmin=in_varmin(1)
  in_ymin=in_varmin(2)
  call getlevs(cstid,in_nz,in_aklev,in_bklev,in_aklay,in_bklay)
  call getgrid(cstid,in_dx,in_dy)
  in_xmax=in_xmin+real(in_nx-1)*in_dx
  in_ymax=in_ymin+real(in_ny-1)*in_dy
  call gettimes(cdfid,in_time,in_ntimes)
  call getstart(cstid,in_idate)
  call getpole(cstid,in_pollon,in_pollat)
  call clscdf(cstid)
  call clscdf(cdfid)

!   Set grid description for output file : height level
  out_vardim(1) = nx
  out_vardim(2) = ny
  out_vardim(3) = in_nz
  out_varmin(1) = xmin
  out_varmin(2) = ymin
  out_varmin(3) = in_varmin(3)
  out_varmax(1) = xmin+real(nx-1)*dx
  out_varmax(2) = ymin+real(ny-1)*dy
  out_varmax(3) = in_varmax(3)
  do i=1,in_nz
     out_aklay(i) = in_aklay(i)
     out_bklay(i) = in_bklay(i)
     out_aklev(i) = in_aklev(i)
     out_bklev(i) = in_bklev(i)
  enddo
  out_dx       = dx
  out_dy       = dy
  out_time     = in_time
  out_ntimes   = in_ntimes
  out_ndim     = 4
  out_mdv      = in_mdv
  out_nx       = out_vardim(1)
  out_ny       = out_vardim(2)
  out_nz       = out_vardim(3)
  out_xmin     = out_varmin(1)
  out_ymin     = out_varmin(2)
  out_pollon   = 0.
  out_pollat   = 90.
  do i=1,5
     out_idate(i) = in_idate(i)
  enddo
  
!   Allocate memory for all fields
  allocate(in_field3(in_nx,in_ny,in_nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array in_field3 ***'
  allocate(in_vect1(in_nx,in_ny,in_nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array in_vect1 ***'
  allocate(in_vect2(in_nx,in_ny,in_nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array in_vect2 ***'
  allocate(out_field3(out_nx,out_ny,out_nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array out_field3 ***'
  allocate(out_vect1(out_nx,out_ny,out_nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array out_vect1 ***'
  allocate(out_vect2(out_nx,out_ny,out_nz),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array out_vect2 ***'
  allocate(out_lat(out_nx,out_ny),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array out_lat ***'
  allocate(out_lon(out_nx,out_ny),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array out_lon ***'
  allocate(out_x(out_nx,out_ny),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array out_x ***'
  allocate(out_y(out_nx,out_ny),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array out_y ***'
  allocate(out_coriol(out_nx,out_ny),stat=stat)
  if (stat.ne.0) print*,'*** error allocating array out_coriol ***'
  
  !   Create constants file and output file (if requested) 
  datar(1)=out_nx
  datar(2)=out_ny
  datar(3)=int(1000.*out_varmax(2))
  datar(4)=int(1000.*out_varmin(1))
  datar(5)=int(1000.*out_varmin(2))
  datar(6)=int(1000.*out_varmax(1))
  datar(7)=int(1000.*out_dx)
  datar(8)=int(1000.*out_dy)
  datar(9)=out_nz
  datar(10)=1
  datar(11)=1      
  datar(12)=0      
  datar(13)=int(1000.*out_pollon) 
  datar(14)=int(1000.*out_pollat) 
  
  cfn=trim(out_cfn)//'_cst'
  call wricst(cfn,datar,out_aklev,out_bklev,out_aklay,out_bklay,out_idate)
      
  call crecdf(trim(out_cfn),cdfid,out_varmin,out_varmax,out_ndim,trim(cfn))
  call clscdf(cdfid)
  
!   -----------------------------------------------------------------
!   Loop through all fields - rotate to new grid
  !   -----------------------------------------------------------------
  
  do ifield=1,nfield
     !      Check if scalar or vectorial field (X for scalar, U.V for vector)
     varname=fieldname(ifield)
     i1=1
     i2=1
     i3=0
     i4=0
100  if (varname(i1:i1).eq.' ') then
        i1=i1+1
        goto 100
     endif
     i2=i1
101  if ((varname(i2:i2).ne.' ').and.(varname(i2:i2).ne.'.')) then
        i2=i2+1
        goto 101
     endif
     if (varname(i2:i2).eq.'.') then
        i3=i2+1
102     if (varname(i3:i3).eq.' ') then
           i3=i3+1
           goto 102
        endif
        i4=i3
104     if (varname(i4:i4).ne.' ') then
           i4=i4+1
           goto 104
        endif
     endif
     if ((i3.ne.0).and.(i4.ne.0)) then
        isvector=1
        i2=i2-1
        varname1=varname(i1:i2)
        i4=i4-1         
        varname2=varname(i3:i4)
        print*,'Rotating vector : ',trim(varname1),' /a ',trim(varname2)
     else
        isvector=0
        i2=i2-1
        varname=varname(i1:i2)
        print*,'Rotating scalar : ',trim(varname)
     endif
     
          !      Rotate a scalar field
     if (isvector.eq.0) then
        
        !         Read input field
        call cdfopn(in_cfn,cdfid)
        call getdef(cdfid,varname,in_ndim,in_mdv,in_vardim,in_varmin,in_varmax,in_stag)
        in_nz=in_vardim(3)
        call getdat(cdfid,varname,in_time,0,in_field3)
        call clscdf(cdfid) 
        
        !         Rotate to new coordinate system
        out_nz=in_nz
        call getenvir (clon,clat,rotation,0,in_field3,in_dx,in_dy,in_xmin,&
             in_ymin,in_nx,in_ny,in_nz,in_mdv,out_field3,out_dx,out_dy,&
             out_xmin,out_ymin,out_nx,out_ny,out_nz)
        !         Write output scalar field
        call cdfwopn(trim(out_cfn),cdfid)
        out_vardim(3)=out_nz
        call putdef(cdfid,varname,4,out_mdv,out_vardim,out_varmin,out_varmax,out_stag)         
        call putdat(cdfid,varname,out_time,0,out_field3)
        call clscdf(cdfid)
        
        !      Rotate vector field
     else if (isvector.eq.1) then
        call cdfopn(in_cfn,cdfid)
        call getdef(cdfid,varname1,in_ndim,in_mdv,in_vardim,in_varmin,in_varmax,in_stag)
        in_nz=in_vardim(3)
        call getdat(cdfid,varname1,in_time,0,in_vect1)
        call getdef(cdfid,varname2,in_ndim,in_mdv,in_vardim,in_varmin,in_varmax,in_stag)
        in_nz=in_vardim(3)
        call getdat(cdfid,varname2,in_time,0,in_vect2)
        call clscdf(cdfid) 
        
        !         Get new vector component in x direction
        out_nz=in_nz
        call getenvir (clon,clat,rotation,1,&
             in_vect1,in_dx,in_dy,&
             in_xmin,in_ymin,in_nx,in_ny,in_nz,in_mdv,&
             out_field3,out_dx,out_dy,&
             out_xmin,out_ymin,out_nx,out_ny,out_nz)
        do i=1,out_nx
           do j=1,out_ny
              do k=1,out_nz
                 out_vect1(i,j,k)=out_field3(i,j,k)
              enddo
           enddo
        enddo
        call getenvir (clon,clat,rotation,2,&
             in_vect2,in_dx,in_dy,&
             in_xmin,in_ymin,in_nx,in_ny,in_nz,in_mdv,&
             out_field3,out_dx,out_dy,&
             out_xmin,out_ymin,out_nx,out_ny,out_nz)
        do i=1,out_nx
           do j=1,out_ny
              do k=1,out_nz
                 if ( (abs(out_vect1 (i,j,k)-out_mdv).gt.eps).and.&
                      (abs(out_field3(i,j,k)-out_mdv).gt.eps) ) then
                    out_vect1(i,j,k)=out_vect1(i,j,k)-out_field3(i,j,k)
                 endif
              enddo
           enddo
        enddo
        !         Get new vector component in y direction
        out_nz=in_nz
        call getenvir (clon,clat,rotation,2,&
             in_vect1,in_dx,in_dy,&
             in_xmin,in_ymin,in_nx,in_ny,in_nz,in_mdv,&
             out_field3,out_dx,out_dy,&
             out_xmin,out_ymin,out_nx,out_ny,out_nz)
        do i=1,out_nx
           do j=1,out_ny
              do k=1,out_nz
                 out_vect2(i,j,k)=out_field3(i,j,k)
              enddo
           enddo
        enddo
        call getenvir (clon,clat,rotation,1,&
             in_vect2,in_dx,in_dy,&
             in_xmin,in_ymin,in_nx,in_ny,in_nz,in_mdv,&
             out_field3,out_dx,out_dy,&
             out_xmin,out_ymin,out_nx,out_ny,out_nz)
        do i=1,out_nx
           do j=1,out_ny
              do k=1,out_nz
                 if ( (abs(out_vect2 (i,j,k)-out_mdv).gt.eps).and.&
                      (abs(out_field3(i,j,k)-out_mdv).gt.eps) ) then
                    out_vect2(i,j,k)=out_vect2(i,j,k)+out_field3(i,j,k)
                 endif
           enddo
        enddo
     enddo
     !         Write output vector field
     call cdfwopn(trim(out_cfn),cdfid)
     out_vardim(3)=out_nz
     call putdef(cdfid,varname1,4,out_mdv,out_vardim,out_varmin,out_varmax,out_stag)         
     call putdat(cdfid,varname1,out_time,0,out_vect1)
     call putdef(cdfid,varname2,4,out_mdv,out_vardim,out_varmin,out_varmax,out_stag)         
     call putdat(cdfid,varname2,out_time,0,out_vect2)
     call clscdf(cdfid)
     


        !         Read input vector field
     end if

     
  enddo

!   -----------------------------------------------------------------
!   Write additional fields: latitude, longitude, Coriolis parameter
!   -----------------------------------------------------------------

!   Open the output file
  call cdfwopn(trim(out_cfn),cdfid)

  !   Geographical latitude
  varname='RLAT'
  print*,'Rotating scalar : ',trim(varname)
  do i=1,in_nx
     do j=1,in_ny
        in_field3(i,j,1)=in_ymin+real(j-1)*in_dy
     enddo
  enddo
  call getenvir (clon,clat,rotation,0,&
       in_field3,in_dx,in_dy,&
       in_xmin,in_ymin,in_nx,in_ny,1,in_mdv,&
       out_lat,out_dx,out_dy,&
       out_xmin,out_ymin,out_nx,out_ny,1)      
  out_vardim(3)=1
  call putdef(cdfid,varname,4,out_mdv,out_vardim,out_varmin,out_varmax,out_stag)         
  call putdatRank2(cdfid,varname,out_time,0,out_lat)
  
  !   Geographical longitude
  varname='RLON'
  print*,'Rotating scalar : ',trim(varname)
  do i=1,in_nx
     do j=1,in_ny
        in_field3(i,j,1)=in_xmin+real(i-1)*in_dx
     enddo
  enddo
  call getenvir (clon,clat,rotation,0,&
       in_field3,in_dx,in_dy,&
       in_xmin,in_ymin,in_nx,in_ny,1,in_mdv,&
       out_lon,out_dx,out_dy,&
       out_xmin,out_ymin,out_nx,out_ny,1)      
  out_vardim(3)=1
  call putdef(cdfid,varname,4,out_mdv,out_vardim,out_varmin,out_varmax,out_stag)         
  call putdatRank2(cdfid,varname,out_time,0,out_lon)
  
  
  
  !   Close output file
  call clscdf(cdfid)

  !   -----------------------------------------------------------------
  !   Exception handling and format specs
  !   -----------------------------------------------------------------
  
  stop

 998  print*,'Z: Problems with input rotated grid'
      stop
  
end PROGRAM rotate_lalo


!   ********************************************************************************
!   * SUBROUTINE: ROTATION TO LATITUDE/LONGITUDE COORDINATE SYSTEM                    *
!   ********************************************************************************

!   --------------------------------------------------------------------------------
!   Subroutine to get environment of strcof
!   --------------------------------------------------------------------------------

SUBROUTINE getenvir (clon,clat,rotation,rotmode,inar, rdx,rdy,rxmin,rymin,rnx,rny,rnz,mdv,outar, dx, dy, xmin, ymin, nx, ny, nz)
  
  !   Rotate from a local quasi-cartesian coordiante system into lat/lon coordinate 
  !   system.
  use kind_parameters,ONLY:&
       sp
  use cdfwriter
  implicit none
  
  !   Declaration of input parameters
  integer     rnx,rny,rnz
  real(kind=sp)        rdx,rdy,rxmin,rymin
  real(kind=sp)        inar(rnx,rny,rnz)
  real(kind=sp)        clon,clat,rotation
  real(kind=sp)        mdv
  integer     rotmode
  
  !   Declaration of output parameters
  integer     nx,ny,nz
  real(kind=sp)        dx,dy,xmin,ymin
  real(kind=sp)        outar(nx,ny,nz)
  
  !   Set numerical and physical constants
  real(kind=sp)	  g2r
  parameter   (g2r=0.0174533)
  real(kind=sp)        pi180
  parameter   (pi180=3.14159265359/180.)
  real(kind=sp)        eps
  parameter   (eps=0.0001)
  
  
!   Flag for test mode
  integer      test
  parameter    (test=0)
  character*80 testfile
  parameter    (testfile='ROTGRID')
  
  !   Auxiliary variables 
  real(kind=sp)         pollon,pollat
  integer      i,j,k,l
!  real(kind=sp)         rlon(nx,ny),rlat(nx,ny)
 ! real(kind=sp)         rlon1(nx,ny),rlat1(nx,ny)
  real(kind=sp),dimension(:,:),allocatable::rlon,rlat,rlon1,rlat1,lon,lat
  real(kind=sp),dimension(:,:),allocatable::rotangle1,rotangle2
  real(kind=sp),dimension(:,:,:),allocatable::sinoutar,cosoutar
  real(kind=sp)         rotangle(nx,ny)

  real(kind=sp)         rxmax,rymax
  real(kind=sp)         xind,yind,pind
  real(kind=sp)         outval
  integer      ix,iy
  real(kind=sp)         ax,ay,az,bx,by,bz,zx,zy,zz
  real(kind=sp)         clon1,clat1,clon2,clat2
  real(kind=sp)         rindx,rindy
  integer      indx,indy,indr,indu
  real(kind=sp)         frac0i,frac0j,frac1i,frac1j
  character*80 cdfname,varname,cstname
  real(kind=sp)         vx1,vx2,vy1,vy2,angle,vx2min
  integer      s
  
  !   Externals
  real(kind=sp)     lmtolms,phtophs
  external lmtolms,phtophs

  allocate(rlon1(nx,ny))
  allocate(rlon(nx,ny))
  allocate(rlat1(nx,ny))
  allocate(rlat(nx,ny))
  allocate(lat(nx,ny))
  allocate(lon(nx,ny))
  allocate(rotangle1(nx,ny))
  allocate(rotangle2(nx,ny))
  allocate(sinoutar(nx,ny,nz))
  allocate(cosoutar(nx,ny,nz))
  !   Get geographical coordinates
  do i=1,nx
     do j=1,ny               
        lon(i,j)=xmin+real(i-1)*dx
        lat(i,j)=ymin+real(j-1)*dy
     enddo
  enddo
  
  !   First rotation
  pollon=clon-180.
  if (pollon.lt.-180.) pollon=pollon+360.
  pollat=90.-clat
  do i=1,nx
     do j=1,ny               
        rlon1(i,j)=lmtolms(lat(i,j),lon(i,j),pollat,pollon)
        rlat1(i,j)=phtophs(lat(i,j),lon(i,j),pollat,pollon)            
     enddo
  enddo
  
  !   Second coordinate transformation 
  pollon=-180.
  pollat=90.+rotation
  do i=1,nx
     do j=1,ny
        rlon(i,j)=90.+lmtolms(rlat1(i,j),rlon1(i,j)-90.,pollat,pollon)
        rlat(i,j)=phtophs(rlat1(i,j),rlon1(i,j)-90.,pollat,pollon)            
     enddo
  enddo
  
  !   Get the angle between the rotated and non-rotated coordinate frame
  if ((rotmode.eq.1).or.(rotmode.eq.2)) then
     do i=2,nx-1
        do j=2,ny-1
           
           !            Angle between latitude circles
           vx1=1.
           vy1=0.
           vx2min=(rlon(i+1,j)-rlon(i-1,j))*cos(pi180*rlat(i,j))
           do s=-1,1,2
              vx2=(rlon(i+1,j)-rlon(i-1,j)+real(s)*360.)*cos(pi180*rlat(i,j))
              if (abs(vx2).lt.abs(vx2min)) vx2min=vx2
           enddo
           vx2=vx2min
           vy2=rlat(i+1,j)-rlat(i-1,j)
           call getangle(vx1,vy1,vx2,vy2,angle)           
           rotangle1(i,j)=angle
           
           !            Angle between longitude circles
           vx1=0.
           vy1=1.
           vx2min=(rlon(i,j+1)-rlon(i,j-1))*cos(pi180*rlat(i,j))
           do s=-1,1,2
              vx2=(rlon(i+1,j)-rlon(i-1,j)+real(s)*360.)*cos(pi180*rlat(i,j))
              if (abs(vx2).lt.abs(vx2min)) vx2min=vx2
           enddo
           vx2=vx2min
           vy2=rlat(i,j+1)-rlat(i,j-1)
           call getangle(vx1,vy1,vx2,vy2,angle)           
           rotangle2(i,j)=angle
           
        enddo
     enddo
     
     !      Set the angle at the boundaries
     do i=1,nx
        rotangle1(i,ny)=2.0*rotangle1(i,ny-1)-rotangle1(i,ny-2)
        rotangle1(i,1) =2.0*rotangle1(i,2)-rotangle1(i,3)
        rotangle2(i,ny)=2.0*rotangle2(i,ny-1)-rotangle2(i,ny-2)
        rotangle2(i,1) =2.0*rotangle2(i,2)-rotangle2(i,3)
     enddo
     do j=1,ny
        rotangle1(nx,j)=2.0*rotangle1(nx-1,j)-rotangle1(nx-2,j)
        rotangle1(1,j) =2.0*rotangle1(2,j)-rotangle1(3,j)
        rotangle2(nx,j)=2.0*rotangle2(nx-1,j)-rotangle2(nx-2,j)
        rotangle2(1,j) =2.0*rotangle2(2,j)-rotangle2(3,j)
     enddo
     
     !      Set the final rotation angle
     do i=1,nx
        do j=1,ny
           rotangle(i,j)=0.5*(rotangle1(i,j)+rotangle2(i,j))
        enddo
     enddo
     
  endif
  
  !   Bring longitude into domain of geographical grid (shift by 360 deg)
  do i=1,nx
     do j=1,ny
100     if (rlon(i,j).lt.rxmin) then
           rlon(i,j)=rlon(i,j)+360.
           goto 100
        endif
102     if (rlon(i,j).gt.(rxmin+real(rnx-1)*rdx)) then
           rlon(i,j)=rlon(i,j)-360.
           goto 102
        endif
     enddo
  enddo
  
!   Rotate the scalar or the vector component
  do i=1,nx
     do j=1,ny
        do k=1,nz   
           
           !            Get index
           rindx=(rlon(i,j)-rxmin)/rdx+1.
           rindy=(rlat(i,j)-rymin)/rdy+1.
           indx=int(rindx)
           indy=int(rindy)
           if ((indx.lt.1).or.(indx.gt.rnx).or.(indy.lt.1).or.(indy.gt.rny)) then
              outar(i,j,k)=mdv
              goto 103
           endif
           
           !            Get inidices of left and upper neighbours
           indr=indx+1
           if (indr.gt.rnx) indr=1
           indu=indy+1
           if (indu.gt.rny) indu=ny
           
           !            Do linear interpolation
           if ( ( abs(inar(indx ,indy, k)-mdv).gt.eps ).and. &
                ( abs(inar(indx ,indu, k)-mdv).gt.eps ).and.&
                ( abs(inar(indr ,indy, k)-mdv).gt.eps ).and.&
                ( abs(inar(indr ,indu, k)-mdv).gt.eps ) ) then
              frac0i=rindx-float(indx)
              frac0j=rindy-float(indy)
              frac1i=1.-frac0i
              frac1j=1.-frac0j
              outval = inar(indx ,indy, k ) * frac1i * frac1j&
                   + inar(indx ,indu, k ) * frac1i * frac0j&
                   + inar(indr ,indy, k ) * frac0i * frac1j&
                   + inar(indr ,indu, k ) * frac0i * frac0j
           else
              outval=mdv
           endif
               
           !            Update output array
           outar(i,j,k)=outval
           
           !            Next
103        continue
           
        enddo
     enddo
  enddo
  
  !   Get components for tests
  if (test.eq.1) then
     do i=1,nx
        do j=1,ny
           do k=1,nz
              cosoutar(i,j,k)=outar(i,j,k)*cos(pi180*rotangle(i,j))
              sinoutar(i,j,k)=-outar(i,j,k)*sin(pi180*rotangle(i,j))
           enddo
        enddo
     enddo
  endif
  
  !   Get correct component of rotated field
  do i=1,nx
     do j=1,ny
        do k=1,nz
           if ( abs(outar(i,j,k)-mdv).gt.eps ) then
              if (rotmode.eq.1) then
                 outar(i,j,k)= outar(i,j,k)*cos(pi180*rotangle(i,j))
              else if (rotmode.eq.2) then
                 outar(i,j,k)=-outar(i,j,k)*sin(pi180*rotangle(i,j))
              else if (rotmode.eq.0) then
                 outar(i,j,k)=outar(i,j,k)
              endif
           endif
        enddo
     enddo
  enddo
  
  !   Test mode: Write grids to cdf file
  if (test.eq.1) then
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='RLON1'
     call  writecdf2D(cdfname,cstname,varname,rlon1,0.,rdx,rdy,rxmin,rymin,rnx,rny,1,1)
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='RLON'
     call  writecdf2D(cdfname,cstname,varname,rlon,0.,rdx,rdy,rxmin,rymin,rnx,rny,0,1)
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='LON'
     call  writecdf2D(cdfname,cstname,varname,lon,0.,rdx,rdy,rxmin,rymin,rnx,rny,0,1)
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='RLAT1'
     call  writecdf2D(cdfname,cstname,varname,rlat1,0.,rdx,rdy,rxmin,rymin,rnx,rny,0,1)
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='RLAT'
     call  writecdf2D(cdfname,cstname,varname,rlat,0.,rdx,rdy,rxmin,rymin,rnx,rny,0,1)
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='LAT'
     call  writecdf2D(cdfname,cstname,varname,lat,0.,rdx,rdy,rxmin,rymin,rnx,rny,0,1)
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='ANGLE1'
     call  writecdf2D(cdfname,cstname,varname,rotangle1,0.,rdx,rdy,rxmin,rymin,rnx,rny,0,1)
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='ANGLE2'
     call  writecdf2D(cdfname,cstname,varname,rotangle2,0.,rdx,rdy,rxmin,rymin,rnx,rny,0,1)
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='U'
     call  writecdf2DRank3(cdfname,cstname,varname,cosoutar,0.,rdx,rdy,rxmin,rymin,rnx,rny,0,1)
     cdfname=testfile
     cstname=trim(testfile)//'_cst'
     varname='V'
     call  writecdf2DRank3(cdfname,cstname,varname,sinoutar,0.,rdx,rdy,rxmin,rymin,rnx,rny,0,1)
     
  endif

  deallocate(rlon1)
  deallocate(rlon)
  deallocate(rlat1)
  deallocate(rlat)
  deallocate(lat)
  deallocate(lon)
  deallocate(rotangle1)
  deallocate(rotangle2)
  deallocate(sinoutar)
  deallocate(cosoutar)

END SUBROUTINE getenvir


!   --------------------------------------------------------------------------------
!   Auxiliary routines: angle between two vectors
!   --------------------------------------------------------------------------------

SUBROUTINE getangle (vx1,vy1,vx2,vy2,angle)
  
  !   Given two vectors <vx1,vy1> and <vx2,vy2>, determine the angle (in deg)
  !   between the two vectors.

  use kind_parameters,ONLY:&
       sp
  
  implicit none
  
  !   Declaration of subroutine parameters
  real(kind=sp) vx1,vy1
  real(kind=sp) vx2,vy2
  real(kind=sp) angle
  
  !   Auxiliary variables and parameters
  real(kind=sp) len1,len2,len3
  real(kind=sp) val1,val2,val3
  real(kind=sp) pi
  parameter (pi=3.14159265359)
  
  len1=sqrt(vx1*vx1+vy1*vy1)
  len2=sqrt(vx2*vx2+vy2*vy2)
  
  if ((len1.gt.0.).and.(len2.gt.0.)) then
     vx1=vx1/len1
     vy1=vy1/len1
     vx2=vx2/len2
     vy2=vy2/len2
     
     val1=vx1*vx2+vy1*vy2
     val2=-vy1*vx2+vx1*vy2
     
     len3=sqrt(val1*val1+val2*val2)
     
     if ( (val1.ge.0.).and.(val2.ge.0.) ) then
        val3=acos(val1/len3)
     else if ( (val1.lt.0.).and.(val2.ge.0.) ) then
        val3=pi-acos(abs(val1)/len3)
     else if ( (val1.ge.0.).and.(val2.le.0.) ) then
        val3=-acos(val1/len3)
     else if ( (val1.lt.0.).and.(val2.le.0.) ) then
        val3=-pi+acos(abs(val1)/len3)
     endif
  else
     val3=0.
  endif
  
  angle=180./pi*val3
  
END SUBROUTINE getangle

!   --------------------------------------------------------------------------------
!   Transformation routine: LMSTOLM and PHSTOPH from library gm2em
!   --------------------------------------------------------------------------------

REAL(KIND=SP) FUNCTION LMTOLMS (PHI, LAM, POLPHI, POLLAM)
  !
  !%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
  !
  !** LMTOLMS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM
  !**                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
  !**                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
  !**                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !   AUFRUF   :   LAM = LMTOLMS (PHI, LAM, POLPHI, POLLAM)
  !   ENTRIES  :   KEINE
  !   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN LAENGE LAM AUF
  !                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
  !                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
  !                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !   VERSIONS-
  !   DATUM    :   03.05.90
  !
  !   EXTERNALS:   KEINE
  !   EINGABE-
  !   PARAMETER:   PHI    REAL(KIND=SP) BREITE DES PUNKTES IM GEOGR. SYSTEM
  !                LAM    REAL(KIND=SP) LAENGE DES PUNKTES IM GEOGR. SYSTEM
  !                POLPHI REAL(KIND=SP) GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
  !                POLLAM REAL(KIND=SP) GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
  !   AUSGABE-
  !   PARAMETER:   WAHRE GEOGRAPHISCHE LAENGE ALS WERT DER FUNKTION
  !                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
  !
  !   COMMON-
  !   BLOECKE  :   KEINE
  !
  !   FEHLERBE-
  !   HANDLUNG :   KEINE
  !   VERFASSER:   G. DE MORSIER
  use kind_parameters,ONLY:&
       sp
  
  REAL(KIND=SP)        LAM,PHI,POLPHI,POLLAM
  
  DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
  
  ZSINPOL = SIN(ZPIR18*POLPHI)
  ZCOSPOL = COS(ZPIR18*POLPHI)
  ZLAMPOL =     ZPIR18*POLLAM
  ZPHI    =     ZPIR18*PHI
  ZLAM    = LAM
  IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
  ZLAM    = ZPIR18*ZLAM
  
  ZARG1   = - SIN(ZLAM-ZLAMPOL)*COS(ZPHI)
  ZARG2   = - ZSINPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL)+ZCOSPOL*SIN(ZPHI)
  IF (ABS(ZARG2).LT.1.E-30) THEN
     IF (ABS(ZARG1).LT.1.E-30) THEN
        LMTOLMS =   0.0
     ELSEIF (ZARG1.GT.0.) THEN
        LMTOLMS =  90.0
     ELSE
        LMTOLMS = -90.0
     ENDIF
  ELSE
     LMTOLMS = ZRPI18*ATAN2(ZARG1,ZARG2)
  ENDIF
  
  RETURN
END FUNCTION LMTOLMS


REAL(KIND=SP) FUNCTION PHTOPHS (PHI, LAM, POLPHI, POLLAM)
  !
  !%Z% Modul %M%, V%I% vom %G%, extrahiert am %H%
  !
  !** PHTOPHS  -   FC:UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI
  !**                 AUF EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS)
  !**                 IM ROTIERTEN SYSTEM. DER NORDPOL DES SYSTEMS HAT
  !**                 DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !   AUFRUF   :   PHI = PHTOPHS (PHI, LAM, POLPHI, POLLAM)
  !   ENTRIES  :   KEINE
  !   ZWECK    :   UMRECHNUNG DER WAHREN GEOGRAPHISCHEN BREITE PHI AUF
  !                EINEM PUNKT MIT DEN KOORDINATEN (PHIS, LAMS) IM
  !                ROTIERTEN SYSTEM. DER NORDPOL DIESES SYSTEMS HAT
  !                DIE WAHREN KOORDINATEN (POLPHI, POLLAM)
  !   VERSIONS-
  !   DATUM    :   03.05.90
  !
  !   EXTERNALS:   KEINE
  !   EINGABE-
  !   PARAMETER:   PHI    REAL(KIND=SP) BREITE DES PUNKTES IM GEOGR. SYSTEM
  !                LAM    REAL(KIND=SP) LAENGE DES PUNKTES IM GEOGR. SYSTEM
  !                POLPHI REAL(KIND=SP) GEOGR.BREITE DES N-POLS DES ROT. SYSTEMS
  !                POLLAM REAL(KIND=SP) GEOGR.LAENGE DES N-POLS DES ROT. SYSTEMS
  !   AUSGABE-
  !   PARAMETER:   ROTIERTE BREITE PHIS ALS WERT DER FUNKTION
  !                ALLE WINKEL IN GRAD (NORDEN>0, OSTEN>0)
  !
  !   COMMON-
  !   BLOECKE  :   KEINE
  !
  !   FEHLERBE-
  !   HANDLUNG :   KEINE
  !   VERFASSER:   G. DE MORSIER
  use kind_parameters,ONLY:&
       sp
  
  REAL(KIND=SP)        LAM,PHI,POLPHI,POLLAM
  
  DATA        ZRPI18 , ZPIR18  / 57.2957795 , 0.0174532925 /
  
  ZSINPOL = SIN(ZPIR18*POLPHI)
  ZCOSPOL = COS(ZPIR18*POLPHI)
  ZLAMPOL = ZPIR18*POLLAM
  ZPHI    = ZPIR18*PHI
  ZLAM    = LAM
  IF(ZLAM.GT.180.0) ZLAM = ZLAM - 360.0
  ZLAM    = ZPIR18*ZLAM
  ZARG    = ZCOSPOL*COS(ZPHI)*COS(ZLAM-ZLAMPOL) + ZSINPOL*SIN(ZPHI)
  
  PHTOPHS = ZRPI18*ASIN(ZARG)
  
  RETURN
END FUNCTION PHTOPHS


!   ---------------------------------------------------------
!   Spherical distance between two lat/lon points
!   ---------------------------------------------------------

real(kind=sp) function sdis(xp,yp,xq,yq)
  !
  !   calculates spherical distance (in km) between two points given
  !   by their spherical coordinates (xp,yp) and (xq,yq), respectively.
  !

  use kind_parameters,ONLY:&
       sp

  real(kind=sp)      re,pi,dconv
  parameter (re=6370.)
  real(kind=sp)      xp,yp,xq,yq,arg
  parameter (pi=3.14159265)
  dconv = pi/180.
  arg=sin(yp*dconv)*sin(yq*dconv)+cos(yp*dconv)*cos(yq*dconv)*cos((xp-xq)*dconv)
  if (arg.lt.-1.) arg=-1.
  if (arg.gt.1.) arg=1.
  sdis=re*acos(arg)
  
end function sdis



!   ----------------------------------------------------------------
!   Check whether variable is found on netcdf file
!   ----------------------------------------------------------------

subroutine check_varok (isok,varname,varlist,nvars)
  
  !   Check whether the variable <varname> is in the list <varlist(nvars)>.
  !   If this is the case, <isok> is incremented by 1. Otherwise <isok>
  !   keeps its value.
  
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
