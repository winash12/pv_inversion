PROGRAM add2p

  !   *********************************************************************
  !   * Add output from the PV inversion to the original P file           *
  !   * Zonal wind (U), meridional wind (V) and temperature (T) are       *
  !   * modified. All other fields are not affected by the programme      *
  !   *                                                                   *
  !   * Michael Sprenger / Summer, Autumn 2006
  !   *********************************************************************
  use kind_parameters,ONLY:&
       sp
  use netcdflibrary
  
  implicit none
  
  !   ---------------------------------------------------------------------
  !   Declaration of parameters and variables
  !   ---------------------------------------------------------------------
  
  !   Input parameters
  character*80 ml_fn
  character*80 zm_fn
  character*80 or_fn
  real(kind=sp)         transxy

!   Variables for input P file : model level
!  real(kind=sp)         ml_varmin(4),ml_varmax(4),ml_stag(4)
  real(kind=sp),allocatable, dimension (:)   ::ml_varmin,ml_varmax,ml_stag
  integer      ml_vardim(4)
  real(kind=sp)         ml_mdv
  integer      ml_ndim
  integer      ml_nx,ml_ny,ml_nz
  real(kind=sp)         ml_xmin,ml_xmax,ml_ymin,ml_ymax,ml_dx,ml_dy
  integer      ml_ntimes
  real(kind=sp)         ml_aklev(500),ml_bklev(500)
  real(kind=sp)         ml_aklay(500),ml_bklay(500)
  real,dimension(:),allocatable::ml_time
  real(kind=sp)         ml_pollon,ml_pollat
  integer      ml_idate(5)
  integer      ml_nvars
  character*80 ml_vnam(80)
  real(kind=sp),allocatable, dimension (:,:)   :: ml_ps,ml_oro
  real(kind=sp),allocatable, dimension (:,:,:) :: ml_t3,ml_u3,ml_v3,ml_z3
  real(kind=sp),allocatable, dimension (:,:,:) :: ml_p3,ml_tv3,ml_q3,ml_p0
  integer,allocatable, dimension (:,:,:) :: flag_ml

!   Variables for GEO.MOD
  !  real(kind=sp)         zm_varmin(4),zm_varmax(4),zm_stag(4)
  real(kind=sp),allocatable, dimension (:)   ::zm_varmin,zm_varmax,zm_stag
  integer      zm_vardim(4)
  real(kind=sp)         zm_mdv
  integer      zm_ndim
  integer      zm_nx,zm_ny,zm_nz
  real(kind=sp)         zm_xmin,zm_xmax,zm_ymin,zm_ymax,zm_dx,zm_dy
  integer      zm_ntimes
  real(kind=sp)         zm_aklev(500),zm_bklev(500)
  real(kind=sp)         zm_aklay(500),zm_bklay(500)
  real(kind=sp),allocatable,dimension(:)::zm_time
  real(kind=sp)         zm_pollon,zm_pollat
  integer      zm_idate(5)
  integer      zm_nvars
  character*80 zm_vnam(80)
  real(kind=sp),allocatable, dimension (:,:,:) :: zm_u3,zm_v3,zm_t3,zm_z3
  real(kind=sp),allocatable, dimension (:,:,:) :: zm_p3
  real,allocatable, dimension (:,:)   :: zm_rlat,zm_rlon
  integer,allocatable, dimension (:,:,:) :: flag_zm
  
!   Variables for GEO.ORG
  real(kind=sp),allocatable, dimension (:)   ::or_varmin,or_varmax,or_stag  
  integer      or_vardim(4)
  real(kind=sp)         or_mdv
  integer      or_ndim
  integer      or_nx,or_ny,or_nz
  real(kind=sp)         or_xmin,or_xmax,or_ymin,or_ymax,or_dx,or_dy
  integer      or_ntimes
  real(kind=sp),allocatable,dimension(:)::or_time
  real(kind=sp)         or_pollon,or_pollat
  integer      or_idate(5)
  integer      or_nvars
  character*80 or_vnam(80)
  real(kind=sp),allocatable, dimension (:,:,:) :: or_p3
  real(kind=sp),allocatable, dimension (:,:,:) :: or_u3,or_v3,or_t3

!   Anomalies
  real(kind=sp),allocatable, dimension (:,:,:) :: an_p3,an_u3,an_v3,an_t3
  
!   Array with the weighting function
  real(kind=sp),allocatable, dimension (:,:)   :: weight,dist
  
  !   Flag for test mode
  integer      test
  parameter    (test=1)
  
  !   Flag for Poisson filling
  integer      poisson
  parameter    (poisson=1)
  
  !   Physical and numerical parameters
  real(kind=sp)         eps                    ! Numerical epsilon
  parameter    (eps=0.001)
  real(kind=sp)         dpextra                ! Allowed range for extrapolation (% of pressure)
  parameter    (dpextra=0.1)
  real(kind=sp)         g                      ! Earth's gravity
  parameter    (g=9.80665)
  real(kind=sp)         tzero                  ! Conversion Celius/Kelvin
  parameter    (tzero=273.15)
  real(kind=sp)         kappa                  ! Kappa (Definition of potential temperature)
  parameter    (kappa=0.6078)
  real(kind=sp)         zerodiv                ! Zero division handling
  parameter    (zerodiv=0.0000000001)
  real(kind=sp)         dpmin                  ! Pressure step for hydrostatic equation
  parameter    (dpmin=10.)
  real(kind=sp)         rdg                    ! Ratio gas constant / Earth's gravity
  parameter    (rdg=29.271)
  
  !   Auxiliray variables
  integer      i,j,k,l
  character*80 varname
  integer      isok
  integer      stat
  integer      cdfid,ierr,cstid
  character*80 cfn
  real(kind=sp)         p1(1000),z1(1000),f1(1000),tv1(1000)
  real(kind=sp)         spline_f1(1000),spline_tv1(1000)
  integer      n1
  real(kind=sp)         hh
  real(kind=sp)         pmin,pmax
  real(kind=sp)         boundlat(100000),boundlon(1000000)
  integer      nbound
  integer      il,ir,ju,jd,im,jm
  real(kind=sp)         lat,lon
  real(kind=sp)         distmin,distpos
  character*80 name
  real(kind=sp)         pu,po,zu,zo,p,z,dp,p0,tvu,tvo,ff
  integer      lmin,n
  real(kind=sp)         mean,nmean
  real(kind=sp)         psfc
  integer      count0,count1
  
  !   Externals
  real(kind=sp)         sdis,kink
  external     sdis,kink

  allocate(ml_time(1))
  allocate(zm_time(1))
  allocate(or_time(1))
  allocate(zm_varmin(4))
  allocate(zm_varmax(4))
  allocate(zm_stag(4))

  allocate(or_varmin(4))
  allocate(or_varmax(4))
  allocate(or_stag(4))
  allocate(ml_varmin(4))
  allocate(ml_varmax(4))
  allocate(ml_stag(4))
  !   -----------------------------------------------------------------
!   Read input fields
!   -----------------------------------------------------------------
  
  print*,'*********************************************************'
  print*,'* add2p                                                 *'
  print*,'*********************************************************'
  print*
  
  !   Read in the parameter file
  open(10,file='fort.10')
  read(10,*) ml_fn
  read(10,*) zm_fn
  read(10,*) or_fn
  read(10,*) name,transxy
  if ( trim(name).ne.'TRANS_XY'  ) stop
  read(10,*) name,psfc
  if ( trim(name).ne.'PS_CHANGE' ) stop
  close(10)
  print*,trim(ml_fn)
  print*,trim(zm_fn)
  print*,trim(or_fn)
  print*
  
  !   Get grid description for P file : model level
  call cdfopn(ml_fn,cdfid)

  call getvars(cdfid,ml_nvars,ml_vnam)

  call getcfn(cdfid,cfn)

  call cdfopn(cfn,cstid)

  isok=0
  varname='T'
  call check_varok(isok,varname,ml_vnam,ml_nvars)
  if (isok.ne.1) goto 998
  call getdef(cdfid,varname,ml_ndim,ml_mdv,ml_vardim,ml_varmin,ml_varmax,ml_stag)

  ml_nx  =ml_vardim(1)
  ml_ny  =ml_vardim(2)
  ml_nz  =ml_vardim(3)
  ml_xmin=ml_varmin(1)
  ml_ymin=ml_varmin(2)
  call getlevs(cstid,ml_nz,ml_aklev,ml_bklev,ml_aklay,ml_bklay)

  call getgrid(cstid,ml_dx,ml_dy)
  ml_xmax=ml_xmin+real(ml_nx-1)*ml_dx
  ml_ymax=ml_ymin+real(ml_ny-1)*ml_dy

  call gettimes(cdfid,ml_time,ml_ntimes)

  call getstart(cstid,ml_idate)

  call getpole(cstid,ml_pollon,ml_pollat)

  call clscdf(cstid)

  call clscdf(cdfid)

  
!   Get grid description for MOD file
  call cdfopn(zm_fn,cdfid)

  call getvars(cdfid,zm_nvars,zm_vnam)

  call getcfn(cdfid,cfn)

  call cdfopn(cfn,cstid)

  isok=0
  varname='T'
  call check_varok(isok,varname,zm_vnam,zm_nvars)
  if (isok.ne.1) goto 997
  call getdef(cdfid,varname,zm_ndim,zm_mdv,zm_vardim,zm_varmin,zm_varmax,zm_stag)

  zm_nx  =zm_vardim(1)
  zm_ny  =zm_vardim(2)
  zm_nz  =zm_vardim(3)
  zm_xmin=zm_varmin(1)
  zm_ymin=zm_varmin(2)
  call getlevs(cstid,zm_nz,zm_aklev,zm_bklev,zm_aklay,zm_bklay)

  call getgrid(cstid,zm_dx,zm_dy)

  zm_xmax=zm_xmin+real(zm_nx-1)*zm_dx
  zm_ymax=zm_ymin+real(zm_ny-1)*zm_dy
  call gettimes(cdfid,zm_time,zm_ntimes)

  call getstart(cstid,zm_idate)

  call getpole(cstid,zm_pollon,zm_pollat)

  call clscdf(cstid)

  call clscdf(cdfid)


  !   Get grid description for ORG file
  call cdfopn(or_fn,cdfid)

  call getvars(cdfid,or_nvars,or_vnam)

  call getcfn(cdfid,cfn)

  call cdfopn(cfn,cstid)

  isok=0
  varname='P'
  call check_varok(isok,varname,or_vnam,or_nvars)
  if (isok.ne.1) goto 997
  call getdef(cdfid,varname,or_ndim,or_mdv,or_vardim,or_varmin,or_varmax,or_stag)

  or_nx  =or_vardim(1)
  or_ny  =or_vardim(2)
  or_nz  =or_vardim(3)
  or_xmin=or_varmin(1)
  or_ymin=or_varmin(2)
  call getgrid(cstid,or_dx,or_dy)

  or_xmax=or_xmin+real(or_nx-1)*or_dx
  or_ymax=or_ymin+real(or_ny-1)*or_dy
  call gettimes(cdfid,or_time,or_ntimes)

  call getstart(cstid,or_idate)

  call getpole(cstid,or_pollon,or_pollat)

  call clscdf(cstid)

  call clscdf(cdfid)

  
  !   Consistency check for the grids
  if ( (abs(ml_xmin-zm_xmin).gt.eps).or.&
          (abs(ml_xmax-zm_xmax).gt.eps).or.&
          (abs(ml_ymin-zm_ymin).gt.eps).or.&
          (abs(ml_ymax-zm_ymax).gt.eps).or.&
          (abs(ml_dx  -zm_dx  ).gt.eps).or.&
          (abs(ml_dy  -zm_dy  ).gt.eps).or.&
          (abs(ml_time(1)-zm_time(1)).gt.eps).or.&
          (abs(ml_xmin-or_xmin).gt.eps).or.&
          (abs(ml_xmax-or_xmax).gt.eps).or.&
          (abs(ml_ymin-or_ymin).gt.eps).or.&
          (abs(ml_ymax-or_ymax).gt.eps).or.&
          (abs(ml_dx  -or_dx  ).gt.eps).or.&
          (abs(ml_dy  -or_dy  ).gt.eps) )&
          then
     print*,'P, GEO.MOD, GEO.ORG  grids not consistent... Stop'
     print*,ml_xmin,zm_xmin,or_xmin
     print*,ml_ymin,zm_ymin,or_ymin
     print*,ml_dx,  zm_dx  ,or_dx
     print*,ml_dy,  zm_dy  ,or_dy
     print*,ml_time,zm_time,or_time
     stop
  endif

!   Allocate memory for input and output P file
  allocate(ml_p3  (ml_nx,ml_ny,ml_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_p3'
  allocate(ml_u3  (ml_nx,ml_ny,ml_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_u3'
  allocate(ml_v3  (ml_nx,ml_ny,ml_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_v3'
  allocate(ml_t3  (ml_nx,ml_ny,ml_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_t3'
  allocate(ml_tv3 (ml_nx,ml_ny,ml_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_tv3'
  allocate(ml_z3  (ml_nx,ml_ny,ml_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_z3'
  allocate(ml_q3  (ml_nx,ml_ny,ml_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_q3'
  allocate(ml_ps  (ml_nx,ml_ny      ),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_ps'
  allocate(ml_oro (ml_nx,ml_ny      ),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_oro'
  allocate(ml_p0  (ml_nx,ml_ny,ml_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating ml_p0'
  
  !   Allocate memory for input GEO.MOD file
  allocate(zm_p3(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_p3'
  allocate(zm_u3(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_u3'
  allocate(zm_v3(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_v3'
  allocate(zm_t3(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_t3'
  allocate(zm_z3(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_z3'
  allocate(zm_rlat(zm_nx,zm_ny),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_rlat'
  allocate(zm_rlon(zm_nx,zm_ny),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_rlon'
  
  !   Allocate memory for input GEO.ORG file
  allocate(or_p3(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_p3'
  allocate(or_t3(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_t3'
  allocate(or_u3(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_u3'
  allocate(or_v3(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating zm_v3'
  
  !   Allocate memory for weighting function, flag and anomalies
  allocate(weight(zm_nx,zm_ny),stat=stat)
  if (stat.ne.0) print*,'error allocating weight'
  allocate(dist  (zm_nx,zm_ny),stat=stat)
  if (stat.ne.0) print*,'error allocating dist'
  allocate(flag_zm(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating flag'
  allocate(flag_ml(zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating flag'
  allocate(an_p3 (zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating an_p3'
  allocate(an_t3 (zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating an_t3'
  allocate(an_u3 (zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating an_u3'
  allocate(an_v3 (zm_nx,zm_ny,zm_nz),stat=stat)
  if (stat.ne.0) print*,'error allocating an_v3'
  
  !   Read variables from P file : model levels
  call cdfopn(ml_fn,cdfid)

  varname='T'
  call getdat(cdfid,varname,ml_time,0,ml_t3)

  print*,'R T    ',trim(ml_fn)
  varname='U'
  call getdat(cdfid,varname,ml_time,0,ml_u3)

  print*,'R U    ',trim(ml_fn)
  varname='V'
  call getdat(cdfid,varname,ml_time,0,ml_v3)

  print*,'R V    ',trim(ml_fn)
  varname='Z'
  call getdat(cdfid,varname,ml_time,0,ml_z3)

  print*,'R Z    ',trim(ml_fn)
  varname='Q'
  call getdat(cdfid,varname,ml_time,0,ml_q3)

  print*,'R Q    ',trim(ml_fn)
  varname='PS'
  call getdatRank2(cdfid,varname,ml_time,0,ml_ps)

  print*,'R PS   ',trim(ml_fn)
  varname='ORO'
  call getdatRank2(cdfid,varname,ml_time,0,ml_oro)

  print*,'R ORO  ',trim(ml_fn)
  call clscdf(cdfid)

  
  !   Read variables from GEO.MOD
  call cdfopn(zm_fn,cdfid)

  
  varname='T'
  call getdat(cdfid,varname,zm_time,0,zm_t3)

  print*,'R T    ',trim(zm_fn)
  varname='U'
  call getdat(cdfid,varname,zm_time,0,zm_u3)

  print*,'R U    ',trim(zm_fn)
  varname='V'
  call getdat(cdfid,varname,zm_time,0,zm_v3)

  print*,'R V    ',trim(zm_fn)
  varname='P'
  call getdat(cdfid,varname,zm_time,0,zm_p3)

  print*,'R P    ',trim(zm_fn)
  varname='RLAT'
  call getdatRank2(cdfid,varname,zm_time,0,zm_rlat)

  print*,'R RLAT ',trim(zm_fn)
  varname='RLON'
  call getdatRank2(cdfid,varname,zm_time,0,zm_rlon)

  print*,'R RLON ',trim(zm_fn)
  
  call clscdf(cdfid)

  
!   Read variables from GEO.ORG
  call cdfopn(or_fn,cdfid)

  varname='P'
  call getdat(cdfid,varname,or_time,0,or_p3)

  print*,'R P    ',trim(or_fn)
  varname='U'
  call getdat(cdfid,varname,or_time,0,or_u3)

  print*,'R U    ',trim(or_fn)
  varname='V'
  call getdat(cdfid,varname,or_time,0,or_v3)

  print*,'R V    ',trim(or_fn)
  varname='T'
  call getdat(cdfid,varname,or_time,0,or_t3)

  print*,'R T    ',trim(or_fn)
  call clscdf(cdfid)
  
  !   Initialize the height levels of the Z file
  do i=1,zm_nx
     do j=1,zm_ny
        do k=1,zm_nz
           zm_z3(i,j,k)=zm_aklay(k)
        enddo
     enddo
  enddo
  
  !   Calculate 3d pressure field on model levels
  print*,'C P (ORIGINAL)'
  do k=1,ml_nz
     do i=1,ml_nx
        do j=1,ml_ny
           if (abs(ml_stag(3)+0.5).lt.eps) then
              ml_p0(i,j,k)=ml_aklay(k)+ml_bklay(k)*ml_ps(i,j)
           else
              ml_p0(i,j,k)=ml_aklev(k)+ml_bklev(k)*ml_ps(i,j)
           endif
        enddo
     enddo
  enddo

  !   Calculate 3d anomalies due to inversion
  print*,'C DP (MODIFIED-ORIGINAL)'
  print*,'C DU (MODIFIED-ORIGINAL)'
  print*,'C DV (MODIFIED-ORIGINAL)'
  print*,'C DT (MODIFIED-ORIGINAL)'
  do k=1,zm_nz
     do i=1,zm_nx
        do j=1,zm_ny
           
!            P
           if ( ( abs(zm_p3(i,j,k)-zm_mdv).gt.eps).and.&
                ( abs(or_p3(i,j,k)-or_mdv).gt.eps) ) &
                then
              an_p3(i,j,k) = zm_p3(i,j,k) - or_p3(i,j,k)
           elseif ( poisson.eq.1 ) then
              an_p3(i,j,k) = zm_mdv
           else
              an_p3(i,j,k) = 0.
           endif
           
!            T
           if ( ( abs(zm_t3(i,j,k)-zm_mdv).gt.eps).and.&
                ( abs(or_t3(i,j,k)-or_mdv).gt.eps) )&
                then
              an_t3(i,j,k) = zm_t3(i,j,k) - or_t3(i,j,k)
           elseif ( poisson.eq.1 ) then
              an_t3(i,j,k) = zm_mdv
           else
              an_t3(i,j,k) = 0.
           endif
           
!            U
           if ( ( abs(zm_u3(i,j,k)-zm_mdv).gt.eps).and.&
                ( abs(or_u3(i,j,k)-or_mdv).gt.eps) )&
                then
              an_u3(i,j,k) = zm_u3(i,j,k) - or_u3(i,j,k)
           elseif ( poisson.eq.1 ) then
              an_u3(i,j,k) = zm_mdv
           else
              an_u3(i,j,k) = 0.
           endif

!            V
           if ( ( abs(zm_v3(i,j,k)-zm_mdv).gt.eps).and. &
                ( abs(or_v3(i,j,k)-or_mdv).gt.eps) )&
                then
              an_v3(i,j,k) = zm_v3(i,j,k) - or_v3(i,j,k)
           elseif ( poisson.eq.1 ) then
              an_v3(i,j,k) = zm_mdv
           else
              an_v3(i,j,k) = 0.
           endif

        enddo
     enddo
  enddo

!   ----------------------------------------------------------------
!   Get the weight function for the inset (from 1 inside to 0 outside)
!   ----------------------------------------------------------------

  print*,'C BOUNDARY FILTER'
  
  !   Init the weight function (1 inside, 0 outside)
  do i=1,ml_nx
     do j=1,ml_ny
            
        if ( (zm_rlat(i,j).lt. -90.).or.&
             (zm_rlat(i,j).gt.  90.).or.&
             (zm_rlon(i,j).lt.-180.).or.&
             (zm_rlon(i,j).gt. 180.) ) then
           weight(i,j)=0.
        else
           weight(i,j)=1.
        endif
        
     enddo
  enddo

!   Get a list of all boundary points
  nbound=0
  do i=1,ml_nx
     do j=1,ml_ny
        
        !         Get neighbouring points
        ir=i+1
        if (ir.gt.ml_nx) ir=ml_nx
        il=i-1
        if (il.lt.    1) il=1
        ju=j+1
        if (ju.gt.ml_ny) ju=ml_ny
        jd=j-1
        if (jd.lt.    1) jd=1
        
!         A boundary point has a 0/1 switch
        if (abs(weight(i,j)-0.).lt.eps) then
               
           if ( (abs(weight(ir, j)-1.).lt.eps).or.&
                (abs(weight(il, j)-1.).lt.eps).or. &
                (abs(weight(i ,ju)-1.).lt.eps).or. &          
                (abs(weight(i ,jd)-1.).lt.eps).or. &
                (abs(weight(ir,ju)-1.).lt.eps).or. &          
                (abs(weight(ir,jd)-1.).lt.eps).or. &
                (abs(weight(il,ju)-1.).lt.eps).or. &          
                (abs(weight(il,jd)-1.).lt.eps).or. &
                (abs(weight(ir,ju)-1.).lt.eps).or. &
                (abs(weight(il,ju)-1.).lt.eps).or. &
                (abs(weight(ir,jd)-1.).lt.eps).or. &
                (abs(weight(il,jd)-1.).lt.eps) ) then
              nbound=nbound+1
              boundlon(nbound)=ml_xmin+real(i-1)*ml_dx
              boundlat(nbound)=ml_ymin+real(j-1)*ml_dy
                  
           endif
               
        endif
            
     enddo
  enddo

!   Get the distance from the subdomain
  do i=1,ml_nx
     do j=1,ml_ny
        
        lon=ml_xmin+real(i-1)*ml_dx
        lat=ml_ymin+real(j-1)*ml_dy
        
        distmin=sdis(lon,lat,boundlon(1),boundlat(1))
        do k=2,nbound
           distpos=sdis(lon,lat,boundlon(k),boundlat(k))
           if (distpos.lt.distmin) distmin=distpos
        enddo
        
        if ( abs(weight(i,j)-1.).lt.eps) then               
           dist(i,j)=distmin
        else
           dist(i,j)=-distmin
        endif
        
     enddo
  enddo

!   Set the new weights
  do i=1,ml_nx
     do j=1,ml_ny
        if (weight(i,j).gt.0.) then
           weight(i,j)=kink(dist(i,j),transxy)
        endif
        
     enddo
  enddo

!   ----------------------------------------------------------------
!   Remove MDV regions in input field
!   ----------------------------------------------------------------

  if ( poisson.ne.1 ) goto 120
  
  !   Define region for mdv filling
  do i=1,zm_nx
     do j=1,zm_ny
        do k=1,zm_nz
           
!          Get neighbour of grid point
           il = i-1
           if (il.lt.1) il=1
           ir = i+1
           if ( ir.gt.zm_nx ) ir=zm_nx
           jd = j-1
           if (jd.lt.1) jd=1
           ju = j+1
           if ( ju.gt.zm_ny ) ju=zm_ny

!          Set flag 2 for boundary and exterior points
           if ( (abs(zm_rlat(il, j)-zm_mdv).lt.eps).or.&
                (abs(zm_rlat(ir, j)-zm_mdv).lt.eps).or.&
                (abs(zm_rlat(i , j)-zm_mdv).lt.eps).or.&
                (abs(zm_rlat(il,ju)-zm_mdv).lt.eps).or.&
                (abs(zm_rlat(ir,ju)-zm_mdv).lt.eps).or.&
                (abs(zm_rlat(i ,ju)-zm_mdv).lt.eps).or.&
                (abs(zm_rlat(il,jd)-zm_mdv).lt.eps).or.&
                (abs(zm_rlat(ir,jd)-zm_mdv).lt.eps).or.&
                (abs(zm_rlat(i ,jd)-zm_mdv).lt.eps) )&
                then
              flag_zm(i,j,k)  = 2
              an_p3(i,j,k) = 0.
              an_u3(i,j,k) = 0.
              an_v3(i,j,k) = 0.
              an_t3(i,j,k) = 0.
              
              !          Set flag 1 for interior MDV points
             elseif ( abs(an_p3(i,j,k)-zm_mdv).le.eps ) then
                flag_zm(i,j,k) = 1

!          Set flag 0 for interior valid points
             else
                flag_zm(i,j,k) = 0
             endif
             
          enddo
       enddo
    enddo

!   Apply mdv filling
      print*,'C DP (POISSON FILLING)'
      call mdvfill(an_p3,an_p3,flag_zm,zm_nx,zm_ny,zm_nz,100)
      print*,'C DT (POISSON FILLING)'
      call mdvfill(an_t3,an_t3,flag_zm,zm_nx,zm_ny,zm_nz,100)
      print*,'C DU (POISSON FILLING)'
      call mdvfill(an_u3,an_u3,flag_zm,zm_nx,zm_ny,zm_nz,100)
      print*,'C DV (POISSON FILLING)'
      call mdvfill(an_v3,an_v3,flag_zm,zm_nx,zm_ny,zm_nz,100)
      
      !   Special treatment: if the number of missing values surpasses 50%
!   on a level, then the anomaly is imported from the level above
      do k=zm_nz-1,1,-1

         count1 = 0
         count0 = 0
         do i=1,zm_nx
            do j= 1,zm_ny
               if ( flag_zm(i,j,k).eq.1 ) count1 = count1 + 1
               if ( flag_zm(i,j,k).eq.0 ) count0 = count0 + 1
          enddo
       enddo

       if ( count1.gt.count0 ) then
          print*,'C P (IMPORTING FROM LEVEL ABOVE)',k
          do i=1,zm_nx
             do j= 1,zm_ny
                an_p3(i,j,k) = an_p3(i,j,k+1)
                flag_zm(i,j,k)  = flag_zm(i,j,k+1)
             enddo
          enddo
          print*,'C U (IMPORTING FROM LEVEL ABOVE)',k
          do i=1,zm_nx
             do j= 1,zm_ny
                an_u3(i,j,k) = an_u3(i,j,k+1)
                flag_zm(i,j,k)  = flag_zm(i,j,k+1)
             enddo
          enddo
          
          print*,'C V (IMPORTING FROM LEVEL ABOVE)',k
          do i=1,zm_nx
             do j= 1,zm_ny
                an_v3(i,j,k) = an_v3(i,j,k+1)
                flag_zm(i,j,k)  = flag_zm(i,j,k+1)
             enddo
          enddo
          
          print*,'C T (IMPORTING FROM LEVEL ABOVE)',k
          do i=1,zm_nx
             do j= 1,zm_ny
                an_t3(i,j,k) = an_t3(i,j,k+1)
                flag_zm(i,j,k)  = flag_zm(i,j,k+1)
             enddo
          enddo
          
       endif

    enddo

!   Write new fields for tests
    if ( test.eq.1 ) then
       
       call cdfwopn(zm_fn,cdfid)

       
       isok=0
       varname='P_ANO'
       call check_varok(isok,varname,zm_vnam,zm_nvars)
       if (isok.eq.0) then
          call putdef(cdfid,varname,zm_ndim,zm_mdv,zm_vardim,zm_varmin,zm_varmax,zm_stag)

       endif
       call putdat(cdfid,varname,zm_time,0,an_p3)

       print*,'W P_ANO ',trim(zm_fn)
       
       isok=0
       varname='T_ANO'
       call check_varok(isok,varname,zm_vnam,zm_nvars)
        if (isok.eq.0) then
         call putdef(cdfid,varname,zm_ndim,zm_mdv,zm_vardim,zm_varmin,zm_varmax,zm_stag)

      endif
      call putdat(cdfid,varname,zm_time,0,an_t3)

      print*,'W T_ANO ',trim(zm_fn)
      
      isok=0
      varname='U_ANO'
      call check_varok(isok,varname,zm_vnam,zm_nvars)
      if (isok.eq.0) then
         call putdef(cdfid,varname,zm_ndim,zm_mdv,zm_vardim,zm_varmin,zm_varmax,zm_stag)

      endif
      call putdat(cdfid,varname,zm_time,0,an_u3)

      print*,'W U_ANO ',trim(zm_fn)
      
      isok=0
      varname='V_ANO'
      call check_varok(isok,varname,zm_vnam,zm_nvars)
      if (isok.eq.0) then
         call putdef(cdfid,varname,zm_ndim,zm_mdv,zm_vardim,zm_varmin,zm_varmax,zm_stag)

      endif
      call putdat(cdfid,varname,zm_time,0,an_v3)

      print*,'W V_ANO ',trim(zm_fn)
      
      call clscdf(cdfid)

      
   endif

   !   Exit point for SOR solver
120 continue
   
   !   ----------------------------------------------------------------
   !   Change surface pressure and get 3d presure field
   !   ----------------------------------------------------------------
   
   print*,'C PS'
   
   !   Change surface pressure: based on PV inversion on GEO
   do i=1,ml_nx
      do j=1,ml_ny
         
         !         Make vertical profile of pressure available
            n1=0
            do k=1,zm_nz
               n1=n1+1
               p1(n1)=an_p3(i,j,k)
               z1(n1)=zm_z3(i,j,k)
            enddo

            if ( n1.ne.0 ) then

!         Keep the original surface pressure (psfc=0)
               if ( abs(psfc).lt.eps ) then
                  ml_ps(i,j) = ml_ps(i,j)

!         Take the full change of surface pressure (psfc=1);
!         Interpolation/extrapolation with a natural cubic spline
               elseif ( abs(psfc-1.).lt.eps ) then
                  call spline (z1,p1,n1,1.e30,1.e30,spline_f1)
                  call splint (z1,p1,spline_f1,n1,ml_oro(i,j),hh)
                  ml_ps(i,j)=ml_ps(i,j) + hh*weight(i,j)
                  
                  !         Only take a fractional change of surface pressure
                  !         Interpolation/extrapolation with a natural cubic spline
            else
               call spline (z1,p1,n1,1.e30,1.e30,spline_f1)
               call splint (z1,p1,spline_f1,n1,ml_oro(i,j),hh)
               ml_ps(i,j)=ml_ps(i,j) + hh*psfc*weight(i,j)
               
            endif

         endif
         
      enddo
   enddo

!   Calculate 3d pressure field on model levels
   print*,'C P'
   do k=1,ml_nz
      do i=1,ml_nx
         do j=1,ml_ny
            if (abs(ml_stag(3)+0.5).lt.eps) then
               ml_p3(i,j,k)=ml_aklay(k)+ml_bklay(k)*ml_ps(i,j)
            else
               ml_p3(i,j,k)=ml_aklev(k)+ml_bklev(k)*ml_ps(i,j)
            endif
         enddo
      enddo
   enddo
   
   !   ----------------------------------------------------------------
!   Get T,U,V at the new pressure levels, based on the original P file
!   ----------------------------------------------------------------

!   Interpolate T from original P file to new pressure levels
   print*,'C T (ORIGINAL)'
   do i=1,ml_nx
      do j=1,ml_ny
         
         n1=0
         do k=ml_nz,1,-1
            n1=n1+1
            f1(n1)=ml_t3(i,j,k)
            p1(n1)=ml_p0(i,j,k)
         enddo
         call spline (p1,f1,n1,1.e30,1.e30,spline_f1)
         do k=1,ml_nz
            call splint (p1,f1,spline_f1,n1,ml_p3(i,j,k),hh)
            ml_t3(i,j,k)=hh
         enddo

      enddo
   enddo
      
!   Interpolate U from original P file to new pressure levels
   print*,'C U (ORIGINAL)'
   do i=1,ml_nx
      do j=1,ml_ny
         n1=0
         do k=ml_nz,1,-1
            n1=n1+1
            f1(n1)=ml_u3(i,j,k)
            p1(n1)=ml_p0(i,j,k)
         enddo
         call spline (p1,f1,n1,1.e30,1.e30,spline_f1)
         do k=1,ml_nz
            call splint (p1,f1,spline_f1,n1,ml_p3(i,j,k),hh)
            ml_u3(i,j,k)=hh
         enddo
         
      enddo
   enddo

!   Interpolate V from original P file to new pressure levels
   print*,'C V (ORIGINAL)'
   do i=1,ml_nx
      do j=1,ml_ny
         
         n1=0
         do k=ml_nz,1,-1
            n1=n1+1
            f1(n1)=ml_v3(i,j,k)
            p1(n1)=ml_p0(i,j,k)
         enddo
         call spline (p1,f1,n1,1.e30,1.e30,spline_f1)
         do k=1,ml_nz
            call splint (p1,f1,spline_f1,n1,ml_p3(i,j,k),hh)
            ml_v3(i,j,k)=hh
         enddo

      enddo
   enddo
   
!   ----------------------------------------------------------------
!   Add T,U,V anomalies at the new pressure levels
!   ----------------------------------------------------------------

!   Change temperature field
   print*,'C T (ANOMALY)'
   do i=1,ml_nx
      do j=1,ml_ny
         
         !         Make vertical profile of temperature available
         n1=0
         pmax=-100.
         pmin=2000.
         do k=zm_nz,1,-1
            if ((abs(an_t3(i,j,k)-zm_mdv).gt.eps).and.&
                 (zm_z3(i,j,k).gt.ml_oro(i,j)).and.&
                 (abs(zm_p3(i,j,k)-zm_mdv).gt.eps)) then
               n1=n1+1
               f1(n1)=an_t3(i,j,k)
               p1(n1)=zm_p3(i,j,k)
               if (p1(n1).lt.pmin) pmin=p1(n1)
               if (p1(n1).gt.pmax) pmax=p1(n1)
            endif
         enddo
         pmin=pmin-dpextra*pmin
         pmax=pmax+dpextra+pmax
         
         !         Cubic spline interpolation
         if (n1.gt.0) then
            call spline (p1,f1,n1,1.e30,1.e30,spline_f1)
            do k=1,ml_nz
               if ( (ml_p3(i,j,k).gt.pmin).and.&
                    (ml_p3(i,j,k).lt.pmax) ) then
                  call splint (p1,f1,spline_f1,n1,ml_p3(i,j,k),hh)
                  ml_t3(i,j,k)=ml_t3(i,j,k) + hh*weight(i,j)
               endif
            enddo
         endif
         
      enddo
   enddo
      
      !   Change zonal velocity field
   print*,'C U'
   do i=1,ml_nx
      do j=1,ml_ny
         
         !         Make vertical profile of zonal velocity available
         n1=0
         pmax=-100.
         pmin=2000.
         do k=zm_nz,1,-1
            if ((abs(an_u3(i,j,k)-zm_mdv).gt.eps).and.&
                 (zm_z3(i,j,k).gt.ml_oro(i,j)).and.&
                 (abs(zm_p3(i,j,k)-zm_mdv).gt.eps)) then
               n1=n1+1
               f1(n1)=an_u3(i,j,k)
               p1(n1)=zm_p3(i,j,k)
               if (p1(n1).lt.pmin) pmin=p1(n1)
               if (p1(n1).gt.pmax) pmax=p1(n1)
            endif
         enddo
         pmin=pmin-dpextra*pmin
         pmax=pmax+dpextra*pmax
         
            !         Cubic spline interpolation
         if (n1.gt.0) then
            call spline (p1,f1,n1,1.e30,1.e30,spline_f1)
            do k=1,ml_nz
               if ( (ml_p3(i,j,k).gt.pmin).and.&
                    (ml_p3(i,j,k).lt.pmax) ) then
                  call splint (p1,f1,spline_f1,n1,ml_p3(i,j,k),hh)
                  ml_u3(i,j,k)=ml_u3(i,j,k) + hh*weight(i,j)
               endif
            enddo
         endif
            
      enddo
   enddo

!   Change meridional velocity field
   print*,'C V'
   do i=1,ml_nx
      do j=1,ml_ny
         
         !         Make vertical profile of zonal velocity available
         n1=0
         pmax=-100.
         pmin=2000.
         do k=zm_nz,1,-1
            if ((abs(an_v3(i,j,k)-zm_mdv).gt.eps).and.&
                 (zm_z3(i,j,k).gt.ml_oro(i,j)).and.&
                 (abs(zm_p3(i,j,k)-zm_mdv).gt.eps)) then
               n1=n1+1
               f1(n1)=an_v3(i,j,k)
               p1(n1)=zm_p3(i,j,k)
               if (p1(n1).lt.pmin) pmin=p1(n1)
               if (p1(n1).gt.pmax) pmax=p1(n1)
            endif
         enddo
         pmin=pmin-dpextra*pmin
         pmax=pmax+dpextra*pmax
         
         !         Cubic spline interpolation
         if (n1.gt.0) then
            call spline (p1,f1,n1,1.e30,1.e30,spline_f1)
            do k=1,ml_nz
               if ( (ml_p3(i,j,k).gt.pmin).and.&
                    (ml_p3(i,j,k).lt.pmax) ) then
                  call splint (p1,f1,spline_f1,n1,ml_p3(i,j,k),hh)
                  ml_v3(i,j,k)=ml_v3(i,j,k) + hh*weight(i,j)
               endif
            enddo
         endif
         
      enddo
   enddo

!   ---------------------------------------------------------------------
!   Change geopotential height
!   ---------------------------------------------------------------------

!   Interpolate Z from original P file to new pressure levels
   print*,'C Z (ORIGINAL)'
   do i=1,ml_nx
      do j=1,ml_ny
         
         n1=0
         do k=ml_nz,1,-1
            n1=n1+1
            f1(n1)=ml_z3(i,j,k)
            p1(n1)=ml_p0(i,j,k)
         enddo
         call spline (p1,f1,n1,1.e30,1.e30,spline_f1)
         do k=1,ml_nz
            call splint (p1,f1,spline_f1,n1,ml_p3(i,j,k),hh)
            ml_z3(i,j,k)=hh
         enddo
            
      enddo
   enddo
   
   !   Calculate 3d virtual temperature
   print*,'C TV'
   do k=1,ml_nz
      do i=1,ml_nx
         do j=1,ml_ny
            ml_tv3(i,j,k) = (ml_t3(i,j,k)+tzero)*&
                 (1.+kappa*ml_q3(i,j,k))
         enddo
      enddo
   enddo
   
!   Add geopotential height anomaly: first, the pressure anomaly is
   !   interpolated to the new position, then it is transformed into
   !   a correction of geopotential height with the hydrostatic equation
   print*,'C DZ (MODIFIED -ORIGINAL)'
   do i=1,ml_nx
      do j=1,ml_ny
         
            !         Make vertical profile of pressure available
         n1=0
         pmax=-100.
         pmin=2000.
         do k=zm_nz,1,-1
            if ((abs(an_p3(i,j,k)-zm_mdv).gt.eps).and.&
                 (zm_z3(i,j,k).gt.ml_oro(i,j)).and.&
                 (abs(zm_p3(i,j,k)-zm_mdv).gt.eps)) then
               n1=n1+1
               f1(n1)=an_p3(i,j,k)
               p1(n1)=zm_p3(i,j,k)
               if (p1(n1).lt.pmin) pmin=p1(n1)
               if (p1(n1).gt.pmax) pmax=p1(n1)
            endif
         enddo
         pmin=pmin-dpextra*pmin
         pmax=pmax+dpextra*pmax
         
         !         Cubic spline interpolation and conversion dp -> dz
         if (n1.gt.0) then
            call spline (p1,f1,n1,1.e30,1.e30,spline_f1)
            do k=1,ml_nz
               if ( (ml_p3(i,j,k).gt.pmin).and.&
                    (ml_p3(i,j,k).lt.pmax) ) then
                  call splint (p1,f1,spline_f1,n1,ml_p3(i,j,k),hh)
                  hh = -rdg * ml_tv3(i,j,k) * hh/ml_p3(i,j,k)
                  ml_z3(i,j,k)=ml_z3(i,j,k) + hh*weight(i,j)
               endif
            enddo
         endif
         
      enddo
   enddo
      
      !   ----------------------------------------------------------------
      !   Remove unrealistic values from final fields
      !   ----------------------------------------------------------------

   if ( poisson.ne.1 ) goto 120
   
   !   Define region for mdv filling
   count1 = 0
   
   do i=1,ml_nx
      do j=1,ml_ny
         do k=1,ml_nz
               
            !          Get neighbour of grid point
            il = i-1
            if (il.lt.1) il=1
            ir = i+1
            if ( ir.gt.ml_nx ) ir=ml_nx
            jd = j-1
            if (jd.lt.1) jd=1
            ju = j+1
            if ( ju.gt.ml_ny ) ju=ml_ny

!          Set flag 2 for boundary and exterior points
            if ( (abs(zm_rlat(il, j)-zm_mdv).lt.eps).or.&
                 (abs(zm_rlat(ir, j)-zm_mdv).lt.eps).or.&
                 (abs(zm_rlat(i , j)-zm_mdv).lt.eps).or.&
                 (abs(zm_rlat(il,ju)-zm_mdv).lt.eps).or.&
                 (abs(zm_rlat(ir,ju)-zm_mdv).lt.eps).or.&
                 (abs(zm_rlat(i ,ju)-zm_mdv).lt.eps).or.&
                 (abs(zm_rlat(il,jd)-zm_mdv).lt.eps).or.&
                 (abs(zm_rlat(ir,jd)-zm_mdv).lt.eps).or.&
                 (abs(zm_rlat(i ,jd)-zm_mdv).lt.eps) )&
                 then
               flag_ml(i,j,k)  = 2
               
               !          Set flag 1 for interior unphysical points
            elseif ( ( abs(ml_t3(i,j,k)).gt.500. ).or.&
                 ( abs(ml_u3(i,j,k)).gt.500. ).or.&
                 ( abs(ml_v3(i,j,k)).gt.500. ) ) &
                 then
               flag_ml(i,j,k) = 1
               count1 = count1 + 1
               
               !          Set flag 0 for interior valid points
            else
               flag_ml(i,j,k) = 0
            endif
            
         enddo
      enddo
   enddo

   print*,'C MASK UNREALISTIC VALUES FOR T,U,V',count1
   
   !   Apply mdv filling
   print*,'C T (POISSON FILLING)'
   call mdvfill(ml_t3,ml_t3,flag_ml,ml_nx,ml_ny,ml_nz,10)
   print*,'C U (POISSON FILLING)'
   call mdvfill(ml_u3,ml_u3,flag_ml,ml_nx,ml_ny,ml_nz,10)
   print*,'C V (POISSON FILLING)'
   call mdvfill(ml_v3,ml_v3,flag_ml,ml_nx,ml_ny,ml_nz,10)
   
!   ----------------------------------------------------------------
!   Write modified fields to P file
!   ----------------------------------------------------------------

!   Open output file
   call cdfwopn(ml_fn,cdfid)

   
   !   Write surface pressure
   isok=0
   varname='PS'
   call check_varok(isok,varname,ml_vnam,ml_nvars)
   if (isok.eq.0) then
      ml_vardim(3)=3
      call putdef(cdfid,varname,ml_ndim,ml_mdv,ml_vardim,ml_varmin,ml_varmax,ml_stag)

      ml_vardim(3)=4
   endif
   call putdatRank2(cdfid,varname,ml_time,0,ml_ps)

   print*,'W PS ',trim(ml_fn)

!   Write temperature
   isok=0
   varname='T'
   call check_varok(isok,varname,ml_vnam,ml_nvars)
   if (isok.eq.0) then
      call putdef(cdfid,varname,ml_ndim,ml_mdv,ml_vardim,ml_varmin,ml_varmax,ml_stag)

   endif
   call putdat(cdfid,varname,ml_time,0,ml_t3)

   print*,'W T  ',trim(ml_fn)
   
   !   Write zonal wind
   isok=0
   varname='U'
   call check_varok(isok,varname,ml_vnam,ml_nvars)
   if (isok.eq.0) then
      call putdef(cdfid,varname,ml_ndim,ml_mdv,ml_vardim,ml_varmin,ml_varmax,ml_stag)

   endif
   call putdat(cdfid,varname,ml_time,0,ml_u3)
   print*,'W U  ',trim(ml_fn)
   
   !   Write meridional wind
   isok=0
   varname='V'
   call check_varok(isok,varname,ml_vnam,ml_nvars)
   if (isok.eq.0) then
      call putdef(cdfid,varname,ml_ndim,ml_mdv,ml_vardim,ml_varmin,ml_varmax,ml_stag)

   endif
   call putdat(cdfid,varname,ml_time,0,ml_v3)

   print*,'W V  ',trim(ml_fn)
   
   !   Write geopotential height
   isok=0
   varname='Z'
   call check_varok(isok,varname,ml_vnam,ml_nvars)
   if (isok.eq.0) then
      call putdef(cdfid,varname,ml_ndim,ml_mdv,ml_vardim,ml_varmin,ml_varmax,ml_stag)

   endif
   call putdat(cdfid,varname,ml_time,0,ml_z3)

   print*,'W Z  ',trim(ml_fn)
   
   !   Filter matrix (only in test mode)
   if (test.eq.1) then
      isok=0
      varname='DIST'
      call check_varok(isok,varname,ml_vnam,ml_nvars)
      if (isok.eq.0) then
         ml_vardim(3)=1
         call putdef(cdfid,varname,ml_ndim,ml_mdv,ml_vardim,ml_varmin,ml_varmax,ml_stag)

      endif
      call putdatRank2(cdfid,varname,ml_time,1,dist)

      print*,'W DIST ',trim(ml_fn)
      
      isok=0
      varname='WEIGHT'
      call check_varok(isok,varname,ml_vnam,ml_nvars)
      if (isok.eq.0) then
         ml_vardim(3)=1
         call putdef(cdfid,varname,ml_ndim,ml_mdv,ml_vardim,ml_varmin,ml_varmax,ml_stag)

      endif
      call putdatRank2(cdfid,varname,ml_time,1,weight)

      print*,'W WEIGHT ',trim(ml_fn)
      
   endif
      
      !   Close output file
   call clscdf(cdfid)

   
      !   ----------------------------------------------------------------
      !   Exception handling
      !   ----------------------------------------------------------------
      
   stop

 998  print*,'Problem with input netcdf file ',trim(ml_fn)
      stop

997   print*,'Problem with input netcdf file ',trim(zm_fn)
      stop

 994  print*,'Problem with output netcdf file ',trim(ml_fn)
      stop

   
 end PROGRAM add2p
 

!   ****************************************************************
!   * SUBROUTINE SECTION: AUXILIARY ROUTINES                       *
!   ****************************************************************

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

!   -------------------------------------------------------------
!   Natural cubic spline /from Numerical Recipes)
!   -------------------------------------------------------------

 SUBROUTINE spline(x,y,n,yp1,ypn,y2)
  use ieee_arithmetic
  use ieee_features
  use kind_parameters,ONLY:&
       sp
  INTEGER n,NMAX
  logical :: underflow_support, gradual, underflow
  REAL(kind=sp) yp1,ypn,x(n),y(n),y2(n)
  PARAMETER (NMAX=500)
  INTEGER i,k
  REAL(kind=sp) p,qn,sig,un,u(NMAX)

  call ieee_set_underflow_mode(.false.)
  call ieee_get_underflow_mode(gradual)
  
  if (yp1.gt..99e30) then
     y2(1)=0.
     u(1)=0.
  else
     y2(1)=-0.5
     u(1)=(3./(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
  endif
  do  i=2,n-1
     sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
     p=sig*y2(i-1)+2.
     y2(i)=(sig-1.)/p
     u(i)=(6.*((y(i+1)-y(i))/(x(i+&
          1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*&
          u(i-1))/p
     call ieee_get_flag(ieee_underflow,underflow)                    
  end do
  if (ypn.gt..99e30) then
     qn=0.
     un=0.
  else
     qn=0.5
     un=(3./(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
  endif
  y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.)
  do  k=n-1,1,-1
     y2(k)=y2(k)*y2(k+1)+u(k)
  end do
     
END SUBROUTINE spline
 
SUBROUTINE splint(xa,ya,y2a,n,x,y)
   use kind_parameters,ONLY:&
       sp

  INTEGER n
  REAL(kind=sp) x,y,xa(n),y2a(n),ya(n)
  INTEGER k,khi,klo
  REAL(kind=sp) a,b,h
  klo=1
  khi=n
1 if (khi-klo.gt.1) then
     k=(khi+klo)/2
     if(xa(k).gt.x)then
        khi=k
     else
        klo=k
     endif
     goto 1
  endif
  h=xa(khi)-xa(klo)
  if (h.eq.0.) pause 'bad xa input in splint'
  a=(xa(khi)-x)/h
  b=(x-xa(klo))/h
  y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.
  return
END SUBROUTINE splint


!   ---------------------------------------------------------
!   Spherical distance between two lat/lon points
!   ---------------------------------------------------------

real function sdis(xp,yp,xq,yq)
!
!   calculates spherical distance (in km) between two points given
!   by their spherical coordinates (xp,yp) and (xq,yq), respectively.
  !
  use kind_parameters,ONLY:&
       sp

  real(kind=sp)      re,pi,dconv
  parameter (re=6370.)
  real(kind=sp)      xp,yp,xq,yq,arg
  parameter(pi=3.14159265)
  dconv=pi/180
  arg=sin(yp*dconv)*sin(yq*dconv)+cos(yp*dconv)*cos(yq*dconv)*cos((xp-xq)*dconv)
  if (arg.lt.-1.) arg=-1.
  if (arg.gt.1.) arg=1.
  sdis=re*acos(arg)
  
end function sdis

!   -----------------------------------------------------------------
!   Rene's KINK function (for smoothing at bounadries)
!   -----------------------------------------------------------------

real function kink (x,a)
  
  implicit none

  !   declaration of parameters
  real   x,a
  
      !   parameters
  real   pi
  parameter  (pi=3.1415926535)
  
  if (x.lt.0.) then
     kink=0.
  elseif (x.gt.a) then
     kink=1.
  else
     kink=(1.+cos(pi*(x-a)/a))/2.
  endif
  
  return
end function kink

!   -----------------------------------------------------------------
!   Rene's KINK function (for smoothing at bounadries)
!   -----------------------------------------------------------------

subroutine mdvfill (out,inp,flag,nx,ny,nz,maxiter)
  
  implicit none
  
  !   Declaration of subroutine parameters
  integer nx,ny,nz
  real    inp (nx,ny,nz)
  real    out (nx,ny,nz)
  integer flag(nx,ny,nz)
  integer maxiter
  
  !   Parameters
  real             omega        ! Omega fopr SOR
  parameter        (omega=1.5)
  
  !   Auxiliary variables
  integer i,j,k
  integer iter
  real    tmp0(nx,ny,nz)
  real    tmp1(nx,ny,nz)
  integer il,ir,ju,jd
  integer count
  real    mean(nz)
  
  !   Calculate mean of variable for all levels
  do k=1,nz
     mean(k) = 0.
     count   = 0
     do i=1,nx
        do j=1,ny
           if ( flag(i,j,k).eq.0 ) then
              count   = count + 1
              mean(k) = mean(k) + inp(i,j,k)
           endif
        enddo
     enddo
     if ( count.ne.0 ) then
        mean(k) = mean(k)/real(count)
     else
        mean(k) = 0.
     endif
  enddo

  !   Create first guess
  do k=1,nz
     do i=1,nx
        do j=1,ny
           if ( flag(i,j,k).eq.0 ) then
              tmp0(i,j,k) = inp(i,j,k)
           else
              tmp0(i,j,k) = mean(k)
        endif
     enddo
  enddo
enddo
   
!   SOR iterations
iter = 0
122 continue

!   Loop over all points
do k=1,nz
   do i=1,nx
      do j=1,ny
         
         !       Apply the updating only for specified points
          if ( flag(i,j,k).ne.1 ) goto 121
          
          !       Get neighbouring points - no handling of periodic domains!
          il = i-1
          if (il.lt.1) il=1
          ir = i+1
          if ( ir.gt.nx ) ir=nx
          jd = j-1
          if (jd.lt.1) jd=1
          ju = j+1
          if ( ju.gt.ny ) ju=ny

!       Update field
          tmp1(i,j,k) = 0.25 * ( tmp0(il,j,k) + tmp0(ir,j,k) +&
               tmp0(i,ju,k) + tmp0(i,jd,k) )
          
          tmp0(i,j,k) = omega * tmp1(i,j,k) +&
               (1. - omega) * tmp1(i,j,k)
          
          !       Exit point for loop
 121      continue
          
       enddo
    enddo
 enddo
 
!   Decide whether further iterations are needed
 iter = iter + 1
 if ( iter.lt.maxiter) goto 122
 
 !   Save output
 do i=1,nx
    do j=1,ny
       do k=1,nz
          if ( flag(i,j,k).eq.1 ) then
             out(i,j,k) = tmp0(i,j,k)
          endif
       enddo
    enddo
 enddo
 
end subroutine mdvfill
      

