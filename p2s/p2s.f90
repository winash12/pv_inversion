program ptos
  use netcdflibrary
  implicit none
!   ******************************************************************

!   NAME:
!   p2s

!   Purpose
!   -------

!   	Calculates secondary data files from primary data files
!   	(based upon IVE-routines).

!   calling:
!   p2s [-m] file variable-list [-s] [-o]

!   example:
!   p2s P911201_00 TH PV CH QX

!   p2s -m #returns man-page


!   ADD YOUR OWN VARIABLES AT ALL PLACES WITH (++)
!   for easy simple calculations, see VEL, M, B
!   for complicated calculations, see VORT, QX and PVR

!   Remarks
!   -------
!   ZB is only read once when needed, i.e. for first time...


  integer,parameter :: ntmax=200,nzmax=200,nvarmax=100

  real,allocatable,dimension(:)::time,time2
  REAL,DIMENSION (:,:),allocatable :: sp1,cl,tl,f,zb,t2m,td2m,vip,u10m,&
       v10m,oro,gradpv
  REAL,DIMENSION (:,:,:),allocatable:: var,th,pv,lpv,the,rh,dhr,&
       tt,qq,uu,vv,ww,rho,alpha,zz,mm,zlay,ug,vg,fl,ipv
  character*80 cdfnam,cstnam,outfnam
  integer	cdfid,cdfid1,cstid,ierr,ndim,vardim(4),stat
  integer	cdfid2,vardim2(4)
  real	dx,dy,mdv
  real,dimension(:),allocatable::varmin,varmax,stag
  real aklev(nzmax),bklev(nzmax),aklay(nzmax),bklay(nzmax),&
       ak(nzmax),bk(nzmax)
  integer	nx,ny,nz,ntimes,ntimes2,i,j,k,n
  integer	stdate(5)
  character*(80) qmode,arg,vnam(nvarmax)
  integer	mode,zdef,nvars
  real	rlat,rlon,lat
  real	pollon,pollat,yphys
  real	phstoph
  real akd,bkd,VEL_calc,VEL_out,FL_calc,FL_out
  logical	prelev
  
  real,parameter ::  pi=3.141592654
  
  integer   PS_out,TH_out,TH_calc,RH_out,RH_calc
  integer   PV_out,PV_calc,THE_out,THE_calc,VIP_out,VIP_calc
  integer   GRADPV_out,GRADPV_calc
  integer	LPV_out,LPV_calc
  integer   CH_out,CH_calc,PVR_out,PVR_calc
  integer   THW_out,THW_calc,Z_outP,Z_out,Z_calc
  integer	DIVQU_out,DIVQU_calc
  integer   NSQ_out,NSQ_calc,RHO_out,RHO_calc,ALPHA_out,ALPHA_calc
  integer   NSQM_out,NSQM_calc,W_out,W_calc,M_out,M_calc
  integer   VORT_out,VORT_calc,UG_out,UG_calc,VG_out,VG_calc
  integer	AVO_out,AVO_calc,CURVO_out,CURVO_calc
  integer	DTHDP_out,DTHDP_calc
  integer	COS_out
  integer   ZLAY_out,ZLAY_calc,UA_out,UA_calc,VA_out,VA_calc
  integer   P_out,P_calc,PLEV_out,PLEV_calc
  integer   QXF_out,QXF_calc,QYF_out,QYF_calc
  integer   QX_out,QX_calc,QY_out,QY_calc,PSRED_calc,PSRED_out
  integer	RI_calc,RI_out,BLH_calc,BLH_out
  integer   GRADTH_out,GRADTH_calc,B_out,B_calc !(++)
  logical   verbose
  logical   ZonP,TonP,PSonP,UonP,VonP,OMEGAonP,ZBonP,ZBneed
  logical	T2MonP,TD2MonP,U10MonP,V10MonP,PVonP,PV3onP
  character*(80) zbfile
  
  !     set defaults:
  allocate(time(1))
  allocate(time2(1))
  allocate(varmin(4))
  allocate(varmax(4))
  allocate(stag(4))
  verbose=.false.
  ZonP=.false.
  TonP=.false.
  PSonP=.false.
  ZBonP=.false.
  UonP=.false.
  VonP=.false.
  OMEGAonP=.false.
  ZBonP=.false.
  T2MonP=.false.
  TD2MonP=.false.
  U10MonP=.false.
  V10MonP=.false.
  PVonP=.false.
  PV3onP=.false.	! Lukas
  ZBneed=.false.
  zbfile=''

  qmode='QNone'
  PS_out=1
  TH_out=0
  TH_calc=0
  RH_out=0
  RH_calc=0
  PV_out=0
  PV_calc=0
  LPV_out=0
  LPV_calc=0
  VIP_out=0
  VIP_calc=0
  THE_out=0
  THE_calc=0
  CH_out=0
  CH_calc=0
  PVR_out=0
  PVR_calc=0
  THW_out=0
  THW_calc=0
  DIVQU_out=0
  DIVQU_calc=0
  Z_calc=0
  Z_out=0
  Z_outP=0
  NSQ_out=0
  NSQ_calc=0
  DTHDP_out=0
  DTHDP_calc=0
  NSQM_out=0
  NSQM_calc=0
  RHO_out=0
  RHO_calc=0
  ALPHA_out=0
  ALPHA_calc=0
  VEL_out=0
  VEL_calc=0
  M_out=0
  M_calc=0
  B_out=0
  B_calc=0
  W_out=0
  W_calc=0
  VORT_out=0
  VORT_calc=0
  AVO_out=0
  AVO_calc=0
  CURVO_out=0
  CURVO_calc=0
  COS_out=0
  UG_out=0
  UG_calc=0
  VG_out=0
  VG_calc=0
  UA_out=0
  UA_calc=0
  VA_out=0
  VA_calc=0
  ZLAY_out=0
  ZLAY_calc=0
  P_out=0
  P_calc=0
  PLEV_out=0
  PLEV_calc=0
  FL_out=0
  FL_calc=0
  PSRED_out=0
  PSRED_calc=0
  RI_out=0
  RI_calc=0
  BLH_out=0
  BLH_calc=0
  QXF_out=0
  QXF_calc=0
  QYF_out=0
  QYF_calc=0
  QX_out=0
  QX_calc=0
  QY_out=0
  QY_calc=0
  GRADTH_out=0
  GRADTH_calc=0
  GRADPV_out=0
  GRADPV_calc=0             !(++)

  
  if (iargc() .lt. 1) then
     print*,'USAGE: p2s [-m] file variable-list [-s] ',&
          '[-o] [-zb file]'
     STOP
  endif
  
  !     REQUESTD INPUT:
  !     ---------------
  !     GET WITH getarg DIRECTLY FROM SHELL:
  call getarg(1,cdfnam)
  if (trim(cdfnam).eq.'-m') then
     print*,' '
     print*,'computes derived variables from primary ones on'
     print*,'the input-file'
     print*,'if the output-file is already present, it will be'
     print*,'updated'
     print*,' '
     print*,'check the source-code itself about details...'
     print*,' '
     print*,'file: netCDF file with basic variables on it'
     print*,'    requested are the variables needed to calculate'
     print*,'    the requested output (variable-list)'
     print*,'    If file is a P-file (starting with P), the output-'
     print*,'    file will be an S-file, otherwise the extension'
     print*,'    _out will be appended, unless -o is used'
     print*,'    if the S-file exists already, it is tried to '
     print*,'    append the new variable'
     print*,'[P]date means that you can either give PYYMMDD_HH'
     print*,'    or YYMMDD_HH alone'
     print*,'variable-list: a list of variables to be calculated'
     print*,'    and written to the S_file, available are:'
     print*,'    TH,PV,LPV,RH,THE,THW,CH,PVR,Z,ZonP,GRADTH,NSQ,NSQM'
     print*,'    M,B,W,RHO,VEL,VORT,AVO,CURVO,UG,VG,ZLAY,UA,VA,P'
     print*,'    PLEV,FL,QX,QY,QXF,QYF,PSRED,RI,BLH,GRADPV,DTHDP'
     print*,'    ALPHA'
     print*,'-s: only small S-file, i.e. TH and PV'
     print*,'-zb: file with ZB (for PSRED)'
     print*,'-o output: filename of the output netCDF file'
     STOP
  endif
!   check, if cdfnam is with or without P:
  if (cdfnam(1:1).eq.'P') then
     outfnam='S'//cdfnam(2:len_trim(cdfnam))
  else        
     outfnam=trim(cdfnam)//'_out'
  endif
  i=2
  do while (iargc().ge.i)
     call getarg(i,arg)
     i=i+1
     if (arg.eq.'TH') TH_out=1
     if (arg.eq.'THE') THE_out=1
     if (arg.eq.'THW') THW_out=1
     if (arg.eq.'DIVQU') DIVQU_out=1
     if (arg.eq.'PV') PV_out=1
     if (arg.eq.'LPV') LPV_out=1
     if (arg.eq.'VIP') VIP_out=1
     if (arg.eq.'RH') RH_out=1
     if (arg.eq.'CH') CH_out=1
     if (arg.eq.'PVR') PVR_out=1
     if (arg.eq.'Z') Z_out=1
     if (arg.eq.'ZonP') Z_outP=1
     if (arg.eq.'GRADTH') GRADTH_out=1
     if (arg.eq.'GRADPV') GRADPV_out=1
     if (arg.eq.'VEL') VEL_out=1
     if (arg.eq.'RHO') RHO_out=1
     if (arg.eq.'ALPHA') ALPHA_out=1
     if (arg.eq.'W') W_out=1
     if (arg.eq.'M') M_out=1
     if (arg.eq.'B') B_out=1
     if (arg.eq.'VORT') VORT_out=1
     if (arg.eq.'AVO') AVO_out=1
     if (arg.eq.'CURVO') CURVO_out=1
     if (arg.eq.'COS') COS_out=1
     if (arg.eq.'UG') UG_out=1
     if (arg.eq.'VG') VG_out=1         
     if (arg.eq.'UA') UA_out=1
     if (arg.eq.'VA') VA_out=1         
     if (arg.eq.'QX') QX_out=1         
     if (arg.eq.'QY') QY_out=1         
     if (arg.eq.'QXF') QXF_out=1         
     if (arg.eq.'QYF') QYF_out=1     
     if (arg.eq.'ZLAY') ZLAY_out=1
     if (arg.eq.'P') P_out=1
     if (arg.eq.'PLEV') PLEV_out=1
     if (arg.eq.'FL') FL_out=1
     if (arg.eq.'PSRED') PSRED_out=1
     if (arg.eq.'RI') RI_out=1
     if (arg.eq.'BLH') BLH_out=1
     if (arg.eq.'NSQ') NSQ_out=1
     if (arg.eq.'DTHDP') DTHDP_out=1
     if (arg.eq.'NSQM') NSQM_out=1 ! (++)      
     if (arg.eq.'-v') verbose=.true.
     if (arg.eq.'-s') then
        TH_out=1
        PV_out=1
     endif
     if (arg.eq.'-o') then
        if (iargc().ge.i) then          
           call getarg(i,outfnam)
           i=i+1
        else
           print*,'option -o requires filename'
           STOP
        endif
     endif
     if (arg.eq.'-zb') then
        if (iargc().ge.i) then          
           call getarg(i,zbfile)
           ZBneed=.true.
           i=i+1
        else
           print*,'option -zb requires filename'
           STOP
        endif
     endif
     if (arg.eq.'-nops') then
        PS_out=0
     endif
  enddo
  
  !   force calculations of requested fields:
  
  if (QX_out.eq.1) QX_calc=1
  if (QY_out.eq.1) QY_calc=1
  if (QXF_out.eq.1) QXF_calc=1
  if (QYF_out.eq.1) QYF_calc=1
  if (ZLAY_out.eq.1) ZLAY_calc=1
  if (P_out.eq.1) P_calc=1
  if (PLEV_out.eq.1) PLEV_calc=1
  if (FL_out.eq.1) FL_calc=1
  if (UG_out.eq.1) UG_calc=1
  if (VG_out.eq.1) VG_calc=1
  if (UA_out.eq.1) UA_calc=1
  if (VA_out.eq.1) VA_calc=1
  if (VORT_out.eq.1) VORT_calc=1
  if (AVO_out.eq.1) AVO_calc=1
  if (CURVO_out.eq.1) CURVO_calc=1
  if (TH_out.eq.1) TH_calc=1
  if (PV_out.eq.1) PV_calc=1
  if (LPV_out.eq.1) LPV_calc=1
  if (VIP_out.eq.1) VIP_calc=1
  if (RH_out.eq.1) RH_calc=1
  if (THE_out.eq.1) THE_calc=1
  if (THW_out.eq.1) THW_calc=1
  if (DIVQU_out.eq.1) DIVQU_calc=1
  if (CH_out.eq.1) CH_calc=1
  if (PVR_out.eq.1) PVR_calc=1
  if (Z_out.eq.1) Z_calc=1
  if (Z_outP.eq.1) Z_calc=1
  if (GRADTH_out.eq.1) GRADTH_calc=1
  if (GRADPV_out.eq.1) GRADPV_calc=1
  if (RHO_out.eq.1) RHO_calc=1
  if (ALPHA_out.eq.1) ALPHA_calc=1
  if (RI_out.eq.1) RI_calc=1
  if (BLH_out.eq.1) BLH_calc=1
  if (VEL_out.eq.1) VEL_calc=1
  if (M_out.eq.1) M_calc=1
  if (B_out.eq.1) B_calc=1
  if (NSQM_out.eq.1) NSQM_calc=1
  if (W_out.eq.1) W_calc=1
  if (DTHDP_out.eq.1) DTHDP_calc=1
  if (NSQ_out.eq.1) NSQ_calc=1 !(++)
  
  !   make dependencies for variable calculations:
  
  if (NSQ_calc.eq.1) TH_calc=1 !(++)
  if (DTHDP_calc.eq.1) TH_calc=1
  if (PSRED_out.eq.1) PSRED_calc=1
  if (PSRED_calc.eq.1) ZBneed=.true.
  if (QX_out.eq.1) UG_calc=1
  if (QY_out.eq.1) VG_calc=1
  if (UA_calc.eq.1) UG_calc=1
  if (VA_calc.eq.1) VG_calc=1
  if (UG_calc.eq.1) ZLAY_calc=1
  if (VG_calc.eq.1) ZLAY_calc=1
  if (ZLAY_calc.eq.1) Z_calc=1
  if (B_calc.eq.1) M_calc=1
  if (M_calc.eq.1) ZLAY_calc=1
  if (NSQ_calc.eq.1) RHO_calc=1
  if (NSQM_calc.eq.1) RHO_calc=1
  if (NSQM_calc.eq.1) THE_calc=1
  if (W_calc.eq.1) RHO_calc=1
  if (GRADTH_calc.eq.1) TH_calc=1
  if (THW_calc.eq.1) THE_calc=1
  if (THE_calc.eq.1) TH_calc=1
  if (VIP_calc.eq.1) PV_calc=1
  if (PV_calc.eq.1) TH_calc=1
  if (LPV_calc.eq.1) PV_calc=1
  if (PVR_calc.eq.1) CH_calc=1
  if (CH_calc.eq.1) TH_calc=1
  if (CH_calc.eq.1) RH_calc=1
  if (RI_calc.eq.1) TH_calc=1
  if (ALPHA_calc.eq.1) RHO_calc=1
  
  print*,'processing: ',trim(cdfnam)

  !   Open files and get infos about data domain
  if (Z_outP.eq.1) then
     call cdfwopn2(trim(cdfnam),cdfid1,ierr)
  else
     call cdfopn2(trim(cdfnam),cdfid1,ierr)
  endif

  if (ierr.ne.0) then
     print*,'ERROR opening input file, stopped'
     stop
  endif

  call getcfn(cdfid1,cstnam)

  call cdfopn(trim(cstnam),cstid)
  
  !   Inquire the variables on the netCDF file:
  
  !   Inquire number of variables and variable names
  call getvars(cdfid1,nvars,vnam)
  !    print*,nvars,vnam(1),vnam(2)
  do i=1,nvars
     if (trim(vnam(i)).eq.'Q') qmode='Q'
     if (trim(vnam(i)).eq.'QD') qmode='QD'
     if (trim(vnam(i)).eq.'Z') ZonP=.true.
     if (trim(vnam(i)).eq.'T') TonP=.true.
     if (trim(vnam(i)).eq.'PS') PSonP=.true.
     if (trim(vnam(i)).eq.'U') UonP=.true.
     if (trim(vnam(i)).eq.'V') VonP=.true.
     if (trim(vnam(i)).eq.'ZB') ZBonP=.true.
     if (trim(vnam(i)).eq.'T2M') T2MonP=.true.
     if (trim(vnam(i)).eq.'TD2M') TD2MonP=.true.
     if (trim(vnam(i)).eq.'U10M') U10MonP=.true.
     if (trim(vnam(i)).eq.'V10M') V10MonP=.true.
     if (trim(vnam(i)).eq.'OMEGA') OMEGAonP=.true.
     if (trim(vnam(i)).eq.'PV') PVonP=.true.
     if (trim(vnam(i)).eq.'PV3') PV3onP=.true.
  enddo
  
  
  if (TonP) then
     call getdef(cdfid1,'T',ndim,mdv,vardim,varmin,varmax,stag)
  else
     call getdef(cdfid1,vnam(2),ndim,mdv,&
          vardim,varmin,varmax,stag)
  endif

  mdv=-999.98999
  
  !   Get the levels, pole, etc.
  
  nx=vardim(1)
  ny=vardim(2)
  nz=vardim(3)
  
  !    print*,nx,ny,nz
  call getgrid(cstid,dx,dy)
  call getlevs(cstid,nz,aklev,bklev,aklay,bklay)
  call getpole(cstid,pollon,pollat)
  call getstart(cstid,stdate)
  
  !   Allocate all arrays:
  !   --------------------
  allocate(sp1(nx,ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating sp1(nx,ny)'
  allocate(oro(nx,ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating oro(nx,ny)'
  allocate(cl(nx,ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating cl(nx,ny)'
  allocate(tl(nx,ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating tl(nx,ny)'
  allocate(f(nx,ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating f(nx,ny)'
  
  
  allocate(var(nx,ny,nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating var(nx,ny,nz)'


  allocate(tt(nx,ny,nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating tt(nx,ny,nz)'
  allocate(qq(nx,ny,nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating qq(nx,ny,nz)'
  allocate(uu(nx,ny,nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating uu(nx,ny,nz)'
  allocate(vv(nx,ny,nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating vv(nx,ny,nz)'
  allocate(ww(nx,ny,nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating ww(nx,ny,nz)'
  
  allocate(t2m(nx,ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating t2m(nx,ny)'
  allocate(td2m(nx,ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating td2m(nx,ny)'
  allocate(u10m(nx,ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating u10m(nx,ny)'
  allocate(v10m(nx,ny),STAT=stat)
  if (stat.ne.0) print*,'error allocating v10m(nx,ny)'
  allocate(ipv(nx,ny,nz),STAT=stat)
  if (stat.ne.0) print*,'error allocating ipv(nx,ny,nz)'
  
  !   Determine if data is on pressure or model levels
  
  prelev=.true.
  do k=1,nz
     if (bklev(k).ne.0.) prelev=.false.
     if (bklay(k).ne.0.) prelev=.false.
  enddo

  print*,'prelev ',prelev
  
  !   Calculate cos(latitude) array and the coriolis parameter
  
  if ((abs(pollon).gt.0.001).or.(abs(pollat-90.).gt.0.001)) then
     do j=1,ny
        rlat=varmin(2)+(j-1)*dy
        do i=1,nx
           rlon=varmin(1)+(i-1)*dx
           yphys=phstoph(rlat,rlon,pollat,pollon)
           !            if I use sind(lat in deg): troubles at the N-pole
           lat=2.*pi*yphys/360.
           cl(i,j)=cos(rlat)
           tl(i,j)=tan(lat)
           f(i,j)=0.000145444*sin(lat)
        enddo
     enddo
  else
     do j=1,ny
        lat=varmin(2)+(j-1)*dy
        lat=2.*pi*lat/360.
        do i=1,nx
           cl(i,j)=cos(lat)
           f(i,j)=0.000145444*sin(lat)
        enddo
     enddo
  endif
  
  !   Determine if data is on levels or layers
  
  if (stag(3).eq.-0.5) then
     ak=aklay
     bk=bklay
  else
     ak=aklev
     bk=bklev
  endif
  
  !   Get all the fields
  !   ------------------
  
  call gettimes(cdfid1,time,ntimes)
  
  !   Loop over all times
  !   -------------------
  
  !    print*,'loop over all times'
  !    print*,ntimes,vardim,nx,ny,nz
  
  do n=1,ntimes
     
     if (.not.prelev) then
        if (PSonP) then
           call getdatRank2(cdfid1,'PS',time,0,sp1)
        else
           if (verbose) print*,'PS not on P-file'
        endif
     endif
     
     if (Z_calc.eq.1) then                        
        allocate(zz(nx,ny,nz),STAT=stat)
        if (stat.ne.0) print*,'error allocating zz'       
        call zlev(zz,tt,qq,oro,sp1,nx,ny,nz,aklev,bklev)            
     endif
     
     if (Z_outP.eq.1) then
        if (n.eq.1) then
           stag(3)=0.0
           call putdef(cdfid1,'Z',4,mdv,&
                vardim,varmin,varmax,stag)
           write(*,*)'*** variable Z created on P-file'
           stag(3)=-0.5
        endif
        call putdat(cdfid1,'Z',time,0,zz)           
     endif

!   Create the secondary data file
     
     if (n.eq.1) then
        print *,outfnam
        call crecdf(trim(outfnam),cdfid,varmin,varmax,&
             3,trim(cstnam))
        write(*,*)'*** NetCDF file ',trim(outfnam),' created'
     endif

        !   Put surface pressure on S-file
     
     if ((.not.prelev).and.(PS_out.eq.1)) then
        vardim(3)=1
        if (n.eq.1) then
           call putdef(cdfid,'PS',4,mdv,vardim,varmin,varmax,stag)
           write(*,*)'*** variable PS created on ',trim(outfnam)
        endif
        call putdatRank2(cdfid,'PS',time,0,sp1)
        vardim(3)=nz
     endif
     
     !   Put geopotential on S-file
     
     if (Z_out.eq.1) then
        if (n.eq.1) then
           stag(3)=0.0
           call putdef(cdfid,'Z',4,mdv,&
                vardim,varmin,varmax,stag)
           write(*,*)'*** variable Z created on ',trim(outfnam)
           stag(3)=-0.5
        endif
        call putdat(cdfid,'Z',time,0,zz)
     endif
     
     !   Calculate the secondary data variables
     !   --------------------------------------
     
     !   Calculation of potential temperature
     
     if (TH_calc.eq.1) then
        allocate(th(nx,ny,nz),STAT=stat)
        if (stat.ne.0) print*,'error allocating th'
        call pottemp(th,tt,sp1,nx,ny,nz,ak,bk)
        if (TH_out.eq.1) then
           if (n.eq.1) then
              call putdef(cdfid,'TH',4,mdv,&
                   vardim,varmin,varmax,stag)
              write(*,*)'*** variable TH created on ',trim(outfnam)
           endif
           call putdat(cdfid,'TH',time,0,th)
        endif
     endif
         
     !   Calculation of potential vorticity (equals PV3 in IVE)
     
     if (PV_calc.eq.1) then
        allocate(pv(nx,ny,nz),STAT=stat)
        if (stat.ne.0) print*,'error allocating pv'
        call potvort(pv,uu,vv,th,sp1,cl,f,nx,ny,nz,ak,bk,varmin,varmax)
        if (PV_out.eq.1) then
           if (n.eq.1) then
              call putdef(cdfid,'PV',4,mdv,vardim,varmin,varmax,stag)
              write(*,*)'*** variable PV created on ',trim(outfnam)
           endif
           call putdat(cdfid,'PV',time,0,pv)
        endif
     endif
     
  end do
!   Close the NetCDF files
  
call clscdf(cdfid)

call clscdf(cdfid1)
call clscdf(cstid)

end program ptos


subroutine calc_rho(var,tt,sp,ie,je,ke,ak,bk)
  
  integer	ie,je,ke
  real,intent(OUT) :: var(ie,je,ke)
  real,intent(IN)  :: tt(ie,je,ke),sp(ie,je)
  real,intent(IN)  :: ak(ke),bk(ke)
  
  integer	i,j,k
  real,parameter ::	tzero=273.15
  
  
  real	pp,psrf
  integer	is
  pp(is)=ak(is)+bk(is)*psrf
  
  
  do i=1,ie
     do j=1,je
        psrf=sp(i,j)
        do k=1,ke
           
           if (tt(i,j,k).lt.100.) then
              var(i,j,k)=pp(k)/(2.87*(tt(i,j,k)+tzero))
           else
              var(i,j,k)=pp(k)/(2.87*tt(i,j,k))                 
           endif
            enddo
         enddo
      enddo
    end subroutine calc_rho
    

    subroutine nsq(var,rho,th,sp,ie,je,ke,ak,bk)
!   ======================================================
      
      
      integer,intent(IN) :: ie,je,ke
      real,intent(OUT)   :: var(ie,je,ke)
      real,intent(IN)    :: rho(ie,je,ke),th(ie,je,ke),sp(ie,je)
      real,intent(IN)    :: ak(ke),bk(ke)
      
      
      integer   stat
      REAL,ALLOCATABLE, DIMENSION (:,:,:) ::  dthdp !3D array
      
      allocate(dthdp(ie,je,ke),STAT=stat)
      IF (stat.ne.0) PRINT*,'nsq: error allocating dthdp(ie,je,ke)'
      
      call ddp(th,dthdp,sp,ie,je,ke,ak,bk)
      
      where(th.ne.0.) var=-96.24*rho/th*dthdp
      
      IF (ALLOCATED(dthdp)) DEALLOCATE(dthdp)
      
    end subroutine nsq
    
    subroutine dthetadp(var,th,sp,ie,je,ke,ak,bk)
      !   ======================================================
      
      
      integer,intent(IN) :: ie,je,ke
      real,intent(OUT)   :: var(ie,je,ke)
      real,intent(IN)    :: th(ie,je,ke),sp(ie,je)
      real,intent(IN)    :: ak(ke),bk(ke)
      

      integer   stat
      REAL,ALLOCATABLE, DIMENSION (:,:,:) ::  dthdp !3D array

      allocate(dthdp(ie,je,ke),STAT=stat)
      IF (stat.ne.0) PRINT*,'nsq: error allocating dthdp(ie,je,ke)'

      call ddp(th,dthdp,sp,ie,je,ke,ak,bk)

      var=-9.80616*100.*dthdp

      IF (ALLOCATED(dthdp)) DEALLOCATE(dthdp)

    end subroutine dthetadp

    subroutine geopot(psi,q,t,oro,sp,ie,je,ke,ak,bk)

      

      integer   ie,je,ke
      real      psi(ie,je,ke),t(ie,je,ke),q(ie,je,ke),oro(ie,je)
      real     sp(ie,je),ak(ke),bk(ke)
      

      integer   i,j,k
      real      r,c,g
      data	r,c,g /287.,0.608,9.80616/
      

      real      pp,psrf
      integer   is
      pp(is)=ak(is)+bk(is)*psrf
      

      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            psi(i,j,1)=1./g*(oro(i,j)&
                 +r*(t(i,j,1)+273.15)*(1.+c*q(i,j,1))*&
                 (psrf-pp(1))/(0.5*(psrf+pp(1))))
         enddo
      enddo
      do j=1,je
         do i=1,ie
            psrf=sp(i,j)
            do k=2,ke
               psi(i,j,k)=psi(i,j,k-1)+r/g*&
                   ((t(i,j,k-1)+273.15)*(1.+c*q(i,j,k-1))+&
                   (t(i,j,k)+273.15)*(1.+c*q(i,j,k)))*&
                   (pp(k-1)-pp(k))/(pp(k-1)+pp(k))
            enddo
         enddo
      enddo
      end
      



      subroutine calc_qx(var,uu,vv,tt,sp,cl,ie,je,ke,ak,bk,vmin,vmax)
!   ===============================================================



      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),tt(ie,je,ke)
      real,intent(IN)  :: sp(ie,je),cl(ie,je)
      real,intent(IN)  :: ak(ke),bk(ke),vmin(4),vmax(4)


      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dudx,dtdx,dvdx,dtdy
      integer	stat
      integer	k
      real	mdv
      logical	prelev


      prelev=.true.
      do k=1,ke
        if (bk(k).ne.0.) prelev=.false.
      enddo
      mdv=-999.98999

      if (.not.prelev) then
        allocate(dspdx(ie,je),STAT=stat)
        if (stat.ne.0) print*,'calc_qx: error allocating dspdx'
        allocate(dspdy(ie,je),STAT=stat)
        if (stat.ne.0) print*,'calc_qx: error allocating dspdy'
      endif

      allocate(dudx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qx: error allocating dudx'
      allocate(dtdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qx: error allocating dtdx'
      allocate(dvdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qx: error allocating dvdx'
      allocate(dtdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qx: error allocating dtdy'

      if (prelev) then
        do k=1,ke
          call ddh2m(uu(1,1,k),dudx(1,1,k),cl,'X',ie,je,1,vmin,vmax,mdv)
          call ddh2m(tt(1,1,k),dtdx(1,1,k),cl,'X',ie,je,1,vmin,vmax,mdv)
          call ddh2m(vv(1,1,k),dvdx(1,1,k),cl,'X',ie,je,1,vmin,vmax,mdv)
          call ddh2m(tt(1,1,k),dtdy(1,1,k),cl,'Y',ie,je,1,vmin,vmax,mdv)
        enddo
      else
        call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
        call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
        call ddh3(uu,dudx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(tt,dtdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(vv,dvdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(tt,dtdy,sp,dspdx,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
      endif

      var=-1.e9*9.80616/273.*(dudx*dtdx+dvdx*dtdy)
      if (prelev) then
        do i=1,ie
        do j=1,je
        do k=1,ke
          if ((dudx(i,j,k).eq.mdv).or.&
            (dtdx(i,j,k).eq.mdv).or.&
            (dvdx(i,j,k).eq.mdv).or.&
            (dtdy(i,j,k).eq.mdv)) then
            var(i,j,k)=mdv
          endif
        enddo
        enddo
        enddo
      endif

      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dudx)) DEALLOCATE(dudx)
      IF (ALLOCATED(dtdx)) DEALLOCATE(dtdx)
      IF (ALLOCATED(dvdx)) DEALLOCATE(dvdx)
      IF (ALLOCATED(dtdy)) DEALLOCATE(dtdy)

    end subroutine calc_qx
    

      subroutine calc_qy(var,uu,vv,tt,sp,cl,tl,ie,je,ke,ak,bk,vmin,vmax)
!   ===============================================================

      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: uu(ie,je,ke),vv(ie,je,ke),tt(ie,je,ke)
      real,intent(IN)  :: sp(ie,je),cl(ie,je),tl(ie,je)
      real,intent(IN)  :: ak(ke),bk(ke),vmin(4),vmax(4)


      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx,dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dudy,dtdx,dvdy,dtdy,tl3d
      integer	k,stat
      real      mdv
      logical   prelev
 

      prelev=.true.
      do k=1,ke
        if (bk(k).ne.0.) prelev=.false.
      enddo
      mdv=-999.98999
      
      if (.not.prelev) then
        allocate(dspdx(ie,je),STAT=stat)
        if (stat.ne.0) print*,'calc_qy: error allocating dspdx'
        allocate(dspdy(ie,je),STAT=stat)
        if (stat.ne.0) print*,'calc_qy: error allocating dspdy'
      endif

      allocate(dudy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating dudy'
      allocate(dtdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating dtdx'
      allocate(dvdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating dvdy'
      allocate(dtdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating dtdy'
      allocate(tl3d(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_qy: error allocating tl3d'

      if (prelev) then
        do k=1,ke
          call ddh2m(uu(1,1,k),dudy(1,1,k),cl,'Y',ie,je,1,vmin,vmax,mdv)
          call ddh2m(tt(1,1,k),dtdx(1,1,k),cl,'X',ie,je,1,vmin,vmax,mdv)
          call ddh2m(vv(1,1,k),dvdy(1,1,k),cl,'Y',ie,je,1,vmin,vmax,mdv)
          call ddh2m(tt(1,1,k),dtdy(1,1,k),cl,'Y',ie,je,1,vmin,vmax,mdv)
        enddo
      else
        call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
        call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
        call ddh3(uu,dudy,sp,dspdx,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(tt,dtdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(vv,dvdy,sp,dspdx,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
        call ddh3(tt,dtdy,sp,dspdx,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)
      endif


      do k=1,ke
         tl3d(1:ie,1:je,k)=tl(1:ie,1:je)
      enddo

      var=-1.e9*9.80616/273.*(dudy*dtdx+dvdy*dtdy&
         -uu/6.37E6*tl3d*dtdx-vv/6.37E6*tl3d*dtdy)
      if (prelev) then
        do i=1,ie
        do j=1,je
        do k=1,ke
          if ((dudy(i,j,k).eq.mdv).or.&
             (dtdx(i,j,k).eq.mdv).or.&
             (dvdy(i,j,k).eq.mdv).or.&
             (dtdy(i,j,k).eq.mdv)) then
            var(i,j,k)=mdv
          endif
        enddo
        enddo
        enddo
      endif

      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dudy)) DEALLOCATE(dudy)
      IF (ALLOCATED(dtdx)) DEALLOCATE(dtdx)
      IF (ALLOCATED(dvdy)) DEALLOCATE(dvdy)
      IF (ALLOCATED(dtdy)) DEALLOCATE(dtdy)
      IF (ALLOCATED(tl3d)) DEALLOCATE(tl3d)

      end


      subroutine calc_ug(var,zz,sp,cl,f,ie,je,ke,ak,bk,vmin,vmax)
!   ===========================================================
!   calculate geostrophic wind
      

      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: zz(ie,je,ke),sp(ie,je),cl(ie,je),f(ie,je)
      real	ak(ke),bk(ke),vmin(4),vmax(4)
      

      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdy
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dzzdy
      integer	k,stat
      
      allocate(dspdy(ie,je),STAT=stat)
      if (stat.ne.0) print*,'calc_ug: error allocating dspdy'
      allocate(dzzdy(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_ug: error allocating dzzdy'
      
      call ddh2(sp,dspdy,cl,'Y',ie,je,1,vmin,vmax)
      call ddh3(zz,dzzdy,sp,dspdy,cl,'Y',ie,je,ke,vmin,vmax,ak,bk)

      do k=1,ke
         var(1:ie,1:je,k)=-dzzdy(1:ie,1:je,k)&
              *9810.*(f(1:ie,1:je)/((ABS(f(1:ie,1:je))+1.e-12)**2))
      enddo  

      IF (ALLOCATED(dspdy)) DEALLOCATE(dspdy)
      IF (ALLOCATED(dzzdy)) DEALLOCATE(dzzdy)

      end


      subroutine calc_vg(var,zz,sp,cl,f,ie,je,ke,ak,bk,vmin,vmax)
!   ===========================================================
!   calculate geostrophic wind
      

      integer   ie,je,ke
      real,intent(OUT) :: var(ie,je,ke)
      real,intent(IN)  :: zz(ie,je,ke),sp(ie,je),cl(ie,je),f(ie,je)
      real	ak(ke),bk(ke),vmin(4),vmax(4)
      

      REAL,ALLOCATABLE, DIMENSION (:,:) :: dspdx
      REAL,ALLOCATABLE, DIMENSION (:,:,:) :: dzzdx
      integer	k,stat
      
      allocate(dspdx(ie,je),STAT=stat)
      if (stat.ne.0) print*,'calc_vg: error allocating dspdx'
      allocate(dzzdx(ie,je,ke),STAT=stat)
      if (stat.ne.0) print*,'calc_vg: error allocating dzzdx'
      
      call ddh2(sp,dspdx,cl,'X',ie,je,1,vmin,vmax)
      call ddh3(zz,dzzdx,sp,dspdx,cl,'X',ie,je,ke,vmin,vmax,ak,bk)

      do k=1,ke
         var(1:ie,1:je,k)=dzzdx(1:ie,1:je,k)&
             *9810.*(f(1:ie,1:je)/((ABS(f(1:ie,1:je))+1.e-12)**2))
      enddo 

      IF (ALLOCATED(dspdx)) DEALLOCATE(dspdx)
      IF (ALLOCATED(dzzdx)) DEALLOCATE(dzzdx)
     
      end

      subroutine zlayer(ap,t,z,sp,ie,je,ke,aklev,bklev,aklay,bklay)

      integer  ie,je,ke
      real,intent(OUT) :: ap(ie,je,ke)
      real,intent(IN)  :: t(ie,je,ke),z(ie,je,ke),sp(ie,je)
      real,intent(IN)  :: aklev(ke),bklev(ke),aklay(ke),bklay(ke)
      

      integer  i,j,k
      real     psrf
      real,parameter :: rdg=29.271,tzero=273.15
 

      real      prlev,prlay
      integer   is
      prlev(is)=aklev(is)+bklev(is)*psrf
      prlay(is)=aklay(is)+bklay(is)*psrf
 

      do i=1,ie
         do j=1,je
            psrf=sp(i,j)
            do k=1,ke-1
               ap(i,j,k) = (z(i,j,k) - rdg*(t(i,j,k)+tzero)*&
                    alog(prlay(k)/prlev(k)))/1000.
            enddo
            ap(i,j,ke) = (z(i,j,ke-1) + rdg*(t(i,j,ke)+tzero)*&
                alog(prlev(ke-1)/prlay(ke)))/1000.
         enddo
      enddo
      end      

      
      subroutine calc_psred(psr,ps,t,zb,ie,je,ke,aklay,bklay)


      integer  ie,je,ke
      real,intent(OUT) :: psr(ie,je)
      real,intent(IN)  :: ps(ie,je),t(ie,je,ke),zb(ie,je)
      real,intent(IN)  :: aklay(ke),bklay(ke)


      integer  i,j
      real     psrf
      real     ztstar,zalpha,zt0
      real,parameter :: rdcp=0.286,tzero=273.15,r=287.05,g=9.80665
     

      real      prlay
      integer   is
      prlay(is)=aklay(is)+bklay(is)*psrf

      do i=1,ie
         do j=1,je
            psrf=ps(i,j)
            if (zb(i,j).lt.1.) then
               psr(i,j)=psrf
            else
               ztstar = (t(i,j,1)+tzero)*(1. + 0.0065*r/g*&
                   (ps(i,j)/prlay(1) - 1.0))
               zalpha = 0.0065*r
               zt0    = ztstar + 0.0065*zb(i,j)
               if (zt0.gt.290.5) then
                  if (ztstar.gt.290.5) then
                     zalpha = 0.0
                     ztstar = 0.5*(ztstar+290.5)
                  else
                     zalpha = r*(290.5-ztstar)/zb(i,j)
                  endif
               else if (ztstar.lt.255.) then
                  ztstar = 0.5*(255.0+ztstar)
               endif
               psr(i,j) = ps(i,j)* exp(g*zb(i,j)/(r*ztstar)*&
                   (1.0 - 0.5*(zalpha*zb(i,j)/(r*ztstar)) +&
                   0.333*(zalpha*zb(i,j)/(r*ztstar))**2))
            endif
         enddo
      enddo
      end

      subroutine calc_blh(blh,ps,t,q,u,v,t2m,td2m,u10m,v10m,&
     			  ie,je,ke,aklay,bklay)

      integer  ie,je,ke
      real,intent(OUT) :: blh(ie,je)
      real,intent(IN)  :: ps(ie,je),t2m(ie,je),td2m(ie,je),&
     			  u10m(ie,je),v10m(ie,je)
      real,intent(IN)  :: t(ie,je,ke),q(ie,je,ke),u(ie,je,ke),&
     		          v(ie,je,ke)
      real,intent(IN)  :: aklay(ke),bklay(ke)
 

      integer  i,j
      real     psrf
      real     ztstar,zalpha,zt0
      real,parameter :: rdcp=0.286,tzero=273.15,r=287.05,g=9.80665
    

      real      prlay
      integer   is
      prlay(is)=aklay(is)+bklay(is)*psrf
 
      do i=1,ie
         do j=1,je
!           call richardson(ps,ust,t(i,j,1),q(i,j,1),u(i,j,1),v(i,j,1),&
 !                         ke,aklay,bklay,hf,t2m,td2m,blh(i,j),wst)
         enddo
      enddo
    end subroutine calc_blh

