      subroutine wricst(cstnam,datar,aklev,bklev,aklay,bklay,stdate)
C------------------------------------------------------------------------

C     Creates the constants file for NetCDF files containing ECMWF
C     data. The constants file is compatible with the one created
C     for EM data (with subroutine writecst).
C
C     Input parameters:
C
C     cstnam    name of constants file
C     datar     array contains all required parameters to write file
C               datar(1):       number of points along x        
C               datar(2):       number of points along y
C               datar(3):       maximum latitude of data region (ymax)
C               datar(4):       minimum longitude of data region (xmin)
C               datar(5):       minimum latitude of data region (ymin)
C               datar(6):       maximum longitude of data region (xmax)
C               datar(7):       grid increment along x
C               datar(8):       grid increment along y
C               datar(9):       number of levels        
C		datar(10):	data type (forecast or analysis)
C		datar(11):	data version
C		datar(12):	constants file version
C		datar(13):	longitude of pole of coordinate system
C		datar(14):	latitude of pole of coordinate system
C     aklev     array contains the aklev values
C     bklev	array contains the bklev values
C     aklay     array contains the aklay values
C     bklay     array contains the bklay values
C     stdate    array contains date (year,month,day,time,step) of first
C               field on file (start-date), dimensionised as stdate(5)
C------------------------------------------------------------------------


      include "netcdf.inc"

      integer   nchar,maxlev

      parameter (nchar=20,maxlev=32)
      real	aklev(maxlev),bklev(maxlev)
      real      aklay(maxlev),bklay(maxlev)
      real	pollat,latmin,latmax
      integer   datar(14)
      integer	stdate(5)
      character*80 cstnam

C     declarations for constants-variables

      integer   nz
      integer   dattyp, datver, cstver

C     further declarations

      integer	ierr			! error flag
      integer	cdfid			! NetCDF id
      integer	xid,yid,zid		! dimension ids
      integer	pollonid, pollatid,	! variable ids
     >		aklevid, bklevid, aklayid, bklayid,
     >		lonminid, lonmaxid, latminid, latmaxid,
     >		dellonid, dellatid,
     >		startyid, startmid, startdid, starthid, startsid,
     >		dattypid, datverid, cstverid

      nz=datar(9)			! number of levels

C     Set data-type and -version, version of cst-file-format

      dattyp=datar(10)
      datver=datar(11)
      cstver=datar(12)

C     Initially set error to false

      ierr=0

C     Create constants file

      cdfid=nccre(trim(cstnam),NCCLOB,ierr)

C     Define the dimensions

      xid = ncddef (cdfid,'nx',datar(1),ierr)
      yid = ncddef (cdfid,'ny',datar(2),ierr)
      zid = ncddef (cdfid,'nz',datar(9),ierr)

C     Define integer constants

      pollonid = ncvdef(cdfid,'pollon', NCFLOAT,0,0,ierr)
      pollatid = ncvdef(cdfid,'pollat', NCFLOAT,0,0,ierr)

      aklevid = ncvdef (cdfid, 'aklev', NCFLOAT, 1, zid, ierr)
      bklevid = ncvdef (cdfid, 'bklev', NCFLOAT, 1, zid, ierr)
      aklayid = ncvdef (cdfid, 'aklay', NCFLOAT, 1, zid, ierr)
      bklayid = ncvdef (cdfid, 'bklay', NCFLOAT, 1, zid, ierr)

      lonminid = ncvdef (cdfid, 'lonmin', NCFLOAT, 0, 0, ierr)
      lonmaxid = ncvdef (cdfid, 'lonmax', NCFLOAT, 0, 0, ierr)
      latminid = ncvdef (cdfid, 'latmin', NCFLOAT, 0, 0, ierr)
      latmaxid = ncvdef (cdfid, 'latmax', NCFLOAT, 0, 0, ierr)
      dellonid = ncvdef (cdfid, 'dellon', NCFLOAT, 0, 0, ierr)
      dellatid = ncvdef (cdfid, 'dellat', NCFLOAT, 0, 0, ierr)
      startyid = ncvdef (cdfid, 'starty', NCLONG, 0, 0, ierr)
      startmid = ncvdef (cdfid, 'startm', NCLONG, 0, 0, ierr)
      startdid = ncvdef (cdfid, 'startd', NCLONG, 0, 0, ierr)
      starthid = ncvdef (cdfid, 'starth', NCLONG, 0, 0, ierr)
      startsid = ncvdef (cdfid, 'starts', NCLONG, 0, 0, ierr)
      dattypid = ncvdef (cdfid, 'dattyp', NCLONG, 0, 0, ierr)
      datverid = ncvdef (cdfid, 'datver', NCLONG, 0, 0, ierr)
      cstverid = ncvdef (cdfid, 'cstver', NCLONG, 0, 0, ierr)

C     Leave define mode

      call ncendf(cdfid,ierr)

C     Store levels
      call ncvpt(cdfid, aklevid, 1, nz, aklev, ierr)
      call ncvpt(cdfid, bklevid, 1, nz, bklev, ierr)
      call ncvpt(cdfid, aklayid, 1, nz, aklay, ierr)
      call ncvpt(cdfid, bklayid, 1, nz, bklay, ierr)

C     Store position of pole (trivial for ECMWF data)
      call ncvpt1(cdfid, pollonid, 1, real(datar(13))/1000., ierr)
      if (datar(14).gt.0) then
        pollat=min(real(datar(14))/1000.,90.)
      else
        pollat=max(real(datar(14))/1000.,-90.)
      endif
      call ncvpt1(cdfid, pollatid, 1, pollat, ierr)

C     Store horizontal data borders and grid increments
      call ncvpt1(cdfid, lonminid, 1, real(datar(4))/1000., ierr)
      call ncvpt1(cdfid, lonmaxid, 1, real(datar(6))/1000., ierr)
      latmin=max(real(datar(5))/1000.,-90.)
      latmax=min(real(datar(3))/1000.,90.)
      call ncvpt1(cdfid, latminid, 1, latmin, ierr)
      call ncvpt1(cdfid, latmaxid, 1, latmax, ierr)
      call ncvpt1(cdfid, dellonid, 1, real(datar(7))/1000., ierr)
      call ncvpt1(cdfid, dellatid, 1, real(datar(8))/1000., ierr)

C     Store date of first field on file (start-date)
      call ncvpt1(cdfid, startyid, 1, stdate(1), ierr)
      call ncvpt1(cdfid, startmid, 1, stdate(2), ierr)
      call ncvpt1(cdfid, startdid, 1, stdate(3), ierr)
      call ncvpt1(cdfid, starthid, 1, stdate(4), ierr)
      call ncvpt1(cdfid, startsid, 1, stdate(5), ierr)

C     Store datatype and version
      call ncvpt1(cdfid, dattypid, 1, dattyp, ierr)
      call ncvpt1(cdfid, datverid, 1, datver, ierr)

C     Store version of the constants file format
      call ncvpt1(cdfid, cstverid, 1, cstver, ierr)

C     Store strings

      call ncclos(cdfid,ierr)
      return

      end
      subroutine writelmcst(cdfid,nx,ny,nz,pollon,pollat,lonmin,
     &lonmax,latmin,latmax,dellon,dellat,dattyp,datver,cstver,
     &psref,tstar,tbeta,pintf,p0top,idate)
c     ------------------------------------------------------------------

      implicit none

      integer   cdfid

c     deklarationen der constants-variablen
      real       pollon,pollat
      real       lonmin,lonmax,latmin,latmax,dellon,dellat
      integer    idate(5)
      integer    nx,ny,nz
      integer    dattyp, datver, cstver
      real       psref, tstar, tbeta, pintf, p0top

      include 'netcdf.inc'

* netcdf declaration
      integer   iret, k
* dimension ids
      integer  nxdim, nydim, nzdim
* variable ids
      integer  startyid, startmid, startdid, starthid
* variable shapes, corners and edge lengths
      integer dims(1), corner(1), edges(1)

* enter define mode
      call ncredf(cdfid, iret)

      startyid = ncvdef (cdfid, 'starty', NCLONG, 0, 0, iret)
      startmid = ncvdef (cdfid, 'startm', NCLONG, 0, 0, iret)
      startdid = ncvdef (cdfid, 'startd', NCLONG, 0, 0, iret)
      starthid = ncvdef (cdfid, 'starth', NCLONG, 0, 0, iret)

* store the rest as global attributes
* store nx,ny,nz
      call ncapt(cdfid,NCGLOBAL,'nx',NCLONG,1,nx,iret)
      call ncapt(cdfid,NCGLOBAL,'ny',NCLONG,1,ny,iret)
      call ncapt(cdfid,NCGLOBAL,'nz',NCLONG,1,nz,iret)

* store pollon, pollat
      call ncapt(cdfid,NCGLOBAL,'pollon',NCFLOAT,1,pollon,iret)
      call ncapt(cdfid,NCGLOBAL,'pollat',NCFLOAT,1,pollat,iret)

* store lonmin, etc
      call ncapt(cdfid,NCGLOBAL,'lonmin',NCFLOAT,1,lonmin,iret)
      call ncapt(cdfid,NCGLOBAL,'lonmax',NCFLOAT,1,lonmax,iret)
      call ncapt(cdfid,NCGLOBAL,'latmin',NCFLOAT,1,latmin,iret)
      call ncapt(cdfid,NCGLOBAL,'latmax',NCFLOAT,1,latmax,iret)
      call ncapt(cdfid,NCGLOBAL,'dellon',NCFLOAT,1,dellon,iret)
      call ncapt(cdfid,NCGLOBAL,'dellat',NCFLOAT,1,dellat,iret)

* store data type and version
      call ncapt(cdfid,NCGLOBAL,'dattyp',NCLONG,1,dattyp,iret)
      call ncapt(cdfid,NCGLOBAL,'datver',NCLONG,1,datver,iret)
      call ncapt(cdfid,NCGLOBAL,'cstver',NCLONG,1,cstver,iret)

* store information of lm model vertical grid
      call ncapt(cdfid,NCGLOBAL,'psref',NCFLOAT,1,psref,iret)
      call ncapt(cdfid,NCGLOBAL,'tstar',NCFLOAT,1,tstar,iret)
      call ncapt(cdfid,NCGLOBAL,'tbeta',NCFLOAT,1,tbeta,iret)
      call ncapt(cdfid,NCGLOBAL,'pintf',NCFLOAT,1,pintf,iret)
      call ncapt(cdfid,NCGLOBAL,'p0top',NCFLOAT,1,p0top,iret)

* leave define mode
      call ncendf(cdfid, iret)

* store starty, etc
      corner(1) = 1
      edges(1) = 1
      call ncvpt(cdfid, startyid, corner, edges, idate(1), iret)
      call ncvpt(cdfid, startmid, corner, edges, idate(2), iret)
      call ncvpt(cdfid, startdid, corner, edges, idate(3), iret)
      call ncvpt(cdfid, starthid, corner, edges, idate(4), iret)

      end
      subroutine globcst(cdfnam,datar,aklev,bklev,aklay,bklay,stdate)
C------------------------------------------------------------------------
C+
C NAME:
C     subroutine globcst
C
C PURPOSE:
C     instead of writing a constants-file (*_cst), the information
C     is added to the netCDF file as global variables
C     the data format is compatible with the one requested by
C     the IVE ETH/MIT version, contact author about details
C
C CATEGORY:
C     model,netCDF
C
C CALLING SEQUENCE:
C     subroutine globcst(cdfnam,datar,aklev,bklev,aklay,bklay,stdate)
C
C INPUTS:
C     cdfnam    name of netCDF file
C               The file needs to exist, otherwise an ERROR occurs,
C               i.e. nothing is done
C     datar     array contains all required parameters to write file
C               datar(1):       number of points along x
C               datar(2):       number of points along y
C               datar(3):       maximum latitude of data region (ymax)
C               datar(4):       minimum longitude of data region (xmin)
C               datar(5):       minimum latitude of data region (ymin)
C               datar(6):       maximum longitude of data region (xmax)
C               datar(7):       grid increment along x
C               datar(8):       grid increment along y
C               datar(9):       number of levels
C               datar(10):      data type (forecast or analysis)
C               datar(11):      data version
C               datar(12):      constants file version
C               datar(13):      longitude of pole of coordinate system
C               datar(14):      latitude of pole of coordinate system
C     aklev     array contains the aklev values
C     bklev     array contains the bklev values
C     aklay     array contains the aklay values
C     bklay     array contains the bklay values
C     stdate    array contains date (year,month,day,time,step) of first
C               field on file (start-date), dimensionised as stdate(5)
C     list    the griblist-ASCII-file
C     varno   the GRIB code number
C
C OUTPUTS:
C     Adds cdf-information to EXISTING netCDF-file
C
C MODIFICATION HISTORY:
C
C     June  93    Christoph Schaer (ETHZ) created
C     Nov   93    Heini Wernli (ETHZ) wricst
C     Nov   98    David N. Bresch (MIT) wricst to globcst
C-
 
C     Sun include statement.
      include "netcdf.inc"
 
      integer   nchar,maxlev
 
      parameter (nchar=20,maxlev=32)
      real      aklev(maxlev),bklev(maxlev)
      real      aklay(maxlev),bklay(maxlev)
      integer   datar(14)
      integer   stdate(5)
      character*80 cdfnam
 
C     declarations for constants-variables
 
      integer   nz
      integer   dattyp, datver, cstver
 
C     further declarations
 
      integer   ierr                    ! error flag
      integer   cdfid                   ! NetCDF id
      integer   xid,yid,zid             ! dimension ids
      integer   pollonid, pollatid,     ! variable ids
     >          aklevid, bklevid, aklayid, bklayid,
     >          lonminid, lonmaxid, latminid, latmaxid,
     >          dellonid, dellatid,
     >          startyid, startmid, startdid, starthid, startsid,
     >          dattypid, datverid, cstverid
 
      nz=datar(9)                       ! number of levels
 
C     Set data-type and -version, version of cst-file-format
 
      dattyp=datar(10)
      datver=datar(11)
      cstver=datar(12)
 
C     Initially set error to false
 
      ierr=0
 
C     open the netCDF-file:
 
      call cdfwopn(cdfnam,cdfid,ierr)
      if (ierr.ne.0) then
         print*,'ERROR opening netCDF-file ',cdfnam
         return
      endif
 
C     Put file into define mode
      call ncredf(cdfid,ierr)
      if (ierr.ne.0) then
         print*,'ERROR switching to netCDF redefine mode'
         return
      endif
 
C     Define the dimensions
 
      xid = ncddef (cdfid,'nx',datar(1),ierr)
      yid = ncddef (cdfid,'ny',datar(2),ierr)
      zid = ncddef (cdfid,'nz',datar(9),ierr)
 
C     Define integer constants
 
      pollonid = ncvdef(cdfid,'pollon', NCFLOAT,0,0,ierr)
      pollatid = ncvdef(cdfid,'pollat', NCFLOAT,0,0,ierr)
 
      aklevid = ncvdef (cdfid, 'aklev', NCFLOAT, 1, zid, ierr)
      bklevid = ncvdef (cdfid, 'bklev', NCFLOAT, 1, zid, ierr)
      aklayid = ncvdef (cdfid, 'aklay', NCFLOAT, 1, zid, ierr)
      bklayid = ncvdef (cdfid, 'bklay', NCFLOAT, 1, zid, ierr)
 
      lonminid = ncvdef (cdfid, 'lonmin', NCFLOAT, 0, 0, ierr)
      lonmaxid = ncvdef (cdfid, 'lonmax', NCFLOAT, 0, 0, ierr)
      latminid = ncvdef (cdfid, 'latmin', NCFLOAT, 0, 0, ierr)
      latmaxid = ncvdef (cdfid, 'latmax', NCFLOAT, 0, 0, ierr)
      dellonid = ncvdef (cdfid, 'dellon', NCFLOAT, 0, 0, ierr)
      dellatid = ncvdef (cdfid, 'dellat', NCFLOAT, 0, 0, ierr)
      startyid = ncvdef (cdfid, 'starty', NCLONG, 0, 0, ierr)
      startmid = ncvdef (cdfid, 'startm', NCLONG, 0, 0, ierr)
      startdid = ncvdef (cdfid, 'startd', NCLONG, 0, 0, ierr)
      starthid = ncvdef (cdfid, 'starth', NCLONG, 0, 0, ierr)
      startsid = ncvdef (cdfid, 'starts', NCLONG, 0, 0, ierr)
      dattypid = ncvdef (cdfid, 'dattyp', NCLONG, 0, 0, ierr)
      datverid = ncvdef (cdfid, 'datver', NCLONG, 0, 0, ierr)
      cstverid = ncvdef (cdfid, 'cstver', NCLONG, 0, 0, ierr)
 
C     Leave define mode
 
      call ncendf(cdfid,ierr)
      if (ierr.ne.0) then
         print*,'ERROR exiting define mode'
         return
      endif
 
C     Store levels
      call ncvpt(cdfid, aklevid, 1, nz, aklev, ierr)
      call ncvpt(cdfid, bklevid, 1, nz, bklev, ierr)
      call ncvpt(cdfid, aklayid, 1, nz, aklay, ierr)
      call ncvpt(cdfid, bklayid, 1, nz, bklay, ierr)
 
C     Store position of pole (trivial for ECMWF data)
      call ncvpt1(cdfid, pollonid, 1, real(datar(13))/1000., ierr)
      call ncvpt1(cdfid, pollatid, 1, real(datar(14))/1000., ierr)
 
C     Store horizontal data borders and grid increments
      call ncvpt1(cdfid, lonminid, 1, real(datar(4))/1000., ierr)
      call ncvpt1(cdfid, lonmaxid, 1, real(datar(6))/1000., ierr)
      call ncvpt1(cdfid, latminid, 1, real(datar(5))/1000., ierr)
      call ncvpt1(cdfid, latmaxid, 1, real(datar(3))/1000., ierr)
      call ncvpt1(cdfid, dellonid, 1, real(datar(7))/1000., ierr)
      call ncvpt1(cdfid, dellatid, 1, real(datar(8))/1000., ierr)
 
C     Store date of first field on file (start-date)
      call ncvpt1(cdfid, startyid, 1, stdate(1), ierr)
      call ncvpt1(cdfid, startmid, 1, stdate(2), ierr)
      call ncvpt1(cdfid, startdid, 1, stdate(3), ierr)
      call ncvpt1(cdfid, starthid, 1, stdate(4), ierr)
      call ncvpt1(cdfid, startsid, 1, stdate(5), ierr)
 
C     Store datatype and version
      call ncvpt1(cdfid, dattypid, 1, dattyp, ierr)
      call ncvpt1(cdfid, datverid, 1, datver, ierr)
 
C     Store version of the constants file format
      call ncvpt1(cdfid, cstverid, 1, cstver, ierr)
 
      if (ierr.ne.0) then
         print*,'ERROR adding cst-date as global variables'
         return
      endif
 
C     Store strings
 
      call ncclos(cdfid,ierr)
      if (ierr.ne.0) then
         print*,'ERROR closing netCDF file'
      endif
 
      return
      end

      subroutine getsdat(cdfid,varnam,time,ix,iy,iz,sx,sy,sz,dat,error)
c-----------------------------------------------------------------------
c     Purpose:
c        This routine is called to read the data within a selected
c	 domain of a variable from an IVE-NetCDF file.
c        Prior to calling this routine, the file must be opened with
c        a call to opncdf (for extension) or crecdf (for creation) or
c        readcdf (for readonly).
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (must be obtained by calling routine
c                              opncdf,readcdf  or crecdf)
c        varnam  char  input   the user-supplied variable name
c        time    real  input   the user-supplied time-level of the
c                              data to be read from the file (the time-
c                              levels stored in the file can be obtained
c                              with a call to gettimes).
c        ix/y/z  int   input   indices of lower left corner of selected
c			       data volume.
c	 sx/y/z  int   input   size of selected data volume
c        dat     real  output  data-array with dimensions (sx,sy,sz).
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not present on
c                                          the file.
c                              error = 2   the value of 'time' is not
c                                          known.to the file.
c			       error = 6,7,8   data volume too large
c                              error =10   another error.
c     History:
c       June  93    Christoph Schaer (ETHZ)  Created getdat
c	Nov   93    Heini Wernli (ETHZ)	     Created getsdat
c-----------------------------------------------------------------------

      include "netcdf.inc"

C     Declaration of local variables
      character*(*) varnam
      character*(20) chars
      integer cdfid

      integer     ix,iy,iz,sx,sy,sz
      real        dat(sx,sy,sz)
      real        misdat,varmin(4),varmax(4),stag(4)
      real        time, timeval

      integer     corner(4),edgeln(4),didtim,vardim(4),ndims
      integer     error, ierr
      integer     ntime
      integer     idtime,idvar,iflag
      integer     i

      call ncpopt(NCVERBOS)

c     access the variable
      call getdef (cdfid, trim(varnam), ndims, misdat,
     &                           vardim, varmin, varmax, stag, ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* in getdef in getdat'
        error=1
        return
      endif
      idvar=ncvid(cdfid,trim(varnam),ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* in ncvid in getsdat'
        error=1
        return
      endif

C     Get times-array
      didtim=ncdid(cdfid,'time',ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* didtim in getsdat'
        error=10
        return
      endif
      call ncdinq(cdfid,didtim,chars,ntime,ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* in ncdinq in getsdat'
        error=10
        return
      endif
      idtime=ncvid(cdfid,'time',ierr)
      if (ierr.ne.0) then
        print *,'*ERROR* in ncvid for time in getsdat'
        error=10
        return
      endif
c     find appropriate time-index
      iflag=0
      do i=1,ntime
        call ncvgt1(cdfid,idtime,i,timeval,ierr)
        if (ierr.ne.0) print *,'*ERROR* in ncvgt1 in getsdat'
        if (time.eq.timeval) iflag=i
      enddo
      if (iflag.eq.0) then
        error=2
        print *,'Error: Unknown time in getsdat'
        print *,time,timeval
        return
      endif

C     Define data volume to be written (index space)
      corner(1)=ix
      corner(2)=iy
      corner(3)=iz
      corner(4)=iflag
      edgeln(1)=sx
      edgeln(2)=sy
      edgeln(3)=sz
      edgeln(4)=1

C     Check if data volume is within data domain

      if (ix+sx-1.gt.vardim(1)) then
        error=7
        print *,'Error: data volume too large in x-direction'
        print *,ix,sx,vardim(1)
        return
      endif
      if (iy+sy-1.gt.vardim(2)) then
        error=8
        print *,'Error: data volume too large in y-direction'
        return
      endif
      if (iz+sz-1.gt.vardim(3)) then
        error=9
        print *,'Error: data volume too large in z-direction'
        return
      endif

C     Read data from NetCDF file

      call ncvgt(cdfid,idvar,corner,edgeln,dat,error)
      if (error.ne.0) then
        print *, 'corner ',corner(1),corner(2),corner(3)
        print *, 'edgeln ',edgeln(1),edgeln(2),edgeln(3)
        print *, '*ERROR* in ncvgt in getsdat'
        error=10
      endif
      end
      subroutine getlevs(cstid,nlev,aklev,bklev,aklay,bklay,error)
c-----------------------------------------------------------------------
c     Purpose:
c     	This routine is called to get the level arrays aklev and
c	bklev from a NetCDF constants file.
c     Arguments:
c	cstid     int	input   identifier for NetCDF constants file
c	nlev	  int	input	number of levels
c	aklev     real	output  array contains all aklev values
c       bklev     real  output  array contains all bklev values
c	aklay	  real  output	array contains all aklay values
c	bklay	  real	output	array contains all bklay values
c	error	  int	output	error flag
c				error = 0   no errors detected
c				error = 1   error detected
c     History:
c	Aug. 93	  Heini Wernli		Created.
c-----------------------------------------------------------------------

      integer   error

      integer   cstid
      integer   ncdid,ncvid		! NetCDF functions
      integer   didz,idak,idbk,idaky,idbky
      integer   nlev
      real      aklev(nlev),bklev(nlev),aklay(nlev),bklay(nlev)
      character*(20) dimnam
      integer   i

      didz	=ncdid(cstid,'nz',error)
      if (error.ne.0) goto 920
      idak	=ncvid(cstid,'aklev',error)
      if (error.ne.0) goto 920
      idbk	=ncvid(cstid,'bklev',error)
      if (error.ne.0) goto 920
      idaky     =ncvid(cstid,'aklay',error)
      if (error.ne.0) goto 920
      idbky     =ncvid(cstid,'bklay',error)
      if (error.ne.0) goto 920

      call ncdinq(cstid,didz,dimnam,nlev,error)	! read number of levels
      if (error.ne.0) goto 920

      do 10 i=1,nlev
        call ncvgt1(cstid,idak,i,aklev(i),error)      ! get aklev
        call ncvgt1(cstid,idbk,i,bklev(i),error)      ! get bklev
        call ncvgt1(cstid,idaky,i,aklay(i),error)      ! get aklay
        call ncvgt1(cstid,idbky,i,bklay(i),error)      ! get bklay
        if (error.ne.0) goto 920
   10 continue

      return

c     Error exits.
  920 write(*,*)'*ERROR*: An error occured in subroutine getlevs'
      return

      end
      subroutine getntim(cdfid,ntimes,ierr)
C------------------------------------------------------------------------
C     Purpose:
C        Get number of times on the specified NetCDF file
C     Arguments:
C        cdfid  int  input   identifier for NetCDF file
C        ntimes int  output  number of times on the file
C        error  int  output  errorflag
C     History:
C        Heini Wernli, ETHZ
C------------------------------------------------------------------------
 
      include "netcdf.inc"
 
      integer   ierr
      integer didtim,ntimes
 
      integer   cdfid,idtime
      integer   ncopts
      character*(20) dimnam
 
c     Get current value of error options, and make sure netCDF-errors do
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt(NCVERBOS)
 
      didtim=ncdid(cdfid,'time',ierr)   ! inquire id for time dimension
      if (ierr.ne.0) goto 900
      idtime=ncvid(cdfid,'time',ierr)   ! inquire id for time array
      if (ierr.ne.0) goto 900
      call ncdinq(cdfid,didtim,dimnam,ntimes,ierr)      ! inquire # of times
      if (ierr.ne.0) goto 900
 
c     normal exit
      call ncpopt (ncopts)
      return
 
c     error exit
 900  ntimes=1
      call ncpopt (ncopts)
      end
      subroutine getstart(cdfid,idate,ierr)
C------------------------------------------------------------------------
C     Purpose:
C	Get start date for fields on specified NetCDF file
C     Arguments:
C	cdfid	int	input	identifier for NetCDF file
C	idate	int	output	array contains date (year,month,day,time,step)
C				dimensioned as idate(5)
C	ierr	int	output	error flag
C------------------------------------------------------------------------

      include "netcdf.inc"

c     variable declarations
      integer   ierr
      integer   idate(5)
      integer   cdfid,ncopts,idvar,nvars
      integer   ndims,ngatts,recdim,i,vartyp,nvatts,vardim(4)
      character*20 vnam(100)

c     Get current value of error options, and make sure NetCDF-errors do
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt (NCVERBOS)

      idvar=ncvid(cdfid,'starty',ierr)
      if (ierr.ne.0) goto 930
      call ncvgt1(cdfid,idvar,1,idate(1),ierr)
      if (ierr.ne.0) goto 920

      idvar=ncvid(cdfid,'startm',ierr)
      if (ierr.ne.0) goto 920
      call ncvgt1(cdfid,idvar,1,idate(2),ierr)
      if (ierr.ne.0) goto 920

      idvar=ncvid(cdfid,'startd',ierr)
      if (ierr.ne.0) goto920
      call ncvgt1(cdfid,idvar,1,idate(3),ierr)
      if (ierr.ne.0) goto 920

      idvar=ncvid(cdfid,'starth',ierr)
      if (ierr.ne.0) goto 920
      call ncvgt1(cdfid,idvar,1,idate(4),ierr)
      if (ierr.ne.0) goto 920

C     Starts is not defined on all files
C     Only ask for it if it exists
C     Inquire number of dimensions, variables and attributes
 
      idate(5)=0
      call ncinq(cdfid,ndims,nvars,ngatts,recdim,ierr)
      do i=1,nvars
        call ncvinq(cdfid,i,vnam(i),vartyp,ndims,vardim,nvatts,ierr)
        if (vnam(i).eq.'starts') then
          idvar=ncvid(cdfid,'starts',ierr)
          call ncvgt1(cdfid,idvar,1,idate(5),ierr)
          if (ierr.ne.0) goto 920
        endif
      enddo

c     normal exit
      call ncpopt (ncopts)
      return

c     error exit
 920  continue
      write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'read the starting-time in subroutine putstart.'
 930  continue
      call ncpopt (ncopts)

      end
      subroutine putstart(cdfid,idate,ierr)
C----------------------------------------------------------------------
C     Purpose:
C        Puts the 'starting-time' on the specified NetCDF file.
C     Arguments:
C        cdfid   int     input   identifier for NetCDF file
C        idate   int     input   array contains date (year,month,day,time,step)
C                                dimensioned as idate(5)
C        ierr    int     output  error flag
C------------------------------------------------------------------------

      include "netcdf.inc"

c     variable declarations
      integer   ierr,idate(5),startid(5),cdfid,ncopts,i

c     Get current value of error options, and make sure NetCDF-errors do
c     not abort execution
      call ncgopt (ncopts)
      call ncpopt (NCVERBOS)

c     define variables
      call ncredf(cdfid,ierr)
      if (ierr.ne.0) goto 920
      startid(1) = ncvdef (cdfid, 'starty', NCLONG, 0, 0, ierr)
      if (ierr.ne.0) goto 920
      startid(2) = ncvdef (cdfid, 'startm', NCLONG, 0, 0, ierr)
      if (ierr.ne.0) goto 920
      startid(3) = ncvdef (cdfid, 'startd', NCLONG, 0, 0, ierr)
      if (ierr.ne.0) goto 920
      startid(4) = ncvdef (cdfid, 'starth', NCLONG, 0, 0, ierr)
      if (ierr.ne.0) goto 920
      startid(5) = ncvdef (cdfid, 'starts', NCLONG, 0, 0, ierr)
      if (ierr.ne.0) goto 920
      call ncendf(cdfid, ierr)
      if (ierr.ne.0) goto 920

c     store variables
      do i=1,5
        call ncvpt1(cdfid,startid(i),1,idate(i),ierr)
        if (ierr.ne.0) goto 920
      enddo

c     synchronyse output to disk, revert to previous error-mode, and exit
      call ncsnc (cdfid,ierr)
      call ncpopt (ncopts)
      return

c     error exit
 920  write (6, *) 'ERROR: An error occurred while attempting to ',
     &             'write the starting-time in subroutine putstart.'
      call ncpopt (ncopts)
      call ncclos (cdfid, ierr)

      end
      subroutine getgrid(cdfid,dx,dy,ierr)
C------------------------------------------------------------------------
C     Purpose:
C       Get grid increments for fields on specified NetCDF file
C     Arguments:
C       cdfid   int     input   identifier for NetCDF file
C	dx	real	output	grid increment along latitude
C	dy	real	output	grid increment along longitude
C       ierr    int     output  error flag
C------------------------------------------------------------------------

      integer   ierr

      integer   cdfid
      integer   ncvid

      integer   idilon,idilat
      real	dx,dy

      idilon    =ncvid(cdfid,'dellon',ierr)
      if (ierr.ne.0) return
      idilat    =ncvid(cdfid,'dellat',ierr)
      if (ierr.ne.0) return

      call ncvgt1(cdfid,idilon,1,dx,ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,idilat,1,dy,ierr)
      if (ierr.ne.0) return

      end
      subroutine getdattyp(cdfid,typ,ierr)
C------------------------------------------------------------------------
C     Purpose:
C       Get data type for specified NetCDF file
C     Arguments:
C       cdfid   int     input   identifier for NetCDF file
C       typ     int     output  data type: 1 (52) for pressure (theta) coord
C       ierr    int     output  error flag
C------------------------------------------------------------------------
 
      integer   ierr
 
      integer   cdfid
      integer   ncvid
 
      integer   idtyp,typ
 
      idtyp    =ncvid(cdfid,'dattyp',ierr)
      if (ierr.ne.0) return
 
      call ncvgt1(cdfid,idtyp,1,typ,ierr)
      if (ierr.ne.0) return
 
      end
      subroutine getpole(cdfid,pollon,pollat,ierr)
C------------------------------------------------------------------------
C     Purpose:
C       Get physical coordinates of pole of coordinate system
C     Arguments:
C       cdfid   int     input   identifier for NetCDF file
C	pollon	real	output	longitude of pole
C	pollat	real	output	latitude of pole
C       ierr    int     output  error flag
C------------------------------------------------------------------------

      integer   ierr

      integer   cdfid
      integer   ncvid

      integer   idplon,idplat
      real      pollon,pollat

      idplon    =ncvid(cdfid,'pollon',ierr)
      if (ierr.ne.0) return
      idplat    =ncvid(cdfid,'pollat',ierr)
      if (ierr.ne.0) return

      call ncvgt1(cdfid,idplon,1,pollon,ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,idplat,1,pollat,ierr)
      if (ierr.ne.0) return

      end
      subroutine getmc2grid(cdfid,polx,poly,delx,shem,phi0,lam0,ierr)
C------------------------------------------------------------------------
C     Purpose:
C       Get physical coordinates of pole of coordinate system
C     Arguments:
C       cdfid   int     input   identifier for NetCDF file
C       ierr    int     output  error flag
C------------------------------------------------------------------------
 
      integer   ierr
 
      integer   cdfid
      integer   ncvid
 
      integer   idpolx,idpoly,iddelx,idshem,idphi0,idlam0
      real      polx,poly,delx,shem,phi0,lam0
 
      idpolx    =ncvid(cdfid,'polx',ierr)
      if (ierr.ne.0) return
      idpoly    =ncvid(cdfid,'poly',ierr)
      if (ierr.ne.0) return
      iddelx	=ncvid(cdfid,'delx',ierr)
      if (ierr.ne.0) return
      idshem	=ncvid(cdfid,'shem',ierr)
      if (ierr.ne.0) return
      idphi0	=ncvid(cdfid,'phi0',ierr)
      if (ierr.ne.0) return
      idlam0	=ncvid(cdfid,'lam0',ierr)
      if (ierr.ne.0) return
 
      call ncvgt1(cdfid,idpolx,1,polx,ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,idpoly,1,poly,ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,iddelx,1,delx,ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,idshem,1,shem,ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,idphi0,1,phi0,ierr)
      if (ierr.ne.0) return
      call ncvgt1(cdfid,idlam0,1,lam0,ierr)
      if (ierr.ne.0) return
 
      end


      subroutine getcfn(cdfid,cfn,ierr)
C------------------------------------------------------------------------
C     Purpose:
C       Get name of constants file
C     Arguments:
C       cdfid   int     input   identifier for NetCDF file
C       cfn     char    output  name of constants file
C       ierr    int     output  error flag
C------------------------------------------------------------------------

      include 'netcdf.inc'

      integer   ierr
      integer   cdfid,lenstr
      character*80 cfn

      lenstr=80
      call ncagtc(cdfid,NCGLOBAL,"constants_file_name",cfn,lenstr,ierr)
      if (ierr.ne.0) write(*,*)'error in SR getcfn'

      end


      subroutine getdim (cdfid, varnam, nx, ny, nz, error)
c-------------------------------------------------------------------------
c     Purpose:
c        This routine is called to get the dimensions of
c        a variable from an IVE-NetCDF file for use with the IVE plotting
c        package. Prior to calling this routine, the file must be opened
c        with a call to opncdf.
c     Arguments:
c        cdfid   int   input   file-identifier
c                              (can be obtained by calling routine
c                              opncdf)
c        varnam  char  input   the user-supplied variable name.
c                              (can be obtained by calling routine
c                              opncdf)
c        nx      int   output  the zonal dimension of the variable.
c        ny      int   output  the meridional dimension of the variable.
c        nz      int   output  the vertical dimension of the variable.
c
c        error   int   output  indicates possible errors found in this
c                              routine.
c                              error = 0   no errors detected.
c                              error = 1   the variable is not on the file.
c                              error =10   other errors.
c     History:
c        March 2000    Heini Wernli (ETHZ)     Created.
c-----------------------------------------------------------------------

      include  "netcdf.inc"

c     Argument declarations.
      character *(*) varnam
      integer        vardim(4), ndim, error, cdfid
      integer        nx,ny,nz

c     Local variable declarations.
      character *(20) dimnam(MAXNCDIM),vnam
      integer         id,i,k
      integer         ndims,nvars,ngatts,recdim,dimsiz(MAXNCDIM)
      integer         vartyp,nvatts, ncopts

c     Get current value of error options.
      call ncgopt (ncopts)

c     make sure NetCDF-errors do not abort execution
      call ncpopt(NCVERBOS)

c     Initially set error to indicate no errors.
      error = 0

c     inquire for number of dimensions
      call ncinq(cdfid,ndims,nvars,ngatts,recdim,error)
      if (error.eq.1) goto 920

c     read dimension-table
      do i=1,ndims
        call ncdinq(cdfid,i,dimnam(i),dimsiz(i),error)
        if (error.gt.0) goto 920
      enddo

c     get id of the variable
      id=ncvid(cdfid,varnam,error)
      if (error.eq.1) goto 910

c     inquire about variable
      call ncvinq(cdfid,id,vnam,vartyp,ndim,vardim,nvatts,error)
      if (vartyp.ne.NCFLOAT) error=1
      if (error.gt.0) goto 920

c     get dimensions from dimension-table
      do k=1,ndim
        vardim(k)=dimsiz(vardim(k))
      enddo

      nx=vardim(1)
      ny=vardim(2)
      nz=vardim(3)

c     normal exit
      call ncpopt (ncopts)
      return

c     Error exits.
 910  write (6, *) '*ERROR*: The selected variable could not be found ',
     &             'in the file by getdim.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      return

 920  write (6, *) '*ERROR*: An error occurred while attempting to ',
     &             'read the data file in subroutine getcdf.'
      call ncpopt (ncopts)
      call ncclos (cdfid, error)
      end


      subroutine new_gettra(cdfid,varnam,ix,ntimes,array,ierr)
C------------------------------------------------------------------------
C
C     Reads the time-evolution for one grid-point of the variable
C     indicated by varnam.
C
C     cdfid     int     input   identifier for NetCDF file
C     varnam    char    input   name of variable
C     ix        int     input   index for trajectory to read
C     ntimes    int     input   number of time-indices to read
C     array     real    output  array contains the readed values
C     ierr      int     output  error flag
C------------------------------------------------------------------------
 
C     Declaration of attributes
 
      integer   cdfid
      character*(*) varnam
      integer   ix
      integer   ntimes
      real      array(ntimes)
 
C     Declaration of local variables
 
      integer   corner(4),edgeln(4)
      integer   idvar,ierr
      integer   ncvid
      integer	strend
 
      corner(1)=ix
      corner(2)=1
      corner(3)=1
      corner(4)=1
      edgeln(1)=1
      edgeln(2)=1
      edgeln(3)=1
      edgeln(4)=ntimes
 
      idvar =ncvid(cdfid,varnam(1:strend(varnam)),ierr)
      call ncvgt(cdfid,idvar,corner,edgeln,array,ierr)
      if (ierr.ne.0) goto 991
 
      return
  991 stop 'Variable not found on NetCDF file in SR new_gettra'
      end
      subroutine getvars(cdfid,nvars,vnam,ierr)
C------------------------------------------------------------------------
 
C     Opens the NetCDF file 'filnam' and returns its identifier cdfid.
 
C     filnam    char    input   name of NetCDF file to open
C     nvars     int     output  number of variables on file
C     vnam	char	output  array with variable names
C     ierr      int     output  error flag
C------------------------------------------------------------------------

      include 'netcdf.inc'
 
      integer   cdfid,ierr,nvars
      character*(*) vnam(*)

      integer	ndims,ngatts,recdim,i,vartyp,nvatts,vardim(4)
 
      call ncpopt(NCVERBOS)

C     Inquire number of dimensions, variables and attributes
 
      call ncinq(cdfid,ndims,nvars,ngatts,recdim,ierr)
 
C     Inquire variable names from NetCDF file
 
      do i=1,nvars
        call ncvinq(cdfid,i,vnam(i),vartyp,ndims,vardim,nvatts,ierr)
      enddo
 
      return
      end
