PROGRAM difference

!  **************************************************************************
!  * Get the difference between two runs                                    *
!  * Michael Sprenger / Autumn 2006                                         *
!  **************************************************************************

  use kind_parameters,ONLY:&
       sp

  implicit none
  
  !  --------------------------------------------------------------------------
  !  Declaration of subroutine parameters
  !  --------------------------------------------------------------------------
  
  !  Physical and numerical parameters
  real(kind=sp)      eps
  parameter (eps=0.01)
  integer   nmax
  parameter (nmax=300*300*200)

  !  Variables for input file 1
  character*80 i1_filename
  real(kind=sp)         i1_varmin(3),i1_varmax(3),i1_stag(3)
  integer      i1_vardim(3)
  real(kind=sp)         i1_mdv
  integer      i1_ndim
  integer      i1_nx,i1_ny,i1_nz
  real(kind=sp)         i1_time
  integer      i1_nvars
  character*80 i1_vnam(100)
  real(kind=sp)         i1_field(nmax)
  
  !  Variables for input file 2
  character*80 i2_filename
  real(kind=sp)         i2_varmin(3),i2_varmax(3),i2_stag(3)
  integer      i2_vardim(3)
  real(kind=sp)         i2_mdv
  integer      i2_ndim
  integer      i2_nx,i2_ny,i2_nz
  real(kind=sp)         i2_time
  integer      i2_nvars
  character*80 i2_vnam(100)
  real(kind=sp)         i2_field(nmax)
  
  !  Variables for output file
  character*80 o_filename
  real(kind=sp)         o_varmin(3),o_varmax(3),o_stag(3)
  integer      o_vardim(3)
  real(kind=sp)         o_mdv
  integer      o_ndim
  integer      o_nx,o_ny,o_nz
  real(kind=sp)         o_time
  real(kind=sp)         o_field(nmax)
  
  !  Auxiliary variables
  integer      i,j,k
  integer      cdfid,ierr
  integer      isok
  character*80 cfn,varname
  
  !  --------------------------------------------------------------------------------
  !  Input
  !  --------------------------------------------------------------------------------
  
  print*,'********************************************************'
  print*,'* DIFFERENCE                                           *'
  print*,'********************************************************'
  
  !  Read parameter file
  open(10,file='fort.10')
  read(10,*) i1_filename
  read(10,*) i2_filename
  read(10,*) o_filename
  close(10)
  print*
  print*,trim(i1_filename)
  print*,trim(i2_filename)
  print*,trim(o_filename)
  print*
  
  !  Get a list of all variables on the two input files
  call cdfopn(i1_filename,cdfid,ierr)

  call getvars(cdfid,i1_nvars,i1_vnam)

  call getdef(cdfid,i1_vnam(i1_nvars),i1_ndim,i1_mdv,i1_vardim,i1_varmin,i1_varmax,i1_stag)
  call clscdf(cdfid)
  print*,(trim(i1_vnam(i))//'  ',i=1,i1_nvars)

  call cdfopn(i2_filename,cdfid)

  call getvars(cdfid,i2_nvars,i2_vnam)

  call clscdf(cdfid)
  print*,(trim(i2_vnam(i))//'  ',i=1,i2_nvars)
  print*

!  Create the output file
  o_ndim=i1_ndim
  do i=1,i1_ndim
     o_varmin(i)=i1_varmin(i)
     o_varmax(i)=i1_varmax(i)
  enddo
  cfn=trim(o_filename)//'_cst'
  call crecdf(trim(o_filename),cdfid,o_varmin,o_varmax,o_ndim,trim(cfn))

  call clscdf(cdfid)

!  --------------------------------------------------------------------------------
!  Loop through all variables
!  --------------------------------------------------------------------------------

  do i=1,i1_nvars
     
     !     Check wether the variable is available on both files
     varname=i1_vnam(i)
     if (varname.eq.'time') goto 100
     isok=0
     call check_varok (isok,varname,i2_vnam,i2_nvars)
     if (isok.eq.0) goto 100
     
     !     Read first input file
     call cdfopn(i1_filename,cdfid)

     call getdef(cdfid,varname,i1_ndim,i1_mdv,i1_vardim,&
          i1_varmin,i1_varmax,i1_stag)
     
     call getdat(cdfid,varname,0.,0,i1_field)

     call clscdf(cdfid)
         
!     Read second input file
     call cdfopn(i2_filename,cdfid)
     call getdef(cdfid,varname,i2_ndim,i2_mdv,i2_vardim,i2_varmin,i2_varmax,i2_stag)
     call getdat(cdfid,varname,0.,0,i2_field)
     
     call clscdf(cdfid)

!     Consistency check
     if (i1_ndim.ne.i2_ndim) then
        print*,'Inconsistent input files... Stop',i1_ndim,i2_ndim
        stop
     endif
     do j=1,3
        if ( (i1_vardim(j).ne.i2_vardim(j)).or.&
             (abs(i1_varmin(j)-i2_varmin(j)).gt.eps).or.&
             (abs(i1_varmax(j)-i2_varmax(j)).gt.eps).or.&
             (abs(i1_stag(j)-i2_stag(j)).gt.eps)) then
           print*,'Inconsistent input files... Stop'
           print*,j,i1_varmin(j),i2_varmin(j)
           print*,j,i1_varmax(j),i2_varmax(j)
           print*,j,i1_vardim(j),i2_vardim(j)
           print*,j,i1_stag(j),  i2_stag(j)
           stop
        endif
     enddo
     
!     Get the difference
     do j=1,i1_vardim(1)*i1_vardim(2)*i1_vardim(3)
        if ( (abs(i1_field(j)-i1_mdv).gt.eps).and.
           >           (abs(i2_field(j)-i2_mdv).gt.eps) ) then
           o_field(j)=i1_field(j)-i2_field(j)
        else
           o_field(j)=i1_mdv
        endif
     enddo
     
!     Write to output file
     o_ndim=i1_ndim
     o_mdv =i1_mdv
     do j=1,i1_ndim
        o_vardim(j)=i1_vardim(j)
        o_varmin(j)=i1_varmin(j)
        o_varmax(j)=i1_varmax(j)
        o_stag(j)  =i1_stag(j)
     enddo
     call cdfwopn(o_filename,cdfid)
     call putdef(cdfid,varname,o_ndim,o_mdv,o_vardim,&
          o_varmin,o_varmax,o_stag)

     call putdat(cdfid,varname,0.,0,o_field)
     print*,'W ',trim(varname),' ',trim(o_filename)
     call clscdf(cdfid)
     
!     Next 
100  continue
     
  enddo
  
  !  -----------------------------------------------------------------
!  Exception handling and format specs
!  -----------------------------------------------------------------

  stop
  
998 print*,'Problems with input file 1  ',trim(i1_filename)
  stop
  
997 print*,'Problems with input file 2  ',trim(i2_filename)
  stop

996 print*,'Problems with output file   ',trim(o_filename)
  stop
  
end PROGRAM difference


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
  
end subroutine check_varok
      
