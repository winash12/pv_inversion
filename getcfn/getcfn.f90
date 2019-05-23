program getcst
!     ==============

!     Get the constants file name of a NetCDF file.
!     The program is invoked by the shell-script getcfn.
!     December 96	H. Wernli

  character*80 cdfnam,cstnam
  integer	cdfid,ierr
  integer	strend
  
  !     Read filename from input file
  
  read(*,'(a)')cdfnam
  
  !     Open file

  call cdfopn(cdfnam(1:strend(cdfnam)),cdfid,ierr)
  
  !     Ask about constants file name and write it's name
  
  call getcfn(cdfid,cstnam,ierr)
  write(*,'(a)') trim(cstnam)
  
end program getcst
