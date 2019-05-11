! Module to determine kinds for different precisions
!------------------------------------------------------------------------------

MODULE kind_parameters

!==============================================================================
!>
!!  Module determines kinds for different precisions.
!!  Number model from which the SELECTED_*\\_KIND are requested: <br>
!! @f{tabular}{{r@{\hspace*{3em}}c@{\hspace*{3em}}c}
!!                     &4 byte REAL     &8 byte REAL        \\\
!!        IEEE:        &precision = 6   &precision =   15   \\\
!!                     &exponent  = 37  &exponent  =  307
!! @f}
!! \\medskip
!!
!!
!! Current Code Owner
!! @author  Ulrich Schaettler, DWD
!!  email:  ulrich.schaettler@dwd.de
!!
!! @par Revision History
!
! History:
! Version      Date       Name
! ----------   ---------- ----
! V5_1         2014-11-28 Ulrich Schaettler
!  Initial Release
!
!!
!! @par Copyright and License
!!
!! This code is subject to the DWD and MPI-M-Software-License-Agreement in
!! its most recent form.
!! Please see the file LICENSE in the root of the source tree for this code.
!! Where software is supplied by third parties, it is indicated in the
!! headers of the routines.
!!
!==============================================================================

  IMPLICIT NONE

  PRIVATE

  PUBLIC  :: sp, dp,i1, i2, i4, i8

!==============================================================================

  ! Floating point section
  ! ----------------------

  INTEGER, PARAMETER :: sp = SELECTED_REAL_KIND( 6, 37) !< single precision
  INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(12,307) !< double precision


  INTEGER, PARAMETER :: wp = sp          !< selected working precision is single precision


  ! Integer section
  ! ---------------

  INTEGER, PARAMETER :: i1 = SELECTED_INT_KIND(  2)   !< at least 1 byte integer
  INTEGER, PARAMETER :: i2 = SELECTED_INT_KIND(  4)   !< at least 2 byte integer
  INTEGER, PARAMETER :: i4 = SELECTED_INT_KIND(  9)   !< at least 4 byte integer
  INTEGER, PARAMETER :: i8 = SELECTED_INT_KIND( 18)   !< at least 8 byte integer

  INTEGER, PARAMETER :: wi = i4                       !< selected working precision

!==============================================================================

END MODULE kind_parameters
