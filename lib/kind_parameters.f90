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
