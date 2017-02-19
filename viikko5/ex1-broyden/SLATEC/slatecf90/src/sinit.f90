subroutine sinit ( n, sa, x, incx )
!
!*******************************************************************************
!
!! SINIT initializes a real vector to a constant.
!
!
!  Modified:
!
!    08 April 1999
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, real SA, the constant value to be used to initialize X.
!
!    Output, real X(*), the vector to be initialized.
!
!    Input, integer INCX, the increment between successive entries of X.
!
  implicit none
!
  integer i
  integer incx
  integer ix
  integer n
  real sa
  real x(*)
!
  if ( n <= 0 ) then

  else if ( incx == 1 ) then

    x(1:n) = sa

  else

    if ( incx >= 0 ) then
      ix = 1
    else
      ix = ( - n + 1 ) * incx + 1
    end if

    do i = 1, n
      x(ix) = sa
      ix = ix + incx
    end do

  end if

  return
end
