subroutine dinit ( n, sa, x, incx )
!
!*******************************************************************************
!
!! DINIT initializes a double precision vector to a constant.
!
!
!  Modified:
!
!    28 October 2002
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, double precision SA, the constant value to be used to initialize X.
!
!    Output, double precision X(*), the vector to be initialized.
!
!    Input, integer INCX, the increment between successive entries of X.
!
  implicit none
!
  integer i
  integer incx
  integer ix
  integer n
  double precision sa
  double precision x(*)
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
