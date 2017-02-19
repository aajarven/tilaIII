function scnrm2 ( n, x, incx )

!*******************************************************************************
!
!! SCNRM2 returns the euclidean norm of a complex vector.
!
!  Discussion:
!
!    SCNRM2 := sqrt ( sum ( conjg ( x(1:n) ) * x(1:n) ) )
!            = sqrt ( dot_product ( x(1:n), x(1:n) ) )
!
!  Reference:
!
!    Lawson, Hanson, Kincaid and Krogh,
!    Basic Linear Algebra Subprograms for FORTRAN usage,
!    ACM Transactions on Mathematical Software,
!    Volume 5, Number 3, pages 308-323, 1979.
!
!  Parameters:
!
!    Input, integer N, the number of entries in the vector.
!
!    Input, complex X(*), the vector.
!
!    Input, integer INCX, the increment between successive entries of X.
!
!    Output, real SCNRM2, the norm of the vector.
!
  implicit none

  integer incx
  integer ix
  integer n
  real norm
  real, parameter :: one = 1.0E+00
  real scale
  real scnrm2
  real ssq
  real temp
  complex x(*)
  real, parameter :: zero = 0.0E+00

  if ( n < 1 .or. incx < 1 ) then

    norm  = zero

  else

    scale = zero
    ssq = one

    do ix = 1, 1 + ( n - 1 ) * incx, incx

      if ( real ( x(ix) ) /= zero ) then
        temp = abs ( real( x(ix) ) )
        if ( scale < temp ) then
          ssq = one + ssq * ( scale / temp )**2
          scale = temp
        else
          ssq = ssq + ( temp / scale )**2
        end if
      end if

      if ( aimag ( x(ix) ) /= zero ) then
        temp = abs ( aimag ( x(ix) ) )
        if ( scale < temp ) then
          ssq = one + ssq * ( scale / temp )**2
          scale = temp
        else
          ssq = ssq + ( temp / scale )**2
        end if

      end if

    end do

    norm  = scale * sqrt ( ssq )

  end if

  scnrm2 = norm

  return
end
