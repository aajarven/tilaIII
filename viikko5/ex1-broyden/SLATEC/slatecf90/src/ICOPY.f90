subroutine ICOPY (N, IX, INCX, IY, INCY)
!
!! ICOPY copies a vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A5
!***TYPE      INTEGER (ICOPY-S, DCOPY-D, CCOPY-C, ICOPY-I)
!***KEYWORDS  BLAS, COPY, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Boland, W. Robert, (LANL)
!           Clemens, Reginald, (PLK)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       IX  integer vector with N elements
!     INCX  storage spacing between elements of IX
!       IY  integer vector with N elements
!     INCY  storage spacing between elements of IY
!
!     --Output--
!       IY  copy of vector IX (unchanged if N  <=  0)
!
!     Copy integer IX to integer IY.
!     For I = 0 to N-1, copy  IX(LX+I*INCX) to IY(LY+I*INCY),
!     where LX = 1 if INCX  >=  0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   930201  DATE WRITTEN
!***END PROLOGUE  ICOPY
  INTEGER IX(*), IY(*)
!***FIRST EXECUTABLE STATEMENT  ICOPY
  if (N  <=  0) RETURN
  if (INCX  ==  INCY) IF (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
    5 IIX = 1
  IIY = 1
  if (INCX  <  0) IIX = (-N+1)*INCX + 1
  if (INCY  <  0) IIY = (-N+1)*INCY + 1
  DO 10 I = 1,N
    IY(IIY) = IX(IIX)
    IIX = IIX + INCX
    IIY = IIY + INCY
   10 CONTINUE
  return
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 7.
!
   20 M = MOD(N,7)
  if (M  ==  0) go to 40
  DO 30 I = 1,M
    IY(I) = IX(I)
   30 CONTINUE
  if (N  <  7) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,7
    IY(I) = IX(I)
    IY(I+1) = IX(I+1)
    IY(I+2) = IX(I+2)
    IY(I+3) = IX(I+3)
    IY(I+4) = IX(I+4)
    IY(I+5) = IX(I+5)
    IY(I+6) = IX(I+6)
   50 CONTINUE
  return
!
!     Code for equal, positive, non-unit increments.
!
   60 NS = N*INCX
  DO 70 I = 1,NS,INCX
    IY(I) = IX(I)
   70 CONTINUE
  return
end
