FUNCTION SDOT (N, SX, INCX, SY, INCY)
!
!! SDOT computes the inner product of two vectors.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      SINGLE PRECISION (SDOT-S, DDOT-D, CDOTU-C)
!***KEYWORDS  BLAS, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Lawson, C. L., (JPL)
!           Hanson, R. J., (SNLA)
!           Kincaid, D. R., (U. of Texas)
!           Krogh, F. T., (JPL)
!***DESCRIPTION
!
!                B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!       SY  single precision vector with N elements
!     INCY  storage spacing between elements of SY
!
!     --Output--
!     SDOT  single precision dot product (zero if N  <=  0)
!
!     Returns the dot product of single precision SX and SY.
!     SDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
!     where LX = 1 if INCX  >=  0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SDOT
  real SDOT
  REAL SX(*), SY(*)
!***FIRST EXECUTABLE STATEMENT  SDOT
  SDOT = 0.0E0
  if (N  <=  0) RETURN
  if (INCX  ==  INCY) IF (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
    5 IX = 1
  IY = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  if (INCY  <  0) IY = (-N+1)*INCY + 1
  DO 10 I = 1,N
    SDOT = SDOT + SX(IX)*SY(IY)
    IX = IX + INCX
    IY = IY + INCY
   10 CONTINUE
  return
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
   20 M = MOD(N,5)
  if (M  ==  0) go to 40
  DO 30 I = 1,M
    SDOT = SDOT + SX(I)*SY(I)
   30 CONTINUE
  if (N  <  5) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,5
  SDOT = SDOT + SX(I)*SY(I) + SX(I+1)*SY(I+1) + SX(I+2)*SY(I+2) + &
                SX(I+3)*SY(I+3) + SX(I+4)*SY(I+4)
   50 CONTINUE
  return
!
!     Code for equal, positive, non-unit increments.
!
   60 NS = N*INCX
  DO 70 I = 1,NS,INCX
    SDOT = SDOT + SX(I)*SY(I)
   70 CONTINUE
  return
end
