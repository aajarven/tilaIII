  DOUBLE PRECISION FUNCTION DSDOT (N, SX, INCX, SY, INCY)
!
!! DSDOT computes the inner product of two vectors with extended ...
!            precision accumulation and result.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      DOUBLE PRECISION (DSDOT-D, DCDOT-C)
!***KEYWORDS  BLAS, COMPLEX VECTORS, DOT PRODUCT, INNER PRODUCT,
!             LINEAR ALGEBRA, VECTOR
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
!    DSDOT  double precision dot product (zero if N <= 0)
!
!     Returns D.P. dot product accumulated in D.P., for S.P. SX and SY
!     DSDOT = sum for I = 0 to N-1 of  SX(LX+I*INCX) * SY(LY+I*INCY),
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
!***END PROLOGUE  DSDOT
  REAL SX(*),SY(*)
!***FIRST EXECUTABLE STATEMENT  DSDOT
  DSDOT = 0.0D0
  if (N  <=  0) RETURN
  if (INCX == INCY .AND. INCX > 0) go to 20
!
!     Code for unequal or nonpositive increments.
!
  KX = 1
  KY = 1
  if (INCX  <  0) KX = 1+(1-N)*INCX
  if (INCY  <  0) KY = 1+(1-N)*INCY
  DO 10 I = 1,N
    DSDOT = DSDOT + DBLE(SX(KX))*DBLE(SY(KY))
    KX = KX + INCX
    KY = KY + INCY
   10 CONTINUE
  return
!
!     Code for equal, positive, non-unit increments.
!
   20 NS = N*INCX
  DO 30 I = 1,NS,INCX
    DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
   30 CONTINUE
  return
end
