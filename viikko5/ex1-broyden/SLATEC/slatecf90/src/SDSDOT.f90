FUNCTION SDSDOT (N, SB, SX, INCX, SY, INCY)
!
!! SDSDOT computes the inner product of two vectors with extended precision.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      SINGLE PRECISION (SDSDOT-S, CDCDOT-C)
!***KEYWORDS  BLAS, DOT PRODUCT, INNER PRODUCT, LINEAR ALGEBRA, VECTOR
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
!       SB  single precision scalar to be added to inner product
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!       SY  single precision vector with N elements
!     INCY  storage spacing between elements of SY
!
!     --Output--
!   SDSDOT  single precision dot product (SB if N  <=  0)
!
!     Returns S.P. result with dot product accumulated in D.P.
!     SDSDOT = SB + sum for I = 0 to N-1 of SX(LX+I*INCX)*SY(LY+I*INCY),
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
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SDSDOT
  real SDSDOT
  REAL SX(*), SY(*), SB
  DOUBLE PRECISION DSDOT
!***FIRST EXECUTABLE STATEMENT  SDSDOT
  DSDOT = SB
  if (N  <=  0) go to 30
  if (INCX == INCY .AND. INCX > 0) go to 40
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
   30 SDSDOT = DSDOT
  return
!
!     Code for equal and positive increments.
!
   40 NS = N*INCX
  DO 50 I = 1,NS,INCX
    DSDOT = DSDOT + DBLE(SX(I))*DBLE(SY(I))
   50 CONTINUE
  SDSDOT = DSDOT
  return
end
