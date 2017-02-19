FUNCTION CDOTC (N, CX, INCX, CY, INCY)
!
!! CDOTC computes the dot product of two complex vectors using the complex ...
!            conjugate of the first vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A4
!***TYPE      COMPLEX (CDOTC-C)
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
!       CX  complex vector with N elements
!     INCX  storage spacing between elements of CX
!       CY  complex vector with N elements
!     INCY  storage spacing between elements of CY
!
!     --Output--
!    CDOTC  complex result (zero if N  <=  0)
!
!     Returns the dot product of complex CX and CY, using CONJUGATE(CX)
!     CDOTC = SUM for I = 0 to N-1 of CONJ(CX(LX+I*INCX))*CY(LY+I*INCY),
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
!***END PROLOGUE  CDOTC
  COMPLEX CDOTC
  COMPLEX CX(*),CY(*)
!***FIRST EXECUTABLE STATEMENT  CDOTC
  CDOTC = (0.0,0.0)
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
    CDOTC = CDOTC + CONJG(CX(KX))*CY(KY)
    KX = KX + INCX
    KY = KY + INCY
   10 CONTINUE
  return
!
!     Code for equal, positive increments.
!
   20 NS = N*INCX
  DO 30 I = 1,NS,INCX
  CDOTC = CDOTC + CONJG(CX(I))*CY(I)
   30 CONTINUE
  return
end
