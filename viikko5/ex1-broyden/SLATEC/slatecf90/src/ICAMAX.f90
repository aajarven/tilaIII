  INTEGER FUNCTION ICAMAX (N, CX, INCX)
!
!! ICAMAX finds the smallest index of the component of a complex ...
!            vector having the maximum sum of magnitudes of real
!            and imaginary parts.
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A2
!***TYPE      COMPLEX (ISAMAX-S, IDAMAX-D, ICAMAX-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, MAXIMUM COMPONENT, VECTOR
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
!
!     --Output--
!   ICAMAX  smallest index (zero if N  <=  0)
!
!     Returns the smallest index of the component of CX having the
!     largest sum of magnitudes of real and imaginary parts.
!     ICAMAX = first I, I = 1 to N, to maximize
!     ABS(REAL(CX(IX+(I-1)*INCX))) + ABS(IMAG(CX(IX+(I-1)*INCX))),
!     where IX = 1 if INCX  >=  0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  ICAMAX
  COMPLEX CX(*)
  REAL SMAX, XMAG
  INTEGER I, INCX, IX, N
  COMPLEX ZDUM
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
!***FIRST EXECUTABLE STATEMENT  ICAMAX
  ICAMAX = 0
  if (N  <=  0) RETURN
  ICAMAX = 1
  if (N  ==  1) RETURN
!
  if (INCX  ==  1) GOTO 20
!
!     Code for increment not equal to 1.
!
  IX = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  SMAX = CABS1(CX(IX))
  IX = IX + INCX
  DO 10 I = 2,N
    XMAG = CABS1(CX(IX))
    if (XMAG  >  SMAX) THEN
      ICAMAX = I
      SMAX = XMAG
    ENDIF
    IX = IX + INCX
   10 CONTINUE
  return
!
!     Code for increment equal to 1.
!
   20 SMAX = CABS1(CX(1))
  DO 30 I = 2,N
    XMAG = CABS1(CX(I))
    if (XMAG  >  SMAX) THEN
      ICAMAX = I
      SMAX = XMAG
    ENDIF
   30 CONTINUE
  return
end
