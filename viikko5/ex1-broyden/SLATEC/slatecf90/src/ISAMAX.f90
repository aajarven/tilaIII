  INTEGER FUNCTION ISAMAX (N, SX, INCX)
!
!! ISAMAX finds the smallest index of that component of a vector
!  having the maximum magnitude.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A2
!***TYPE      SINGLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
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
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!
!     --Output--
!   ISAMAX  smallest index (zero if N  <=  0)
!
!     Find smallest index of maximum magnitude of single precision SX.
!     ISAMAX = first I, I = 1 to N, to maximize  ABS(SX(IX+(I-1)*INCX)),
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
!   920618  Slight restructuring of code.  (RWC, WRB)
!***END PROLOGUE  ISAMAX
  REAL SX(*), SMAX, XMAG
  INTEGER I, INCX, IX, N
!***FIRST EXECUTABLE STATEMENT  ISAMAX
  ISAMAX = 0
  if (N  <=  0) RETURN
  ISAMAX = 1
  if (N  ==  1) RETURN
!
  if (INCX  ==  1) GOTO 20
!
!     Code for increment not equal to 1.
!
  IX = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  SMAX = ABS(SX(IX))
  IX = IX + INCX
  DO 10 I = 2,N
    XMAG = ABS(SX(IX))
    if (XMAG  >  SMAX) THEN
      ISAMAX = I
      SMAX = XMAG
    ENDIF
    IX = IX + INCX
   10 CONTINUE
  return
!
!     Code for increments equal to 1.
!
   20 SMAX = ABS(SX(1))
  DO 30 I = 2,N
    XMAG = ABS(SX(I))
    if (XMAG  >  SMAX) THEN
      ISAMAX = I
      SMAX = XMAG
    ENDIF
   30 CONTINUE
  return
end
