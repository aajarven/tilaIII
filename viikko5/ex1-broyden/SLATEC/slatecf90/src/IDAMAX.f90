  INTEGER FUNCTION IDAMAX (N, DX, INCX)
!
!! IDAMAX finds the smallest index of that component of a vector ...
!            having the maximum magnitude.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A2
!***TYPE      DOUBLE PRECISION (ISAMAX-S, IDAMAX-D, ICAMAX-C)
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
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!   IDAMAX  smallest index (zero if N  <=  0)
!
!     Find smallest index of maximum magnitude of double precision DX.
!     IDAMAX = first I, I = 1 to N, to maximize ABS(DX(IX+(I-1)*INCX)),
!     where IX = 1 if INCX  >=  0, else IX = 1+(1-N)*INCX.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  IDAMAX
  DOUBLE PRECISION DX(*), DMAX, XMAG
  INTEGER I, INCX, IX, N
!***FIRST EXECUTABLE STATEMENT  IDAMAX
  IDAMAX = 0
  if (N  <=  0) RETURN
  IDAMAX = 1
  if (N  ==  1) RETURN
!
  if (INCX  ==  1) GOTO 20
!
!     Code for increments not equal to 1.
!
  IX = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  DMAX = ABS(DX(IX))
  IX = IX + INCX
  DO 10 I = 2,N
    XMAG = ABS(DX(IX))
    if (XMAG  >  DMAX) THEN
      IDAMAX = I
      DMAX = XMAG
    ENDIF
    IX = IX + INCX
   10 CONTINUE
  return
!
!     Code for increments equal to 1.
!
   20 DMAX = ABS(DX(1))
  DO 30 I = 2,N
    XMAG = ABS(DX(I))
    if (XMAG  >  DMAX) THEN
      IDAMAX = I
      DMAX = XMAG
    ENDIF
   30 CONTINUE
  return
end
