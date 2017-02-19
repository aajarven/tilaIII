FUNCTION SASUM (N, SX, INCX)
!
!! SASUM compute the sum of the magnitudes of the elements of a vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A3A
!***TYPE      SINGLE PRECISION (SASUM-S, DASUM-D, SCASUM-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SUM OF MAGNITUDES OF A VECTOR
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
!        N  number of elements in input vector(S)
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!
!     --Output--
!    SASUM  single precision result (zero if N  <=  0)
!
!     Returns sum of magnitudes of single precision SX.
!     SASUM = sum from 0 to N-1 of ABS(SX(IX+I*INCX)),
!     where IX = 1 if INCX  >=  0, else IX = 1+(1-N)*INCX.
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
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SASUM
  real SASUM
  REAL SX(*)
  INTEGER I, INCX, IX, M, MP1, N
!***FIRST EXECUTABLE STATEMENT  SASUM
  SASUM = 0.0E0
  if (N  <=  0) RETURN
!
  if (INCX  ==  1) GOTO 20
!
!     Code for increment not equal to 1.
!
  IX = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  DO 10 I = 1,N
    SASUM = SASUM + ABS(SX(IX))
    IX = IX + INCX
   10 CONTINUE
  return
!
!     Code for increment equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 6.
!
   20 M = MOD(N,6)
  if (M  ==  0) GOTO 40
  DO 30 I = 1,M
    SASUM = SASUM + ABS(SX(I))
   30 CONTINUE
  if (N  <  6) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,6
    SASUM = SASUM + ABS(SX(I)) + ABS(SX(I+1)) + ABS(SX(I+2)) + &
            ABS(SX(I+3)) + ABS(SX(I+4)) + ABS(SX(I+5))
   50 CONTINUE
  return
end
