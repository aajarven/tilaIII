DOUBLE PRECISION FUNCTION DASUM (N, DX, INCX)
!
!! DASUM sums the magnitudes of the elements of a vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A3A
!***TYPE      DOUBLE PRECISION (SASUM-S, DASUM-D, SCASUM-C)
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
!        N  number of elements in input vector(s)
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!
!     --Output--
!    DASUM  double precision result (zero if N  <=  0)
!
!     Returns sum of magnitudes of double precision DX.
!     DASUM = sum from 0 to N-1 of ABS(DX(IX+I*INCX)),
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
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900821  Modified to correct problem with a negative increment.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DASUM
  DOUBLE PRECISION DX(*)
  INTEGER I, INCX, IX, M, MP1, N
!***FIRST EXECUTABLE STATEMENT  DASUM
  DASUM = 0.0D0
  if (N  <=  0) RETURN
!
  if (INCX  ==  1) GOTO 20
!
!     Code for increment not equal to 1.
!
  IX = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  DO 10 I = 1,N
    DASUM = DASUM + ABS(DX(IX))
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
    DASUM = DASUM + ABS(DX(I))
   30 CONTINUE
  if (N  <  6) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,6
    DASUM = DASUM + ABS(DX(I)) + ABS(DX(I+1)) + ABS(DX(I+2)) + &
            ABS(DX(I+3)) + ABS(DX(I+4)) + ABS(DX(I+5))
   50 CONTINUE
  return
end
