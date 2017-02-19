subroutine SSCAL (N, SA, SX, INCX)
!
!! SSCAL multiplies a vector by a constant.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A6
!***TYPE      SINGLE PRECISION (SSCAL-S, DSCAL-D, CSCAL-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, SCALE, VECTOR
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
!       SA  single precision scale factor
!       SX  single precision vector with N elements
!     INCX  storage spacing between elements of SX
!
!     --Output--
!       SX  single precision result (unchanged if N  <=  0)
!
!     Replace single precision SX by single precision SA*SX.
!     For I = 0 to N-1, replace SX(IX+I*INCX) with  SA * SX(IX+I*INCX),
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
!***END PROLOGUE  SSCAL
  REAL SA, SX(*)
  INTEGER I, INCX, IX, M, MP1, N
!***FIRST EXECUTABLE STATEMENT  SSCAL
  if (N  <=  0) RETURN
  if (INCX  ==  1) GOTO 20
!
!     Code for increment not equal to 1.
!
  IX = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  DO 10 I = 1,N
    SX(IX) = SA*SX(IX)
    IX = IX + INCX
   10 CONTINUE
  return
!
!     Code for increment equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 5.
!
   20 M = MOD(N,5)
  if (M  ==  0) GOTO 40
  DO 30 I = 1,M
    SX(I) = SA*SX(I)
   30 CONTINUE
  if (N  <  5) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,5
    SX(I) = SA*SX(I)
    SX(I+1) = SA*SX(I+1)
    SX(I+2) = SA*SX(I+2)
    SX(I+3) = SA*SX(I+3)
    SX(I+4) = SA*SX(I+4)
   50 CONTINUE
  return
end
