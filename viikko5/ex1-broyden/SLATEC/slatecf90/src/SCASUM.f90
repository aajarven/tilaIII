function SCASUM (N, CX, INCX)
!
!! SCASUM computes the sum of the magnitudes of the real and ...
!  imaginary elements of a complex vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A3A
!***TYPE      COMPLEX (SASUM-S, DASUM-D, SCASUM-C)
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
!       CX  complex vector with N elements
!     INCX  storage spacing between elements of CX
!
!     --Output--
!   SCASUM  single precision result (zero if N  <=  0)
!
!     Returns sums of magnitudes of real and imaginary parts of
!     components of CX.  Note that this is not the L1 norm of CX.
!     CASUM = sum from 0 to N-1 of ABS(REAL(CX(IX+I*INCX))) +
!             ABS(IMAG(CX(IX+I*INCX))),
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
!***END PROLOGUE  SCASUM
  COMPLEX CX(*)
  INTEGER I, INCX, IX, N
!***FIRST EXECUTABLE STATEMENT  SCASUM
  SCASUM = 0.0E0
  if (N  <=  0) RETURN
!
  if (INCX  ==  1) GOTO 20
!
!     Code for increment not equal to 1.
!
  IX = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  DO 10 I = 1,N
    SCASUM = SCASUM + ABS(REAL(CX(IX))) + ABS(AIMAG(CX(IX)))
    IX = IX + INCX
   10 CONTINUE
  return
!
!     Code for increment equal to 1.
!
   20 DO 30 I = 1,N
    SCASUM = SCASUM + ABS(REAL(CX(I))) + ABS(AIMAG(CX(I)))
   30 CONTINUE
  return
end
