subroutine CSSCAL (N, SA, CX, INCX)
!
!! CSSCAL scales a complex vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A6
!***TYPE      COMPLEX (CSSCAL-C)
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
!       CX  complex vector with N elements
!     INCX  storage spacing between elements of CX
!
!     --Output--
!       CX  scaled result (unchanged if N  <=  0)
!
!     Replace complex CX by (single precision SA) * (complex CX)
!     For I = 0 to N-1, replace CX(IX+I*INCX) with  SA * CX(IX+I*INCX),
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
!***END PROLOGUE  CSSCAL
  COMPLEX CX(*)
  REAL SA
  INTEGER I, INCX, IX, N
!***FIRST EXECUTABLE STATEMENT  CSSCAL
  if (N  <=  0) RETURN
!
  if (INCX  ==  1) GOTO 20
!
!     Code for increment not equal to 1.
!
  IX = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  DO 10 I = 1,N
    CX(IX) = SA*CX(IX)
    IX = IX + INCX
   10 CONTINUE
  return
!
!     Code for increment equal to 1.
!
   20 DO 30 I = 1,N
    CX(I) = SA*CX(I)
   30 CONTINUE
  return
end
