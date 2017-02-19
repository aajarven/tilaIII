subroutine DAXPY (N, DA, DX, INCX, DY, INCY)
!
!! DAXPY computes a constant times a vector plus a vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A7
!***TYPE      DOUBLE PRECISION (SAXPY-S, DAXPY-D, CAXPY-C)
!***KEYWORDS  BLAS, LINEAR ALGEBRA, TRIAD, VECTOR
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
!       DA  double precision scalar multiplier
!       DX  double precision vector with N elements
!     INCX  storage spacing between elements of DX
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DY  double precision result (unchanged if N  <=  0)
!
!     Overwrite double precision DY with double precision DA*DX + DY.
!     For I = 0 to N-1, replace  DY(LY+I*INCY) with DA*DX(LX+I*INCX) +
!       DY(LY+I*INCY),
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
!***END PROLOGUE  DAXPY
  DOUBLE PRECISION DX(*), DY(*), DA
!***FIRST EXECUTABLE STATEMENT  DAXPY
  if (N <= 0 .OR. DA == 0.0D0) RETURN
  if (INCX  ==  INCY) IF (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
    5 IX = 1
  IY = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  if (INCY  <  0) IY = (-N+1)*INCY + 1
  DO 10 I = 1,N
    DY(IY) = DY(IY) + DA*DX(IX)
    IX = IX + INCX
    IY = IY + INCY
   10 CONTINUE
  return
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 4.
!
   20 M = MOD(N,4)
  if (M  ==  0) go to 40
  DO 30 I = 1,M
    DY(I) = DY(I) + DA*DX(I)
   30 CONTINUE
  if (N  <  4) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,4
    DY(I) = DY(I) + DA*DX(I)
    DY(I+1) = DY(I+1) + DA*DX(I+1)
    DY(I+2) = DY(I+2) + DA*DX(I+2)
    DY(I+3) = DY(I+3) + DA*DX(I+3)
   50 CONTINUE
  return
!
!     Code for equal, positive, non-unit increments.
!
   60 NS = N*INCX
  DO 70 I = 1,NS,INCX
    DY(I) = DA*DX(I) + DY(I)
   70 CONTINUE
  return
end
