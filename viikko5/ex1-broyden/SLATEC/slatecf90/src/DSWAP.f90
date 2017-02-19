subroutine DSWAP (N, DX, INCX, DY, INCY)
!
!! DSWAP interchanges two vectors.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A5
!***TYPE      DOUBLE PRECISION (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
!***KEYWORDS  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR
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
!       DY  double precision vector with N elements
!     INCY  storage spacing between elements of DY
!
!     --Output--
!       DX  input vector DY (unchanged if N  <=  0)
!       DY  input vector DX (unchanged if N  <=  0)
!
!     Interchange double precision DX and double precision DY.
!     For I = 0 to N-1, interchange  DX(LX+I*INCX) and DY(LY+I*INCY),
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
!***END PROLOGUE  DSWAP
  DOUBLE PRECISION DX(*), DY(*), DTEMP1, DTEMP2, DTEMP3
!***FIRST EXECUTABLE STATEMENT  DSWAP
  if (N  <=  0) RETURN
  if (INCX  ==  INCY) IF (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
    5 IX = 1
  IY = 1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  if (INCY  <  0) IY = (-N+1)*INCY + 1
  DO 10 I = 1,N
    DTEMP1 = DX(IX)
    DX(IX) = DY(IY)
    DY(IY) = DTEMP1
    IX = IX + INCX
    IY = IY + INCY
   10 CONTINUE
  return
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 3.
!
   20 M = MOD(N,3)
  if (M  ==  0) go to 40
  DO 30 I = 1,M
    DTEMP1 = DX(I)
    DX(I) = DY(I)
    DY(I) = DTEMP1
   30 CONTINUE
  if (N  <  3) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,3
    DTEMP1 = DX(I)
    DTEMP2 = DX(I+1)
    DTEMP3 = DX(I+2)
    DX(I) = DY(I)
    DX(I+1) = DY(I+1)
    DX(I+2) = DY(I+2)
    DY(I) = DTEMP1
    DY(I+1) = DTEMP2
    DY(I+2) = DTEMP3
   50 CONTINUE
  return
!
!     Code for equal, positive, non-unit increments.
!
   60 NS = N*INCX
  DO 70 I = 1,NS,INCX
    DTEMP1 = DX(I)
    DX(I) = DY(I)
    DY(I) = DTEMP1
   70 CONTINUE
  return
end
