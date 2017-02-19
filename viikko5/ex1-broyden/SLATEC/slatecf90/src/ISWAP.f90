subroutine ISWAP (N, IX, INCX, IY, INCY)
!
!! ISWAP interchanges two vectors.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A5
!***TYPE      INTEGER (SSWAP-S, DSWAP-D, CSWAP-C, ISWAP-I)
!***KEYWORDS  BLAS, INTERCHANGE, LINEAR ALGEBRA, VECTOR
!***AUTHOR  Vandevender, W. H., (SNLA)
!***DESCRIPTION
!
!                Extended B L A S  Subprogram
!    Description of Parameters
!
!     --Input--
!        N  number of elements in input vector(s)
!       IX  integer vector with N elements
!     INCX  storage spacing between elements of IX
!       IY  integer vector with N elements
!     INCY  storage spacing between elements of IY
!
!     --Output--
!       IX  input vector IY (unchanged if N  <=  0)
!       IY  input vector IX (unchanged if N  <=  0)
!
!     Interchange integer IX and integer IY.
!     For I = 0 to N-1, interchange  IX(LX+I*INCX) and IY(LY+I*INCY),
!     where LX = 1 if INCX  >=  0, else LX = 1+(1-N)*INCX, and LY is
!     defined in a similar way using INCY.
!
!***REFERENCES  C. L. Lawson, R. J. Hanson, D. R. Kincaid and F. T.
!                 Krogh, Basic linear algebra subprograms for Fortran
!                 usage, Algorithm No. 539, Transactions on Mathematical
!                 Software 5, 3 (September 1979), pp. 308-323.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   850601  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  ISWAP
  INTEGER IX(*), IY(*), ITEMP1, ITEMP2, ITEMP3
!***FIRST EXECUTABLE STATEMENT  ISWAP
  if (N  <=  0) RETURN
  if (INCX  /=  INCY) go to 5
  if (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
    5 IIX = 1
  IIY = 1
  if (INCX  <  0) IIX = (1-N)*INCX + 1
  if (INCY  <  0) IIY = (1-N)*INCY + 1
  DO 10 I = 1,N
    ITEMP1 = IX(IIX)
    IX(IIX) = IY(IIY)
    IY(IIY) = ITEMP1
    IIX = IIX + INCX
    IIY = IIY + INCY
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
    ITEMP1 = IX(I)
    IX(I) = IY(I)
    IY(I) = ITEMP1
   30 CONTINUE
  if (N  <  3) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,3
    ITEMP1 = IX(I)
    ITEMP2 = IX(I+1)
    ITEMP3 = IX(I+2)
    IX(I) = IY(I)
    IX(I+1) = IY(I+1)
    IX(I+2) = IY(I+2)
    IY(I) = ITEMP1
    IY(I+1) = ITEMP2
    IY(I+2) = ITEMP3
   50 CONTINUE
  return
!
!     Code for equal, positive, non-unit increments.
!
   60 NS = N*INCX
  DO 70 I = 1,NS,INCX
    ITEMP1 = IX(I)
    IX(I) = IY(I)
    IY(I) = ITEMP1
   70 CONTINUE
  return
end
