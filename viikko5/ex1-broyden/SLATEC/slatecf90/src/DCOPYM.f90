subroutine DCOPYM (N, DX, INCX, DY, INCY)
!
!! DCOPYM copies the negative of a vector to a vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A5
!***TYPE      DOUBLE PRECISION (SCOPYM-S, DCOPYM-D)
!***KEYWORDS  BLAS, COPY, VECTOR
!***AUTHOR  Kahaner, D. K., (NBS)
!***DESCRIPTION
!
!       Description of Parameters
!           The * Flags Output Variables
!
!       N   Number of elements in vector(s)
!      DX   Double precision vector with N elements
!    INCX   Storage spacing between elements of DX
!      DY*  Double precision negative copy of DX
!    INCY   Storage spacing between elements of DY
!
!      ***  Note that DY = -DX  ***
!
!     Copy negative of d.p. DX to d.p. DY.  For I=0 to N-1,
!     copy  -DX(LX+I*INCX) to DY(LY+I*INCY), where LX=1 if
!     INCX  >=  0, else LX = 1+(1-N)*INCX, and LY is defined
!     in a similar way using INCY.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!***END PROLOGUE  DCOPYM
  DOUBLE PRECISION DX(*), DY(*)
!***FIRST EXECUTABLE STATEMENT  DCOPYM
  if (N  <=  0) RETURN
  if (INCX  ==  INCY) IF (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
   5  IX=1
  IY=1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  if (INCY  <  0) IY = (-N+1)*INCY + 1
  DO 10 I = 1,N
    DY(IY) = -DX(IX)
    IX = IX + INCX
    IY = IY + INCY
   10 CONTINUE
  return
!
!     Code for both increments equal to 1.
!
!     Clean-up loop so remaining vector length is a multiple of 7.
!
   20 M = MOD(N,7)
  if (M  ==  0) go to 40
  DO 30 I = 1,M
    DY(I) = -DX(I)
   30 CONTINUE
  if (N  <  7) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,7
    DY(I) = -DX(I)
    DY(I+1) = -DX(I+1)
    DY(I+2) = -DX(I+2)
    DY(I+3) = -DX(I+3)
    DY(I+4) = -DX(I+4)
    DY(I+5) = -DX(I+5)
    DY(I+6) = -DX(I+6)
   50 CONTINUE
  return
!
!     Code for equal, positive, non-unit increments.
!
   60 NS = N*INCX
  DO 70 I = 1,NS,INCX
    DY(I) = -DX(I)
   70 CONTINUE
  return
end
