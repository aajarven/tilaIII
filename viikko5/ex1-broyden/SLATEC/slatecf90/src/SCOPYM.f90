subroutine SCOPYM (N, SX, INCX, SY, INCY)
!
!! SCOPYM copies the negative of a vector to a vector.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1A5
!***TYPE      SINGLE PRECISION (SCOPYM-S, DCOPYM-D)
!***KEYWORDS  BLAS, COPY, VECTOR
!***AUTHOR  Kahaner, D. K., (NBS)
!***DESCRIPTION
!
!       Description of Parameters
!           The * Flags Output Variables
!
!       N   Number of elements in vector(s)
!      SX   Real vector with N elements
!    INCX   Storage spacing between elements of SX
!      SY*  Real negative copy of SX
!    INCY   Storage spacing between elements of SY
!
!      ***  Note that SY = -SX  ***
!
!     Copy negative of real SX to real SY.  For I=0 to N-1,
!     copy  -SX(LX+I*INCX) to SY(LY+I*INCY), where LX=1 if
!     INCX  >=  0, else LX = 1+(1-N)*INCX, and LY is defined
!     in a similar way using INCY.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920310  Corrected definition of LX in DESCRIPTION.  (WRB)
!***END PROLOGUE  SCOPYM
  REAL SX(*),SY(*)
!***FIRST EXECUTABLE STATEMENT  SCOPYM
  if (N  <=  0) RETURN
  if (INCX  ==  INCY) IF (INCX-1) 5,20,60
!
!     Code for unequal or nonpositive increments.
!
    5 IX=1
  IY=1
  if (INCX  <  0) IX = (-N+1)*INCX + 1
  if (INCY  <  0) IY = (-N+1)*INCY + 1
  DO 10 I = 1,N
    SY(IY) = -SX(IX)
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
    SY(I) = -SX(I)
   30 CONTINUE
  if (N  <  7) RETURN
   40 MP1 = M + 1
  DO 50 I = MP1,N,7
    SY(I) = -SX(I)
    SY(I+1) = -SX(I+1)
    SY(I+2) = -SX(I+2)
    SY(I+3) = -SX(I+3)
    SY(I+4) = -SX(I+4)
    SY(I+5) = -SX(I+5)
    SY(I+6) = -SX(I+6)
   50 CONTINUE
  return
!
!     Code for equal, positive, non-unit increments.
!
   60 NS = N*INCX
  DO 70 I = 1,NS,INCX
    SY(I) = -SX(I)
   70 CONTINUE
  return
end
