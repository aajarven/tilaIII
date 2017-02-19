subroutine D9UPAK (X, Y, N)
!
!! D9UPAK unpacks a floating point number X so that X = Y*2**N.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  A6B
!***TYPE      DOUBLE PRECISION (R9UPAK-S, D9UPAK-D)
!***KEYWORDS  FNLIB, UNPACK
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
!   Unpack a floating point number X so that X = Y*2.0**N, where
!   0.5  <=  ABS(Y)  <  1.0.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   780701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900820  Corrected code to find Y between 0.5 and 1.0 rather than
!           between 0.05 and 1.0.  (WRB)
!***END PROLOGUE  D9UPAK
  DOUBLE PRECISION X,Y,ABSX
!***FIRST EXECUTABLE STATEMENT  D9UPAK
  ABSX = ABS(X)
  N = 0
  if (X == 0.0D0) go to 30
!
   10 if (ABSX >= 0.5D0) go to 20
  N = N-1
  ABSX = ABSX*2.0D0
  go to 10
!
   20 if (ABSX < 1.0D0) go to 30
  N = N+1
  ABSX = ABSX*0.5D0
  go to 20
!
   30 Y = SIGN(ABSX,X)
  return
!
end
