subroutine R9UPAK (X, Y, N)
!
!! R9UPAK unpacks a floating point number X so that X = Y*2**N.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  A6B
!***TYPE      SINGLE PRECISION (R9UPAK-S, D9UPAK-D)
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
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  R9UPAK
!***FIRST EXECUTABLE STATEMENT  R9UPAK
  ABSX = ABS(X)
  N = 0
  if (X == 0.0E0) go to 30
!
   10 if (ABSX >= 0.5E0) go to 20
  N = N-1
  ABSX = ABSX*2.0E0
  go to 10
!
   20 if (ABSX < 1.0E0) go to 30
  N = N+1
  ABSX = ABSX*0.5E0
  go to 20
!
   30 Y = SIGN(ABSX,X)
  return
!
end
