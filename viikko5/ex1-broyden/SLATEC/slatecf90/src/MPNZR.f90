subroutine MPNZR (RS, RE, Z, TRUNC)
!
!! MPNZR is subsidiary to DQDOTA and DQDOTI.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (MPNZR-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  Modified for use with BLAS.  Blank COMMON changed to named COMMON.
!  Assumes long (i.e. (t+4)-DIGIT) fraction in R, sign = RS, exponent
!  = RE.  Normalizes, and returns 'mp' result in Z. Integer arguments
!  RS and RE are not preserved. R*-rounding is used if TRUNC == 0
!
!  The argument Z(*) and the variable R in COMMON are INTEGER arrays
!  of size 30.  See the comments in the routine MPBLAS for the reason
!  for this choice.
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***ROUTINES CALLED  MPERR, MPOVFL, MPUNFL
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPNZR
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, Z(*), RE, RS, TRUNC, B2
!***FIRST EXECUTABLE STATEMENT  MPNZR
  I2 = T + 4
  if (RS /= 0) go to 20
! STORE ZERO IN Z
   10 Z(1) = 0
  return
! CHECK THAT SIGN = +-1
   20 if (ABS(RS) <= 1) go to 40
  WRITE (LUN, 30)
   30 FORMAT (' *** SIGN NOT 0, +1 OR -1 IN call TO MPNZR,', &
          ' POSSIBLE OVERWRITING PROBLEM ***')
  call MPERR
  go to 10
! LOOK FOR FIRST NONZERO DIGIT
   40 DO 50 I = 1, I2
  IS = I - 1
  if (R(I) > 0) go to 60
   50 CONTINUE
! FRACTION ZERO
  go to 10
   60 if (IS == 0) go to 90
! NORMALIZE
  RE = RE - IS
  I2M = I2 - IS
  DO 70 J = 1, I2M
  K = J + IS
   70 R(J) = R(K)
  I2P = I2M + 1
  DO 80 J = I2P, I2
   80 R(J) = 0
! CHECK TO SEE if TRUNCATION IS DESIRED
   90 if (TRUNC /= 0) go to 150
! SEE if ROUNDING NECESSARY
! TREAT EVEN AND ODD BASES DIFFERENTLY
  B2 = B/2
  if ((2*B2) /= B) go to 130
! B EVEN.  ROUND if R(T+1) >= B2 UNLESS R(T) ODD AND ALL ZEROS
! AFTER R(T+2).
  if (R(T+1) - B2) 150, 100, 110
  100 if (MOD(R(T),2) == 0) go to 110
  if ((R(T+2)+R(T+3)+R(T+4)) == 0) go to 150
! ROUND
  110 DO 120 J = 1, T
  I = T + 1 - J
  R(I) = R(I) + 1
  if (R(I) < B) go to 150
  120 R(I) = 0
! EXCEPTIONAL CASE, ROUNDED UP TO .10000...
  RE = RE + 1
  R(1) = 1
  go to 150
! ODD BASE, ROUND if R(T+1)...  >  1/2
  130 DO 140 I = 1, 4
  IT = T + I
  if (R(IT) - B2) 150, 140, 110
  140 CONTINUE
! CHECK FOR OVERFLOW
  150 if (RE <= M) go to 170
  WRITE (LUN, 160)
  160 FORMAT (' *** OVERFLOW OCCURRED IN MPNZR ***')
  call MPOVFL (Z)
  return
! CHECK FOR UNDERFLOW
  170 if (RE < (-M)) go to 190
! STORE RESULT IN Z
  Z(1) = RS
  Z(2) = RE
  DO 180 I = 1, T
  180 Z(I+2) = R(I)
  return
! UNDERFLOW HERE
  190 call MPUNFL (Z)
  return
end
