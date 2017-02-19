subroutine MPDIVI (X, IY, Z)
!
!! MPDIVI is subsidiary to DQDOTA and DQDOTI.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (MPDIVI-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!  Divides 'mp' X by the single-precision integer IY giving 'mp' Z.
!  This is much faster than division by an 'mp' number.
!
!  The arguments X(*) and Z(*), and the variable R in COMMON are all
!  INTEGER arrays of size 30.  See the comments in the routine MPBLAS
!  for the reason for this choice.
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***ROUTINES CALLED  MPCHK, MPERR, MPNZR, MPSTR, MPUNFL
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   ??????  Modified for use with BLAS.  Blank COMMON changed to named
!           COMMON.  R given dimension 12.
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPDIVI
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, X(*), Z(*), RS, RE, R1, C, C2, B2
!***FIRST EXECUTABLE STATEMENT  MPDIVI
  RS = X(1)
  J = IY
  if (J) 30, 10, 40
   10 WRITE (LUN, 20)
   20 FORMAT (' *** ATTEMPTED DIVISION BY ZERO IN call TO MPDIVI ***')
  go to 230
   30 J = -J
  RS = -RS
   40 RE = X(2)
! CHECK FOR ZERO DIVIDEND
  if (RS == 0) go to 120
! CHECK FOR DIVISION BY B
  if (J /= B) go to 50
  call MPSTR (X, Z)
  if (RE <= (-M)) go to 240
  Z(1) = RS
  Z(2) = RE - 1
  return
! CHECK FOR DIVISION BY 1 OR -1
   50 if (J /= 1) go to 60
  call MPSTR (X, Z)
  Z(1) = RS
  return
   60 C = 0
  I2 = T + 4
  I = 0
! if J*B NOT REPRESENTABLE AS AN INTEGER HAVE TO SIMULATE
! LONG DIVISION.   ASSUME AT LEAST 16-BIT WORD.
  B2 = MAX(8*B,32767/B)
  if (J >= B2) go to 130
! LOOK FOR FIRST NONZERO DIGIT IN QUOTIENT
   70 I = I + 1
  C = B*C
  if (I <= T) C = C + X(I+2)
  R1 = C/J
  if (R1) 210, 70, 80
! ADJUST EXPONENT AND GET T+4 DIGITS IN QUOTIENT
   80 RE = RE + 1 - I
  R(1) = R1
  C = B*(C - J*R1)
  KH = 2
  if (I >= T) go to 100
  KH = 1 + T - I
  DO 90 K = 2, KH
  I = I + 1
  C = C + X(I+2)
  R(K) = C/J
   90 C = B*(C - J*R(K))
  if (C < 0) go to 210
  KH = KH + 1
  100 DO 110 K = KH, I2
  R(K) = C/J
  110 C = B*(C - J*R(K))
  if (C < 0) go to 210
! NORMALIZE AND ROUND RESULT
  120 call MPNZR (RS, RE, Z, 0)
  return
! HERE NEED SIMULATED DOUBLE-PRECISION DIVISION
  130 C2 = 0
  J1 = J/B
  J2 = J - J1*B
  J11 = J1 + 1
! LOOK FOR FIRST NONZERO DIGIT
  140 I = I + 1
  C = B*C + C2
  C2 = 0
  if (I <= T) C2 = X(I+2)
  if (C-J1) 140, 150, 160
  150 if (C2 < J2) go to 140
! COMPUTE T+4 QUOTIENT DIGITS
  160 RE = RE + 1 - I
  K = 1
  go to 180
! MAIN LOOP FOR LARGE ABS(IY) CASE
  170 K = K + 1
  if (K > I2) go to 120
  I = I + 1
! GET APPROXIMATE QUOTIENT FIRST
  180 IR = C/J11
! NOW REDUCE SO OVERFLOW DOES NOT OCCUR
  IQ = C - IR*J1
  if (IQ < B2) go to 190
! HERE IQ*B WOULD POSSIBLY OVERFLOW SO INCREASE IR
  IR = IR + 1
  IQ = IQ - J1
  190 IQ = IQ*B - IR*J2
  if (IQ >= 0) go to 200
! HERE IQ NEGATIVE SO IR WAS TOO LARGE
  IR = IR - 1
  IQ = IQ + J
  200 if (I <= T) IQ = IQ + X(I+2)
  IQJ = IQ/J
! R(K) = QUOTIENT, C = REMAINDER
  R(K) = IQJ + IR
  C = IQ - J*IQJ
  if (C >= 0) go to 170
! CARRY NEGATIVE SO OVERFLOW MUST HAVE OCCURRED
  210 call MPCHK (1, 4)
  WRITE (LUN, 220)
  220 FORMAT (' *** INTEGER OVERFLOW IN MPDIVI, B TOO LARGE ***')
  230 call MPERR
  Z(1) = 0
  return
! UNDERFLOW HERE
  240 call MPUNFL(Z)
  return
end
