subroutine MPADD3 (X, Y, S, MED, RE)
!
!! MPADD3 is subsidiary to DQDOTA and DQDOTI.
!
!***LIBRARY   SLATEC
!***TYPE      ALL (MPADD3-A)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   Called by MPADD2; does inner loops of addition
!
!   The arguments X(*) and Y(*) and the variable R in COMMON are all
!   INTEGER arrays of size 30.  See the comments in the routine MPBLAS
!   for the reason for this choice.
!
!***SEE ALSO  DQDOTA, DQDOTI, MPBLAS
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    MPCOM
!***REVISION HISTORY  (YYMMDD)
!   791001  DATE WRITTEN
!   ??????  Modified for use with BLAS.  Blank COMMON changed to named
!           COMMON.  R given dimension 12.
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   930124  Increased Array size in MPCON for SUN -r8.  (RWC)
!***END PROLOGUE  MPADD3
  COMMON /MPCOM/ B, T, M, LUN, MXR, R(30)
  INTEGER B, T, R, X(*), Y(*), S, RE, C, TED
!***FIRST EXECUTABLE STATEMENT  MPADD3
  TED = T + MED
  I2 = T + 4
  I = I2
  C = 0
! CLEAR GUARD DIGITS TO RIGHT OF X DIGITS
   10 if (I <= TED) go to 20
  R(I) = 0
  I = I - 1
  go to 10
   20 if (S < 0) go to 130
! HERE DO ADDITION, EXPONENT(Y)  >=  EXPONENT(X)
  if (I < T) go to 40
   30 J = I - MED
  R(I) = X(J+2)
  I = I - 1
  if (I > T) go to 30
   40 if (I <= MED) go to 60
  J = I - MED
  C = Y(I+2) + X(J+2) + C
  if (C < B) go to 50
! CARRY GENERATED HERE
  R(I) = C - B
  C = 1
  I = I - 1
  go to 40
! NO CARRY GENERATED HERE
   50 R(I) = C
  C = 0
  I = I - 1
  go to 40
   60 if (I <= 0) go to 90
  C = Y(I+2) + C
  if (C < B) go to 70
  R(I) = 0
  C = 1
  I = I - 1
  go to 60
   70 R(I) = C
  I = I - 1
! NO CARRY POSSIBLE HERE
   80 if (I <= 0) RETURN
  R(I) = Y(I+2)
  I = I - 1
  go to 80
   90 if (C == 0) RETURN
! MUST SHIFT RIGHT HERE AS CARRY OFF END
  I2P = I2 + 1
  DO 100 J = 2, I2
  I = I2P - J
  100 R(I+1) = R(I)
  R(1) = 1
  RE = RE + 1
  return
! HERE DO SUBTRACTION, ABS(Y)  >  ABS(X)
  110 J = I - MED
  R(I) = C - X(J+2)
  C = 0
  if (R(I) >= 0) go to 120
! BORROW GENERATED HERE
  C = -1
  R(I) = R(I) + B
  120 I = I - 1
  130 if (I > T) go to 110
  140 if (I <= MED) go to 160
  J = I - MED
  C = Y(I+2) + C - X(J+2)
  if (C >= 0) go to 150
! BORROW GENERATED HERE
  R(I) = C + B
  C = -1
  I = I - 1
  go to 140
! NO BORROW GENERATED HERE
  150 R(I) = C
  C = 0
  I = I - 1
  go to 140
  160 if (I <= 0) RETURN
  C = Y(I+2) + C
  if (C >= 0) go to 70
  R(I) = C + B
  C = -1
  I = I - 1
  go to 160
end
