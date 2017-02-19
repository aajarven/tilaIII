subroutine XADD (X, IX, Y, IY, Z, IZ, IERROR)
!
!! XADD provides single-precision floating-point arithmetic ...
!  with an extended exponent range.
!
!***LIBRARY   SLATEC
!***CATEGORY  A3D
!***TYPE      SINGLE PRECISION (XADD-S, DXADD-D)
!***KEYWORDS  EXTENDED-RANGE SINGLE-PRECISION ARITHMETIC
!***AUTHOR  Lozier, Daniel W., (National Bureau of Standards)
!           Smith, John M., (NBS and George Mason University)
!***DESCRIPTION
!     REAL X, Y, Z
!     INTEGER IX, IY, IZ
!
!                  FORMS THE EXTENDED-RANGE SUM  (Z,IZ) =
!                  (X,IX) + (Y,IY).  (Z,IZ) IS ADJUSTED
!                  BEFORE RETURNING. THE INPUT OPERANDS
!                  NEED NOT BE IN ADJUSTED FORM, BUT THEIR
!                  PRINCIPAL PARTS MUST SATISFY
!                  RADIX**(-2L) <= ABS(X) <= RADIX**(2L),
!                  RADIX**(-2L) <= ABS(Y) <= RADIX**(2L).
!
!***SEE ALSO  XSET
!***REFERENCES  (NONE)
!***ROUTINES CALLED  XADJ
!***COMMON BLOCKS    XBLK2
!***REVISION HISTORY  (YYMMDD)
!   820712  DATE WRITTEN
!   881020  Revised to meet SLATEC CML recommendations.  (DWL and JMS)
!   901019  Revisions to prologue.  (DWL and WRB)
!   901106  Changed all specific intrinsics to generic.  (WRB)
!           Corrected order of sections in prologue and added TYPE
!           section.  (WRB)
!   920127  Revised PURPOSE section of prologue.  (DWL)
!***END PROLOGUE  XADD
  REAL X, Y, Z
  INTEGER IX, IY, IZ
  REAL RADIX, RADIXL, RAD2L, DLG10R
  INTEGER L, L2, KMAX
  COMMON /XBLK2/ RADIX, RADIXL, RAD2L, DLG10R, L, L2, KMAX
  SAVE /XBLK2/
!
!
!   THE CONDITIONS IMPOSED ON L AND KMAX BY THIS SUBROUTINE
! ARE
!     (1) 1  <  L  <=  0.5*LOGR(0.5*DZERO)
!
!     (2) NRADPL  <  L  <=  KMAX/6
!
!     (3) KMAX  <=  (2**NBITS - 4*L - 1)/2
!
! THESE CONDITIONS MUST BE MET BY APPROPRIATE CODING
! IN SUBROUTINE XSET.
!
!***FIRST EXECUTABLE STATEMENT  XADD
  IERROR=0
  if (X /= 0.0) go to 10
  Z = Y
  IZ = IY
  go to 220
   10 if (Y /= 0.0) go to 20
  Z = X
  IZ = IX
  go to 220
   20 CONTINUE
  if (IX >= 0 .AND. IY >= 0) go to 40
  if (IX < 0 .AND. IY < 0) go to 40
  if (ABS(IX) <= 6*L .AND. ABS(IY) <= 6*L) go to 40
  if (IX >= 0) go to 30
  Z = Y
  IZ = IY
  go to 220
   30 CONTINUE
  Z = X
  IZ = IX
  go to 220
   40 I = IX - IY
  if (I) 80, 50, 90
   50 if (ABS(X) > 1.0 .AND. ABS(Y) > 1.0) go to 60
  if (ABS(X) < 1.0 .AND. ABS(Y) < 1.0) go to 70
  Z = X + Y
  IZ = IX
  go to 220
   60 S = X/RADIXL
  T = Y/RADIXL
  Z = S + T
  IZ = IX + L
  go to 220
   70 S = X*RADIXL
  T = Y*RADIXL
  Z = S + T
  IZ = IX - L
  go to 220
   80 S = Y
  IS = IY
  T = X
  go to 100
   90 S = X
  IS = IX
  T = Y
  100 CONTINUE
!
!  AT THIS POINT, THE ONE OF (X,IX) OR (Y,IY) THAT HAS THE
! LARGER AUXILIARY INDEX IS STORED IN (S,IS). THE PRINCIPAL
! PART OF THE OTHER INPUT IS STORED IN T.
!
  I1 = ABS(I)/L
  I2 = MOD(ABS(I),L)
  if (ABS(T) >= RADIXL) go to 130
  if (ABS(T) >= 1.0) go to 120
  if (RADIXL*ABS(T) >= 1.0) go to 110
  J = I1 + 1
  T = T*RADIX**(L-I2)
  go to 140
  110 J = I1
  T = T*RADIX**(-I2)
  go to 140
  120 J = I1 - 1
  if (J < 0) go to 110
  T = T*RADIX**(-I2)/RADIXL
  go to 140
  130 J = I1 - 2
  if (J < 0) go to 120
  T = T*RADIX**(-I2)/RAD2L
  140 CONTINUE
!
!  AT THIS POINT, SOME OR ALL OF THE DIFFERENCE IN THE
! AUXILIARY INDICES HAS BEEN USED TO EFFECT A LEFT SHIFT
! OF T.  THE SHIFTED VALUE OF T SATISFIES
!
!       RADIX**(-2*L)  <=  ABS(T)  <=  1.0
!
! AND, if J=0, NO FURTHER SHIFTING REMAINS TO BE DONE.
!
  if (J == 0) go to 190
  if (ABS(S) >= RADIXL .OR. J > 3) go to 150
  if (ABS(S) >= 1.0) go to (180, 150, 150), J
  if (RADIXL*ABS(S) >= 1.0) go to (180, 170, 150), J
  go to (180, 170, 160), J
  150 Z = S
  IZ = IS
  go to 220
  160 S = S*RADIXL
  170 S = S*RADIXL
  180 S = S*RADIXL
  190 CONTINUE
!
!   AT THIS POINT, THE REMAINING DIFFERENCE IN THE
! AUXILIARY INDICES HAS BEEN USED TO EFFECT A RIGHT SHIFT
! OF S.  if THE SHIFTED VALUE OF S WOULD HAVE EXCEEDED
! RADIX**L, THEN (S,IS) IS RETURNED AS THE VALUE OF THE
! SUM.
!
  if (ABS(S) > 1.0 .AND. ABS(T) > 1.0) go to 200
  if (ABS(S) < 1.0 .AND. ABS(T) < 1.0) go to 210
  Z = S + T
  IZ = IS - J*L
  go to 220
  200 S = S/RADIXL
  T = T/RADIXL
  Z = S + T
  IZ = IS - J*L + L
  go to 220
  210 S = S*RADIXL
  T = T*RADIXL
  Z = S + T
  IZ = IS - J*L - L
  220 call XADJ(Z, IZ,IERROR)
  return
end
