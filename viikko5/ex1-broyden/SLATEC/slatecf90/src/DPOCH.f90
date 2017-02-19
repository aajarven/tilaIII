  DOUBLE PRECISION FUNCTION DPOCH (A, X)
!
!! DPOCH evaluates a generalization of Pochhammer's symbol.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C1, C7A
!***TYPE      DOUBLE PRECISION (POCH-S, DPOCH-D)
!***KEYWORDS  FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate a double precision generalization of Pochhammer's symbol
! (A)-sub-X = GAMMA(A+X)/GAMMA(A) for double precision A and X.
! For X a non-negative integer, POCH(A,X) is just Pochhammer's symbol.
! This is a preliminary version that does not handle wrong arguments
! properly and may not properly handle the case when the result is
! computed to less than half of double precision.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D9LGMC, DFAC, DGAMMA, DGAMR, DLGAMS, DLNREL, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  DPOCH
  DOUBLE PRECISION A, X, ABSA, ABSAX, ALNGA, ALNGAX, AX, B, PI, &
    SGNGA, SGNGAX, DFAC, DLNREL, D9LGMC, DGAMMA, DGAMR, DCOT
  EXTERNAL DGAMMA
  SAVE PI
  DATA PI / 3.141592653589793238462643383279503D0 /
!***FIRST EXECUTABLE STATEMENT  DPOCH
  AX = A + X
  if (AX > 0.0D0) go to 30
  if (AINT(AX) /= AX) go to 30
!
  if (A  >  0.0D0 .OR. AINT(A)  /=  A) call XERMSG ('SLATEC', &
     'DPOCH', 'A+X IS NON-POSITIVE INTEGER BUT A IS NOT', 2, 2)
!
! WE KNOW HERE THAT BOTH A+X AND A ARE NON-POSITIVE INTEGERS.
!
  DPOCH = 1.0D0
  if (X == 0.D0) RETURN
!
  N = X
  if (MIN(A+X,A) < (-20.0D0)) go to 20
!
  IA = A
  DPOCH = (-1.0D0)**N * DFAC(-IA)/DFAC(-IA-N)
  return
!
 20   DPOCH = (-1.0D0)**N * EXP ((A-0.5D0)*DLNREL(X/(A-1.0D0)) &
    + X*LOG(-A+1.0D0-X) - X + D9LGMC(-A+1.0D0) - D9LGMC(-A-X+1.D0))
  return
!
! A+X IS NOT ZERO OR A NEGATIVE INTEGER.
!
 30   DPOCH = 0.0D0
  if (A <= 0.0D0 .AND. AINT(A) == A) RETURN
!
  N = ABS(X)
  if (DBLE(N) /= X .OR. N > 20) go to 50
!
! X IS A SMALL NON-POSITIVE INTEGER, PRESUMMABLY A COMMON CASE.
!
  DPOCH = 1.0D0
  if (N == 0) RETURN
  DO 40 I=1,N
    DPOCH = DPOCH * (A+I-1)
 40   CONTINUE
  return
!
 50   ABSAX = ABS(A+X)
  ABSA = ABS(A)
  if (MAX(ABSAX,ABSA) > 20.0D0) go to 60
  DPOCH = DGAMMA(A+X) * DGAMR(A)
  return
!
 60   if (ABS(X) > 0.5D0*ABSA) go to 70
!
! ABS(X) IS SMALL AND BOTH ABS(A+X) AND ABS(A) ARE LARGE.  THUS,
! A+X AND A MUST HAVE THE SAME SIGN.  FOR NEGATIVE A, WE USE
! GAMMA(A+X)/GAMMA(A) = GAMMA(-A+1)/GAMMA(-A-X+1) *
! SIN(PI*A)/SIN(PI*(A+X))
!
  B = A
  if (B < 0.0D0) B = -A - X + 1.0D0
  DPOCH = EXP ((B-0.5D0)*DLNREL(X/B) + X*LOG(B+X) - X &
    + D9LGMC(B+X) - D9LGMC(B) )
  if (A < 0.0D0 .AND. DPOCH /= 0.0D0) DPOCH = &
    DPOCH/(COS(PI*X) + DCOT(PI*A)*SIN(PI*X) )
  return
!
 70   call DLGAMS (A+X, ALNGAX, SGNGAX)
  call DLGAMS (A, ALNGA, SGNGA)
  DPOCH = SGNGAX * SGNGA * EXP(ALNGAX-ALNGA)
!
  return
end
