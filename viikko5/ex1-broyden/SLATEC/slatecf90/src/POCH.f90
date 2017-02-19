function POCH (A, X)
!
!! POCH evaluates a generalization of Pochhammer's symbol.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C1, C7A
!***TYPE      SINGLE PRECISION (POCH-S, DPOCH-D)
!***KEYWORDS  FNLIB, POCHHAMMER, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate a generalization of Pochhammer's symbol
! (A)-sub-X = GAMMA(A+X)/GAMMA(A).  For X a non-negative integer,
! POCH(A,X) is just Pochhammer's symbol.  A and X are single precision.
! This is a preliminary version.  Error handling when POCH(A,X) is
! less than half precision is probably incorrect.  Grossly incorrect
! arguments are not handled properly.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  ALGAMS, ALNREL, FAC, GAMMA, GAMR, R9LGMC, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900727  Added EXTERNAL statement.  (WRB)
!***END PROLOGUE  POCH
  EXTERNAL GAMMA
  SAVE PI
  DATA PI / 3.141592653589793238E0 /
!***FIRST EXECUTABLE STATEMENT  POCH
  AX = A + X
  if (AX > 0.0) go to 30
  if (AINT(AX) /= AX) go to 30
!
  if (A  >  0.0 .OR. AINT(A)  /=  A) call XERMSG ('SLATEC', 'POCH', &
     'A+X IS NON-POSITIVE INTEGER BUT A IS NOT', 2, 2)
!
! WE KNOW HERE THAT BOTH A+X AND A ARE NON-POSITIVE INTEGERS.
!
  POCH = 1.0
  if (X == 0.0) RETURN
!
  N = X
  if (MIN(A+X,A) < (-20.0)) go to 20
!
  POCH = (-1.0)**N * FAC(-INT(A))/FAC(-INT(A)-N)
  return
!
 20   POCH = (-1.0)**N * EXP ((A-0.5)*ALNREL(X/(A-1.0)) &
    + X*LOG(-A+1.0-X) - X + R9LGMC(-A+1.) - R9LGMC(-A-X+1.) )
  return
!
! HERE WE KNOW A+X IS NOT ZERO OR A NEGATIVE INTEGER.
!
 30   POCH = 0.0
  if (A <= 0.0 .AND. AINT(A) == A) RETURN
!
  N = ABS(X)
  if (REAL(N) /= X .OR. N > 20) go to 50
!
! X IS A SMALL NON-POSITIVE INTEGER, PRESUMMABLY A COMMON CASE.
!
  POCH = 1.0
  if (N == 0) RETURN
  DO 40 I=1,N
    POCH = POCH * (A+I-1)
 40   CONTINUE
  return
!
 50   ABSAX = ABS(A+X)
  ABSA = ABS(A)
  if (MAX(ABSAX,ABSA) > 20.0) go to 60
  POCH = GAMMA(A+X)*GAMR(A)
  return
!
 60   if (ABS(X) > 0.5*ABSA) go to 70
!
! HERE ABS(X) IS SMALL AND BOTH ABS(A+X) AND ABS(A) ARE LARGE.  THUS,
! A+X AND A MUST HAVE THE SAME SIGN.  FOR NEGATIVE A, WE USE
! GAMMA(A+X)/GAMMA(A) = GAMMA(-A+1)/GAMMA(-A-X+1) *
! SIN(PI*A)/SIN(PI*(A+X))
!
  B = A
  if (B < 0.0) B = -A - X + 1.0
  POCH = EXP ((B-0.5)*ALNREL(X/B) + X*LOG(B+X) - X + &
    R9LGMC(B+X) - R9LGMC(B) )
  if (A < 0.0 .AND. POCH /= 0.0) POCH = POCH/(COS(PI*X) + &
    COT(PI*A)*SIN(PI*X))
  return
!
 70   call ALGAMS (A+X, ALNGAX, SGNGAX)
  call ALGAMS (A, ALNGA, SGNGA)
  POCH = SGNGAX * SGNGA * EXP(ALNGAX-ALNGA)
!
  return
end
