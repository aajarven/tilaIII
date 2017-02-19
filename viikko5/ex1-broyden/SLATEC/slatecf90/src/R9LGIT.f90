function R9LGIT (A, X, ALGAP1)
!
!! R9LGIT computes the logarithm of Tricomi's incomplete Gamma ...
!            function with Perron's continued fraction for large X and ...
!            A  >=  X.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7E
!***TYPE      SINGLE PRECISION (R9LGIT-S, D9LGIT-D)
!***KEYWORDS  FNLIB, INCOMPLETE GAMMA FUNCTION, LOGARITHM,
!             PERRON'S CONTINUED FRACTION, SPECIAL FUNCTIONS, TRICOMI
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the log of Tricomi's incomplete gamma function with Perron's
! continued fraction for large X and for A  >=  X.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770701  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  R9LGIT
  SAVE EPS, SQEPS
  DATA EPS, SQEPS / 2*0.0 /
!***FIRST EXECUTABLE STATEMENT  R9LGIT
  if (EPS == 0.0) EPS = 0.5*R1MACH(3)
  if (SQEPS == 0.0) SQEPS = SQRT(R1MACH(4))
!
  if (X  <=  0.0 .OR. A  <  X) call XERMSG ('SLATEC', 'R9LGIT', &
     'X SHOULD BE GT 0.0 AND LE A', 2, 2)
!
  AX = A + X
  A1X = AX + 1.0
  R = 0.0
  P = 1.0
  S = P
  DO 20 K=1,200
    FK = K
    T = (A+FK)*X*(1.0+R)
    R = T/((AX+FK)*(A1X+FK)-T)
    P = R*P
    S = S + P
    if (ABS(P) < EPS*S) go to 30
 20   CONTINUE
  call XERMSG ('SLATEC', 'R9LGIT', &
     'NO CONVERGENCE IN 200 TERMS OF CONTINUED FRACTION', 3, 2)
!
 30   HSTAR = 1.0 - X*S/A1X
  if (HSTAR  <  SQEPS) call XERMSG ('SLATEC', 'R9LGIT', &
     'RESULT LESS THAN HALF PRECISION', 1, 1)
!
  R9LGIT = -X - ALGAP1 - LOG(HSTAR)
!
  return
end
