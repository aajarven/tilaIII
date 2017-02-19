function R9CHU (A, B, Z)
!
!! R9CHU evaluates for large Z  Z**A * U(A,B,Z) where U is the logarithmic ...
!  confluent hypergeometric function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C11
!***TYPE      SINGLE PRECISION (R9CHU-S, D9CHU-D)
!***KEYWORDS  FNLIB, LOGARITHMIC CONFLUENT HYPERGEOMETRIC FUNCTION,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Evaluate for large Z  Z**A * U(A,B,Z)  where U is the logarithmic
! confluent hypergeometric function.  A rational approximation due to Y.
! L. Luke is used.  When U is not in the asymptotic region, i.e., when A
! or B is large compared with Z, considerable significance loss occurs.
! A warning is provided when the computed result is less than half
! precision.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770801  DATE WRITTEN
!   890206  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  R9CHU
  DIMENSION AA(4), BB(4)
  LOGICAL FIRST
  SAVE EPS, SQEPS, FIRST
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  R9CHU
  if (FIRST) THEN
     EPS = 4.0*R1MACH(4)
     SQEPS = SQRT (R1MACH(4))
  end if
  FIRST = .FALSE.
!
  BP = 1.0 + A - B
  AB = A*BP
  CT2 = 2.0*(Z-AB)
  SAB = A + BP
!
  BB(1) = 1.0
  AA(1) = 1.0
!
  CT3 = SAB + 1.0 + AB
  BB(2) = 1.0 + 2.0*Z/CT3
  AA(2) = 1.0 + CT2/CT3
!
  ANBN = CT3 + SAB + 3.0
  CT1 = 1.0 + 2.0*Z/ANBN
  BB(3) = 1.0 + 6.0*CT1*Z/CT3
  AA(3) = 1.0 + 6.0*AB/ANBN + 3.0*CT1*CT2/CT3
!
  DO 30 I=4,300
    X2I1 = 2*I - 3
    CT1 = X2I1/(X2I1-2.0)
    ANBN = ANBN + X2I1 + SAB
    CT2 = (X2I1 - 1.0) / ANBN
    C2 = X2I1*CT2 - 1.0
    D1Z = X2I1*2.0*Z/ANBN
!
    CT3 = SAB*CT2
    G1 = D1Z + CT1*(C2+CT3)
    G2 = D1Z - C2
    G3 = CT1*(1.0 - CT3 - 2.0*CT2)
!
    BB(4) = G1*BB(3) + G2*BB(2) + G3*BB(1)
    AA(4) = G1*AA(3) + G2*AA(2) + G3*AA(1)
    if (ABS(AA(4)*BB(1)-AA(1)*BB(4)) < EPS*ABS(BB(4)*BB(1))) &
      go to 40
!
! if OVERFLOWS OR UNDERFLOWS PROVE TO BE A PROBLEM, THE STATEMENTS
! BELOW COULD BE ALTERED TO INCORPORATE A DYNAMICALLY ADJUSTED SCALE
! FACTOR.
!
    DO 20 J=1,3
      BB(J) = BB(J+1)
      AA(J) = AA(J+1)
 20     CONTINUE
 30   CONTINUE
  call XERMSG ('SLATEC', 'R9CHU', 'NO CONVERGENCE IN 300 TERMS', 1, &
     2)
!
 40   R9CHU = AA(4)/BB(4)
!
  if (R9CHU  <  SQEPS .OR. R9CHU  >  1.0/SQEPS) call XERMSG &
     ('SLATEC', 'R9CHU', 'ANSWER LESS THAN HALF PRECISION', 2, 1)
!
  return
end
