function COT (X)
!
!! COT computes the cotangent.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C4A
!***TYPE      SINGLE PRECISION (COT-S, DCOT-D, CCOT-C)
!***KEYWORDS  COTANGENT, ELEMENTARY FUNCTIONS, FNLIB, TRIGONOMETRIC
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! COT(X) calculates the cotangent of the real argument X.  X is in
! units of radians.
!
! Series for COT        on the interval  0.          to  6.25000D-02
!                                        with weighted error   3.76E-17
!                                         log weighted error  16.42
!                               significant figures required  15.51
!                                    decimal places required  16.88
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  COT
  DIMENSION COTCS(8)
  LOGICAL FIRST
  SAVE COTCS, PI2REC, NTERMS, XMAX, XSML, XMIN, SQEPS, FIRST
  DATA COTCS( 1) /    .24025916098295630E0 /
  DATA COTCS( 2) /   -.016533031601500228E0 /
  DATA COTCS( 3) /   -.000042998391931724E0 /
  DATA COTCS( 4) /   -.000000159283223327E0 /
  DATA COTCS( 5) /   -.000000000619109313E0 /
  DATA COTCS( 6) /   -.000000000002430197E0 /
  DATA COTCS( 7) /   -.000000000000009560E0 /
  DATA COTCS( 8) /   -.000000000000000037E0 /
  DATA PI2REC / .0116197723675813430E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  COT
  if (FIRST) THEN
     NTERMS = INITS (COTCS, 8, 0.1*R1MACH(3))
     XMAX = 1.0/R1MACH(4)
     XSML = SQRT (3.0*R1MACH(3))
     XMIN = EXP ( MAX(LOG(R1MACH(1)), -LOG(R1MACH(2))) + 0.01)
     SQEPS = SQRT (R1MACH(4))
  end if
  FIRST = .FALSE.
!
  Y = ABS(X)
  if (ABS(X)  <  XMIN) call XERMSG ('SLATEC', 'COT', &
     'ABS(X) IS ZERO OR SO SMALL COT OVERFLOWS', 2, 2)
  if (Y  >  XMAX) call XERMSG ('SLATEC', 'COT', &
     'NO PRECISION BECAUSE ABS(X) IS TOO BIG', 3, 2)
!
! CAREFULLY COMPUTE Y * (2/PI) = (AINT(Y) + REM(Y)) * (.625 + PI2REC)
! = AINT(.625*Y) + REM(.625*Y) + Y*PI2REC  =  AINT(.625*Y) + Z
! = AINT(.625*Y) + AINT(Z) + REM(Z)
!
  AINTY = AINT (Y)
  YREM = Y - AINTY
  PRODBG = 0.625*AINTY
  AINTY = AINT (PRODBG)
  Y = (PRODBG-AINTY) + 0.625*YREM + Y*PI2REC
  AINTY2 = AINT (Y)
  AINTY = AINTY + AINTY2
  Y = Y - AINTY2
!
  IFN = MOD (AINTY, 2.)
  if (IFN == 1) Y = 1.0 - Y
!
  if (ABS(X)  >  0.5 .AND. Y  <  ABS(X)*SQEPS) call XERMSG &
     ('SLATEC', 'COT', &
     'ANSWER LT HALF PRECISION, ABS(X) TOO BIG OR X NEAR N*PI ' // &
     '(N /= 0)' , 1, 1)
!
  if (Y > 0.25) go to 20
  COT = 1.0/X
  if (Y > XSML) COT = (0.5 + CSEVL (32.0*Y*Y-1., COTCS, NTERMS)) /Y
  go to 40
!
 20   if (Y > 0.5) go to 30
  COT = (0.5 + CSEVL (8.0*Y*Y-1., COTCS, NTERMS)) / (0.5*Y)
  COT = (COT**2 - 1.0) * 0.5 / COT
  go to 40
!
 30   COT = (0.5 + CSEVL (2.0*Y*Y-1., COTCS, NTERMS)) / (0.25*Y)
  COT = (COT**2 - 1.0) * 0.5 / COT
  COT = (COT**2 - 1.0) * 0.5 / COT
!
 40   if (X /= 0.) COT = SIGN (COT, X)
  if (IFN == 1) COT = -COT
!
  return
end
