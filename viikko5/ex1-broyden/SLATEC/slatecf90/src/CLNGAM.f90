FUNCTION CLNGAM (ZIN)
!
!! CLNGAM computes the logarithm of the absolute value of the Gamma function.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      COMPLEX (ALNGAM-S, DLNGAM-D, CLNGAM-C)
!***KEYWORDS  ABSOLUTE VALUE, COMPLETE GAMMA FUNCTION, FNLIB, LOGARITHM,
!             SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! CLNGAM computes the natural log of the complex valued gamma function
! at ZIN, where ZIN is a complex number.  This is a preliminary version,
! which is not accurate.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  C9LGMC, CARG, CLNREL, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   780401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  CLNGAM
  COMPLEX CLNGAM
  COMPLEX ZIN, Z, CORR, CLNREL, C9LGMC
  LOGICAL FIRST
  SAVE PI, SQ2PIL, BOUND, DXREL, FIRST
  DATA PI / 3.14159265358979324E0 /
  DATA SQ2PIL / 0.91893853320467274E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  CLNGAM
  if (FIRST) THEN
     N = -0.30*LOG(R1MACH(3))
! BOUND = N*(0.1*EPS)**(-1/(2*N-1))/(PI*EXP(1))
     BOUND = 0.1171*N*(0.1*R1MACH(3))**(-1./(2*N-1))
     DXREL = SQRT (R1MACH(4))
  end if
  FIRST = .FALSE.
!
  Z = ZIN
  X = REAL(ZIN)
  Y = AIMAG(ZIN)
!
  CORR = (0.0, 0.0)
  CABSZ = ABS(Z)
  if (X >= 0.0 .AND. CABSZ > BOUND) go to 50
  if (X < 0.0 .AND. ABS(Y) > BOUND) go to 50
!
  if (CABSZ < BOUND) go to 20
!
! USE THE REFLECTION FORMULA FOR REAL(Z) NEGATIVE, ABS(Z) LARGE, AND
! ABS(AIMAG(Y)) SMALL.
!
  if (Y > 0.0) Z = CONJG (Z)
  CORR = EXP (-CMPLX(0.0,2.0*PI)*Z)
  if (REAL(CORR)  ==  1.0 .AND. AIMAG(CORR)  ==  0.0) call XERMSG &
     ('SLATEC', 'CLNGAM', 'Z IS A NEGATIVE INTEGER', 3, 2)
!
  CLNGAM = SQ2PIL + 1.0 - CMPLX(0.0,PI)*(Z-0.5) - CLNREL(-CORR) &
    + (Z-0.5)*LOG(1.0-Z) - Z - C9LGMC(1.0-Z)
  if (Y > 0.0) CLNGAM = CONJG (CLNGAM)
  return
!
! USE THE RECURSION RELATION FOR ABS(Z) SMALL.
!
 20   if (X >= (-0.5) .OR. ABS(Y) > DXREL) go to 30
  if (ABS((Z-AINT(X-0.5))/X)  <  DXREL) call XERMSG ('SLATEC', &
     'CLNGAM', &
     'ANSWER LT HALF PRECISION BECAUSE Z TOO NEAR NEGATIVE INTEGER', &
     1, 1)
!
 30   N = SQRT (BOUND**2 - Y**2) - X + 1.0
  ARGSUM = 0.0
  CORR = (1.0, 0.0)
  DO 40 I=1,N
    ARGSUM = ARGSUM + CARG(Z)
    CORR = Z*CORR
    Z = 1.0 + Z
 40   CONTINUE
!
  if (REAL(CORR)  ==  0.0 .AND. AIMAG(CORR)  ==  0.0) call XERMSG &
     ('SLATEC', 'CLNGAM', 'Z IS A NEGATIVE INTEGER', 3, 2)
  CORR = -CMPLX (LOG(ABS(CORR)), ARGSUM)
!
! USE STIRLING-S APPROXIMATION FOR LARGE Z.
!
 50   CLNGAM = SQ2PIL + (Z-0.5)*LOG(Z) - Z + CORR + C9LGMC(Z)
  return
!
end
