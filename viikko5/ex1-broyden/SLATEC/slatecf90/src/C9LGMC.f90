FUNCTION C9LGMC (ZIN)
!
!! C9LGMC computes the log gamma correction factor so that ...
!  LOG(CGAMMA(Z)) = 0.5*LOG(2.*PI) + (Z-0.5)*LOG(Z) - Z + C9LGMC(Z).
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C7A
!***TYPE      COMPLEX (R9LGMC-S, D9LGMC-D, C9LGMC-C)
!***KEYWORDS  COMPLETE GAMMA FUNCTION, CORRECTION TERM, FNLIB,
!             LOG GAMMA, LOGARITHM, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute the LOG GAMMA correction term for large ABS(Z) when REAL(Z)
!  >=  0.0 and for large ABS(AIMAG(Y)) when REAL(Z)  <  0.0.  We find
! C9LGMC so that
!   LOG(Z) = 0.5*LOG(2.*PI) + (Z-0.5)*LOG(Z) - Z + C9LGMC(Z)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   780401  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!***END PROLOGUE  C9LGMC
  COMPLEX C9LGMC
  COMPLEX ZIN, Z, Z2INV
  DIMENSION BERN(11)
  LOGICAL FIRST
  SAVE BERN, NTERM, BOUND, XBIG, XMAX, FIRST
  DATA BERN( 1) /    .083333333333333333E0   /
  DATA BERN( 2) /   -.0027777777777777778E0  /
  DATA BERN( 3) /    .00079365079365079365E0 /
  DATA BERN( 4) /   -.00059523809523809524E0 /
  DATA BERN( 5) /    .00084175084175084175E0 /
  DATA BERN( 6) /   -.0019175269175269175E0  /
  DATA BERN( 7) /    .0064102564102564103E0  /
  DATA BERN( 8) /   -.029550653594771242E0   /
  DATA BERN( 9) /    .17964437236883057E0    /
  DATA BERN(10) /  -1.3924322169059011E0     /
  DATA BERN(11) /  13.402864044168392E0      /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  C9LGMC
  if (FIRST) THEN
     NTERM = -0.30*LOG(R1MACH(3))
     BOUND = 0.1170*NTERM*(0.1*R1MACH(3))**(-1./(2*NTERM-1))
     XBIG = 1.0/SQRT(R1MACH(3))
     XMAX = EXP (MIN(LOG(R1MACH(2)/12.0), -LOG(12.*R1MACH(1))) )
  end if
  FIRST = .FALSE.
!
  Z = ZIN
  X = REAL (Z)
  Y = AIMAG(Z)
  CABSZ = ABS(Z)

  if (X  <  0.0 .AND. ABS(Y)  <  BOUND) call XERMSG ('SLATEC', &
     'C9LGMC', 'NOT VALID FOR NEGATIVE REAL(Z) AND SMALL ' // &
     'ABS(AIMAG(Z))', 2, 2)
  if (CABSZ  <  BOUND) call XERMSG ('SLATEC', 'C9LGMC', &
     'NOT VALID FOR SMALL ABS(Z)', 3, 2)

  if (CABSZ >= XMAX) go to 50

  if (CABSZ >= XBIG) C9LGMC = 1.0/(12.0*Z)
  if (CABSZ >= XBIG) RETURN

  Z2INV = 1.0/Z**2
  C9LGMC = (0.0, 0.0)
  DO I=1,NTERM
    NDX = NTERM + 1 - I
    C9LGMC = BERN(NDX) + C9LGMC*Z2INV
  end do

  C9LGMC = C9LGMC/Z
  return

 50   C9LGMC = (0.0, 0.0)
  call XERMSG ('SLATEC', 'C9LGMC', 'Z SO BIG C9LGMC UNDERFLOWS', 1, &
     1)
  return
end
