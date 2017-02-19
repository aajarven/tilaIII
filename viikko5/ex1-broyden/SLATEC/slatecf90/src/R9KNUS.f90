subroutine R9KNUS (XNU, X, BKNU, BKNU1, ISWTCH)
!
!! R9KNUS: Bessel functions EXP(X)*K-SUB-XNU(X) and EXP(X)*K-SUB-XNU+1(X) ...
!  for 0.0  <=  XNU  <  1.0.
!
!***LIBRARY   SLATEC (FNLIB)
!***CATEGORY  C10B3
!***TYPE      SINGLE PRECISION (R9KNUS-S, D9KNUS-D)
!***KEYWORDS  BESSEL FUNCTION, FNLIB, SPECIAL FUNCTIONS
!***AUTHOR  Fullerton, W., (LANL)
!***DESCRIPTION
!
! Compute Bessel functions EXP(X) * K-sub-XNU (X)  and
! EXP(X) * K-sub-XNU+1 (X) for 0.0  <=  XNU  <  1.0 .
!
! Series for C0K        on the interval  0.          to  2.50000D-01
!                                        with weighted error   1.60E-17
!                                         log weighted error  16.79
!                               significant figures required  15.99
!                                    decimal places required  17.40
!
! Series for ZNU1       on the interval -7.00000D-01 to  0.
!                                        with weighted error   1.43E-17
!                                         log weighted error  16.85
!                               significant figures required  16.08
!                                    decimal places required  17.38
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  CSEVL, GAMMA, INITS, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   770601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900720  Routine changed from user-callable to subsidiary.  (WRB)
!   900727  Added EXTERNAL statement.  (WRB)
!   920618  Removed space from variable names.  (RWC, WRB)
!***END PROLOGUE  R9KNUS
  DIMENSION ALPHA(15), BETA(15), A(15), C0KCS(16), ZNU1CS(12)
  LOGICAL FIRST
  EXTERNAL GAMMA
  SAVE C0KCS, ZNU1CS, EULER, SQPI2, ALN2, NTC0K, NTZNU1, &
   XNUSML, XSML, ALNSML, ALNBIG, ALNEPS, FIRST
  DATA C0KCS( 1) /    .060183057242626108E0 /
  DATA C0KCS( 2) /   -.15364871433017286E0 /
  DATA C0KCS( 3) /   -.011751176008210492E0 /
  DATA C0KCS( 4) /   -.000852487888919795E0 /
  DATA C0KCS( 5) /   -.000061329838767496E0 /
  DATA C0KCS( 6) /   -.000004405228124551E0 /
  DATA C0KCS( 7) /   -.000000316312467283E0 /
  DATA C0KCS( 8) /   -.000000022710719382E0 /
  DATA C0KCS( 9) /   -.000000001630564460E0 /
  DATA C0KCS(10) /   -.000000000117069392E0 /
  DATA C0KCS(11) /   -.000000000008405206E0 /
  DATA C0KCS(12) /   -.000000000000603466E0 /
  DATA C0KCS(13) /   -.000000000000043326E0 /
  DATA C0KCS(14) /   -.000000000000003110E0 /
  DATA C0KCS(15) /   -.000000000000000223E0 /
  DATA C0KCS(16) /   -.000000000000000016E0 /
  DATA ZNU1CS( 1) /    .20330675699419173E0 /
  DATA ZNU1CS( 2) /    .14007793341321977E0 /
  DATA ZNU1CS( 3) /    .007916796961001613E0 /
  DATA ZNU1CS( 4) /    .000339801182532104E0 /
  DATA ZNU1CS( 5) /    .000011741975688989E0 /
  DATA ZNU1CS( 6) /    .000000339357570612E0 /
  DATA ZNU1CS( 7) /    .000000008425941769E0 /
  DATA ZNU1CS( 8) /    .000000000183336677E0 /
  DATA ZNU1CS( 9) /    .000000000003549698E0 /
  DATA ZNU1CS(10) /    .000000000000061903E0 /
  DATA ZNU1CS(11) /    .000000000000000981E0 /
  DATA ZNU1CS(12) /    .000000000000000014E0 /
  DATA EULER / 0.57721566490153286E0 /
  DATA SQPI2 / 1.2533141373155003E0 /
  DATA ALN2 / 0.69314718055994531E0 /
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  R9KNUS
  if (FIRST) THEN
     NTC0K = INITS (C0KCS, 16, 0.1*R1MACH(3))
     NTZNU1 = INITS (ZNU1CS, 12, 0.1*R1MACH(3))
!
     XNUSML = SQRT (R1MACH(3)/8.0)
     XSML = 0.1*R1MACH(3)
     ALNSML = LOG (R1MACH(1))
     ALNBIG = LOG (R1MACH(2))
     ALNEPS = LOG (0.1*R1MACH(3))
  end if
  FIRST = .FALSE.
!
  if (XNU  <  0. .OR. XNU  >=  1.0) call XERMSG ('SLATEC', &
     'R9KNUS', 'XNU MUST BE GE 0 AND LT 1', 1, 2)
  if (X  <=  0.) call XERMSG ('SLATEC', 'R9KNUS', 'X MUST BE GT 0', &
     2, 2)
!
  ISWTCH = 0
  if (X > 2.0) go to 50
!
! X IS SMALL.  COMPUTE K-SUB-XNU (X) AND THE DERIVATIVE OF K-SUB-XNU (X)
! THEN FIND K-SUB-XNU+1 (X).  XNU IS REDUCED TO THE INTERVAL (-.5,+.5)
! THEN TO (0., .5), BECAUSE K OF NEGATIVE ORDER (-NU) = K OF POSITIVE
! ORDER (+NU).
!
  V = XNU
  if (XNU > 0.5) V = 1.0 - XNU
!
! CAREFULLY FIND (X/2)**XNU AND Z**XNU WHERE Z = X*X/4.
  ALNZ = 2.0 * (LOG(X) - ALN2)
!
  if (X > XNU) go to 20
  if (-0.5*XNU*ALNZ-ALN2-LOG(XNU)  >  ALNBIG) call XERMSG &
     ('SLATEC', 'R9KNUS', 'X SO SMALL BESSEL K-SUB-XNU OVERFLOWS', &
     3, 2)
!
 20   VLNZ = V*ALNZ
  X2TOV = EXP (0.5*VLNZ)
  ZTOV = 0.0
  if (VLNZ > ALNSML) ZTOV = X2TOV**2
!
  A0 = 0.5*GAMMA(1.0+V)
  B0 = 0.5*GAMMA(1.0-V)
  C0 = -EULER
  if (ZTOV > 0.5 .AND. V > XNUSML) C0 = -0.75 + &
    CSEVL ((8.0*V)*V-1., C0KCS, NTC0K)
!
  if (ZTOV <= 0.5) ALPHA(1) = (A0-ZTOV*B0)/V
  if (ZTOV > 0.5) ALPHA(1) = C0 - ALNZ*(0.75 + &
    CSEVL (VLNZ/0.35+1.0, ZNU1CS, NTZNU1))*B0
  BETA(1) = -0.5*(A0+ZTOV*B0)
!
  Z = 0.0
  if (X > XSML) Z = 0.25*X*X
  NTERMS = MAX (2.0, 11.0+(8.*ALNZ-25.19-ALNEPS)/(4.28-ALNZ))
  DO 30 I=2,NTERMS
    XI = I - 1
    A0 = A0/(XI*(XI-V))
    B0 = B0/(XI*(XI+V))
    ALPHA(I) = (ALPHA(I-1)+2.0*XI*A0)/(XI*(XI+V))
    BETA(I) = (XI-0.5*V)*ALPHA(I) - ZTOV*B0
 30   CONTINUE
!
  BKNU = ALPHA(NTERMS)
  BKNUD = BETA(NTERMS)
  DO 40 II=2,NTERMS
    I = NTERMS + 1 - II
    BKNU = ALPHA(I) + BKNU*Z
    BKNUD = BETA(I) + BKNUD*Z
 40   CONTINUE
!
  EXPX = EXP(X)
  BKNU = EXPX*BKNU/X2TOV
!
  if (-0.5*(XNU+1.)*ALNZ-2.0*ALN2 > ALNBIG) ISWTCH = 1
  if (ISWTCH == 1) RETURN
  BKNUD = EXPX*BKNUD*2.0/(X2TOV*X)
!
  if (XNU <= 0.5) BKNU1 = V*BKNU/X - BKNUD
  if (XNU <= 0.5) RETURN
!
  BKNU0 = BKNU
  BKNU = -V*BKNU/X - BKNUD
  BKNU1 = 2.0*XNU*BKNU/X + BKNU0
  return
!
! X IS LARGE.  FIND K-SUB-XNU (X) AND K-SUB-XNU+1 (X) WITH Y. L. LUKE-S
! RATIONAL EXPANSION.
!
 50   SQRTX = SQRT(X)
  if (X > 1.0/XSML) go to 90
  AN = -1.56 + 4.0/X
  BN = -0.29 - 0.22/X
  NTERMS = MIN (15, MAX1 (3.0, AN+BN*ALNEPS))
!
  DO 80 INU=1,2
    XMU = 0.
    if (INU == 1 .AND. XNU > XNUSML) XMU = (4.0*XNU)*XNU
    if (INU == 2) XMU = 4.0*(ABS(XNU)+1.)**2
!
    A(1) = 1.0 - XMU
    A(2) = 9.0 - XMU
    A(3) = 25.0 - XMU
    if (A(2) == 0.) RESULT = SQPI2*(16.*X+XMU+7.)/(16.*X*SQRTX)
    if (A(2) == 0.) go to 70
!
    ALPHA(1) = 1.0
    ALPHA(2) = (16.*X+A(2))/A(2)
    ALPHA(3) = ((768.*X+48.*A(3))*X + A(2)*A(3))/(A(2)*A(3))
!
    BETA(1) = 1.0
    BETA(2) = (16.*X+(XMU+7.))/A(2)
    BETA(3) = ((768.*X+48.*(XMU+23.))*X + ((XMU+62.)*XMU+129.)) &
      / (A(2)*A(3))
!
    if (NTERMS < 4) go to 65
    DO 60 I=4,NTERMS
      N = I - 1
      X2N = 2*N - 1
!
      A(I) = (X2N+2.)**2 - XMU
      QQ = 16.*X2N/A(I)
      P1 = -X2N*(12*N*N-20*N-A(1))/((X2N-2.)*A(I)) - QQ*X
      P2 = (12*N*N-28*N+8-A(1))/A(I) - QQ*X
      P3 = -X2N*A(I-3)/((X2N-2.)*A(I))
!
      ALPHA(I) = -P1*ALPHA(I-1) - P2*ALPHA(I-2) - P3*ALPHA(I-3)
      BETA(I) = -P1*BETA(I-1) - P2*BETA(I-2) - P3*BETA(I-3)
 60     CONTINUE
!
 65     RESULT = SQPI2*BETA(NTERMS)/(SQRTX*ALPHA(NTERMS))
!
 70     if (INU == 1) BKNU = RESULT
    if (INU == 2) BKNU1 = RESULT
 80   CONTINUE
  return
!
 90   BKNU = SQPI2/SQRTX
  BKNU1 = BKNU
  return
!
end
