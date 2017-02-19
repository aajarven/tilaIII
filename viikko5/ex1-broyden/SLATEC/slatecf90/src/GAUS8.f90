subroutine GAUS8 (FUN, A, B, ERR, ANS, IERR)
!
!! GAUS8 integrates a real function of one variable over a finite interval ...
!  using an adaptive 8-point Legendre-Gauss algorithm.  Intended primarily
!  for high accuracy integration or integration of smooth functions.
!
!***LIBRARY   SLATEC
!***CATEGORY  H2A1A1
!***TYPE      SINGLE PRECISION (GAUS8-S, DGAUS8-D)
!***KEYWORDS  ADAPTIVE QUADRATURE, AUTOMATIC INTEGRATOR,
!             GAUSS QUADRATURE, NUMERICAL INTEGRATION
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!        GAUS8 integrates real functions of one variable over finite
!        intervals using an adaptive 8-point Legendre-Gauss algorithm.
!        GAUS8 is intended primarily for high accuracy integration
!        or integration of smooth functions.
!
!     Description of Arguments
!
!        Input--
!        FUN - name of external function to be integrated.  This name
!              must be in an EXTERNAL statement in the calling program.
!              FUN must be a REAL function of one REAL argument.  The
!              value of the argument to FUN is the variable of
!              integration which ranges from A to B.
!        A   - lower limit of integration
!        B   - upper limit of integration (may be less than A)
!        ERR - is a requested pseudorelative error tolerance.  Normally
!              pick a value of ABS(ERR) so that STOL  <  ABS(ERR)  <=
!              1.0E-3 where STOL is the single precision unit roundoff
!              R1MACH(4).  ANS will normally have no more error than
!              ABS(ERR) times the integral of the absolute value of
!              FUN(X).  Usually, smaller values for ERR yield more
!              accuracy and require more function evaluations.
!
!              A negative value for ERR causes an estimate of the
!              absolute error in ANS to be returned in ERR.  Note that
!              ERR must be a variable (not a constant) in this case.
!              Note also that the user must reset the value of ERR
!              before making any more calls that use the variable ERR.
!
!        Output--
!        ERR - will be an estimate of the absolute error in ANS if the
!              input value of ERR was negative.  (ERR is unchanged if
!              the input value of ERR was non-negative.)  The estimated
!              error is solely for information to the user and should
!              not be used as a correction to the computed integral.
!        ANS - computed value of integral
!        IERR- a status code
!            --Normal codes
!               1 ANS most likely meets requested error tolerance,
!                 or A=B.
!              -1 A and B are too nearly equal to allow normal
!                 integration.  ANS is set to zero.
!            --Abnormal code
!               2 ANS probably does not meet requested error tolerance.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  I1MACH, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   810223  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!***END PROLOGUE  GAUS8
  INTEGER IERR, K, KML, KMX, L, LMN, LMX, LR, MXL, NBITS, &
   NIB, NLMN, NLMX
  INTEGER I1MACH
  REAL A, AA, AE, ANIB, ANS, AREA, B, C, CE, EE, EF, EPS, ERR, EST, &
   GL, GLR, GR, HH, SQ2, TOL, VL, VR, W1, W2, W3, W4, X1, X2, X3, &
   X4, X, H
  REAL R1MACH, G8, FUN
  DIMENSION AA(30), HH(30), LR(30), VL(30), GR(30)
  SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2, &
   NLMN, KMX, KML
  DATA X1, X2, X3, X4/ &
       1.83434642495649805E-01,     5.25532409916328986E-01, &
       7.96666477413626740E-01,     9.60289856497536232E-01/
  DATA W1, W2, W3, W4/ &
       3.62683783378361983E-01,     3.13706645877887287E-01, &
       2.22381034453374471E-01,     1.01228536290376259E-01/
  DATA SQ2/1.41421356E0/
  DATA NLMN/1/,KMX/5000/,KML/6/
  G8(X,H)=H*((W1*(FUN(X-X1*H) + FUN(X+X1*H)) &
             +W2*(FUN(X-X2*H) + FUN(X+X2*H))) &
            +(W3*(FUN(X-X3*H) + FUN(X+X3*H)) &
             +W4*(FUN(X-X4*H) + FUN(X+X4*H))))
!***FIRST EXECUTABLE STATEMENT  GAUS8
!
!     Initialize
!
  K = I1MACH(11)
  ANIB = R1MACH(5)*K/0.30102000E0
  NBITS = ANIB
  NLMX = MIN(30,(NBITS*5)/8)
  ANS = 0.0E0
  IERR = 1
  CE = 0.0E0
  if (A  ==  B) go to 140
  LMX = NLMX
  LMN = NLMN
  if (B  ==  0.0E0) go to 10
  if (SIGN(1.0E0,B)*A  <=  0.0E0) go to 10
  C = ABS(1.0E0-A/B)
  if (C  >  0.1E0) go to 10
  if (C  <=  0.0E0) go to 140
  ANIB = 0.5E0 - LOG(C)/0.69314718E0
  NIB = ANIB
  LMX = MIN(NLMX,NBITS-NIB-7)
  if (LMX  <  1) go to 130
  LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0E0**(5-NBITS))/2.0E0
  if (ERR  ==  0.0E0) TOL = SQRT(R1MACH(4))
  EPS = TOL
  HH(1) = (B-A)/4.0E0
  AA(1) = A
  LR(1) = 1
  L = 1
  EST = G8(AA(L)+2.0E0*HH(L),2.0E0*HH(L))
  K = 8
  AREA = ABS(EST)
  EF = 0.5E0
  MXL = 0
!
!     Compute refined estimates, estimate the error, etc.
!
   20 GL = G8(AA(L)+HH(L),HH(L))
  GR(L) = G8(AA(L)+3.0E0*HH(L),HH(L))
  K = K + 16
  AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
!     if (L  <  LMN) go to 11
  GLR = GL + GR(L)
  EE = ABS(EST-GLR)*EF
  AE = MAX(EPS*AREA,TOL*ABS(GLR))
  if (EE-AE) 40, 40, 50
   30 MXL = 1
   40 CE = CE + (EST-GLR)
  if (LR(L)) 60, 60, 80
!
!     Consider the left half of this level
!
   50 if (K  >  KMX) LMX = KML
  if (L  >=  LMX) go to 30
  L = L + 1
  EPS = EPS*0.5E0
  EF = EF/SQ2
  HH(L) = HH(L-1)*0.5E0
  LR(L) = -1
  AA(L) = AA(L-1)
  EST = GL
  go to 20
!
!     Proceed to right half at this level
!
   60 VL(L) = GLR
   70 EST = GR(L-1)
  LR(L) = 1
  AA(L) = AA(L) + 4.0E0*HH(L)
  go to 20
!
!     Return one level
!
   80 VR = GLR
   90 if (L  <=  1) go to 120
  L = L - 1
  EPS = EPS*2.0E0
  EF = EF*SQ2
  if (LR(L)) 100, 100, 110
  100 VL(L) = VL(L+1) + VR
  go to 70
  110 VR = VL(L+1) + VR
  go to 90
!
!     Exit
!
  120 ANS = VR
  if ((MXL == 0) .OR. (ABS(CE) <= 2.0E0*TOL*AREA)) go to 140
  IERR = 2
  call XERMSG ('SLATEC', 'GAUS8', &
     'ANS is probably insufficiently accurate.', 3, 1)
  go to 140
  130 IERR = -1
  call XERMSG ('SLATEC', 'GAUS8', &
     'A and B are too nearly equal to allow normal integration. $$' &
     // 'ANS is set to zero and IERR to -1.', 1, -1)
  140 if (ERR  <  0.0E0) ERR = CE
  return
end
