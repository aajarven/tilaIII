subroutine DPPGQ8 (FUN, LDC, C, XI, LXI, KK, ID, A, B, INPPV, ERR, &
     ANS, IERR)
!
!! DPPGQ8 is subsidiary to DPFQAD.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (PPGQ8-S, DPPGQ8-D)
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract    **** A DOUBLE PRECISION routine ****
!
!        DPPGQ8, a modification of GAUS8, integrates the
!        product of FUN(X) by the ID-th derivative of a spline
!        DPPVAL(LDC,C,XI,LXI,KK,ID,X,INPPV)  between limits A and B.
!
!     Description of Arguments
!
!      Input-- FUN,C,XI,A,B,ERR are DOUBLE PRECISION
!        FUN - Name of external function of one argument which
!              multiplies DPPVAL.
!        LDC - Leading dimension of matrix C, LDC  >=  KK
!        C   - Matrix of Taylor derivatives of dimension at least
!              (K,LXI)
!        XI  - Breakpoint vector of length LXI+1
!        LXI - Number of polynomial pieces
!        KK  - Order of the spline, KK  >=  1
!        ID  - Order of the spline derivative, 0  <=  ID  <=  KK-1
!        A   - Lower limit of integral
!        B   - Upper limit of integral (may be less than A)
!        INPPV- Initialization parameter for DPPVAL
!        ERR - Is a requested pseudorelative error tolerance.  Normally
!              pick a value of ABS(ERR)  <  1D-3.  ANS will normally
!              have no more error than ABS(ERR) times the integral of
!              the absolute value of FUN(X)*DPPVAL(LDC,C,XI,LXI,KK,ID,X,
!              INPPV).
!
!
!      Output-- ERR,ANS are DOUBLE PRECISION
!        ERR - Will be an estimate of the absolute error in ANS if the
!              input value of ERR was negative.  (ERR Is unchanged if
!              the input value of ERR was nonnegative.)  The estimated
!              error is solely for information to the user and should
!              not be used as a correction to the computed integral.
!        ANS - Computed value of integral
!        IERR- A status code
!            --Normal Codes
!               1 ANS most likely meets requested error tolerance,
!                 or A=B.
!              -1 A and B are too nearly equal to allow normal
!                 integration.  ANS is set to zero.
!            --Abnormal Code
!               2 ANS probably does not meet requested error tolerance.
!
!***SEE ALSO  DPFQAD
!***ROUTINES CALLED  D1MACH, DPPVAL, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated the AUTHOR section.  (WRB)
!***END PROLOGUE  DPPGQ8
!
  INTEGER ID,IERR,INPPV,K,KK,KML,KMX,L,LDC,LMN,LMX,LR,LXI,MXL, &
   NBITS, NIB, NLMN, NLMX
  INTEGER I1MACH
  DOUBLE PRECISION A,AA,AE,ANIB,ANS,AREA,B,BE,C,CC,EE,EF,EPS,ERR, &
   EST,GL,GLR,GR,HH,SQ2,TOL,VL,VR,W1, W2, W3, W4, XI, X1, &
   X2, X3, X4, X, H
  DOUBLE PRECISION D1MACH, DPPVAL, G8, FUN
  DIMENSION XI(*), C(LDC,*)
  DIMENSION AA(60), HH(60), LR(60), VL(60), GR(60)
  SAVE X1, X2, X3, X4, W1, W2, W3, W4, SQ2, NLMN, KMX, KML
  DATA X1, X2, X3, X4/ &
       1.83434642495649805D-01,     5.25532409916328986D-01, &
       7.96666477413626740D-01,     9.60289856497536232D-01/
  DATA W1, W2, W3, W4/ &
       3.62683783378361983D-01,     3.13706645877887287D-01, &
       2.22381034453374471D-01,     1.01228536290376259D-01/
  DATA SQ2/1.41421356D0/
  DATA NLMN/1/,KMX/5000/,KML/6/
  G8(X,H)= &
     H*((W1*(FUN(X-X1*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X-X1*H,INPPV) &
            +FUN(X+X1*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X+X1*H,INPPV)) &
        +W2*(FUN(X-X2*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X-X2*H,INPPV) &
            +FUN(X+X2*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X+X2*H,INPPV))) &
       +(W3*(FUN(X-X3*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X-X3*H,INPPV) &
            +FUN(X+X3*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X+X3*H,INPPV)) &
        +W4*(FUN(X-X4*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X-X4*H,INPPV) &
            +FUN(X+X4*H)*DPPVAL(LDC,C,XI,LXI,KK,ID,X+X4*H,INPPV))))
!
!     INITIALIZE
!
!***FIRST EXECUTABLE STATEMENT  DPPGQ8
  K = I1MACH(14)
  ANIB = D1MACH(5)*K/0.30102000D0
  NBITS = INT(ANIB)
  NLMX = MIN((NBITS*5)/8,60)
  ANS = 0.0D0
  IERR = 1
  BE = 0.0D0
  if (A == B) go to 140
  LMX = NLMX
  LMN = NLMN
  if (B == 0.0D0) go to 10
  if (SIGN(1.0D0,B)*A <= 0.0D0) go to 10
  CC = ABS(1.0D0-A/B)
  if (CC > 0.1D0) go to 10
  if (CC <= 0.0D0) go to 140
  ANIB = 0.5D0 - LOG(CC)/0.69314718D0
  NIB = INT(ANIB)
  LMX = MIN(NLMX,NBITS-NIB-7)
  if (LMX < 1) go to 130
  LMN = MIN(LMN,LMX)
   10 TOL = MAX(ABS(ERR),2.0D0**(5-NBITS))/2.0D0
  if (ERR == 0.0D0) TOL = SQRT(D1MACH(4))
  EPS = TOL
  HH(1) = (B-A)/4.0D0
  AA(1) = A
  LR(1) = 1
  L = 1
  EST = G8(AA(L)+2.0D0*HH(L),2.0D0*HH(L))
  K = 8
  AREA = ABS(EST)
  EF = 0.5D0
  MXL = 0
!
!     COMPUTE REFINED ESTIMATES, ESTIMATE THE ERROR, ETC.
!
   20 GL = G8(AA(L)+HH(L),HH(L))
  GR(L) = G8(AA(L)+3.0D0*HH(L),HH(L))
  K = K + 16
  AREA = AREA + (ABS(GL)+ABS(GR(L))-ABS(EST))
  GLR = GL + GR(L)
  EE = ABS(EST-GLR)*EF
  AE = MAX(EPS*AREA,TOL*ABS(GLR))
  if (EE-AE) 40, 40, 50
   30 MXL = 1
   40 BE = BE + (EST-GLR)
  if (LR(L)) 60, 60, 80
!
!     CONSIDER THE LEFT HALF OF THIS LEVEL
!
   50 if (K > KMX) LMX = KML
  if (L >= LMX) go to 30
  L = L + 1
  EPS = EPS*0.5D0
  EF = EF/SQ2
  HH(L) = HH(L-1)*0.5D0
  LR(L) = -1
  AA(L) = AA(L-1)
  EST = GL
  go to 20
!
!     PROCEED TO RIGHT HALF AT THIS LEVEL
!
   60 VL(L) = GLR
   70 EST = GR(L-1)
  LR(L) = 1
  AA(L) = AA(L) + 4.0D0*HH(L)
  go to 20
!
!     return ONE LEVEL
!
   80 VR = GLR
   90 if (L <= 1) go to 120
  L = L - 1
  EPS = EPS*2.0D0
  EF = EF*SQ2
  if (LR(L)) 100, 100, 110
  100 VL(L) = VL(L+1) + VR
  go to 70
  110 VR = VL(L+1) + VR
  go to 90
!
!      EXIT
!
  120 ANS = VR
  if ((MXL == 0) .OR. (ABS(BE) <= 2.0D0*TOL*AREA)) go to 140
  IERR = 2
  call XERMSG ('SLATEC', 'DPPGQ8', &
     'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.', 3, 1)
  go to 140
  130 IERR = -1
  call XERMSG ('SLATEC', 'DPPGQ8', &
     'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL ' // &
     'INTEGRATION.  ANSWER IS SET TO ZERO, AND IERR=-1.', 1, -1)
  140 CONTINUE
  if (ERR < 0.0D0) ERR = BE
  return
end
