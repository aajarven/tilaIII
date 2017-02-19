subroutine DBSGQ8 (FUN, XT, BC, N, KK, ID, A, B, INBV, ERR, ANS, &
     IERR, WORK)
!
!! DBSGQ8 is subsidiary to DBFQAD.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (BSGQ8-S, DBSGQ8-D)
!***AUTHOR  Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract    **** A DOUBLE PRECISION routine ****
!
!        DBSGQ8, a modification of GAUS8, integrates the
!        product of FUN(X) by the ID-th derivative of a spline
!        DBVALU(XT,BC,N,KK,ID,X,INBV,WORK)  between limits A and B.
!
!     Description of Arguments
!
!        INPUT-- FUN,XT,BC,A,B,ERR are DOUBLE PRECISION
!        FUN - Name of external function of one argument which
!              multiplies DBVALU.
!        XT  - Knot array for DBVALU
!        BC  - B-coefficient array for DBVALU
!        N   - Number of B-coefficients for DBVALU
!        KK  - Order of the spline, KK >= 1
!        ID  - Order of the spline derivative, 0 <= ID <= KK-1
!        A   - Lower limit of integral
!        B   - Upper limit of integral (may be less than A)
!        INBV- Initialization parameter for DBVALU
!        ERR - Is a requested pseudorelative error tolerance.  Normally
!              pick a value of ABS(ERR) < 1D-3.  ANS will normally
!              have no more error than ABS(ERR) times the integral of
!              the absolute value of FUN(X)*DBVALU(XT,BC,N,KK,X,ID,
!              INBV,WORK).
!
!
!        OUTPUT-- ERR,ANS,WORK are DOUBLE PRECISION
!        ERR - Will be an estimate of the absolute error in ANS if the
!              input value of ERR was negative.  (ERR is unchanged if
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
!        WORK- Work vector of length 3*K for DBVALU
!
!***SEE ALSO  DBFQAD
!***ROUTINES CALLED  D1MACH, DBVALU, I1MACH, XERMSG
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
!***END PROLOGUE  DBSGQ8
!
  INTEGER ID, IERR, INBV, K, KK, KML, KMX, L, LMN, LMX, LR, MXL, &
   N, NBITS, NIB, NLMN, NLMX
  INTEGER I1MACH
  DOUBLE PRECISION A,AA,AE,ANIB,ANS,AREA,B,BC,C,CE,EE,EF,EPS,ERR, &
   EST,GL,GLR,GR,HH,SQ2,TOL,VL,VR,WORK,W1, W2, W3, W4, XT, X1, &
   X2, X3, X4, X, H
  DOUBLE PRECISION D1MACH, DBVALU, G8, FUN
  DIMENSION XT(*), BC(*), WORK(*)
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
  G8(X,H)=H*((W1*(FUN(X-X1*H)*DBVALU(XT,BC,N,KK,ID,X-X1*H,INBV,WORK) &
              + FUN(X+X1*H)*DBVALU(XT,BC,N,KK,ID,X+X1*H,INBV,WORK)) &
            +W2*(FUN(X-X2*H)*DBVALU(XT,BC,N,KK,ID,X-X2*H,INBV,WORK)+ &
                FUN(X+X2*H)*DBVALU(XT,BC,N,KK,ID,X+X2*H,INBV,WORK))) &
           +(W3*(FUN(X-X3*H)*DBVALU(XT,BC,N,KK,ID,X-X3*H,INBV,WORK)+ &
                 FUN(X+X3*H)*DBVALU(XT,BC,N,KK,ID,X+X3*H,INBV,WORK)) &
            +W4*(FUN(X-X4*H)*DBVALU(XT,BC,N,KK,ID,X-X4*H,INBV,WORK)+ &
               FUN(X+X4*H)*DBVALU(XT,BC,N,KK,ID,X+X4*H,INBV,WORK))))
!
!     INITIALIZE
!
!***FIRST EXECUTABLE STATEMENT  DBSGQ8
  K = I1MACH(14)
  ANIB = D1MACH(5)*K/0.30102000D0
  NBITS = INT(ANIB)
  NLMX = MIN((NBITS*5)/8,60)
  ANS = 0.0D0
  IERR = 1
  CE = 0.0D0
  if (A == B) go to 140
  LMX = NLMX
  LMN = NLMN
  if (B == 0.0D0) go to 10
  if (SIGN(1.0D0,B)*A <= 0.0D0) go to 10
  C = ABS(1.0D0-A/B)
  if (C > 0.1D0) go to 10
  if (C <= 0.0D0) go to 140
  ANIB = 0.5D0 - LOG(C)/0.69314718D0
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
   40 CE = CE + (EST-GLR)
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
  if ((MXL == 0) .OR. (ABS(CE) <= 2.0D0*TOL*AREA)) go to 140
  IERR = 2
  call XERMSG ('SLATEC', 'DBSGQ8', &
     'ANS IS PROBABLY INSUFFICIENTLY ACCURATE.', 3, 1)
  go to 140
  130 IERR = -1
  call XERMSG ('SLATEC', 'DBSGQ8', &
     'A AND B ARE TOO NEARLY EQUAL TO ALLOW NORMAL INTEGRATION. ' // &
     ' ANS IS SET TO ZERO AND IERR TO -1.', 1, -1)
  140 CONTINUE
  if (ERR < 0.0D0) ERR = CE
  return
end
