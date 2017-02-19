subroutine STEPS (F, NEQN, Y, X, H, EPS, WT, START, HOLD, K, KOLD, &
     CRASH, PHI, P, YP, PSI, ALPHA, BETA, SIG, V, W, G, PHASE1, NS, &
     NORND, KSTEPS, TWOU, FOURU, XOLD, KPREV, IVC, IV, KGI, GI, &
     RPAR, IPAR)
!
!! STEPS integrates a system of ordinary differential equations one step.
!
!***LIBRARY   SLATEC (DEPAC)
!***CATEGORY  I1A1B
!***TYPE      SINGLE PRECISION (STEPS-S, DSTEPS-D)
!***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
!             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR
!***AUTHOR  Shampine, L. F., (SNLA)
!           Gordon, M. K., (SNLA)
!             MODIFIED BY H.A. WATTS
!***DESCRIPTION
!
!   Written by L. F. Shampine and M. K. Gordon
!
!   Abstract
!
!   Subroutine  STEPS  is normally used indirectly through subroutine
!   DEABM .  Because  DEABM  suffices for most problems and is much
!   easier to use, using it should be considered before using  STEPS
!   alone.
!
!   Subroutine STEPS integrates a system of  NEQN  first order ordinary
!   differential equations one step, normally from X to X+H, using a
!   modified divided difference form of the Adams Pece formulas.  Local
!   extrapolation is used to improve absolute stability and accuracy.
!   The code adjusts its order and step size to control the local error
!   per unit step in a generalized sense.  Special devices are included
!   to control roundoff error and to detect when the user is requesting
!   too much accuracy.
!
!   This code is completely explained and documented in the text,
!   Computer Solution of Ordinary Differential Equations, The Initial
!   Value Problem  by L. F. Shampine and M. K. Gordon.
!   Further details on use of this code are available in "Solving
!   Ordinary Differential Equations with ODE, STEP, and INTRP",
!   by L. F. Shampine and M. K. Gordon, SLA-73-1060.
!
!
!   The parameters represent --
!      F -- subroutine to evaluate derivatives
!      NEQN -- number of equations to be integrated
!      Y(*) -- solution vector at X
!      X -- independent variable
!      H -- appropriate step size for next step.  Normally determined by
!           code
!      EPS -- local error tolerance
!      WT(*) -- vector of weights for error criterion
!      START -- logical variable set .TRUE. for first step,  .FALSE.
!           otherwise
!      HOLD -- step size used for last successful step
!      K -- appropriate order for next step (determined by code)
!      KOLD -- order used for last successful step
!      CRASH -- logical variable set .TRUE. when no step can be taken,
!           .FALSE. otherwise.
!      YP(*) -- derivative of solution vector at  X  after successful
!           step
!      KSTEPS -- counter on attempted steps
!      TWOU -- 2.*U where U is machine unit roundoff quantity
!      FOURU -- 4.*U where U is machine unit roundoff quantity
!      RPAR,IPAR -- parameter arrays which you may choose to use
!            for communication between your program and subroutine F.
!            They are not altered or used by STEPS.
!   The variables X,XOLD,KOLD,KGI and IVC and the arrays Y,PHI,ALPHA,G,
!   W,P,IV and GI are required for the interpolation subroutine SINTRP.
!   The remaining variables and arrays are included in the call list
!   only to eliminate local retention of variables between calls.
!
!   Input to STEPS
!
!      First call --
!
!   The user must provide storage in his calling program for all arrays
!   in the call list, namely
!
!     DIMENSION Y(NEQN),WT(NEQN),PHI(NEQN,16),P(NEQN),YP(NEQN),PSI(12),
!    1  ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10),
!    2  RPAR(*),IPAR(*)
!
!    **Note**
!
!   The user must also declare  START ,  CRASH ,  PHASE1  and  NORND
!   logical variables and  F  an EXTERNAL subroutine, supply the
!   subroutine  F(X,Y,YP)  to evaluate
!      DY(I)/DX = YP(I) = F(X,Y(1),Y(2),...,Y(NEQN))
!   and initialize only the following parameters.
!      NEQN -- number of equations to be integrated
!      Y(*) -- vector of initial values of dependent variables
!      X -- initial value of the independent variable
!      H -- nominal step size indicating direction of integration
!           and maximum size of step.  Must be variable
!      EPS -- local error tolerance per step.  Must be variable
!      WT(*) -- vector of non-zero weights for error criterion
!      START -- .TRUE.
!      YP(*) -- vector of initial derivative values
!      KSTEPS -- set KSTEPS to zero
!      TWOU -- 2.*U where U is machine unit roundoff quantity
!      FOURU -- 4.*U where U is machine unit roundoff quantity
!   Define U to be the machine unit roundoff quantity by calling
!   the function routine  R1MACH,  U = R1MACH(4), or by
!   computing U so that U is the smallest positive number such
!   that 1.0+U  >  1.0.
!
!   STEPS  requires that the L2 norm of the vector with components
!   LOCAL ERROR(L)/WT(L)  be less than  EPS  for a successful step.  The
!   array  WT  allows the user to specify an error test appropriate
!   for his problem.  For example,
!      WT(L) = 1.0  specifies absolute error,
!            = ABS(Y(L))  error relative to the most recent value of the
!                 L-th component of the solution,
!            = ABS(YP(L))  error relative to the most recent value of
!                 the L-th component of the derivative,
!            = MAX(WT(L),ABS(Y(L)))  error relative to the largest
!                 magnitude of L-th component obtained so far,
!            = ABS(Y(L))*RELERR/EPS + ABSERR/EPS  specifies a mixed
!                 relative-absolute test where  RELERR  is relative
!                 error,  ABSERR  is absolute error and  EPS =
!                 MAX(RELERR,ABSERR) .
!
!      Subsequent calls --
!
!   Subroutine  STEPS  is designed so that all information needed to
!   continue the integration, including the step size  H  and the order
!   K , is returned with each step.  With the exception of the step
!   size, the error tolerance, and the weights, none of the parameters
!   should be altered.  The array  WT  must be updated after each step
!   to maintain relative error tests like those above.  Normally the
!   integration is continued just beyond the desired endpoint and the
!   solution interpolated there with subroutine  SINTRP .  If it is
!   impossible to integrate beyond the endpoint, the step size may be
!   reduced to hit the endpoint since the code will not take a step
!   larger than the  H  input.  Changing the direction of integration,
!   i.e., the sign of  H , requires the user set  START = .TRUE. before
!   calling  STEPS  again.  This is the only situation in which  START
!   should be altered.
!
!   Output from STEPS
!
!      Successful Step --
!
!   The subroutine returns after each successful step with  START  and
!   CRASH  set .FALSE. .  X  represents the independent variable
!   advanced one step of length  HOLD  from its value on input and  Y
!   the solution vector at the new value of  X .  All other parameters
!   represent information corresponding to the new  X  needed to
!   continue the integration.
!
!      Unsuccessful Step --
!
!   When the error tolerance is too small for the machine precision,
!   the subroutine returns without taking a step and  CRASH = .TRUE. .
!   An appropriate step size and error tolerance for continuing are
!   estimated and all other information is restored as upon input
!   before returning.  To continue with the larger tolerance, the user
!   just calls the code again.  A restart is neither required nor
!   desirable.
!
!***REFERENCES  L. F. Shampine and M. K. Gordon, Solving ordinary
!                 differential equations with ODE, STEP, and INTRP,
!                 Report SLA-73-1060, Sandia Laboratories, 1973.
!***ROUTINES CALLED  HSTART, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   740101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  STEPS
!
  LOGICAL START,CRASH,PHASE1,NORND
  DIMENSION Y(*),WT(*),PHI(NEQN,16),P(*),YP(*),PSI(12), &
    ALPHA(12),BETA(12),SIG(13),V(12),W(12),G(13),GI(11),IV(10), &
    RPAR(*),IPAR(*)
  DIMENSION TWO(13),GSTR(13)
  EXTERNAL F
  SAVE TWO, GSTR
!
  DATA  TWO(1),TWO(2),TWO(3),TWO(4),TWO(5),TWO(6),TWO(7),TWO(8), &
        TWO(9),TWO(10),TWO(11),TWO(12),TWO(13) /2.0,4.0,8.0,16.0, &
        32.0,64.0,128.0,256.0,512.0,1024.0,2048.0,4096.0,8192.0/
  DATA  GSTR(1),GSTR(2),GSTR(3),GSTR(4),GSTR(5),GSTR(6),GSTR(7), &
        GSTR(8),GSTR(9),GSTR(10),GSTR(11),GSTR(12),GSTR(13)/0.500, &
        0.0833,0.0417,0.0264,0.0188,0.0143,0.0114,0.00936,0.00789, &
        0.00679,0.00592,0.00524,0.00468/
!
!
!       ***     BEGIN BLOCK 0     ***
!   CHECK if STEP SIZE OR ERROR TOLERANCE IS TOO SMALL FOR MACHINE
!   PRECISION.  if FIRST STEP, INITIALIZE PHI ARRAY AND ESTIMATE A
!   STARTING STEP SIZE.
!                   ***
!
!   if STEP SIZE IS TOO SMALL, DETERMINE AN ACCEPTABLE ONE
!
!***FIRST EXECUTABLE STATEMENT  STEPS
  CRASH = .TRUE.
  if ( ABS(H)  >=  FOURU*ABS(X)) go to 5
  H = SIGN(FOURU*ABS(X),H)
  return
 5    P5EPS = 0.5*EPS
!
!   if ERROR TOLERANCE IS TOO SMALL, INCREASE IT TO AN ACCEPTABLE VALUE
!
  ROUND = 0.0
  DO 10 L = 1,NEQN
 10     ROUND = ROUND + (Y(L)/WT(L))**2
  ROUND = TWOU*SQRT(ROUND)
  if ( P5EPS  >=  ROUND) go to 15
  EPS = 2.0*ROUND*(1.0 + FOURU)
  return
 15   CRASH = .FALSE.
  G(1) = 1.0
  G(2) = 0.5
  SIG(1) = 1.0
  if ( .NOT.START) go to 99
!
!   INITIALIZE.  COMPUTE APPROPRIATE STEP SIZE FOR FIRST STEP
!
!     call F(X,Y,YP,RPAR,IPAR)
!     SUM = 0.0
  DO 20 L = 1,NEQN
    PHI(L,1) = YP(L)
   20   PHI(L,2) = 0.0
!20     SUM = SUM + (YP(L)/WT(L))**2
!     SUM = SQRT(SUM)
!     ABSH = ABS(H)
!     if ( EPS  <  16.0*SUM*H*H) ABSH = 0.25*SQRT(EPS/SUM)
!     H = SIGN(MAX(ABSH,FOURU*ABS(X)),H)
!
  U = R1MACH(4)
  BIG = SQRT(R1MACH(2))
  call HSTART (F,NEQN,X,X+H,Y,YP,WT,1,U,BIG, &
               PHI(1,3),PHI(1,4),PHI(1,5),PHI(1,6),RPAR,IPAR,H)
!
  HOLD = 0.0
  K = 1
  KOLD = 0
  KPREV = 0
  START = .FALSE.
  PHASE1 = .TRUE.
  NORND = .TRUE.
  if ( P5EPS  >  100.0*ROUND) go to 99
  NORND = .FALSE.
  DO 25 L = 1,NEQN
 25     PHI(L,15) = 0.0
 99   IFAIL = 0
!       ***     END BLOCK 0     ***
!
!       ***     BEGIN BLOCK 1     ***
!   COMPUTE COEFFICIENTS OF FORMULAS FOR THIS STEP.  AVOID COMPUTING
!   THOSE QUANTITIES NOT CHANGED WHEN STEP SIZE IS NOT CHANGED.
!                   ***
!
 100  KP1 = K+1
  KP2 = K+2
  KM1 = K-1
  KM2 = K-2
!
!   NS IS THE NUMBER OF STEPS TAKEN WITH SIZE H, INCLUDING THE CURRENT
!   ONE.  WHEN K < NS, NO COEFFICIENTS CHANGE
!
  if ( H  /=  HOLD) NS = 0
  if (NS <= KOLD) NS = NS+1
  NSP1 = NS+1
  if (K  <  NS) go to 199
!
!   COMPUTE THOSE COMPONENTS OF ALPHA(*),BETA(*),PSI(*),SIG(*) WHICH
!   ARE CHANGED
!
  BETA(NS) = 1.0
  REALNS = NS
  ALPHA(NS) = 1.0/REALNS
  TEMP1 = H*REALNS
  SIG(NSP1) = 1.0
  if ( K  <  NSP1) go to 110
  DO 105 I = NSP1,K
    IM1 = I-1
    TEMP2 = PSI(IM1)
    PSI(IM1) = TEMP1
    BETA(I) = BETA(IM1)*PSI(IM1)/TEMP2
    TEMP1 = TEMP2 + H
    ALPHA(I) = H/TEMP1
    REALI = I
 105    SIG(I+1) = REALI*ALPHA(I)*SIG(I)
 110  PSI(K) = TEMP1
!
!   COMPUTE COEFFICIENTS G(*)
!
!   INITIALIZE V(*) AND SET W(*).
!
  if ( NS  >  1) go to 120
  DO 115 IQ = 1,K
    TEMP3 = IQ*(IQ+1)
    V(IQ) = 1.0/TEMP3
 115    W(IQ) = V(IQ)
  IVC = 0
  KGI = 0
  if (K  ==  1) go to 140
  KGI = 1
  GI(1) = W(2)
  go to 140
!
!   if ORDER WAS RAISED, UPDATE DIAGONAL PART OF V(*)
!
 120  if ( K  <=  KPREV) go to 130
  if (IVC  ==  0) go to 122
  JV = KP1 - IV(IVC)
  IVC = IVC - 1
  go to 123
 122  JV = 1
  TEMP4 = K*KP1
  V(K) = 1.0/TEMP4
  W(K) = V(K)
  if (K  /=  2) go to 123
  KGI = 1
  GI(1) = W(2)
 123  NSM2 = NS-2
  if ( NSM2  <  JV) go to 130
  DO 125 J = JV,NSM2
    I = K-J
    V(I) = V(I) - ALPHA(J+1)*V(I+1)
 125    W(I) = V(I)
  if (I  /=  2) go to 130
  KGI = NS - 1
  GI(KGI) = W(2)
!
!   UPDATE V(*) AND SET W(*)
!
 130  LIMIT1 = KP1 - NS
  TEMP5 = ALPHA(NS)
  DO 135 IQ = 1,LIMIT1
    V(IQ) = V(IQ) - TEMP5*V(IQ+1)
 135    W(IQ) = V(IQ)
  G(NSP1) = W(1)
  if (LIMIT1  ==  1) go to 137
  KGI = NS
  GI(KGI) = W(2)
 137  W(LIMIT1+1) = V(LIMIT1+1)
  if (K  >=  KOLD) go to 140
  IVC = IVC + 1
  IV(IVC) = LIMIT1 + 2
!
!   COMPUTE THE G(*) IN THE WORK VECTOR W(*)
!
 140  NSP2 = NS + 2
  KPREV = K
  if ( KP1  <  NSP2) go to 199
  DO 150 I = NSP2,KP1
    LIMIT2 = KP2 - I
    TEMP6 = ALPHA(I-1)
    DO 145 IQ = 1,LIMIT2
 145      W(IQ) = W(IQ) - TEMP6*W(IQ+1)
 150    G(I) = W(1)
 199    CONTINUE
!       ***     END BLOCK 1     ***
!
!       ***     BEGIN BLOCK 2     ***
!   PREDICT A SOLUTION P(*), EVALUATE DERIVATIVES USING PREDICTED
!   SOLUTION, ESTIMATE LOCAL ERROR AT ORDER K AND ERRORS AT ORDERS K,
!   K-1, K-2 AS if CONSTANT STEP SIZE WERE USED.
!                   ***
!
!   INCREMENT COUNTER ON ATTEMPTED STEPS
!
  KSTEPS = KSTEPS + 1
!
!   CHANGE PHI TO PHI STAR
!
  if ( K  <  NSP1) go to 215
  DO 210 I = NSP1,K
    TEMP1 = BETA(I)
    DO 205 L = 1,NEQN
 205      PHI(L,I) = TEMP1*PHI(L,I)
 210    CONTINUE
!
!   PREDICT SOLUTION AND DIFFERENCES
!
 215  DO 220 L = 1,NEQN
    PHI(L,KP2) = PHI(L,KP1)
    PHI(L,KP1) = 0.0
 220    P(L) = 0.0
  DO 230 J = 1,K
    I = KP1 - J
    IP1 = I+1
    TEMP2 = G(I)
    DO 225 L = 1,NEQN
      P(L) = P(L) + TEMP2*PHI(L,I)
 225      PHI(L,I) = PHI(L,I) + PHI(L,IP1)
 230    CONTINUE
  if ( NORND) go to 240
  DO 235 L = 1,NEQN
    TAU = H*P(L) - PHI(L,15)
    P(L) = Y(L) + TAU
 235    PHI(L,16) = (P(L) - Y(L)) - TAU
  go to 250
 240  DO 245 L = 1,NEQN
 245    P(L) = Y(L) + H*P(L)
 250  XOLD = X
  X = X + H
  ABSH = ABS(H)
  call F(X,P,YP,RPAR,IPAR)
!
!   ESTIMATE ERRORS AT ORDERS K,K-1,K-2
!
  ERKM2 = 0.0
  ERKM1 = 0.0
  ERK = 0.0
  DO 265 L = 1,NEQN
    TEMP3 = 1.0/WT(L)
    TEMP4 = YP(L) - PHI(L,1)
    if ( KM2)265,260,255
 255    ERKM2 = ERKM2 + ((PHI(L,KM1)+TEMP4)*TEMP3)**2
 260    ERKM1 = ERKM1 + ((PHI(L,K)+TEMP4)*TEMP3)**2
 265    ERK = ERK + (TEMP4*TEMP3)**2
  if ( KM2)280,275,270
 270  ERKM2 = ABSH*SIG(KM1)*GSTR(KM2)*SQRT(ERKM2)
 275  ERKM1 = ABSH*SIG(K)*GSTR(KM1)*SQRT(ERKM1)
 280  TEMP5 = ABSH*SQRT(ERK)
  ERR = TEMP5*(G(K)-G(KP1))
  ERK = TEMP5*SIG(KP1)*GSTR(K)
  KNEW = K
!
!   TEST if ORDER SHOULD BE LOWERED
!
  if ( KM2)299,290,285
 285  if ( MAX(ERKM1,ERKM2)  <=  ERK) KNEW = KM1
  go to 299
 290  if ( ERKM1  <=  0.5*ERK) KNEW = KM1
!
!   TEST if STEP SUCCESSFUL
!
 299  if ( ERR  <=  EPS) go to 400
!       ***     END BLOCK 2     ***
!
!       ***     BEGIN BLOCK 3     ***
!   THE STEP IS UNSUCCESSFUL.  RESTORE  X, PHI(*,*), PSI(*) .
!   if THIRD CONSECUTIVE FAILURE, SET ORDER TO ONE.  IF STEP FAILS MORE
!   THAN THREE TIMES, CONSIDER AN OPTIMAL STEP SIZE.  DOUBLE ERROR
!   TOLERANCE AND RETURN if ESTIMATED STEP SIZE IS TOO SMALL FOR MACHINE
!   PRECISION.
!                   ***
!
!   RESTORE X, PHI(*,*) AND PSI(*)
!
  PHASE1 = .FALSE.
  X = XOLD
  DO 310 I = 1,K
    TEMP1 = 1.0/BETA(I)
    IP1 = I+1
    DO 305 L = 1,NEQN
 305      PHI(L,I) = TEMP1*(PHI(L,I) - PHI(L,IP1))
 310    CONTINUE
  if ( K  <  2) go to 320
  DO 315 I = 2,K
 315    PSI(I-1) = PSI(I) - H
!
!   ON THIRD FAILURE, SET ORDER TO ONE.  THEREAFTER, USE OPTIMAL STEP
!   SIZE
!
 320  IFAIL = IFAIL + 1
  TEMP2 = 0.5
  if ( IFAIL - 3) 335,330,325
 325  if ( P5EPS  <  0.25*ERK) TEMP2 = SQRT(P5EPS/ERK)
 330  KNEW = 1
 335  H = TEMP2*H
  K = KNEW
  NS = 0
  if ( ABS(H)  >=  FOURU*ABS(X)) go to 340
  CRASH = .TRUE.
  H = SIGN(FOURU*ABS(X),H)
  EPS = EPS + EPS
  return
 340  go to 100
!       ***     END BLOCK 3     ***
!
!       ***     BEGIN BLOCK 4     ***
!   THE STEP IS SUCCESSFUL.  CORRECT THE PREDICTED SOLUTION, EVALUATE
!   THE DERIVATIVES USING THE CORRECTED SOLUTION AND UPDATE THE
!   DIFFERENCES.  DETERMINE BEST ORDER AND STEP SIZE FOR NEXT STEP.
!                   ***
 400  KOLD = K
  HOLD = H
!
!   CORRECT AND EVALUATE
!
  TEMP1 = H*G(KP1)
  if ( NORND) go to 410
  DO 405 L = 1,NEQN
    TEMP3 = Y(L)
    RHO = TEMP1*(YP(L) - PHI(L,1)) - PHI(L,16)
    Y(L) = P(L) + RHO
    PHI(L,15) = (Y(L) - P(L)) - RHO
 405    P(L) = TEMP3
  go to 420
 410  DO 415 L = 1,NEQN
    TEMP3 = Y(L)
    Y(L) = P(L) + TEMP1*(YP(L) - PHI(L,1))
 415    P(L) = TEMP3
 420  call F(X,Y,YP,RPAR,IPAR)
!
!   UPDATE DIFFERENCES FOR NEXT STEP
!
  DO 425 L = 1,NEQN
    PHI(L,KP1) = YP(L) - PHI(L,1)
 425    PHI(L,KP2) = PHI(L,KP1) - PHI(L,KP2)
  DO 435 I = 1,K
    DO 430 L = 1,NEQN
 430      PHI(L,I) = PHI(L,I) + PHI(L,KP1)
 435    CONTINUE
!
!   ESTIMATE ERROR AT ORDER K+1 UNLESS:
!     IN FIRST PHASE WHEN ALWAYS RAISE ORDER,
!     ALREADY DECIDED TO LOWER ORDER,
!     STEP SIZE NOT CONSTANT SO ESTIMATE UNRELIABLE
!
  ERKP1 = 0.0
  if ( KNEW  ==  KM1  .OR.  K  ==  12) PHASE1 = .FALSE.
  if ( PHASE1) go to 450
  if ( KNEW  ==  KM1) go to 455
  if ( KP1  >  NS) go to 460
  DO 440 L = 1,NEQN
 440    ERKP1 = ERKP1 + (PHI(L,KP2)/WT(L))**2
  ERKP1 = ABSH*GSTR(KP1)*SQRT(ERKP1)
!
!   USING ESTIMATED ERROR AT ORDER K+1, DETERMINE APPROPRIATE ORDER
!   FOR NEXT STEP
!
  if ( K  >  1) go to 445
  if ( ERKP1  >=  0.5*ERK) go to 460
  go to 450
 445  if ( ERKM1  <=  MIN(ERK,ERKP1)) go to 455
  if ( ERKP1  >=  ERK  .OR.  K  ==  12) go to 460
!
!   HERE ERKP1  <  ERK  <  MAX(ERKM1,ERKM2) ELSE ORDER WOULD HAVE
!   BEEN LOWERED IN BLOCK 2.  THUS ORDER IS TO BE RAISED
!
!   RAISE ORDER
!
 450  K = KP1
  ERK = ERKP1
  go to 460
!
!   LOWER ORDER
!
 455  K = KM1
  ERK = ERKM1
!
!   WITH NEW ORDER DETERMINE APPROPRIATE STEP SIZE FOR NEXT STEP
!
 460  HNEW = H + H
  if ( PHASE1) go to 465
  if ( P5EPS  >=  ERK*TWO(K+1)) go to 465
  HNEW = H
  if ( P5EPS  >=  ERK) go to 465
  TEMP2 = K+1
  R = (P5EPS/ERK)**(1.0/TEMP2)
  HNEW = ABSH*MAX(0.5,MIN(0.9,R))
  HNEW = SIGN(MAX(HNEW,FOURU*ABS(X)),H)
 465  H = HNEW
  return
!       ***     END BLOCK 4     ***
end
