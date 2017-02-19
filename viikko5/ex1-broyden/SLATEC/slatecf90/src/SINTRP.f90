subroutine SINTRP (X, Y, XOUT, YOUT, YPOUT, NEQN, KOLD, PHI, IVC, &
     IV, KGI, GI, ALPHA, OG, OW, OX, OY)
!
!! SINTRP approximates the solution at XOUT by evaluating the polynomial ...
!  computed in STEPS at XOUT.  Must be used in conjunction with STEPS.
!
!***LIBRARY   SLATEC (DEPAC)
!***CATEGORY  I1A1B
!***TYPE      SINGLE PRECISION (SINTRP-S, DINTP-D)
!***KEYWORDS  ADAMS METHOD, DEPAC, INITIAL VALUE PROBLEMS, ODE,
!             ORDINARY DIFFERENTIAL EQUATIONS, PREDICTOR-CORRECTOR,
!             SMOOTH INTERPOLANT
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   The methods in subroutine  STEPS  approximate the solution near  X
!   by a polynomial.  Subroutine  SINTRP  approximates the solution at
!   XOUT  by evaluating the polynomial there.  Information defining this
!   polynomial is passed from  STEPS  so  SINTRP  cannot be used alone.
!
!   Subroutine STEPS is completely explained and documented in the text,
!   "Computer Solution of Ordinary Differential Equations, the Initial
!   Value Problem"  by L. F. Shampine and M. K. Gordon.
!
!   Input to SINTRP --
!
!   The user provides storage in the calling program for the arrays in
!   the call list
!      DIMENSION Y(NEQN),YOUT(NEQN),YPOUT(NEQN),PHI(NEQN,16),OY(NEQN)
!                AND ALPHA(12),OG(13),OW(12),GI(11),IV(10)
!   and defines
!      XOUT -- point at which solution is desired.
!   The remaining parameters are defined in  STEPS  and passed to
!   SINTRP  from that subroutine
!
!   Output from  SINTRP --
!
!      YOUT(*) -- solution at  XOUT
!      YPOUT(*) -- derivative of solution at  XOUT
!   The remaining parameters are returned unaltered from their input
!   values.  Integration with  STEPS  may be continued.
!
!***REFERENCES  H. A. Watts, A smoother interpolant for DE/STEP, INTRP
!                 II, Report SAND84-0293, Sandia Laboratories, 1984.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   840201  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SINTRP
!
  DIMENSION Y(*),YOUT(*),YPOUT(*),PHI(NEQN,16),OY(*)
  DIMENSION G(13),C(13),W(13),OG(13),OW(12),ALPHA(12),GI(11),IV(10)
!
!***FIRST EXECUTABLE STATEMENT  SINTRP
  KP1 = KOLD + 1
  KP2 = KOLD + 2
!
  HI = XOUT - OX
  H = X - OX
  XI = HI/H
  XIM1 = XI - 1.
!
!   INITIALIZE W(*) FOR COMPUTING G(*)
!
  XIQ = XI
  DO 10 IQ = 1,KP1
    XIQ = XI*XIQ
    TEMP1 = IQ*(IQ+1)
 10     W(IQ) = XIQ/TEMP1
!
!   COMPUTE THE DOUBLE INTEGRAL TERM GDI
!
  if (KOLD  <=  KGI) go to 50
  if (IVC  >  0) go to 20
  GDI = 1.0/TEMP1
  M = 2
  go to 30
 20   IW = IV(IVC)
  GDI = OW(IW)
  M = KOLD - IW + 3
 30   if (M  >  KOLD) go to 60
  DO 40 I = M,KOLD
 40     GDI = OW(KP2-I) - ALPHA(I)*GDI
  go to 60
 50   GDI = GI(KOLD)
!
!   COMPUTE G(*) AND C(*)
!
 60   G(1) = XI
  G(2) = 0.5*XI*XI
  C(1) = 1.0
  C(2) = XI
  if (KOLD  <  2) go to 90
  DO 80 I = 2,KOLD
    ALP = ALPHA(I)
    GAMMA = 1.0 + XIM1*ALP
    L = KP2 - I
    DO 70 JQ = 1,L
 70       W(JQ) = GAMMA*W(JQ) - ALP*W(JQ+1)
    G(I+1) = W(1)
 80     C(I+1) = GAMMA*C(I)
!
!   DEFINE INTERPOLATION PARAMETERS
!
 90   SIGMA = (W(2) - XIM1*W(1))/GDI
  RMU = XIM1*C(KP1)/GDI
  HMU = RMU/H
!
!   INTERPOLATE FOR THE SOLUTION -- YOUT
!   AND FOR THE DERIVATIVE OF THE SOLUTION -- YPOUT
!
  DO 100 L = 1,NEQN
    YOUT(L) = 0.0
 100    YPOUT(L) = 0.0
  DO 120 J = 1,KOLD
    I = KP2 - J
    GDIF = OG(I) - OG(I-1)
    TEMP2 = (G(I) - G(I-1)) - SIGMA*GDIF
    TEMP3 = (C(I) - C(I-1)) + RMU*GDIF
    DO 110 L = 1,NEQN
      YOUT(L) = YOUT(L) + TEMP2*PHI(L,I)
 110      YPOUT(L) = YPOUT(L) + TEMP3*PHI(L,I)
 120    CONTINUE
  DO 130 L = 1,NEQN
    YOUT(L) = ((1.0 - SIGMA)*OY(L) + SIGMA*Y(L)) + &
               H*(YOUT(L) + (G(1) - SIGMA*OG(1))*PHI(L,1))
 130    YPOUT(L) = HMU*(OY(L) - Y(L)) + &
                  (YPOUT(L) + (C(1) + RMU*OG(1))*PHI(L,1))
!
  return
end
