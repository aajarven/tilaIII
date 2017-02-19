subroutine HSTART (F, NEQ, A, B, Y, YPRIME, ETOL, MORDER, SMALL, &
     BIG, SPY, PV, YP, SF, RPAR, IPAR, H)
!
!! HSTART is subsidiary to DEABM, DEBDF and DERKF.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (HSTART-S, DHSTRT-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   HSTART computes a starting step size to be used in solving initial
!   value problems in ordinary differential equations.
! **********************************************************************
!  Abstract
!
!     Subroutine HSTART computes a starting step size to be used by an
!     initial value method in solving ordinary differential equations.
!     It is based on an estimate of the local Lipschitz constant for the
!     differential equation (lower bound on a norm of the Jacobian),
!     a bound on the differential equation (first derivative), and
!     a bound on the partial derivative of the equation with respect to
!     the independent variable.
!     (All approximated near the initial point A.)
!
!     Subroutine HSTART uses a function subprogram HVNRM for computing
!     a vector norm.  The maximum norm is presently utilized though it
!     can easily be replaced by any other vector norm.  It is presumed
!     that any replacement norm routine would be carefully coded to
!     prevent unnecessary underflows or overflows from occurring, and
!     also, would not alter the vector or number of components.
!
! **********************************************************************
!  On Input you must provide the following
!
!      F -- This is a subroutine of the form
!                               F(X,U,UPRIME,RPAR,IPAR)
!             which defines the system of first order differential
!             equations to be solved.  For the given values of X and the
!             vector  U(*)=(U(1),U(2),...,U(NEQ)) , the subroutine must
!             evaluate the NEQ components of the system of differential
!             equations  dU/DX=F(X,U)  and store the derivatives in the
!             array UPRIME(*), that is,  UPRIME(I) = * dU(I)/DX *  for
!             equations I=1,...,NEQ.
!
!             Subroutine F must not alter X or U(*).  You must declare
!             the name F in an EXTERNAL statement in your program that
!             calls HSTART.  You must dimension U and UPRIME in F.
!
!             RPAR and IPAR are real and integer parameter arrays which
!             you can use for communication between your program and
!             subroutine F.  They are not used or altered by HSTART.  If
!             you do not need RPAR or IPAR, ignore these parameters by
!             treating them as dummy arguments.  If you do choose to use
!             them, dimension them in your program and in F as arrays
!             of appropriate length.
!
!      NEQ -- This is the number of (first order) differential equations
!             to be integrated.
!
!      A -- This is the initial point of integration.
!
!      B -- This is a value of the independent variable used to define
!             the direction of integration.  A reasonable choice is to
!             set  B  to the first point at which a solution is desired.
!             You can also use  B, if necessary, to restrict the length
!             of the first integration step because the algorithm will
!             not compute a starting step length which is bigger than
!             ABS(B-A), unless  B  has been chosen too close to  A.
!             (It is presumed that HSTART has been called with  B
!             different from  A  on the machine being used.  Also see
!             the discussion about the parameter  SMALL.)
!
!      Y(*) -- This is the vector of initial values of the NEQ solution
!             components at the initial point  A.
!
!      YPRIME(*) -- This is the vector of derivatives of the NEQ
!             solution components at the initial point  A.
!             (defined by the differential equations in subroutine F)
!
!      ETOL -- This is the vector of error tolerances corresponding to
!             the NEQ solution components.  It is assumed that all
!             elements are positive.  Following the first integration
!             step, the tolerances are expected to be used by the
!             integrator in an error test which roughly requires that
!                        ABS(local error)  <=  ETOL
!             for each vector component.
!
!      MORDER -- This is the order of the formula which will be used by
!             the initial value method for taking the first integration
!             step.
!
!      SMALL -- This is a small positive machine dependent constant
!             which is used for protecting against computations with
!             numbers which are too small relative to the precision of
!             floating point arithmetic.  SMALL  should be set to
!             (approximately) the smallest positive real number such
!             that  (1.+SMALL)  >  1.  on the machine being used. the
!             quantity  SMALL**(3/8)  is used in computing increments of
!             variables for approximating derivatives by differences.
!             also the algorithm will not compute a starting step length
!             which is smaller than  100*SMALL*ABS(A).
!
!      BIG -- This is a large positive machine dependent constant which
!             is used for preventing machine overflows.  A reasonable
!             choice is to set big to (approximately) the square root of
!             the largest real number which can be held in the machine.
!
!      SPY(*),PV(*),YP(*),SF(*) -- These are real work arrays of length
!             NEQ which provide the routine with needed storage space.
!
!      RPAR,IPAR -- These are parameter arrays, of real and integer
!             type, respectively, which can be used for communication
!             between your program and the F subroutine.  They are not
!             used or altered by HSTART.
!
! **********************************************************************
!  On Output  (after the return from HSTART),
!
!      H -- Is an appropriate starting step size to be attempted by the
!             differential equation method.
!
!           All parameters in the call list remain unchanged except for
!           the working arrays SPY(*),PV(*),YP(*) and SF(*).
!
! **********************************************************************
!
!***SEE ALSO  DEABM, DEBDF, DERKF
!***ROUTINES CALLED  HVNRM
!***REVISION HISTORY  (YYMMDD)
!   800501  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891024  Changed references from VNORM to HVNRM.  (WRB)
!   891024  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  HSTART
!
  DIMENSION Y(*),YPRIME(*),ETOL(*),SPY(*),PV(*),YP(*),SF(*), &
     RPAR(*),IPAR(*)
  EXTERNAL F
!
!.......................................................................
!
!***FIRST EXECUTABLE STATEMENT  HSTART
  DX = B - A
  ABSDX = ABS(DX)
  RELPER = SMALL**0.375
  YNORM = HVNRM(Y,NEQ)
!
!.......................................................................
!
!     COMPUTE A WEIGHTED APPROXIMATE BOUND (DFDXB) ON THE PARTIAL
!     DERIVATIVE OF THE EQUATION WITH RESPECT TO THE
!     INDEPENDENT VARIABLE. PROTECT AGAINST AN OVERFLOW. ALSO
!     COMPUTE A WEIGHTED BOUND (FBND) ON THE FIRST DERIVATIVE LOCALLY.
!
  DA = SIGN(MAX(MIN(RELPER*ABS(A),ABSDX),100.*SMALL*ABS(A)),DX)
  if (DA  ==  0.) DA = RELPER*DX
  call F(A+DA,Y,SF,RPAR,IPAR)
!
  if (MORDER  ==  1) go to 20
  POWER = 2./(MORDER+1)
  DO 10 J=1,NEQ
    WTJ = ETOL(J)**POWER
    SPY(J) = SF(J)/WTJ
    YP(J) = YPRIME(J)/WTJ
   10   PV(J) = SPY(J) - YP(J)
  go to 40
!
   20 DO 30 J=1,NEQ
    SPY(J) = SF(J)/ETOL(J)
    YP(J) = YPRIME(J)/ETOL(J)
   30   PV(J) = SPY(J) - YP(J)
!
   40 DELF = HVNRM(PV,NEQ)
  DFDXB = BIG
  if (DELF  <  BIG*ABS(DA)) DFDXB = DELF/ABS(DA)
  YPNORM = HVNRM(YP,NEQ)
  FBND = MAX(HVNRM(SPY,NEQ),YPNORM)
!
!.......................................................................
!
!     COMPUTE AN ESTIMATE (DFDUB) OF THE LOCAL LIPSCHITZ CONSTANT FOR
!     THE SYSTEM OF DIFFERENTIAL EQUATIONS. THIS ALSO REPRESENTS AN
!     ESTIMATE OF THE NORM OF THE JACOBIAN LOCALLY.
!     THREE ITERATIONS (TWO WHEN NEQ=1) ARE USED TO ESTIMATE THE
!     LIPSCHITZ CONSTANT BY NUMERICAL DIFFERENCES. THE FIRST
!     PERTURBATION VECTOR IS BASED ON THE INITIAL DERIVATIVES AND
!     DIRECTION OF INTEGRATION. THE SECOND PERTURBATION VECTOR IS
!     FORMED USING ANOTHER EVALUATION OF THE DIFFERENTIAL EQUATION.
!     THE THIRD PERTURBATION VECTOR IS FORMED USING PERTURBATIONS BASED
!     ONLY ON THE INITIAL VALUES. COMPONENTS THAT ARE ZERO ARE ALWAYS
!     CHANGED TO NON-ZERO VALUES (EXCEPT ON THE FIRST ITERATION). WHEN
!     INFORMATION IS AVAILABLE, CARE IS TAKEN TO ENSURE THAT COMPONENTS
!     OF THE PERTURBATION VECTOR HAVE SIGNS WHICH ARE CONSISTENT WITH
!     THE SLOPES OF LOCAL SOLUTION CURVES.
!     ALSO CHOOSE THE LARGEST BOUND (FBND) FOR THE FIRST DERIVATIVE.
!     NO ATTEMPT IS MADE TO KEEP THE PERTURBATION VECTOR SIZE CONSTANT.
!
  if (YPNORM  ==  0.) go to 60
!                       USE INITIAL DERIVATIVES FOR FIRST PERTURBATION
  ICASE = 1
  DO 50 J=1,NEQ
    SPY(J) = YPRIME(J)
   50   YP(J) = YPRIME(J)
  go to 80
!                       CANNOT HAVE A NULL PERTURBATION VECTOR
   60 ICASE = 2
  DO 70 J=1,NEQ
    SPY(J) = YPRIME(J)
   70   YP(J) = ETOL(J)
!
   80 DFDUB = 0.
  LK = MIN(NEQ+1,3)
  DO 260 K=1,LK
!                       SET YPNORM AND DELX
    YPNORM = HVNRM(YP,NEQ)
    if (ICASE  ==  1  .OR.  ICASE  ==  3) go to 90
    DELX = SIGN(1.0,DX)
    go to 120
!                       TRY TO ENFORCE MEANINGFUL PERTURBATION VALUES
   90   DELX = DX
    if (ABS(DELX)*YPNORM  >=  RELPER*YNORM) go to 100
    DELXB = BIG
    if (RELPER*YNORM  <  BIG*YPNORM) DELXB = RELPER*YNORM/YPNORM
    DELX = SIGN(DELXB,DX)
  100   DO 110 J=1,NEQ
      if (ABS(DELX*YP(J))  >  ETOL(J)) DELX=SIGN(ETOL(J)/YP(J),DX)
  110     CONTINUE
!                       DEFINE PERTURBED VECTOR OF INITIAL VALUES
  120   DO 130 J=1,NEQ
  130     PV(J) = Y(J) + DELX*YP(J)
    if (K  ==  2) go to 150
!                       EVALUATE DERIVATIVES ASSOCIATED WITH PERTURBED
!                       VECTOR  AND  COMPUTE CORRESPONDING DIFFERENCES
    call F(A,PV,YP,RPAR,IPAR)
    DO 140 J=1,NEQ
  140     PV(J) = YP(J) - YPRIME(J)
    go to 170
!                       USE A SHIFTED VALUE OF THE INDEPENDENT VARIABLE
!                                             IN COMPUTING ONE ESTIMATE
  150   call F(A+DA,PV,YP,RPAR,IPAR)
    DO 160 J=1,NEQ
  160     PV(J) = YP(J) - SF(J)
!                       CHOOSE LARGEST BOUND ON THE WEIGHTED FIRST
!                                                   DERIVATIVE
  170   if (MORDER  ==  1) go to 190
    DO 180 J=1,NEQ
  180     YP(J) = YP(J)/ETOL(J)**POWER
    go to 210
  190   DO 200 J=1,NEQ
  200     YP(J) = YP(J)/ETOL(J)
  210   FBND = MAX(FBND,HVNRM(YP,NEQ))
!                       COMPUTE BOUND ON A LOCAL LIPSCHITZ CONSTANT
    DELF = HVNRM(PV,NEQ)
    if (DELF  ==  0.) go to 220
    DELY = ABS(DELX)*YPNORM
    if (DELF  >=  BIG*DELY) go to 270
    DFDUB = MAX(DFDUB,DELF/DELY)
!
  220   if (K  ==  LK) go to 280
!                       CHOOSE NEXT PERTURBATION VECTOR
    DO 250 J=1,NEQ
      if (K  ==  LK-1) go to 230
      ICASE = 3
      DY = ABS(PV(J))
      if (DY  ==  0.) DY = MAX(DELF,ETOL(J))
      go to 240
  230     ICASE = 4
      DY = MAX(RELPER*ABS(Y(J)),ETOL(J))
  240     if (SPY(J)  ==  0.) SPY(J) = YP(J)
      if (SPY(J)  /=  0.) DY = SIGN(DY,SPY(J))
  250     YP(J) = DY
  260   CONTINUE
!
!                       PROTECT AGAINST AN OVERFLOW
  270 DFDUB = BIG
!
!.......................................................................
!
!     COMPUTE A BOUND (YDPB) ON THE NORM OF THE SECOND DERIVATIVE
!
  280 YDPB = DFDXB + DFDUB*FBND
!
!.......................................................................
!
!     COMPUTE A STARTING STEP SIZE BASED ON THE ABOVE FIRST AND SECOND
!     DERIVATIVE INFORMATION
!
!                       RESTRICT THE STEP LENGTH TO BE NOT BIGGER THAN
!                       ABS(B-A).   (UNLESS  B  IS TOO CLOSE TO  A)
  H = ABSDX
!
  if (YDPB  /=  0.  .OR.  FBND  /=  0.) go to 290
!
!                       BOTH FIRST DERIVATIVE TERM (FBND) AND SECOND
!                                    DERIVATIVE TERM (YDPB) ARE ZERO
  go to 310
!
  290 if (YDPB  /=  0.) go to 300
!
!                       ONLY SECOND DERIVATIVE TERM (YDPB) IS ZERO
  if (1.0  <  FBND*ABSDX) H = 1./FBND
  go to 310
!
!                       SECOND DERIVATIVE TERM (YDPB) IS NON-ZERO
  300 SRYDPB = SQRT(0.5*YDPB)
  if (1.0  <  SRYDPB*ABSDX) H = 1./SRYDPB
!
!                       FURTHER RESTRICT THE STEP LENGTH TO BE NOT
!                                                 BIGGER THAN  1/DFDUB
  310 if (H*DFDUB  >  1.) H = 1./DFDUB
!
!                       FINALLY, RESTRICT THE STEP LENGTH TO BE NOT
!                       SMALLER THAN  100*SMALL*ABS(A).  HOWEVER, IF
!                       A=0. AND THE COMPUTED H UNDERFLOWED TO ZERO,
!                       THE ALGORITHM RETURNS  SMALL*ABS(B)  FOR THE
!                                                       STEP LENGTH.
  H = MAX(H,100.*SMALL*ABS(A))
  if (H  ==  0.) H = SMALL*ABS(B)
!
!                       NOW SET DIRECTION OF INTEGRATION
  H = SIGN(H,DX)
!
  return
end
