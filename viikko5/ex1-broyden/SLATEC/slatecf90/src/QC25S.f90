subroutine QC25S (F, A, B, BL, BR, ALFA, BETA, RI, RJ, RG, RH, &
     RESULT, ABSERR, RESASC, INTEGR, NEV)
!
!! QC25S computes I = Integral of F*W over (BL,BR), with error estimate, ...
!  where the weight function W has a singular
!            behaviour of ALGEBRAICO-LOGARITHMIC type at the points
!            A and/or B. (BL,BR) is a part of (A,B).
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A2A2
!***TYPE      SINGLE PRECISION (QC25S-S, DQC25S-D)
!***KEYWORDS  25-POINT CLENSHAW-CURTIS INTEGRATION, QUADPACK, QUADRATURE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Integration rules for integrands having ALGEBRAICO-LOGARITHMIC
!        end point singularities
!        Standard fortran subroutine
!        Real version
!
!        PARAMETERS
!           F      - Real
!                    Function subprogram defining the integrand
!                    F(X). The actual name for F needs to be declared
!                    E X T E R N A L  in the driver program.
!
!           A      - Real
!                    Left end point of the original interval
!
!           B      - Real
!                    Right end point of the original interval, B > A
!
!           BL     - Real
!                    Lower limit of integration, BL >= A
!
!           BR     - Real
!                    Upper limit of integration, BR <= B
!
!           ALFA   - Real
!                    PARAMETER IN THE WEIGHT FUNCTION
!
!           BETA   - Real
!                    Parameter in the weight function
!
!           RI,RJ,RG,RH - Real
!                    Modified CHEBYSHEV moments for the application
!                    of the generalized CLENSHAW-CURTIS
!                    method (computed in subroutine DQMOMO)
!
!           RESULT - Real
!                    Approximation to the integral
!                    RESULT is computed by using a generalized
!                    CLENSHAW-CURTIS method if B1 = A or BR = B.
!                    in all other cases the 15-POINT KRONROD
!                    RULE is applied, obtained by optimal addition of
!                    Abscissae to the 7-POINT GAUSS RULE.
!
!           ABSERR - Real
!                    Estimate of the modulus of the absolute error,
!                    which should equal or exceed ABS(I-RESULT)
!
!           RESASC - Real
!                    Approximation to the integral of ABS(F*W-I/(B-A))
!
!           INTEGR - Integer
!                    Which determines the weight function
!                    = 1   W(X) = (X-A)**ALFA*(B-X)**BETA
!                    = 2   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
!                    = 3   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
!                    = 4   W(X) = (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*
!                                 LOG(B-X)
!
!           NEV    - Integer
!                    Number of integrand evaluations
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  QCHEB, QK15W, QWGTS
!***REVISION HISTORY  (YYMMDD)
!   810101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  QC25S
!
  REAL A,ABSERR,ALFA,B,BETA,BL,BR,CENTR,CHEB12,CHEB24, &
    DC,F,FACTOR,FIX,FVAL,HLGTH,RESABS,RESASC, &
    RESULT,RES12,RES24,RG,RH,RI,RJ,U,QWGTS,X
  INTEGER I,INTEGR,ISYM,NEV
!
  DIMENSION CHEB12(13),CHEB24(25),FVAL(25),RG(25),RH(25),RI(25), &
    RJ(25),X(11)
!
  EXTERNAL F, QWGTS
!
!           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24)
!           K = 1, ..., 11, TO BE USED FOR THE COMPUTATION OF THE
!           CHEBYSHEV SERIES EXPANSION OF F.
!
  SAVE X
  DATA X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9),X(10), &
    X(11)/ &
       0.9914448613738104E+00,     0.9659258262890683E+00, &
       0.9238795325112868E+00,     0.8660254037844386E+00, &
       0.7933533402912352E+00,     0.7071067811865475E+00, &
       0.6087614290087206E+00,     0.5000000000000000E+00, &
       0.3826834323650898E+00,     0.2588190451025208E+00, &
       0.1305261922200516E+00/
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
!                    (BR-BL)*0.5*COS(K*PI/24)+(BR+BL)*0.5
!                    K = 0, ..., 24
!           CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
!                    OF DEGREE 12, FOR THE FUNCTION F, IN THE
!                    INTERVAL (BL,BR)
!           CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
!                    OF DEGREE 24, FOR THE FUNCTION F, IN THE
!                    INTERVAL (BL,BR)
!           RES12  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB12
!           RES24  - APPROXIMATION TO THE INTEGRAL OBTAINED FROM CHEB24
!           QWGTS - EXTERNAL FUNCTION SUBPROGRAM DEFINING
!                    THE FOUR POSSIBLE WEIGHT FUNCTIONS
!           HLGTH  - HALF-LENGTH OF THE INTERVAL (BL,BR)
!           CENTR  - MID POINT OF THE INTERVAL (BL,BR)
!
!***FIRST EXECUTABLE STATEMENT  QC25S
  NEV = 25
  if ( BL == A.AND.(ALFA /= 0.0E+00.OR.INTEGR == 2.OR.INTEGR == 4)) &
   go to 10
  if ( BR == B.AND.(BETA /= 0.0E+00.OR.INTEGR == 3.OR.INTEGR == 4)) &
   go to 140
!
!           if A > BL AND B < BR, APPLY THE 15-POINT GAUSS-KRONROD
!           SCHEME.
!
!
  call QK15W(F,QWGTS,A,B,ALFA,BETA,INTEGR,BL,BR, &
      RESULT,ABSERR,RESABS,RESASC)
  NEV = 15
  go to 270
!
!           THIS PART OF THE PROGRAM IS EXECUTED ONLY if A = BL.
!           ----------------------------------------------------
!
!           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
!           FOLLOWING FUNCTION
!           F1 = (0.5*(B+B-BR-A)-0.5*(BR-A)*X)**BETA
!                  *F(0.5*(BR-A)*X+0.5*(BR+A))
!
   10 HLGTH = 0.5E+00*(BR-BL)
  CENTR = 0.5E+00*(BR+BL)
  FIX = B-CENTR
  FVAL(1) = 0.5E+00*F(HLGTH+CENTR)*(FIX-HLGTH)**BETA
  FVAL(13) = F(CENTR)*(FIX**BETA)
  FVAL(25) = 0.5E+00*F(CENTR-HLGTH)*(FIX+HLGTH)**BETA
  DO 20 I=2,12
    U = HLGTH*X(I-1)
    ISYM = 26-I
    FVAL(I) = F(U+CENTR)*(FIX-U)**BETA
    FVAL(ISYM) = F(CENTR-U)*(FIX+U)**BETA
   20 CONTINUE
  FACTOR = HLGTH**(ALFA+0.1E+01)
  RESULT = 0.0E+00
  ABSERR = 0.0E+00
  RES12 = 0.0E+00
  RES24 = 0.0E+00
  if ( INTEGR > 2) go to 70
  call QCHEB(X,FVAL,CHEB12,CHEB24)
!
!           INTEGR = 1  (OR 2)
!
  DO 30 I=1,13
    RES12 = RES12+CHEB12(I)*RI(I)
    RES24 = RES24+CHEB24(I)*RI(I)
   30 CONTINUE
  DO 40 I=14,25
    RES24 = RES24+CHEB24(I)*RI(I)
   40 CONTINUE
  if ( INTEGR == 1) go to 130
!
!           INTEGR = 2
!
  DC = LOG(BR-BL)
  RESULT = RES24*DC
  ABSERR = ABS((RES24-RES12)*DC)
  RES12 = 0.0E+00
  RES24 = 0.0E+00
  DO 50 I=1,13
    RES12 = RES12+CHEB12(I)*RG(I)
    RES24 = RES12+CHEB24(I)*RG(I)
   50 CONTINUE
  DO 60 I=14,25
    RES24 = RES24+CHEB24(I)*RG(I)
   60 CONTINUE
  go to 130
!
!           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
!           FOLLOWING FUNCTION
!           F4 = F1*LOG(0.5*(B+B-BR-A)-0.5*(BR-A)*X)
!
   70 FVAL(1) = FVAL(1)*LOG(FIX-HLGTH)
  FVAL(13) = FVAL(13)*LOG(FIX)
  FVAL(25) = FVAL(25)*LOG(FIX+HLGTH)
  DO 80 I=2,12
    U = HLGTH*X(I-1)
    ISYM = 26-I
    FVAL(I) = FVAL(I)*LOG(FIX-U)
    FVAL(ISYM) = FVAL(ISYM)*LOG(FIX+U)
   80 CONTINUE
  call QCHEB(X,FVAL,CHEB12,CHEB24)
!
!           INTEGR = 3  (OR 4)
!
  DO 90 I=1,13
    RES12 = RES12+CHEB12(I)*RI(I)
    RES24 = RES24+CHEB24(I)*RI(I)
   90 CONTINUE
  DO 100 I=14,25
    RES24 = RES24+CHEB24(I)*RI(I)
  100 CONTINUE
  if ( INTEGR == 3) go to 130
!
!           INTEGR = 4
!
  DC = LOG(BR-BL)
  RESULT = RES24*DC
  ABSERR = ABS((RES24-RES12)*DC)
  RES12 = 0.0E+00
  RES24 = 0.0E+00
  DO 110 I=1,13
    RES12 = RES12+CHEB12(I)*RG(I)
    RES24 = RES24+CHEB24(I)*RG(I)
  110 CONTINUE
  DO 120 I=14,25
    RES24 = RES24+CHEB24(I)*RG(I)
  120 CONTINUE
  130 RESULT = (RESULT+RES24)*FACTOR
  ABSERR = (ABSERR+ABS(RES24-RES12))*FACTOR
  go to 270
!
!           THIS PART OF THE PROGRAM IS EXECUTED ONLY if B = BR.
!           ----------------------------------------------------
!
!           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
!           FOLLOWING FUNCTION
!           F2 = (0.5*(B+BL-A-A)+0.5*(B-BL)*X)**ALFA
!                *F(0.5*(B-BL)*X+0.5*(B+BL))
!
  140 HLGTH = 0.5E+00*(BR-BL)
  CENTR = 0.5E+00*(BR+BL)
  FIX = CENTR-A
  FVAL(1) = 0.5E+00*F(HLGTH+CENTR)*(FIX+HLGTH)**ALFA
  FVAL(13) = F(CENTR)*(FIX**ALFA)
  FVAL(25) = 0.5E+00*F(CENTR-HLGTH)*(FIX-HLGTH)**ALFA
  DO 150 I=2,12
    U = HLGTH*X(I-1)
    ISYM = 26-I
    FVAL(I) = F(U+CENTR)*(FIX+U)**ALFA
    FVAL(ISYM) = F(CENTR-U)*(FIX-U)**ALFA
  150 CONTINUE
  FACTOR = HLGTH**(BETA+0.1E+01)
  RESULT = 0.0E+00
  ABSERR = 0.0E+00
  RES12 = 0.0E+00
  RES24 = 0.0E+00
  if ( INTEGR == 2.OR.INTEGR == 4) go to 200
!
!           INTEGR = 1  (OR 3)
!
  call QCHEB(X,FVAL,CHEB12,CHEB24)
  DO 160 I=1,13
    RES12 = RES12+CHEB12(I)*RJ(I)
    RES24 = RES24+CHEB24(I)*RJ(I)
  160 CONTINUE
  DO 170 I=14,25
    RES24 = RES24+CHEB24(I)*RJ(I)
  170 CONTINUE
  if ( INTEGR == 1) go to 260
!
!           INTEGR = 3
!
  DC = LOG(BR-BL)
  RESULT = RES24*DC
  ABSERR = ABS((RES24-RES12)*DC)
  RES12 = 0.0E+00
  RES24 = 0.0E+00
  DO 180 I=1,13
    RES12 = RES12+CHEB12(I)*RH(I)
    RES24 = RES24+CHEB24(I)*RH(I)
  180 CONTINUE
  DO 190 I=14,25
    RES24 = RES24+CHEB24(I)*RH(I)
  190 CONTINUE
  go to 260
!
!           COMPUTE THE CHEBYSHEV SERIES EXPANSION OF THE
!           FOLLOWING FUNCTION
!           F3 = F2*LOG(0.5*(B-BL)*X+0.5*(B+BL-A-A))
!
  200 FVAL(1) = FVAL(1)*LOG(FIX+HLGTH)
  FVAL(13) = FVAL(13)*LOG(FIX)
  FVAL(25) = FVAL(25)*LOG(FIX-HLGTH)
  DO 210 I=2,12
    U = HLGTH*X(I-1)
    ISYM = 26-I
    FVAL(I) = FVAL(I)*LOG(FIX+U)
    FVAL(ISYM) = FVAL(ISYM)*LOG(FIX-U)
  210 CONTINUE
  call QCHEB(X,FVAL,CHEB12,CHEB24)
!
!           INTEGR = 2  (OR 4)
!
  DO 220 I=1,13
    RES12 = RES12+CHEB12(I)*RJ(I)
    RES24 = RES24+CHEB24(I)*RJ(I)
  220 CONTINUE
  DO 230 I=14,25
    RES24 = RES24+CHEB24(I)*RJ(I)
  230 CONTINUE
  if ( INTEGR == 2) go to 260
  DC = LOG(BR-BL)
  RESULT = RES24*DC
  ABSERR = ABS((RES24-RES12)*DC)
  RES12 = 0.0E+00
  RES24 = 0.0E+00
!
!           INTEGR = 4
!
  DO 240 I=1,13
    RES12 = RES12+CHEB12(I)*RH(I)
    RES24 = RES24+CHEB24(I)*RH(I)
  240 CONTINUE
  DO 250 I=14,25
    RES24 = RES24+CHEB24(I)*RH(I)
  250 CONTINUE
  260 RESULT = (RESULT+RES24)*FACTOR
  ABSERR = (ABSERR+ABS(RES24-RES12))*FACTOR
  270 RETURN
end
