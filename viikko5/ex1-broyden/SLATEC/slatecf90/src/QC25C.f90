subroutine QC25C (F, A, B, C, RESULT, ABSERR, KRUL, NEVAL)
!
!! QC25C computes I = Integral of F*W over (A,B) with error estimate, ...
!  where W(X) = 1/(X-C)
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A2A2, J4
!***TYPE      SINGLE PRECISION (QC25C-S, DQC25C-D)
!***KEYWORDS  25-POINT CLENSHAW-CURTIS INTEGRATION, QUADPACK, QUADRATURE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Integration rules for the computation of CAUCHY
!        PRINCIPAL VALUE integrals
!        Standard fortran subroutine
!        Real version
!
!        PARAMETERS
!           F      - Real
!                    Function subprogram defining the integrand function
!                    F(X). The actual name for F needs to be declared
!                    E X T E R N A L  in the driver program.
!
!           A      - Real
!                    Left end point of the integration interval
!
!           B      - Real
!                    Right end point of the integration interval, B > A
!
!           C      - Real
!                    Parameter in the WEIGHT function
!
!           RESULT - Real
!                    Approximation to the integral
!                    result is computed by using a generalized
!                    Clenshaw-Curtis method if C lies within ten percent
!                    of the integration interval. In the other case the
!                    15-point Kronrod rule obtained by optimal addition
!                    of abscissae to the 7-point Gauss rule, is applied.
!
!           ABSERR - Real
!                    Estimate of the modulus of the absolute error,
!                    which should equal or exceed ABS(I-RESULT)
!
!           KRUL   - Integer
!                    Key which is decreased by 1 if the 15-point
!                    Gauss-Kronrod scheme has been used
!
!           NEVAL  - Integer
!                    Number of integrand evaluations
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  QCHEB, QK15W, QWGTC
!***REVISION HISTORY  (YYMMDD)
!   810101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  QC25C
!
  REAL A,ABSERR,AK22,AMOM0,AMOM1,AMOM2,B,C,CC, &
    CENTR,CHEB12,CHEB24,QWGTC,F,FVAL,HLGTH,P2,P3,P4, &
    RESABS,RESASC,RESULT,RES12,RES24,U,X
  INTEGER I,ISYM,K,KP,KRUL,NEVAL
!
  DIMENSION X(11),FVAL(25),CHEB12(13),CHEB24(25)
!
  EXTERNAL F, QWGTC
!
!           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24),
!           K = 1, ..., 11, TO BE USED FOR THE CHEBYSHEV SERIES
!           EXPANSION OF F
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
!           ----------------------
!           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
!                    COS(K*PI/24),  K = 0, ..., 24
!           CHEB12 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS,
!                    FOR THE FUNCTION F, OF DEGREE 12
!           CHEB24 - CHEBYSHEV SERIES EXPANSION COEFFICIENTS,
!                    FOR THE FUNCTION F, OF DEGREE 24
!           RES12  - APPROXIMATION TO THE INTEGRAL CORRESPONDING
!                    TO THE USE OF CHEB12
!           RES24  - APPROXIMATION TO THE INTEGRAL CORRESPONDING
!                    TO THE USE OF CHEB24
!           QWGTC - EXTERNAL FUNCTION SUBPROGRAM DEFINING
!                    THE WEIGHT FUNCTION
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           CENTR  - MID POINT OF THE INTERVAL
!
!
!           CHECK THE POSITION OF C.
!
!***FIRST EXECUTABLE STATEMENT  QC25C
  CC = (0.2E+01*C-B-A)/(B-A)
  if ( ABS(CC) < 0.11E+01) go to 10
!
!           APPLY THE 15-POINT GAUSS-KRONROD SCHEME.
!
  KRUL = KRUL-1
  call QK15W(F,QWGTC,C,P2,P3,P4,KP,A,B,RESULT,ABSERR, &
    RESABS,RESASC)
  NEVAL = 15
  if (RESASC == ABSERR) KRUL = KRUL+1
  go to 50
!
!           USE THE GENERALIZED CLENSHAW-CURTIS METHOD.
!
   10 HLGTH = 0.5E+00*(B-A)
  CENTR = 0.5E+00*(B+A)
  NEVAL = 25
  FVAL(1) = 0.5E+00*F(HLGTH+CENTR)
  FVAL(13) = F(CENTR)
  FVAL(25) = 0.5E+00*F(CENTR-HLGTH)
  DO 20 I=2,12
    U = HLGTH*X(I-1)
    ISYM = 26-I
    FVAL(I) = F(U+CENTR)
    FVAL(ISYM) = F(CENTR-U)
   20 CONTINUE
!
!           COMPUTE THE CHEBYSHEV SERIES EXPANSION.
!
  call QCHEB(X,FVAL,CHEB12,CHEB24)
!
!           THE MODIFIED CHEBYSHEV MOMENTS ARE COMPUTED
!           BY FORWARD RECURSION, USING AMOM0 AND AMOM1
!           AS STARTING VALUES.
!
  AMOM0 = LOG(ABS((0.1E+01-CC)/(0.1E+01+CC)))
  AMOM1 = 0.2E+01+CC*AMOM0
  RES12 = CHEB12(1)*AMOM0+CHEB12(2)*AMOM1
  RES24 = CHEB24(1)*AMOM0+CHEB24(2)*AMOM1
  DO 30 K=3,13
    AMOM2 = 0.2E+01*CC*AMOM1-AMOM0
    AK22 = (K-2)*(K-2)
    if ( (K/2)*2 == K) AMOM2 = AMOM2-0.4E+01/(AK22-0.1E+01)
    RES12 = RES12+CHEB12(K)*AMOM2
    RES24 = RES24+CHEB24(K)*AMOM2
    AMOM0 = AMOM1
    AMOM1 = AMOM2
   30 CONTINUE
  DO 40 K=14,25
    AMOM2 = 0.2E+01*CC*AMOM1-AMOM0
    AK22 = (K-2)*(K-2)
    if ( (K/2)*2 == K) AMOM2 = AMOM2-0.4E+01/ &
    (AK22-0.1E+01)
    RES24 = RES24+CHEB24(K)*AMOM2
    AMOM0 = AMOM1
    AMOM1 = AMOM2
   40 CONTINUE
  RESULT = RES24
  ABSERR = ABS(RES24-RES12)
   50 RETURN
end
