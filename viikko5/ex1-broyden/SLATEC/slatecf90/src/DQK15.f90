subroutine DQK15 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
!
!! DQK15 computes Integral of F over (A,B), with error estimate.
!                       J = integral of ABS(F) over (A,B)
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A2
!***TYPE      DOUBLE PRECISION (QK15-S, DQK15-D)
!***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!           Integration rules
!           Standard fortran subroutine
!           Double precision version
!
!           PARAMETERS
!            ON ENTRY
!              F      - Double precision
!                       Function subprogram defining the integrand
!                       FUNCTION F(X). The actual name for F needs to be
!                       Declared E X T E R N A L in the calling program.
!
!              A      - Double precision
!                       Lower limit of integration
!
!              B      - Double precision
!                       Upper limit of integration
!
!            ON RETURN
!              RESULT - Double precision
!                       Approximation to the integral I
!                       Result is computed by applying the 15-POINT
!                       KRONROD RULE (RESK) obtained by optimal addition
!                       of abscissae to the 7-POINT GAUSS RULE(RESG).
!
!              ABSERR - Double precision
!                       Estimate of the modulus of the absolute error,
!                       which should not exceed ABS(I-RESULT)
!
!              RESABS - Double precision
!                       Approximation to the integral J
!
!              RESASC - Double precision
!                       Approximation to the integral of ABS(F-I/(B-A))
!                       over (A,B)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DQK15
!
  DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DHLGTH, &
    D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC, &
    RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
  INTEGER J,JTW,JTWM1
  EXTERNAL F
!
  DIMENSION FV1(7),FV2(7),WG(4),WGK(8),XGK(8)
!
!           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
!           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!           CORRESPONDING WEIGHTS ARE GIVEN.
!
!           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
!                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 7-POINT
!                    GAUSS RULE
!                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!                    ADDED TO THE 7-POINT GAUSS RULE
!
!           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
!
!           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
!
!
! GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
! BELL LABS, NOV. 1981.
!
  SAVE WG, XGK, WGK
  DATA WG  (  1) / 0.129484966168869693270611432679082D0 /
  DATA WG  (  2) / 0.279705391489276667901467771423780D0 /
  DATA WG  (  3) / 0.381830050505118944950369775488975D0 /
  DATA WG  (  4) / 0.417959183673469387755102040816327D0 /
!
  DATA XGK (  1) / 0.991455371120812639206854697526329D0 /
  DATA XGK (  2) / 0.949107912342758524526189684047851D0 /
  DATA XGK (  3) / 0.864864423359769072789712788640926D0 /
  DATA XGK (  4) / 0.741531185599394439863864773280788D0 /
  DATA XGK (  5) / 0.586087235467691130294144838258730D0 /
  DATA XGK (  6) / 0.405845151377397166906606412076961D0 /
  DATA XGK (  7) / 0.207784955007898467600689403773245D0 /
  DATA XGK (  8) / 0.000000000000000000000000000000000D0 /
!
  DATA WGK (  1) / 0.022935322010529224963732008058970D0 /
  DATA WGK (  2) / 0.063092092629978553290700663189204D0 /
  DATA WGK (  3) / 0.104790010322250183839876322541518D0 /
  DATA WGK (  4) / 0.140653259715525918745189590510238D0 /
  DATA WGK (  5) / 0.169004726639267902826583426598550D0 /
  DATA WGK (  6) / 0.190350578064785409913256402421014D0 /
  DATA WGK (  7) / 0.204432940075298892414161999234649D0 /
  DATA WGK (  8) / 0.209482141084727828012999174891714D0 /
!
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           CENTR  - MID POINT OF THE INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           ABSC   - ABSCISSA
!           FVAL*  - FUNCTION VALUE
!           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
!           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
!           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
!                    I.E. TO I/(B-A)
!
!           MACHINE DEPENDENT CONSTANTS
!           ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  DQK15
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
!
  CENTR = 0.5D+00*(A+B)
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
!
!           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
!           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
  FC = F(CENTR)
  RESG = FC*WG(4)
  RESK = FC*WGK(8)
  RESABS = ABS(RESK)
  DO 10 J=1,3
    JTW = J*2
    ABSC = HLGTH*XGK(JTW)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTW) = FVAL1
    FV2(JTW) = FVAL2
    FSUM = FVAL1+FVAL2
    RESG = RESG+WG(J)*FSUM
    RESK = RESK+WGK(JTW)*FSUM
    RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
  DO 15 J = 1,4
    JTWM1 = J*2-1
    ABSC = HLGTH*XGK(JTWM1)
    FVAL1 = F(CENTR-ABSC)
    FVAL2 = F(CENTR+ABSC)
    FV1(JTWM1) = FVAL1
    FV2(JTWM1) = FVAL2
    FSUM = FVAL1+FVAL2
    RESK = RESK+WGK(JTWM1)*FSUM
    RESABS = RESABS+WGK(JTWM1)*(ABS(FVAL1)+ABS(FVAL2))
   15 CONTINUE
  RESKH = RESK*0.5D+00
  RESASC = WGK(8)*ABS(FC-RESKH)
  DO 20 J=1,7
    RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
  RESULT = RESK*HLGTH
  RESABS = RESABS*DHLGTH
  RESASC = RESASC*DHLGTH
  ABSERR = ABS((RESK-RESG)*HLGTH)
  if ( RESASC /= 0.0D+00.AND.ABSERR /= 0.0D+00) &
    ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  if ( RESABS > UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
    ((EPMACH*0.5D+02)*RESABS,ABSERR)
  return
end
