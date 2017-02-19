subroutine DQK21 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
!
!! DQK21 computes Integral of F over (A,B), with error estimate.
!                       J = Integral of ABS(F) over (A,B)
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A2
!***TYPE      DOUBLE PRECISION (QK21-S, DQK21-D)
!***KEYWORDS  21-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
!                       Declared E X T E R N A L in the driver program.
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
!                       RESULT is computed by applying the 21-POINT
!                       KRONROD RULE (RESK) obtained by optimal addition
!                       of abscissae to the 10-POINT GAUSS RULE (RESG).
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
!***END PROLOGUE  DQK21
!
  DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DHLGTH, &
    D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH,RESABS,RESASC, &
    RESG,RESK,RESKH,RESULT,UFLOW,WG,WGK,XGK
  INTEGER J,JTW,JTWM1
  EXTERNAL F
!
  DIMENSION FV1(10),FV2(10),WG(5),WGK(11),XGK(11)
!
!           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
!           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!           CORRESPONDING WEIGHTS ARE GIVEN.
!
!           XGK    - ABSCISSAE OF THE 21-POINT KRONROD RULE
!                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 10-POINT
!                    GAUSS RULE
!                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!                    ADDED TO THE 10-POINT GAUSS RULE
!
!           WGK    - WEIGHTS OF THE 21-POINT KRONROD RULE
!
!           WG     - WEIGHTS OF THE 10-POINT GAUSS RULE
!
!
! GAUSS QUADRATURE WEIGHTS AND KRONROD QUADRATURE ABSCISSAE AND WEIGHTS
! AS EVALUATED WITH 80 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON,
! BELL LABS, NOV. 1981.
!
  SAVE WG, XGK, WGK
  DATA WG  (  1) / 0.066671344308688137593568809893332D0 /
  DATA WG  (  2) / 0.149451349150580593145776339657697D0 /
  DATA WG  (  3) / 0.219086362515982043995534934228163D0 /
  DATA WG  (  4) / 0.269266719309996355091226921569469D0 /
  DATA WG  (  5) / 0.295524224714752870173892994651338D0 /
!
  DATA XGK (  1) / 0.995657163025808080735527280689003D0 /
  DATA XGK (  2) / 0.973906528517171720077964012084452D0 /
  DATA XGK (  3) / 0.930157491355708226001207180059508D0 /
  DATA XGK (  4) / 0.865063366688984510732096688423493D0 /
  DATA XGK (  5) / 0.780817726586416897063717578345042D0 /
  DATA XGK (  6) / 0.679409568299024406234327365114874D0 /
  DATA XGK (  7) / 0.562757134668604683339000099272694D0 /
  DATA XGK (  8) / 0.433395394129247190799265943165784D0 /
  DATA XGK (  9) / 0.294392862701460198131126603103866D0 /
  DATA XGK ( 10) / 0.148874338981631210884826001129720D0 /
  DATA XGK ( 11) / 0.000000000000000000000000000000000D0 /
!
  DATA WGK (  1) / 0.011694638867371874278064396062192D0 /
  DATA WGK (  2) / 0.032558162307964727478818972459390D0 /
  DATA WGK (  3) / 0.054755896574351996031381300244580D0 /
  DATA WGK (  4) / 0.075039674810919952767043140916190D0 /
  DATA WGK (  5) / 0.093125454583697605535065465083366D0 /
  DATA WGK (  6) / 0.109387158802297641899210590325805D0 /
  DATA WGK (  7) / 0.123491976262065851077958109831074D0 /
  DATA WGK (  8) / 0.134709217311473325928054001771707D0 /
  DATA WGK (  9) / 0.142775938577060080797094273138717D0 /
  DATA WGK ( 10) / 0.147739104901338491374841515972068D0 /
  DATA WGK ( 11) / 0.149445554002916905664936468389821D0 /
!
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           CENTR  - MID POINT OF THE INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           ABSC   - ABSCISSA
!           FVAL*  - FUNCTION VALUE
!           RESG   - RESULT OF THE 10-POINT GAUSS FORMULA
!           RESK   - RESULT OF THE 21-POINT KRONROD FORMULA
!           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
!                    I.E. TO I/(B-A)
!
!
!           MACHINE DEPENDENT CONSTANTS
!           ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  DQK21
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
!
  CENTR = 0.5D+00*(A+B)
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
!
!           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
!           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
  RESG = 0.0D+00
  FC = F(CENTR)
  RESK = WGK(11)*FC
  RESABS = ABS(RESK)
  DO 10 J=1,5
    JTW = 2*J
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
  DO 15 J = 1,5
    JTWM1 = 2*J-1
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
  RESASC = WGK(11)*ABS(FC-RESKH)
  DO 20 J=1,10
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
