subroutine QK21 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
!
!! QK21 estimates an integral with a 21 point Gauss Kronrod rule.
!
!***PURPOSE  To compute I = Integral of F over (A,B), with error
!                           estimate
!                       J = Integral of ABS(F) over (A,B)
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A2
!***TYPE      SINGLE PRECISION (QK21-S, DQK21-D)
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
!           Real version
!
!           PARAMETERS
!            ON ENTRY
!              F      - Real
!                       Function subprogram defining the integrand
!                       FUNCTION F(X). The actual name for F needs to be
!                       Declared E X T E R N A L in the driver program.
!
!              A      - Real
!                       Lower limit of integration
!
!              B      - Real
!                       Upper limit of integration
!
!            ON RETURN
!              RESULT - Real
!                       Approximation to the integral I
!                       RESULT is computed by applying the 21-POINT
!                       KRONROD RULE (RESK) obtained by optimal addition
!                       of abscissae to the 10-POINT GAUSS RULE (RESG).
!
!              ABSERR - Real
!                       Estimate of the modulus of the absolute error,
!                       which should not exceed ABS(I-RESULT)
!
!              RESABS - Real
!                       Approximation to the integral J
!
!              RESASC - Real
!                       Approximation to the integral of ABS(F-I/(B-A))
!                       over (A,B)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  QK21
!
  REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2, &
    FV1,FV2,HLGTH,RESABS,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW,WG,WGK, &
    XGK
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
  SAVE XGK, WGK, WG
  DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7), &
    XGK(8),XGK(9),XGK(10),XGK(11)/ &
           0.9956571630258081E+00,     0.9739065285171717E+00, &
       0.9301574913557082E+00,     0.8650633666889845E+00, &
       0.7808177265864169E+00,     0.6794095682990244E+00, &
       0.5627571346686047E+00,     0.4333953941292472E+00, &
       0.2943928627014602E+00,     0.1488743389816312E+00, &
       0.0000000000000000E+00/
!
  DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7), &
    WGK(8),WGK(9),WGK(10),WGK(11)/ &
       0.1169463886737187E-01,     0.3255816230796473E-01, &
       0.5475589657435200E-01,     0.7503967481091995E-01, &
       0.9312545458369761E-01,     0.1093871588022976E+00, &
       0.1234919762620659E+00,     0.1347092173114733E+00, &
       0.1427759385770601E+00,     0.1477391049013385E+00, &
       0.1494455540029169E+00/
!
  DATA WG(1),WG(2),WG(3),WG(4),WG(5)/ &
       0.6667134430868814E-01,     0.1494513491505806E+00, &
       0.2190863625159820E+00,     0.2692667193099964E+00, &
       0.2955242247147529E+00/
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
!***FIRST EXECUTABLE STATEMENT  QK21
  EPMACH = R1MACH(4)
  UFLOW = R1MACH(1)
!
  CENTR = 0.5E+00*(A+B)
  HLGTH = 0.5E+00*(B-A)
  DHLGTH = ABS(HLGTH)
!
!           COMPUTE THE 21-POINT KRONROD APPROXIMATION TO
!           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
  RESG = 0.0E+00
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
  RESKH = RESK*0.5E+00
  RESASC = WGK(11)*ABS(FC-RESKH)
  DO 20 J=1,10
    RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
  RESULT = RESK*HLGTH
  RESABS = RESABS*DHLGTH
  RESASC = RESASC*DHLGTH
  ABSERR = ABS((RESK-RESG)*HLGTH)
  if ( RESASC /= 0.0E+00.AND.ABSERR /= 0.0E+00) &
    ABSERR = RESASC*MIN(0.1E+01, &
    (0.2E+03*ABSERR/RESASC)**1.5E+00)
  if ( RESABS > UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX &
    ((EPMACH*0.5E+02)*RESABS,ABSERR)
  return
end
