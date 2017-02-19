subroutine QK31 (F, A, B, RESULT, ABSERR, RESABS, RESASC)
!
!! QK31 estimates an integral with a 31 point Gauss-Kronrod rule.
!
!***PURPOSE  To compute I = Integral of F over (A,B) with error
!                           estimate
!                       J = Integral of ABS(F) over (A,B)
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A2
!***TYPE      SINGLE PRECISION (QK31-S, DQK31-D)
!***KEYWORDS  31-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
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
!                       Declared E X T E R N A L in the calling program.
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
!                       RESULT is computed by applying the 31-POINT
!                       GAUSS-KRONROD RULE (RESK), obtained by optimal
!                       addition of abscissae to the 15-POINT GAUSS
!                       RULE (RESG).
!
!              ABSERR - Real
!                       Estimate of the modulus of the modulus,
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
!***END PROLOGUE  QK31
  REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,F,FC,FSUM,FVAL1,FVAL2, &
    FV1,FV2,HLGTH,RESABS,RESASC,RESG,RESK,RESKH,RESULT,R1MACH,UFLOW, &
    WG,WGK,XGK
  INTEGER J,JTW,JTWM1
  EXTERNAL F
!
  DIMENSION FV1(15),FV2(15),XGK(16),WGK(16),WG(8)
!
!           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
!           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!           CORRESPONDING WEIGHTS ARE GIVEN.
!
!           XGK    - ABSCISSAE OF THE 31-POINT KRONROD RULE
!                    XGK(2), XGK(4), ...  ABSCISSAE OF THE 15-POINT
!                    GAUSS RULE
!                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!                    ADDED TO THE 15-POINT GAUSS RULE
!
!           WGK    - WEIGHTS OF THE 31-POINT KRONROD RULE
!
!           WG     - WEIGHTS OF THE 15-POINT GAUSS RULE
!
  SAVE XGK, WGK, WG
  DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8), &
    XGK(9),XGK(10),XGK(11),XGK(12),XGK(13),XGK(14),XGK(15), &
    XGK(16)/ &
       0.9980022986933971E+00,   0.9879925180204854E+00, &
       0.9677390756791391E+00,   0.9372733924007059E+00, &
       0.8972645323440819E+00,   0.8482065834104272E+00, &
       0.7904185014424659E+00,   0.7244177313601700E+00, &
       0.6509967412974170E+00,   0.5709721726085388E+00, &
       0.4850818636402397E+00,   0.3941513470775634E+00, &
       0.2991800071531688E+00,   0.2011940939974345E+00, &
       0.1011420669187175E+00,   0.0E+00               /
  DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8), &
    WGK(9),WGK(10),WGK(11),WGK(12),WGK(13),WGK(14),WGK(15), &
    WGK(16)/ &
       0.5377479872923349E-02,   0.1500794732931612E-01, &
       0.2546084732671532E-01,   0.3534636079137585E-01, &
       0.4458975132476488E-01,   0.5348152469092809E-01, &
       0.6200956780067064E-01,   0.6985412131872826E-01, &
       0.7684968075772038E-01,   0.8308050282313302E-01, &
       0.8856444305621177E-01,   0.9312659817082532E-01, &
       0.9664272698362368E-01,   0.9917359872179196E-01, &
       0.1007698455238756E+00,   0.1013300070147915E+00/
  DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/ &
       0.3075324199611727E-01,   0.7036604748810812E-01, &
       0.1071592204671719E+00,   0.1395706779261543E+00, &
       0.1662692058169939E+00,   0.1861610000155622E+00, &
       0.1984314853271116E+00,   0.2025782419255613E+00/
!
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!           CENTR  - MID POINT OF THE INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           ABSC   - ABSCISSA
!           FVAL*  - FUNCTION VALUE
!           RESG   - RESULT OF THE 15-POINT GAUSS FORMULA
!           RESK   - RESULT OF THE 31-POINT KRONROD FORMULA
!           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F OVER (A,B),
!                    I.E. TO I/(B-A)
!
!           MACHINE DEPENDENT CONSTANTS
!           ---------------------------
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  QK31
  EPMACH = R1MACH(4)
  UFLOW = R1MACH(1)
!
  CENTR = 0.5E+00*(A+B)
  HLGTH = 0.5E+00*(B-A)
  DHLGTH = ABS(HLGTH)
!
!           COMPUTE THE 31-POINT KRONROD APPROXIMATION TO
!           THE INTEGRAL, AND ESTIMATE THE ABSOLUTE ERROR.
!
  FC = F(CENTR)
  RESG = WG(8)*FC
  RESK = WGK(16)*FC
  RESABS = ABS(RESK)
  DO 10 J=1,7
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
  DO 15 J = 1,8
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
  RESKH = RESK*0.5E+00
  RESASC = WGK(16)*ABS(FC-RESKH)
  DO 20 J=1,15
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
