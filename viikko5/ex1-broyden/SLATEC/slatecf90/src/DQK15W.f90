subroutine DQK15W (F, W, P1, P2, P3, P4, KP, A, B, RESULT, ABSERR, &
     RESABS, RESASC)
!
!! DQK15W computes Integral of F*W over (A,B), with error estimate.
!                       J = Integral of ABS(F*W) over (A,B)
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A2A2
!***TYPE      DOUBLE PRECISION (QK15W-S, DQK15W-D)
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
!             ON ENTRY
!              F      - Double precision
!                       Function subprogram defining the integrand
!                       function F(X). The actual name for F needs to be
!                       declared E X T E R N A L in the driver program.
!
!              W      - Double precision
!                       Function subprogram defining the integrand
!                       WEIGHT function W(X). The actual name for W
!                       needs to be declared E X T E R N A L in the
!                       calling program.
!
!              P1, P2, P3, P4 - Double precision
!                       Parameters in the WEIGHT function
!
!              KP     - Integer
!                       Key for indicating the type of WEIGHT function
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
!                       RESULT is computed by applying the 15-point
!                       Kronrod rule (RESK) obtained by optimal addition
!                       of abscissae to the 7-point Gauss rule (RESG).
!
!              ABSERR - Double precision
!                       Estimate of the modulus of the absolute error,
!                       which should equal or exceed ABS(I-RESULT)
!
!              RESABS - Double precision
!                       Approximation to the integral of ABS(F)
!
!              RESASC - Double precision
!                       Approximation to the integral of ABS(F-I/(B-A))
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   810101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DQK15W
!
  DOUBLE PRECISION A,ABSC,ABSC1,ABSC2,ABSERR,B,CENTR,DHLGTH, &
    D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH, &
    P1,P2,P3,P4,RESABS,RESASC,RESG,RESK,RESKH,RESULT,UFLOW,W,WG,WGK, &
    XGK
  INTEGER J,JTW,JTWM1,KP
  EXTERNAL F, W
!
  DIMENSION FV1(7),FV2(7),XGK(8),WGK(8),WG(4)
!
!           THE ABSCISSAE AND WEIGHTS ARE GIVEN FOR THE INTERVAL (-1,1).
!           BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND THEIR
!           CORRESPONDING WEIGHTS ARE GIVEN.
!
!           XGK    - ABSCISSAE OF THE 15-POINT GAUSS-KRONROD RULE
!                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
!                    GAUSS RULE
!                    XGK(1), XGK(3), ... ABSCISSAE WHICH ARE OPTIMALLY
!                    ADDED TO THE 7-POINT GAUSS RULE
!
!           WGK    - WEIGHTS OF THE 15-POINT GAUSS-KRONROD RULE
!
!           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE
!
  SAVE XGK, WGK, WG
  DATA XGK(1),XGK(2),XGK(3),XGK(4),XGK(5),XGK(6),XGK(7),XGK(8)/ &
       0.9914553711208126D+00,     0.9491079123427585D+00, &
       0.8648644233597691D+00,     0.7415311855993944D+00, &
       0.5860872354676911D+00,     0.4058451513773972D+00, &
       0.2077849550078985D+00,     0.0000000000000000D+00/
!
  DATA WGK(1),WGK(2),WGK(3),WGK(4),WGK(5),WGK(6),WGK(7),WGK(8)/ &
       0.2293532201052922D-01,     0.6309209262997855D-01, &
       0.1047900103222502D+00,     0.1406532597155259D+00, &
       0.1690047266392679D+00,     0.1903505780647854D+00, &
       0.2044329400752989D+00,     0.2094821410847278D+00/
!
  DATA WG(1),WG(2),WG(3),WG(4)/ &
       0.1294849661688697D+00,    0.2797053914892767D+00, &
       0.3818300505051889D+00,    0.4179591836734694D+00/
!
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           CENTR  - MID POINT OF THE INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           ABSC*  - ABSCISSA
!           FVAL*  - FUNCTION VALUE
!           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
!           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
!           RESKH  - APPROXIMATION TO THE MEAN VALUE OF F*W OVER (A,B),
!                    I.E. TO I/(B-A)
!
!           MACHINE DEPENDENT CONSTANTS
!           ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  DQK15W
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
!
  CENTR = 0.5D+00*(A+B)
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
!
!           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO THE
!           INTEGRAL, AND ESTIMATE THE ERROR.
!
  FC = F(CENTR)*W(CENTR,P1,P2,P3,P4,KP)
  RESG = WG(4)*FC
  RESK = WGK(8)*FC
  RESABS = ABS(RESK)
  DO 10 J=1,3
    JTW = J*2
    ABSC = HLGTH*XGK(JTW)
    ABSC1 = CENTR-ABSC
    ABSC2 = CENTR+ABSC
    FVAL1 = F(ABSC1)*W(ABSC1,P1,P2,P3,P4,KP)
    FVAL2 = F(ABSC2)*W(ABSC2,P1,P2,P3,P4,KP)
    FV1(JTW) = FVAL1
    FV2(JTW) = FVAL2
    FSUM = FVAL1+FVAL2
    RESG = RESG+WG(J)*FSUM
    RESK = RESK+WGK(JTW)*FSUM
    RESABS = RESABS+WGK(JTW)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
  DO 15 J=1,4
    JTWM1 = J*2-1
    ABSC = HLGTH*XGK(JTWM1)
    ABSC1 = CENTR-ABSC
    ABSC2 = CENTR+ABSC
    FVAL1 = F(ABSC1)*W(ABSC1,P1,P2,P3,P4,KP)
    FVAL2 = F(ABSC2)*W(ABSC2,P1,P2,P3,P4,KP)
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
  if ( RESABS > UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX((EPMACH* &
    0.5D+02)*RESABS,ABSERR)
  return
end
