subroutine DQK15I (F, BOUN, INF, A, B, RESULT, ABSERR, RESABS, &
     RESASC)
!
!! DQK15I computes an integral over an infinite range.
!
!***PURPOSE  The original (infinite integration range is mapped
!            onto the interval (0,1) and (A,B) is a part of (0,1).
!            it is the purpose to compute
!            I = Integral of transformed integrand over (A,B),
!            J = Integral of ABS(Transformed Integrand) over (A,B).
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A3A2, H2A4A2
!***TYPE      DOUBLE PRECISION (QK15I-S, DQK15I-D)
!***KEYWORDS  15-POINT GAUSS-KRONROD RULES, QUADPACK, QUADRATURE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!           Integration Rule
!           Standard Fortran subroutine
!           Double precision version
!
!           PARAMETERS
!            ON ENTRY
!              F      - Double precision
!                       Function subprogram defining the integrand
!                       FUNCTION F(X). The actual name for F needs to be
!                       Declared E X T E R N A L in the calling program.
!
!              BOUN   - Double precision
!                       Finite bound of original integration
!                       Range (SET TO ZERO if INF = +2)
!
!              INF    - Integer
!                       If INF = -1, the original interval is
!                                   (-INFINITY,BOUND),
!                       If INF = +1, the original interval is
!                                   (BOUND,+INFINITY),
!                       If INF = +2, the original interval is
!                                   (-INFINITY,+INFINITY) AND
!                       The integral is computed as the sum of two
!                       integrals, one over (-INFINITY,0) and one over
!                       (0,+INFINITY).
!
!              A      - Double precision
!                       Lower limit for integration over subrange
!                       of (0,1)
!
!              B      - Double precision
!                       Upper limit for integration over subrange
!                       of (0,1)
!
!            ON RETURN
!              RESULT - Double precision
!                       Approximation to the integral I
!                       Result is computed by applying the 15-POINT
!                       KRONROD RULE(RESK) obtained by optimal addition
!                       of abscissae to the 7-POINT GAUSS RULE(RESG).
!
!              ABSERR - Double precision
!                       Estimate of the modulus of the absolute error,
!                       WHICH SHOULD EQUAL or EXCEED ABS(I-RESULT)
!
!              RESABS - Double precision
!                       Approximation to the integral J
!
!              RESASC - Double precision
!                       Approximation to the integral of
!                       ABS((TRANSFORMED INTEGRAND)-I/(B-A)) over (A,B)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DQK15I
!
  DOUBLE PRECISION A,ABSC,ABSC1,ABSC2,ABSERR,B,BOUN,CENTR,DINF, &
    D1MACH,EPMACH,F,FC,FSUM,FVAL1,FVAL2,FV1,FV2,HLGTH, &
    RESABS,RESASC,RESG,RESK,RESKH,RESULT,TABSC1,TABSC2,UFLOW,WG,WGK, &
    XGK
  INTEGER INF,J
  EXTERNAL F
!
  DIMENSION FV1(7),FV2(7),XGK(8),WGK(8),WG(8)
!
!           THE ABSCISSAE AND WEIGHTS ARE SUPPLIED FOR THE INTERVAL
!           (-1,1).  BECAUSE OF SYMMETRY ONLY THE POSITIVE ABSCISSAE AND
!           THEIR CORRESPONDING WEIGHTS ARE GIVEN.
!
!           XGK    - ABSCISSAE OF THE 15-POINT KRONROD RULE
!                    XGK(2), XGK(4), ... ABSCISSAE OF THE 7-POINT
!                    GAUSS RULE
!                    XGK(1), XGK(3), ...  ABSCISSAE WHICH ARE OPTIMALLY
!                    ADDED TO THE 7-POINT GAUSS RULE
!
!           WGK    - WEIGHTS OF THE 15-POINT KRONROD RULE
!
!           WG     - WEIGHTS OF THE 7-POINT GAUSS RULE, CORRESPONDING
!                    TO THE ABSCISSAE XGK(2), XGK(4), ...
!                    WG(1), WG(3), ... ARE SET TO ZERO.
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
  DATA WG(1),WG(2),WG(3),WG(4),WG(5),WG(6),WG(7),WG(8)/ &
       0.0000000000000000D+00,     0.1294849661688697D+00, &
       0.0000000000000000D+00,     0.2797053914892767D+00, &
       0.0000000000000000D+00,     0.3818300505051189D+00, &
       0.0000000000000000D+00,     0.4179591836734694D+00/
!
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           CENTR  - MID POINT OF THE INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTERVAL
!           ABSC*  - ABSCISSA
!           TABSC* - TRANSFORMED ABSCISSA
!           FVAL*  - FUNCTION VALUE
!           RESG   - RESULT OF THE 7-POINT GAUSS FORMULA
!           RESK   - RESULT OF THE 15-POINT KRONROD FORMULA
!           RESKH  - APPROXIMATION TO THE MEAN VALUE OF THE TRANSFORMED
!                    INTEGRAND OVER (A,B), I.E. TO I/(B-A)
!
!           MACHINE DEPENDENT CONSTANTS
!           ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  DQK15I
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
  DINF = MIN(1,INF)
!
  CENTR = 0.5D+00*(A+B)
  HLGTH = 0.5D+00*(B-A)
  TABSC1 = BOUN+DINF*(0.1D+01-CENTR)/CENTR
  FVAL1 = F(TABSC1)
  if ( INF == 2) FVAL1 = FVAL1+F(-TABSC1)
  FC = (FVAL1/CENTR)/CENTR
!
!           COMPUTE THE 15-POINT KRONROD APPROXIMATION TO
!           THE INTEGRAL, AND ESTIMATE THE ERROR.
!
  RESG = WG(8)*FC
  RESK = WGK(8)*FC
  RESABS = ABS(RESK)
  DO 10 J=1,7
    ABSC = HLGTH*XGK(J)
    ABSC1 = CENTR-ABSC
    ABSC2 = CENTR+ABSC
    TABSC1 = BOUN+DINF*(0.1D+01-ABSC1)/ABSC1
    TABSC2 = BOUN+DINF*(0.1D+01-ABSC2)/ABSC2
    FVAL1 = F(TABSC1)
    FVAL2 = F(TABSC2)
    if ( INF == 2) FVAL1 = FVAL1+F(-TABSC1)
    if ( INF == 2) FVAL2 = FVAL2+F(-TABSC2)
    FVAL1 = (FVAL1/ABSC1)/ABSC1
    FVAL2 = (FVAL2/ABSC2)/ABSC2
    FV1(J) = FVAL1
    FV2(J) = FVAL2
    FSUM = FVAL1+FVAL2
    RESG = RESG+WG(J)*FSUM
    RESK = RESK+WGK(J)*FSUM
    RESABS = RESABS+WGK(J)*(ABS(FVAL1)+ABS(FVAL2))
   10 CONTINUE
  RESKH = RESK*0.5D+00
  RESASC = WGK(8)*ABS(FC-RESKH)
  DO 20 J=1,7
    RESASC = RESASC+WGK(J)*(ABS(FV1(J)-RESKH)+ABS(FV2(J)-RESKH))
   20 CONTINUE
  RESULT = RESK*HLGTH
  RESASC = RESASC*HLGTH
  RESABS = RESABS*HLGTH
  ABSERR = ABS((RESK-RESG)*HLGTH)
  if ( RESASC /= 0.0D+00.AND.ABSERR /= 0.D0) ABSERR = RESASC* &
   MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  if ( RESABS > UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
   ((EPMACH*0.5D+02)*RESABS,ABSERR)
  return
end
