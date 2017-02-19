subroutine QNG (F, A, B, EPSABS, EPSREL, RESULT, ABSERR, NEVAL, &
     IER)
!
!! QNG calculates an approximation RESULT to an integral
!  I = integral of F over (A,B),
!            hopefully satisfying following claim for accuracy
!            ABS(I-RESULT) <= MAX(EPSABS,EPSREL*ABS(I)).
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A1
!***TYPE      SINGLE PRECISION (QNG-S, DQNG-D)
!***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD(PATTERSON) RULES,
!             NONADAPTIVE, QUADPACK, QUADRATURE, SMOOTH INTEGRAND
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
! NON-ADAPTIVE INTEGRATION
! STANDARD FORTRAN SUBROUTINE
! REAL VERSION
!
!           F      - Real version
!                    Function subprogram defining the integrand function
!                    F(X). The actual name for F needs to be declared
!                    E X T E R N A L in the driver program.
!
!           A      - Real version
!                    Lower limit of integration
!
!           B      - Real version
!                    Upper limit of integration
!
!           EPSABS - Real
!                    Absolute accuracy requested
!           EPSREL - Real
!                    Relative accuracy requested
!                    If  EPSABS <= 0
!                    And EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28),
!                    The routine will end with IER = 6.
!
!         ON RETURN
!           RESULT - Real
!                    Approximation to the integral I
!                    Result is obtained by applying the 21-POINT
!                    GAUSS-KRONROD RULE (RES21) obtained by optimal
!                    addition of abscissae to the 10-POINT GAUSS RULE
!                    (RES10), or by applying the 43-POINT RULE (RES43)
!                    obtained by optimal addition of abscissae to the
!                    21-POINT GAUSS-KRONROD RULE, or by applying the
!                    87-POINT RULE (RES87) obtained by optimal addition
!                    of abscissae to the 43-POINT RULE.
!
!           ABSERR - Real
!                    Estimate of the modulus of the absolute error,
!                    which should EQUAL or EXCEED ABS(I-RESULT)
!
!           NEVAL  - Integer
!                    Number of integrand evaluations
!
!           IER    - IER = 0 normal and reliable termination of the
!                            routine. It is assumed that the requested
!                            accuracy has been achieved.
!                    IER > 0 Abnormal termination of the routine. It is
!                            assumed that the requested accuracy has
!                            not been achieved.
!           ERROR MESSAGES
!                    IER = 1 The maximum number of steps has been
!                            executed. The integral is probably too
!                            difficult to be calculated by DQNG.
!                        = 6 The input is invalid, because
!                            EPSABS <= 0 AND
!                            EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28).
!                            RESULT, ABSERR and NEVAL are set to zero.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  QNG
!
  REAL A,ABSC,ABSERR,B,CENTR,DHLGTH,EPMACH,EPSABS,EPSREL,F,FCENTR, &
    FVAL,FVAL1,FVAL2,FV1,FV2,FV3,FV4,HLGTH,RESULT,RES10,RES21,RES43, &
    RES87,RESABS,RESASC,RESKH,R1MACH,SAVFUN,UFLOW,W10,W21A,W43A, &
    W43B,W87A,W87B,X1,X2,X3,X4
  INTEGER IER,IPX,K,L,NEVAL
  EXTERNAL F
!
  DIMENSION FV1(5),FV2(5),FV3(5),FV4(5),X1(5),X2(5),X3(11),X4(22), &
    W10(5),W21A(5),W21B(6),W43A(10),W43B(12),W87A(21),W87B(23), &
    SAVFUN(21)
!
!           THE FOLLOWING DATA STATEMENTS CONTAIN THE
!           ABSCISSAE AND WEIGHTS OF THE INTEGRATION RULES USED.
!
!           X1      ABSCISSAE COMMON TO THE 10-, 21-, 43-
!                   AND 87-POINT RULE
!           X2      ABSCISSAE COMMON TO THE 21-, 43- AND
!                   87-POINT RULE
!           X3      ABSCISSAE COMMON TO THE 43- AND 87-POINT
!                   RULE
!           X4      ABSCISSAE OF THE 87-POINT RULE
!           W10     WEIGHTS OF THE 10-POINT FORMULA
!           W21A    WEIGHTS OF THE 21-POINT FORMULA FOR
!                   ABSCISSAE X1
!           W21B    WEIGHTS OF THE 21-POINT FORMULA FOR
!                   ABSCISSAE X2
!           W43A    WEIGHTS OF THE 43-POINT FORMULA FOR
!                   ABSCISSAE X1, X3
!           W43B    WEIGHTS OF THE 43-POINT FORMULA FOR
!                   ABSCISSAE X3
!           W87A    WEIGHTS OF THE 87-POINT FORMULA FOR
!                   ABSCISSAE X1, X2, X3
!           W87B    WEIGHTS OF THE 87-POINT FORMULA FOR
!                   ABSCISSAE X4
!
  SAVE X1, X2, X3, X4, W10, W21A, W21B, W43A, W43B, W87A, W87B
  DATA X1(1),X1(2),X1(3),X1(4),X1(5)/ &
       0.9739065285171717E+00,     0.8650633666889845E+00, &
       0.6794095682990244E+00,     0.4333953941292472E+00, &
       0.1488743389816312E+00/
  DATA X2(1),X2(2),X2(3),X2(4),X2(5)/ &
       0.9956571630258081E+00,     0.9301574913557082E+00, &
       0.7808177265864169E+00,     0.5627571346686047E+00, &
       0.2943928627014602E+00/
  DATA X3(1),X3(2),X3(3),X3(4),X3(5),X3(6),X3(7),X3(8), &
    X3(9),X3(10),X3(11)/ &
       0.9993333609019321E+00,     0.9874334029080889E+00, &
       0.9548079348142663E+00,     0.9001486957483283E+00, &
       0.8251983149831142E+00,     0.7321483889893050E+00, &
       0.6228479705377252E+00,     0.4994795740710565E+00, &
       0.3649016613465808E+00,     0.2222549197766013E+00, &
       0.7465061746138332E-01/
  DATA X4(1),X4(2),X4(3),X4(4),X4(5),X4(6),X4(7),X4(8),X4(9), &
    X4(10),X4(11),X4(12),X4(13),X4(14),X4(15),X4(16),X4(17),X4(18), &
    X4(19),X4(20),X4(21),X4(22)/   0.9999029772627292E+00, &
       0.9979898959866787E+00,     0.9921754978606872E+00, &
       0.9813581635727128E+00,     0.9650576238583846E+00, &
       0.9431676131336706E+00,     0.9158064146855072E+00, &
       0.8832216577713165E+00,     0.8457107484624157E+00, &
       0.8035576580352310E+00,     0.7570057306854956E+00, &
       0.7062732097873218E+00,     0.6515894665011779E+00, &
       0.5932233740579611E+00,     0.5314936059708319E+00, &
       0.4667636230420228E+00,     0.3994248478592188E+00, &
       0.3298748771061883E+00,     0.2585035592021616E+00, &
       0.1856953965683467E+00,     0.1118422131799075E+00, &
       0.3735212339461987E-01/
  DATA W10(1),W10(2),W10(3),W10(4),W10(5)/ &
       0.6667134430868814E-01,     0.1494513491505806E+00, &
       0.2190863625159820E+00,     0.2692667193099964E+00, &
       0.2955242247147529E+00/
  DATA W21A(1),W21A(2),W21A(3),W21A(4),W21A(5)/ &
       0.3255816230796473E-01,     0.7503967481091995E-01, &
       0.1093871588022976E+00,     0.1347092173114733E+00, &
       0.1477391049013385E+00/
  DATA W21B(1),W21B(2),W21B(3),W21B(4),W21B(5),W21B(6)/ &
       0.1169463886737187E-01,     0.5475589657435200E-01, &
       0.9312545458369761E-01,     0.1234919762620659E+00, &
       0.1427759385770601E+00,     0.1494455540029169E+00/
  DATA W43A(1),W43A(2),W43A(3),W43A(4),W43A(5),W43A(6),W43A(7), &
    W43A(8),W43A(9),W43A(10)/      0.1629673428966656E-01, &
       0.3752287612086950E-01,     0.5469490205825544E-01, &
       0.6735541460947809E-01,     0.7387019963239395E-01, &
       0.5768556059769796E-02,     0.2737189059324884E-01, &
       0.4656082691042883E-01,     0.6174499520144256E-01, &
       0.7138726726869340E-01/
  DATA W43B(1),W43B(2),W43B(3),W43B(4),W43B(5),W43B(6), &
    W43B(7),W43B(8),W43B(9),W43B(10),W43B(11),W43B(12)/ &
       0.1844477640212414E-02,     0.1079868958589165E-01, &
       0.2189536386779543E-01,     0.3259746397534569E-01, &
       0.4216313793519181E-01,     0.5074193960018458E-01, &
       0.5837939554261925E-01,     0.6474640495144589E-01, &
       0.6956619791235648E-01,     0.7282444147183321E-01, &
       0.7450775101417512E-01,     0.7472214751740301E-01/
  DATA W87A(1),W87A(2),W87A(3),W87A(4),W87A(5),W87A(6), &
    W87A(7),W87A(8),W87A(9),W87A(10),W87A(11),W87A(12), &
    W87A(13),W87A(14),W87A(15),W87A(16),W87A(17),W87A(18), &
    W87A(19),W87A(20),W87A(21)/ &
       0.8148377384149173E-02,     0.1876143820156282E-01, &
       0.2734745105005229E-01,     0.3367770731163793E-01, &
       0.3693509982042791E-01,     0.2884872430211531E-02, &
       0.1368594602271270E-01,     0.2328041350288831E-01, &
       0.3087249761171336E-01,     0.3569363363941877E-01, &
       0.9152833452022414E-03,     0.5399280219300471E-02, &
       0.1094767960111893E-01,     0.1629873169678734E-01, &
       0.2108156888920384E-01,     0.2537096976925383E-01, &
       0.2918969775647575E-01,     0.3237320246720279E-01, &
       0.3478309895036514E-01,     0.3641222073135179E-01, &
       0.3725387550304771E-01/
  DATA W87B(1),W87B(2),W87B(3),W87B(4),W87B(5),W87B(6),W87B(7), &
    W87B(8),W87B(9),W87B(10),W87B(11),W87B(12),W87B(13),W87B(14), &
    W87B(15),W87B(16),W87B(17),W87B(18),W87B(19),W87B(20), &
    W87B(21),W87B(22),W87B(23)/    0.2741455637620724E-03, &
       0.1807124155057943E-02,     0.4096869282759165E-02, &
       0.6758290051847379E-02,     0.9549957672201647E-02, &
       0.1232944765224485E-01,     0.1501044734638895E-01, &
       0.1754896798624319E-01,     0.1993803778644089E-01, &
       0.2219493596101229E-01,     0.2433914712600081E-01, &
       0.2637450541483921E-01,     0.2828691078877120E-01, &
       0.3005258112809270E-01,     0.3164675137143993E-01, &
       0.3305041341997850E-01,     0.3425509970422606E-01, &
       0.3526241266015668E-01,     0.3607698962288870E-01, &
       0.3669860449845609E-01,     0.3712054926983258E-01, &
       0.3733422875193504E-01,     0.3736107376267902E-01/
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           CENTR  - MID POINT OF THE INTEGRATION INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL
!           FCENTR - FUNCTION VALUE AT MID POINT
!           ABSC   - ABSCISSA
!           FVAL   - FUNCTION VALUE
!           SAVFUN - ARRAY OF FUNCTION VALUES WHICH
!                    HAVE ALREADY BEEN COMPUTED
!           RES10  - 10-POINT GAUSS RESULT
!           RES21  - 21-POINT KRONROD RESULT
!           RES43  - 43-POINT RESULT
!           RES87  - 87-POINT RESULT
!           RESABS - APPROXIMATION TO THE INTEGRAL OF ABS(F)
!           RESASC - APPROXIMATION TO THE INTEGRAL OF ABS(F-I/(B-A))
!
!           MACHINE DEPENDENT CONSTANTS
!           ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  QNG
  EPMACH = R1MACH(4)
  UFLOW = R1MACH(1)
!
!           TEST ON VALIDITY OF PARAMETERS
!           ------------------------------
!
  RESULT = 0.0E+00
  ABSERR = 0.0E+00
  NEVAL = 0
  IER = 6
  if ( EPSABS <= 0.0E+00.AND.EPSREL < MAX(0.5E-14,0.5E+02*EPMACH)) &
    go to 80
  HLGTH = 0.5E+00*(B-A)
  DHLGTH = ABS(HLGTH)
  CENTR = 0.5E+00*(B+A)
  FCENTR = F(CENTR)
  NEVAL = 21
  IER = 1
!
!          COMPUTE THE INTEGRAL USING THE 10- AND 21-POINT FORMULA.
!
  DO 70 L = 1,3
  go to (5,25,45),L
    5 RES10 = 0.0E+00
  RES21 = W21B(6)*FCENTR
  RESABS = W21B(6)*ABS(FCENTR)
  DO 10 K=1,5
    ABSC = HLGTH*X1(K)
    FVAL1 = F(CENTR+ABSC)
    FVAL2 = F(CENTR-ABSC)
    FVAL = FVAL1+FVAL2
    RES10 = RES10+W10(K)*FVAL
    RES21 = RES21+W21A(K)*FVAL
    RESABS = RESABS+W21A(K)*(ABS(FVAL1)+ABS(FVAL2))
    SAVFUN(K) = FVAL
    FV1(K) = FVAL1
    FV2(K) = FVAL2
   10 CONTINUE
  IPX = 5
  DO 15 K=1,5
    IPX = IPX+1
    ABSC = HLGTH*X2(K)
    FVAL1 = F(CENTR+ABSC)
    FVAL2 = F(CENTR-ABSC)
    FVAL = FVAL1+FVAL2
    RES21 = RES21+W21B(K)*FVAL
    RESABS = RESABS+W21B(K)*(ABS(FVAL1)+ABS(FVAL2))
    SAVFUN(IPX) = FVAL
    FV3(K) = FVAL1
    FV4(K) = FVAL2
   15 CONTINUE
!
!          TEST FOR CONVERGENCE.
!
  RESULT = RES21*HLGTH
  RESABS = RESABS*DHLGTH
  RESKH = 0.5E+00*RES21
  RESASC = W21B(6)*ABS(FCENTR-RESKH)
  DO 20 K = 1,5
    RESASC = RESASC+W21A(K)*(ABS(FV1(K)-RESKH)+ABS(FV2(K)-RESKH)) &
                    +W21B(K)*(ABS(FV3(K)-RESKH)+ABS(FV4(K)-RESKH))
   20 CONTINUE
  ABSERR = ABS((RES21-RES10)*HLGTH)
  RESASC = RESASC*DHLGTH
  go to 65
!
!          COMPUTE THE INTEGRAL USING THE 43-POINT FORMULA.
!
   25 RES43 = W43B(12)*FCENTR
  NEVAL = 43
  DO 30 K=1,10
    RES43 = RES43+SAVFUN(K)*W43A(K)
   30 CONTINUE
  DO 40 K=1,11
    IPX = IPX+1
    ABSC = HLGTH*X3(K)
    FVAL = F(ABSC+CENTR)+F(CENTR-ABSC)
    RES43 = RES43+FVAL*W43B(K)
    SAVFUN(IPX) = FVAL
   40 CONTINUE
!
!          TEST FOR CONVERGENCE.
!
  RESULT = RES43*HLGTH
  ABSERR = ABS((RES43-RES21)*HLGTH)
  go to 65
!
!          COMPUTE THE INTEGRAL USING THE 87-POINT FORMULA.
!
   45 RES87 = W87B(23)*FCENTR
  NEVAL = 87
  DO 50 K=1,21
    RES87 = RES87+SAVFUN(K)*W87A(K)
   50 CONTINUE
  DO 60 K=1,22
    ABSC = HLGTH*X4(K)
    RES87 = RES87+W87B(K)*(F(ABSC+CENTR)+F(CENTR-ABSC))
   60 CONTINUE
  RESULT = RES87*HLGTH
  ABSERR = ABS((RES87-RES43)*HLGTH)
   65 if ( RESASC /= 0.0E+00.AND.ABSERR /= 0.0E+00) &
    ABSERR = RESASC*MIN(0.1E+01, &
    (0.2E+03*ABSERR/RESASC)**1.5E+00)
  if (RESABS > UFLOW/(0.5E+02*EPMACH)) ABSERR = MAX &
    ((EPMACH*0.5E+02)*RESABS,ABSERR)
  if (ABSERR <= MAX(EPSABS,EPSREL*ABS(RESULT))) IER = 0
! ***JUMP OUT OF DO-LOOP
  if (IER == 0) go to 999
   70 CONTINUE
   80 call XERMSG ('SLATEC', 'QNG', 'ABNORMAL RETURN', IER, 0)
  999 RETURN
end
