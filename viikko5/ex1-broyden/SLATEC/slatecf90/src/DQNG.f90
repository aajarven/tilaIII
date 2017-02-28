subroutine DQNG (F, A, B, EPSABS, EPSREL, RESULT, ABSERR, NEVAL, IER)
!
!! DQNG approximates the integral of a function over a finite interval.
!
!***PURPOSE  The routine calculates an approximation result to a
!            given definite integral I = integral of F over (A,B),
!            hopefully satisfying following claim for accuracy
!            ABS(I-RESULT) <= MAX(EPSABS,EPSREL*ABS(I)).
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A1
!***TYPE      DOUBLE PRECISION (QNG-S, DQNG-D)
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
! DOUBLE PRECISION VERSION
!
!           F      - Double precision
!                    Function subprogram defining the integrand function
!                    F(X). The actual name for F needs to be declared
!                    E X T E R N A L in the driver program.
!
!           A      - Double precision
!                    Lower limit of integration
!
!           B      - Double precision
!                    Upper limit of integration
!
!           EPSABS - Double precision
!                    Absolute accuracy requested
!           EPSREL - Double precision
!                    Relative accuracy requested
!                    If  EPSABS <= 0
!                    And EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28),
!                    The routine will end with IER = 6.
!
!         ON RETURN
!           RESULT - Double precision
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
!           ABSERR - Double precision
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
!***ROUTINES CALLED  D1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  DQNG
!
  DOUBLE PRECISION A,ABSC,ABSERR,B,CENTR,DHLGTH, &
    D1MACH,EPMACH,EPSABS,EPSREL,F,FCENTR,FVAL,FVAL1,FVAL2,FV1,FV2, &
    FV3,FV4,HLGTH,RESULT,RES10,RES21,RES43,RES87,RESABS,RESASC, &
    RESKH,SAVFUN,UFLOW,W10,W21A,W21B,W43A,W43B,W87A,W87B,X1,X2,X3,X4
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
!           X1      ABSCISSAE COMMON TO THE 10-, 21-, 43- AND 87-
!                   POINT RULE
!           X2      ABSCISSAE COMMON TO THE 21-, 43- AND 87-POINT RULE
!           X3      ABSCISSAE COMMON TO THE 43- AND 87-POINT RULE
!           X4      ABSCISSAE OF THE 87-POINT RULE
!           W10     WEIGHTS OF THE 10-POINT FORMULA
!           W21A    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X1
!           W21B    WEIGHTS OF THE 21-POINT FORMULA FOR ABSCISSAE X2
!           W43A    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X1, X3
!           W43B    WEIGHTS OF THE 43-POINT FORMULA FOR ABSCISSAE X3
!           W87A    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X1,
!                   X2, X3
!           W87B    WEIGHTS OF THE 87-POINT FORMULA FOR ABSCISSAE X4
!
!
! GAUSS-KRONROD-PATTERSON QUADRATURE COEFFICIENTS FOR USE IN
! QUADPACK ROUTINE QNG.  THESE COEFFICIENTS WERE CALCULATED WITH
! 101 DECIMAL DIGIT ARITHMETIC BY L. W. FULLERTON, BELL LABS, NOV 1981.
!
  SAVE X1, W10, X2, W21A, W21B, X3, W43A, W43B, X4, W87A, W87B
  DATA X1    (  1) / 0.973906528517171720077964012084452D0 /
  DATA X1    (  2) / 0.865063366688984510732096688423493D0 /
  DATA X1    (  3) / 0.679409568299024406234327365114874D0 /
  DATA X1    (  4) / 0.433395394129247190799265943165784D0 /
  DATA X1    (  5) / 0.148874338981631210884826001129720D0 /
  DATA W10   (  1) / 0.066671344308688137593568809893332D0 /
  DATA W10   (  2) / 0.149451349150580593145776339657697D0 /
  DATA W10   (  3) / 0.219086362515982043995534934228163D0 /
  DATA W10   (  4) / 0.269266719309996355091226921569469D0 /
  DATA W10   (  5) / 0.295524224714752870173892994651338D0 /
!
  DATA X2    (  1) / 0.995657163025808080735527280689003D0 /
  DATA X2    (  2) / 0.930157491355708226001207180059508D0 /
  DATA X2    (  3) / 0.780817726586416897063717578345042D0 /
  DATA X2    (  4) / 0.562757134668604683339000099272694D0 /
  DATA X2    (  5) / 0.294392862701460198131126603103866D0 /
  DATA W21A  (  1) / 0.032558162307964727478818972459390D0 /
  DATA W21A  (  2) / 0.075039674810919952767043140916190D0 /
  DATA W21A  (  3) / 0.109387158802297641899210590325805D0 /
  DATA W21A  (  4) / 0.134709217311473325928054001771707D0 /
  DATA W21A  (  5) / 0.147739104901338491374841515972068D0 /
  DATA W21B  (  1) / 0.011694638867371874278064396062192D0 /
  DATA W21B  (  2) / 0.054755896574351996031381300244580D0 /
  DATA W21B  (  3) / 0.093125454583697605535065465083366D0 /
  DATA W21B  (  4) / 0.123491976262065851077958109831074D0 /
  DATA W21B  (  5) / 0.142775938577060080797094273138717D0 /
  DATA W21B  (  6) / 0.149445554002916905664936468389821D0 /
!
  DATA X3    (  1) / 0.999333360901932081394099323919911D0 /
  DATA X3    (  2) / 0.987433402908088869795961478381209D0 /
  DATA X3    (  3) / 0.954807934814266299257919200290473D0 /
  DATA X3    (  4) / 0.900148695748328293625099494069092D0 /
  DATA X3    (  5) / 0.825198314983114150847066732588520D0 /
  DATA X3    (  6) / 0.732148388989304982612354848755461D0 /
  DATA X3    (  7) / 0.622847970537725238641159120344323D0 /
  DATA X3    (  8) / 0.499479574071056499952214885499755D0 /
  DATA X3    (  9) / 0.364901661346580768043989548502644D0 /
  DATA X3    ( 10) / 0.222254919776601296498260928066212D0 /
  DATA X3    ( 11) / 0.074650617461383322043914435796506D0 /
  DATA W43A  (  1) / 0.016296734289666564924281974617663D0 /
  DATA W43A  (  2) / 0.037522876120869501461613795898115D0 /
  DATA W43A  (  3) / 0.054694902058255442147212685465005D0 /
  DATA W43A  (  4) / 0.067355414609478086075553166302174D0 /
  DATA W43A  (  5) / 0.073870199632393953432140695251367D0 /
  DATA W43A  (  6) / 0.005768556059769796184184327908655D0 /
  DATA W43A  (  7) / 0.027371890593248842081276069289151D0 /
  DATA W43A  (  8) / 0.046560826910428830743339154433824D0 /
  DATA W43A  (  9) / 0.061744995201442564496240336030883D0 /
  DATA W43A  ( 10) / 0.071387267268693397768559114425516D0 /
  DATA W43B  (  1) / 0.001844477640212414100389106552965D0 /
  DATA W43B  (  2) / 0.010798689585891651740465406741293D0 /
  DATA W43B  (  3) / 0.021895363867795428102523123075149D0 /
  DATA W43B  (  4) / 0.032597463975345689443882222526137D0 /
  DATA W43B  (  5) / 0.042163137935191811847627924327955D0 /
  DATA W43B  (  6) / 0.050741939600184577780189020092084D0 /
  DATA W43B  (  7) / 0.058379395542619248375475369330206D0 /
  DATA W43B  (  8) / 0.064746404951445885544689259517511D0 /
  DATA W43B  (  9) / 0.069566197912356484528633315038405D0 /
  DATA W43B  ( 10) / 0.072824441471833208150939535192842D0 /
  DATA W43B  ( 11) / 0.074507751014175118273571813842889D0 /
  DATA W43B  ( 12) / 0.074722147517403005594425168280423D0 /
!
  DATA X4    (  1) / 0.999902977262729234490529830591582D0 /
  DATA X4    (  2) / 0.997989895986678745427496322365960D0 /
  DATA X4    (  3) / 0.992175497860687222808523352251425D0 /
  DATA X4    (  4) / 0.981358163572712773571916941623894D0 /
  DATA X4    (  5) / 0.965057623858384619128284110607926D0 /
  DATA X4    (  6) / 0.943167613133670596816416634507426D0 /
  DATA X4    (  7) / 0.915806414685507209591826430720050D0 /
  DATA X4    (  8) / 0.883221657771316501372117548744163D0 /
  DATA X4    (  9) / 0.845710748462415666605902011504855D0 /
  DATA X4    ( 10) / 0.803557658035230982788739474980964D0 /
  DATA X4    ( 11) / 0.757005730685495558328942793432020D0 /
  DATA X4    ( 12) / 0.706273209787321819824094274740840D0 /
  DATA X4    ( 13) / 0.651589466501177922534422205016736D0 /
  DATA X4    ( 14) / 0.593223374057961088875273770349144D0 /
  DATA X4    ( 15) / 0.531493605970831932285268948562671D0 /
  DATA X4    ( 16) / 0.466763623042022844871966781659270D0 /
  DATA X4    ( 17) / 0.399424847859218804732101665817923D0 /
  DATA X4    ( 18) / 0.329874877106188288265053371824597D0 /
  DATA X4    ( 19) / 0.258503559202161551802280975429025D0 /
  DATA X4    ( 20) / 0.185695396568346652015917141167606D0 /
  DATA X4    ( 21) / 0.111842213179907468172398359241362D0 /
  DATA X4    ( 22) / 0.037352123394619870814998165437704D0 /
  DATA W87A  (  1) / 0.008148377384149172900002878448190D0 /
  DATA W87A  (  2) / 0.018761438201562822243935059003794D0 /
  DATA W87A  (  3) / 0.027347451050052286161582829741283D0 /
  DATA W87A  (  4) / 0.033677707311637930046581056957588D0 /
  DATA W87A  (  5) / 0.036935099820427907614589586742499D0 /
  DATA W87A  (  6) / 0.002884872430211530501334156248695D0 /
  DATA W87A  (  7) / 0.013685946022712701888950035273128D0 /
  DATA W87A  (  8) / 0.023280413502888311123409291030404D0 /
  DATA W87A  (  9) / 0.030872497611713358675466394126442D0 /
  DATA W87A  ( 10) / 0.035693633639418770719351355457044D0 /
  DATA W87A  ( 11) / 0.000915283345202241360843392549948D0 /
  DATA W87A  ( 12) / 0.005399280219300471367738743391053D0 /
  DATA W87A  ( 13) / 0.010947679601118931134327826856808D0 /
  DATA W87A  ( 14) / 0.016298731696787335262665703223280D0 /
  DATA W87A  ( 15) / 0.021081568889203835112433060188190D0 /
  DATA W87A  ( 16) / 0.025370969769253827243467999831710D0 /
  DATA W87A  ( 17) / 0.029189697756475752501446154084920D0 /
  DATA W87A  ( 18) / 0.032373202467202789685788194889595D0 /
  DATA W87A  ( 19) / 0.034783098950365142750781997949596D0 /
  DATA W87A  ( 20) / 0.036412220731351787562801163687577D0 /
  DATA W87A  ( 21) / 0.037253875503047708539592001191226D0 /
  DATA W87B  (  1) / 0.000274145563762072350016527092881D0 /
  DATA W87B  (  2) / 0.001807124155057942948341311753254D0 /
  DATA W87B  (  3) / 0.004096869282759164864458070683480D0 /
  DATA W87B  (  4) / 0.006758290051847378699816577897424D0 /
  DATA W87B  (  5) / 0.009549957672201646536053581325377D0 /
  DATA W87B  (  6) / 0.012329447652244853694626639963780D0 /
  DATA W87B  (  7) / 0.015010447346388952376697286041943D0 /
  DATA W87B  (  8) / 0.017548967986243191099665352925900D0 /
  DATA W87B  (  9) / 0.019938037786440888202278192730714D0 /
  DATA W87B  ( 10) / 0.022194935961012286796332102959499D0 /
  DATA W87B  ( 11) / 0.024339147126000805470360647041454D0 /
  DATA W87B  ( 12) / 0.026374505414839207241503786552615D0 /
  DATA W87B  ( 13) / 0.028286910788771200659968002987960D0 /
  DATA W87B  ( 14) / 0.030052581128092695322521110347341D0 /
  DATA W87B  ( 15) / 0.031646751371439929404586051078883D0 /
  DATA W87B  ( 16) / 0.033050413419978503290785944862689D0 /
  DATA W87B  ( 17) / 0.034255099704226061787082821046821D0 /
  DATA W87B  ( 18) / 0.035262412660156681033782717998428D0 /
  DATA W87B  ( 19) / 0.036076989622888701185500318003895D0 /
  DATA W87B  ( 20) / 0.036698604498456094498018047441094D0 /
  DATA W87B  ( 21) / 0.037120549269832576114119958413599D0 /
  DATA W87B  ( 22) / 0.037334228751935040321235449094698D0 /
  DATA W87B  ( 23) / 0.037361073762679023410321241766599D0 /
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           CENTR  - MID POINT OF THE INTEGRATION INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL
!           FCENTR - FUNCTION VALUE AT MID POINT
!           ABSC   - ABSCISSA
!           FVAL   - FUNCTION VALUE
!           SAVFUN - ARRAY OF FUNCTION VALUES WHICH HAVE ALREADY BEEN
!                    COMPUTED
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
!***FIRST EXECUTABLE STATEMENT  DQNG
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
!
!           TEST ON VALIDITY OF PARAMETERS
!           ------------------------------
!
  RESULT = 0.0D+00
  ABSERR = 0.0D+00
  NEVAL = 0
  IER = 6
  if ( EPSABS <= 0.0D+00.AND.EPSREL < MAX(0.5D+02*EPMACH,0.5D-28)) &
    go to 80
  HLGTH = 0.5D+00*(B-A)
  DHLGTH = ABS(HLGTH)
  CENTR = 0.5D+00*(B+A)
  FCENTR = F(CENTR)
  NEVAL = 21
  IER = 1
!
!          COMPUTE THE INTEGRAL USING THE 10- AND 21-POINT FORMULA.
!
  DO 70 L = 1,3
  go to (5,25,45),L
    5 RES10 = 0.0D+00
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
  RESKH = 0.5D+00*RES21
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
   65 if ( RESASC /= 0.0D+00.AND.ABSERR /= 0.0D+00) &
    ABSERR = RESASC*MIN(0.1D+01,(0.2D+03*ABSERR/RESASC)**1.5D+00)
  if (RESABS > UFLOW/(0.5D+02*EPMACH)) ABSERR = MAX &
    ((EPMACH*0.5D+02)*RESABS,ABSERR)
  if (ABSERR <= MAX(EPSABS,EPSREL*ABS(RESULT))) IER = 0
! ***JUMP OUT OF DO-LOOP
  if (IER == 0) go to 999
   70 CONTINUE
   80 call XERMSG ('SLATEC', 'DQNG', 'ABNORMAL RETURN', IER, 0)
  999 RETURN
end