subroutine DQAGPE (F, A, B, NPTS2, POINTS, EPSABS, EPSREL, LIMIT, &
     RESULT, ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, PTS, &
     IORD, LEVEL, NDIN, LAST)
!
!! DQAGPE approximates the integral of a function with singularities.
!
!  Approximate a given definite integral I = Integral of F
!            over (A,B), hopefully satisfying the accuracy claim:
!                 ABS(I-RESULT) <= MAX(EPSABS,EPSREL*ABS(I)).
!            Break points of the integration interval, where local
!            difficulties of the integrand may occur (e.g. singularities
!            or discontinuities) are provided by the user.
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A2A1
!***TYPE      DOUBLE PRECISION (QAGPE-S, DQAGPE-D)
!***KEYWORDS  AUTOMATIC INTEGRATOR, EXTRAPOLATION, GENERAL-PURPOSE,
!             GLOBALLY ADAPTIVE, QUADPACK, QUADRATURE,
!             SINGULARITIES AT USER SPECIFIED POINTS
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Computation of a definite integral
!        Standard fortran subroutine
!        Double precision version
!
!        PARAMETERS
!         ON ENTRY
!            F      - Double precision
!                     Function subprogram defining the integrand
!                     function F(X). The actual name for F needs to be
!                     declared E X T E R N A L in the driver program.
!
!            A      - Double precision
!                     Lower limit of integration
!
!            B      - Double precision
!                     Upper limit of integration
!
!            NPTS2  - Integer
!                     Number equal to two more than the number of
!                     user-supplied break points within the integration
!                     range, NPTS2 >= 2.
!                     If NPTS2 < 2, the routine will end with IER = 6.
!
!            POINTS - Double precision
!                     Vector of dimension NPTS2, the first (NPTS2-2)
!                     elements of which are the user provided break
!                     POINTS. If these POINTS do not constitute an
!                     ascending sequence there will be an automatic
!                     sorting.
!
!            EPSABS - Double precision
!                     Absolute accuracy requested
!            EPSREL - Double precision
!                     Relative accuracy requested
!                     If  EPSABS <= 0
!                     and EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28),
!                     the routine will end with IER = 6.
!
!            LIMIT  - Integer
!                     Gives an upper bound on the number of subintervals
!                     in the partition of (A,B), LIMIT >= NPTS2
!                     If LIMIT < NPTS2, the routine will end with
!                     IER = 6.
!
!         ON RETURN
!            RESULT - Double precision
!                     Approximation to the integral
!
!            ABSERR - Double precision
!                     Estimate of the modulus of the absolute error,
!                     which should equal or exceed ABS(I-RESULT)
!
!            NEVAL  - Integer
!                     Number of integrand evaluations
!
!            IER    - Integer
!                     IER = 0 Normal and reliable termination of the
!                             routine. It is assumed that the requested
!                             accuracy has been achieved.
!                     IER > 0 Abnormal termination of the routine.
!                             The estimates for integral and error are
!                             less reliable. It is assumed that the
!                             requested accuracy has not been achieved.
!            ERROR MESSAGES
!                     IER = 1 Maximum number of subdivisions allowed
!                             has been achieved. One can allow more
!                             subdivisions by increasing the value of
!                             LIMIT (and taking the according dimension
!                             adjustments into account). However, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. If
!                             the position of a local difficulty can be
!                             determined (i.e. SINGULARITY,
!                             DISCONTINUITY within the interval), it
!                             should be supplied to the routine as an
!                             element of the vector points. If necessary
!                             an appropriate special-purpose integrator
!                             must be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 The occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             The error may be under-estimated.
!                         = 3 Extremely bad integrand behaviour occurs
!                             At some points of the integration
!                             interval.
!                         = 4 The algorithm does not converge.
!                             Roundoff error is detected in the
!                             extrapolation table. It is presumed that
!                             the requested tolerance cannot be
!                             achieved, and that the returned result is
!                             the best which can be obtained.
!                         = 5 The integral is probably divergent, or
!                             slowly convergent. It must be noted that
!                             divergence can occur with any other value
!                             of IER > 0.
!                         = 6 The input is invalid because
!                             NPTS2 < 2 or
!                             Break points are specified outside
!                             the integration range or
!                             (EPSABS <= 0 and
!                              EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28))
!                             or LIMIT < NPTS2.
!                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
!                             and ELIST(1) are set to zero. ALIST(1) and
!                             BLIST(1) are set to A and B respectively.
!
!            ALIST  - Double precision
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the left end points
!                     of the subintervals in the partition of the given
!                     integration range (A,B)
!
!            BLIST  - Double precision
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the right end points
!                     of the subintervals in the partition of the given
!                     integration range (A,B)
!
!            RLIST  - Double precision
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the integral
!                     approximations on the subintervals
!
!            ELIST  - Double precision
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the moduli of the
!                     absolute error estimates on the subintervals
!
!            PTS    - Double precision
!                     Vector of dimension at least NPTS2, containing the
!                     integration limits and the break points of the
!                     interval in ascending sequence.
!
!            LEVEL  - Integer
!                     Vector of dimension at least LIMIT, containing the
!                     subdivision levels of the subinterval, i.e. if
!                     (AA,BB) is a subinterval of (P1,P2) where P1 as
!                     well as P2 is a user-provided break point or
!                     integration limit, then (AA,BB) has level L if
!                     ABS(BB-AA) = ABS(P2-P1)*2**(-L).
!
!            NDIN   - Integer
!                     Vector of dimension at least NPTS2, after first
!                     integration over the intervals (PTS(I)),PTS(I+1),
!                     I = 0,1, ..., NPTS2-2, the error estimates over
!                     some of the intervals may have been increased
!                     artificially, in order to put their subdivision
!                     forward. If this happens for the subinterval
!                     numbered K, NDIN(K) is put to 1, otherwise
!                     NDIN(K) = 0.
!
!            IORD   - Integer
!                     Vector of dimension at least LIMIT, the first K
!                     elements of which are pointers to the
!                     error estimates over the subintervals,
!                     such that ELIST(IORD(1)), ..., ELIST(IORD(K))
!                     form a decreasing sequence, with K = LAST
!                     If LAST <= (LIMIT/2+2), and K = LIMIT+1-LAST
!                     otherwise
!
!            LAST   - Integer
!                     Number of subintervals actually produced in the
!                     subdivisions process
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DQELG, DQK21, DQPSRT
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DQAGPE
  DOUBLE PRECISION A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1, &
    A2,B,BLIST,B1,B2,CORREC,DEFABS,DEFAB1,DEFAB2, &
    DRES,D1MACH,ELIST,EPMACH,EPSABS,EPSREL,ERLARG,ERLAST,ERRBND, &
    ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,ERTEST,F,OFLOW,POINTS,PTS, &
    RESA,RESABS,RESEPS,RESULT,RES3LA,RLIST,RLIST2,SIGN,TEMP,UFLOW
  INTEGER I,ID,IER,IERRO,IND1,IND2,IORD,IP1,IROFF1,IROFF2,IROFF3,J, &
    JLOW,JUPBND,K,KSGN,KTMIN,LAST,LEVCUR,LEVEL,LEVMAX,LIMIT,MAXERR, &
    NDIN,NEVAL,NINT,NINTP1,NPTS,NPTS2,NRES,NRMAX,NUMRL2
  LOGICAL EXTRAP,NOEXT
!
!
  DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*), &
    LEVEL(*),NDIN(*),POINTS(*),PTS(*),RES3LA(3), &
    RLIST(*),RLIST2(52)
!
  EXTERNAL F
!
!            THE DIMENSION OF RLIST2 IS DETERMINED BY THE VALUE OF
!            LIMEXP IN SUBROUTINE EPSALG (RLIST2 SHOULD BE OF DIMENSION
!            (LIMEXP+2) AT LEAST).
!
!
!            LIST OF MAJOR VARIABLES
!            -----------------------
!
!           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
!                       CONSIDERED UP TO NOW
!           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
!                       CONSIDERED UP TO NOW
!           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
!                       (ALIST(I),BLIST(I))
!           RLIST2    - ARRAY OF DIMENSION AT LEAST LIMEXP+2
!                       CONTAINING THE PART OF THE EPSILON TABLE WHICH
!                       IS STILL NEEDED FOR FURTHER COMPUTATIONS
!           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
!           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST ERROR
!                       ESTIMATE
!           ERRMAX    - ELIST(MAXERR)
!           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
!                       (BEFORE THAT SUBDIVISION HAS TAKEN PLACE)
!           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
!           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
!           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
!                       ABS(RESULT))
!           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
!           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
!           LAST      - INDEX FOR SUBDIVISION
!           NRES      - NUMBER OF CALLS TO THE EXTRAPOLATION ROUTINE
!           NUMRL2    - NUMBER OF ELEMENTS IN RLIST2. if AN APPROPRIATE
!                       APPROXIMATION TO THE COMPOUNDED INTEGRAL HAS
!                       BEEN OBTAINED, IT IS PUT IN RLIST2(NUMRL2) AFTER
!                       NUMRL2 HAS BEEN INCREASED BY ONE.
!           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
!                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
!           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE
!                       IS ATTEMPTING TO PERFORM EXTRAPOLATION. I.E.
!                       BEFORE SUBDIVIDING THE SMALLEST INTERVAL WE
!                       TRY TO DECREASE THE VALUE OF ERLARG.
!           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION IS
!                       NO LONGER ALLOWED (TRUE-VALUE)
!
!            MACHINE DEPENDENT CONSTANTS
!            ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  DQAGPE
  EPMACH = D1MACH(4)
!
!            TEST ON VALIDITY OF PARAMETERS
!            -----------------------------
!
  IER = 0
  NEVAL = 0
  LAST = 0
  RESULT = 0.0D+00
  ABSERR = 0.0D+00
  ALIST(1) = A
  BLIST(1) = B
  RLIST(1) = 0.0D+00
  ELIST(1) = 0.0D+00
  IORD(1) = 0
  LEVEL(1) = 0
  NPTS = NPTS2-2
  if ( NPTS2 < 2.OR.LIMIT <= NPTS.OR.(EPSABS <= 0.0D+00.AND. &
    EPSREL < MAX(0.5D+02*EPMACH,0.5D-28))) IER = 6
  if ( IER == 6) go to 999
!
!            if ANY BREAK POINTS ARE PROVIDED, SORT THEM INTO AN
!            ASCENDING SEQUENCE.
!
  SIGN = 1.0D+00
  if ( A > B) SIGN = -1.0D+00
  PTS(1) = MIN(A,B)
  if ( NPTS == 0) go to 15
  DO 10 I = 1,NPTS
    PTS(I+1) = POINTS(I)
   10 CONTINUE
   15 PTS(NPTS+2) = MAX(A,B)
  NINT = NPTS+1
  A1 = PTS(1)
  if ( NPTS == 0) go to 40
  NINTP1 = NINT+1
  DO 20 I = 1,NINT
    IP1 = I+1
    DO 20 J = IP1,NINTP1
      if ( PTS(I) <= PTS(J)) go to 20
      TEMP = PTS(I)
      PTS(I) = PTS(J)
      PTS(J) = TEMP
   20 CONTINUE
  if ( PTS(1) /= MIN(A,B).OR.PTS(NINTP1) /= MAX(A,B)) IER = 6
  if ( IER == 6) go to 999
!
!            COMPUTE FIRST INTEGRAL AND ERROR APPROXIMATIONS.
!            ------------------------------------------------
!
   40 RESABS = 0.0D+00
  DO 50 I = 1,NINT
    B1 = PTS(I+1)
    call DQK21(F,A1,B1,AREA1,ERROR1,DEFABS,RESA)
    ABSERR = ABSERR+ERROR1
    RESULT = RESULT+AREA1
    NDIN(I) = 0
    if ( ERROR1 == RESA.AND.ERROR1 /= 0.0D+00) NDIN(I) = 1
    RESABS = RESABS+DEFABS
    LEVEL(I) = 0
    ELIST(I) = ERROR1
    ALIST(I) = A1
    BLIST(I) = B1
    RLIST(I) = AREA1
    IORD(I) = I
    A1 = B1
   50 CONTINUE
  ERRSUM = 0.0D+00
  DO 55 I = 1,NINT
    if ( NDIN(I) == 1) ELIST(I) = ABSERR
    ERRSUM = ERRSUM+ELIST(I)
   55 CONTINUE
!
!           TEST ON ACCURACY.
!
  LAST = NINT
  NEVAL = 21*NINT
  DRES = ABS(RESULT)
  ERRBND = MAX(EPSABS,EPSREL*DRES)
  if ( ABSERR <= 0.1D+03*EPMACH*RESABS.AND.ABSERR > ERRBND) IER = 2
  if ( NINT == 1) go to 80
  DO 70 I = 1,NPTS
    JLOW = I+1
    IND1 = IORD(I)
    DO 60 J = JLOW,NINT
      IND2 = IORD(J)
      if ( ELIST(IND1) > ELIST(IND2)) go to 60
      IND1 = IND2
      K = J
   60   CONTINUE
    if ( IND1 == IORD(I)) go to 70
    IORD(K) = IORD(I)
    IORD(I) = IND1
   70 CONTINUE
  if ( LIMIT < NPTS2) IER = 1
   80 if ( IER /= 0.OR.ABSERR <= ERRBND) go to 999
!
!           INITIALIZATION
!           --------------
!
  RLIST2(1) = RESULT
  MAXERR = IORD(1)
  ERRMAX = ELIST(MAXERR)
  AREA = RESULT
  NRMAX = 1
  NRES = 0
  NUMRL2 = 1
  KTMIN = 0
  EXTRAP = .FALSE.
  NOEXT = .FALSE.
  ERLARG = ERRSUM
  ERTEST = ERRBND
  LEVMAX = 1
  IROFF1 = 0
  IROFF2 = 0
  IROFF3 = 0
  IERRO = 0
  UFLOW = D1MACH(1)
  OFLOW = D1MACH(2)
  ABSERR = OFLOW
  KSGN = -1
  if ( DRES >= (0.1D+01-0.5D+02*EPMACH)*RESABS) KSGN = 1
!
!           MAIN DO-LOOP
!           ------------
!
  DO 160 LAST = NPTS2,LIMIT
!
!           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST ERROR
!           ESTIMATE.
!
    LEVCUR = LEVEL(MAXERR)+1
    A1 = ALIST(MAXERR)
    B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
    A2 = B1
    B2 = BLIST(MAXERR)
    ERLAST = ERRMAX
    call DQK21(F,A1,B1,AREA1,ERROR1,RESA,DEFAB1)
    call DQK21(F,A2,B2,AREA2,ERROR2,RESA,DEFAB2)
!
!           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
!           AND ERROR AND TEST FOR ACCURACY.
!
    NEVAL = NEVAL+42
    AREA12 = AREA1+AREA2
    ERRO12 = ERROR1+ERROR2
    ERRSUM = ERRSUM+ERRO12-ERRMAX
    AREA = AREA+AREA12-RLIST(MAXERR)
    if ( DEFAB1 == ERROR1.OR.DEFAB2 == ERROR2) go to 95
    if ( ABS(RLIST(MAXERR)-AREA12) > 0.1D-04*ABS(AREA12) &
    .OR.ERRO12 < 0.99D+00*ERRMAX) go to 90
    if ( EXTRAP) IROFF2 = IROFF2+1
    if ( .NOT.EXTRAP) IROFF1 = IROFF1+1
   90   if ( LAST > 10.AND.ERRO12 > ERRMAX) IROFF3 = IROFF3+1
   95   LEVEL(MAXERR) = LEVCUR
    LEVEL(LAST) = LEVCUR
    RLIST(MAXERR) = AREA1
    RLIST(LAST) = AREA2
    ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
!
!           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
!
    if ( IROFF1+IROFF2 >= 10.OR.IROFF3 >= 20) IER = 2
    if ( IROFF2 >= 5) IERRO = 3
!
!           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
!           SUBINTERVALS EQUALS LIMIT.
!
    if ( LAST == LIMIT) IER = 1
!
!           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
!           AT A POINT OF THE INTEGRATION RANGE
!
    if ( MAX(ABS(A1),ABS(B2)) <= (0.1D+01+0.1D+03*EPMACH)* &
    (ABS(A2)+0.1D+04*UFLOW)) IER = 4
!
!           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
!
    if ( ERROR2 > ERROR1) go to 100
    ALIST(LAST) = A2
    BLIST(MAXERR) = B1
    BLIST(LAST) = B2
    ELIST(MAXERR) = ERROR1
    ELIST(LAST) = ERROR2
    go to 110
  100   ALIST(MAXERR) = A2
    ALIST(LAST) = A1
    BLIST(LAST) = B1
    RLIST(MAXERR) = AREA2
    RLIST(LAST) = AREA1
    ELIST(MAXERR) = ERROR2
    ELIST(LAST) = ERROR1
!
!           call SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
!           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
!           WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
!
  110   call DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
! ***JUMP OUT OF DO-LOOP
    if ( ERRSUM <= ERRBND) go to 190
! ***JUMP OUT OF DO-LOOP
    if ( IER /= 0) go to 170
    if ( NOEXT) go to 160
    ERLARG = ERLARG-ERLAST
    if ( LEVCUR+1 <= LEVMAX) ERLARG = ERLARG+ERRO12
    if ( EXTRAP) go to 120
!
!           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
!           SMALLEST INTERVAL.
!
    if ( LEVEL(MAXERR)+1 <= LEVMAX) go to 160
    EXTRAP = .TRUE.
    NRMAX = 2
  120   if ( IERRO == 3.OR.ERLARG <= ERTEST) go to 140
!
!           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
!           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS OVER
!           THE LARGER INTERVALS (ERLARG) AND PERFORM EXTRAPOLATION.
!
    ID = NRMAX
    JUPBND = LAST
    if ( LAST > (2+LIMIT/2)) JUPBND = LIMIT+3-LAST
    DO 130 K = ID,JUPBND
      MAXERR = IORD(NRMAX)
      ERRMAX = ELIST(MAXERR)
! ***JUMP OUT OF DO-LOOP
      if ( LEVEL(MAXERR)+1 <= LEVMAX) go to 160
      NRMAX = NRMAX+1
  130   CONTINUE
!
!           PERFORM EXTRAPOLATION.
!
  140   NUMRL2 = NUMRL2+1
    RLIST2(NUMRL2) = AREA
    if ( NUMRL2 <= 2) go to 155
    call DQELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
    KTMIN = KTMIN+1
    if ( KTMIN > 5.AND.ABSERR < 0.1D-02*ERRSUM) IER = 5
    if ( ABSEPS >= ABSERR) go to 150
    KTMIN = 0
    ABSERR = ABSEPS
    RESULT = RESEPS
    CORREC = ERLARG
    ERTEST = MAX(EPSABS,EPSREL*ABS(RESEPS))
! ***JUMP OUT OF DO-LOOP
    if ( ABSERR < ERTEST) go to 170
!
!           PREPARE BISECTION OF THE SMALLEST INTERVAL.
!
  150   if ( NUMRL2 == 1) NOEXT = .TRUE.
    if ( IER >= 5) go to 170
  155   MAXERR = IORD(1)
    ERRMAX = ELIST(MAXERR)
    NRMAX = 1
    EXTRAP = .FALSE.
    LEVMAX = LEVMAX+1
    ERLARG = ERRSUM
  160 CONTINUE
!
!           SET THE FINAL RESULT.
!           ---------------------
!
!
  170 if ( ABSERR == OFLOW) go to 190
  if ( (IER+IERRO) == 0) go to 180
  if ( IERRO == 3) ABSERR = ABSERR+CORREC
  if ( IER == 0) IER = 3
  if ( RESULT /= 0.0D+00.AND.AREA /= 0.0D+00)go to 175
  if ( ABSERR > ERRSUM)go to 190
  if ( AREA == 0.0D+00) go to 210
  go to 180
  175 if ( ABSERR/ABS(RESULT) > ERRSUM/ABS(AREA))go to 190
!
!           TEST ON DIVERGENCE.
!
  180 if ( KSGN == (-1).AND.MAX(ABS(RESULT),ABS(AREA)) <=  &
    DEFABS*0.1D-01) go to 210
  if ( 0.1D-01 > (RESULT/AREA).OR.(RESULT/AREA) > 0.1D+03.OR. &
    ERRSUM > ABS(AREA)) IER = 6
  go to 210
!
!           COMPUTE GLOBAL INTEGRAL SUM.
!
  190 RESULT = 0.0D+00
  DO 200 K = 1,LAST
    RESULT = RESULT+RLIST(K)
  200 CONTINUE
  ABSERR = ERRSUM
  210 if ( IER > 2) IER = IER-1
  RESULT = RESULT*SIGN
  999 RETURN
end
