subroutine QAWOE (F, A, B, OMEGA, INTEGR, EPSABS, EPSREL, LIMIT, &
     ICALL, MAXP1, RESULT, ABSERR, NEVAL, IER, LAST, ALIST, BLIST, &
     RLIST, ELIST, IORD, NNLOG, MOMCOM, CHEBMO)
!
!! QAWOE calculates an approximation to a given definite integral
!               I = Integral of F(X)*W(X) over (A,B), where
!                  W(X) = COS(OMEGA*X)
!               or W(X) = SIN(OMEGA*X),
!            hopefully satisfying the following claim for accuracy
!               ABS(I-RESULT) <= MAX(EPSABS,EPSREL*ABS(I)).
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A2A1
!***TYPE      SINGLE PRECISION (QAWOE-S, DQAWOE-D)
!***KEYWORDS  AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD,
!             EXTRAPOLATION, GLOBALLY ADAPTIVE,
!             INTEGRAND WITH OSCILLATORY COS OR SIN FACTOR, QUADPACK,
!             QUADRATURE, SPECIAL-PURPOSE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Computation of Oscillatory integrals
!        Standard fortran subroutine
!        Real version
!
!        PARAMETERS
!         ON ENTRY
!            F      - Real
!                     Function subprogram defining the integrand
!                     function F(X). The actual name for F needs to be
!                     declared E X T E R N A L in the driver program.
!
!            A      - Real
!                     Lower limit of integration
!
!            B      - Real
!                     Upper limit of integration
!
!            OMEGA  - Real
!                     Parameter in the integrand weight function
!
!            INTEGR - Integer
!                     Indicates which of the WEIGHT functions is to be
!                     used
!                     INTEGR = 1      W(X) = COS(OMEGA*X)
!                     INTEGR = 2      W(X) = SIN(OMEGA*X)
!                     If INTEGR /= 1 and INTEGR /= 2, the routine
!                     will end with IER = 6.
!
!            EPSABS - Real
!                     Absolute accuracy requested
!            EPSREL - Real
!                     Relative accuracy requested
!                     If  EPSABS <= 0
!                     and EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28),
!                     the routine will end with IER = 6.
!
!            LIMIT  - Integer
!                     Gives an upper bound on the number of subdivisions
!                     in the partition of (A,B), LIMIT >= 1.
!
!            ICALL  - Integer
!                     If QAWOE is to be used only once, ICALL must
!                     be set to 1.  Assume that during this call, the
!                     Chebyshev moments (for CLENSHAW-CURTIS integration
!                     of degree 24) have been computed for intervals of
!                     lengths (ABS(B-A))*2**(-L), L=0,1,2,...MOMCOM-1.
!                     If ICALL > 1 this means that QAWOE has been
!                     called twice or more on intervals of the same
!                     length ABS(B-A). The Chebyshev moments already
!                     computed are then re-used in subsequent calls.
!                     If ICALL < 1, the routine will end with IER = 6.
!
!            MAXP1  - Integer
!                     Gives an upper bound on the number of Chebyshev
!                     moments which can be stored, i.e. for the
!                     intervals of lengths ABS(B-A)*2**(-L),
!                     L=0,1, ..., MAXP1-2, MAXP1 >= 1.
!                     If MAXP1 < 1, the routine will end with IER = 6.
!
!         ON RETURN
!            RESULT - Real
!                     Approximation to the integral
!
!            ABSERR - Real
!                     Estimate of the modulus of the absolute error,
!                     which should equal or exceed ABS(I-RESULT)
!
!            NEVAL  - Integer
!                     Number of integrand evaluations
!
!            IER    - Integer
!                     IER = 0 Normal and reliable termination of the
!                             routine. It is assumed that the
!                             requested accuracy has been achieved.
!                   - IER > 0 Abnormal termination of the routine.
!                             The estimates for integral and error are
!                             less reliable. It is assumed that the
!                             requested accuracy has not been achieved.
!            ERROR MESSAGES
!                     IER = 1 Maximum number of subdivisions allowed
!                             has been achieved. One can allow more
!                             subdivisions by increasing the value of
!                             LIMIT (and taking according dimension
!                             adjustments into account). However, if
!                             this yields no improvement it is advised
!                             to analyze the integrand, in order to
!                             determine the integration difficulties.
!                             If the position of a local difficulty can
!                             be determined (e.g. SINGULARITY,
!                             DISCONTINUITY within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. If possible,
!                             an appropriate special-purpose integrator
!                             should be used which is designed for
!                             handling the type of difficulty involved.
!                         = 2 The occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                             The error may be under-estimated.
!                         = 3 Extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 4 The algorithm does not converge.
!                             Roundoff error is detected in the
!                             extrapolation table.
!                             It is presumed that the requested
!                             tolerance cannot be achieved due to
!                             roundoff in the extrapolation table,
!                             and that the returned result is the
!                             best which can be obtained.
!                         = 5 The integral is probably divergent, or
!                             slowly convergent. It must be noted that
!                             divergence can occur with any other value
!                             of IER > 0.
!                         = 6 The input is invalid, because
!                             (EPSABS <= 0 and
!                              EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28))
!                             or (INTEGR /= 1 and INTEGR /= 2) or
!                             ICALL < 1 or MAXP1 < 1.
!                             RESULT, ABSERR, NEVAL, LAST, RLIST(1),
!                             ELIST(1), IORD(1) and NNLOG(1) are set
!                             to ZERO. ALIST(1) and BLIST(1) are set
!                             to A and B respectively.
!
!            LAST  -  Integer
!                     On return, LAST equals the number of
!                     subintervals produces in the subdivision
!                     process, which determines the number of
!                     significant elements actually in the
!                     WORK ARRAYS.
!            ALIST  - Real
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the left
!                     end points of the subintervals in the partition
!                     of the given integration range (A,B)
!
!            BLIST  - Real
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the right
!                     end points of the subintervals in the partition
!                     of the given integration range (A,B)
!
!            RLIST  - Real
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the integral
!                     approximations on the subintervals
!
!            ELIST  - Real
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the moduli of the
!                     absolute error estimates on the subintervals
!
!            IORD   - Integer
!                     Vector of dimension at least LIMIT, the first K
!                     elements of which are pointers to the error
!                     estimates over the subintervals,
!                     such that ELIST(IORD(1)), ...,
!                     ELIST(IORD(K)) form a decreasing sequence, with
!                     K = LAST if LAST <= (LIMIT/2+2), and
!                     K = LIMIT+1-LAST otherwise.
!
!            NNLOG  - Integer
!                     Vector of dimension at least LIMIT, containing the
!                     subdivision levels of the subintervals, i.e.
!                     IWORK(I) = L means that the subinterval
!                     numbered I is of length ABS(B-A)*2**(1-L)
!
!         ON ENTRY AND RETURN
!            MOMCOM - Integer
!                     Indicating that the Chebyshev moments
!                     have been computed for intervals of lengths
!                     (ABS(B-A))*2**(-L), L=0,1,2, ..., MOMCOM-1,
!                     MOMCOM < MAXP1
!
!            CHEBMO - Real
!                     Array of dimension (MAXP1,25) containing the
!                     Chebyshev moments
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  QC25F, QELG, QPSRT, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  QAWOE
!
  REAL A,ABSEPS,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1, &
    A2,B,BLIST,B1,B2,CHEBMO,CORREC,DEFAB1,DEFAB2,DEFABS, &
    DOMEGA,R1MACH,DRES,ELIST,EPMACH,EPSABS,EPSREL,ERLARG, &
    ERLAST,ERRBND,ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,ERTEST,F,OFLOW, &
    OMEGA,RESABS,RESEPS,RESULT,RES3LA,RLIST,RLIST2,SMALL,UFLOW,WIDTH
  INTEGER ICALL,ID,IER,IERRO,INTEGR,IORD,IROFF1,IROFF2,IROFF3, &
    JUPBND,K,KSGN,KTMIN,LAST,LIMIT,MAXERR,MAXP1,MOMCOM,NEV, &
    NEVAL,NNLOG,NRES,NRMAX,NRMOM,NUMRL2
  LOGICAL EXTRAP,NOEXT,EXTALL
!
  DIMENSION ALIST(*),BLIST(*),RLIST(*),ELIST(*), &
    IORD(*),RLIST2(52),RES3LA(3),CHEBMO(MAXP1,25),NNLOG(*)
!
  EXTERNAL F
!
!            THE DIMENSION OF RLIST2 IS DETERMINED BY  THE VALUE OF
!            LIMEXP IN SUBROUTINE QELG (RLIST2 SHOULD BE OF
!            DIMENSION (LIMEXP+2) AT LEAST).
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
!                       CONTAINING THE PART OF THE EPSILON TABLE
!                       WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS
!           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
!           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
!                       ERROR ESTIMATE
!           ERRMAX    - ELIST(MAXERR)
!           ERLAST    - ERROR ON THE INTERVAL CURRENTLY SUBDIVIDED
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
!                       BEEN OBTAINED IT IS PUT IN RLIST2(NUMRL2) AFTER
!                       NUMRL2 HAS BEEN INCREASED BY ONE
!           SMALL     - LENGTH OF THE SMALLEST INTERVAL CONSIDERED
!                       UP TO NOW, MULTIPLIED BY 1.5
!           ERLARG    - SUM OF THE ERRORS OVER THE INTERVALS LARGER
!                       THAN THE SMALLEST INTERVAL CONSIDERED UP TO NOW
!           EXTRAP    - LOGICAL VARIABLE DENOTING THAT THE ROUTINE IS
!                       ATTEMPTING TO PERFORM EXTRAPOLATION, I.E. BEFORE
!                       SUBDIVIDING THE SMALLEST INTERVAL WE TRY TO
!                       DECREASE THE VALUE OF ERLARG
!           NOEXT     - LOGICAL VARIABLE DENOTING THAT EXTRAPOLATION
!                       IS NO LONGER ALLOWED (TRUE VALUE)
!
!            MACHINE DEPENDENT CONSTANTS
!            ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  QAWOE
  EPMACH = R1MACH(4)
!
!         TEST ON VALIDITY OF PARAMETERS
!         ------------------------------
!
  IER = 0
  NEVAL = 0
  LAST = 0
  RESULT = 0.0E+00
  ABSERR = 0.0E+00
  ALIST(1) = A
  BLIST(1) = B
  RLIST(1) = 0.0E+00
  ELIST(1) = 0.0E+00
  IORD(1) = 0
  NNLOG(1) = 0
  if ( (INTEGR /= 1.AND.INTEGR /= 2).OR.(EPSABS <= 0.0E+00.AND. &
    EPSREL < MAX(0.5E+02*EPMACH,0.5E-14)).OR.ICALL < 1.OR. &
    MAXP1 < 1) IER = 6
  if ( IER == 6) go to 999
!
!           FIRST APPROXIMATION TO THE INTEGRAL
!           -----------------------------------
!
  DOMEGA = ABS(OMEGA)
  NRMOM = 0
  if (ICALL > 1) go to 5
  MOMCOM = 0
    5 call QC25F(F,A,B,DOMEGA,INTEGR,NRMOM,MAXP1,0,RESULT,ABSERR, &
    NEVAL,DEFABS,RESABS,MOMCOM,CHEBMO)
!
!           TEST ON ACCURACY.
!
  DRES = ABS(RESULT)
  ERRBND = MAX(EPSABS,EPSREL*DRES)
  RLIST(1) = RESULT
  ELIST(1) = ABSERR
  IORD(1) = 1
  if ( ABSERR <= 0.1E+03*EPMACH*DEFABS.AND.ABSERR >  &
    ERRBND) IER = 2
  if ( LIMIT == 1) IER = 1
  if ( IER /= 0.OR.ABSERR <= ERRBND) go to 200
!
!           INITIALIZATIONS
!           ---------------
!
  UFLOW = R1MACH(1)
  OFLOW = R1MACH(2)
  ERRMAX = ABSERR
  MAXERR = 1
  AREA = RESULT
  ERRSUM = ABSERR
  ABSERR = OFLOW
  NRMAX = 1
  EXTRAP = .FALSE.
  NOEXT = .FALSE.
  IERRO = 0
  IROFF1 = 0
  IROFF2 = 0
  IROFF3 = 0
  KTMIN = 0
  SMALL = ABS(B-A)*0.75E+00
  NRES = 0
  NUMRL2 = 0
  EXTALL = .FALSE.
  if ( 0.5E+00*ABS(B-A)*DOMEGA > 0.2E+01) go to 10
  NUMRL2 = 1
  EXTALL = .TRUE.
  RLIST2(1) = RESULT
   10 if ( 0.25E+00*ABS(B-A)*DOMEGA <= 0.2E+01) EXTALL = .TRUE.
  KSGN = -1
  if ( DRES >= (0.1E+01-0.5E+02*EPMACH)*DEFABS) KSGN = 1
!
!           MAIN DO-LOOP
!           ------------
!
  DO 140 LAST = 2,LIMIT
!
!           BISECT THE SUBINTERVAL WITH THE NRMAX-TH LARGEST
!           ERROR ESTIMATE.
!
    NRMOM = NNLOG(MAXERR)+1
    A1 = ALIST(MAXERR)
    B1 = 0.5E+00*(ALIST(MAXERR)+BLIST(MAXERR))
    A2 = B1
    B2 = BLIST(MAXERR)
    ERLAST = ERRMAX
    call QC25F(F,A1,B1,DOMEGA,INTEGR,NRMOM,MAXP1,0, &
    AREA1,ERROR1,NEV,RESABS,DEFAB1,MOMCOM,CHEBMO)
    NEVAL = NEVAL+NEV
    call QC25F(F,A2,B2,DOMEGA,INTEGR,NRMOM,MAXP1,1, &
    AREA2,ERROR2,NEV,RESABS,DEFAB2,MOMCOM,CHEBMO)
    NEVAL = NEVAL+NEV
!
!           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
!           AND ERROR AND TEST FOR ACCURACY.
!
    AREA12 = AREA1+AREA2
    ERRO12 = ERROR1+ERROR2
    ERRSUM = ERRSUM+ERRO12-ERRMAX
    AREA = AREA+AREA12-RLIST(MAXERR)
    if ( DEFAB1 == ERROR1.OR.DEFAB2 == ERROR2) go to 25
    if ( ABS(RLIST(MAXERR)-AREA12) > 0.1E-04*ABS(AREA12) &
    .OR.ERRO12 < 0.99E+00*ERRMAX) go to 20
    if ( EXTRAP) IROFF2 = IROFF2+1
    if ( .NOT.EXTRAP) IROFF1 = IROFF1+1
   20   if ( LAST > 10.AND.ERRO12 > ERRMAX) IROFF3 = IROFF3+1
   25   RLIST(MAXERR) = AREA1
    RLIST(LAST) = AREA2
    NNLOG(MAXERR) = NRMOM
    NNLOG(LAST) = NRMOM
    ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
!
!           TEST FOR ROUNDOFF ERROR AND EVENTUALLY
!           SET ERROR FLAG
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
!           AT A POINT OF THE INTEGRATION RANGE.
!
    if ( MAX(ABS(A1),ABS(B2)) <= (0.1E+01+0.1E+03*EPMACH) &
    *(ABS(A2)+0.1E+04*UFLOW)) IER = 4
!
!           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
!
    if ( ERROR2 > ERROR1) go to 30
    ALIST(LAST) = A2
    BLIST(MAXERR) = B1
    BLIST(LAST) = B2
    ELIST(MAXERR) = ERROR1
    ELIST(LAST) = ERROR2
    go to 40
   30   ALIST(MAXERR) = A2
    ALIST(LAST) = A1
    BLIST(LAST) = B1
    RLIST(MAXERR) = AREA2
    RLIST(LAST) = AREA1
    ELIST(MAXERR) = ERROR2
    ELIST(LAST) = ERROR1
!
!           call SUBROUTINE QPSRT TO MAINTAIN THE DESCENDING ORDERING
!           IN THE LIST OF ERROR ESTIMATES AND SELECT THE
!           SUBINTERVAL WITH NRMAX-TH LARGEST ERROR ESTIMATE (TO BE
!           BISECTED NEXT).
!
   40   call QPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
! ***JUMP OUT OF DO-LOOP
  if ( ERRSUM <= ERRBND) go to 170
  if ( IER /= 0) go to 150
    if ( LAST == 2.AND.EXTALL) go to 120
    if ( NOEXT) go to 140
    if ( .NOT.EXTALL) go to 50
    ERLARG = ERLARG-ERLAST
    if ( ABS(B1-A1) > SMALL) ERLARG = ERLARG+ERRO12
    if ( EXTRAP) go to 70
!
!           TEST WHETHER THE INTERVAL TO BE BISECTED NEXT IS THE
!           SMALLEST INTERVAL.
!
   50   WIDTH = ABS(BLIST(MAXERR)-ALIST(MAXERR))
    if ( WIDTH > SMALL) go to 140
    if ( EXTALL) go to 60
!
!           TEST WHETHER WE CAN START WITH THE EXTRAPOLATION
!           PROCEDURE (WE DO THIS if WE INTEGRATE OVER THE
!           NEXT INTERVAL WITH USE OF A GAUSS-KRONROD RULE - SEE
!           SUBROUTINE QC25F).
!
    SMALL = SMALL*0.5E+00
    if ( 0.25E+00*WIDTH*DOMEGA > 0.2E+01) go to 140
    EXTALL = .TRUE.
    go to 130
   60   EXTRAP = .TRUE.
    NRMAX = 2
   70   if ( IERRO == 3.OR.ERLARG <= ERTEST) go to 90
!
!           THE SMALLEST INTERVAL HAS THE LARGEST ERROR.
!           BEFORE BISECTING DECREASE THE SUM OF THE ERRORS
!           OVER THE LARGER INTERVALS (ERLARG) AND PERFORM
!           EXTRAPOLATION.
!
    JUPBND = LAST
    if (LAST > (LIMIT/2+2)) JUPBND = LIMIT+3-LAST
    ID = NRMAX
    DO 80 K = ID,JUPBND
      MAXERR = IORD(NRMAX)
      ERRMAX = ELIST(MAXERR)
      if ( ABS(BLIST(MAXERR)-ALIST(MAXERR)) > SMALL) go to 140
      NRMAX = NRMAX+1
   80   CONTINUE
!
!           PERFORM EXTRAPOLATION.
!
   90   NUMRL2 = NUMRL2+1
    RLIST2(NUMRL2) = AREA
    if ( NUMRL2 < 3) go to 110
    call QELG(NUMRL2,RLIST2,RESEPS,ABSEPS,RES3LA,NRES)
    KTMIN = KTMIN+1
    if ( KTMIN > 5.AND.ABSERR < 0.1E-02*ERRSUM) IER = 5
    if ( ABSEPS >= ABSERR) go to 100
    KTMIN = 0
    ABSERR = ABSEPS
    RESULT = RESEPS
    CORREC = ERLARG
    ERTEST = MAX(EPSABS,EPSREL*ABS(RESEPS))
! ***JUMP OUT OF DO-LOOP
    if ( ABSERR <= ERTEST) go to 150
!
!           PREPARE BISECTION OF THE SMALLEST INTERVAL.
!
  100   if ( NUMRL2 == 1) NOEXT = .TRUE.
    if ( IER == 5) go to 150
  110   MAXERR = IORD(1)
    ERRMAX = ELIST(MAXERR)
    NRMAX = 1
    EXTRAP = .FALSE.
    SMALL = SMALL*0.5E+00
    ERLARG = ERRSUM
    go to 140
  120   SMALL = SMALL*0.5E+00
    NUMRL2 = NUMRL2+1
    RLIST2(NUMRL2) = AREA
  130   ERTEST = ERRBND
    ERLARG = ERRSUM
  140 CONTINUE
!
!           SET THE FINAL RESULT.
!           ---------------------
!
  150 if ( ABSERR == OFLOW.OR.NRES == 0) go to 170
  if ( IER+IERRO == 0) go to 165
  if ( IERRO == 3) ABSERR = ABSERR+CORREC
  if ( IER == 0) IER = 3
  if ( RESULT /= 0.0E+00.AND.AREA /= 0.0E+00) go to 160
  if ( ABSERR > ERRSUM) go to 170
  if ( AREA == 0.0E+00) go to 190
  go to 165
  160 if ( ABSERR/ABS(RESULT) > ERRSUM/ABS(AREA)) go to 170
!
!           TEST ON DIVERGENCE.
!
  165 if ( KSGN == (-1).AND.MAX(ABS(RESULT),ABS(AREA)) <=  &
   DEFABS*0.1E-01) go to 190
  if ( 0.1E-01 > (RESULT/AREA).OR.(RESULT/AREA) > 0.1E+03 &
   .OR.ERRSUM >= ABS(AREA)) IER = 6
  go to 190
!
!           COMPUTE GLOBAL INTEGRAL SUM.
!
  170 RESULT = 0.0E+00
  DO 180 K=1,LAST
    RESULT = RESULT+RLIST(K)
  180 CONTINUE
  ABSERR = ERRSUM
  190 if (IER > 2) IER=IER-1
  200 if (INTEGR == 2.AND.OMEGA < 0.0E+00) RESULT=-RESULT
  999 RETURN
end
