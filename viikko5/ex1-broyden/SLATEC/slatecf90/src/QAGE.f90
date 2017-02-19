subroutine QAGE (F, A, B, EPSABS, EPSREL, KEY, LIMIT, RESULT, &
     ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
!
!! QAGE calculates an approximation RESULT to a given definite integral
!  I = Integral of F over (A,B),
!            hopefully satisfying following claim for accuracy
!            ABS(I-RESLT) <= MAX(EPSABS,EPSREL*ABS(I)).
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A1
!***TYPE      SINGLE PRECISION (QAGE-S, DQAGE-D)
!***KEYWORDS  AUTOMATIC INTEGRATOR, GAUSS-KRONROD RULES,
!             GENERAL-PURPOSE, GLOBALLY ADAPTIVE, INTEGRAND EXAMINATOR,
!             QUADPACK, QUADRATURE
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
!            EPSABS - Real
!                     Absolute accuracy requested
!            EPSREL - Real
!                     Relative accuracy requested
!                     If  EPSABS <= 0
!                     and EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28),
!                     the routine will end with IER = 6.
!
!            KEY    - Integer
!                     Key for choice of local integration rule
!                     A Gauss-Kronrod pair is used with
!                          7 - 15 points if KEY < 2,
!                         10 - 21 points if KEY = 2,
!                         15 - 31 points if KEY = 3,
!                         20 - 41 points if KEY = 4,
!                         25 - 51 points if KEY = 5,
!                         30 - 61 points if KEY > 5.
!
!            LIMIT  - Integer
!                     Gives an upper bound on the number of subintervals
!                     in the partition of (A,B), LIMIT >= 1.
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
!                             routine. It is assumed that the requested
!                             accuracy has been achieved.
!                     IER > 0 Abnormal termination of the routine
!                             The estimates for result and error are
!                             less reliable. It is assumed that the
!                             requested accuracy has not been achieved.
!            ERROR MESSAGES
!                     IER = 1 Maximum number of subdivisions allowed
!                             has been achieved. One can allow more
!                             subdivisions by increasing the value
!                             of LIMIT.
!                             However, if this yields no improvement it
!                             is rather advised to analyze the integrand
!                             in order to determine the integration
!                             difficulties. If the position of a local
!                             difficulty can be determined(e.g.
!                             SINGULARITY, DISCONTINUITY within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling the integrator on the
!                             subranges. If possible, an appropriate
!                             special-purpose integrator should be used
!                             which is designed for handling the type of
!                             difficulty involved.
!                         = 2 The occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 Extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 The input is invalid, because
!                             (EPSABS <= 0 and
!                              EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28),
!                             RESULT, ABSERR, NEVAL, LAST, RLIST(1) ,
!                             ELIST(1) and IORD(1) are set to zero.
!                             ALIST(1) and BLIST(1) are set to A and B
!                             respectively.
!
!            ALIST   - Real
!                      Vector of dimension at least LIMIT, the first
!                       LAST  elements of which are the left
!                      end points of the subintervals in the partition
!                      of the given integration range (A,B)
!
!            BLIST   - Real
!                      Vector of dimension at least LIMIT, the first
!                       LAST  elements of which are the right
!                      end points of the subintervals in the partition
!                      of the given integration range (A,B)
!
!            RLIST   - Real
!                      Vector of dimension at least LIMIT, the first
!                       LAST  elements of which are the
!                      integral approximations on the subintervals
!
!            ELIST   - Real
!                      Vector of dimension at least LIMIT, the first
!                       LAST  elements of which are the moduli of the
!                      absolute error estimates on the subintervals
!
!            IORD    - Integer
!                      Vector of dimension at least LIMIT, the first K
!                      elements of which are pointers to the
!                      error estimates over the subintervals,
!                      such that ELIST(IORD(1)), ...,
!                      ELIST(IORD(K)) form a decreasing sequence,
!                      with K = LAST if LAST <= (LIMIT/2+2), and
!                      K = LIMIT+1-LAST otherwise
!
!            LAST    - Integer
!                      Number of subintervals actually produced in the
!                      subdivision process
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  QK15, QK21, QK31, QK41, QK51, QK61, QPSRT, R1MACH
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  QAGE
!
  REAL A,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2,B,BLIST, &
    B1,B2,DEFABS,DEFAB1,DEFAB2,R1MACH,ELIST,EPMACH, &
    EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERROR2,ERRO12,ERRSUM,F, &
    RESABS,RESULT,RLIST,UFLOW
  INTEGER IER,IORD,IROFF1,IROFF2,K,KEY,KEYF,LAST, &
    LIMIT,MAXERR,NEVAL,NRMAX
!
  DIMENSION ALIST(*),BLIST(*),ELIST(*),IORD(*), &
    RLIST(*)
!
  EXTERNAL F
!
!            LIST OF MAJOR VARIABLES
!            -----------------------
!
!           ALIST     - LIST OF LEFT END POINTS OF ALL SUBINTERVALS
!                       CONSIDERED UP TO NOW
!           BLIST     - LIST OF RIGHT END POINTS OF ALL SUBINTERVALS
!                       CONSIDERED UP TO NOW
!           RLIST(I)  - APPROXIMATION TO THE INTEGRAL OVER
!                      (ALIST(I),BLIST(I))
!           ELIST(I)  - ERROR ESTIMATE APPLYING TO RLIST(I)
!           MAXERR    - POINTER TO THE INTERVAL WITH LARGEST
!                       ERROR ESTIMATE
!           ERRMAX    - ELIST(MAXERR)
!           AREA      - SUM OF THE INTEGRALS OVER THE SUBINTERVALS
!           ERRSUM    - SUM OF THE ERRORS OVER THE SUBINTERVALS
!           ERRBND    - REQUESTED ACCURACY MAX(EPSABS,EPSREL*
!                       ABS(RESULT))
!           *****1    - VARIABLE FOR THE LEFT SUBINTERVAL
!           *****2    - VARIABLE FOR THE RIGHT SUBINTERVAL
!           LAST      - INDEX FOR SUBDIVISION
!
!
!           MACHINE DEPENDENT CONSTANTS
!           ---------------------------
!
!           EPMACH  IS THE LARGEST RELATIVE SPACING.
!           UFLOW  IS THE SMALLEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  QAGE
  EPMACH = R1MACH(4)
  UFLOW = R1MACH(1)
!
!           TEST ON VALIDITY OF PARAMETERS
!           ------------------------------
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
  if ( EPSABS <= 0.0E+00.AND. &
    EPSREL < MAX(0.5E+02*EPMACH,0.5E-14)) IER = 6
  if ( IER == 6) go to 999
!
!           FIRST APPROXIMATION TO THE INTEGRAL
!           -----------------------------------
!
  KEYF = KEY
  if ( KEY <= 0) KEYF = 1
  if ( KEY >= 7) KEYF = 6
  NEVAL = 0
  if ( KEYF == 1) call QK15(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  if ( KEYF == 2) call QK21(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  if ( KEYF == 3) call QK31(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  if ( KEYF == 4) call QK41(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  if ( KEYF == 5) call QK51(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  if ( KEYF == 6) call QK61(F,A,B,RESULT,ABSERR,DEFABS,RESABS)
  LAST = 1
  RLIST(1) = RESULT
  ELIST(1) = ABSERR
  IORD(1) = 1
!
!           TEST ON ACCURACY.
!
  ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
  if ( ABSERR <= 0.5E+02*EPMACH*DEFABS.AND.ABSERR >  &
    ERRBND) IER = 2
  if ( LIMIT == 1) IER = 1
  if ( IER /= 0.OR.(ABSERR <= ERRBND.AND.ABSERR /= RESABS) &
    .OR.ABSERR == 0.0E+00) go to 60
!
!           INITIALIZATION
!           --------------
!
!
  ERRMAX = ABSERR
  MAXERR = 1
  AREA = RESULT
  ERRSUM = ABSERR
  NRMAX = 1
  IROFF1 = 0
  IROFF2 = 0
!
!           MAIN DO-LOOP
!           ------------
!
  DO 30 LAST = 2,LIMIT
!
!           BISECT THE SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE.
!
    A1 = ALIST(MAXERR)
    B1 = 0.5E+00*(ALIST(MAXERR)+BLIST(MAXERR))
    A2 = B1
    B2 = BLIST(MAXERR)
    if ( KEYF == 1) call QK15(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    if ( KEYF == 2) call QK21(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    if ( KEYF == 3) call QK31(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    if ( KEYF == 4) call QK41(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    if ( KEYF == 5) call QK51(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    if ( KEYF == 6) call QK61(F,A1,B1,AREA1,ERROR1,RESABS,DEFAB1)
    if ( KEYF == 1) call QK15(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    if ( KEYF == 2) call QK21(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    if ( KEYF == 3) call QK31(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    if ( KEYF == 4) call QK41(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    if ( KEYF == 5) call QK51(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
    if ( KEYF == 6) call QK61(F,A2,B2,AREA2,ERROR2,RESABS,DEFAB2)
!
!           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
!           AND ERROR AND TEST FOR ACCURACY.
!
    NEVAL = NEVAL+1
    AREA12 = AREA1+AREA2
    ERRO12 = ERROR1+ERROR2
    ERRSUM = ERRSUM+ERRO12-ERRMAX
    AREA = AREA+AREA12-RLIST(MAXERR)
    if ( DEFAB1 == ERROR1.OR.DEFAB2 == ERROR2) go to 5
    if ( ABS(RLIST(MAXERR)-AREA12) <= 0.1E-04*ABS(AREA12) &
    .AND.ERRO12 >= 0.99E+00*ERRMAX) IROFF1 = IROFF1+1
    if ( LAST > 10.AND.ERRO12 > ERRMAX) IROFF2 = IROFF2+1
    5   RLIST(MAXERR) = AREA1
    RLIST(LAST) = AREA2
    ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
    if ( ERRSUM <= ERRBND) go to 8
!
!           TEST FOR ROUNDOFF ERROR AND EVENTUALLY
!           SET ERROR FLAG.
!
    if ( IROFF1 >= 6.OR.IROFF2 >= 20) IER = 2
!
!           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF
!           SUBINTERVALS EQUALS LIMIT.
!
    if ( LAST == LIMIT) IER = 1
!
!           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
!           AT A POINT OF THE INTEGRATION RANGE.
!
    if ( MAX(ABS(A1),ABS(B2)) <= (0.1E+01+0.1E+03* &
    EPMACH)*(ABS(A2)+0.1E+04*UFLOW)) IER = 3
!
!           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
!
    8   if ( ERROR2 > ERROR1) go to 10
    ALIST(LAST) = A2
    BLIST(MAXERR) = B1
    BLIST(LAST) = B2
    ELIST(MAXERR) = ERROR1
    ELIST(LAST) = ERROR2
    go to 20
   10   ALIST(MAXERR) = A2
    ALIST(LAST) = A1
    BLIST(LAST) = B1
    RLIST(MAXERR) = AREA2
    RLIST(LAST) = AREA1
    ELIST(MAXERR) = ERROR2
    ELIST(LAST) = ERROR1
!
!           call SUBROUTINE QPSRT TO MAINTAIN THE DESCENDING ORDERING
!           IN THE LIST OF ERROR ESTIMATES AND SELECT THE
!           SUBINTERVAL WITH THE LARGEST ERROR ESTIMATE (TO BE
!           BISECTED NEXT).
!
   20   call QPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
! ***JUMP OUT OF DO-LOOP
    if ( IER /= 0.OR.ERRSUM <= ERRBND) go to 40
   30 CONTINUE
!
!           COMPUTE FINAL RESULT.
!           ---------------------
!
   40 RESULT = 0.0E+00
  DO 50 K=1,LAST
    RESULT = RESULT+RLIST(K)
   50 CONTINUE
  ABSERR = ERRSUM
   60 if ( KEYF /= 1) NEVAL = (10*KEYF+1)*(2*NEVAL+1)
  if ( KEYF == 1) NEVAL = 30*NEVAL+15
  999 RETURN
end
