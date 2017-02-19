subroutine DQAWCE (F, A, B, C, EPSABS, EPSREL, LIMIT, RESULT, &
     ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, IORD, LAST)
!
!! DQAWCE approximates the Cauchy principal value of the integral of F(X)/(X-C).
!
!  The routine calculates an approximation result to a
!            CAUCHY PRINCIPAL VALUE I = Integral of F*W over (A,B)
!            (W(X) = 1/(X-C), (C /= A, C /= B), hopefully satisfying
!            following claim for accuracy
!            ABS(I-RESULT) <= MAX(EPSABS,EPSREL*ABS(I))
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A2A1, J4
!***TYPE      DOUBLE PRECISION (QAWCE-S, DQAWCE-D)
!***KEYWORDS  AUTOMATIC INTEGRATOR, CAUCHY PRINCIPAL VALUE,
!             CLENSHAW-CURTIS METHOD, QUADPACK, QUADRATURE,
!             SPECIAL-PURPOSE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Computation of a CAUCHY PRINCIPAL VALUE
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
!            C      - Double precision
!                     Parameter in the WEIGHT function, C /= A, C /= B
!                     If C = A OR C = B, the routine will end with
!                     IER = 6.
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
!                     in the partition of (A,B), LIMIT >= 1
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
!                     IER > 0 Abnormal termination of the routine
!                             the estimates for integral and error are
!                             less reliable. It is assumed that the
!                             requested accuracy has not been achieved.
!            ERROR MESSAGES
!                     IER = 1 Maximum number of subdivisions allowed
!                             has been achieved. One can allow more sub-
!                             divisions by increasing the value of
!                             LIMIT. However, if this yields no
!                             improvement it is advised to analyze the
!                             the integrand, in order to determine the
!                             the integration difficulties. If the
!                             position of a local difficulty can be
!                             determined (e.g. SINGULARITY,
!                             DISCONTINUITY within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling
!                             appropriate integrators on the subranges.
!                         = 2 The occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 Extremely bad integrand behaviour
!                             occurs at some interior points of
!                             the integration interval.
!                         = 6 The input is invalid, because
!                             C = A or C = B or
!                             (EPSABS <= 0 and
!                              EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28))
!                             or LIMIT < 1.
!                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
!                             IORD(1) and LAST are set to zero. ALIST(1)
!                             and BLIST(1) are set to A and B
!                             respectively.
!
!            ALIST   - Double precision
!                      Vector of dimension at least LIMIT, the first
!                       LAST  elements of which are the left
!                      end points of the subintervals in the partition
!                      of the given integration range (A,B)
!
!            BLIST   - Double precision
!                      Vector of dimension at least LIMIT, the first
!                       LAST  elements of which are the right
!                      end points of the subintervals in the partition
!                      of the given integration range (A,B)
!
!            RLIST   - Double precision
!                      Vector of dimension at least LIMIT, the first
!                       LAST  elements of which are the integral
!                      approximations on the subintervals
!
!            ELIST   - Double precision
!                      Vector of dimension LIMIT, the first  LAST
!                      elements of which are the moduli of the absolute
!                      error estimates on the subintervals
!
!            IORD    - Integer
!                      Vector of dimension at least LIMIT, the first K
!                      elements of which are pointers to the error
!                      estimates over the subintervals, so that
!                      ELIST(IORD(1)), ..., ELIST(IORD(K)) with K = LAST
!                      If LAST <= (LIMIT/2+2), and K = LIMIT+1-LAST
!                      otherwise, form a decreasing sequence
!
!            LAST    - Integer
!                      Number of subintervals actually produced in
!                      the subdivision process
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DQC25C, DQPSRT
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DQAWCE
!
  DOUBLE PRECISION A,AA,ABSERR,ALIST,AREA,AREA1,AREA12,AREA2,A1,A2, &
    B,BB,BLIST,B1,B2,C,D1MACH,ELIST,EPMACH,EPSABS,EPSREL, &
    ERRBND,ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,F,RESULT,RLIST,UFLOW
  INTEGER IER,IORD,IROFF1,IROFF2,K,KRULE,LAST,LIMIT,MAXERR,NEV, &
    NEVAL,NRMAX
!
  DIMENSION ALIST(*),BLIST(*),RLIST(*),ELIST(*), &
    IORD(*)
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
!                       (ALIST(I),BLIST(I))
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
!            MACHINE DEPENDENT CONSTANTS
!            ---------------------------
!
!           EPMACH IS THE LARGEST RELATIVE SPACING.
!           UFLOW IS THE SMALLEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  DQAWCE
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
!
!
!           TEST ON VALIDITY OF PARAMETERS
!           ------------------------------
!
  IER = 6
  NEVAL = 0
  LAST = 0
  ALIST(1) = A
  BLIST(1) = B
  RLIST(1) = 0.0D+00
  ELIST(1) = 0.0D+00
  IORD(1) = 0
  RESULT = 0.0D+00
  ABSERR = 0.0D+00
  if (C == A.OR.C == B.OR.(EPSABS <= 0.0D+00.AND. &
      EPSREL < MAX(0.5D+02*EPMACH,0.5D-28))) go to 999
!
!           FIRST APPROXIMATION TO THE INTEGRAL
!           -----------------------------------
!
  AA=A
  BB=B
  if (A <= B) go to 10
  AA=B
  BB=A
10    IER=0
  KRULE = 1
  call DQC25C(F,AA,BB,C,RESULT,ABSERR,KRULE,NEVAL)
  LAST = 1
  RLIST(1) = RESULT
  ELIST(1) = ABSERR
  IORD(1) = 1
  ALIST(1) = A
  BLIST(1) = B
!
!           TEST ON ACCURACY
!
  ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
  if ( LIMIT == 1) IER = 1
  if ( ABSERR < MIN(0.1D-01*ABS(RESULT),ERRBND) &
    .OR.IER == 1) go to 70
!
!           INITIALIZATION
!           --------------
!
  ALIST(1) = AA
  BLIST(1) = BB
  RLIST(1) = RESULT
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
  DO 40 LAST = 2,LIMIT
!
!           BISECT THE SUBINTERVAL WITH NRMAX-TH LARGEST
!           ERROR ESTIMATE.
!
    A1 = ALIST(MAXERR)
    B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
    B2 = BLIST(MAXERR)
    if ( C <= B1.AND.C > A1) B1 = 0.5D+00*(C+B2)
    if ( C > B1.AND.C < B2) B1 = 0.5D+00*(A1+C)
    A2 = B1
    KRULE = 2
    call DQC25C(F,A1,B1,C,AREA1,ERROR1,KRULE,NEV)
    NEVAL = NEVAL+NEV
    call DQC25C(F,A2,B2,C,AREA2,ERROR2,KRULE,NEV)
    NEVAL = NEVAL+NEV
!
!           IMPROVE PREVIOUS APPROXIMATIONS TO INTEGRAL
!           AND ERROR AND TEST FOR ACCURACY.
!
    AREA12 = AREA1+AREA2
    ERRO12 = ERROR1+ERROR2
    ERRSUM = ERRSUM+ERRO12-ERRMAX
    AREA = AREA+AREA12-RLIST(MAXERR)
    if ( ABS(RLIST(MAXERR)-AREA12) < 0.1D-04*ABS(AREA12) &
      .AND.ERRO12 >= 0.99D+00*ERRMAX.AND.KRULE == 0) &
      IROFF1 = IROFF1+1
    if ( LAST > 10.AND.ERRO12 > ERRMAX.AND.KRULE == 0) &
      IROFF2 = IROFF2+1
    RLIST(MAXERR) = AREA1
    RLIST(LAST) = AREA2
    ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
    if ( ERRSUM <= ERRBND) go to 15
!
!           TEST FOR ROUNDOFF ERROR AND EVENTUALLY SET ERROR FLAG.
!
    if ( IROFF1 >= 6.AND.IROFF2 > 20) IER = 2
!
!           SET ERROR FLAG IN THE CASE THAT NUMBER OF INTERVAL
!           BISECTIONS EXCEEDS LIMIT.
!
    if ( LAST == LIMIT) IER = 1
!
!           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
!           AT A POINT OF THE INTEGRATION RANGE.
!
    if ( MAX(ABS(A1),ABS(B2)) <= (0.1D+01+0.1D+03*EPMACH) &
      *(ABS(A2)+0.1D+04*UFLOW)) IER = 3
!
!           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
!
   15   if ( ERROR2 > ERROR1) go to 20
    ALIST(LAST) = A2
    BLIST(MAXERR) = B1
    BLIST(LAST) = B2
    ELIST(MAXERR) = ERROR1
    ELIST(LAST) = ERROR2
    go to 30
   20   ALIST(MAXERR) = A2
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
   30    call DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
! ***JUMP OUT OF DO-LOOP
    if ( IER /= 0.OR.ERRSUM <= ERRBND) go to 50
   40 CONTINUE
!
!           COMPUTE FINAL RESULT.
!           ---------------------
!
   50 RESULT = 0.0D+00
  DO 60 K=1,LAST
    RESULT = RESULT+RLIST(K)
   60 CONTINUE
  ABSERR = ERRSUM
   70 if (AA == B) RESULT=-RESULT
  999 RETURN
end
