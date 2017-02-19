subroutine DQAWSE (F, A, B, ALFA, BETA, INTEGR, EPSABS, EPSREL, &
     LIMIT, RESULT, ABSERR, NEVAL, IER, ALIST, BLIST, RLIST, ELIST, &
     IORD, LAST)
!
!! DQAWSE approimates the integral F(X)*W(X); W(X) has endpoint singularities.
!
!***PURPOSE  The routine calculates an approximation result to a given
!            definite integral I = Integral of F*W over (A,B),
!            (where W shows a singular behaviour at the end points,
!            see parameter INTEGR).
!            Hopefully satisfying following claim for accuracy
!            ABS(I-RESULT) <= MAX(EPSABS,EPSREL*ABS(I)).
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A2A1
!***TYPE      DOUBLE PRECISION (QAWSE-S, DQAWSE-D)
!***KEYWORDS  ALGEBRAIC-LOGARITHMIC END POINT SINGULARITIES,
!             AUTOMATIC INTEGRATOR, CLENSHAW-CURTIS METHOD, QUADPACK,
!             QUADRATURE, SPECIAL-PURPOSE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Integration of functions having algebraico-logarithmic
!        end point singularities
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
!                     Upper limit of integration, B > A
!                     If B <= A, the routine will end with IER = 6.
!
!            ALFA   - Double precision
!                     Parameter in the WEIGHT function, ALFA > (-1)
!                     If ALFA <= (-1), the routine will end with
!                     IER = 6.
!
!            BETA   - Double precision
!                     Parameter in the WEIGHT function, BETA > (-1)
!                     If BETA <= (-1), the routine will end with
!                     IER = 6.
!
!            INTEGR - Integer
!                     Indicates which WEIGHT function is to be used
!                     = 1  (X-A)**ALFA*(B-X)**BETA
!                     = 2  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)
!                     = 3  (X-A)**ALFA*(B-X)**BETA*LOG(B-X)
!                     = 4  (X-A)**ALFA*(B-X)**BETA*LOG(X-A)*LOG(B-X)
!                     If INTEGR < 1 or INTEGR > 4, the routine
!                     will end with IER = 6.
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
!                     in the partition of (A,B), LIMIT >= 2
!                     If LIMIT < 2, the routine will end with IER = 6.
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
!                             the estimates for the integral and error
!                             are less reliable. It is assumed that the
!                             requested accuracy has not been achieved.
!            ERROR MESSAGES
!                         = 1 Maximum number of subdivisions allowed
!                             has been achieved. One can allow more
!                             subdivisions by increasing the value of
!                             LIMIT. However, if this yields no
!                             improvement, it is advised to analyze the
!                             integrand in order to determine the
!                             integration difficulties which prevent the
!                             requested tolerance from being achieved.
!                             In case of a jump DISCONTINUITY or a local
!                             SINGULARITY of algebraico-logarithmic type
!                             at one or more interior points of the
!                             integration range, one should proceed by
!                             splitting up the interval at these
!                             points and calling the integrator on the
!                             subranges.
!                         = 2 The occurrence of roundoff error is
!                             detected, which prevents the requested
!                             tolerance from being achieved.
!                         = 3 Extremely bad integrand behaviour occurs
!                             at some points of the integration
!                             interval.
!                         = 6 The input is invalid, because
!                             B <= A or ALFA <= (-1) or BETA <= (-1), or
!                             INTEGR < 1 or INTEGR > 4, or
!                             (EPSABS <= 0 and
!                              EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28),
!                             or LIMIT < 2.
!                             RESULT, ABSERR, NEVAL, RLIST(1), ELIST(1),
!                             IORD(1) and LAST are set to zero. ALIST(1)
!                             and BLIST(1) are set to A and B
!                             respectively.
!
!            ALIST  - Double precision
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the left
!                     end points of the subintervals in the partition
!                     of the given integration range (A,B)
!
!            BLIST  - Double precision
!                     Vector of dimension at least LIMIT, the first
!                      LAST  elements of which are the right
!                     end points of the subintervals in the partition
!                     of the given integration range (A,B)
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
!            IORD   - Integer
!                     Vector of dimension at least LIMIT, the first K
!                     of which are pointers to the error
!                     estimates over the subintervals, so that
!                     ELIST(IORD(1)), ..., ELIST(IORD(K)) with K = LAST
!                     If LAST <= (LIMIT/2+2), and K = LIMIT+1-LAST
!                     otherwise form a decreasing sequence
!
!            LAST   - Integer
!                     Number of subintervals actually produced in
!                     the subdivision process
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DQC25S, DQMOMO, DQPSRT
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DQAWSE
!
  DOUBLE PRECISION A,ABSERR,ALFA,ALIST,AREA,AREA1,AREA12,AREA2,A1, &
    A2,B,BETA,BLIST,B1,B2,CENTRE,D1MACH,ELIST,EPMACH, &
    EPSABS,EPSREL,ERRBND,ERRMAX,ERROR1,ERRO12,ERROR2,ERRSUM,F, &
    RESAS1,RESAS2,RESULT,RG,RH,RI,RJ,RLIST,UFLOW
  INTEGER IER,INTEGR,IORD,IROFF1,IROFF2,K,LAST,LIMIT,MAXERR,NEV, &
    NEVAL,NRMAX
!
  EXTERNAL F
!
  DIMENSION ALIST(*),BLIST(*),RLIST(*),ELIST(*), &
    IORD(*),RI(25),RJ(25),RH(25),RG(25)
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
!***FIRST EXECUTABLE STATEMENT  DQAWSE
  EPMACH = D1MACH(4)
  UFLOW = D1MACH(1)
!
!           TEST ON VALIDITY OF PARAMETERS
!           ------------------------------
!
  IER = 6
  NEVAL = 0
  LAST = 0
  RLIST(1) = 0.0D+00
  ELIST(1) = 0.0D+00
  IORD(1) = 0
  RESULT = 0.0D+00
  ABSERR = 0.0D+00
  if (B <= A.OR.(EPSABS == 0.0D+00.AND. &
      EPSREL < MAX(0.5D+02*EPMACH,0.5D-28)).OR.ALFA <= (-0.1D+01) &
      .OR.BETA <= (-0.1D+01).OR.INTEGR < 1.OR.INTEGR > 4.OR. &
      LIMIT < 2) go to 999
  IER = 0
!
!           COMPUTE THE MODIFIED CHEBYSHEV MOMENTS.
!
  call DQMOMO(ALFA,BETA,RI,RJ,RG,RH,INTEGR)
!
!           INTEGRATE OVER THE INTERVALS (A,(A+B)/2) AND ((A+B)/2,B).
!
  CENTRE = 0.5D+00*(B+A)
  call DQC25S(F,A,B,A,CENTRE,ALFA,BETA,RI,RJ,RG,RH,AREA1, &
    ERROR1,RESAS1,INTEGR,NEV)
  NEVAL = NEV
  call DQC25S(F,A,B,CENTRE,B,ALFA,BETA,RI,RJ,RG,RH,AREA2, &
    ERROR2,RESAS2,INTEGR,NEV)
  LAST = 2
  NEVAL = NEVAL+NEV
  RESULT = AREA1+AREA2
  ABSERR = ERROR1+ERROR2
!
!           TEST ON ACCURACY.
!
  ERRBND = MAX(EPSABS,EPSREL*ABS(RESULT))
!
!           INITIALIZATION
!           --------------
!
  if ( ERROR2 > ERROR1) go to 10
  ALIST(1) = A
  ALIST(2) = CENTRE
  BLIST(1) = CENTRE
  BLIST(2) = B
  RLIST(1) = AREA1
  RLIST(2) = AREA2
  ELIST(1) = ERROR1
  ELIST(2) = ERROR2
  go to 20
   10 ALIST(1) = CENTRE
  ALIST(2) = A
  BLIST(1) = B
  BLIST(2) = CENTRE
  RLIST(1) = AREA2
  RLIST(2) = AREA1
  ELIST(1) = ERROR2
  ELIST(2) = ERROR1
   20 IORD(1) = 1
  IORD(2) = 2
  if ( LIMIT == 2) IER = 1
  if ( ABSERR <= ERRBND.OR.IER == 1) go to 999
  ERRMAX = ELIST(1)
  MAXERR = 1
  NRMAX = 1
  AREA = RESULT
  ERRSUM = ABSERR
  IROFF1 = 0
  IROFF2 = 0
!
!            MAIN DO-LOOP
!            ------------
!
  DO 60 LAST = 3,LIMIT
!
!           BISECT THE SUBINTERVAL WITH LARGEST ERROR ESTIMATE.
!
    A1 = ALIST(MAXERR)
    B1 = 0.5D+00*(ALIST(MAXERR)+BLIST(MAXERR))
    A2 = B1
    B2 = BLIST(MAXERR)
!
    call DQC25S(F,A,B,A1,B1,ALFA,BETA,RI,RJ,RG,RH,AREA1, &
    ERROR1,RESAS1,INTEGR,NEV)
    NEVAL = NEVAL+NEV
    call DQC25S(F,A,B,A2,B2,ALFA,BETA,RI,RJ,RG,RH,AREA2, &
    ERROR2,RESAS2,INTEGR,NEV)
    NEVAL = NEVAL+NEV
!
!           IMPROVE PREVIOUS APPROXIMATIONS INTEGRAL AND ERROR
!           AND TEST FOR ACCURACY.
!
    AREA12 = AREA1+AREA2
    ERRO12 = ERROR1+ERROR2
    ERRSUM = ERRSUM+ERRO12-ERRMAX
    AREA = AREA+AREA12-RLIST(MAXERR)
    if ( A == A1.OR.B == B2) go to 30
    if ( RESAS1 == ERROR1.OR.RESAS2 == ERROR2) go to 30
!
!           TEST FOR ROUNDOFF ERROR.
!
    if ( ABS(RLIST(MAXERR)-AREA12) < 0.1D-04*ABS(AREA12) &
    .AND.ERRO12 >= 0.99D+00*ERRMAX) IROFF1 = IROFF1+1
    if ( LAST > 10.AND.ERRO12 > ERRMAX) IROFF2 = IROFF2+1
   30   RLIST(MAXERR) = AREA1
    RLIST(LAST) = AREA2
!
!           TEST ON ACCURACY.
!
    ERRBND = MAX(EPSABS,EPSREL*ABS(AREA))
    if ( ERRSUM <= ERRBND) go to 35
!
!           SET ERROR FLAG IN THE CASE THAT THE NUMBER OF INTERVAL
!           BISECTIONS EXCEEDS LIMIT.
!
    if ( LAST == LIMIT) IER = 1
!
!
!           SET ERROR FLAG IN THE CASE OF ROUNDOFF ERROR.
!
    if ( IROFF1 >= 6.OR.IROFF2 >= 20) IER = 2
!
!           SET ERROR FLAG IN THE CASE OF BAD INTEGRAND BEHAVIOUR
!           AT INTERIOR POINTS OF INTEGRATION RANGE.
!
    if ( MAX(ABS(A1),ABS(B2)) <= (0.1D+01+0.1D+03*EPMACH)* &
    (ABS(A2)+0.1D+04*UFLOW)) IER = 3
!
!           APPEND THE NEWLY-CREATED INTERVALS TO THE LIST.
!
   35   if ( ERROR2 > ERROR1) go to 40
    ALIST(LAST) = A2
    BLIST(MAXERR) = B1
    BLIST(LAST) = B2
    ELIST(MAXERR) = ERROR1
    ELIST(LAST) = ERROR2
    go to 50
   40   ALIST(MAXERR) = A2
    ALIST(LAST) = A1
    BLIST(LAST) = B1
    RLIST(MAXERR) = AREA2
    RLIST(LAST) = AREA1
    ELIST(MAXERR) = ERROR2
    ELIST(LAST) = ERROR1
!
!           call SUBROUTINE DQPSRT TO MAINTAIN THE DESCENDING ORDERING
!           IN THE LIST OF ERROR ESTIMATES AND SELECT THE SUBINTERVAL
!           WITH LARGEST ERROR ESTIMATE (TO BE BISECTED NEXT).
!
   50   call DQPSRT(LIMIT,LAST,MAXERR,ERRMAX,ELIST,IORD,NRMAX)
! ***JUMP OUT OF DO-LOOP
    if (IER /= 0.OR.ERRSUM <= ERRBND) go to 70
   60 CONTINUE
!
!           COMPUTE FINAL RESULT.
!           ---------------------
!
   70 RESULT = 0.0D+00
  DO 80 K=1,LAST
    RESULT = RESULT+RLIST(K)
   80 CONTINUE
  ABSERR = ERRSUM
  999 RETURN
end
