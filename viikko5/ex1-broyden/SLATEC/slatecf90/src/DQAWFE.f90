subroutine DQAWFE (F, A, OMEGA, INTEGR, EPSABS, LIMLST, LIMIT, &
     MAXP1, RESULT, ABSERR, NEVAL, IER, RSLST, ERLST, IERLST, LST, &
     ALIST, BLIST, RLIST, ELIST, IORD, NNLOG, CHEBMO)
!
!! DQAWFE approximates a given Fourier integral.
!
!  The routine calculates an approximation result to a
!            given Fourier integral
!            I = Integral of F(X)*W(X) over (A,INFINITY)
!            where W(X)=COS(OMEGA*X) or W(X)=SIN(OMEGA*X),
!            hopefully satisfying following claim for accuracy
!            ABS(I-RESULT) <= EPSABS.
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A3A1
!***TYPE      DOUBLE PRECISION (QAWFE-S, DQAWFE-D)
!***KEYWORDS  AUTOMATIC INTEGRATOR, CONVERGENCE ACCELERATION,
!             FOURIER INTEGRALS, INTEGRATION BETWEEN ZEROS, QUADPACK,
!             QUADRATURE, SPECIAL-PURPOSE INTEGRAL
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Computation of Fourier integrals
!        Standard fortran subroutine
!        Double precision version
!
!        PARAMETERS
!         ON ENTRY
!            F      - Double precision
!                     Function subprogram defining the integrand
!                     Function F(X). The actual name for F needs to
!                     be declared E X T E R N A L in the driver program.
!
!            A      - Double precision
!                     Lower limit of integration
!
!            OMEGA  - Double precision
!                     Parameter in the WEIGHT function
!
!            INTEGR - Integer
!                     Indicates which WEIGHT function is used
!                     INTEGR = 1      W(X) = COS(OMEGA*X)
!                     INTEGR = 2      W(X) = SIN(OMEGA*X)
!                     If INTEGR /= 1.AND.INTEGR /= 2, the routine will
!                     end with IER = 6.
!
!            EPSABS - Double precision
!                     absolute accuracy requested, EPSABS > 0
!                     If EPSABS <= 0, the routine will end with IER = 6.
!
!            LIMLST - Integer
!                     LIMLST gives an upper bound on the number of
!                     cycles, LIMLST >= 1.
!                     If LIMLST < 3, the routine will end with IER = 6.
!
!            LIMIT  - Integer
!                     Gives an upper bound on the number of subintervals
!                     allowed in the partition of each cycle, LIMIT >= 1
!                     each cycle, LIMIT >= 1.
!
!            MAXP1  - Integer
!                     Gives an upper bound on the number of
!                     Chebyshev moments which can be stored, I.E.
!                     for the intervals of lengths ABS(B-A)*2**(-L),
!                     L=0,1, ..., MAXP1-2, MAXP1 >= 1
!
!         ON RETURN
!            RESULT - Double precision
!                     Approximation to the integral X
!
!            ABSERR - Double precision
!                     Estimate of the modulus of the absolute error,
!                     which should equal or exceed ABS(I-RESULT)
!
!            NEVAL  - Integer
!                     Number of integrand evaluations
!
!            IER    - IER = 0 Normal and reliable termination of
!                             the routine. It is assumed that the
!                             requested accuracy has been achieved.
!                     IER > 0 Abnormal termination of the routine. The
!                             estimates for integral and error are less
!                             reliable. It is assumed that the requested
!                             accuracy has not been achieved.
!            ERROR MESSAGES
!                    If OMEGA /= 0
!                     IER = 1 Maximum number of  cycles  allowed
!                             Has been achieved., i.e. of subintervals
!                             (A+(K-1)C,A+KC) where
!                             C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
!                             for K = 1, 2, ..., LST.
!                             One can allow more cycles by increasing
!                             the value of LIMLST (and taking the
!                             according dimension adjustments into
!                             account).
!                             Examine the array IWORK which contains
!                             the error flags on the cycles, in order to
!                             look for eventual local integration
!                             difficulties. If the position of a local
!                             difficulty can be determined (e.g.
!                             SINGULARITY, DISCONTINUITY within the
!                             interval) one will probably gain from
!                             splitting up the interval at this point
!                             and calling appropriate integrators on
!                             the subranges.
!                         = 4 The extrapolation table constructed for
!                             convergence acceleration of the series
!                             formed by the integral contributions over
!                             the cycles, does not converge to within
!                             the requested accuracy. As in the case of
!                             IER = 1, it is advised to examine the
!                             array IWORK which contains the error
!                             flags on the cycles.
!                         = 6 The input is invalid because
!                             (INTEGR /= 1 AND INTEGR /= 2) or
!                              EPSABS <= 0 or LIMLST < 3.
!                              RESULT, ABSERR, NEVAL, LST are set
!                              to zero.
!                         = 7 Bad integrand behaviour occurs within one
!                             or more of the cycles. Location and type
!                             of the difficulty involved can be
!                             determined from the vector IERLST. Here
!                             LST is the number of cycles actually
!                             needed (see below).
!                             IERLST(K) = 1 The maximum number of
!                                           subdivisions (= LIMIT) has
!                                           been achieved on the K th
!                                           cycle.
!                                       = 2 Occurrence of roundoff error
!                                           is detected and prevents the
!                                           tolerance imposed on the
!                                           K th cycle, from being
!                                           achieved.
!                                       = 3 Extremely bad integrand
!                                           behaviour occurs at some
!                                           points of the K th cycle.
!                                       = 4 The integration procedure
!                                           over the K th cycle does
!                                           not converge (to within the
!                                           required accuracy) due to
!                                           roundoff in the
!                                           extrapolation procedure
!                                           invoked on this cycle. It
!                                           is assumed that the result
!                                           on this interval is the
!                                           best which can be obtained.
!                                       = 5 The integral over the K th
!                                           cycle is probably divergent
!                                           or slowly convergent. It
!                                           must be noted that
!                                           divergence can occur with
!                                           any other value of
!                                           IERLST(K).
!                    If OMEGA = 0 and INTEGR = 1,
!                    The integral is calculated by means of DQAGIE
!                    and IER = IERLST(1) (with meaning as described
!                    for IERLST(K), K = 1).
!
!            RSLST  - Double precision
!                     Vector of dimension at least LIMLST
!                     RSLST(K) contains the integral contribution
!                     over the interval (A+(K-1)C,A+KC) where
!                     C = (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA),
!                     K = 1, 2, ..., LST.
!                     Note that, if OMEGA = 0, RSLST(1) contains
!                     the value of the integral over (A,INFINITY).
!
!            ERLST  - Double precision
!                     Vector of dimension at least LIMLST
!                     ERLST(K) contains the error estimate corresponding
!                     with RSLST(K).
!
!            IERLST - Integer
!                     Vector of dimension at least LIMLST
!                     IERLST(K) contains the error flag corresponding
!                     with RSLST(K). For the meaning of the local error
!                     flags see description of output parameter IER.
!
!            LST    - Integer
!                     Number of subintervals needed for the integration
!                     If OMEGA = 0 then LST is set to 1.
!
!            ALIST, BLIST, RLIST, ELIST - Double precision
!                     vector of dimension at least LIMIT,
!
!            IORD, NNLOG - Integer
!                     Vector of dimension at least LIMIT, providing
!                     space for the quantities needed in the subdivision
!                     process of each cycle
!
!            CHEBMO - Double precision
!                     Array of dimension at least (MAXP1,25), providing
!                     space for the Chebyshev moments needed within the
!                     cycles
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, DQAGIE, DQAWOE, DQELG
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced variable.  (WRB)
!   891009  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  DQAWFE
!
  DOUBLE PRECISION A,ABSEPS,ABSERR,ALIST,BLIST,CHEBMO,CORREC,CYCLE, &
    C1,C2,DL,DRL,D1MACH,ELIST,ERLST,EP,EPS,EPSA, &
    EPSABS,ERRSUM,F,FACT,OMEGA,P,PI,P1,PSUM,RESEPS,RESULT,RES3LA, &
    RLIST,RSLST,UFLOW
  INTEGER IER,IERLST,INTEGR,IORD,KTMIN,L,LAST,LST,LIMIT,LIMLST,LL, &
      MAXP1,MOMCOM,NEV,NEVAL,NNLOG,NRES,NUMRL2
!
  DIMENSION ALIST(*),BLIST(*),CHEBMO(MAXP1,25),ELIST(*), &
    ERLST(*),IERLST(*),IORD(*),NNLOG(*),PSUM(52), &
    RES3LA(3),RLIST(*),RSLST(*)
!
  EXTERNAL F
!
!
!            THE DIMENSION OF  PSUM  IS DETERMINED BY THE VALUE OF
!            LIMEXP IN SUBROUTINE DQELG (PSUM MUST BE OF DIMENSION
!            (LIMEXP+2) AT LEAST).
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           C1, C2    - END POINTS OF SUBINTERVAL (OF LENGTH CYCLE)
!           CYCLE     - (2*INT(ABS(OMEGA))+1)*PI/ABS(OMEGA)
!           PSUM      - VECTOR OF DIMENSION AT LEAST (LIMEXP+2)
!                       (SEE ROUTINE DQELG)
!                       PSUM CONTAINS THE PART OF THE EPSILON TABLE
!                       WHICH IS STILL NEEDED FOR FURTHER COMPUTATIONS.
!                       EACH ELEMENT OF PSUM IS A PARTIAL SUM OF THE
!                       SERIES WHICH SHOULD SUM TO THE VALUE OF THE
!                       INTEGRAL.
!           ERRSUM    - SUM OF ERROR ESTIMATES OVER THE SUBINTERVALS,
!                       CALCULATED CUMULATIVELY
!           EPSA      - ABSOLUTE TOLERANCE REQUESTED OVER CURRENT
!                       SUBINTERVAL
!           CHEBMO    - ARRAY CONTAINING THE MODIFIED CHEBYSHEV
!                       MOMENTS (SEE ALSO ROUTINE DQC25F)
!
  SAVE P, PI
  DATA P/0.9D+00/
  DATA PI / 3.14159265358979323846264338327950D0 /
!
!           TEST ON VALIDITY OF PARAMETERS
!           ------------------------------
!
!***FIRST EXECUTABLE STATEMENT  DQAWFE
  RESULT = 0.0D+00
  ABSERR = 0.0D+00
  NEVAL = 0
  LST = 0
  IER = 0
  if ( (INTEGR /= 1.AND.INTEGR /= 2).OR.EPSABS <= 0.0D+00.OR. &
    LIMLST < 3) IER = 6
  if ( IER == 6) go to 999
  if ( OMEGA /= 0.0D+00) go to 10
!
!           INTEGRATION BY DQAGIE if OMEGA IS ZERO
!           --------------------------------------
!
  if ( INTEGR == 1) call DQAGIE(F,A,1,EPSABS,0.0D+00,LIMIT, &
    RESULT,ABSERR,NEVAL,IER,ALIST,BLIST,RLIST,ELIST,IORD,LAST)
  RSLST(1) = RESULT
  ERLST(1) = ABSERR
  IERLST(1) = IER
  LST = 1
  go to 999
!
!           INITIALIZATIONS
!           ---------------
!
   10 L = ABS(OMEGA)
  DL = 2*L+1
  CYCLE = DL*PI/ABS(OMEGA)
  IER = 0
  KTMIN = 0
  NEVAL = 0
  NUMRL2 = 0
  NRES = 0
  C1 = A
  C2 = CYCLE+A
  P1 = 0.1D+01-P
  UFLOW = D1MACH(1)
  EPS = EPSABS
  if ( EPSABS > UFLOW/P1) EPS = EPSABS*P1
  EP = EPS
  FACT = 0.1D+01
  CORREC = 0.0D+00
  ABSERR = 0.0D+00
  ERRSUM = 0.0D+00
!
!           MAIN DO-LOOP
!           ------------
!
  DO 50 LST = 1,LIMLST
!
!           INTEGRATE OVER CURRENT SUBINTERVAL.
!
    EPSA = EPS*FACT
    call DQAWOE(F,C1,C2,OMEGA,INTEGR,EPSA,0.0D+00,LIMIT,LST,MAXP1, &
    RSLST(LST),ERLST(LST),NEV,IERLST(LST),LAST,ALIST,BLIST,RLIST, &
    ELIST,IORD,NNLOG,MOMCOM,CHEBMO)
    NEVAL = NEVAL+NEV
    FACT = FACT*P
    ERRSUM = ERRSUM+ERLST(LST)
    DRL = 0.5D+02*ABS(RSLST(LST))
!
!           TEST ON ACCURACY WITH PARTIAL SUM
!
    if ( (ERRSUM+DRL) <= EPSABS.AND.LST >= 6) go to 80
    CORREC = MAX(CORREC,ERLST(LST))
    if ( IERLST(LST) /= 0) EPS = MAX(EP,CORREC*P1)
    if ( IERLST(LST) /= 0) IER = 7
    if ( IER == 7.AND.(ERRSUM+DRL) <= CORREC*0.1D+02.AND. &
    LST > 5) go to 80
    NUMRL2 = NUMRL2+1
    if ( LST > 1) go to 20
    PSUM(1) = RSLST(1)
    go to 40
   20   PSUM(NUMRL2) = PSUM(LL)+RSLST(LST)
    if ( LST == 2) go to 40
!
!           TEST ON MAXIMUM NUMBER OF SUBINTERVALS
!
    if ( LST == LIMLST) IER = 1
!
!           PERFORM NEW EXTRAPOLATION
!
    call DQELG(NUMRL2,PSUM,RESEPS,ABSEPS,RES3LA,NRES)
!
!           TEST WHETHER EXTRAPOLATED RESULT IS INFLUENCED BY ROUNDOFF
!
    KTMIN = KTMIN+1
    if ( KTMIN >= 15.AND.ABSERR <= 0.1D-02*(ERRSUM+DRL)) IER = 4
    if ( ABSEPS > ABSERR.AND.LST /= 3) go to 30
    ABSERR = ABSEPS
    RESULT = RESEPS
    KTMIN = 0
!
!           if IER IS NOT 0, CHECK WHETHER DIRECT RESULT (PARTIAL SUM)
!           OR EXTRAPOLATED RESULT YIELDS THE BEST INTEGRAL
!           APPROXIMATION
!
    if ( (ABSERR+0.1D+02*CORREC) <= EPSABS.OR. &
    (ABSERR <= EPSABS.AND.0.1D+02*CORREC >= EPSABS)) go to 60
   30   if ( IER /= 0.AND.IER /= 7) go to 60
   40   LL = NUMRL2
    C1 = C2
    C2 = C2+CYCLE
   50 CONTINUE
!
!         SET FINAL RESULT AND ERROR ESTIMATE
!         -----------------------------------
!
   60 ABSERR = ABSERR+0.1D+02*CORREC
  if ( IER == 0) go to 999
  if ( RESULT /= 0.0D+00.AND.PSUM(NUMRL2) /= 0.0D+00) go to 70
  if ( ABSERR > ERRSUM) go to 80
  if ( PSUM(NUMRL2) == 0.0D+00) go to 999
   70 if ( ABSERR/ABS(RESULT) > (ERRSUM+DRL)/ABS(PSUM(NUMRL2))) &
    go to 80
  if ( IER >= 1.AND.IER /= 7) ABSERR = ABSERR+DRL
  go to 999
   80 RESULT = PSUM(NUMRL2)
  ABSERR = ERRSUM+DRL
  999 RETURN
end
