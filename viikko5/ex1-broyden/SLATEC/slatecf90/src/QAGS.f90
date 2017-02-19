subroutine QAGS (F, A, B, EPSABS, EPSREL, RESULT, ABSERR, NEVAL, &
     IER, LIMIT, LENW, LAST, IWORK, WORK)
!
!! QAGS approximates the definite integral of F(X) over (A,B).
!
!  Hopefully satisfying following claim for accuracy
!  ABS(I-RESULT) <= MAX(EPSABS,EPSREL*ABS(I)).
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A1A1
!***TYPE      SINGLE PRECISION (QAGS-S, DQAGS-D)
!***KEYWORDS  AUTOMATIC INTEGRATOR, END POINT SINGULARITIES,
!             EXTRAPOLATION, GENERAL-PURPOSE, GLOBALLY ADAPTIVE,
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
!
!        PARAMETERS
!         ON ENTRY
!            F      - Real
!                     Function subprogram defining the integrand
!                     Function F(X). The actual name for F needs to be
!                     Declared E X T E R N A L in the driver program.
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
!                     And EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28),
!                     The routine will end with IER = 6.
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
!                             The estimates for integral and error are
!                             less reliable. It is assumed that the
!                             requested accuracy has not been achieved.
!            ERROR MESSAGES
!                     IER = 1 Maximum number of subdivisions allowed
!                             has been achieved. One can allow more sub-
!                             divisions by increasing the value of LIMIT
!                             (and taking the according dimension
!                             adjustments into account. However, if
!                             this yields no improvement it is advised
!                             to analyze the integrand in order to
!                             determine the integration difficulties. If
!                             the position of a local difficulty can be
!                             determined (e.g. SINGULARITY,
!                             DISCONTINUITY within the interval) one
!                             will probably gain from splitting up the
!                             interval at this point and calling the
!                             integrator on the subranges. If possible,
!                             an appropriate special-purpose integrator
!                             should be used, which is designed for
!                             handling the type of difficulty involved.
!                         = 2 The occurrence of roundoff error is detec-
!                             ted, which prevents the requested
!                             tolerance from being achieved.
!                             The error may be under-estimated.
!                         = 3 Extremely bad integrand behaviour
!                             occurs at some points of the integration
!                             interval.
!                         = 4 The algorithm does not converge.
!                             Roundoff error is detected in the
!                             Extrapolation table. It is presumed that
!                             the requested tolerance cannot be
!                             achieved, and that the returned result is
!                             the best which can be obtained.
!                         = 5 The integral is probably divergent, or
!                             slowly convergent. It must be noted that
!                             divergence can occur with any other value
!                             of IER.
!                         = 6 The input is invalid, because
!                             (EPSABS <= 0 AND
!                              EPSREL < MAX(50*REL.MACH.ACC.,0.5D-28)
!                             OR LIMIT < 1 OR LENW < LIMIT*4.
!                             RESULT, ABSERR, NEVAL, LAST are set to
!                             zero.  Except when LIMIT or LENW is
!                             invalid, IWORK(1), WORK(LIMIT*2+1) and
!                             WORK(LIMIT*3+1) are set to zero, WORK(1)
!                             is set to A and WORK(LIMIT+1) TO B.
!
!         DIMENSIONING PARAMETERS
!            LIMIT - Integer
!                    Dimensioning parameter for IWORK
!                    LIMIT determines the maximum number of subintervals
!                    in the partition of the given integration interval
!                    (A,B), LIMIT >= 1.
!                    if LIMIT < 1, the routine will end with IER = 6.
!
!            LENW  - Integer
!                    Dimensioning parameter for WORK
!                    LENW must be at least LIMIT*4.
!                    If LENW < LIMIT*4, the routine will end
!                    with IER = 6.
!
!            LAST  - Integer
!                    On return, LAST equals the number of subintervals
!                    produced in the subdivision process, determines the
!                    number of significant elements actually in the WORK
!                    Arrays.
!
!         WORK ARRAYS
!            IWORK - Integer
!                    Vector of dimension at least LIMIT, the first K
!                    elements of which contain pointers
!                    to the error estimates over the subintervals
!                    such that WORK(LIMIT*3+IWORK(1)),... ,
!                    WORK(LIMIT*3+IWORK(K)) form a decreasing
!                    sequence, with K = LAST if LAST <= (LIMIT/2+2),
!                    and K = LIMIT+1-LAST otherwise
!
!            WORK  - Real
!                    Vector of dimension at least LENW
!                    on return
!                    WORK(1), ..., WORK(LAST) contain the left
!                     end-points of the subintervals in the
!                     partition of (A,B),
!                    WORK(LIMIT+1), ..., WORK(LIMIT+LAST) contain
!                     the right end-points,
!                    WORK(LIMIT*2+1), ..., WORK(LIMIT*2+LAST) contain
!                     the integral approximations over the subintervals,
!                    WORK(LIMIT*3+1), ..., WORK(LIMIT*3+LAST)
!                     contain the error estimates.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  QAGSE, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!***END PROLOGUE  QAGS
!
!
  REAL A,ABSERR,B,EPSABS,EPSREL,F,RESULT,WORK
  INTEGER IER,IWORK,LENW,LIMIT,LVL,L1,L2,L3,NEVAL
!
  DIMENSION IWORK(*),WORK(*)
!
  EXTERNAL F
!
!         CHECK VALIDITY OF LIMIT AND LENW.
!
!***FIRST EXECUTABLE STATEMENT  QAGS
  IER = 6
  NEVAL = 0
  LAST = 0
  RESULT = 0.0E+00
  ABSERR = 0.0E+00
  if ( LIMIT < 1.OR.LENW < LIMIT*4) go to 10
!
!         PREPARE call FOR QAGSE.
!
  L1 = LIMIT+1
  L2 = LIMIT+L1
  L3 = LIMIT+L2
!
  call QAGSE(F,A,B,EPSABS,EPSREL,LIMIT,RESULT,ABSERR,NEVAL, &
    IER,WORK(1),WORK(L1),WORK(L2),WORK(L3),IWORK,LAST)
!
!         call ERROR HANDLER if NECESSARY.
!
  LVL = 0
10    if ( IER == 6) LVL = 1
  if (IER  /=  0) call XERMSG ('SLATEC', 'QAGS', &
     'ABNORMAL RETURN', IER, LVL)
  return
end
