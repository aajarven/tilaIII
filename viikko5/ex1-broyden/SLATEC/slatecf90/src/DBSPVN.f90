subroutine DBSPVN (T, JHIGH, K, INDEX, X, ILEFT, VNIKX, WORK, &
     IWORK)
!
!! DBSPVN calculates the value of all (possibly) nonzero basis functions at X.
!
!***LIBRARY   SLATEC
!***CATEGORY  E3, K6
!***TYPE      DOUBLE PRECISION (BSPVN-S, DBSPVN-D)
!***KEYWORDS  EVALUATION OF B-SPLINE
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Written by Carl de Boor and modified by D. E. Amos
!
!     Abstract    **** a double precision routine ****
!         DBSPVN is the BSPLVN routine of the reference.
!
!         DBSPVN calculates the value of all (possibly) nonzero basis
!         functions at X of order MAX(JHIGH,(J+1)*(INDEX-1)), where T(K)
!          <=  X  <=  T(N+1) and J=IWORK is set inside the routine on
!         the first call when INDEX=1.  ILEFT is such that T(ILEFT)  <=
!         X  <  T(ILEFT+1).  A call to DINTRV(T,N+1,X,ILO,ILEFT,MFLAG)
!         produces the proper ILEFT.  DBSPVN calculates using the basic
!         algorithm needed in DBSPVD.  If only basis functions are
!         desired, setting JHIGH=K and INDEX=1 can be faster than
!         calling DBSPVD, but extra coding is required for derivatives
!         (INDEX=2) and DBSPVD is set up for this purpose.
!
!         Left limiting values are set up as described in DBSPVD.
!
!     Description of Arguments
!
!         Input      T,X are double precision
!          T       - knot vector of length N+K, where
!                    N = number of B-spline basis functions
!                    N = sum of knot multiplicities-K
!          JHIGH   - order of B-spline, 1  <=  JHIGH  <=  K
!          K       - highest possible order
!          INDEX   - INDEX = 1 gives basis functions of order JHIGH
!                          = 2 denotes previous entry with work, IWORK
!                              values saved for subsequent calls to
!                              DBSPVN.
!          X       - argument of basis functions,
!                    T(K)  <=  X  <=  T(N+1)
!          ILEFT   - largest integer such that
!                    T(ILEFT)  <=  X  <   T(ILEFT+1)
!
!         Output     VNIKX, WORK are double precision
!          VNIKX   - vector of length K for spline values.
!          WORK    - a work vector of length 2*K
!          IWORK   - a work parameter.  Both WORK and IWORK contain
!                    information necessary to continue for INDEX = 2.
!                    When INDEX = 1 exclusively, these are scratch
!                    variables and can be used for other purposes.
!
!     Error Conditions
!         Improper input is a fatal error.
!
!***REFERENCES  Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!***ROUTINES CALLED  XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBSPVN
!
  INTEGER ILEFT, IMJP1, INDEX, IPJ, IWORK, JHIGH, JP1, JP1ML, K, L
  DOUBLE PRECISION T, VM, VMPREV, VNIKX, WORK, X
!     DIMENSION T(ILEFT+JHIGH)
  DIMENSION T(*), VNIKX(*), WORK(*)
!     CONTENT OF J, DELTAM, DELTAP IS EXPECTED UNCHANGED BETWEEN CALLS.
!     WORK(I) = DELTAP(I), WORK(K+I) = DELTAM(I), I = 1,K
!***FIRST EXECUTABLE STATEMENT  DBSPVN
  if ( K < 1) go to 90
  if ( JHIGH > K .OR. JHIGH < 1) go to 100
  if ( INDEX < 1 .OR. INDEX > 2) go to 105
  if ( X < T(ILEFT) .OR. X > T(ILEFT+1)) go to 110
  go to (10, 20), INDEX
   10 IWORK = 1
  VNIKX(1) = 1.0D0
  if (IWORK >= JHIGH) go to 40
!
   20 IPJ = ILEFT + IWORK
  WORK(IWORK) = T(IPJ) - X
  IMJP1 = ILEFT - IWORK + 1
  WORK(K+IWORK) = X - T(IMJP1)
  VMPREV = 0.0D0
  JP1 = IWORK + 1
  DO 30 L=1,IWORK
    JP1ML = JP1 - L
    VM = VNIKX(L)/(WORK(L)+WORK(K+JP1ML))
    VNIKX(L) = VM*WORK(L) + VMPREV
    VMPREV = VM*WORK(K+JP1ML)
   30 CONTINUE
  VNIKX(JP1) = VMPREV
  IWORK = JP1
  if (IWORK < JHIGH) go to 20
!
   40 RETURN
!
!
   90 CONTINUE
  call XERMSG ('SLATEC', 'DBSPVN', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
  100 CONTINUE
  call XERMSG ('SLATEC', 'DBSPVN', &
     'JHIGH DOES NOT SATISFY 1 <= JHIGH <= K', 2, 1)
  return
  105 CONTINUE
  call XERMSG ('SLATEC', 'DBSPVN', 'INDEX IS NOT 1 OR 2', 2, 1)
  return
  110 CONTINUE
  call XERMSG ('SLATEC', 'DBSPVN', &
     'X DOES NOT SATISFY T(ILEFT) <= X <= T(ILEFT+1)', 2, 1)
  return
end
