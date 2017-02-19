subroutine BFQAD (F, T, BCOEF, N, K, ID, X1, X2, TOL, QUAD, IERR, WORK)
!
!! BFQAD computes the integral of a product of a function and a ...
!  derivative of a B-spline.
!
!***LIBRARY   SLATEC
!***CATEGORY  H2A2A1, E3, K6
!***TYPE      SINGLE PRECISION (BFQAD-S, DBFQAD-D)
!***KEYWORDS  INTEGRAL OF B-SPLINE, QUADRATURE
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!         BFQAD computes the integral on (X1,X2) of a product of a
!         function F and the ID-th derivative of a K-th order B-spline,
!         using the B-representation (T,BCOEF,N,K).  (X1,X2) must be
!         a subinterval of T(K)  <=  X .le. T(N+1).  An integration
!         routine BSGQ8 (a modification
!         of GAUS8), integrates the product on sub-
!         intervals of (X1,X2) formed by included (distinct) knots.
!
!     Description of Arguments
!         Input
!           F      - external function of one argument for the
!                    integrand BF(X)=F(X)*BVALU(T,BCOEF,N,K,ID,X,INBV,
!                    WORK)
!           T      - knot array of length N+K
!           BCOEF  - coefficient array of length N
!           N      - length of coefficient array
!           K      - order of B-spline, K  >=  1
!           ID     - order of the spline derivative, 0  <=  ID  <=  K-1
!                    ID=0 gives the spline function
!           X1,X2  - end points of quadrature interval in
!                    T(K)  <=  X  <=  T(N+1)
!           TOL    - desired accuracy for the quadrature, suggest
!                    10.*STOL  <  TOL  <=  0.1 where STOL is the single
!                    precision unit roundoff for the machine = R1MACH(4)
!
!         Output
!           QUAD   - integral of BF(X) on (X1,X2)
!           IERR   - a status code
!                    IERR=1  normal return
!                         2  some quadrature on (X1,X2) does not meet
!                            the requested tolerance.
!           WORK   - work vector of length 3*K
!
!     Error Conditions
!         X1 or X2 not in T(K)  <=  X  <=  T(N+1) is a fatal error.
!         TOL not greater than the single precision unit roundoff or
!         less than 0.1 is a fatal error.
!         Some quadrature fails to meet the requested tolerance.
!
!***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
!                 B-splines, Report SAND79-1825, Sandia Laboratories,
!                 December 1979.
!***ROUTINES CALLED  BSGQ8, INTRV, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BFQAD
!
!
  INTEGER ID, IERR, IFLG, ILO, IL1, IL2, K, LEFT, MFLAG, N, NPK, NP1
  REAL A,AA,ANS,B,BB,BCOEF,Q,QUAD,T,TA,TB,TOL,WORK,WTOL, X1, &
   X2
  REAL R1MACH, F
  DIMENSION T(*), BCOEF(*), WORK(*)
  EXTERNAL F
!***FIRST EXECUTABLE STATEMENT  BFQAD
  IERR = 1
  QUAD = 0.0E0
  if ( K < 1) go to 100
  if ( N < K) go to 105
  if ( ID < 0 .OR. ID >= K) go to 110
  WTOL = R1MACH(4)
  if (TOL < WTOL .OR. TOL > 0.1E0) go to 30
  AA = MIN(X1,X2)
  BB = MAX(X1,X2)
  if (AA < T(K)) go to 20
  NP1 = N + 1
  if (BB > T(NP1)) go to 20
  if (AA == BB) RETURN
  NPK = N + K
!
  ILO = 1
  call INTRV(T, NPK, AA, ILO, IL1, MFLAG)
  call INTRV(T, NPK, BB, ILO, IL2, MFLAG)
  if (IL2 >= NP1) IL2 = N
  INBV = 1
  Q = 0.0E0
  DO 10 LEFT=IL1,IL2
    TA = T(LEFT)
    TB = T(LEFT+1)
    if (TA == TB) go to 10
    A = MAX(AA,TA)
    B = MIN(BB,TB)
    call BSGQ8(F,T,BCOEF,N,K,ID,A,B,INBV,TOL,ANS,IFLG,WORK)
    if (IFLG > 1) IERR = 2
    Q = Q + ANS
   10 CONTINUE
  if (X1 > X2) Q = -Q
  QUAD = Q
  return
!
!
   20 CONTINUE
  call XERMSG ('SLATEC', 'BFQAD', &
     'X1 OR X2 OR BOTH DO NOT SATISFY T(K) <= X <= T(N+1)', 2, 1)
  return
   30 CONTINUE
  call XERMSG ('SLATEC', 'BFQAD', &
     'TOL IS LESS THAN THE SINGLE PRECISION TOLERANCE OR ' // &
     'GREATER THAN 0.1', 2, 1)
  return
  100 CONTINUE
  call XERMSG ('SLATEC', 'BFQAD', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
  105 CONTINUE
  call XERMSG ('SLATEC', 'BFQAD', 'N DOES NOT SATISFY N >= K', 2, &
     1)
  return
  110 CONTINUE
  call XERMSG ('SLATEC', 'BFQAD', &
     'ID DOES NOT SATISFY 0  <=  ID  <  K', 2, 1)
  return
end
