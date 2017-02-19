subroutine PFQAD (F, LDC, C, XI, LXI, K, ID, X1, X2, TOL, QUAD, &
     IERR)
!
!! PFQAD computes the integral on (X1,X2) of a product of a function F and ...
!  the ID-th derivative of a B-spline, (PP-representation).
!
!***LIBRARY   SLATEC
!***CATEGORY  H2A2A1, E3, K6
!***TYPE      SINGLE PRECISION (PFQAD-S, DPFQAD-D)
!***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract
!         PFQAD computes the integral on (X1,X2) of a product of a
!         function F and the ID-th derivative of a B-spline, using the
!         PP-representation (C,XI,LXI,K). (X1,X2) is normally a sub-
!         interval of XI(1)  <=  X  <=  XI(LXI+1).  An integration rou-
!         tine, PPGQ8(a modification of GAUS8), integrates the product
!         on sub-intervals of (X1,X2) formed by the included break
!         points.  Integration outside of (XI(1),XI(LXI+1)) is permitted
!         provided F is defined.
!
!     Description of Arguments
!         Input
!           F      - external function of one argument for the
!                    integrand PF(X)=F(X)*PPVAL(LDC,C,XI,LXI,K,ID,X,
!                    INPPV)
!           LDC    - leading dimension of matrix C, LDC  >=  K
!           C(I,J) - right Taylor derivatives at XI(J), I=1,K , J=1,LXI
!           XI(*)  - break point array of length LXI+1
!           LXI    - number of polynomial pieces
!           K      - order of B-spline, K  >=  1
!           ID     - order of the spline derivative, 0  <=  ID  <=  K-1
!                    ID=0 gives the spline function
!           X1,X2  - end points of quadrature interval, normally in
!                    XI(1)  <=  X  <=  XI(LXI+1)
!           TOL    - desired accuracy for the quadrature, suggest
!                    10.*STOL  <  TOL  <=  0.1 where STOL is the single
!                    precision unit roundoff for the machine = R1MACH(4)
!
!         Output
!           QUAD   - integral of PF(X) on (X1,X2)
!           IERR   - a status code
!                    IERR=1 normal return
!                         2 some quadrature does not meet the
!                           requested tolerance
!
!     Error Conditions
!         TOL not greater than the single precision unit roundoff or
!         less than 0.1 is a fatal error.
!         Some quadrature does not meet the requested tolerance.
!
!***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
!                 B-splines, Report SAND79-1825, Sandia Laboratories,
!                 December 1979.
!***ROUTINES CALLED  INTRV, PPGQ8, R1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  PFQAD
!
  INTEGER ID,IERR,IFLG,ILO,IL1,IL2,INPPV,K,LDC,LEFT,LXI,MF1,MF2
  REAL A, AA, ANS, B, BB, C, Q, QUAD, TA, TB, TOL, WTOL, XI, X1, X2
  REAL R1MACH, F
  DIMENSION XI(*), C(LDC,*)
  EXTERNAL F
!
!***FIRST EXECUTABLE STATEMENT  PFQAD
  IERR = 1
  QUAD = 0.0E0
  if ( K < 1) go to 100
  if ( LDC < K) go to 105
  if ( ID < 0 .OR. ID >= K) go to 110
  if ( LXI < 1) go to 115
  WTOL = R1MACH(4)
  if (TOL < WTOL .OR. TOL > 0.1E0) go to 20
  AA = MIN(X1,X2)
  BB = MAX(X1,X2)
  if (AA == BB) RETURN
  ILO = 1
  call INTRV(XI, LXI, AA, ILO, IL1, MF1)
  call INTRV(XI, LXI, BB, ILO, IL2, MF2)
  Q = 0.0E0
  INPPV = 1
  DO 10 LEFT=IL1,IL2
    TA = XI(LEFT)
    A = MAX(AA,TA)
    if (LEFT == 1) A = AA
    TB = BB
    if (LEFT < LXI) TB = XI(LEFT+1)
    B = MIN(BB,TB)
    call PPGQ8(F,LDC,C,XI,LXI,K,ID,A,B,INPPV,TOL,ANS,IFLG)
    if (IFLG > 1) IERR = 2
    Q = Q + ANS
   10 CONTINUE
  if (X1 > X2) Q = -Q
  QUAD = Q
  return
!
   20 CONTINUE
  call XERMSG ('SLATEC', 'PFQAD', &
     'TOL IS LESS THAN THE SINGLE PRECISION TOLERANCE OR ' // &
     'GREATER THAN 0.1', 2, 1)
  return
  100 CONTINUE
  call XERMSG ('SLATEC', 'PFQAD', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
  105 CONTINUE
  call XERMSG ('SLATEC', 'PFQAD', 'LDC DOES NOT SATISFY LDC >= K', &
     2, 1)
  return
  110 CONTINUE
  call XERMSG ('SLATEC', 'PFQAD', &
     'ID DOES NOT SATISFY 0 <= ID < K', 2, 1)
  return
  115 CONTINUE
  call XERMSG ('SLATEC', 'PFQAD', 'LXI DOES NOT SATISFY LXI >= 1', &
     2, 1)
  return
end
