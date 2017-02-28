subroutine DBSQAD (T, BCOEF, N, K, X1, X2, BQUAD, WORK)
!
!! DBSQAD computes the integral of a K-th order B-spline using the ...
!  B-representation.
!
!***LIBRARY   SLATEC
!***CATEGORY  H2A2A1, E3, K6
!***TYPE      DOUBLE PRECISION (BSQAD-S, DBSQAD-D)
!***KEYWORDS  INTEGRAL OF B-SPLINES, QUADRATURE
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract    **** a double precision routine ****
!
!         DBSQAD computes the integral on (X1,X2) of a K-th order
!         B-spline using the B-representation (T,BCOEF,N,K).  Orders
!         K as high as 20 are permitted by applying a 2, 6, or 10
!         point Gauss formula on subintervals of (X1,X2) which are
!         formed by included (distinct) knots.
!
!         If orders K greater than 20 are needed, use DBFQAD with
!         F(X) = 1.
!
!         The maximum number of significant digits obtainable in
!         DBSQAD is the smaller of 18 and the number of digits
!         carried in double precision arithmetic.
!
!     Description of Arguments
!         Input      T,BCOEF,X1,X2 are double precision
!           T      - knot array of length N+K
!           BCOEF  - B-spline coefficient array of length N
!           N      - length of coefficient array
!           K      - order of B-spline, 1  <=  K  <=  20
!           X1,X2  - end points of quadrature interval in
!                    T(K)  <=  X  <=  T(N+1)
!
!         Output     BQUAD,WORK are double precision
!           BQUAD  - integral of the B-spline over (X1,X2)
!           WORK   - work vector of length 3*K
!
!     Error Conditions
!         Improper input is a fatal error
!
!***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
!                 B-splines, Report SAND79-1825, Sandia Laboratories,
!                 December 1979.
!***ROUTINES CALLED  DBVALU, DINTRV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBSQAD
!
  INTEGER I,IL1,IL2,ILO,INBV, JF,K,LEFT,M,MF,MFLAG,N, NPK, NP1
  DOUBLE PRECISION A,AA,B,BB,BCOEF,BMA,BPA,BQUAD,C1,GPTS,GWTS,GX,Q, &
   SUM, T, TA, TB, WORK, X1, X2, Y1, Y2
  DOUBLE PRECISION DBVALU
  DIMENSION T(*), BCOEF(*), GPTS(9), GWTS(9), SUM(5), WORK(*)
!
  SAVE GPTS, GWTS
  DATA GPTS(1), GPTS(2), GPTS(3), GPTS(4), GPTS(5), GPTS(6), &
       GPTS(7), GPTS(8), GPTS(9)/ &
       5.77350269189625764D-01,     2.38619186083196909D-01, &
       6.61209386466264514D-01,     9.32469514203152028D-01, &
       1.48874338981631211D-01,     4.33395394129247191D-01, &
       6.79409568299024406D-01,     8.65063366688984511D-01, &
       9.73906528517171720D-01/
  DATA GWTS(1), GWTS(2), GWTS(3), GWTS(4), GWTS(5), GWTS(6), &
       GWTS(7), GWTS(8), GWTS(9)/ &
       1.00000000000000000D+00,     4.67913934572691047D-01, &
       3.60761573048138608D-01,     1.71324492379170345D-01, &
       2.95524224714752870D-01,     2.69266719309996355D-01, &
       2.19086362515982044D-01,     1.49451349150580593D-01, &
       6.66713443086881376D-02/
!
!***FIRST EXECUTABLE STATEMENT  DBSQAD
  BQUAD = 0.0D0
  if ( K < 1 .OR. K > 20) go to 65
  if ( N < K) go to 70
  AA = MIN(X1,X2)
  BB = MAX(X1,X2)
  if (AA < T(K)) go to 60
  NP1 = N + 1
  if (BB > T(NP1)) go to 60
  if (AA == BB) RETURN
  NPK = N + K
!     SELECTION OF 2, 6, OR 10 POINT GAUSS FORMULA
  JF = 0
  MF = 1
  if (K <= 4) go to 10
  JF = 1
  MF = 3
  if (K <= 12) go to 10
  JF = 4
  MF = 5
   10 CONTINUE
!
  DO 20 I=1,MF
    SUM(I) = 0.0D0
   20 CONTINUE
  ILO = 1
  INBV = 1
  call DINTRV(T, NPK, AA, ILO, IL1, MFLAG)
  call DINTRV(T, NPK, BB, ILO, IL2, MFLAG)
  if (IL2 >= NP1) IL2 = N
  DO 40 LEFT=IL1,IL2
    TA = T(LEFT)
    TB = T(LEFT+1)
    if (TA == TB) go to 40
    A = MAX(AA,TA)
    B = MIN(BB,TB)
    BMA = 0.5D0*(B-A)
    BPA = 0.5D0*(B+A)
    DO 30 M=1,MF
      C1 = BMA*GPTS(JF+M)
      GX = -C1 + BPA
      Y2 = DBVALU(T,BCOEF,N,K,0,GX,INBV,WORK)
      GX = C1 + BPA
      Y1 = DBVALU(T,BCOEF,N,K,0,GX,INBV,WORK)
      SUM(M) = SUM(M) + (Y1+Y2)*BMA
   30   CONTINUE
   40 CONTINUE
  Q = 0.0D0
  DO 50 M=1,MF
    Q = Q + GWTS(JF+M)*SUM(M)
   50 CONTINUE
  if (X1 > X2) Q = -Q
  BQUAD = Q
  return
!
!
   60 CONTINUE
  call XERMSG ('SLATEC', 'DBSQAD', &
     'X1 OR X2 OR BOTH DO NOT SATISFY T(K) <= X <= T(N+1)', 2, 1)
  return
   65 CONTINUE
  call XERMSG ('SLATEC', 'DBSQAD', &
     'K DOES NOT SATISFY 1 <= K <= 20', 2, 1)
  return
   70 CONTINUE
  call XERMSG ('SLATEC', 'DBSQAD', 'N DOES NOT SATISFY N >= K', 2, &
     1)
  return
end