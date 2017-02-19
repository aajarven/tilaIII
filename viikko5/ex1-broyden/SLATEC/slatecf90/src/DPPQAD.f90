subroutine DPPQAD (LDC, C, XI, LXI, K, X1, X2, PQUAD)
!
!! DPPQAD computes the integral on (X1,X2) of a K-th order B-spline ...
!  using the piecewise polynomial (PP) representation.
!
!***LIBRARY   SLATEC
!***CATEGORY  H2A2A1, E3, K6
!***TYPE      DOUBLE PRECISION (PPQAD-S, DPPQAD-D)
!***KEYWORDS  B-SPLINE, DATA FITTING, INTERPOLATION, QUADRATURE, SPLINES
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract    **** a double precision routine ****
!         DPPQAD computes the integral on (X1,X2) of a K-th order
!         B-spline using the piecewise polynomial representation
!         (C,XI,LXI,K).  Here the Taylor expansion about the left
!         end point XI(J) of the J-th interval is integrated and
!         evaluated on subintervals of (X1,X2) which are formed by
!         included break points.  Integration outside (XI(1),XI(LXI+1))
!         is permitted.
!
!     Description of Arguments
!         Input      C,XI,X1,X2 are double precision
!           LDC    - leading dimension of matrix C, LDC  >=  K
!           C(I,J) - right Taylor derivatives at XI(J), I=1,K , J=1,LXI
!           XI(*)  - break point array of length LXI+1
!           LXI    - number of polynomial pieces
!           K      - order of B-spline, K  >=  1
!           X1,X2  - end points of quadrature interval, normally in
!                    XI(1)  <=  X  <=  XI(LXI+1)
!
!         Output     PQUAD is double precision
!           PQUAD  - integral of the PP representation over (X1,X2)
!
!     Error Conditions
!         Improper input is a fatal error
!
!***REFERENCES  D. E. Amos, Quadrature subroutines for splines and
!                 B-splines, Report SAND79-1825, Sandia Laboratories,
!                 December 1979.
!***ROUTINES CALLED  DINTRV, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPPQAD
!
  INTEGER I, II, IL, ILO, IL1, IL2, IM, K, LDC, LEFT, LXI, MF1, MF2
  DOUBLE PRECISION A,AA,BB,C,DX,FLK,PQUAD,Q,S,SS,TA,TB,X,XI,X1,X2
  DIMENSION XI(*), C(LDC,*), SS(2)
!
!***FIRST EXECUTABLE STATEMENT  DPPQAD
  PQUAD = 0.0D0
  if ( K < 1) go to 100
  if ( LXI < 1) go to 105
  if ( LDC < K) go to 110
  AA = MIN(X1,X2)
  BB = MAX(X1,X2)
  if (AA == BB) RETURN
  ILO = 1
  call DINTRV(XI, LXI, AA, ILO, IL1, MF1)
  call DINTRV(XI, LXI, BB, ILO, IL2, MF2)
  Q = 0.0D0
  DO 40 LEFT=IL1,IL2
    TA = XI(LEFT)
    A = MAX(AA,TA)
    if (LEFT == 1) A = AA
    TB = BB
    if (LEFT < LXI) TB = XI(LEFT+1)
    X = MIN(BB,TB)
    DO 30 II=1,2
      SS(II) = 0.0D0
      DX = X - XI(LEFT)
      if (DX == 0.0D0) go to 20
      S = C(K,LEFT)
      FLK = K
      IM = K - 1
      IL = IM
      DO 10 I=1,IL
        S = S*DX/FLK + C(IM,LEFT)
        IM = IM - 1
        FLK = FLK - 1.0D0
   10     CONTINUE
      SS(II) = S*DX
   20     CONTINUE
      X = A
   30   CONTINUE
    Q = Q + (SS(1)-SS(2))
   40 CONTINUE
  if (X1 > X2) Q = -Q
  PQUAD = Q
  return
!
!
  100 CONTINUE
  call XERMSG ('SLATEC', 'DPPQAD', 'K DOES NOT SATISFY K >= 1', 2, &
     1)
  return
  105 CONTINUE
  call XERMSG ('SLATEC', 'DPPQAD', 'LXI DOES NOT SATISFY LXI >= 1', &
     2, 1)
  return
  110 CONTINUE
  call XERMSG ('SLATEC', 'DPPQAD', 'LDC DOES NOT SATISFY LDC >= K', &
     2, 1)
  return
end
