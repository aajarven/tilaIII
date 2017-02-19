subroutine DBINT4 (X, Y, NDATA, IBCL, IBCR, FBCL, FBCR, KNTOPT, T, &
     BCOEF, N, K, W)
!
!! DBINT4 computes the B-representation of a cubic spline ...
!  which interpolates given data.
!
!***LIBRARY   SLATEC
!***CATEGORY  E1A
!***TYPE      DOUBLE PRECISION (BINT4-S, DBINT4-D)
!***KEYWORDS  B-SPLINE, CUBIC SPLINES, DATA FITTING, INTERPOLATION
!***AUTHOR  Amos, D. E., (SNLA)
!***DESCRIPTION
!
!     Abstract    **** a double precision routine ****
!
!         DBINT4 computes the B representation (T,BCOEF,N,K) of a
!         cubic spline (K=4) which interpolates data (X(I),Y(I)),
!         I=1,NDATA.  Parameters IBCL, IBCR, FBCL, FBCR allow the
!         specification of the spline first or second derivative at
!         both X(1) and X(NDATA).  When this data is not specified
!         by the problem, it is common practice to use a natural
!         spline by setting second derivatives at X(1) and X(NDATA)
!         to zero (IBCL=IBCR=2,FBCL=FBCR=0.0).  The spline is defined
!         on T(4)  <=  X  <=  T(N+1) with (ordered) interior knots at
!         X(I) values where N=NDATA+2.  The knots T(1),T(2),T(3) lie to
!         the left of T(4)=X(1) and the knots T(N+2), T(N+3), T(N+4)
!         lie to the right of T(N+1)=X(NDATA) in increasing order.  If
!         no extrapolation outside (X(1),X(NDATA)) is anticipated, the
!         knots T(1)=T(2)=T(3)=T(4)=X(1) and T(N+2)=T(N+3)=T(N+4)=
!         T(N+1)=X(NDATA) can be specified by KNTOPT=1.  KNTOPT=2
!         selects a knot placement for T(1), T(2), T(3) to make the
!         first 7 knots symmetric about T(4)=X(1) and similarly for
!         T(N+2), T(N+3), T(N+4) about T(N+1)=X(NDATA).  KNTOPT=3
!         allows the user to make his own selection, in increasing
!         order, for T(1), T(2), T(3) to the left of X(1) and T(N+2),
!         T(N+3), T(N+4) to the right of X(NDATA) in the work array
!         W(1) through W(6).  In any case, the interpolation on
!         T(4)  <=  X  <=  T(N+1) by using function DBVALU is unique
!         for given boundary conditions.
!
!     Description of Arguments
!
!         Input      X,Y,FBCL,FBCR,W are double precision
!           X      - X vector of abscissae of length NDATA, distinct
!                    and in increasing order
!           Y      - Y vector of ordinates of length NDATA
!           NDATA  - number of data points, NDATA  >=  2
!           IBCL   - selection parameter for left boundary condition
!                    IBCL = 1 constrain the first derivative at
!                             X(1) to FBCL
!                         = 2 constrain the second derivative at
!                             X(1) to FBCL
!           IBCR   - selection parameter for right boundary condition
!                    IBCR = 1 constrain first derivative at
!                             X(NDATA) to FBCR
!                    IBCR = 2 constrain second derivative at
!                             X(NDATA) to FBCR
!           FBCL   - left boundary values governed by IBCL
!           FBCR   - right boundary values governed by IBCR
!           KNTOPT - knot selection parameter
!                    KNTOPT = 1 sets knot multiplicity at T(4) and
!                               T(N+1) to 4
!                           = 2 sets a symmetric placement of knots
!                               about T(4) and T(N+1)
!                           = 3 sets T(I)=W(I) and T(N+1+I)=W(3+I),I=1,3
!                               where W(I),I=1,6 is supplied by the user
!           W      - work array of dimension at least 5*(NDATA+2)
!                    If KNTOPT=3, then W(1),W(2),W(3) are knot values to
!                    the left of X(1) and W(4),W(5),W(6) are knot
!                    values to the right of X(NDATA) in increasing
!                    order to be supplied by the user
!
!         Output     T,BCOEF are double precision
!           T      - knot array of length N+4
!           BCOEF  - B spline coefficient array of length N
!           N      - number of coefficients, N=NDATA+2
!           K      - order of spline, K=4
!
!     Error Conditions
!         Improper  input is a fatal error
!         Singular system of equations is a fatal error
!
!***REFERENCES  D. E. Amos, Computation with splines and B-splines,
!                 Report SAND78-1968, Sandia Laboratories, March 1979.
!               Carl de Boor, Package for calculating with B-splines,
!                 SIAM Journal on Numerical Analysis 14, 3 (June 1977),
!                 pp. 441-472.
!               Carl de Boor, A Practical Guide to Splines, Applied
!                 Mathematics Series 27, Springer-Verlag, New York,
!                 1978.
!***ROUTINES CALLED  D1MACH, DBNFAC, DBNSLV, DBSPVD, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBINT4
!
  INTEGER I, IBCL, IBCR, IFLAG, ILB, ILEFT, IT, IUB, IW, IWP, J, &
   JW, K, KNTOPT, N, NDATA, NDM, NP, NWROW
  DOUBLE PRECISION BCOEF,FBCL,FBCR,T,TOL,TXN,TX1,VNIKX,W,WDTOL, &
   WORK,X,XL,Y
  DOUBLE PRECISION D1MACH
  DIMENSION X(*), Y(*), T(*), BCOEF(*), W(5,*), VNIKX(4,4), WORK(15)
!***FIRST EXECUTABLE STATEMENT  DBINT4
  WDTOL = D1MACH(4)
  TOL = SQRT(WDTOL)
  if (NDATA < 2) go to 200
  NDM = NDATA - 1
  DO 10 I=1,NDM
    if (X(I) >= X(I+1)) go to 210
   10 CONTINUE
  if (IBCL < 1 .OR. IBCL > 2) go to 220
  if (IBCR < 1 .OR. IBCR > 2) go to 230
  if (KNTOPT < 1 .OR. KNTOPT > 3) go to 240
  K = 4
  N = NDATA + 2
  NP = N + 1
  DO 20 I=1,NDATA
    T(I+3) = X(I)
   20 CONTINUE
  go to (30, 50, 90), KNTOPT
!     SET UP KNOT ARRAY WITH MULTIPLICITY 4 AT X(1) AND X(NDATA)
   30 CONTINUE
  DO 40 I=1,3
    T(4-I) = X(1)
    T(NP+I) = X(NDATA)
   40 CONTINUE
  go to 110
!     SET UP KNOT ARRAY WITH SYMMETRIC PLACEMENT ABOUT END POINTS
   50 CONTINUE
  if (NDATA > 3) go to 70
  XL = (X(NDATA)-X(1))/3.0D0
  DO 60 I=1,3
    T(4-I) = T(5-I) - XL
    T(NP+I) = T(NP+I-1) + XL
   60 CONTINUE
  go to 110
   70 CONTINUE
  TX1 = X(1) + X(1)
  TXN = X(NDATA) + X(NDATA)
  DO 80 I=1,3
    T(4-I) = TX1 - X(I+1)
    T(NP+I) = TXN - X(NDATA-I)
   80 CONTINUE
  go to 110
!     SET UP KNOT ARRAY LESS THAN X(1) AND GREATER THAN X(NDATA) TO BE
!     SUPPLIED BY USER IN WORK LOCATIONS W(1) THROUGH W(6) WHEN KNTOPT=3
   90 CONTINUE
  DO 100 I=1,3
    T(4-I) = W(4-I,1)
    JW = MAX(1,I-1)
    IW = MOD(I+2,5)+1
    T(NP+I) = W(IW,JW)
    if (T(4-I) > T(5-I)) go to 250
    if (T(NP+I) < T(NP+I-1)) go to 250
  100 CONTINUE
  110 CONTINUE
!
  DO 130 I=1,5
    DO 120 J=1,N
      W(I,J) = 0.0D0
  120   CONTINUE
  130 CONTINUE
!     SET UP LEFT INTERPOLATION POINT AND LEFT BOUNDARY CONDITION FOR
!     RIGHT LIMITS
  IT = IBCL + 1
  call DBSPVD(T, K, IT, X(1), K, 4, VNIKX, WORK)
  IW = 0
  if (ABS(VNIKX(3,1)) < TOL) IW = 1
  DO 140 J=1,3
    W(J+1,4-J) = VNIKX(4-J,IT)
    W(J,4-J) = VNIKX(4-J,1)
  140 CONTINUE
  BCOEF(1) = Y(1)
  BCOEF(2) = FBCL
!     SET UP INTERPOLATION EQUATIONS FOR POINTS I=2 TO I=NDATA-1
  ILEFT = 4
  if (NDM < 2) go to 170
  DO 160 I=2,NDM
    ILEFT = ILEFT + 1
    call DBSPVD(T, K, 1, X(I), ILEFT, 4, VNIKX, WORK)
    DO 150 J=1,3
      W(J+1,3+I-J) = VNIKX(4-J,1)
  150   CONTINUE
    BCOEF(I+1) = Y(I)
  160 CONTINUE
!     SET UP RIGHT INTERPOLATION POINT AND RIGHT BOUNDARY CONDITION FOR
!     LEFT LIMITS(ILEFT IS ASSOCIATED WITH T(N)=X(NDATA-1))
  170 CONTINUE
  IT = IBCR + 1
  call DBSPVD(T, K, IT, X(NDATA), ILEFT, 4, VNIKX, WORK)
  JW = 0
  if (ABS(VNIKX(2,1)) < TOL) JW = 1
  DO 180 J=1,3
    W(J+1,3+NDATA-J) = VNIKX(5-J,IT)
    W(J+2,3+NDATA-J) = VNIKX(5-J,1)
  180 CONTINUE
  BCOEF(N-1) = FBCR
  BCOEF(N) = Y(NDATA)
!     SOLVE SYSTEM OF EQUATIONS
  ILB = 2 - JW
  IUB = 2 - IW
  NWROW = 5
  IWP = IW + 1
  call DBNFAC(W(IWP,1), NWROW, N, ILB, IUB, IFLAG)
  if (IFLAG == 2) go to 190
  call DBNSLV(W(IWP,1), NWROW, N, ILB, IUB, BCOEF)
  return
!
!
  190 CONTINUE
  call XERMSG ('SLATEC', 'DBINT4', &
     'THE SYSTEM OF EQUATIONS IS SINGULAR', 2, 1)
  return
  200 CONTINUE
  call XERMSG ('SLATEC', 'DBINT4', 'NDATA IS LESS THAN 2', 2, 1)
  return
  210 CONTINUE
  call XERMSG ('SLATEC', 'DBINT4', &
     'X VALUES ARE NOT DISTINCT OR NOT ORDERED', 2, 1)
  return
  220 CONTINUE
  call XERMSG ('SLATEC', 'DBINT4', 'IBCL IS NOT 1 OR 2', 2, 1)
  return
  230 CONTINUE
  call XERMSG ('SLATEC', 'DBINT4', 'IBCR IS NOT 1 OR 2', 2, 1)
  return
  240 CONTINUE
  call XERMSG ('SLATEC', 'DBINT4', 'KNTOPT IS NOT 1, 2, OR 3', 2, &
     1)
  return
  250 CONTINUE
  call XERMSG ('SLATEC', 'DBINT4', &
     'KNOT INPUT THROUGH W ARRAY IS NOT ORDERED PROPERLY', 2, 1)
  return
end
