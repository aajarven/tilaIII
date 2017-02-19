subroutine LPDP (A, MDA, M, N1, N2, PRGOPT, X, WNORM, MODE, WS, IS)
!
!! LPDP is subsidiary to LSEI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (LPDP-S, DLPDP-D)
!***AUTHOR  Hanson, R. J., (SNLA)
!           Haskell, K. H., (SNLA)
!***DESCRIPTION
!
!     DIMENSION A(MDA,N+1),PRGOPT(*),X(N),WS((M+2)*(N+7)),IS(M+N+1),
!     where N=N1+N2.  This is a slight overestimate for WS(*).
!
!     Determine an N1-vector W, and
!               an N2-vector Z
!     which minimizes the Euclidean length of W
!     subject to G*W+H*Z  >=  Y.
!     This is the least projected distance problem, LPDP.
!     The matrices G and H are of respective
!     dimensions M by N1 and M by N2.
!
!     Called by subprogram LSI( ).
!
!     The matrix
!                (G H Y)
!
!     occupies rows 1,...,M and cols 1,...,N1+N2+1 of A(*,*).
!
!     The solution (W) is returned in X(*).
!                  (Z)
!
!     The value of MODE indicates the status of
!     the computation after returning to the user.
!
!          MODE=1  The solution was successfully obtained.
!
!          MODE=2  The inequalities are inconsistent.
!
!***SEE ALSO  LSEI
!***ROUTINES CALLED  SCOPY, SDOT, SNRM2, SSCAL, WNNLS
!***REVISION HISTORY  (YYMMDD)
!   790701  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910408  Updated the AUTHOR section.  (WRB)
!***END PROLOGUE  LPDP
!
!     SUBROUTINES CALLED
!
!     WNNLS         SOLVES A NONNEGATIVELY CONSTRAINED LINEAR LEAST
!                   SQUARES PROBLEM WITH LINEAR EQUALITY CONSTRAINTS.
!                   PART OF THIS PACKAGE.
!
!++
!     SDOT,         SUBROUTINES FROM THE BLAS PACKAGE.
!     SSCAL,SNRM2,  SEE TRANS. MATH. SOFT., VOL. 5, NO. 3, P. 308.
!     SCOPY
!
  REAL             A(MDA,*), PRGOPT(*), WS(*), WNORM, X(*)
  INTEGER IS(*)
  REAL             FAC, ONE, RNORM, SC, YNORM, ZERO
  REAL             SDOT, SNRM2
  SAVE ZERO, ONE, FAC
  DATA ZERO, ONE /0.E0,1.E0/, FAC /0.1E0/
!***FIRST EXECUTABLE STATEMENT  LPDP
  N = N1 + N2
  MODE = 1
  if (.NOT.(M <= 0)) go to 20
  if (.NOT.(N > 0)) go to 10
  X(1) = ZERO
  call SCOPY(N, X, 0, X, 1)
   10 WNORM = ZERO
  return
   20 NP1 = N + 1
!
!     SCALE NONZERO ROWS OF INEQUALITY MATRIX TO HAVE LENGTH ONE.
  DO 40 I=1,M
    SC = SNRM2(N,A(I,1),MDA)
    if (.NOT.(SC /= ZERO)) go to 30
    SC = ONE/SC
    call SSCAL(NP1, SC, A(I,1), MDA)
   30   CONTINUE
   40 CONTINUE
!
!     SCALE RT.-SIDE VECTOR TO HAVE LENGTH ONE (OR ZERO).
  YNORM = SNRM2(M,A(1,NP1),1)
  if (.NOT.(YNORM /= ZERO)) go to 50
  SC = ONE/YNORM
  call SSCAL(M, SC, A(1,NP1), 1)
!
!     SCALE COLS OF MATRIX H.
   50 J = N1 + 1
   60 if (.NOT.(J <= N)) go to 70
  SC = SNRM2(M,A(1,J),1)
  if (SC /= ZERO) SC = ONE/SC
  call SSCAL(M, SC, A(1,J), 1)
  X(J) = SC
  J = J + 1
  go to 60
   70 if (.NOT.(N1 > 0)) go to 130
!
!     COPY TRANSPOSE OF (H G Y) TO WORK ARRAY WS(*).
  IW = 0
  DO 80 I=1,M
!
!     MOVE COL OF TRANSPOSE OF H INTO WORK ARRAY.
    call SCOPY(N2, A(I,N1+1), MDA, WS(IW+1), 1)
    IW = IW + N2
!
!     MOVE COL OF TRANSPOSE OF G INTO WORK ARRAY.
    call SCOPY(N1, A(I,1), MDA, WS(IW+1), 1)
    IW = IW + N1
!
!     MOVE COMPONENT OF VECTOR Y INTO WORK ARRAY.
    WS(IW+1) = A(I,NP1)
    IW = IW + 1
   80 CONTINUE
  WS(IW+1) = ZERO
  call SCOPY(N, WS(IW+1), 0, WS(IW+1), 1)
  IW = IW + N
  WS(IW+1) = ONE
  IW = IW + 1
!
!     SOLVE EU=F SUBJECT TO (TRANSPOSE OF H)U=0, U >= 0.  THE
!     MATRIX E = TRANSPOSE OF (G Y), AND THE (N+1)-VECTOR
!     F = TRANSPOSE OF (0,...,0,1).
  IX = IW + 1
  IW = IW + M
!
!     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLS( ).
  IS(1) = 0
  IS(2) = 0
  call WNNLS(WS, NP1, N2, NP1-N2, M, 0, PRGOPT, WS(IX), RNORM, &
   MODEW, IS, WS(IW+1))
!
!     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY W.
  SC = ONE - SDOT(M,A(1,NP1),1,WS(IX),1)
  if (.NOT.(ONE+FAC*ABS(SC) /= ONE .AND. RNORM > ZERO)) go to 110
  SC = ONE/SC
  DO 90 J=1,N1
    X(J) = SC*SDOT(M,A(1,J),1,WS(IX),1)
   90 CONTINUE
!
!     COMPUTE THE VECTOR Q=Y-GW.  OVERWRITE Y WITH THIS VECTOR.
  DO 100 I=1,M
    A(I,NP1) = A(I,NP1) - SDOT(N1,A(I,1),MDA,X,1)
  100 CONTINUE
  go to 120
  110 MODE = 2
  return
  120 CONTINUE
  130 if (.NOT.(N2 > 0)) go to 180
!
!     COPY TRANSPOSE OF (H Q) TO WORK ARRAY WS(*).
  IW = 0
  DO 140 I=1,M
    call SCOPY(N2, A(I,N1+1), MDA, WS(IW+1), 1)
    IW = IW + N2
    WS(IW+1) = A(I,NP1)
    IW = IW + 1
  140 CONTINUE
  WS(IW+1) = ZERO
  call SCOPY(N2, WS(IW+1), 0, WS(IW+1), 1)
  IW = IW + N2
  WS(IW+1) = ONE
  IW = IW + 1
  IX = IW + 1
  IW = IW + M
!
!     SOLVE RV=S SUBJECT TO V >= 0.  THE MATRIX R =(TRANSPOSE
!     OF (H Q)), WHERE Q=Y-GW.  THE (N2+1)-VECTOR S =(TRANSPOSE
!     OF (0,...,0,1)).
!
!     DO NOT CHECK LENGTHS OF WORK ARRAYS IN THIS USAGE OF WNNLS( ).
  IS(1) = 0
  IS(2) = 0
  call WNNLS(WS, N2+1, 0, N2+1, M, 0, PRGOPT, WS(IX), RNORM, MODEW, &
   IS, WS(IW+1))
!
!     COMPUTE THE COMPONENTS OF THE SOLN DENOTED ABOVE BY Z.
  SC = ONE - SDOT(M,A(1,NP1),1,WS(IX),1)
  if (.NOT.(ONE+FAC*ABS(SC) /= ONE .AND. RNORM > ZERO)) go to 160
  SC = ONE/SC
  DO 150 J=1,N2
    L = N1 + J
    X(L) = SC*SDOT(M,A(1,L),1,WS(IX),1)*X(L)
  150 CONTINUE
  go to 170
  160 MODE = 2
  return
  170 CONTINUE
!
!     ACCOUNT FOR SCALING OF RT.-SIDE VECTOR IN SOLUTION.
  180 call SSCAL(N, YNORM, X, 1)
  WNORM = SNRM2(N1,X,1)
  return
end
