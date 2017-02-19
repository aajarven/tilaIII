subroutine D1UPDT (M, N, S, LS, U, V, W, SING)
!
!! D1UPDT is subsidiary to DNSQ and DNSQE.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (R1UPDT-S, D1UPDT-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an M by N lower trapezoidal matrix S, an M-vector U,
!     and an N-vector V, the problem is to determine an
!     orthogonal matrix Q such that
!
!                   t
!           (S + U*V )*Q
!
!     is again lower trapezoidal.
!
!     This subroutine determines Q as the product of 2*(N - 1)
!     transformations
!
!           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!     where GV(I), GW(I) are Givens rotations in the (I,N) plane
!     which eliminate elements in the I-th and N-th planes,
!     respectively. Q itself is not accumulated, rather the
!     information to recover the GV, GW rotations is returned.
!
!     The SUBROUTINE statement is
!
!       SUBROUTINE D1UPDT(M,N,S,LS,U,V,W,SING)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of S.
!
!       N is a positive integer input variable set to the number
!         of columns of S. N must not exceed M.
!
!       S is an array of length LS. On input S must contain the lower
!         trapezoidal matrix S stored by columns. On output S contains
!         the lower trapezoidal matrix produced as described above.
!
!       LS is a positive integer input variable not less than
!         (N*(2*M-N+1))/2.
!
!       U is an input array of length M which must contain the
!         vector U.
!
!       V is an array of length N. On input V must contain the vector
!         V. On output V(I) contains the information necessary to
!         recover the Givens rotation GV(I) described above.
!
!       W is an output array of length M. W(I) contains information
!         necessary to recover the Givens rotation GW(I) described
!         above.
!
!       SING is a LOGICAL output variable. SING is set TRUE if any
!         of the diagonal elements of the output S are zero. Otherwise
!         SING is set FALSE.
!
!***SEE ALSO  DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  D1UPDT
  DOUBLE PRECISION D1MACH
  INTEGER I, J, JJ, L, LS, M, N, NM1, NMJ
  DOUBLE PRECISION COS, COTAN, GIANT, ONE, P25, P5, S(*), &
       SIN, TAN, TAU, TEMP, U(*), V(*), W(*), ZERO
  LOGICAL SING
  SAVE ONE, P5, P25, ZERO
  DATA ONE,P5,P25,ZERO /1.0D0,5.0D-1,2.5D-1,0.0D0/
!
!     GIANT IS THE LARGEST MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  D1UPDT
  GIANT = D1MACH(2)
!
!     INITIALIZE THE DIAGONAL ELEMENT POINTER.
!
  JJ = (N*(2*M - N + 1))/2 - (M - N)
!
!     MOVE THE NONTRIVIAL PART OF THE LAST COLUMN OF S INTO W.
!
  L = JJ
  DO 10 I = N, M
     W(I) = S(L)
     L = L + 1
   10    CONTINUE
!
!     ROTATE THE VECTOR V INTO A MULTIPLE OF THE N-TH UNIT VECTOR
!     IN SUCH A WAY THAT A SPIKE IS INTRODUCED INTO W.
!
  NM1 = N - 1
  if (NM1  <  1) go to 70
  DO 60 NMJ = 1, NM1
     J = N - NMJ
     JJ = JJ - (M - J + 1)
     W(J) = ZERO
     if (V(J)  ==  ZERO) go to 50
!
!        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!        J-TH ELEMENT OF V.
!
     if (ABS(V(N))  >=  ABS(V(J))) go to 20
        COTAN = V(N)/V(J)
        SIN = P5/SQRT(P25+P25*COTAN**2)
        COS = SIN*COTAN
        TAU = ONE
        if (ABS(COS)*GIANT  >  ONE) TAU = ONE/COS
        go to 30
   20    CONTINUE
        TAN = V(J)/V(N)
        COS = P5/SQRT(P25+P25*TAN**2)
        SIN = COS*TAN
        TAU = SIN
   30    CONTINUE
!
!        APPLY THE TRANSFORMATION TO V AND STORE THE INFORMATION
!        NECESSARY TO RECOVER THE GIVENS ROTATION.
!
     V(N) = SIN*V(J) + COS*V(N)
     V(J) = TAU
!
!        APPLY THE TRANSFORMATION TO S AND EXTEND THE SPIKE IN W.
!
     L = JJ
     DO 40 I = J, M
        TEMP = COS*S(L) - SIN*W(I)
        W(I) = SIN*S(L) + COS*W(I)
        S(L) = TEMP
        L = L + 1
   40       CONTINUE
   50    CONTINUE
   60    CONTINUE
   70 CONTINUE
!
!     ADD THE SPIKE FROM THE RANK 1 UPDATE TO W.
!
  DO 80 I = 1, M
     W(I) = W(I) + V(N)*U(I)
   80    CONTINUE
!
!     ELIMINATE THE SPIKE.
!
  SING = .FALSE.
  if (NM1  <  1) go to 140
  DO 130 J = 1, NM1
     if (W(J)  ==  ZERO) go to 120
!
!        DETERMINE A GIVENS ROTATION WHICH ELIMINATES THE
!        J-TH ELEMENT OF THE SPIKE.
!
     if (ABS(S(JJ))  >=  ABS(W(J))) go to 90
        COTAN = S(JJ)/W(J)
        SIN = P5/SQRT(P25+P25*COTAN**2)
        COS = SIN*COTAN
        TAU = ONE
        if (ABS(COS)*GIANT  >  ONE) TAU = ONE/COS
        go to 100
   90    CONTINUE
        TAN = W(J)/S(JJ)
        COS = P5/SQRT(P25+P25*TAN**2)
        SIN = COS*TAN
        TAU = SIN
  100    CONTINUE
!
!        APPLY THE TRANSFORMATION TO S AND REDUCE THE SPIKE IN W.
!
     L = JJ
     DO 110 I = J, M
        TEMP = COS*S(L) + SIN*W(I)
        W(I) = -SIN*S(L) + COS*W(I)
        S(L) = TEMP
        L = L + 1
  110       CONTINUE
!
!        STORE THE INFORMATION NECESSARY TO RECOVER THE
!        GIVENS ROTATION.
!
     W(J) = TAU
  120    CONTINUE
!
!        TEST FOR ZERO DIAGONAL ELEMENTS IN THE OUTPUT S.
!
     if (S(JJ)  ==  ZERO) SING = .TRUE.
     JJ = JJ + (M - J + 1)
  130    CONTINUE
  140 CONTINUE
!
!     MOVE W BACK INTO THE LAST COLUMN OF THE OUTPUT S.
!
  L = JJ
  DO 150 I = N, M
     S(L) = W(I)
     L = L + 1
  150    CONTINUE
  if (S(JJ)  ==  ZERO) SING = .TRUE.
  return
!
!     LAST CARD OF SUBROUTINE D1UPDT.
!
end
