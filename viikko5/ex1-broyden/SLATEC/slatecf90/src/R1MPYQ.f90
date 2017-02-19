subroutine R1MPYQ (M, N, A, LDA, V, W)
!
!! R1MPYQ is subsidiary to SNSQ and SNSQE.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (R1MPYQ-S, D1MPYQ-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Given an M by N matrix A, this subroutine computes A*Q where
!     Q is the product of 2*(N - 1) transformations
!
!           GV(N-1)*...*GV(1)*GW(1)*...*GW(N-1)
!
!     and GV(I), GW(I) are Givens rotations in the (I,N) plane which
!     eliminate elements in the I-th and N-th planes, respectively.
!     Q itself is not given, rather the information to recover the
!     GV, GW rotations is supplied.
!
!     The subroutine statement is
!
!       SUBROUTINE R1MPYQ(M,N,A,LDA,V,W)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       A is an M by N ARRAY. On input A must contain the matrix
!         to be postmultiplied by the orthogonal matrix Q
!         described above. On output A*Q has replaced A.
!
!       LDA is a positive integer input variable not less than M
!         which specifies the leading dimension of the array A.
!
!       V is an input array of length N. V(I) must contain the
!         information necessary to recover the Givens rotation GV(I)
!         described above.
!
!       W is an input array of length N. W(I) must contain the
!         information necessary to recover the Givens rotation GW(I)
!         described above.
!
!***SEE ALSO  SNSQ, SNSQE
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  R1MPYQ
  INTEGER M,N,LDA
  REAL A(LDA,*),V(*),W(*)
  INTEGER I,J,NMJ,NM1
  REAL COS,ONE,SIN,TEMP
  SAVE ONE
  DATA ONE /1.0E0/
!***FIRST EXECUTABLE STATEMENT  R1MPYQ
  NM1 = N - 1
  if (NM1  <  1) go to 50
  DO 20 NMJ = 1, NM1
     J = N - NMJ
     if (ABS(V(J))  >  ONE) COS = ONE/V(J)
     if (ABS(V(J))  >  ONE) SIN = SQRT(ONE-COS**2)
     if (ABS(V(J))  <=  ONE) SIN = V(J)
     if (ABS(V(J))  <=  ONE) COS = SQRT(ONE-SIN**2)
     DO 10 I = 1, M
        TEMP = COS*A(I,J) - SIN*A(I,N)
        A(I,N) = SIN*A(I,J) + COS*A(I,N)
        A(I,J) = TEMP
   10       CONTINUE
   20    CONTINUE
!
!     APPLY THE SECOND SET OF GIVENS ROTATIONS TO A.
!
  DO 40 J = 1, NM1
     if (ABS(W(J))  >  ONE) COS = ONE/W(J)
     if (ABS(W(J))  >  ONE) SIN = SQRT(ONE-COS**2)
     if (ABS(W(J))  <=  ONE) SIN = W(J)
     if (ABS(W(J))  <=  ONE) COS = SQRT(ONE-SIN**2)
     DO 30 I = 1, M
        TEMP = COS*A(I,J) + SIN*A(I,N)
        A(I,N) = -SIN*A(I,J) + COS*A(I,N)
        A(I,J) = TEMP
   30       CONTINUE
   40    CONTINUE
   50 CONTINUE
  return
!
!     LAST CARD OF SUBROUTINE R1MPYQ.
!
end
