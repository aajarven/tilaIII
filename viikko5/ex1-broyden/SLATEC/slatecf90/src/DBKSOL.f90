subroutine DBKSOL (N, A, X)
!
!! DBKSOL is subsidiary to DBVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (BKSOL-S, DBKSOL-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!     Solution of an upper triangular linear system by
!     back-substitution
!
!     The matrix A is assumed to be stored in a linear
!     array proceeding in a row-wise manner. The
!     vector X contains the given constant vector on input
!     and contains the solution on return.
!     The actual diagonal of A is unity while a diagonal
!     scaling matrix is stored there.
! **********************************************************************
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DDOT
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DBKSOL
!
  DOUBLE PRECISION DDOT
  INTEGER J, K, M, N, NM1
  DOUBLE PRECISION A(*), X(*)
!
!***FIRST EXECUTABLE STATEMENT  DBKSOL
  M = (N*(N + 1))/2
  X(N) = X(N)*A(M)
  NM1 = N - 1
  if (NM1  <  1) go to 20
  DO 10 K = 1, NM1
     J = N - K
     M = M - K - 1
     X(J) = X(J)*A(M) - DDOT(K,A(M+1),1,X(J+1),1)
   10 CONTINUE
   20 CONTINUE
!
  return
end
