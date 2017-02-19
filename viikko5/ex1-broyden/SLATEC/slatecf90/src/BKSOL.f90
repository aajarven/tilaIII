subroutine BKSOL (N, A, X)
!
!! BKSOL is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BKSOL-S, DBKSOL-D)
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
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  SDOT
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  BKSOL
!
  DIMENSION A(*),X(*)
!
!***FIRST EXECUTABLE STATEMENT  BKSOL
  M=(N*(N+1))/2
  X(N)=X(N)*A(M)

  NM1=N-1
  DO K=1,NM1
    J=N-K
    M=M-K-1
    X(J)=X(J)*A(M) - SDOT(K,A(M+1),1,X(J+1),1)
  end do

  RETURN
end
