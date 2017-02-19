subroutine RSGBA (NM, N, A, B, W, MATZ, Z, FV1, FV2, IERR)
!
!! RSGBA: eigenvalues and eigenvectors of a symmetric generalized eigenproblem.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4B1
!***TYPE      SINGLE PRECISION (RSGBA-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (EISPACK)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the REAL SYMMETRIC generalized eigenproblem  BAx = (LAMBDA)x.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A, B, and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrices A and B.  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        A contains a real symmetric matrix.  A is a two-dimensional
!          REAL array, dimensioned A(NM,N).
!
!        B contains a positive definite real symmetric matrix.  B is a
!          two-dimensional REAL array, dimensioned B(NM,N).
!
!        MATZ is an INTEGER variable set equal to zero if only
!          eigenvalues are desired.  Otherwise, it is set to any
!          non-zero integer for both eigenvalues and eigenvectors.
!
!     On Output
!
!        W contains the eigenvalues in ascending order.  W is a
!          one-dimensional REAL array, dimensioned W(N).
!
!        Z contains the eigenvectors if MATZ is not zero.  Z is a
!          two-dimensional REAL array, dimensioned Z(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          10*N       if N is greater than NM,
!          7*N+1      if B is not positive definite,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!                     The eigenvalues should be correct for indices
!                     1, 2, ..., IERR-1, but no eigenvectors are
!                     computed.
!
!        FV1 and FV2 are one-dimensional REAL arrays used for temporary
!          storage, dimensioned FV1(N) and FV2(N).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  REBAKB, REDUC2, TQL2, TQLRAT, TRED1, TRED2
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RSGBA
!
  INTEGER N,NM,IERR,MATZ
  REAL A(NM,*),B(NM,*),W(*),Z(NM,*),FV1(*),FV2(*)
!
!***FIRST EXECUTABLE STATEMENT  RSGBA
  if (N  <=  NM) go to 10
  IERR = 10 * N
  go to 50
!
   10 call  REDUC2(NM,N,A,B,FV2,IERR)
  if (IERR  /=  0) go to 50
  if (MATZ  /=  0) go to 20
!     .......... FIND EIGENVALUES ONLY ..........
  call  TRED1(NM,N,A,W,FV1,FV2)
  call  TQLRAT(N,W,FV2,IERR)
  go to 50
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 call  TRED2(NM,N,A,W,FV1,Z)
  call  TQL2(NM,N,W,FV1,Z,IERR)
  if (IERR  /=  0) go to 50
  call  REBAKB(NM,N,B,FV2,N,Z)
   50 RETURN
end
