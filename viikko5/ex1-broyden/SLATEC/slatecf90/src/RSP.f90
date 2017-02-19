subroutine RSP (NM, N, NV, A, W, MATZ, Z, FV1, FV2, IERR)
!
!! RSP eigenvalues and eigenvectors of real symmetric packed matrix.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A1
!***TYPE      SINGLE PRECISION (RSP-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (EISPACK)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a REAL SYMMETRIC PACKED matrix.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, Z, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        NV is an INTEGER variable set equal to the dimension of the
!          array A as specified in the calling program.  NV must not
!          be less than  N*(N+1)/2.
!
!        A contains the lower triangle, stored row-wise, of the real
!          symmetric packed matrix.  A is a one-dimensional REAL
!          array, dimensioned A(NV).
!
!        MATZ is an INTEGER variable set equal to zero if only
!          eigenvalues are desired.  Otherwise, it is set to any
!          non-zero integer for both eigenvalues and eigenvectors.
!
!     On Output
!
!        A has been destroyed.
!
!        W contains the eigenvalues in ascending order.  W is a
!          one-dimensional REAL array, dimensioned W(N).
!
!        Z contains the eigenvectors if MATZ is not zero.  The eigen-
!          vectors are orthonormal.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          10*N       if N is greater than NM,
!          20*N       if NV is less than N*(N+1)/2,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!                     The eigenvalues and eigenvectors in the W and Z
!                     arrays should be correct for indices
!                     1, 2, ..., IERR-1.
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
!***ROUTINES CALLED  TQL2, TQLRAT, TRBAK3, TRED3
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RSP
!
  INTEGER I,J,N,NM,NV,IERR,MATZ
  REAL A(*),W(*),Z(NM,*),FV1(*),FV2(*)
!
!***FIRST EXECUTABLE STATEMENT  RSP
  if (N  <=  NM) go to 5
  IERR = 10 * N
  go to 50
    5 if (NV  >=  (N * (N + 1)) / 2) go to 10
  IERR = 20 * N
  go to 50
!
   10 call  TRED3(N,NV,A,W,FV1,FV2)
  if (MATZ  /=  0) go to 20
!     .......... FIND EIGENVALUES ONLY ..........
  call  TQLRAT(N,W,FV2,IERR)
  go to 50
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 DO 40 I = 1, N
!
     DO 30 J = 1, N
        Z(J,I) = 0.0E0
   30    CONTINUE
!
     Z(I,I) = 1.0E0
   40 CONTINUE
!
  call  TQL2(NM,N,W,FV1,Z,IERR)
  if (IERR  /=  0) go to 50
  call  TRBAK3(NM,N,NV,A,N,Z)
   50 RETURN
end
