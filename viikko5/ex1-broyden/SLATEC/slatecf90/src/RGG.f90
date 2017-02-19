subroutine RGG (NM, N, A, B, ALFR, ALFI, BETA, MATZ, Z, IERR)
!
!! RGG computes eigenvalues and eigenvectors for real generalized eigenproblem.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4B2
!***TYPE      SINGLE PRECISION (RGG-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (EISPACK)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the REAL GENERAL GENERALIZED eigenproblem  Ax = (LAMBDA)Bx.
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
!        A contains a real general matrix.  A is a two-dimensional
!          REAL array, dimensioned A(NM,N).
!
!        B contains a real general matrix.  B is a two-dimensional
!          REAL array, dimensioned B(NM,N).
!
!        MATZ is an INTEGER variable set equal to zero if only
!          eigenvalues are desired.  Otherwise, it is set to any
!          non-zero integer for both eigenvalues and eigenvectors.
!
!     On Output
!
!        A and B have been destroyed.
!
!        ALFR and ALFI contain the real and imaginary parts,
!          respectively, of the numerators of the eigenvalues.
!          ALFR and ALFI are one-dimensional REAL arrays,
!          dimensioned ALFR(N) and ALFI(N).
!
!        BETA contains the denominators of the eigenvalues,
!          which are thus given by the ratios  (ALFR+I*ALFI)/BETA.
!          Complex conjugate pairs of eigenvalues appear consecutively
!          with the eigenvalue having the positive imaginary part first.
!          BETA is a one-dimensional REAL array, dimensioned BETA(N).
!
!        Z contains the real and imaginary parts of the eigenvectors
!          if MATZ is not zero.  If the J-th eigenvalue is real, the
!          J-th column of  Z  contains its eigenvector.  If the J-th
!          eigenvalue is complex with positive imaginary part, the
!          J-th and (J+1)-th columns of  Z  contain the real and
!          imaginary parts of its eigenvector.  The conjugate of this
!          vector is the eigenvector for the conjugate eigenvalue.
!          Z is a two-dimensional REAL array, dimensioned Z(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          10*N       if N is greater than NM,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30*N iterations.
!                     The eigenvalues should be correct for indices
!                     IERR+1, IERR+2, ..., N, but no eigenvectors are
!                     computed.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  QZHES, QZIT, QZVAL, QZVEC
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RGG
!
  INTEGER N,NM,IERR,MATZ
  REAL A(NM,*),B(NM,*),ALFR(*),ALFI(*),BETA(*),Z(NM,*)
  LOGICAL TF
!
!***FIRST EXECUTABLE STATEMENT  RGG
  if (N  <=  NM) go to 10
  IERR = 10 * N
  go to 50
!
   10 if (MATZ  /=  0) go to 20
!     .......... FIND EIGENVALUES ONLY ..........
  TF = .FALSE.
  call  QZHES(NM,N,A,B,TF,Z)
  call  QZIT(NM,N,A,B,0.0E0,TF,Z,IERR)
  call  QZVAL(NM,N,A,B,ALFR,ALFI,BETA,TF,Z)
  go to 50
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 TF = .TRUE.
  call  QZHES(NM,N,A,B,TF,Z)
  call  QZIT(NM,N,A,B,0.0E0,TF,Z,IERR)
  call  QZVAL(NM,N,A,B,ALFR,ALFI,BETA,TF,Z)
  if (IERR  /=  0) go to 50
  call  QZVEC(NM,N,A,B,ALFR,ALFI,BETA,Z)
   50 RETURN
end
