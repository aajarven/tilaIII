subroutine RG (NM, N, A, WR, WI, MATZ, Z, IV1, FV1, IERR)
!
!! RG computes the eigenvalues and eigenvectors of a real general matrix.
!!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A2
!***TYPE      SINGLE PRECISION (RG-S, CG-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (EISPACK)
!     To find the eigenvalues and eigenvectors (if desired)
!     of a REAL GENERAL matrix.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains the real general matrix.  A is a two-dimensional
!          REAL array, dimensioned A(NM,N).
!
!        MATZ is an INTEGER variable set equal to zero if only
!          eigenvalues are desired.  Otherwise, it is set to any
!          non-zero integer for both eigenvalues and eigenvectors.
!
!     On Output
!
!        A has been destroyed.
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues.  The eigenvalues are unordered except
!          that complex conjugate pairs of eigenvalues appear consecu-
!          tively with the eigenvalue having the positive imaginary part
!          first.  If an error exit is made, the eigenvalues should be
!          correct for indices IERR+1, IERR+2, ..., N.  WR and WI are
!          one-dimensional REAL arrays, dimensioned WR(N) and WI(N).
!
!        Z contains the real and imaginary parts of the eigenvectors
!          if MATZ is not zero.  If the J-th eigenvalue is real, the
!          J-th column of Z contains its eigenvector.  If the J-th
!          eigenvalue is complex with positive imaginary part, the
!          J-th and (J+1)-th columns of Z contain the real and
!          imaginary parts of its eigenvector.  The conjugate of this
!          vector is the eigenvector for the conjugate eigenvalue.
!          Z is a two-dimensional REAL array, dimensioned Z(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          10*N       if N is greater than NM,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30 iterations.
!                     The eigenvalues should be correct for indices
!                     IERR+1, IERR+2, ..., N, but no eigenvectors are
!                     computed.
!
!        IV1 and FV1 are one-dimensional temporary storage arrays of
!          dimension N.  IV1 is of type INTEGER and FV1 of type REAL.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  BALANC, BALBAK, ELMHES, ELTRAN, HQR, HQR2
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   921103  Corrected description of IV1.  (DWL, FNF and WRB)
!***END PROLOGUE  RG
!
  INTEGER N,NM,IS1,IS2,IERR,MATZ
  REAL A(NM,*),WR(*),WI(*),Z(NM,*),FV1(*)
  INTEGER IV1(*)
!
!***FIRST EXECUTABLE STATEMENT  RG
  if (N  <=  NM) go to 10
  IERR = 10 * N
  go to 50
!
   10 call  BALANC(NM,N,A,IS1,IS2,FV1)
  call  ELMHES(NM,N,IS1,IS2,A,IV1)
  if (MATZ  /=  0) go to 20
!     .......... FIND EIGENVALUES ONLY ..........
  call  HQR(NM,N,IS1,IS2,A,WR,WI,IERR)
  go to 50
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 call  ELTRAN(NM,N,IS1,IS2,A,IV1,Z)
  call  HQR2(NM,N,IS1,IS2,A,WR,WI,Z,IERR)
  if (IERR  /=  0) go to 50
  call  BALBAK(NM,N,IS1,IS2,FV1,N,Z)
   50 RETURN
end
