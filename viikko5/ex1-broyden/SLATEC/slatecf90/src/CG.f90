subroutine CG (NM, N, AR, AI, WR, WI, MATZ, ZR, ZI, FV1, FV2, FV3, &
     IERR)
!
!! CG computes the eigenvalues and, optionally, the eigenvectors ...
!            of a complex general matrix.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A4
!***TYPE      COMPLEX (RG-S, CG-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (EISPACK)
!     to find the eigenvalues and eigenvectors (if desired)
!     of a COMPLEX GENERAL matrix.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, AR, AI, ZR and ZI, as declared in the
!          calling program dimension statement.  NM is an INTEGER
!          variable.
!
!        N is the order of the matrix A=(AR,AI).  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        AR and AI contain the real and imaginary parts, respectively,
!          of the complex general matrix.  AR and AI are two-dimensional
!          REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
!
!        MATZ is an INTEGER variable set equal to zero if only
!          eigenvalues are desired.  Otherwise, it is set to any
!          non-zero integer for both eigenvalues and eigenvectors.
!
!     On OUTPUT
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues.  WR and WI are one-dimensional REAL
!          arrays, dimensioned WR(N) and WI(N).
!
!        ZR and ZI contain the real and imaginary parts, respectively,
!          of the eigenvectors if MATZ is not zero.  ZR and ZI are
!          two-dimensional REAL arrays, dimensioned ZR(NM,N) and
!          ZI(NM,N).
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
!        FV1, FV2, and FV3 are one-dimensional REAL arrays used for
!          temporary storage, dimensioned FV1(N), FV2(N), and FV3(N).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  CBABK2, CBAL, COMQR, COMQR2, CORTH
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CG
!
  INTEGER N,NM,IS1,IS2,IERR,MATZ
  REAL AR(NM,*),AI(NM,*),WR(*),WI(*),ZR(NM,*),ZI(NM,*)
  REAL FV1(*),FV2(*),FV3(*)
!
!***FIRST EXECUTABLE STATEMENT  CG
  if (N  <=  NM) go to 10
  IERR = 10 * N
  go to 50
!
   10 call  CBAL(NM,N,AR,AI,IS1,IS2,FV1)
  call  CORTH(NM,N,IS1,IS2,AR,AI,FV2,FV3)
  if (MATZ  /=  0) go to 20
!     .......... FIND EIGENVALUES ONLY ..........
  call  COMQR(NM,N,IS1,IS2,AR,AI,WR,WI,IERR)
  go to 50
!     .......... FIND BOTH EIGENVALUES AND EIGENVECTORS ..........
   20 call  COMQR2(NM,N,IS1,IS2,FV2,FV3,AR,AI,WR,WI,ZR,ZI,IERR)
  if (IERR  /=  0) go to 50
  call  CBABK2(NM,N,IS1,IS2,FV1,N,ZR,ZI)
   50 RETURN
end
