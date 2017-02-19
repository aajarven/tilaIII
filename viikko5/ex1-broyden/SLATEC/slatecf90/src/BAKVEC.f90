subroutine BAKVEC (NM, N, T, E, M, Z, IERR)
!
!! BAKVEC forms the eigenvectors of a certain real non-symmetric tridiagonal
!  matrix from a symmetric tridiagonal matrix output from FIGI.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (BAKVEC-S)
!***KEYWORDS  EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine forms the eigenvectors of a NONSYMMETRIC
!     TRIDIAGONAL matrix by back transforming those of the
!     corresponding symmetric matrix determined by  FIGI.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, T and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix T.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        T contains the nonsymmetric matrix.  Its subdiagonal is
!          stored in the last N-1 positions of the first column,
!          its diagonal in the N positions of the second column,
!          and its superdiagonal in the first N-1 positions of
!          the third column.  T(1,1) and T(N,3) are arbitrary.
!          T is a two-dimensional REAL array, dimensioned T(NM,3).
!
!        E contains the subdiagonal elements of the symmetric
!          matrix in its last N-1 positions.  E(1) is arbitrary.
!          E is a one-dimensional REAL array, dimensioned E(N).
!
!        M is the number of eigenvectors to be back transformed.
!          M is an INTEGER variable.
!
!        Z contains the eigenvectors to be back transformed
!          in its first M columns.  Z is a two-dimensional REAL
!          array, dimensioned Z(NM,M).
!
!     On OUTPUT
!
!        T is unaltered.
!
!        E is destroyed.
!
!        Z contains the transformed eigenvectors in its first M columns.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          2*N+I      if E(I) is zero with T(I,1) or T(I-1,3) non-zero.
!                     In this case, the symmetric matrix is not similar
!                     to the original matrix, and the eigenvectors
!                     cannot be found by this program.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BAKVEC
!
  INTEGER I,J,M,N,NM,IERR
  REAL T(NM,3),E(*),Z(NM,*)
!
!***FIRST EXECUTABLE STATEMENT  BAKVEC
  IERR = 0
  if (M  ==  0) go to 1001
  E(1) = 1.0E0
  if (N  ==  1) go to 1001
!
  DO 100 I = 2, N
     if (E(I)  /=  0.0E0) go to 80
     if (T(I,1)  /=  0.0E0 .OR. T(I-1,3)  /=  0.0E0) go to 1000
     E(I) = 1.0E0
     go to 100
   80    E(I) = E(I-1) * E(I) / T(I-1,3)
  100 CONTINUE
!
  DO 120 J = 1, M
!
     DO 120 I = 2, N
     Z(I,J) = Z(I,J) * E(I)
  120 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- EIGENVECTORS CANNOT BE
!                FOUND BY THIS PROGRAM ..........
 1000 IERR = 2 * N + I
 1001 RETURN
end
