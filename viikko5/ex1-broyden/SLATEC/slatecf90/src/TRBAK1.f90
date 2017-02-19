subroutine TRBAK1 (NM, N, A, E, M, Z)
!
!! TRBAK1 forms the eigenvectors of real symmetric matrix ...
!  from the eigenvectors of a symmetric tridiagonal matrix formed by TRED1.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (TRBAK1-S)
!***KEYWORDS  EIGENVECTORS OF A REAL SYMMETRIC MATRIX, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TRBAK1,
!     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine forms the eigenvectors of a REAL SYMMETRIC
!     matrix by back transforming those of the corresponding
!     symmetric tridiagonal matrix determined by  TRED1.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains information about the orthogonal transformations
!          used in the reduction by  TRED1  in its strict lower
!          triangle.  A is a two-dimensional REAL array, dimensioned
!          A(NM,N).
!
!        E contains the subdiagonal elements of the tridiagonal matrix
!          in its last N-1 positions.  E(1) is arbitrary.  These
!          elements provide the remaining information about the
!          orthogonal transformations.  E is a one-dimensional REAL
!          array, dimensioned E(N).
!
!        M is the number of columns of Z to be back transformed.
!          M is an INTEGER variable.
!
!        Z contains the eigenvectors to be back transformed in its
!          first M columns.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,M).
!
!     On Output
!
!        Z contains the transformed eigenvectors in its first M columns.
!
!     Note that TRBAK1 preserves vector Euclidean norms.
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
!***END PROLOGUE  TRBAK1
!
  INTEGER I,J,K,L,M,N,NM
  REAL A(NM,*),E(*),Z(NM,*)
  REAL S
!
!***FIRST EXECUTABLE STATEMENT  TRBAK1
  if (M  ==  0) go to 200
  if (N  ==  1) go to 200
!
  DO 140 I = 2, N
     L = I - 1
     if (E(I)  ==  0.0E0) go to 140
!
     DO 130 J = 1, M
        S = 0.0E0
!
        DO 110 K = 1, L
  110       S = S + A(I,K) * Z(K,J)
!     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN TRED1.
!                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
        S = (S / A(I,L)) / E(I)
!
        DO 120 K = 1, L
  120       Z(K,J) = Z(K,J) + S * A(I,K)
!
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
end
