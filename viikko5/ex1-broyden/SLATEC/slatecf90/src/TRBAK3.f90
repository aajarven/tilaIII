subroutine TRBAK3 (NM, N, NV, A, M, Z)
!
!! TRBAK3 forms the eigenvectors of a real symmetric matrix from the ...
!  eigenvectors of a symmetric tridiagonal matrix formed by TRED3.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (TRBAK3-S)
!***KEYWORDS  EIGENVECTORS OF A REAL SYMMETRIC MATRIX, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TRBAK3,
!     NUM. MATH. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine forms the eigenvectors of a REAL SYMMETRIC
!     matrix by back transforming those of the corresponding
!     symmetric tridiagonal matrix determined by  TRED3.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, Z, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        NV is an INTEGER variable set equal to the dimension of the
!          array A as specified in the calling program.  NV must not
!          be less than  N*(N+1)/2.
!
!        A contains information about the orthogonal transformations
!          used in the reduction by  TRED3  in its first N*(N+1)/2
!          positions.  A is a one-dimensional REAL array, dimensioned
!          A(NV).
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
!     Note that TRBAK3 preserves vector Euclidean norms.
!
!     Questions and comments should be directed to b. s. Garbow,
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
!***END PROLOGUE  TRBAK3
!
  INTEGER I,J,K,L,M,N,IK,IZ,NM,NV
  REAL A(*),Z(NM,*)
  REAL H,S
!
!***FIRST EXECUTABLE STATEMENT  TRBAK3
  if (M  ==  0) go to 200
  if (N  ==  1) go to 200
!
  DO 140 I = 2, N
     L = I - 1
     IZ = (I * L) / 2
     IK = IZ + I
     H = A(IK)
     if (H  ==  0.0E0) go to 140
!
     DO 130 J = 1, M
        S = 0.0E0
        IK = IZ
!
        DO 110 K = 1, L
           IK = IK + 1
           S = S + A(IK) * Z(K,J)
  110       CONTINUE
!     .......... DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
        S = (S / H) / H
        IK = IZ
!
        DO 120 K = 1, L
           IK = IK + 1
           Z(K,J) = Z(K,J) - S * A(IK)
  120       CONTINUE
!
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
end
