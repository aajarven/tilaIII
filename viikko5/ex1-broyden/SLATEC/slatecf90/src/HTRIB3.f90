subroutine HTRIB3 (NM, N, A, TAU, M, ZR, ZI)
!
!! HTRIB3 computes the eigenvectors of a complex Hermitian matrix from ...
!  the eigenvectors of a real symmetric tridiagonal matrix output from HTRID3.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (HTRIB3-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of a complex analogue of
!     the ALGOL procedure TRBAK3, NUM. MATH. 11, 181-195(1968)
!     by Martin, Reinsch, and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 212-226(1971).
!
!     This subroutine forms the eigenvectors of a COMPLEX HERMITIAN
!     matrix by back transforming those of the corresponding
!     real symmetric tridiagonal matrix determined by  HTRID3.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A, ZR, and ZI, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains some information about the unitary transformations
!          used in the reduction by  HTRID3.  A is a two-dimensional
!          REAL array, dimensioned A(NM,N).
!
!        TAU contains further information about the transformations.
!          TAU is a one-dimensional REAL array, dimensioned TAU(2,N).
!
!        M is the number of eigenvectors to be back transformed.
!          M is an INTEGER variable.
!
!        ZR contains the eigenvectors to be back transformed in its
!          first M columns.  The contents of ZI are immaterial.  ZR and
!          ZI are two-dimensional REAL arrays, dimensioned ZR(NM,M) and
!          ZI(NM,M).
!
!     On OUTPUT
!
!        ZR and ZI contain the real and imaginary parts, respectively,
!          of the transformed eigenvectors in their first M columns.
!
!     NOTE that the last component of each returned vector
!     is real and that vector Euclidean norms are preserved.
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
!***END PROLOGUE  HTRIB3
!
  INTEGER I,J,K,L,M,N,NM
  REAL A(NM,*),TAU(2,*),ZR(NM,*),ZI(NM,*)
  REAL H,S,SI
!
!***FIRST EXECUTABLE STATEMENT  HTRIB3
  if (M  ==  0) go to 200
!     .......... TRANSFORM THE EIGENVECTORS OF THE REAL SYMMETRIC
!                TRIDIAGONAL MATRIX TO THOSE OF THE HERMITIAN
!                TRIDIAGONAL MATRIX. ..........
  DO 50 K = 1, N
!
     DO 50 J = 1, M
        ZI(K,J) = -ZR(K,J) * TAU(2,K)
        ZR(K,J) = ZR(K,J) * TAU(1,K)
   50 CONTINUE
!
  if (N  ==  1) go to 200
!     .......... RECOVER AND APPLY THE HOUSEHOLDER MATRICES ..........
  DO 140 I = 2, N
     L = I - 1
     H = A(I,I)
     if (H  ==  0.0E0) go to 140
!
     DO 130 J = 1, M
        S = 0.0E0
        SI = 0.0E0

        DO K = 1, L
           S = S + A(I,K) * ZR(K,J) - A(K,I) * ZI(K,J)
           SI = SI + A(I,K) * ZI(K,J) + A(K,I) * ZR(K,J)
        end do
!
!  DOUBLE DIVISIONS AVOID POSSIBLE UNDERFLOW.
!
        S = (S / H) / H
        SI = (SI / H) / H

        DO K = 1, L
           ZR(K,J) = ZR(K,J) - S * A(I,K) - SI * A(K,I)
           ZI(K,J) = ZI(K,J) - SI * A(I,K) + S * A(K,I)
        end do

  130    CONTINUE

  140 CONTINUE

  200 RETURN
end
