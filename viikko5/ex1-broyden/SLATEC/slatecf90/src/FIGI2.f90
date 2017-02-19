subroutine FIGI2 (NM, N, T, D, E, Z, IERR)
!
!! FIGI2 transforms certain real non-symmetric tridiagonal matrix ...
!            to symmetric tridiagonal matrix.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1C
!***TYPE      SINGLE PRECISION (FIGI2-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     Given a NONSYMMETRIC TRIDIAGONAL matrix such that the products
!     of corresponding pairs of off-diagonal elements are all
!     non-negative, and zero only when both factors are zero, this
!     subroutine reduces it to a SYMMETRIC TRIDIAGONAL matrix
!     using and accumulating diagonal similarity transformations.
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
!     On OUTPUT
!
!        T is unaltered.
!
!        D contains the diagonal elements of the tridiagonal symmetric
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the tridiagonal
!          symmetric matrix in its last N-1 positions.  E(1) is not set.
!          E is a one-dimensional REAL array, dimensioned E(N).
!
!        Z contains the diagonal transformation matrix produced in the
!          symmetrization.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          N+I        if T(I,1)*T(I-1,3) is negative,
!          2*N+I      if T(I,1)*T(I-1,3) is zero with one factor
!                     non-zero.  In these cases, there does not exist
!                     a symmetrizing similarity transformation which
!                     is essential for the validity of the later
!                     eigenvector computation.
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
!***END PROLOGUE  FIGI2
!
  INTEGER I,J,N,NM,IERR
  REAL T(NM,3),D(*),E(*),Z(NM,*)
  REAL H
!
!***FIRST EXECUTABLE STATEMENT  FIGI2
  IERR = 0
!
  DO 100 I = 1, N
!
     DO 50 J = 1, N
   50    Z(I,J) = 0.0E0
!
     if (I  ==  1) go to 70
     H = T(I,1) * T(I-1,3)
     if (H) 900, 60, 80
   60    if (T(I,1)  /=  0.0E0 .OR. T(I-1,3)  /=  0.0E0) go to 1000
     E(I) = 0.0E0
   70    Z(I,I) = 1.0E0
     go to 90
   80    E(I) = SQRT(H)
     Z(I,I) = Z(I-1,I-1) * E(I) / T(I-1,3)
   90    D(I) = T(I,2)
  100 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
!                ELEMENTS IS NEGATIVE ..........
  900 IERR = N + I
  go to 1001
!     .......... SET ERROR -- PRODUCT OF SOME PAIR OF OFF-DIAGONAL
!                ELEMENTS IS ZERO WITH ONE MEMBER NON-ZERO ..........
 1000 IERR = 2 * N + I
 1001 RETURN
end
