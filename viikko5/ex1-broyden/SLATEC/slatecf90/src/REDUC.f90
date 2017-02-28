subroutine REDUC (NM, N, A, B, DL, IERR)
!
!! REDUC reduces a generalized symmetric eigenproblem to a standard ...
!            symmetric eigenproblem using Cholesky factorization.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1C
!***TYPE      SINGLE PRECISION (REDUC-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure REDUC1,
!     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
!
!     This subroutine reduces the generalized SYMMETRIC eigenproblem
!     Ax=(LAMBDA)Bx, where B is POSITIVE DEFINITE, to the standard
!     symmetric eigenproblem using the Cholesky factorization of B.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and B, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrices A and B.  If the Cholesky
!          factor L of B is already available, N should be prefixed
!          with a minus sign.  N is an INTEGER variable.
!
!        A and B contain the real symmetric input matrices.  Only
!          the full upper triangles of the matrices need be supplied.
!          If N is negative, the strict lower triangle of B contains,
!          instead, the strict lower triangle of its Cholesky factor L.
!          A and B are two-dimensional REAL arrays, dimensioned A(NM,N)
!          and B(NM,N).
!
!       DL contains, if N is negative, the diagonal elements of L.
!          DL is a one-dimensional REAL array, dimensioned DL(N).
!
!     On Output
!
!        A contains in its full lower triangle the full lower triangle
!          of the symmetric matrix derived from the reduction to the
!          standard form.  The strict upper triangle of A is unaltered.
!
!        B contains in its strict lower triangle the strict lower
!          triangle of its Cholesky factor L.  The full upper triangle
!          of B is unaltered.
!
!        DL contains the diagonal elements of L.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          7*N+1      if B is not positive definite.
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
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  REDUC
!
  INTEGER I,J,K,N,I1,J1,NM,NN,IERR
  REAL A(NM,*),B(NM,*),DL(*)
  REAL X,Y
!
!***FIRST EXECUTABLE STATEMENT  REDUC
  IERR = 0
  NN = ABS(N)
  if (N  <  0) go to 100
!     .......... FORM L IN THE ARRAYS B AND DL ..........
  DO 80 I = 1, N
     I1 = I - 1
!
     DO 80 J = I, N
        X = B(I,J)
        if (I  ==  1) go to 40
!
        DO 20 K = 1, I1
   20       X = X - B(I,K) * B(J,K)
!
   40       if (J  /=  I) go to 60
        if (X  <=  0.0E0) go to 1000
        Y = SQRT(X)
        DL(I) = Y
        go to 80
   60       B(J,I) = X / Y
   80 CONTINUE
!     .......... FORM THE TRANSPOSE OF THE UPPER TRIANGLE OF INV(L)*A
!                IN THE LOWER TRIANGLE OF THE ARRAY A ..........
  100 DO 200 I = 1, NN
     I1 = I - 1
     Y = DL(I)
!
     DO 200 J = I, NN
        X = A(I,J)
        if (I  ==  1) go to 180
!
        DO 160 K = 1, I1
  160       X = X - B(I,K) * A(J,K)
!
  180       A(J,I) = X / Y
  200 CONTINUE
!     .......... PRE-MULTIPLY BY INV(L) AND OVERWRITE ..........
  DO 300 J = 1, NN
     J1 = J - 1
!
     DO 300 I = J, NN
        X = A(I,J)
        if (I  ==  J) go to 240
        I1 = I - 1
!
        DO 220 K = J, I1
  220       X = X - A(K,J) * B(I,K)
!
  240       if (J  ==  1) go to 280
!
        DO 260 K = 1, J1
  260       X = X - A(J,K) * B(I,K)
!
  280       A(I,J) = X / DL(I)
  300 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- B IS NOT POSITIVE DEFINITE ..........
 1000 IERR = 7 * N + 1
 1001 RETURN
end