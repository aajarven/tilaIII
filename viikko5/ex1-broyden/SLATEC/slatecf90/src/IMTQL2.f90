subroutine IMTQL2 (NM, N, D, E, Z, IERR)
!
!! IMTQL2 computes eigenvalues and eigenvectors of a symmetric tridiagonal ...
!  matrix using the implicit QL method.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (IMTQL2-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure IMTQL2,
!     NUM. MATH. 12, 377-383(1968) by Martin and Wilkinson,
!     as modified in NUM. MATH. 15, 450(1970) by Dubrulle.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 241-248(1971).
!
!     This subroutine finds the eigenvalues and eigenvectors
!     of a SYMMETRIC TRIDIAGONAL matrix by the implicit QL method.
!     The eigenvectors of a FULL SYMMETRIC matrix can also
!     be found if  TRED2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, Z, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is
!          arbitrary.  E is a one-dimensional REAL array, dimensioned
!          E(N).
!
!        Z contains the transformation matrix produced in the reduction
!          by  TRED2,  if performed.  This transformation matrix is
!          necessary if you want to obtain the eigenvectors of the full
!          symmetric matrix.  If the eigenvectors of the symmetric
!          tridiagonal matrix are desired, Z must contain the identity
!          matrix.  Z is a two-dimensional REAL array, dimensioned
!          Z(NM,N).
!
!      On OUTPUT
!
!        D contains the eigenvalues in ascending order.  If an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1, 2, ..., IERR-1.
!
!        E has been destroyed.
!
!        Z contains orthonormal eigenvectors of the full symmetric
!          or symmetric tridiagonal matrix, depending on what it
!          contained on input.  If an error exit is made,  Z contains
!          the eigenvectors associated with the stored eigenvalues.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!                     The eigenvalues and eigenvectors should be correct
!                     for indices 1, 2, ..., IERR-1, but the eigenvalues
!                     are not ordered.
!
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  IMTQL2
!
  INTEGER I,J,K,L,M,N,II,NM,MML,IERR
  REAL D(*),E(*),Z(NM,*)
  REAL B,C,F,G,P,R,S,S1,S2
  REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  IMTQL2
  IERR = 0
  if (N  ==  1) go to 1001
!
  DO 100 I = 2, N
  100 E(I-1) = E(I)
!
  E(N) = 0.0E0
!
  DO 240 L = 1, N
     J = 0
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
  105    DO 110 M = L, N
        if (M  ==  N) go to 120
        S1 = ABS(D(M)) + ABS(D(M+1))
        S2 = S1 + ABS(E(M))
        if (S2  ==  S1) go to 120
  110    CONTINUE
!
  120    P = D(L)
     if (M  ==  L) go to 240
     if (J  ==  30) go to 1000
     J = J + 1
!     .......... FORM SHIFT ..........
     G = (D(L+1) - P) / (2.0E0 * E(L))
     R = PYTHAG(G,1.0E0)
     G = D(M) - P + E(L) / (G + SIGN(R,G))
     S = 1.0E0
     C = 1.0E0
     P = 0.0E0
     MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
     DO 200 II = 1, MML
        I = M - II
        F = S * E(I)
        B = C * E(I)
        if (ABS(F)  <  ABS(G)) go to 150
        C = G / F
        R = SQRT(C*C+1.0E0)
        E(I+1) = F * R
        S = 1.0E0 / R
        C = C * S
        go to 160
  150       S = F / G
        R = SQRT(S*S+1.0E0)
        E(I+1) = G * R
        C = 1.0E0 / R
        S = S * C
  160       G = D(I+1) - P
        R = (D(I) - G) * S + 2.0E0 * C * B
        P = S * R
        D(I+1) = G + P
        G = C * R - B
!     .......... FORM VECTOR ..........
        DO 180 K = 1, N
           F = Z(K,I+1)
           Z(K,I+1) = S * Z(K,I) + C * F
           Z(K,I) = C * Z(K,I) - S * F
  180       CONTINUE
!
  200    CONTINUE
!
     D(L) = D(L) - P
     E(L) = G
     E(M) = 0.0E0
     go to 105
  240 CONTINUE
!     .......... ORDER EIGENVALUES AND EIGENVECTORS ..........
  DO 300 II = 2, N
     I = II - 1
     K = I
     P = D(I)
!
     DO 260 J = II, N
        if (D(J)  >=  P) go to 260
        K = J
        P = D(J)
  260    CONTINUE
!
     if (K  ==  I) go to 300
     D(K) = D(I)
     D(I) = P
!
     DO 280 J = 1, N
        P = Z(J,I)
        Z(J,I) = Z(J,K)
        Z(J,K) = P
  280    CONTINUE
!
  300 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS ..........
 1000 IERR = L
 1001 RETURN
end
