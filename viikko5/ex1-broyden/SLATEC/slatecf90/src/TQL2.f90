subroutine TQL2 (NM, N, D, E, Z, IERR)
!
!! TQL2 computes eigenvalues and eigenvectors of symmetric tridiagonal matrix.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4A5, D4C2A
!***TYPE      SINGLE PRECISION (TQL2-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure TQL2,
!     NUM. MATH. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
!     Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 227-240(1971).
!
!     This subroutine finds the eigenvalues and eigenvectors
!     of a SYMMETRIC TRIDIAGONAL matrix by the QL method.
!     The eigenvectors of a FULL SYMMETRIC matrix can also
!     be found if  TRED2  has been used to reduce this
!     full matrix to tridiagonal form.
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
!        D contains the diagonal elements of the symmetric tridiagonal
!          matrix.  D is a one-dimensional REAL array, dimensioned D(N).
!
!        E contains the subdiagonal elements of the symmetric
!          tridiagonal matrix in its last N-1 positions.  E(1) is
!          arbitrary.  E is a one-dimensional REAL array, dimensioned
!          E(N).
!
!        Z contains the transformation matrix produced in the
!          reduction by  TRED2, if performed.  If the eigenvectors
!          of the tridiagonal matrix are desired, Z must contain
!          the identity matrix.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,N).
!
!      On Output
!
!        D contains the eigenvalues in ascending order.  If an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1, 2, ..., IERR-1.
!
!        E has been destroyed.
!
!        Z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  If an error exit is made,
!          Z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
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
!***END PROLOGUE  TQL2
!
  INTEGER I,J,K,L,M,N,II,L1,L2,NM,MML,IERR
  REAL D(*),E(*),Z(NM,*)
  REAL B,C,C2,C3,DL1,EL1,F,G,H,P,R,S,S2
  REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  TQL2
  IERR = 0
  if (N  ==  1) go to 1001
!
  DO 100 I = 2, N
  100 E(I-1) = E(I)
!
  F = 0.0E0
  B = 0.0E0
  E(N) = 0.0E0
!
  DO 240 L = 1, N
     J = 0
     H = ABS(D(L)) + ABS(E(L))
     if (B  <  H) B = H
!     .......... LOOK FOR SMALL SUB-DIAGONAL ELEMENT ..........
     DO 110 M = L, N
        if (B + ABS(E(M))  ==  B) go to 120
!     .......... E(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP ..........
  110    CONTINUE
!
  120    if (M  ==  L) go to 220
  130    if (J  ==  30) go to 1000
     J = J + 1
!     .......... FORM SHIFT ..........
     L1 = L + 1
     L2 = L1 + 1
     G = D(L)
     P = (D(L1) - G) / (2.0E0 * E(L))
     R = PYTHAG(P,1.0E0)
     D(L) = E(L) / (P + SIGN(R,P))
     D(L1) = E(L) * (P + SIGN(R,P))
     DL1 = D(L1)
     H = G - D(L)
     if (L2  >  N) go to 145
!
     DO 140 I = L2, N
  140    D(I) = D(I) - H
!
  145    F = F + H
!     .......... QL TRANSFORMATION ..........
     P = D(M)
     C = 1.0E0
     C2 = C
     EL1 = E(L1)
     S = 0.0E0
     MML = M - L
!     .......... FOR I=M-1 STEP -1 UNTIL L DO -- ..........
     DO 200 II = 1, MML
        C3 = C2
        C2 = C
        S2 = S
        I = M - II
        G = C * E(I)
        H = C * P
        if (ABS(P)  <  ABS(E(I))) go to 150
        C = E(I) / P
        R = SQRT(C*C+1.0E0)
        E(I+1) = S * P * R
        S = C / R
        C = 1.0E0 / R
        go to 160
  150       C = P / E(I)
        R = SQRT(C*C+1.0E0)
        E(I+1) = S * E(I) * R
        S = 1.0E0 / R
        C = C * S
  160       P = C * D(I) - S * G
        D(I+1) = H + S * (C * G + S * D(I))
!     .......... FORM VECTOR ..........
        DO 180 K = 1, N
           H = Z(K,I+1)
           Z(K,I+1) = S * Z(K,I) + C * H
           Z(K,I) = C * Z(K,I) - S * H
  180       CONTINUE
!
  200    CONTINUE
!
     P = -S * S2 * C3 * EL1 * E(L) / DL1
     E(L) = S * P
     D(L) = C * P
     if (B + ABS(E(L))  >  B) go to 130
  220    D(L) = D(L) + F
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
