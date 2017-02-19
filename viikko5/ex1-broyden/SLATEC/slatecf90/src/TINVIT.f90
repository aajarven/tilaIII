subroutine TINVIT (NM, N, D, E, E2, M, W, IND, Z, IERR, RV1, RV2, &
     RV3, RV4, RV6)
!
!! TINVIT computes eigenvectors of symmetric tridiagonal matrix ...
!            corresponding to specified eigenvalues, using inverse
!            iteration.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C3
!***TYPE      SINGLE PRECISION (TINVIT-S)
!***KEYWORDS  EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the inverse iteration tech-
!     nique in the ALGOL procedure TRISTURM by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 418-439(1971).
!
!     This subroutine finds those eigenvectors of a TRIDIAGONAL
!     SYMMETRIC matrix corresponding to specified eigenvalues,
!     using inverse iteration.
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
!        E2 contains the squares of the corresponding elements of E,
!          with zeros corresponding to negligible elements of E.
!          E(I) is considered negligible if it is not larger than
!          the product of the relative machine precision and the sum
!          of the magnitudes of D(I) and D(I-1).  E2(1) must contain
!          0.0e0 if the eigenvalues are in ascending order, or 2.0e0
!          if the eigenvalues are in descending order.  If  BISECT,
!          TRIDIB, or  IMTQLV  has been used to find the eigenvalues,
!          their output E2 array is exactly what is expected here.
!          E2 is a one-dimensional REAL array, dimensioned E2(N).
!
!        M is the number of specified eigenvalues for which eigenvectors
!          are to be determined.  M is an INTEGER variable.
!
!        W contains the M eigenvalues in ascending or descending order.
!          W is a one-dimensional REAL array, dimensioned W(M).
!
!        IND contains in its first M positions the submatrix indices
!          associated with the corresponding eigenvalues in W --
!          1 for eigenvalues belonging to the first submatrix from
!          the top, 2 for those belonging to the second submatrix, etc.
!          If  BISECT  or  TRIDIB  has been used to determine the
!          eigenvalues, their output IND array is suitable for input
!          to TINVIT.  IND is a one-dimensional INTEGER array,
!          dimensioned IND(M).
!
!     On Output
!
!       ** All input arrays are unaltered.**
!
!        Z contains the associated set of orthonormal eigenvectors.
!          Any vector which fails to converge is set to zero.
!          Z is a two-dimensional REAL array, dimensioned Z(NM,M).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          -J         if the eigenvector corresponding to the J-th
!                     eigenvalue fails to converge in 5 iterations.
!
!        RV1, RV2 and RV3 are one-dimensional REAL arrays used for
!          temporary storage.  They are used to store the main diagonal
!          and the two adjacent diagonals of the triangular matrix
!          produced in the inverse iteration process.  RV1, RV2 and
!          RV3 are dimensioned RV1(N), RV2(N) and RV3(N).
!
!        RV4 and RV6 are one-dimensional REAL arrays used for temporary
!          storage.  RV4 holds the multipliers of the Gaussian
!          elimination process.  RV6 holds the approximate eigenvectors
!          in this process.  RV4 and RV6 are dimensioned RV4(N) and
!          RV6(N).
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
!***END PROLOGUE  TINVIT
!
  INTEGER I,J,M,N,P,Q,R,S,II,IP,JJ,NM,ITS,TAG,IERR,GROUP
  INTEGER IND(*)
  REAL D(*),E(*),E2(*),W(*),Z(NM,*)
  REAL RV1(*),RV2(*),RV3(*),RV4(*),RV6(*)
  REAL U,V,UK,XU,X0,X1,EPS2,EPS3,EPS4,NORM,ORDER
!
!***FIRST EXECUTABLE STATEMENT  TINVIT
  IERR = 0
  if (M  ==  0) go to 1001
  TAG = 0
  ORDER = 1.0E0 - E2(1)
  Q = 0
!     .......... ESTABLISH AND PROCESS NEXT SUBMATRIX ..........
  100 P = Q + 1
!
  DO 120 Q = P, N
     if (Q  ==  N) go to 140
     if (E2(Q+1)  ==  0.0E0) go to 140
  120 CONTINUE
!     .......... FIND VECTORS BY INVERSE ITERATION ..........
  140 TAG = TAG + 1
  S = 0
!
  DO 920 R = 1, M
     if (IND(R)  /=  TAG) go to 920
     ITS = 1
     X1 = W(R)
     if (S  /=  0) go to 510
!     .......... CHECK FOR ISOLATED ROOT ..........
     XU = 1.0E0
     if (P  /=  Q) go to 490
     RV6(P) = 1.0E0
     go to 870
  490    NORM = ABS(D(P))
     IP = P + 1
!
     DO 500 I = IP, Q
  500    NORM = MAX(NORM, ABS(D(I)) + ABS(E(I)))
!     .......... EPS2 IS THE CRITERION FOR GROUPING,
!                EPS3 REPLACES ZERO PIVOTS AND EQUAL
!                ROOTS ARE MODIFIED BY EPS3,
!                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
     EPS2 = 1.0E-3 * NORM
     EPS3 = NORM
  502    EPS3 = 0.5E0*EPS3
     if (NORM + EPS3  >  NORM) go to 502
     UK = SQRT(REAL(Q-P+5))
     EPS3 = UK * EPS3
     EPS4 = UK * EPS3
     UK = EPS4 / UK
     S = P
  505    GROUP = 0
     go to 520
!     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  510    if (ABS(X1-X0)  >=  EPS2) go to 505
     GROUP = GROUP + 1
     if (ORDER * (X1 - X0)  <=  0.0E0) X1 = X0 + ORDER * EPS3
!     .......... ELIMINATION WITH INTERCHANGES AND
!                INITIALIZATION OF VECTOR ..........
  520    V = 0.0E0
!
     DO 580 I = P, Q
        RV6(I) = UK
        if (I  ==  P) go to 560
        if (ABS(E(I))  <  ABS(U)) go to 540
!     .......... WARNING -- A DIVIDE CHECK MAY OCCUR HERE IF
!                E2 ARRAY HAS NOT BEEN SPECIFIED CORRECTLY ..........
        XU = U / E(I)
        RV4(I) = XU
        RV1(I-1) = E(I)
        RV2(I-1) = D(I) - X1
        RV3(I-1) = 0.0E0
        if (I  /=  Q) RV3(I-1) = E(I+1)
        U = V - XU * RV2(I-1)
        V = -XU * RV3(I-1)
        go to 580
  540       XU = E(I) / U
        RV4(I) = XU
        RV1(I-1) = U
        RV2(I-1) = V
        RV3(I-1) = 0.0E0
  560       U = D(I) - X1 - XU * V
        if (I  /=  Q) V = E(I+1)
  580    CONTINUE
!
     if (U  ==  0.0E0) U = EPS3
     RV1(Q) = U
     RV2(Q) = 0.0E0
     RV3(Q) = 0.0E0
!     .......... BACK SUBSTITUTION
!                FOR I=Q STEP -1 UNTIL P DO -- ..........
  600    DO 620 II = P, Q
        I = P + Q - II
        RV6(I) = (RV6(I) - U * RV2(I) - V * RV3(I)) / RV1(I)
        V = U
        U = RV6(I)
  620    CONTINUE
!     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
!                MEMBERS OF GROUP ..........
     if (GROUP  ==  0) go to 700
     J = R
!
     DO 680 JJ = 1, GROUP
  630       J = J - 1
        if (IND(J)  /=  TAG) go to 630
        XU = 0.0E0
!
        DO 640 I = P, Q
  640       XU = XU + RV6(I) * Z(I,J)
!
        DO 660 I = P, Q
  660       RV6(I) = RV6(I) - XU * Z(I,J)
!
  680    CONTINUE
!
  700    NORM = 0.0E0
!
     DO 720 I = P, Q
  720    NORM = NORM + ABS(RV6(I))
!
     if (NORM  >=  1.0E0) go to 840
!     .......... FORWARD SUBSTITUTION ..........
     if (ITS  ==  5) go to 830
     if (NORM  /=  0.0E0) go to 740
     RV6(S) = EPS4
     S = S + 1
     if (S  >  Q) S = P
     go to 780
  740    XU = EPS4 / NORM
!
     DO 760 I = P, Q
  760    RV6(I) = RV6(I) * XU
!     .......... ELIMINATION OPERATIONS ON NEXT VECTOR
!                ITERATE ..........
  780    DO 820 I = IP, Q
        U = RV6(I)
!     .......... if RV1(I-1)  ==  E(I), A ROW INTERCHANGE
!                WAS PERFORMED EARLIER IN THE
!                TRIANGULARIZATION PROCESS ..........
        if (RV1(I-1)  /=  E(I)) go to 800
        U = RV6(I-1)
        RV6(I-1) = RV6(I)
  800       RV6(I) = U - RV4(I) * RV6(I-1)
  820    CONTINUE
!
     ITS = ITS + 1
     go to 600
!     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  830    IERR = -R
     XU = 0.0E0
     go to 870
!     .......... NORMALIZE SO THAT SUM OF SQUARES IS
!                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
!
     DO 860 I = P, Q
  860    U = U + RV6(I)**2
!
     XU = 1.0E0 / SQRT(U)
!
  870    DO 880 I = 1, N
  880    Z(I,R) = 0.0E0
!
     DO 900 I = P, Q
  900    Z(I,R) = RV6(I) * XU
!
     X0 = X1
  920 CONTINUE
!
  if (Q  <  N) go to 100
 1001 RETURN
end
