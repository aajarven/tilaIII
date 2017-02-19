subroutine BANDV (NM, N, MBW, A, E21, M, W, Z, IERR, NV, RV, RV6)
!
!! BANDV forms the eigenvectors of a real symmetric band matrix ...
!  associated with a set of ordered approximate eigenvalues
!  by inverse iteration.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C3
!***TYPE      SINGLE PRECISION (BANDV-S)
!***KEYWORDS  EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine finds those eigenvectors of a REAL SYMMETRIC
!     BAND matrix corresponding to specified eigenvalues, using inverse
!     iteration.  The subroutine may also be used to solve systems
!     of linear equations with a symmetric or non-symmetric band
!     coefficient matrix.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        MBW is the number of columns of the array A used to store the
!          band matrix.  If the matrix is symmetric, MBW is its (half)
!          band width, denoted MB and defined as the number of adjacent
!          diagonals, including the principal diagonal, required to
!          specify the non-zero portion of the lower triangle of the
!          matrix.  If the subroutine is being used to solve systems
!          of linear equations and the coefficient matrix is not
!          symmetric, it must however have the same number of adjacent
!          diagonals above the main diagonal as below, and in this
!          case, MBW=2*MB-1.  MBW is an INTEGER variable.  MB must not
!          be greater than N.
!
!        A contains the lower triangle of the symmetric band input
!          matrix stored as an N by MB array.  Its lowest subdiagonal
!          is stored in the last N+1-MB positions of the first column,
!          its next subdiagonal in the last N+2-MB positions of the
!          second column, further subdiagonals similarly, and finally
!          its principal diagonal in the N positions of column MB.
!          If the subroutine is being used to solve systems of linear
!          equations and the coefficient matrix is not symmetric, A is
!          N by 2*MB-1 instead with lower triangle as above and with
!          its first superdiagonal stored in the first N-1 positions of
!          column MB+1, its second superdiagonal in the first N-2
!          positions of column MB+2, further superdiagonals similarly,
!          and finally its highest superdiagonal in the first N+1-MB
!          positions of the last column.  Contents of storage locations
!          not part of the matrix are arbitrary.  A is a two-dimensional
!          REAL array, dimensioned A(NM,MBW).
!
!        E21 specifies the ordering of the eigenvalues and contains
!            0.0E0 if the eigenvalues are in ascending order, or
!            2.0E0 if the eigenvalues are in descending order.
!          If the subroutine is being used to solve systems of linear
!          equations, E21 should be set to 1.0E0 if the coefficient
!          matrix is symmetric and to -1.0E0 if not.  E21 is a REAL
!          variable.
!
!        M is the number of specified eigenvalues or the number of
!          systems of linear equations.  M is an INTEGER variable.
!
!        W contains the M eigenvalues in ascending or descending order.
!          If the subroutine is being used to solve systems of linear
!          equations (A-W(J)*I)*X(J)=B(J), where I is the identity
!          matrix, W(J) should be set accordingly, for J=1,2,...,M.
!          W is a one-dimensional REAL array, dimensioned W(M).
!
!        Z contains the constant matrix columns (B(J),J=1,2,...,M), if
!          the subroutine is used to solve systems of linear equations.
!          Z is a two-dimensional REAL array, dimensioned Z(NM,M).
!
!        NV must be set to the dimension of the array parameter RV
!          as declared in the calling program dimension statement.
!          NV is an INTEGER variable.
!
!     On OUTPUT
!
!        A and W are unaltered.
!
!        Z contains the associated set of orthogonal eigenvectors.
!          Any vector which fails to converge is set to zero.  If the
!          subroutine is used to solve systems of linear equations,
!          Z contains the solution matrix columns (X(J),J=1,2,...,M).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          -J         if the eigenvector corresponding to the J-th
!                     eigenvalue fails to converge, or if the J-th
!                     system of linear equations is nearly singular.
!
!        RV and RV6 are temporary storage arrays.  If the subroutine
!          is being used to solve systems of linear equations, the
!          determinant (up to sign) of A-W(M)*I is available, upon
!          return, as the product of the first N elements of RV.
!          RV and RV6 are one-dimensional REAL arrays.  Note that RV
!          is dimensioned RV(NV), where NV must be at least N*(2*MB-1).
!          RV6 is dimensioned RV6(N).
!
!     Questions and comments should be directed to B. S. Garbow,
!     Applied Mathematics Division, ARGONNE NATIONAL LABORATORY
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
!***END PROLOGUE  BANDV
!
  INTEGER I,J,K,M,N,R,II,IJ,JJ,KJ,MB,M1,NM,NV,IJ1,ITS,KJ1,MBW,M21
  INTEGER IERR,MAXJ,MAXK,GROUP
  REAL A(NM,*),W(*),Z(NM,*),RV(*),RV6(*)
  REAL U,V,UK,XU,X0,X1,E21,EPS2,EPS3,EPS4,NORM,ORDER,S
!
!***FIRST EXECUTABLE STATEMENT  BANDV
  IERR = 0
  if (M  ==  0) go to 1001
  MB = MBW
  if (E21  <  0.0E0) MB = (MBW + 1) / 2
  M1 = MB - 1
  M21 = M1 + MB
  ORDER = 1.0E0 - ABS(E21)
!     .......... FIND VECTORS BY INVERSE ITERATION ..........
  DO 920 R = 1, M
     ITS = 1
     X1 = W(R)
     if (R  /=  1) go to 100
!     .......... COMPUTE NORM OF MATRIX ..........
     NORM = 0.0E0
!
     DO 60 J = 1, MB
        JJ = MB + 1 - J
        KJ = JJ + M1
        IJ = 1
        S = 0.0E0
!
        DO 40 I = JJ, N
           S = S + ABS(A(I,J))
           if (E21  >=  0.0E0) go to 40
           S = S + ABS(A(IJ,KJ))
           IJ = IJ + 1
   40       CONTINUE
!
        NORM = MAX(NORM,S)
   60    CONTINUE
!
     if (E21  <  0.0E0) NORM = 0.5E0 * NORM
!     .......... EPS2 IS THE CRITERION FOR GROUPING,
!                EPS3 REPLACES ZERO PIVOTS AND EQUAL
!                ROOTS ARE MODIFIED BY EPS3,
!                EPS4 IS TAKEN VERY SMALL TO AVOID OVERFLOW ..........
     if (NORM  ==  0.0E0) NORM = 1.0E0
     EPS2 = 1.0E-3 * NORM * ABS(ORDER)
     EPS3 = NORM
   70    EPS3 = 0.5E0*EPS3
     if (NORM + EPS3  >  NORM) go to 70
     UK = SQRT(REAL(N))
     EPS3 = UK * EPS3
     EPS4 = UK * EPS3
   80    GROUP = 0
     go to 120
!     .......... LOOK FOR CLOSE OR COINCIDENT ROOTS ..........
  100    if (ABS(X1-X0)  >=  EPS2) go to 80
     GROUP = GROUP + 1
     if (ORDER * (X1 - X0)  <=  0.0E0) X1 = X0 + ORDER * EPS3
!     .......... EXPAND MATRIX, SUBTRACT EIGENVALUE,
!                AND INITIALIZE VECTOR ..........
  120    DO 200 I = 1, N
        IJ = I + MIN(0,I-M1) * N
        KJ = IJ + MB * N
        IJ1 = KJ + M1 * N
        if (M1  ==  0) go to 180
!
        DO 150 J = 1, M1
           if (IJ  >  M1) go to 125
           if (IJ  >  0) go to 130
           RV(IJ1) = 0.0E0
           IJ1 = IJ1 + N
           go to 130
  125          RV(IJ) = A(I,J)
  130          IJ = IJ + N
           II = I + J
           if (II  >  N) go to 150
           JJ = MB - J
           if (E21  >=  0.0E0) go to 140
           II = I
           JJ = MB + J
  140          RV(KJ) = A(II,JJ)
           KJ = KJ + N
  150       CONTINUE
!
  180       RV(IJ) = A(I,MB) - X1
        RV6(I) = EPS4
        if (ORDER  ==  0.0E0) RV6(I) = Z(I,R)
  200    CONTINUE
!
     if (M1  ==  0) go to 600
!     .......... ELIMINATION WITH INTERCHANGES ..........
     DO 580 I = 1, N
        II = I + 1
        MAXK = MIN(I+M1-1,N)
        MAXJ = MIN(N-I,M21-2) * N
!
        DO 360 K = I, MAXK
           KJ1 = K
           J = KJ1 + N
           JJ = J + MAXJ
!
           DO 340 KJ = J, JJ, N
              RV(KJ1) = RV(KJ)
              KJ1 = KJ
  340          CONTINUE
!
           RV(KJ1) = 0.0E0
  360       CONTINUE
!
        if (I  ==  N) go to 580
        U = 0.0E0
        MAXK = MIN(I+M1,N)
        MAXJ = MIN(N-II,M21-2) * N
!
        DO 450 J = I, MAXK
           if (ABS(RV(J))  <  ABS(U)) go to 450
           U = RV(J)
           K = J
  450       CONTINUE
!
        J = I + N
        JJ = J + MAXJ
        if (K  ==  I) go to 520
        KJ = K
!
        DO 500 IJ = I, JJ, N
           V = RV(IJ)
           RV(IJ) = RV(KJ)
           RV(KJ) = V
           KJ = KJ + N
  500       CONTINUE
!
        if (ORDER  /=  0.0E0) go to 520
        V = RV6(I)
        RV6(I) = RV6(K)
        RV6(K) = V
  520       if (U  ==  0.0E0) go to 580
!
        DO 560 K = II, MAXK
           V = RV(K) / U
           KJ = K
!
           DO 540 IJ = J, JJ, N
              KJ = KJ + N
              RV(KJ) = RV(KJ) - V * RV(IJ)
  540          CONTINUE
!
           if (ORDER  ==  0.0E0) RV6(K) = RV6(K) - V * RV6(I)
  560       CONTINUE
!
  580    CONTINUE
!     .......... BACK SUBSTITUTION
!                FOR I=N STEP -1 UNTIL 1 DO -- ..........
  600    DO 630 II = 1, N
        I = N + 1 - II
        MAXJ = MIN(II,M21)
        if (MAXJ  ==  1) go to 620
        IJ1 = I
        J = IJ1 + N
        JJ = J + (MAXJ - 2) * N
!
        DO 610 IJ = J, JJ, N
           IJ1 = IJ1 + 1
           RV6(I) = RV6(I) - RV(IJ) * RV6(IJ1)
  610       CONTINUE
!
  620       V = RV(I)
        if (ABS(V)  >=  EPS3) go to 625
!     .......... SET ERROR -- NEARLY SINGULAR LINEAR SYSTEM ..........
        if (ORDER  ==  0.0E0) IERR = -R
        V = SIGN(EPS3,V)
  625       RV6(I) = RV6(I) / V
  630    CONTINUE
!
     XU = 1.0E0
     if (ORDER  ==  0.0E0) go to 870
!     .......... ORTHOGONALIZE WITH RESPECT TO PREVIOUS
!                MEMBERS OF GROUP ..........
     if (GROUP  ==  0) go to 700
!
     DO 680 JJ = 1, GROUP
        J = R - GROUP - 1 + JJ
        XU = 0.0E0
!
        DO 640 I = 1, N
  640       XU = XU + RV6(I) * Z(I,J)
!
        DO 660 I = 1, N
  660       RV6(I) = RV6(I) - XU * Z(I,J)
!
  680    CONTINUE
!
  700    NORM = 0.0E0
!
     DO 720 I = 1, N
  720    NORM = NORM + ABS(RV6(I))
!
     if (NORM  >=  0.1E0) go to 840
!     .......... IN-LINE PROCEDURE FOR CHOOSING
!                A NEW STARTING VECTOR ..........
     if (ITS  >=  N) go to 830
     ITS = ITS + 1
     XU = EPS4 / (UK + 1.0E0)
     RV6(1) = EPS4
!
     DO 760 I = 2, N
  760    RV6(I) = XU
!
     RV6(ITS) = RV6(ITS) - EPS4 * UK
     go to 600
!     .......... SET ERROR -- NON-CONVERGED EIGENVECTOR ..........
  830    IERR = -R
     XU = 0.0E0
     go to 870
!     .......... NORMALIZE SO THAT SUM OF SQUARES IS
!                1 AND EXPAND TO FULL ORDER ..........
  840    U = 0.0E0
!
     DO 860 I = 1, N
  860    U = U + RV6(I)**2
!
     XU = 1.0E0 / SQRT(U)
!
  870    DO 900 I = 1, N
  900    Z(I,R) = RV6(I) * XU
!
     X0 = X1
  920 CONTINUE
!
 1001 RETURN
end
