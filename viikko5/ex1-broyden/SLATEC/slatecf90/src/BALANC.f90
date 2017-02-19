subroutine BALANC (NM, N, A, LOW, IGH, SCALE)
!
!! BALANC balances a real general matrix and isolates eigenvalues when possible.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1A
!***TYPE      SINGLE PRECISION (BALANC-S, CBAL-C)
!***KEYWORDS  EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure BALANCE,
!     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
!     HANDBOOK FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971).
!
!     This subroutine balances a REAL matrix and isolates
!     eigenvalues whenever possible.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, A, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        A contains the input matrix to be balanced.  A is a
!          two-dimensional REAL array, dimensioned A(NM,N).
!
!     On OUTPUT
!
!        A contains the balanced matrix.
!
!        LOW and IGH are two INTEGER variables such that A(I,J)
!          is equal to zero if
!           (1) I is greater than J and
!           (2) J=1,...,LOW-1 or I=IGH+1,...,N.
!
!        SCALE contains information determining the permutations and
!          scaling factors used.  SCALE is a one-dimensional REAL array,
!          dimensioned SCALE(N).
!
!     Suppose that the principal submatrix in rows LOW through IGH
!     has been balanced, that P(J) denotes the index interchanged
!     with J during the permutation step, and that the elements
!     of the diagonal matrix used are denoted by D(I,J).  Then
!        SCALE(J) = P(J),    for J = 1,...,LOW-1
!                 = D(J,J),      J = LOW,...,IGH
!                 = P(J)         J = IGH+1,...,N.
!     The order in which the interchanges are made is N to IGH+1,
!     then 1 TO LOW-1.
!
!     Note that 1 is returned for IGH if IGH is zero formally.
!
!     The ALGOL procedure EXC contained in BALANCE appears in
!     BALANC  in line.  (Note that the ALGOL roles of identifiers
!     K,L have been reversed.)
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
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  BALANC
!
  INTEGER I,J,K,L,M,N,JJ,NM,IGH,LOW,IEXC
  REAL A(NM,*),SCALE(*)
  REAL C,F,G,R,S,B2,RADIX
  LOGICAL NOCONV
!
!***FIRST EXECUTABLE STATEMENT  BALANC
  RADIX = 16
!
  B2 = RADIX * RADIX
  K = 1
  L = N
  go to 100
!     .......... IN-LINE PROCEDURE FOR ROW AND
!                COLUMN EXCHANGE ..........
   20 SCALE(M) = J
  if (J  ==  M) go to 50
!
  DO 30 I = 1, L
     F = A(I,J)
     A(I,J) = A(I,M)
     A(I,M) = F
   30 CONTINUE
!
  DO 40 I = K, N
     F = A(J,I)
     A(J,I) = A(M,I)
     A(M,I) = F
   40 CONTINUE
!
   50 go to (80,130), IEXC
!     .......... SEARCH FOR ROWS ISOLATING AN EIGENVALUE
!                AND PUSH THEM DOWN ..........
   80 if (L  ==  1) go to 280
  L = L - 1
!     .......... FOR J=L STEP -1 UNTIL 1 DO -- ..........
  100 DO 120 JJ = 1, L
     J = L + 1 - JJ
!
     DO 110 I = 1, L
        if (I  ==  J) go to 110
        if (A(J,I)  /=  0.0E0) go to 120
  110    CONTINUE
!
     M = L
     IEXC = 1
     go to 20
  120 CONTINUE
!
  go to 140
!     .......... SEARCH FOR COLUMNS ISOLATING AN EIGENVALUE
!                AND PUSH THEM LEFT ..........
  130 K = K + 1
!
  140 DO 170 J = K, L
!
     DO 150 I = K, L
        if (I  ==  J) go to 150
        if (A(I,J)  /=  0.0E0) go to 170
  150    CONTINUE
!
     M = K
     IEXC = 2
     go to 20
  170 CONTINUE
!     .......... NOW BALANCE THE SUBMATRIX IN ROWS K TO L ..........
  DO 180 I = K, L
  180 SCALE(I) = 1.0E0
!     .......... ITERATIVE LOOP FOR NORM REDUCTION ..........
  190 NOCONV = .FALSE.
!
  DO 270 I = K, L
     C = 0.0E0
     R = 0.0E0
!
     DO 200 J = K, L
        if (J  ==  I) go to 200
        C = C + ABS(A(J,I))
        R = R + ABS(A(I,J))
  200    CONTINUE
!     .......... GUARD AGAINST ZERO C OR R DUE TO UNDERFLOW ..........
     if (C  ==  0.0E0 .OR. R  ==  0.0E0) go to 270
     G = R / RADIX
     F = 1.0E0
     S = C + R
  210    if (C  >=  G) go to 220
     F = F * RADIX
     C = C * B2
     go to 210
  220    G = R * RADIX
  230    if (C  <  G) go to 240
     F = F / RADIX
     C = C / B2
     go to 230
!     .......... NOW BALANCE ..........
  240    if ((C + R) / F  >=  0.95E0 * S) go to 270
     G = 1.0E0 / F
     SCALE(I) = SCALE(I) * F
     NOCONV = .TRUE.
!
     DO 250 J = K, N
  250    A(I,J) = A(I,J) * G
!
     DO 260 J = 1, L
  260    A(J,I) = A(J,I) * F
!
  270 CONTINUE
!
  if (NOCONV) go to 190
!
  280 LOW = K
  IGH = L
  return
end
