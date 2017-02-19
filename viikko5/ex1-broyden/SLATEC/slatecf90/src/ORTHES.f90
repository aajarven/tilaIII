subroutine ORTHES (NM, N, LOW, IGH, A, ORT)
!
!! ORTHES reduces a real general matrix to upper Hessenberg form ...
!            using orthogonal similarity transformations.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B2
!***TYPE      SINGLE PRECISION (ORTHES-S, CORTH-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure ORTHES,
!     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!     Given a REAL GENERAL matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     LOW through IGH to upper Hessenberg form by
!     orthogonal similarity transformations.
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
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
!
!        A contains the general matrix to be reduced to upper
!          Hessenberg form.  A is a two-dimensional REAL array,
!          dimensioned A(NM,N).
!
!     On OUTPUT
!
!        A contains the upper Hessenberg matrix.  Some information about
!          the orthogonal transformations used in the reduction
!          is stored in the remaining triangle under the Hessenberg
!          matrix.
!
!        ORT contains further information about the orthogonal trans-
!          formations used in the reduction.  Only elements LOW+1
!          through IGH are used.  ORT is a one-dimensional REAL array,
!          dimensioned ORT(IGH).
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
!***END PROLOGUE  ORTHES
!
  INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
  REAL A(NM,*),ORT(*)
  REAL F,G,H,SCALE
!
!***FIRST EXECUTABLE STATEMENT  ORTHES
  LA = IGH - 1
  KP1 = LOW + 1
  if (LA  <  KP1) go to 200
!
  DO 180 M = KP1, LA
     H = 0.0E0
     ORT(M) = 0.0E0
     SCALE = 0.0E0
!     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
     DO 90 I = M, IGH
   90    SCALE = SCALE + ABS(A(I,M-1))
!
     if (SCALE  ==  0.0E0) go to 180
     MP = M + IGH
!     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
     DO 100 II = M, IGH
        I = MP - II
        ORT(I) = A(I,M-1) / SCALE
        H = H + ORT(I) * ORT(I)
  100    CONTINUE
!
     G = -SIGN(SQRT(H),ORT(M))
     H = H - ORT(M) * G
     ORT(M) = ORT(M) - G
!     .......... FORM (I-(U*UT)/H) * A ..........
     DO 130 J = M, N
        F = 0.0E0
!     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
        DO 110 II = M, IGH
           I = MP - II
           F = F + ORT(I) * A(I,J)
  110       CONTINUE
!
        F = F / H
!
        DO 120 I = M, IGH
  120       A(I,J) = A(I,J) - F * ORT(I)
!
  130    CONTINUE
!     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
     DO 160 I = 1, IGH
        F = 0.0E0
!     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
        DO 140 JJ = M, IGH
           J = MP - JJ
           F = F + ORT(J) * A(I,J)
  140       CONTINUE
!
        F = F / H
!
        DO 150 J = M, IGH
  150       A(I,J) = A(I,J) - F * ORT(J)
!
  160    CONTINUE
!
     ORT(M) = SCALE * ORT(M)
     A(M,M-1) = SCALE * G
  180 CONTINUE
!
  200 RETURN
end
