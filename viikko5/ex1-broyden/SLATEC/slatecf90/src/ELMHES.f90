subroutine ELMHES (NM, N, LOW, IGH, A, INT)
!
!! ELMHES reduces a real general matrix to upper Hessenberg form ...
!  using stabilized elementary similarity transformations.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B2
!***TYPE      SINGLE PRECISION (ELMHES-S, COMHES-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure ELMHES,
!     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!     Given a REAL GENERAL matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     LOW through IGH to upper Hessenberg form by
!     stabilized elementary similarity transformations.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, A, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix, A.  N is an INTEGER variable.
!          N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
!
!        A contains the input matrix.  A is a two-dimensional REAL
!          array, dimensioned A(NM,N).
!
!     On OUTPUT
!
!        A contains the upper Hessenberg matrix.  The multipliers which
!          were used in the reduction are stored in the remaining
!          triangle under the Hessenberg matrix.
!
!        INT contains information on the rows and columns interchanged
!          in the reduction.  Only elements LOW through IGH are used.
!          INT is a one-dimensional INTEGER array, dimensioned INT(IGH).
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
!***END PROLOGUE  ELMHES
!
  INTEGER I,J,M,N,LA,NM,IGH,KP1,LOW,MM1,MP1
  REAL A(NM,*)
  REAL X,Y
  INTEGER INT(*)
!
!***FIRST EXECUTABLE STATEMENT  ELMHES
  LA = IGH - 1
  KP1 = LOW + 1
  if (LA  <  KP1) go to 200
!
  DO 180 M = KP1, LA
     MM1 = M - 1
     X = 0.0E0
     I = M
!
     DO 100 J = M, IGH
        if (ABS(A(J,MM1))  <=  ABS(X)) go to 100
        X = A(J,MM1)
        I = J
  100    CONTINUE
!
     INT(M) = I
     if (I  ==  M) go to 130
!    .......... INTERCHANGE ROWS AND COLUMNS OF A ..........
     DO 110 J = MM1, N
        Y = A(I,J)
        A(I,J) = A(M,J)
        A(M,J) = Y
  110    CONTINUE
!
     DO 120 J = 1, IGH
        Y = A(J,I)
        A(J,I) = A(J,M)
        A(J,M) = Y
  120    CONTINUE
!    .......... END INTERCHANGE ..........
  130    if (X  ==  0.0E0) go to 180
     MP1 = M + 1
!
     DO 160 I = MP1, IGH
        Y = A(I,MM1)
        if (Y  ==  0.0E0) go to 160
        Y = Y / X
        A(I,MM1) = Y
!
        DO 140 J = M, N
  140       A(I,J) = A(I,J) - Y * A(M,J)
!
        DO 150 J = 1, IGH
  150       A(J,M) = A(J,M) + Y * A(J,I)
!
  160    CONTINUE
!
  180 CONTINUE
!
  200 RETURN
end
