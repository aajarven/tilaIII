subroutine CORTH (NM, N, LOW, IGH, AR, AI, ORTR, ORTI)
!
!! CORTH reduces a complex general matrix to complex upper Hessenberg ...
!  form using unitary similarity transformations.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B2
!***TYPE      COMPLEX (ORTHES-S, CORTH-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of a complex analogue of
!     the ALGOL procedure ORTHES, NUM. MATH. 12, 349-368(1968)
!     by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!     Given a COMPLEX GENERAL matrix, this subroutine
!     reduces a submatrix situated in rows and columns
!     LOW through IGH to upper Hessenberg form by
!     unitary similarity transformations.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, AR and AI, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix A=(AR,AI).  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  CBAL.  If  CBAL  has not been used,
!          set LOW=1 and IGH equal to the order of the matrix, N.
!
!        AR and AI contain the real and imaginary parts, respectively,
!          of the complex input matrix.  AR and AI are two-dimensional
!          REAL arrays, dimensioned AR(NM,N) and AI(NM,N).
!
!     On OUTPUT
!
!        AR and AI contain the real and imaginary parts, respectively,
!          of the Hessenberg matrix.  Information about the unitary
!          transformations used in the reduction is stored in the
!          remaining triangles under the Hessenberg matrix.
!
!        ORTR and ORTI contain further information about the unitary
!          transformations.  Only elements LOW through IGH are used.
!          ORTR and ORTI are one-dimensional REAL arrays, dimensioned
!          ORTR(IGH) and ORTI(IGH).
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
!***END PROLOGUE  CORTH
!
  INTEGER I,J,M,N,II,JJ,LA,MP,NM,IGH,KP1,LOW
  REAL AR(NM,*),AI(NM,*),ORTR(*),ORTI(*)
  REAL F,G,H,FI,FR,SCALE
  REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  CORTH
  LA = IGH - 1
  KP1 = LOW + 1
  if (LA  <  KP1) go to 200
!
  DO 180 M = KP1, LA
     H = 0.0E0
     ORTR(M) = 0.0E0
     ORTI(M) = 0.0E0
     SCALE = 0.0E0
!     .......... SCALE COLUMN (ALGOL TOL THEN NOT NEEDED) ..........
     DO 90 I = M, IGH
   90    SCALE = SCALE + ABS(AR(I,M-1)) + ABS(AI(I,M-1))
!
     if (SCALE  ==  0.0E0) go to 180
     MP = M + IGH
!     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
     DO 100 II = M, IGH
        I = MP - II
        ORTR(I) = AR(I,M-1) / SCALE
        ORTI(I) = AI(I,M-1) / SCALE
        H = H + ORTR(I) * ORTR(I) + ORTI(I) * ORTI(I)
  100    CONTINUE
!
     G = SQRT(H)
     F = PYTHAG(ORTR(M),ORTI(M))
     if (F  ==  0.0E0) go to 103
     H = H + F * G
     G = G / F
     ORTR(M) = (1.0E0 + G) * ORTR(M)
     ORTI(M) = (1.0E0 + G) * ORTI(M)
     go to 105
!
  103    ORTR(M) = G
     AR(M,M-1) = SCALE
!     .......... FORM (I-(U*UT)/H) * A ..........
  105    DO 130 J = M, N
        FR = 0.0E0
        FI = 0.0E0
!     .......... FOR I=IGH STEP -1 UNTIL M DO -- ..........
        DO 110 II = M, IGH
           I = MP - II
           FR = FR + ORTR(I) * AR(I,J) + ORTI(I) * AI(I,J)
           FI = FI + ORTR(I) * AI(I,J) - ORTI(I) * AR(I,J)
  110       CONTINUE
!
        FR = FR / H
        FI = FI / H
!
        DO 120 I = M, IGH
           AR(I,J) = AR(I,J) - FR * ORTR(I) + FI * ORTI(I)
           AI(I,J) = AI(I,J) - FR * ORTI(I) - FI * ORTR(I)
  120       CONTINUE
!
  130    CONTINUE
!     .......... FORM (I-(U*UT)/H)*A*(I-(U*UT)/H) ..........
     DO 160 I = 1, IGH
        FR = 0.0E0
        FI = 0.0E0
!     .......... FOR J=IGH STEP -1 UNTIL M DO -- ..........
        DO 140 JJ = M, IGH
           J = MP - JJ
           FR = FR + ORTR(J) * AR(I,J) - ORTI(J) * AI(I,J)
           FI = FI + ORTR(J) * AI(I,J) + ORTI(J) * AR(I,J)
  140       CONTINUE
!
        FR = FR / H
        FI = FI / H
!
        DO 150 J = M, IGH
           AR(I,J) = AR(I,J) - FR * ORTR(J) - FI * ORTI(J)
           AI(I,J) = AI(I,J) + FR * ORTI(J) - FI * ORTR(J)
  150       CONTINUE
!
  160    CONTINUE
!
     ORTR(M) = SCALE * ORTR(M)
     ORTI(M) = SCALE * ORTI(M)
     AR(M,M-1) = -G * AR(M,M-1)
     AI(M,M-1) = -G * AI(M,M-1)
  180 CONTINUE
!
  200 RETURN
end
