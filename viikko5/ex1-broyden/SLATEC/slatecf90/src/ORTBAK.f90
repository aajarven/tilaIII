subroutine ORTBAK (NM, LOW, IGH, A, ORT, M, Z)
!
!! ORTBAK forms the eigenvectors of a general real matrix from the ...
!  eigenvectors of the upper Hessenberg matrix output from ORTHES.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (ORTBAK-S, CORTB-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure ORTBAK,
!     NUM. MATH. 12, 349-368(1968) by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 339-358(1971).
!
!     This subroutine forms the eigenvectors of a REAL GENERAL
!     matrix by back transforming those of the corresponding
!     upper Hessenberg matrix determined by  ORTHES.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix.
!
!        A contains some information about the orthogonal trans-
!          formations used in the reduction to Hessenberg form by
!          ORTHES  in its strict lower triangle.  A is a two-dimensional
!          REAL array, dimensioned A(NM,IGH).
!
!        ORT contains further information about the orthogonal trans-
!          formations used in the reduction by  ORTHES.  Only elements
!          LOW through IGH are used.  ORT is a one-dimensional REAL
!          array, dimensioned ORT(IGH).
!
!        M is the number of columns of Z to be back transformed.
!          M is an INTEGER variable.
!
!        Z contains the real and imaginary parts of the eigenvectors to
!          be back transformed in its first M columns.  Z is a two-
!          dimensional REAL array, dimensioned Z(NM,M).
!
!     On OUTPUT
!
!        Z contains the real and imaginary parts of the transformed
!          eigenvectors in its first M columns.
!
!        ORT has been used for temporary storage as is not restored.
!
!     NOTE that ORTBAK preserves vector Euclidean norms.
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
!***END PROLOGUE  ORTBAK
!
  INTEGER I,J,M,LA,MM,MP,NM,IGH,KP1,LOW,MP1
  REAL A(NM,*),ORT(*),Z(NM,*)
  REAL G
!
!***FIRST EXECUTABLE STATEMENT  ORTBAK
  if (M  ==  0) go to 200
  LA = IGH - 1
  KP1 = LOW + 1
  if (LA  <  KP1) go to 200
!     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  DO 140 MM = KP1, LA
     MP = LOW + IGH - MM
     if (A(MP,MP-1)  ==  0.0E0) go to 140
     MP1 = MP + 1
!
     DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
!
     DO 130 J = 1, M
        G = 0.0E0
!
        DO 110 I = MP, IGH
  110       G = G + ORT(I) * Z(I,J)
!     .......... DIVISOR BELOW IS NEGATIVE OF H FORMED IN ORTHES.
!                DOUBLE DIVISION AVOIDS POSSIBLE UNDERFLOW ..........
        G = (G / ORT(MP)) / A(MP,MP-1)
!
        DO 120 I = MP, IGH
  120       Z(I,J) = Z(I,J) + G * ORT(I)
!
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
end
