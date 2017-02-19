subroutine ORTRAN (NM, N, LOW, IGH, A, ORT, Z)
!
!! ORTRAN accumulates orthogonal similarity transformations in the ...
!            reduction of real general matrix by ORTHES.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (ORTRAN-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure ORTRANS,
!     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!     This subroutine accumulates the orthogonal similarity
!     transformations used in the reduction of a REAL GENERAL
!     matrix to upper Hessenberg form by  ORTHES.
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
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  BALANC.  If  BALANC  has not been
!          used, set LOW=1 and IGH equal to the order of the matrix, N.
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
!     On OUTPUT
!
!        Z contains the transformation matrix produced in the reduction
!          by  ORTHES  to the upper Hessenberg form.  Z is a two-
!          dimensional REAL array, dimensioned Z(NM,N).
!
!        ORT has been used for temporary storage as is not restored.
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
!***END PROLOGUE  ORTRAN
!
  INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
  REAL A(NM,*),ORT(*),Z(NM,*)
  REAL G
!
!     .......... INITIALIZE Z TO IDENTITY MATRIX ..........
!***FIRST EXECUTABLE STATEMENT  ORTRAN
  DO 80 I = 1, N
!
     DO 60 J = 1, N
   60    Z(I,J) = 0.0E0
!
     Z(I,I) = 1.0E0
   80 CONTINUE
!
  KL = IGH - LOW - 1
  if (KL  <  1) go to 200
!     .......... FOR MP=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  DO 140 MM = 1, KL
     MP = IGH - MM
     if (A(MP,MP-1)  ==  0.0E0) go to 140
     MP1 = MP + 1
!
     DO 100 I = MP1, IGH
  100    ORT(I) = A(I,MP-1)
!
     DO 130 J = MP, IGH
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
