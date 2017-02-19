subroutine ELTRAN (NM, N, LOW, IGH, A, INT, Z)
!
!! ELTRAN accumulates the stabilized elementary similarity ...
!            transformations used in the reduction of a real general
!            matrix to upper Hessenberg form by ELMHES.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (ELTRAN-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure ELMTRANS,
!     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!     This subroutine accumulates the stabilized elementary
!     similarity transformations used in the reduction of a
!     REAL GENERAL matrix to upper Hessenberg form by  ELMHES.
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
!        A contains the multipliers which were used in the reduction
!          by  ELMHES  in its lower triangle below the subdiagonal.
!          A is a two-dimensional REAL array, dimensioned A(NM,IGH).
!
!        INT contains information on the rows and columns interchanged
!          in the reduction by  ELMHES.  Only elements LOW through IGH
!          are used.  INT is a one-dimensional INTEGER array,
!          dimensioned INT(IGH).
!
!     On OUTPUT
!
!        Z contains the transformation matrix produced in the reduction
!          by  ELMHES.  Z is a two-dimensional REAL array, dimensioned
!          Z(NM,N).
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
!***END PROLOGUE  ELTRAN
!
  INTEGER I,J,N,KL,MM,MP,NM,IGH,LOW,MP1
  REAL A(NM,*),Z(NM,*)
  INTEGER INT(*)
!
!***FIRST EXECUTABLE STATEMENT  ELTRAN
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
     MP1 = MP + 1
!
     DO 100 I = MP1, IGH
  100    Z(I,MP) = A(I,MP-1)
!
     I = INT(MP)
     if (I  ==  MP) go to 140
!
     DO 130 J = MP, IGH
        Z(MP,J) = Z(I,J)
        Z(I,J) = 0.0E0
  130    CONTINUE
!
     Z(I,MP) = 1.0E0
  140 CONTINUE
!
  200 RETURN
end
