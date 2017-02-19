subroutine BALBAK (NM, N, LOW, IGH, SCALE, M, Z)
!
!! BALBAK forms the eigenvectors of a real general matrix from the ...
!  eigenvectors of matrix output from BALANC.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (BALBAK-S, CBABK2-C)
!***KEYWORDS  EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure BALBAK,
!     NUM. MATH. 13, 293-304(1969) by Parlett and Reinsch.
!     HANDBOOK FOR AUTO. COMP., Vol.II-LINEAR ALGEBRA, 315-326(1971).
!
!     This subroutine forms the eigenvectors of a REAL GENERAL
!     matrix by back transforming those of the corresponding
!     balanced matrix determined by  BALANC.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameter, Z, as declared in the calling program
!          dimension statement.  NM is an INTEGER variable.
!
!        N is the number of components of the vectors in matrix Z.
!          N is an INTEGER variable.  N must be less than or equal
!          to NM.
!
!        LOW and IGH are INTEGER variables determined by  BALANC.
!
!        SCALE contains information determining the permutations and
!          scaling factors used by  BALANC.  SCALE is a one-dimensional
!          REAL array, dimensioned SCALE(N).
!
!        M is the number of columns of Z to be back transformed.
!          M is an INTEGER variable.
!
!        Z contains the real and imaginary parts of the eigen-
!          vectors to be back transformed in its first M columns.
!          Z is a two-dimensional REAL array, dimensioned Z(NM,M).
!
!     On OUTPUT
!
!        Z contains the real and imaginary parts of the
!          transformed eigenvectors in its first M columns.
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
!***END PROLOGUE  BALBAK
!
  INTEGER I,J,K,M,N,II,NM,IGH,LOW
  REAL SCALE(*),Z(NM,*)
  REAL S
!
!***FIRST EXECUTABLE STATEMENT  BALBAK
  if (M  ==  0) go to 200
  if (IGH  ==  LOW) go to 120
!
  DO 110 I = LOW, IGH
     S = SCALE(I)
!     .......... LEFT HAND EIGENVECTORS ARE BACK TRANSFORMED
!                if THE FOREGOING STATEMENT IS REPLACED BY
!                S=1.0E0/SCALE(I). ..........
     DO 100 J = 1, M
  100    Z(I,J) = Z(I,J) * S
!
  110 CONTINUE
!     ......... FOR I=LOW-1 STEP -1 UNTIL 1,
!               IGH+1 STEP 1 UNTIL N DO -- ..........
  120 DO 140 II = 1, N
     I = II
     if (I  >=  LOW .AND. I  <=  IGH) go to 140
     if (I  <  LOW) I = LOW - II
     K = SCALE(I)
     if (K  ==  I) go to 140
!
     DO 130 J = 1, M
        S = Z(I,J)
        Z(I,J) = Z(K,J)
        Z(K,J) = S
  130    CONTINUE
!
  140 CONTINUE
!
  200 RETURN
end
