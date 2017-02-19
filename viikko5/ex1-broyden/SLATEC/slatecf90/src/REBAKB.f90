subroutine REBAKB (NM, N, B, DL, M, Z)
!
!! REBAKB forms the eigenvectors of a generalized symmetric ...
!            eigensystem from the eigenvectors of derived matrix output
!            from REDUC2.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C4
!***TYPE      SINGLE PRECISION (REBAKB-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure REBAKB,
!     NUM. MATH. 11, 99-110(1968) by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 303-314(1971).
!
!     This subroutine forms the eigenvectors of a generalized
!     SYMMETRIC eigensystem by back transforming those of the
!     derived symmetric matrix determined by  REDUC2.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, B and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix system.  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        B contains information about the similarity transformation
!          (Cholesky decomposition) used in the reduction by  REDUC2
!          in its strict lower triangle.  B is a two-dimensional
!          REAL array, dimensioned B(NM,N).
!
!        DL contains further information about the transformation.
!          DL is a one-dimensional REAL array, dimensioned DL(N).
!
!        M is the number of eigenvectors to be back transformed.
!          M is an INTEGER variable.
!
!        Z contains the eigenvectors to be back transformed in its
!          first M columns.  Z is a two-dimensional REAL array
!          dimensioned Z(NM,M).
!
!     On Output
!
!        Z contains the transformed eigenvectors in its first
!          M columns.
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
!***END PROLOGUE  REBAKB
!
  INTEGER I,J,K,M,N,I1,II,NM
  REAL B(NM,*),DL(*),Z(NM,*)
  REAL X
!
!***FIRST EXECUTABLE STATEMENT  REBAKB
  if (M  ==  0) go to 200
!
  DO 100 J = 1, M
!     .......... FOR I=N STEP -1 UNTIL 1 DO -- ..........
     DO 100 II = 1, N
        I1 = N - II
        I = I1 + 1
        X = DL(I) * Z(I,J)
        if (I  ==  1) go to 80
!
        DO 60 K = 1, I1
   60       X = X + B(I,K) * Z(K,J)
!
   80       Z(I,J) = X
  100 CONTINUE
!
  200 RETURN
end
