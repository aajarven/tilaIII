subroutine DPBDI (ABD, LDA, N, M, DET)
!
!! DPBDI computes the determinant of a symmetric positive definite
!  band matrix using the factors computed by DPBCO or DPBFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D3B2
!***TYPE      DOUBLE PRECISION (SPBDI-S, DPBDI-D, CPBDI-C)
!***KEYWORDS  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
!             MATRIX, POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     DPBDI computes the determinant
!     of a double precision symmetric positive definite band matrix
!     using the factors computed by DPBCO or DPBFA.
!     If the inverse is needed, use DPBSL  N  times.
!
!     On Entry
!
!        ABD     DOUBLE PRECISION(LDA, N)
!                the output from DPBCO or DPBFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        M       INTEGER
!                the number of diagonals above the main diagonal.
!
!     On Return
!
!        DET     DOUBLE PRECISION(2)
!                determinant of original matrix in the form
!                DETERMINANT = DET(1) * 10.0**DET(2)
!                with  1.0  <=  DET(1)  <  10.0
!                or  DET(1)  ==  0.0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DPBDI
  INTEGER LDA,N,M
  DOUBLE PRECISION ABD(LDA,*)
  DOUBLE PRECISION DET(2)
!
  DOUBLE PRECISION S
  INTEGER I
!***FIRST EXECUTABLE STATEMENT  DPBDI
!
!     COMPUTE DETERMINANT
!
  DET(1) = 1.0D0
  DET(2) = 0.0D0
  S = 10.0D0
  DO 50 I = 1, N
     DET(1) = ABD(M+1,I)**2*DET(1)
     if (DET(1)  ==  0.0D0) go to 60
   10    if (DET(1)  >=  1.0D0) go to 20
        DET(1) = S*DET(1)
        DET(2) = DET(2) - 1.0D0
     go to 10
   20    CONTINUE
   30    if (DET(1)  <  S) go to 40
        DET(1) = DET(1)/S
        DET(2) = DET(2) + 1.0D0
     go to 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
  return
end
