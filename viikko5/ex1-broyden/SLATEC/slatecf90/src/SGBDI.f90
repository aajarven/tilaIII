subroutine SGBDI (ABD, LDA, N, ML, MU, IPVT, DET)
!
!! SGBDI computes the determinant of a band matrix using the factors ...
!            computed by SGBCO or SGBFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D3A2
!***TYPE      SINGLE PRECISION (SGBDI-S, DGBDI-D, CGBDI-C)
!***KEYWORDS  BANDED, DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK,
!             MATRIX
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     SGBDI computes the determinant of a band matrix
!     using the factors computed by SBGCO or SGBFA.
!     If the inverse is needed, use SGBSL  N  times.
!
!     On Entry
!
!        ABD     REAL(LDA, N)
!                the output from SBGCO or SGBFA.
!
!        LDA     INTEGER
!                the leading dimension of the array  ABD .
!
!        N       INTEGER
!                the order of the original matrix.
!
!        ML      INTEGER
!                number of diagonals below the main diagonal.
!
!        MU      INTEGER
!                number of diagonals above the main diagonal.
!
!        IPVT    INTEGER(N)
!                the pivot vector from SBGCO or SGBFA.
!
!     On Return
!
!        DET     REAL(2)
!                determinant of original matrix.
!                Determinant = DET(1) * 10.0**DET(2)
!                with  1.0  <=  ABS(DET(1))  <  10.0
!                or  DET(1) = 0.0 .
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
!***END PROLOGUE  SGBDI
  INTEGER LDA,N,ML,MU,IPVT(*)
  REAL ABD(LDA,*),DET(2)
!
  REAL TEN
  INTEGER I,M
!***FIRST EXECUTABLE STATEMENT  SGBDI
  M = ML + MU + 1
  DET(1) = 1.0E0
  DET(2) = 0.0E0
  TEN = 10.0E0
  DO 50 I = 1, N
     if (IPVT(I)  /=  I) DET(1) = -DET(1)
     DET(1) = ABD(M,I)*DET(1)
     if (DET(1)  ==  0.0E0) go to 60
   10    if (ABS(DET(1))  >=  1.0E0) go to 20
        DET(1) = TEN*DET(1)
        DET(2) = DET(2) - 1.0E0
     go to 10
   20    CONTINUE
   30    if (ABS(DET(1))  <  TEN) go to 40
        DET(1) = DET(1)/TEN
        DET(2) = DET(2) + 1.0E0
     go to 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
  return
end
