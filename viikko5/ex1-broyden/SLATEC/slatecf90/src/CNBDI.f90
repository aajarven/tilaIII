subroutine CNBDI (ABE, LDA, N, ML, MU, IPVT, DET)
!
!! CNBDI computes the determinant of a band matrix using the factors
!  computed by CNBCO or CNBFA.
!
!***LIBRARY   SLATEC
!***CATEGORY  D3C2
!***TYPE      COMPLEX (SNBDI-S, DNBDI-D, CNBDI-C)
!***KEYWORDS  BANDED, DETERMINANT, LINEAR EQUATIONS, NONSYMMETRIC
!***AUTHOR  Voorhees, E. A., (LANL)
!***DESCRIPTION
!
!     CNBDI computes the determinant of a band matrix
!     using the factors computed by CNBCO or CNBFA.
!     If the inverse is needed, use CNBSL  N  times.
!
!     On Entry
!
!        ABE     COMPLEX(LDA, NC)
!                the output from CNBCO or CNBFA.
!                NC must be  >=  2*ML+MU+1 .
!
!        LDA     INTEGER
!                the leading dimension of the array  ABE .
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
!                the pivot vector from CNBCO or CNBFA.
!
!     On Return
!
!        DET     COMPLEX(2)
!                determinant of original matrix.
!                Determinant = DET(1) * 10.0**DET(2)
!                with  1.0  <=  CABS1(DET(1))  <  10.0
!                or  DET(1) = 0.0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800730  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CNBDI
  INTEGER LDA,N,ML,MU,IPVT(*)
  COMPLEX ABE(LDA,*),DET(2)
!
  REAL TEN
  INTEGER I
  COMPLEX ZDUM
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
!
!***FIRST EXECUTABLE STATEMENT  CNBDI
  DET(1) = (1.0E0,0.0E0)
  DET(2) = (0.0E0,0.0E0)
  TEN = 10.0E0
  DO 50 I = 1, N
     if (IPVT(I)  /=  I) DET(1) = -DET(1)
     DET(1) = ABE(I,ML+1)*DET(1)
     if (CABS1(DET(1))  ==  0.0E0) go to 60
   10    if (CABS1(DET(1))  >=  1.0E0) go to 20
        DET(1) = CMPLX(TEN,0.0E0)*DET(1)
        DET(2) = DET(2) - (1.0E0,0.0E0)
     go to 10
   20    CONTINUE
   30    if (CABS1(DET(1))  <  TEN) go to 40
        DET(1) = DET(1)/CMPLX(TEN,0.0E0)
        DET(2) = DET(2) + (1.0E0,0.0E0)
     go to 30
   40    CONTINUE
   50 CONTINUE
   60 CONTINUE
  return
end
