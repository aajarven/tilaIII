subroutine CROTG (CA, CB, C, S)
!
!! CROTG constructs a Givens transformation.
!
!***LIBRARY   SLATEC (BLAS)
!***CATEGORY  D1B10
!***TYPE      COMPLEX (SROTG-S, DROTG-D, CROTG-C)
!***KEYWORDS  BLAS, GIVENS ROTATION, GIVENS TRANSFORMATION,
!             LINEAR ALGEBRA, VECTOR
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!    Complex Givens transformation
!
!    Construct the Givens transformation
!
!             (C    S)
!       G  =  (      ),  C**2 + ABS(S)**2 =1,
!             (-S   C)
!
!    which zeros the second entry of the complex 2-vector (CA,CB)**T
!
!    The quantity CA/ABS(CA)*NORM(CA,CB) overwrites CA in storage.
!
!    Input:
!        CA (Complex)
!        CB (Complex)
!
!    Output:
!        CA (Complex)      CA/ABS(CA)*NORM(CA,CB)
!        C  (Real)
!        S  (Complex)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  CROTG
  COMPLEX CA, CB, S
  REAL C
  REAL NORM, SCALE
  COMPLEX ALPHA
!***FIRST EXECUTABLE STATEMENT  CROTG
  if (ABS(CA)  ==  0.0) THEN
    C = 0.0
    S = (1.0,0.0)
    CA = CB
  ELSE
    SCALE = ABS(CA) + ABS(CB)
    NORM = SCALE * SQRT((ABS(CA/SCALE))**2 + (ABS(CB/SCALE))**2)
    ALPHA = CA /ABS(CA)
    C = ABS(CA) / NORM
    S = ALPHA * CONJG(CB) / NORM
    CA = ALPHA * NORM
  end if
  return
end
