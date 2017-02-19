FUNCTION VNWRMS (N, V, W)
!
!! VNWRMS is subsidiary to DEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (VNWRMS-S, DVNRMS-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   VNWRMS computes a weighted root-mean-square vector norm for the
!   integrator package DEBDF.
!
!***SEE ALSO  DEBDF
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  VNWRMS
!
!
!LLL. OPTIMIZE
!-----------------------------------------------------------------------
! THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED ROOT-MEAN-SQUARE NORM
! OF THE VECTOR OF LENGTH N CONTAINED IN THE ARRAY V, WITH WEIGHTS
! CONTAINED IN THE ARRAY W OF LENGTH N..
!   VNWRMS = SQRT( (1/N) * SUM( V(I)/W(I) )**2 )
!-----------------------------------------------------------------------
  INTEGER N, I
  real VNWRMS
  REAL V, W, SUM
  DIMENSION V(*), W(*)
!***FIRST EXECUTABLE STATEMENT  VNWRMS
  SUM = 0.0E0
  DO 10 I = 1,N
 10     SUM = SUM + (V(I)/W(I))**2
  VNWRMS = SQRT(SUM/N)
  return
!----------------------- END OF FUNCTION VNWRMS ------------------------
end
