  DOUBLE PRECISION FUNCTION DVNRMS (N, V, W)
!
!! DVNRMS computes a weighted root-mean-square vector norm for DDEBDF.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DDEBDF
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (VNWRMS-S, DVNRMS-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   DVNRMS computes a weighted root-mean-square vector norm for the
!   integrator package DDEBDF.
!
!***SEE ALSO  DDEBDF
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DVNRMS
  INTEGER I, N
  DOUBLE PRECISION SUM, V, W
  DIMENSION V(*),W(*)
!***FIRST EXECUTABLE STATEMENT  DVNRMS
  SUM = 0.0D0
  DO 10 I = 1, N
     SUM = SUM + (V(I)/W(I))**2
   10 CONTINUE
  DVNRMS = SQRT(SUM/N)
  return
!     ----------------------- END OF FUNCTION DVNRMS
!     ------------------------
end
