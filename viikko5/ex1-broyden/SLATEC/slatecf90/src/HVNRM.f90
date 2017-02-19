function HVNRM (V, NCOMP)
!
!! HVNRM is subsidiary to DEABM, DEBDF and DERKF.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (HVNRM-S, DHVNRM-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     Compute the maximum norm of the vector V(*) of length NCOMP and
!     return the result as HVNRM.
!
!***SEE ALSO  DEABM, DEBDF, DERKF
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800501  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891024  Changed routine name from VNORM to HVNRM.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  HVNRM
  DIMENSION V(*)
!***FIRST EXECUTABLE STATEMENT  HVNRM
  HVNRM=0.
  DO 10 K=1,NCOMP
   10   HVNRM=MAX(HVNRM,ABS(V(K)))
  return
end
