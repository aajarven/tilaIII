subroutine CPEVLR (N, M, A, X, C)
!
!! CPEVLR is subsidiary to CPZERO.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CPEVLR-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  CPZERO
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   810223  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CPEVLR
  REAL A(*),C(*)
!***FIRST EXECUTABLE STATEMENT  CPEVLR
  NP1=N+1
  DO J=1,NP1
        CI=0.0
        CIM1=A(J)
        MINI=MIN(M+1,N+2-J)
        DO I=1,MINI
           if ( J  /=  1) CI=C(I)
           if ( I  /=  1) CIM1=C(I-1)
           C(I)=CIM1+X*CI
        end do
  end do

  return
end
