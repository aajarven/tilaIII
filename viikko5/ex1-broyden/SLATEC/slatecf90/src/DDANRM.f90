  DOUBLE PRECISION FUNCTION DDANRM (NEQ, V, WT, RPAR, IPAR)
!
!! DDANRM computes vector norms for DDASSL.
!
!***LIBRARY   SLATEC (DASSL)
!***TYPE      DOUBLE PRECISION (SDANRM-S, DDANRM-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
!     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
!     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
!     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
!        DDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
!-----------------------------------------------------------------------
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!***END PROLOGUE  DDANRM
!
  INTEGER  NEQ, IPAR(*)
  DOUBLE PRECISION  V(NEQ), WT(NEQ), RPAR(*)
!
  INTEGER  I
  DOUBLE PRECISION  SUM, VMAX
!
!***FIRST EXECUTABLE STATEMENT  DDANRM
  DDANRM = 0.0D0
  VMAX = 0.0D0

  DO I = 1, NEQ
    if ( ABS(V(I)/WT(I))  >  VMAX) VMAX = ABS(V(I)/WT(I))
  end do

  if ( 0.0D0 < VMAX ) then

    SUM = 0.0D0
    DO I = 1,NEQ
      SUM = SUM + ((V(I)/WT(I))/VMAX)**2
    end do

    DDANRM = VMAX*SQRT(SUM/NEQ)

  end if

  return
end
