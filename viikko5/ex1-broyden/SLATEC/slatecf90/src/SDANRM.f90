FUNCTION SDANRM (NEQ, V, WT, RPAR, IPAR)
!
!! SDANRM computes vector norms for SDASSL.
!
!***LIBRARY   SLATEC (DASSL)
!***TYPE      SINGLE PRECISION (SDANRM-S, DDANRM-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     THIS FUNCTION ROUTINE COMPUTES THE WEIGHTED
!     ROOT-MEAN-SQUARE NORM OF THE VECTOR OF LENGTH
!     NEQ CONTAINED IN THE ARRAY V,WITH WEIGHTS
!     CONTAINED IN THE ARRAY WT OF LENGTH NEQ.
!        SDANRM=SQRT((1/NEQ)*SUM(V(I)/WT(I))**2)
!-----------------------------------------------------------------------
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!***END PROLOGUE  SDANRM
!
  real SDANRM
  INTEGER  NEQ, IPAR(*)
  REAL  V(NEQ), WT(NEQ), RPAR(*)
!
  INTEGER  I
  REAL  SUM, VMAX
!
!***FIRST EXECUTABLE STATEMENT  SDANRM
  SDANRM = 0.0E0
  VMAX = 0.0E0
  DO 10 I = 1,NEQ
    if ( ABS(V(I)/WT(I))  >  VMAX) VMAX = ABS(V(I)/WT(I))
10      CONTINUE
  if ( VMAX  <=  0.0E0) go to 30
  SUM = 0.0E0
  DO 20 I = 1,NEQ
20      SUM = SUM + ((V(I)/WT(I))/VMAX)**2
  SDANRM = VMAX*SQRT(SUM/NEQ)
30    CONTINUE
  return
!------END OF FUNCTION SDANRM------
end
