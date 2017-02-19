subroutine DEFEHL (F, NEQ, T, Y, H, YP, F1, F2, F3, F4, F5, YS, &
     RPAR, IPAR)
!
!! DEFEHL is subsidiary to DERKF.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (DEFEHL-S, DFEHL-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     Fehlberg Fourth-Fifth order Runge-Kutta Method
! **********************************************************************
!
!    DEFEHL integrates a system of NEQ first order
!    ordinary differential equations of the form
!               dU/DX = F(X,U)
!    over one step when the vector Y(*) of initial values for U(*) and
!    the vector YP(*) of initial derivatives, satisfying  YP = F(T,Y),
!    are given at the starting point X=T.
!
!    DEFEHL advances the solution over the fixed step H and returns
!    the fifth order (sixth order accurate locally) solution
!    approximation at T+H in the array YS(*).
!    F1,---,F5 are arrays of dimension NEQ which are needed
!    for internal storage.
!    The formulas have been grouped to control loss of significance.
!    DEFEHL should be called with an H not smaller than 13 units of
!    roundoff in T so that the various independent arguments can be
!    distinguished.
!
!    This subroutine has been written with all variables and statement
!    numbers entirely compatible with DERKFS. For greater efficiency,
!    the call to DEFEHL can be replaced by the module beginning with
!    line 222 and extending to the last line just before the return
!    statement.
!
! **********************************************************************
!
!***SEE ALSO  DERKF
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   800501  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891009  Removed unreferenced statement label.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DEFEHL
!
!
  DIMENSION Y(*),YP(*),F1(*),F2(*),F3(*),F4(*),F5(*), &
            YS(*),RPAR(*),IPAR(*)
!
!***FIRST EXECUTABLE STATEMENT  DEFEHL
  CH=H/4.
  DO 230 K=1,NEQ
  230   YS(K)=Y(K)+CH*YP(K)
  call F(T+CH,YS,F1,RPAR,IPAR)
!
  CH=3.*H/32.
  DO 240 K=1,NEQ
  240   YS(K)=Y(K)+CH*(YP(K)+3.*F1(K))
  call F(T+3.*H/8.,YS,F2,RPAR,IPAR)
!
  CH=H/2197.
  DO 250 K=1,NEQ
  250   YS(K)=Y(K)+CH*(1932.*YP(K)+(7296.*F2(K)-7200.*F1(K)))
  call F(T+12.*H/13.,YS,F3,RPAR,IPAR)
!
  CH=H/4104.
  DO 260 K=1,NEQ
  260   YS(K)=Y(K)+CH*((8341.*YP(K)-845.*F3(K))+ &
                              (29440.*F2(K)-32832.*F1(K)))
  call F(T+H,YS,F4,RPAR,IPAR)
!
  CH=H/20520.
  DO 270 K=1,NEQ
  270   YS(K)=Y(K)+CH*((-6080.*YP(K)+(9295.*F3(K)-5643.*F4(K)))+ &
                               (41040.*F1(K)-28352.*F2(K)))
  call F(T+H/2.,YS,F5,RPAR,IPAR)
!
!     COMPUTE APPROXIMATE SOLUTION AT T+H
!
  CH=H/7618050.
  DO 290 K=1,NEQ
  290   YS(K)=Y(K)+CH*((902880.*YP(K)+(3855735.*F3(K)-1371249.*F4(K)))+ &
                  (3953664.*F2(K)+277020.*F5(K)))
!
  return
end
