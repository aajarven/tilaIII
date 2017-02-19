subroutine DFEHL (DF, NEQ, T, Y, H, YP, F1, F2, F3, F4, F5, YS, &
     RPAR, IPAR)
!
!! DFEHL implements a (4,5) order Runge-Kutta-Fehlberg ODE method.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DDERKF
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (DEFEHL-S, DFEHL-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!     Fehlberg Fourth-Fifth Order Runge-Kutta Method
! **********************************************************************
!
!    DFEHL integrates a system of NEQ first order
!    ordinary differential equations of the form
!               DU/DX = DF(X,U)
!    over one step when the vector Y(*) of initial values for U(*) and
!    the vector YP(*) of initial derivatives, satisfying  YP = DF(T,Y),
!    are given at the starting point X=T.
!
!    DFEHL advances the solution over the fixed step H and returns
!    the fifth order (sixth order accurate locally) solution
!    approximation at T+H in the array YS(*).
!    F1,---,F5 are arrays of dimension NEQ which are needed
!    for internal storage.
!    The formulas have been grouped to control loss of significance.
!    DFEHL should be called with an H not smaller than 13 units of
!    roundoff in T so that the various independent arguments can be
!    distinguished.
!
!    This subroutine has been written with all variables and statement
!    numbers entirely compatible with DRKFS. For greater efficiency,
!    the call to DFEHL can be replaced by the module beginning with
!    line 222 and extending to the last line just before the return
!    statement.
!
! **********************************************************************
!
!***SEE ALSO  DDERKF
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DFEHL
!
  INTEGER IPAR, K, NEQ
  DOUBLE PRECISION CH, F1, F2, F3, F4, F5, H, RPAR, T, Y, YP, YS
  DIMENSION Y(*),YP(*),F1(*),F2(*),F3(*),F4(*),F5(*), &
            YS(*),RPAR(*),IPAR(*)
!
!***FIRST EXECUTABLE STATEMENT  DFEHL
  CH = H/4.0D0
  DO 10 K = 1, NEQ
     YS(K) = Y(K) + CH*YP(K)
   10 CONTINUE
  call DF(T+CH,YS,F1,RPAR,IPAR)
!
  CH = 3.0D0*H/32.0D0
  DO 20 K = 1, NEQ
     YS(K) = Y(K) + CH*(YP(K) + 3.0D0*F1(K))
   20 CONTINUE
  call DF(T+3.0D0*H/8.0D0,YS,F2,RPAR,IPAR)
!
  CH = H/2197.0D0
  DO 30 K = 1, NEQ
     YS(K) = Y(K) &
             + CH &
               *(1932.0D0*YP(K) + (7296.0D0*F2(K) - 7200.0D0*F1(K)))
   30 CONTINUE
  call DF(T+12.0D0*H/13.0D0,YS,F3,RPAR,IPAR)
!
  CH = H/4104.0D0
  DO 40 K = 1, NEQ
     YS(K) = Y(K) &
             + CH &
               *((8341.0D0*YP(K) - 845.0D0*F3(K)) &
                 + (29440.0D0*F2(K) - 32832.0D0*F1(K)))
   40 CONTINUE
  call DF(T+H,YS,F4,RPAR,IPAR)
!
  CH = H/20520.0D0
  DO 50 K = 1, NEQ
     YS(K) = Y(K) &
             + CH &
               *((-6080.0D0*YP(K) &
                  + (9295.0D0*F3(K) - 5643.0D0*F4(K))) &
                 + (41040.0D0*F1(K) - 28352.0D0*F2(K)))
   50 CONTINUE
  call DF(T+H/2.0D0,YS,F5,RPAR,IPAR)
!
!     COMPUTE APPROXIMATE SOLUTION AT T+H
!
  CH = H/7618050.0D0
  DO 60 K = 1, NEQ
     YS(K) = Y(K) &
             + CH &
               *((902880.0D0*YP(K) &
                  + (3855735.0D0*F3(K) - 1371249.0D0*F4(K))) &
                 + (3953664.0D0*F2(K) + 277020.0D0*F5(K)))
   60 CONTINUE
!
  return
end
