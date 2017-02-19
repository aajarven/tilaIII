subroutine DDAWTS (NEQ, IWT, RTOL, ATOL, Y, WT, RPAR, IPAR)
!
!! DDAWTS sets the error weight vector for DDASSL.
!
!***LIBRARY   SLATEC (DASSL)
!***TYPE      DOUBLE PRECISION (SDAWTS-S, DDAWTS-D)
!***AUTHOR  Petzold, Linda R., (LLNL)
!***DESCRIPTION
!-----------------------------------------------------------------------
!     THIS SUBROUTINE SETS THE ERROR WEIGHT VECTOR
!     WT ACCORDING TO WT(I)=RTOL(I)*ABS(Y(I))+ATOL(I),
!     I=1,-,N.
!     RTOL AND ATOL ARE SCALARS if IWT = 0,
!     AND VECTORS if IWT = 1.
!-----------------------------------------------------------------------
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   830315  DATE WRITTEN
!   901009  Finished conversion to SLATEC 4.0 format (F.N.Fritsch)
!   901019  Merged changes made by C. Ulrich with SLATEC 4.0 format.
!   901026  Added explicit declarations for all variables and minor
!           cosmetic changes to prologue.  (FNF)
!***END PROLOGUE  DDAWTS
!
  INTEGER  NEQ, IWT, IPAR(*)
  DOUBLE PRECISION  RTOL(*), ATOL(*), Y(*), WT(*), RPAR(*)
!
  INTEGER  I
  DOUBLE PRECISION  ATOLI, RTOLI
!
!***FIRST EXECUTABLE STATEMENT  DDAWTS
  RTOLI=RTOL(1)
  ATOLI=ATOL(1)
  DO 20 I=1,NEQ
     if (IWT  == 0) go to 10
       RTOLI=RTOL(I)
       ATOLI=ATOL(I)
10         WT(I)=RTOLI*ABS(Y(I))+ATOLI
20         CONTINUE
  return
!-----------END OF SUBROUTINE DDAWTS------------------------------------
end
