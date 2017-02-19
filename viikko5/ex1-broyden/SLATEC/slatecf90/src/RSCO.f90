subroutine RSCO (RSAV, ISAV)
!
!! RSCO is subsidiary to DEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (RSCO-S, DRSCO-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   RSCO transfers data from arrays to a common block within the
!   integrator package DEBDF.
!
!***SEE ALSO  DEBDF
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DEBDF1
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  RSCO
!
!
!-----------------------------------------------------------------------
! THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON
! BLOCK DEBDF1  , WHICH IS USED INTERNALLY IN THE DEBDF
! PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS
! OF SUBROUTINE SVCO OR THE EQUIVALENT.
!-----------------------------------------------------------------------
  INTEGER ISAV, I,      ILS, LENILS, LENRLS
  REAL RSAV, RLS
  DIMENSION RSAV(*), ISAV(*)
  COMMON /DEBDF1/ RLS(218), ILS(33)
  SAVE LENRLS, LENILS
  DATA LENRLS/218/, LENILS/33/
!
!***FIRST EXECUTABLE STATEMENT  RSCO
  DO 10 I = 1,LENRLS
 10     RLS(I) = RSAV(I)
  DO 20 I = 1,LENILS
 20     ILS(I) = ISAV(I)
  return
!----------------------- END OF SUBROUTINE RSCO -----------------------
end
