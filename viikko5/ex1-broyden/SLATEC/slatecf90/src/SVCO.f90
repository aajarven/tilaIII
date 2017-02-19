subroutine SVCO (RSAV, ISAV)
!
!! SVCO transfers data from a common block to arrays for DEBDF.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DEBDF
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SVCO-S, DSVCO-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   SVCO transfers data from a common block to arrays within the
!   integrator package DEBDF.
!
!***SEE ALSO  DEBDF
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DEBDF1
!***REVISION HISTORY  (YYMMDD)
!   800901  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  SVCO
!
!
!-----------------------------------------------------------------------
! THIS ROUTINE STORES IN RSAV AND ISAV THE CONTENTS OF COMMON BLOCK
! DEBDF1  , WHICH IS USED INTERNALLY IN THE DEBDF PACKAGE.
!
! RSAV = REAL ARRAY OF LENGTH 218 OR MORE.
! ISAV = INTEGER ARRAY OF LENGTH 33 OR MORE.
!-----------------------------------------------------------------------
  INTEGER ISAV, I,      ILS, LENILS, LENRLS
  REAL RSAV, RLS
  DIMENSION RSAV(*), ISAV(*)
  COMMON /DEBDF1/ RLS(218), ILS(33)
  SAVE LENRLS, LENILS
  DATA LENRLS/218/, LENILS/33/
!
!***FIRST EXECUTABLE STATEMENT  SVCO
  DO 10 I = 1,LENRLS
 10     RSAV(I) = RLS(I)
  DO 20 I = 1,LENILS
 20     ISAV(I) = ILS(I)
  return
!----------------------- END OF SUBROUTINE SVCO -----------------------
end
