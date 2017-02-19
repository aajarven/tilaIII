subroutine DSVCO (RSAV, ISAV)
!
!! DSVCO transfers data from a common block to arrays for DDEBDF.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DDEBDF
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SVCO-S, DSVCO-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   DSVCO transfers data from a common block to arrays within the
!   integrator package DDEBDF.
!
!***SEE ALSO  DDEBDF
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DDEBD1
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DSVCO
!-----------------------------------------------------------------------
! THIS ROUTINE STORES IN RSAV AND ISAV THE CONTENTS OF COMMON BLOCK
! DDEBD1  , WHICH IS USED INTERNALLY IN THE DDEBDF PACKAGE.
!
! RSAV = DOUBLE PRECISION ARRAY OF LENGTH 218 OR MORE.
! ISAV = INTEGER ARRAY OF LENGTH 33 OR MORE.
!-----------------------------------------------------------------------
  INTEGER I, ILS, ISAV, LENILS, LENRLS
  DOUBLE PRECISION RLS, RSAV
  DIMENSION RSAV(*),ISAV(*)
  SAVE LENRLS, LENILS
  COMMON /DDEBD1/ RLS(218),ILS(33)
  DATA LENRLS /218/, LENILS /33/
!
!***FIRST EXECUTABLE STATEMENT  DSVCO
  DO 10 I = 1, LENRLS
     RSAV(I) = RLS(I)
   10 CONTINUE
  DO 20 I = 1, LENILS
     ISAV(I) = ILS(I)
   20 CONTINUE
  return
!     ----------------------- END OF SUBROUTINE DSVCO
!     -----------------------
end
