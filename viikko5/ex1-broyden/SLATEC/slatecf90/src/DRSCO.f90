subroutine DRSCO (RSAV, ISAV)
!
!! DRSCO transfers data from arrays to common blocks for DDEBDF.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DDEBDF
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (RSCO-S, DRSCO-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   DRSCO transfers data from arrays to a common block within the
!   integrator package DDEBDF.
!
!***SEE ALSO  DDEBDF
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DDEBD1
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DRSCO
!-----------------------------------------------------------------------
! THIS ROUTINE RESTORES FROM RSAV AND ISAV THE CONTENTS OF COMMON
! BLOCK DDEBD1  , WHICH IS USED INTERNALLY IN THE DDEBDF
! PACKAGE.  THIS PRESUMES THAT RSAV AND ISAV WERE LOADED BY MEANS
! OF SUBROUTINE DSVCO OR THE EQUIVALENT.
!-----------------------------------------------------------------------
!
  INTEGER I, ILS, ISAV, LENILS, LENRLS
  DOUBLE PRECISION RLS, RSAV
  DIMENSION RSAV(*),ISAV(*)
  SAVE LENRLS, LENILS
  COMMON /DDEBD1/ RLS(218),ILS(33)
  DATA LENRLS /218/, LENILS /33/
!
!***FIRST EXECUTABLE STATEMENT  DRSCO
  DO 10 I = 1, LENRLS
     RLS(I) = RSAV(I)
   10 CONTINUE
  DO 20 I = 1, LENILS
     ILS(I) = ISAV(I)
   20 CONTINUE
  return
!     ----------------------- END OF SUBROUTINE DRSCO
!     -----------------------
end
