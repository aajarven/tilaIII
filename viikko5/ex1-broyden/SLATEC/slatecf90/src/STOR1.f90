subroutine STOR1 (U, YH, V, YP, NTEMP, NDISK, NTAPE)
!
!! STOR1 is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (STOR1-S, DSTOR1-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!             0 -- Storage at output points.
!     NTEMP =
!             1 -- Temporary storage
! **********************************************************************
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    ML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  STOR1
  DIMENSION U(*),YH(*),V(*),YP(*)
!
! **********************************************************************
!
  COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
!
! **********************************************************************
!
!***FIRST EXECUTABLE STATEMENT  STOR1
  NCTNF = NCOMP * NFC
  DO 10 J = 1,NCTNF
   10 U(J) = YH(J)
  if (INHOMO  ==  1)  go to 30
!
!   ZERO PARTICULAR SOLUTION
!
  if (NTEMP  ==  1)  return
  DO 20 J = 1,NCOMP
   20 V(J) = 0.
  go to 70
!
!   NONZERO PARTICULAR SOLUTION
!
   30 if (NTEMP  ==  0)  go to 50
!
  DO 40 J = 1,NCOMP
   40 V(J) = YP(J)
  return
!
   50 DO 60 J = 1,NCOMP
   60 V(J) = C * YP(J)
!
!  IS OUTPUT INFORMATION TO BE WRITTEN TO DISK
!
   70 if (NDISK  ==  1)  WRITE (NTAPE) (V(J),J=1,NCOMP),(U(J),J=1,NCTNF)
!
  return
end
