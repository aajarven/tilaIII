subroutine DSTOR1 (U, YH, V, YP, NTEMP, NDISK, NTAPE)
!
!! DSTOR1 is subsidiary to DBVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (STOR1-S, DSTOR1-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!             0 -- storage at output points.
!     NTEMP =
!             1 -- temporary storage
! **********************************************************************
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DSTOR1
  INTEGER IGOFX, INHOMO, IVP, J, NCOMP, NCTNF, NDISK, NFC, NTAPE, &
       NTEMP
  DOUBLE PRECISION C, U(*), V(*), XSAV, YH(*), YP(*)
!
!     ******************************************************************
!
  COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
!
!      *****************************************************************
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 80
!***FIRST EXECUTABLE STATEMENT  DSTOR1
     NCTNF = NCOMP*NFC
     DO 10 J = 1, NCTNF
        U(J) = YH(J)
   10    CONTINUE
     if (INHOMO  ==  1) go to 30
!
!           ZERO PARTICULAR SOLUTION
!
!     ......EXIT
        if (NTEMP  ==  1) go to 80
        DO 20 J = 1, NCOMP
           V(J) = 0.0D0
   20       CONTINUE
     go to 70
   30    CONTINUE
!
!           NONZERO PARTICULAR SOLUTION
!
        if (NTEMP  ==  0) go to 50
!
           DO 40 J = 1, NCOMP
              V(J) = YP(J)
   40          CONTINUE
!     .........EXIT
           go to 80
   50       CONTINUE
!
        DO 60 J = 1, NCOMP
           V(J) = C*YP(J)
   60       CONTINUE
   70    CONTINUE
!
!        IS OUTPUT INFORMATION TO BE WRITTEN TO DISK
!
     if (NDISK  ==  1) &
        WRITE (NTAPE) (V(J), J = 1, NCOMP),(U(J), J = 1, NCTNF)
   80 CONTINUE
!
  return
end
