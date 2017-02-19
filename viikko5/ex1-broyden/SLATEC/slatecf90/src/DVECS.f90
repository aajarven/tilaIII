subroutine DVECS (NCOMP, LNFC, YHP, WORK, IWORK, INHOMO, IFLAG)
!
!! DVECS is subsidiary to DBVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SVECS-S, DVECS-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!  This subroutine is used for the special structure of COMPLEX*16
!  valued problems. DMGSBV is called upon to obtain LNFC vectors from an
!  original set of 2*LNFC independent vectors so that the resulting
!  LNFC vectors together with their imaginary product or mate vectors
!  form an independent set.
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DMGSBV
!***COMMON BLOCKS    DML18J
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891009  Removed unreferenced statement label.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DVECS
!
  INTEGER ICOCO, IDP, IFLAG, INDPVT, INHOMO, INTEG, IWORK(*), K, &
       KP, LNFC, LNFCC, MXNON, NCOMP, NDISK, NEQ, NEQIVP, NIC, NIV, &
       NOPG, NPS, NTAPE, NTP, NUMORT, NXPTS
  DOUBLE PRECISION AE, DUM, RE, TOL, WORK(*), YHP(NCOMP,*)
  COMMON /DML18J/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,LNFCC, &
                  ICOCO
!***FIRST EXECUTABLE STATEMENT  DVECS
     if (LNFC  /=  1) go to 20
        DO 10 K = 1, NCOMP
           YHP(K,LNFC+1) = YHP(K,LNFCC+1)
   10       CONTINUE
        IFLAG = 1
     go to 60
   20    CONTINUE
        NIV = LNFC
        LNFC = 2*LNFC
        LNFCC = 2*LNFCC
        KP = LNFC + 2 + LNFCC
        IDP = INDPVT
        INDPVT = 0
        call DMGSBV(NCOMP,LNFC,YHP,NCOMP,NIV,IFLAG,WORK(1),WORK(KP), &
                    IWORK(1),INHOMO,YHP(1,LNFC+1),WORK(LNFC+2),DUM)
        LNFC = LNFC/2
        LNFCC = LNFCC/2
        INDPVT = IDP
        if (IFLAG  /=  0 .OR. NIV  /=  LNFC) go to 40
           DO 30 K = 1, NCOMP
              YHP(K,LNFC+1) = YHP(K,LNFCC+1)
   30          CONTINUE
           IFLAG = 1
        go to 50
   40       CONTINUE
           IFLAG = 99
   50       CONTINUE
   60    CONTINUE
  CONTINUE
  return
end
