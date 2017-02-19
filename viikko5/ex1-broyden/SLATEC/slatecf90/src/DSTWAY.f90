subroutine DSTWAY (U, V, YHP, INOUT, STOWA)
!
!! DSTWAY is subsidiary to DBVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (STWAY-S, DSTWAY-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!  This subroutine stores (recalls) integration data in the event
!  that a restart is needed (the homogeneous solution vectors become
!  too dependent to continue).
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DSTOR1
!***COMMON BLOCKS    DML15T, DML18J, DML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DSTWAY
!
  INTEGER ICOCO, IGOFX, INDPVT, INFO, INHOMO, INOUT, INTEG, ISTKOP, &
       IVP, J, K, KNSWOT, KO, KOP, KS, KSJ, LOTJP, MNSWOT, MXNON, &
       NCOMP, NDISK, NEQ, NEQIVP, NFC, NFCC, NIC, NOPG, NPS, NSWOT, &
       NTAPE, NTP, NUMORT, NXPTS
  DOUBLE PRECISION AE, C, PWCND, PX, RE, STOWA(*), TND, TOL, U(*), &
       V(*), X, XBEG, XEND, XOP, XOT, XSAV, YHP(*)
!
  COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
  COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
  COMMON /DML18J/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC, &
                  ICOCO
!
!***FIRST EXECUTABLE STATEMENT  DSTWAY
  if (INOUT  ==  1) go to 30
!
!        SAVE IN STOWA ARRAY AND ISTKOP
!
     KS = NFC*NCOMP
     call DSTOR1(STOWA,U,STOWA(KS+1),V,1,0,0)
     KS = KS + NCOMP
     if (NEQIVP  <  1) go to 20
     DO 10 J = 1, NEQIVP
        KSJ = KS + J
        STOWA(KSJ) = YHP(KSJ)
   10    CONTINUE
   20    CONTINUE
     KS = KS + NEQIVP
     STOWA(KS+1) = X
     ISTKOP = KOP
     if (XOP  ==  X) ISTKOP = KOP + 1
  go to 80
   30 CONTINUE
!
!        RECALL FROM STOWA ARRAY AND ISTKOP
!
     KS = NFC*NCOMP
     call DSTOR1(YHP,STOWA,YHP(KS+1),STOWA(KS+1),1,0,0)
     KS = KS + NCOMP
     if (NEQIVP  <  1) go to 50
     DO 40 J = 1, NEQIVP
        KSJ = KS + J
        YHP(KSJ) = STOWA(KSJ)
   40    CONTINUE
   50    CONTINUE
     KS = KS + NEQIVP
     X = STOWA(KS+1)
     INFO(1) = 0
     KO = KOP - ISTKOP
     KOP = ISTKOP
     if (NDISK  ==  0 .OR. KO  ==  0) go to 70
        DO 60 K = 1, KO
           BACKSPACE NTAPE
   60       CONTINUE
   70    CONTINUE
   80 CONTINUE
  return
end
