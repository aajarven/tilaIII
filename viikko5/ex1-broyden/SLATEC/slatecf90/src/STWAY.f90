subroutine STWAY (U, V, YHP, INOUT, STOWA)
!
!! STWAY stores or recalls integration data for a restart of BVSUP.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BVSUP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (STWAY-S, DSTWAY-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!  This subroutine stores (recalls) integration data in the event
!  that a restart is needed (the homogeneous solution vectors become
!  too dependent to continue)
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  STOR1
!***COMMON BLOCKS    ML15TO, ML18JR, ML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  STWAY
!
  DIMENSION U(*),V(*),YHP(*),STOWA(*)
!
  COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
  COMMON /ML15TO/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
  COMMON /ML18JR/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC, &
                  ICOCO
!
!***FIRST EXECUTABLE STATEMENT  STWAY
  if (INOUT  ==  1) go to 100
!
!     SAVE IN STOWA ARRAY AND ISTKOP
!
  KS=NFC*NCOMP
  call STOR1(STOWA,U,STOWA(KS+1),V,1,0,0)
  KS=KS+NCOMP
  if (NEQIVP  ==  0) go to 50
  DO 25 J=1,NEQIVP
  KSJ=KS+J
   25 STOWA(KSJ)=YHP(KSJ)
   50 KS=KS+NEQIVP
  STOWA(KS+1)=X
  ISTKOP=KOP
  if (XOP  ==  X) ISTKOP=KOP+1
  return
!
!     RECALL FROM STOWA ARRAY AND ISTKOP
!
  100 KS=NFC*NCOMP
  call STOR1(YHP,STOWA,YHP(KS+1),STOWA(KS+1),1,0,0)
  KS=KS+NCOMP
  if (NEQIVP  ==  0) go to 150
  DO 125 J=1,NEQIVP
  KSJ=KS+J
  125 YHP(KSJ)=STOWA(KSJ)
  150 KS=KS+NEQIVP
  X=STOWA(KS+1)
  INFO(1)=0
  KO=KOP-ISTKOP
  KOP=ISTKOP
  if (NDISK  ==  0  .OR.  KO  ==  0) RETURN
  DO 175 K=1,KO
  175 BACKSPACE NTAPE
  return
end
