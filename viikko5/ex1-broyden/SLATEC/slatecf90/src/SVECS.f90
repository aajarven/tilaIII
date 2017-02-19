subroutine SVECS (NCOMP, LNFC, YHP, WORK, IWORK, INHOMO, IFLAG)
!
!! SVECS is subsidiary to BVSUP.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BVSUP
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SVECS-S, DVECS-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!  This subroutine is used for the special structure of complex valued
!  problems. MGSBV is called upon to obtain LNFC vectors from an
!  original set of 2*LNFC independent vectors so that the resulting
!  LNFC vectors together with their imaginary product or mate vectors
!  form an independent set.
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  MGSBV
!***COMMON BLOCKS    ML18JR
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  SVECS
!
  DIMENSION YHP(NCOMP,*),WORK(*),IWORK(*)
  COMMON /ML18JR/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,LNFCC, &
                  ICOCO
!***FIRST EXECUTABLE STATEMENT  SVECS
  if (LNFC  ==  1) go to 5
  NIV=LNFC
  LNFC=2*LNFC
  LNFCC=2*LNFCC
  KP=LNFC+2+LNFCC
  IDP=INDPVT
  INDPVT=0
  call MGSBV(NCOMP,LNFC,YHP,NCOMP,NIV,IFLAG,WORK(1),WORK(KP), &
           IWORK(1),INHOMO,YHP(1,LNFC+1),WORK(LNFC+2),DUM)
  LNFC=LNFC/2
  LNFCC=LNFCC/2
  INDPVT=IDP
  if (IFLAG  ==  0  .AND.  NIV  ==  LNFC) go to 5
  IFLAG=99
  return
    5 DO 6 K=1,NCOMP
    6 YHP(K,LNFC+1)=YHP(K,LNFCC+1)
  IFLAG=1
  return
end
