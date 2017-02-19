subroutine RKFAB (NCOMP, XPTS, NXPTS, NFC, IFLAG, Z, MXNON, P, &
     NTP, IP, YHP, NIV, U, V, W, S, STOWA, G, WORK, IWORK, NFCC)
!
!! RKFAB integrates an initial value problem for BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (RKFAB-S, DRKFAB-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!
!     Subroutine RKFAB integrates the initial value equations using
!     the variable-step RUNGE-KUTTA-FEHLBERG integration scheme or
!     the variable-order ADAMS method and orthonormalization
!     determined by a linear dependence test.
!
! **********************************************************************
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  BVDER, DEABM, DERKF, REORT, STOR1
!***COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  RKFAB
!
  DIMENSION P(NTP,*),IP(NFCC,*),U(NCOMP,NFC,*), &
            V(NCOMP,*),W(NFCC,*),Z(*),YHP(NCOMP,*), &
            XPTS(*),S(*),STOWA(*),WORK(*),IWORK(*), &
            G(*)
!
! **********************************************************************
!
  COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFCD
  COMMON /ML15TO/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
  COMMON /ML18JR/ AE,RE,TOL,NXPTSD,NIC,NOPG,MXNOND,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTPD,NEQIVP,NUMORT,NFCCD, &
                  ICOCO
  COMMON /ML17BW/ KKKZPW,NEEDW,NEEDIW,K1,K2,K3,K4,K5,K6,K7,K8,K9, &
                  K10,K11,L1,L2,KKKINT,LLLINT
!
  EXTERNAL BVDER
!
! **********************************************************************
!  INITIALIZATION OF COUNTERS AND VARIABLES.
!
!***FIRST EXECUTABLE STATEMENT  RKFAB
  KOD = 1
  NON = 1
  X = XBEG
  JON = 1
  INFO(1) = 0
  INFO(2) = 0
  INFO(3) = 1
  INFO(4) = 1
  WORK(1) = XEND
  if (NOPG  ==  0)  go to 1
  INFO(3) = 0
  if (X  ==  Z(1))  JON = 2
    1 NFCP1 = NFC + 1
!
! **********************************************************************
! *****BEGINNING OF INTEGRATION LOOP AT OUTPUT POINTS.******************
! **********************************************************************
!
  DO 110 KOPP = 2,NXPTS
  KOP=KOPP
!
    5 XOP = XPTS(KOP)
  if (NDISK  ==  0)  KOD = KOP
!
!     STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS.
!
   10 XXOP = XOP
  if (NOPG  ==  0)   go to 15
  if (XEND > XBEG.AND.XOP > Z(JON)) XXOP=Z(JON)
  if (XEND < XBEG.AND.XOP < Z(JON)) XXOP=Z(JON)
!
! **********************************************************************
   15 go to (20,25),INTEG
!     DERKF INTEGRATOR
!
   20 call DERKF(BVDER,NEQ,X,YHP,XXOP,INFO,RE,AE,IDID,WORK,KKKINT, &
             IWORK,LLLINT,G,IPAR)
  go to 28
!     DEABM INTEGRATOR
!
   25 call DEABM(BVDER,NEQ,X,YHP,XXOP,INFO,RE,AE,IDID,WORK,KKKINT, &
             IWORK,LLLINT,G,IPAR)
   28 if ( IDID  >=  1) go to 30
  INFO(1) = 1
  if ( IDID  ==  -1) go to 15
  IFLAG = 20 - IDID
  return
!
! **********************************************************************
!     GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR ORTHONORMALIZATION
!     (TEMPORARILY USING U AND V IN THE TEST)
!
   30 if (NOPG  ==  0)  go to 35
  if (XXOP  /=  Z(JON))  go to 100
  JFLAG=2
  go to 40
   35 JFLAG=1
  if (INHOMO  ==  3  .AND.  X  ==  XEND) JFLAG=3
!
   40 if (NDISK  ==  0) NON=NUMORT+1
  call REORT(NCOMP,U(1,1,KOD),V(1,KOD),YHP,NIV, &
             W(1,NON),S,P(1,NON),IP(1,NON),STOWA,JFLAG)
!
  if (JFLAG  /=  30) go to 45
  IFLAG=30
  return
!
   45 if (JFLAG  ==  10) go to 5
!
  if (JFLAG  /=  0)  go to 100
!
! **********************************************************************
!     STORE ORTHONORMALIZED VECTORS INTO SOLUTION VECTORS.
!
  if (NUMORT  <  MXNON)  go to 65
  if (X  ==  XEND) go to 65
  IFLAG = 13
  return
!
   65 NUMORT = NUMORT + 1
  call STOR1(YHP,U(1,1,KOD),YHP(1,NFCP1),V(1,KOD),1, &
             NDISK,NTAPE)
!
! **********************************************************************
!     STORE ORTHONORMALIZATION INFORMATION, INITIALIZE
!     INTEGRATION FLAG, AND CONTINUE INTEGRATION TO THE NEXT
!     ORTHONORMALIZATION POINT OR OUTPUT POINT.
!
  Z(NUMORT) = X
  if (INHOMO  ==  1  .AND.  NPS  ==  0)  C = S(NFCP1) * C
  if (NDISK  ==  0)  go to 90
  if (INHOMO  ==  1)  WRITE (NTAPE) (W(J,1), J = 1,NFCC)
  WRITE(NTAPE) (IP(J,1), J = 1,NFCC),(P(J,1), J = 1,NTP)
   90 INFO(1) = 0
  JON = JON + 1
  if (NOPG  ==  1  .AND.  X  /=  XOP)  go to 10
!
! **********************************************************************
!     CONTINUE INTEGRATION if WE ARE NOT AT AN OUTPUT POINT.
!
  100 if (IDID  ==  1)  go to 15
!
!     STORAGE OF HOMOGENEOUS SOLUTIONS IN U AND THE PARTICULAR
!     SOLUTION IN V AT THE OUTPUT POINTS.
!
  call STOR1(U(1,1,KOD),YHP,V(1,KOD),YHP(1,NFCP1),0,NDISK,NTAPE)
  110 CONTINUE
! **********************************************************************
! **********************************************************************
!
  IFLAG = 0
  return
end
