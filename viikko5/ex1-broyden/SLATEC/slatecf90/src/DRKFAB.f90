subroutine DRKFAB (NCOMP, XPTS, NXPTS, NFC, IFLAG, Z, MXNON, P, &
     NTP, IP, YHP, NIV, U, V, W, S, STOWA, G, WORK, IWORK, NFCC)
!
!! DRKFAB integrates an initial value problem for DBVSUP.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DBVSUP
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (RKFAB-S, DRKFAB-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!
!     Subroutine DRKFAB integrates the initial value equations using
!     the variable-step Runge-Kutta-Fehlberg integration scheme or
!     the variable-order Adams method and orthonormalization
!     determined by a linear dependence test.
!
! **********************************************************************
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DBVDER, DDEABM, DDERKF, DREORT, DSTOR1
!***COMMON BLOCKS    DML15T, DML17B, DML18J, DML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DRKFAB
!
  INTEGER ICOCO, IDID, IFLAG, IGOFX, INDPVT, INFO, INHOMO, INTEG, &
       IPAR, ISTKOP, IVP, J, JFLAG, JON, &
       K1, K10, K11, K2, K3, K4, K5, K6, K7, K8, K9, KKKINT, &
       KKKZPW, KNSWOT, KOD, KOP, KOPP, L1, L2, LLLINT, LOTJP, &
       MNSWOT, MXNON, MXNOND, NCOMP, NCOMPD, NDISK, NEEDIW, NEEDW, &
       NEQ, NEQIVP, NFC, NFCC, NFCCD, NFCD, NFCP1, NIC, NIV, NON, &
       NOPG, NPS, NSWOT, NTAPE, NTP, NTPD, NUMORT, NXPTS, NXPTSD, &
       IP(NFCC,*), IWORK(*)
  DOUBLE PRECISION AE, C, G(*), P(NTP,*), PWCND, PX, RE, &
       S(*), STOWA(*), TND, TOL, U(NCOMP,NFC,*), &
       V(NCOMP,*), W(NFCC,*), WORK(*), X, XBEG, XEND, XOP, &
       XOT, XPTS(*), XSAV, XXOP, YHP(NCOMP,*), Z(*)
!
!     ******************************************************************
!
  COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFCD
  COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
  COMMON /DML18J/ AE,RE,TOL,NXPTSD,NIC,NOPG,MXNOND,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTPD,NEQIVP,NUMORT,NFCCD, &
                  ICOCO
  COMMON /DML17B/ KKKZPW,NEEDW,NEEDIW,K1,K2,K3,K4,K5,K6,K7,K8,K9, &
                  K10,K11,L1,L2,KKKINT,LLLINT
!
  EXTERNAL DBVDER
!
!      *****************************************************************
!       INITIALIZATION OF COUNTERS AND VARIABLES.
!
!     BEGIN BLOCK PERMITTING ...EXITS TO 220
!        BEGIN BLOCK PERMITTING ...EXITS TO 10
!***FIRST EXECUTABLE STATEMENT  DRKFAB
        KOD = 1
        NON = 1
        X = XBEG
        JON = 1
        INFO(1) = 0
        INFO(2) = 0
        INFO(3) = 1
        INFO(4) = 1
        WORK(1) = XEND
!        ...EXIT
        if (NOPG  ==  0) go to 10
        INFO(3) = 0
        if (X  ==  Z(1)) JON = 2
   10    CONTINUE
     NFCP1 = NFC + 1
!
!        ***************************************************************
!        *****BEGINNING OF INTEGRATION LOOP AT OUTPUT
!        POINTS.******************
!        ***************************************************************
!
     DO 210 KOPP = 2, NXPTS
        KOP = KOPP
        XOP = XPTS(KOP)
        if (NDISK  ==  0) KOD = KOP
!
   20       CONTINUE
!
!              STEP BY STEP INTEGRATION LOOP BETWEEN OUTPUT POINTS.
!
!              BEGIN BLOCK PERMITTING ...EXITS TO 190
!                 BEGIN BLOCK PERMITTING ...EXITS TO 30
                 XXOP = XOP
!                 ...EXIT
                 if (NOPG  ==  0) go to 30
                 if (XEND  >  XBEG .AND. XOP  >  Z(JON)) &
                    XXOP = Z(JON)
                 if (XEND  <  XBEG .AND. XOP  <  Z(JON)) &
                    XXOP = Z(JON)
   30             CONTINUE
!
!                 ******************************************************
   40             CONTINUE
!                    BEGIN BLOCK PERMITTING ...EXITS TO 170
                    go to (50,60), INTEG
!                       DDERKF INTEGRATOR
!
   50                   CONTINUE
                       call DDERKF(DBVDER,NEQ,X,YHP,XXOP,INFO,RE,AE, &
                                   IDID,WORK,KKKINT,IWORK,LLLINT,G, &
                                   IPAR)
                    go to 70
!                       DDEABM INTEGRATOR
!
   60                   CONTINUE
                       call DDEABM(DBVDER,NEQ,X,YHP,XXOP,INFO,RE,AE, &
                                   IDID,WORK,KKKINT,IWORK,LLLINT,G, &
                                   IPAR)
   70                   CONTINUE
                    if (IDID  >=  1) go to 80
                       INFO(1) = 1
!                    ......EXIT
                       if (IDID  ==  -1) go to 170
                       IFLAG = 20 - IDID
!     .....................EXIT
                       go to 220
   80                   CONTINUE
!
!                       ************************************************
!                           GRAM-SCHMIDT ORTHOGONALIZATION TEST FOR
!                           ORTHONORMALIZATION (TEMPORARILY USING U AND
!                           V IN THE TEST)
!
                    if (NOPG  ==  0) go to 100
                       if (XXOP  ==  Z(JON)) go to 90
!
!                             ******************************************
!                                 CONTINUE INTEGRATION if WE ARE NOT AT
!                                 AN OUTPUT POINT.
!
!           ..................EXIT
                          if (IDID  /=  1) go to 200
!                    .........EXIT
                          go to 170
   90                      CONTINUE
                       JFLAG = 2
                    go to 110
  100                   CONTINUE
                       JFLAG = 1
                       if (INHOMO  ==  3 .AND. X  ==  XEND) &
                          JFLAG = 3
  110                   CONTINUE
!
                    if (NDISK  ==  0) NON = NUMORT + 1
                    call DREORT(NCOMP,U(1,1,KOD),V(1,KOD),YHP,NIV, &
                                W(1,NON),S,P(1,NON),IP(1,NON),STOWA, &
                                JFLAG)
!
                    if (JFLAG  /=  30) go to 120
                       IFLAG = 30
!     .....................EXIT
                       go to 220
  120                   CONTINUE
!
                    if (JFLAG  /=  10) go to 130
                       XOP = XPTS(KOP)
                       if (NDISK  ==  0) KOD = KOP
!              ............EXIT
                       go to 190
  130                   CONTINUE
!
                    if (JFLAG  ==  0) go to 140
!
!                          *********************************************
!                              CONTINUE INTEGRATION if WE ARE NOT AT AN
!                              OUTPUT POINT.
!
!           ...............EXIT
                       if (IDID  /=  1) go to 200
!                    ......EXIT
                       go to 170
  140                   CONTINUE
!
!                       ************************************************
!                           STORE ORTHONORMALIZED VECTORS INTO SOLUTION
!                           VECTORS.
!
                    if (NUMORT  <  MXNON) go to 150
                    if (X  ==  XEND) go to 150
                       IFLAG = 13
!     .....................EXIT
                       go to 220
  150                   CONTINUE
!
                    NUMORT = NUMORT + 1
                    call DSTOR1(YHP,U(1,1,KOD),YHP(1,NFCP1), &
                                V(1,KOD),1,NDISK,NTAPE)
!
!                       ************************************************
!                           STORE ORTHONORMALIZATION INFORMATION,
!                           INITIALIZE INTEGRATION FLAG, AND CONTINUE
!                           INTEGRATION TO THE NEXT ORTHONORMALIZATION
!                           POINT OR OUTPUT POINT.
!
                    Z(NUMORT) = X
                    if (INHOMO  ==  1 .AND. NPS  ==  0) &
                       C = S(NFCP1)*C
                    if (NDISK  ==  0) go to 160
                       if (INHOMO  ==  1) &
                          WRITE (NTAPE) (W(J,1), J = 1, NFCC)
                       WRITE (NTAPE) &
                             (IP(J,1), J = 1, NFCC), &
                             (P(J,1), J = 1, NTP)
  160                   CONTINUE
                    INFO(1) = 0
                    JON = JON + 1
!                 ......EXIT
                    if (NOPG  ==  1 .AND. X  /=  XOP) go to 180
!
!                       ************************************************
!                           CONTINUE INTEGRATION if WE ARE NOT AT AN
!                           OUTPUT POINT.
!
!           ............EXIT
                    if (IDID  /=  1) go to 200
  170                CONTINUE
              go to 40
  180             CONTINUE
  190          CONTINUE
        go to 20
  200       CONTINUE
!
!           STORAGE OF HOMOGENEOUS SOLUTIONS IN U AND THE PARTICULAR
!           SOLUTION IN V AT THE OUTPUT POINTS.
!
        call DSTOR1(U(1,1,KOD),YHP,V(1,KOD),YHP(1,NFCP1),0,NDISK, &
                    NTAPE)
  210    CONTINUE
!        ***************************************************************
!        ***************************************************************
!
     IFLAG = 0
  220 CONTINUE
  return
end
