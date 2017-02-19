subroutine DREORT (NCOMP, Y, YP, YHP, NIV, W, S, P, IP, STOWA, &
     IFLAG)
!
!! DREORT orthonormalizes the solution vector of a homogeneous system.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DBVSUP
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (REORT-S, DREORT-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!   INPUT
! *********
!     Y, YP and YHP = homogeneous solution matrix and particular
!                     solution vector to be orthonormalized.
!     IFLAG = 1 --  store YHP into Y and YP, test for
!                   reorthonormalization, orthonormalize if needed,
!                   save restart data.
!             2 --  store YHP into Y and YP, reorthonormalization,
!                   no restarts.
!                   (preset orthonormalization mode)
!             3 --  store YHP into Y and YP, reorthonormalization
!                   (when INHOMO=3 and X=XEND).
! **********************************************************************
!   OUTPUT
! *********
!     Y, YP = orthonormalized solutions.
!     NIV = number of independent vectors returned from DMGSBV.
!     IFLAG = 0 --  reorthonormalization was performed.
!            10 --  solution process must be restarted at the last
!                   orthonormalization point.
!            30 --  solutions are linearly dependent, problem must
!                   be restarted from the beginning.
!     W, P, IP = orthonormalization information.
! **********************************************************************
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DDOT, DMGSBV, DSTOR1, DSTWAY
!***COMMON BLOCKS    DML15T, DML18J, DML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DREORT
!
  DOUBLE PRECISION DDOT
  INTEGER ICOCO, IFLAG, IGOFX, IJK, INDPVT, INFO, INHOMO, INTEG, &
       IP(*), ISTKOP, IVP, J, K, KK, KNSWOT, KOP, L, LOTJP, MFLAG, &
       MNSWOT, MXNON, NCOMP, NCOMPD, NDISK, NEQ, NEQIVP, NFC, &
       NFCC, NFCP, NIC, NIV, NOPG, NPS, NSWOT, NTAPE, NTP, NUMORT, &
       NXPTS
  DOUBLE PRECISION AE, C, DND, DNDT, DX, P(*), PWCND, PX, RE, S(*), &
       SRP, STOWA(*), TND, TOL, VNORM, W(*), WCND, X, XBEG, XEND, &
       XOP, XOT, XSAV, Y(NCOMP,*), YHP(NCOMP,*), YP(*), YPNM
!
!     ******************************************************************
!
  COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFC
  COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
  COMMON /DML18J/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC, &
                  ICOCO
!
! **********************************************************************
!     BEGIN BLOCK PERMITTING ...EXITS TO 210
!        BEGIN BLOCK PERMITTING ...EXITS TO 10
!***FIRST EXECUTABLE STATEMENT  DREORT
        NFCP = NFC + 1
!
!           CHECK TO SEE if ORTHONORMALIZATION TEST IS TO BE PERFORMED
!
!        ...EXIT
        if (IFLAG  /=  1) go to 10
        KNSWOT = KNSWOT + 1
!        ...EXIT
        if (KNSWOT  >=  NSWOT) go to 10
!     ......EXIT
        if ((XEND - X)*(X - XOT)  <  0.0D0) go to 210
   10    CONTINUE
     call DSTOR1(Y,YHP,YP,YHP(1,NFCP),1,0,0)
!
!        ***************************************************************
!
!        ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y
!        AND PARTICULAR SOLUTION YP.
!
     NIV = NFC
     call DMGSBV(NCOMP,NFC,Y,NCOMP,NIV,MFLAG,S,P,IP,INHOMO,YP,W, &
                 WCND)
!
!           ************************************************************
!
!        CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS.
!
     if (MFLAG  ==  0) go to 50
!           BEGIN BLOCK PERMITTING ...EXITS TO 40
           if (IFLAG  ==  2) go to 30
              if (NSWOT  <=  1 .AND. LOTJP  /=  0) go to 20
!
!                    RETRIEVE DATA FOR A RESTART AT LAST
!                    ORTHONORMALIZATION POINT
!
                 call DSTWAY(Y,YP,YHP,1,STOWA)
                 LOTJP = 1
                 NSWOT = 1
                 KNSWOT = 0
                 MNSWOT = MNSWOT/2
                 TND = TND + 1.0D0
                 IFLAG = 10
!           .........EXIT
                 go to 40
   20             CONTINUE
   30          CONTINUE
           IFLAG = 30
   40       CONTINUE
     go to 200
   50    CONTINUE
!           BEGIN BLOCK PERMITTING ...EXITS TO 190
!              BEGIN BLOCK PERMITTING ...EXITS TO 110
!
!                 ******************************************************
!
!              ...EXIT
              if (IFLAG  /=  1) go to 110
!
!                 TEST FOR ORTHONORMALIZATION
!
!              ...EXIT
              if (WCND  <  50.0D0*TOL) go to 110
              DO 60 IJK = 1, NFCP
!              ......EXIT
                 if (S(IJK)  >  1.0D20) go to 110
   60             CONTINUE
!
!                 USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE
!                 NORM DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION
!                 CHECKPOINT.  OTHER CONTROLS ON THE NUMBER OF STEPS TO
!                 THE NEXT CHECKPOINT ARE ADDED FOR SAFETY PURPOSES.
!
              NSWOT = KNSWOT
              KNSWOT = 0
              LOTJP = 0
              WCND = LOG10(WCND)
              if (WCND  >  TND + 3.0D0) NSWOT = 2*NSWOT
              if (WCND  <  PWCND) go to 70
                 XOT = XEND
                 NSWOT = MIN(MNSWOT,NSWOT)
                 PWCND = WCND
                 PX = X
              go to 100
   70             CONTINUE
                 DX = X - PX
                 DND = PWCND - WCND
                 if (DND  >=  4) NSWOT = NSWOT/2
                 DNDT = WCND - TND
                 if (ABS(DX*DNDT)  <=  DND*ABS(XEND-X)) go to 80
                    XOT = XEND
                    NSWOT = MIN(MNSWOT,NSWOT)
                    PWCND = WCND
                    PX = X
                 go to 90
   80                CONTINUE
                    XOT = X + DX*DNDT/DND
                    NSWOT = MIN(MNSWOT,NSWOT)
                    PWCND = WCND
                    PX = X
   90                CONTINUE
  100             CONTINUE
!           ......EXIT
              go to 190
  110          CONTINUE
!
!              *********************************************************
!
!              ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE
!              HOMOGENEOUS SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
!
           NSWOT = 1
           KNSWOT = 0
           LOTJP = 1
           KK = 1
           L = 1
           DO 150 K = 1, NFCC
!                 BEGIN BLOCK PERMITTING ...EXITS TO 140
                 SRP = SQRT(P(KK))
                 if (INHOMO  ==  1) W(K) = SRP*W(K)
                 VNORM = 1.0D0/SRP
                 P(KK) = VNORM
                 KK = KK + NFCC + 1 - K
                 if (NFC  ==  NFCC) go to 120
!                 ......EXIT
                    if (L  /=  K/2) go to 140
  120                CONTINUE
                 DO 130 J = 1, NCOMP
                    Y(J,L) = Y(J,L)*VNORM
  130                CONTINUE
                 L = L + 1
  140             CONTINUE
  150          CONTINUE
!
           if (INHOMO  /=  1 .OR. NPS  ==  1) go to 180
!
!                 NORMALIZE THE PARTICULAR SOLUTION
!
              YPNM = DDOT(NCOMP,YP,1,YP,1)
              if (YPNM  ==  0.0D0) YPNM = 1.0D0
              YPNM = SQRT(YPNM)
              S(NFCP) = YPNM
              DO 160 J = 1, NCOMP
                 YP(J) = YP(J)/YPNM
  160             CONTINUE
              DO 170 J = 1, NFCC
                 W(J) = C*W(J)
  170             CONTINUE
  180          CONTINUE
!
           if (IFLAG  ==  1) call DSTWAY(Y,YP,YHP,0,STOWA)
           IFLAG = 0
  190       CONTINUE
  200    CONTINUE
  210 CONTINUE
  return
end
