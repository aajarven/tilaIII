subroutine REORT (NCOMP, Y, YP, YHP, NIV, W, S, P, IP, STOWA, &
     IFLAG)
!
!! REORT is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (REORT-S, DREORT-D)
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
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  MGSBV, SDOT, STOR1, STWAY
!***COMMON BLOCKS    ML15TO, ML18JR, ML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  REORT
!
  DIMENSION Y(NCOMP,*),YP(*),W(*),S(*),P(*),IP(*), &
            STOWA(*),YHP(NCOMP,*)
!
! **********************************************************************
!
  COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFC
  COMMON /ML15TO/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
  COMMON /ML18JR/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC, &
                  ICOCO
!
! **********************************************************************
!***FIRST EXECUTABLE STATEMENT  REORT
  NFCP=NFC+1
!
!     CHECK TO SEE if ORTHONORMALIZATION TEST IS TO BE PERFORMED
!
  if (IFLAG  /=  1) go to 5
  KNSWOT=KNSWOT+1
  if (KNSWOT  >=  NSWOT) go to 5
  if ((XEND-X)*(X-XOT)  <  0.) RETURN
    5 call STOR1(Y,YHP,YP,YHP(1,NFCP),1,0,0)
!
!     ****************************************
!
!     ORTHOGONALIZE THE HOMOGENEOUS SOLUTIONS Y
!     AND PARTICULAR SOLUTION YP.
!
  NIV=NFC
  call MGSBV(NCOMP,NFC,Y,NCOMP,NIV,MFLAG,S,P,IP,INHOMO,YP,W,WCND)
!
!     ****************************************
!
!  CHECK FOR LINEAR DEPENDENCE OF THE SOLUTIONS.
!
  if (MFLAG  ==  0)  go to 25
  if (IFLAG  ==  2) go to 15
  if (NSWOT  >  1  .OR.  LOTJP  ==  0) go to 20
   15 IFLAG=30
  return
!
!     RETRIEVE DATA FOR A RESTART AT LAST ORTHONORMALIZATION POINT
!
   20 call STWAY(Y,YP,YHP,1,STOWA)
  LOTJP=1
  NSWOT=1
  KNSWOT=0
  MNSWOT=MNSWOT/2
  TND=TND+1.
  IFLAG=10
  return
!
!     ****************************************
!
   25 if (IFLAG  /=  1) go to 60
!
!     TEST FOR ORTHONORMALIZATION
!
  if (WCND  <  50.*TOL) go to 60
  DO 30 IJK=1,NFCP
  if (S(IJK)  >  1.0E+20) go to 60
   30 CONTINUE
!
!     USE LINEAR EXTRAPOLATION ON LOGARITHMIC VALUES OF THE NORM
!     DECREMENTS TO DETERMINE NEXT ORTHONORMALIZATION CHECKPOINT.
!     OTHER CONTROLS ON THE NUMBER OF STEPS TO THE NEXT CHECKPOINT
!     ARE ADDED FOR SAFETY PURPOSES.
!
  NSWOT=KNSWOT
  KNSWOT=0
  LOTJP=0
  WCND=LOG10(WCND)
  if (WCND  >  TND+3.) NSWOT=2*NSWOT
  if (WCND  >=  PWCND) go to 40
  DX=X-PX
  DND=PWCND-WCND
  if (DND  >=  4) NSWOT=NSWOT/2
  DNDT=WCND-TND
  if (ABS(DX*DNDT)  >  DND*ABS(XEND-X)) go to 40
  XOT=X+DX*DNDT/DND
  go to 50
   40 XOT=XEND
   50 NSWOT=MIN(MNSWOT,NSWOT)
  PWCND=WCND
  PX=X
  return
!
!     ****************************************
!
!     ORTHONORMALIZATION NECESSARY SO WE NORMALIZE THE HOMOGENEOUS
!     SOLUTION VECTORS AND CHANGE W ACCORDINGLY.
!
   60 NSWOT=1
  KNSWOT=0
  LOTJP=1
  KK = 1
  L=1
  DO 70 K = 1,NFCC
  SRP=SQRT(P(KK))
  if (INHOMO  ==  1) W(K)=SRP*W(K)
  VNORM=1./SRP
  P(KK)=VNORM
  KK = KK + NFCC + 1 - K
  if (NFC  ==  NFCC) go to 63
  if (L  /=  K/2) go to 70
   63 DO 65 J = 1,NCOMP
   65 Y(J,L) = Y(J,L)*VNORM
  L=L+1
   70 CONTINUE
!
  if (INHOMO  /=  1  .OR.  NPS  ==  1)  go to 100
!
!     NORMALIZE THE PARTICULAR SOLUTION
!
  YPNM=SDOT(NCOMP,YP,1,YP,1)
  if (YPNM  ==  0.0)  YPNM = 1.0
  YPNM = SQRT(YPNM)
  S(NFCP) = YPNM
  DO 80 J = 1,NCOMP
   80 YP(J) = YP(J) / YPNM
  DO 90 J = 1,NFCC
   90 W(J) = C * W(J)
!
  100 if (IFLAG  ==  1) call STWAY(Y,YP,YHP,0,STOWA)
  IFLAG=0
  return
end
