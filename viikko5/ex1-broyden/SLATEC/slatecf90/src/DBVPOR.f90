subroutine DBVPOR (Y, NROWY, NCOMP, XPTS, NXPTS, A, NROWA, ALPHA, &
     NIC, B, NROWB, BETA, NFC, IFLAG, Z, MXNON, P, NTP, IP, W, NIV, &
     YHP, U, V, COEF, S, STOWA, G, WORK, IWORK, NFCC)
!
!! DBVPOR is subsidiary to DBVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (BVPOR-S, DBVPOR-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!     INPUT to DBVPOR    (items not defined in DBVSUP comments)
! **********************************************************************
!
!     NOPG = 0 -- orthonormalization points not pre-assigned
!          = 1 -- orthonormalization points pre-assigned
!
!     MXNON = maximum number of orthogonalizations allowed.
!
!     NDISK = 0 -- in-core storage
!           = 1 -- disk storage.  Value of NTAPE in data statement
!                  is set to 13.  If another value is desired,
!                  the data statement must be changed.
!
!     INTEG = type of integrator and associated test to be used
!             to determine when to orthonormalize.
!
!             1 -- use GRAM-SCHMIDT test and DDERKF
!             2 -- use GRAM-SCHMIDT test and DDEABM
!
!     TOL = tolerance for allowable error in orthogonalization test.
!
!     NPS = 0 normalize particular solution to unit length at each
!             point of orthonormalization.
!         = 1 do not normalize particular solution.
!
!     NTP = must be  >=  NFC*(NFC+1)/2.
!
!     NFCC = 2*NFC for special treatment of a COMPLEX*16 valued problem
!
!     ICOCO = 0 skip final computations (superposition coefficients
!               and, hence, boundary problem solution)
!           = 1 calculate superposition coefficients and obtain
!               solution to the boundary value problem
!
! **********************************************************************
!     OUTPUT from DBVPOR
! **********************************************************************
!
!     Y(NROWY,NXPTS) = solution at specified output points.
!
!     MXNON = number of orthonormalizations performed by DBVPOR.
!
!     Z(MXNON+1) = locations of orthonormalizations performed by DBVPOR.
!
!     NIV = number of independent vectors returned from DMGSBV. Normally
!           this parameter will be meaningful only when DMGSBV returns
!           with MFLAG = 2.
!
! **********************************************************************
!
!     The following variables are in the argument list because of
!     variable dimensioning.  In general, they contain no information of
!     use to the user.  The amount of storage set aside by the user must
!     be greater than or equal to that indicated by the dimension
!     statements.  For the disk storage mode, NON = 0 and KPTS = 1,
!     while for the in-core storage mode, NON = MXNON and KPTS = NXPTS.
!
!     P(NTP,NON+1)
!     IP(NFCC,NON+1)
!     YHP(NCOMP,NFC+1)  plus an additional column of the length  NEQIVP
!     U(NCOMP,NFC,KPTS)
!     V(NCOMP,KPTS)
!     W(NFCC,NON+1)
!     COEF(NFCC)
!     S(NFC+1)
!     STOWA(NCOMP*(NFC+1)+NEQIVP+1)
!     G(NCOMP)
!     WORK(KKKWS)
!     IWORK(LLLIWS)
!
! **********************************************************************
!     SUBROUTINES used by DBVPOR
!         DLSSUD -- solves an underdetermined system of linear
!                   equations.  This routine is used to get a full
!                   set of initial conditions for integration.
!                   Called by DBVPOR.
!
!         DVECS -- obtains starting vectors for special treatment
!                   of COMPLEX*16 valued problems, called by DBVPOR.
!
!         DRKFAB -- routine which conducts integration using DDERKF or
!                   DDEABM.
!
!         DSTWAY -- storage for backup capability, called by
!                   DBVPOR and DREORT.
!
!         DSTOR1 -- storage at output points, called by DBVPOR,
!                   DRKFAB, DREORT and DSTWAY.
!
!         DDOT -- single precision vector inner product routine,
!                   called by DBVPOR, DCOEF, DLSSUD, DMGSBV,
!                   DBKSOL, DREORT and DPRVEC.
!         ** NOTE **
!         a considerable improvement in speed can be achieved if a
!         machine language version is used for DDOT.
!
!         DCOEF -- computes the superposition constants from the
!                   boundary conditions at XFINAL.
!
!         DBKSOL -- solves an upper triangular set of linear equations.
!
! **********************************************************************
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DBKSOL, DCOEF, DDOT, DLSSUD, DRKFAB, DSTOR1,
!                    DSTWAY, DVECS
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
!***END PROLOGUE  DBVPOR
!
  DOUBLE PRECISION DDOT
  INTEGER I, I1, I2, IC, ICOCO, IFLAG, IGOFX, INDPVT, INFO, INHOMO, &
       INTEG, IRA, ISFLG, ISTKOP, IVP, J, &
       K, KNSWOT, KOD, KOP, KPTS, KWC, KWD, KWS, KWT, L, LOTJP, M, &
       MNSWOT, MXNON, MXNOND, N, NCOMP, NCOMP2, NCOMPD, NDISK, NDW, &
       NEQ, NEQIVP, NFC, NFCC, NFCCD, NFCD, NFCP1, NFCP2, NIC, &
       NICD, NIV, NN, NON, NOPG, NPS, NROWA, NROWB, NROWY, NSWOT, &
       NTAPE, NTP, NTPD, NUMORT, NXPTS, NXPTSD, &
       IP(NFCC,*), IWORK(*)
  DOUBLE PRECISION A(NROWA,*), AE, ALPHA(*), B(NROWB,*), &
       BETA(*), C, COEF(*), G(*), P(NTP,*), PWCND, PX, &
       RE, S(*), STOWA(*), TND, TOL, U(NCOMP,NFC,*), &
       V(NCOMP,*), W(NFCC,*), WORK(*), X, XBEG, XEND, XOP, &
       XOT, XPTS(*), XSAV, Y(NROWY,*), YHP(NCOMP,*), &
       Z(*)
!
!     ******************************************************************
!
  COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFCD
  COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
  COMMON /DML18J/ AE,RE,TOL,NXPTSD,NICD,NOPG,MXNOND,NDISK,NTAPE, &
                  NEQ,INDPVT,INTEG,NPS,NTPD,NEQIVP,NUMORT,NFCCD, &
                  ICOCO
!
!      *****************************************************************
!
!***FIRST EXECUTABLE STATEMENT  DBVPOR
  NFCP1 = NFC + 1
  NUMORT = 0
  C = 1.0D0
!
!     ******************************************************************
!         CALCULATE INITIAL CONDITIONS WHICH SATISFY
!                       A*YH(XINITIAL)=0  AND  A*YP(XINITIAL)=ALPHA.
!         WHEN NFC  /=  NFCC DLSSUD DEFINES VALUES YHP IN A MATRIX OF
!         SIZE (NFCC+1)*NCOMP AND ,HENCE, OVERFLOWS THE STORAGE
!         ALLOCATION INTO THE U ARRAY. HOWEVER, THIS IS OKAY SINCE
!         PLENTY OF SPACE IS AVAILABLE IN U AND IT HAS NOT YET BEEN
!         USED.
!
  NDW = NROWA*NCOMP
  KWS = NDW + NIC + 1
  KWD = KWS + NIC
  KWT = KWD + NIC
  KWC = KWT + NIC
  IFLAG = 0
  call DLSSUD(A,YHP(1,NFCC+1),ALPHA,NIC,NCOMP,NROWA,YHP,NCOMP,IFLAG, &
              1,IRA,0,WORK(1),WORK(NDW+1),IWORK,WORK(KWS),WORK(KWD), &
              WORK(KWT),ISFLG,WORK(KWC))
  if (IFLAG  ==  1) go to 10
     IFLAG = -4
  go to 200
   10 CONTINUE
     if (NFC  /=  NFCC) &
        call DVECS(NCOMP,NFC,YHP,WORK,IWORK,INHOMO,IFLAG)
     if (IFLAG  ==  1) go to 20
        IFLAG = -5
     go to 190
   20    CONTINUE
!
!           ************************************************************
!               DETERMINE THE NUMBER OF DIFFERENTIAL EQUATIONS TO BE
!               INTEGRATED, INITIALIZE VARIABLES FOR AUXILIARY INITIAL
!               VALUE PROBLEM AND STORE INITIAL CONDITIONS.
!
        NEQ = NCOMP*NFC
        if (INHOMO  ==  1) NEQ = NEQ + NCOMP
        IVP = 0
        if (NEQIVP  ==  0) go to 40
           IVP = NEQ
           NEQ = NEQ + NEQIVP
           NFCP2 = NFCP1
           if (INHOMO  ==  1) NFCP2 = NFCP1 + 1
           DO 30 K = 1, NEQIVP
              YHP(K,NFCP2) = ALPHA(NIC+K)
   30          CONTINUE
   40       CONTINUE
        call DSTOR1(U,YHP,V,YHP(1,NFCP1),0,NDISK,NTAPE)
!
!           ************************************************************
!               SET UP DATA FOR THE ORTHONORMALIZATION TESTING PROCEDURE
!               AND SAVE INITIAL CONDITIONS IN CASE A RESTART IS
!               NECESSARY.
!
        NSWOT = 1
        KNSWOT = 0
        LOTJP = 1
        TND = LOG10(10.0D0*TOL)
        PWCND = LOG10(SQRT(TOL))
        X = XBEG
        PX = X
        XOT = XEND
        XOP = X
        KOP = 1
        call DSTWAY(U,V,YHP,0,STOWA)
!
!           ************************************************************
!           ******** FORWARD INTEGRATION OF ALL INITIAL VALUE EQUATIONS
!           **********
!           ************************************************************
!
        call DRKFAB(NCOMP,XPTS,NXPTS,NFC,IFLAG,Z,MXNON,P,NTP,IP,YHP, &
                    NIV,U,V,W,S,STOWA,G,WORK,IWORK,NFCC)
        if (IFLAG  /=  0 .OR. ICOCO  ==  0) go to 180
!
!              *********************************************************
!              **************** BACKWARD SWEEP TO OBTAIN SOLUTION
!              *******************
!              *********************************************************
!
!                  CALCULATE SUPERPOSITION COEFFICIENTS AT XFINAL.
!
!                FOR THE DISK STORAGE VERSION, IT IS NOT NECESSARY TO
!                READ  U  AND  V AT THE LAST OUTPUT POINT, SINCE THE
!                LOCAL COPY OF EACH STILL EXISTS.
!
           KOD = 1
           if (NDISK  ==  0) KOD = NXPTS
           I1 = 1 + NFCC*NFCC
           I2 = I1 + NFCC
           call DCOEF(U(1,1,KOD),V(1,KOD),NCOMP,NROWB,NFC,NIC,B, &
                       BETA,COEF,INHOMO,RE,AE,WORK,WORK(I1), &
                       WORK(I2),IWORK,IFLAG,NFCC)
!
!              *********************************************************
!                  CALCULATE SOLUTION AT OUTPUT POINTS BY RECURRING
!                  BACKWARDS.  AS WE RECUR BACKWARDS FROM XFINAL TO
!                  XINITIAL WE MUST CALCULATE NEW SUPERPOSITION
!                  COEFFICIENTS EACH TIME WE CROSS A POINT OF
!                  ORTHONORMALIZATION.
!
           K = NUMORT
           NCOMP2 = NCOMP/2
           IC = 1
           if (NFC  /=  NFCC) IC = 2
           DO 170 J = 1, NXPTS
              KPTS = NXPTS - J + 1
              KOD = KPTS
              if (NDISK  ==  1) KOD = 1
   50             CONTINUE
!                 ...EXIT
                 if (K  ==  0) go to 120
!                 ...EXIT
                 if (XEND  >  XBEG .AND. XPTS(KPTS)  >=  Z(K)) &
                    go to 120
!                 ...EXIT
                 if (XEND  <  XBEG .AND. XPTS(KPTS)  <=  Z(K)) &
                    go to 120
                 NON = K
                 if (NDISK  ==  0) go to 60
                    NON = 1
                    BACKSPACE NTAPE
                    READ (NTAPE) &
                         (IP(I,1), I = 1, NFCC),(P(I,1), I = 1, NTP)
                    BACKSPACE NTAPE
   60                CONTINUE
                 if (INHOMO  /=  1) go to 90
                    if (NDISK  ==  0) go to 70
                       BACKSPACE NTAPE
                       READ (NTAPE) (W(I,1), I = 1, NFCC)
                       BACKSPACE NTAPE
   70                   CONTINUE
                    DO 80 N = 1, NFCC
                       COEF(N) = COEF(N) - W(N,NON)
   80                   CONTINUE
   90                CONTINUE
                 call DBKSOL(NFCC,P(1,NON),COEF)
                 DO 100 M = 1, NFCC
                    WORK(M) = COEF(M)
  100                CONTINUE
                 DO 110 M = 1, NFCC
                    L = IP(M,NON)
                    COEF(L) = WORK(M)
  110                CONTINUE
                 K = K - 1
              go to 50
  120             CONTINUE
              if (NDISK  ==  0) go to 130
                 BACKSPACE NTAPE
                 READ (NTAPE) &
                      (V(I,1), I = 1, NCOMP), &
                      ((U(I,M,1), I = 1, NCOMP), M = 1, NFC)
                 BACKSPACE NTAPE
  130             CONTINUE
              DO 140 N = 1, NCOMP
                 Y(N,KPTS) = V(N,KOD) &
                             + DDOT(NFC,U(N,1,KOD),NCOMP,COEF,IC)
  140             CONTINUE
              if (NFC  ==  NFCC) go to 160
                 DO 150 N = 1, NCOMP2
                    NN = NCOMP2 + N
                    Y(N,KPTS) = Y(N,KPTS) &
                                - DDOT(NFC,U(NN,1,KOD),NCOMP, &
                                       COEF(2),2)
                    Y(NN,KPTS) = Y(NN,KPTS) &
                                 + DDOT(NFC,U(N,1,KOD),NCOMP, &
                                        COEF(2),2)
  150                CONTINUE
  160             CONTINUE
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
!
!     ******************************************************************
!
  MXNON = NUMORT
  return
end
