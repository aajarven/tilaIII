subroutine BVPOR (Y, NROWY, NCOMP, XPTS, NXPTS, A, NROWA, ALPHA, &
     NIC, B, NROWB, BETA, NFC, IFLAG, Z, MXNON, P, NTP, IP, W, NIV, &
     YHP, U, V, COEF, S, STOWA, G, WORK, IWORK, NFCC)
!
!! BVPOR is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (BVPOR-S, DBVPOR-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!     INPUT to BVPOR    (items not defined in BVSUP comments)
! **********************************************************************
!
!     NOPG = 0 -- Orthonormalization points not pre-assigned
!          = 1 -- Orthonormalization points pre-assigned
!
!     MXNON = Maximum number of orthogonalizations allowed.
!
!     NDISK = 0 -- IN-CORE storage
!           = 1 -- DISK storage.  Value of NTAPE in data statement
!                  is set to 13.  If another value is desired,
!                  the data statement must be changed.
!
!     INTEG = Type of integrator and associated test to be used
!             to determine when to orthonormalize.
!
!             1 -- Use GRAM-SCHMIDT test and DERKF
!             2 -- Use GRAM-SCHMIDT test and DEABM
!
!     TOL = Tolerance for allowable error in orthogonalization test.
!
!     NPS = 0 Normalize particular solution to unit length at each
!             point of orthonormalization.
!         = 1 Do not normalize particular solution.
!
!     NTP = Must be  >=  NFC*(NFC+1)/2.
!
!
!     NFCC = 2*NFC for special treatment of a complex valued problem
!
!     ICOCO = 0 Skip final computations (superposition coefficients
!               and ,hence, boundary problem solution)
!           = 1 Calculate superposition coefficients and obtain
!               solution to the boundary value problem
!
! **********************************************************************
!     OUTPUT from BVPOR
! **********************************************************************
!
!     Y(NROWY,NXPTS) = Solution at specified output points.
!
!     MXNON = Number of orthonormalizations performed by BVPOR.
!
!     Z(MXNON+1) = Locations of orthonormalizations performed by BVPOR.
!
!     NIV = Number of independent vectors returned from MGSBV. Normally
!        this parameter will be meaningful only when MGSBV returns with
!           MFLAG = 2.
!
! **********************************************************************
!
!     The following variables are in the argument list because of
!     variable dimensioning. In general, they contain no information of
!     use to the user.  The amount of storage set aside by the user must
!     be greater than or equal to that indicated by the dimension
!     statements.   For the DISK storage mode, NON = 0 and KPTS = 1,
!     while for the IN-CORE storage mode, NON = MXNON and KPTS = NXPTS.
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
!     Subroutines used by BVPOR
!         LSSUDS -- Solves an underdetermined system of linear
!                   equations.  This routine is used to get a full
!                   set of initial conditions for integration.
!                   Called by BVPOR
!
!         SVECS -- Obtains starting vectors for special treatment
!                  of complex valued problems , called by BVPOR
!
!         RKFAB -- Routine which conducts integration using DERKF or
!                   DEABM
!
!         STWAY -- Storage for backup capability, called by
!                   BVPOR and REORT
!
!         STOR1 -- Storage at output points, called by BVPOR,
!                  RKFAB, REORT and STWAY.
!
!         SDOT -- Single precision vector inner product routine,
!                   called by BVPOR, SCOEF, LSSUDS, MGSBV,
!                   BKSOL, REORT and PRVEC.
!         ** NOTE **
!         A considerable improvement in speed can be achieved if a
!         machine language version is used for SDOT.
!
!         SCOEF -- Computes the superposition constants from the
!                  boundary conditions at Xfinal.
!
!         BKSOL -- Solves an upper triangular set of linear equations.
!
! **********************************************************************
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  BKSOL, LSSUDS, RKFAB, SCOEF, SDOT, STOR1, STWAY,
!                    SVECS
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
!***END PROLOGUE  BVPOR
!
  DIMENSION Y(NROWY,*),A(NROWA,*),ALPHA(*),B(NROWB,*), &
            BETA(*),P(NTP,*),IP(NFCC,*), &
            U(NCOMP,NFC,*),V(NCOMP,*),W(NFCC,*), &
            COEF(*),Z(*),YHP(NCOMP,*),XPTS(*),S(*), &
            WORK(*),IWORK(*),STOWA(*),G(*)
!
! **********************************************************************
!
  COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMPD,NFCD
  COMMON /ML15TO/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
  COMMON /ML18JR/ AE,RE,TOL,NXPTSD,NICD,NOPG,MXNOND,NDISK,NTAPE, &
                  NEQ,INDPVT,INTEG,NPS,NTPD,NEQIVP,NUMORT,NFCCD, &
                  ICOCO
!
! **********************************************************************
!
!***FIRST EXECUTABLE STATEMENT  BVPOR
  NFCP1 = NFC + 1
  NUMORT = 0
  C = 1.0
!
! **********************************************************************
!     CALCULATE INITIAL CONDITIONS WHICH SATISFY
!                   A*YH(XINITIAL)=0  AND  A*YP(XINITIAL)=ALPHA.
!     WHEN NFC  /=  NFCC LSSUDS DEFINES VALUES YHP IN A MATRIX OF SIZE
!     (NFCC+1)*NCOMP AND ,HENCE, OVERFLOWS THE STORAGE ALLOCATION INTO
!     THE U ARRAY. HOWEVER, THIS IS OKAY SINCE PLENTY OF SPACE IS
!     AVAILABLE IN U AND IT HAS NOT YET BEEN USED.
!
  NDW = NROWA * NCOMP
  KWS = NDW + NIC + 1
  KWD = KWS + NIC
  KWT = KWD + NIC
  KWC = KWT + NIC
  IFLAG = 0
  call LSSUDS(A,YHP(1,NFCC+1),ALPHA,NIC,NCOMP,NROWA,YHP,NCOMP, &
              IFLAG,1,IRA,0,WORK(1),WORK(NDW+1),IWORK,WORK(KWS), &
              WORK(KWD),WORK(KWT),ISFLG,WORK(KWC))
  if (IFLAG  ==  1) go to 3
  IFLAG=-4
  go to 250
    3 if (NFC  /=  NFCC) call SVECS(NCOMP,NFC,YHP,WORK,IWORK, &
                     INHOMO,IFLAG)
  if (IFLAG  ==  1)  go to 5
  IFLAG=-5
  go to 250
!
! **********************************************************************
!     DETERMINE THE NUMBER OF DIFFERENTIAL EQUATIONS TO BE INTEGRATED,
!     INITIALIZE VARIABLES FOR AUXILIARY INITIAL VALUE PROBLEM AND
!     STORE INITIAL CONDITIONS.
!
    5 NEQ = NCOMP * NFC
  if (INHOMO  ==  1)  NEQ = NEQ + NCOMP
  IVP = 0
  if (NEQIVP  ==  0)  go to 10
  IVP = NEQ
  NEQ = NEQ + NEQIVP
  NFCP2 = NFCP1
  if (INHOMO  ==  1)  NFCP2 = NFCP1 + 1
  DO 7 K = 1,NEQIVP
    7 YHP(K,NFCP2) = ALPHA(NIC+K)
   10 call STOR1(U,YHP,V,YHP(1,NFCP1),0,NDISK,NTAPE)
!
! **********************************************************************
!     SET UP DATA FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND
!     SAVE INITIAL CONDITIONS IN CASE A RESTART IS NECESSARY.
!
  NSWOT=1
  KNSWOT=0
  LOTJP=1
  TND=LOG10(10.*TOL)
  PWCND=LOG10(SQRT(TOL))
  X=XBEG
  PX=X
  XOT=XEND
  XOP=X
  KOP=1
  call STWAY(U,V,YHP,0,STOWA)
!
! **********************************************************************
! ******** FORWARD INTEGRATION OF ALL INITIAL VALUE EQUATIONS **********
! **********************************************************************
!
  call RKFAB(NCOMP,XPTS,NXPTS,NFC,IFLAG,Z,MXNON,P,NTP,IP, &
              YHP,NIV,U,V,W,S,STOWA,G,WORK,IWORK,NFCC)
  if (IFLAG  /=  0  .OR.  ICOCO  ==  0)  go to 250
!
! **********************************************************************
! **************** BACKWARD SWEEP TO OBTAIN SOLUTION *******************
! **********************************************************************
!
!     CALCULATE SUPERPOSITION COEFFICIENTS AT XFINAL.
!
!   FOR THE DISK STORAGE VERSION, IT IS NOT NECESSARY TO READ  U  AND  V
!   AT THE LAST OUTPUT POINT, SINCE THE LOCAL COPY OF EACH STILL EXISTS.
!
  KOD = 1
  if (NDISK  ==  0)  KOD = NXPTS
  I1=1+NFCC*NFCC
  I2=I1+NFCC
  call SCOEF(U(1,1,KOD),V(1,KOD),NCOMP,NROWB,NFC,NIC,B,BETA,COEF, &
             INHOMO,RE,AE,WORK,WORK(I1),WORK(I2),IWORK,IFLAG,NFCC)
!
! **********************************************************************
!     CALCULATE SOLUTION AT OUTPUT POINTS BY RECURRING BACKWARDS.
!     AS WE RECUR BACKWARDS FROM XFINAL TO XINITIAL WE MUST CALCULATE
!     NEW SUPERPOSITION COEFFICIENTS EACH TIME WE CROSS A POINT OF
!     ORTHONORMALIZATION.
!
  K = NUMORT
  NCOMP2=NCOMP/2
  IC=1
  if (NFC  /=  NFCC) IC=2
  DO 200 J = 1,NXPTS
  KPTS = NXPTS - J + 1
  KOD = KPTS
  if (NDISK  ==  1)  KOD = 1
  135 if (K  ==  0)  go to 170
  if (XEND > XBEG .AND. XPTS(KPTS) >= Z(K))  go to 170
  if (XEND < XBEG .AND. XPTS(KPTS) <= Z(K))  go to 170
  NON = K
  if (NDISK  ==  0)  go to 136
  NON = 1
  BACKSPACE NTAPE
  READ (NTAPE) (IP(I,1), I = 1,NFCC),(P(I,1), I = 1,NTP)
  BACKSPACE NTAPE
  136 if (INHOMO  /=  1)  go to 150
  if (NDISK  ==  0)  go to 138
  BACKSPACE NTAPE
  READ (NTAPE) (W(I,1), I = 1,NFCC)
  BACKSPACE NTAPE
  138 DO 140 N = 1,NFCC
  140 COEF(N) = COEF(N) - W(N,NON)
  150 call BKSOL(NFCC,P(1,NON),COEF)
  DO 155 M = 1,NFCC
  155 WORK(M) = COEF(M)
  DO 160 M = 1,NFCC
  L = IP(M,NON)
  160 COEF(L) = WORK(M)
  K = K - 1
  go to 135
  170 if (NDISK  ==  0)  go to 175
  BACKSPACE NTAPE
  READ (NTAPE) (V(I,1), I = 1,NCOMP), &
               ((U(I,M,1), I = 1,NCOMP), M = 1,NFC)
  BACKSPACE NTAPE
  175 DO 180 N = 1,NCOMP
  180 Y(N,KPTS) = V(N,KOD) + SDOT(NFC,U(N,1,KOD),NCOMP,COEF,IC)
  if (NFC  ==  NFCC) go to 200
  DO 190 N=1,NCOMP2
  NN=NCOMP2+N
  Y(N,KPTS)=Y(N,KPTS) - SDOT(NFC,U(NN,1,KOD),NCOMP,COEF(2),2)
  190 Y(NN,KPTS)=Y(NN,KPTS) + SDOT(NFC,U(N,1,KOD),NCOMP,COEF(2),2)
  200 CONTINUE
!
! **********************************************************************
!
  250 MXNON = NUMORT
  return
end
