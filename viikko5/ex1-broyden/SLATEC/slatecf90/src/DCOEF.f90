subroutine DCOEF (YH, YP, NCOMP, NROWB, NFC, NIC, B, BETA, COEF, &
     INHOMO, RE, AE, BY, CVEC, WORK, IWORK, IFLAG, NFCC)
!
!! DCOEF is subsidiary to DBVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SCOEF-S, DCOEF-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
! INPUT to DCOEF
! **********************************************************************
!
!     YH = matrix of homogeneous solutions.
!     YP = vector containing particular solution.
!     NCOMP = number of components per solution vector.
!     NROWB = first dimension of B in calling program.
!     NFC = number of base solution vectors.
!     NFCC = 2*NFC for the special treatment of COMPLEX*16 valued
!            equations. Otherwise, NFCC=NFC.
!     NIC = number of specified initial conditions.
!     B = boundary condition matrix at X = XFINAL.
!     BETA = vector of nonhomogeneous boundary conditions at X = XFINAL.
!              1 - nonzero particular solution
!     INHOMO = 2 - zero particular solution
!              3 - eigenvalue problem
!     RE = relative error tolerance.
!     AE = absolute error tolerance.
!     BY = storage space for the matrix  B*YH
!     CVEC = storage space for the vector  BETA-B*YP
!     WORK = double precision array of internal storage. Dimension must
!     be GE
!            NFCC*(NFCC+4)
!     IWORK = integer array of internal storage. Dimension must be GE
!             3+NFCC
!
! **********************************************************************
! OUTPUT from DCOEF
! **********************************************************************
!
!     COEF = array containing superposition constants.
!     IFLAG = indicator of success from DSUDS in solving the
!             boundary equations.
!           = 0 boundary equations are solved.
!           = 1 boundary equations appear to have many solutions.
!           = 2 boundary equations appear to be inconsistent.
!           = 3 for this value of an eigenparameter, the boundary
!               equations have only the zero solution.
!
! **********************************************************************
!
!     Subroutine DCOEF solves for the superposition constants from the
!     linear equations defined by the boundary conditions at X = XFINAL.
!
!                          B*YP + B*YH*COEF = BETA
!
! **********************************************************************
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  DDOT, DSUDS, XGETF, XSETF
!***COMMON BLOCKS    DML5MC
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   890921  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  DCOEF
!
  DOUBLE PRECISION DDOT
  INTEGER I, IFLAG, INHOMO, IWORK(*), J, K, KFLAG, KI, L, LPAR, &
       MLSO, NCOMP, NCOMP2, NF, NFC, NFCC, NFCCM1, NIC, &
       NROWB
  DOUBLE PRECISION AE, B(NROWB,*), BBN, BETA(*), BN, BRN, &
       BY(NFCC,*), BYKL, BYS, COEF(*), CONS, CVEC(*), EPS, &
       FOURU, GAM, RE, SQOVFL, SRU, TWOU, UN, URO, WORK(*), &
       YH(NCOMP,*), YP(*), YPN
!
  COMMON /DML5MC/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
!***FIRST EXECUTABLE STATEMENT  DCOEF
!
!     SET UP MATRIX  B*YH  AND VECTOR  BETA - B*YP
!
  NCOMP2 = NCOMP/2
  DO 80 K = 1, NFCC
     DO 10 J = 1, NFC
        L = J
        if (NFC  /=  NFCC) L = 2*J - 1
        BY(K,L) = DDOT(NCOMP,B(K,1),NROWB,YH(1,J),1)
   10    CONTINUE
     if (NFC  ==  NFCC) go to 30
        DO 20 J = 1, NFC
           L = 2*J
           BYKL = DDOT(NCOMP2,B(K,1),NROWB,YH(NCOMP2+1,J),1)
           BY(K,L) = DDOT(NCOMP2,B(K,NCOMP2+1),NROWB,YH(1,J),1) &
                     - BYKL
   20       CONTINUE
   30    CONTINUE
     go to (40,50,60), INHOMO
!        CASE 1
   40    CONTINUE
        CVEC(K) = BETA(K) - DDOT(NCOMP,B(K,1),NROWB,YP,1)
     go to 70
!        CASE 2
   50    CONTINUE
        CVEC(K) = BETA(K)
     go to 70
!        CASE 3
   60    CONTINUE
        CVEC(K) = 0.0D0
   70    CONTINUE
   80 CONTINUE
  CONS = ABS(CVEC(1))
  BYS = ABS(BY(1,1))
!
!         SOLVE LINEAR SYSTEM
!
  IFLAG = 0
  MLSO = 0
  if (INHOMO  ==  3) MLSO = 1
  KFLAG = 0.5D0 * LOG10(EPS)
  call XGETF(NF)
  call XSETF(0)
   90 CONTINUE
     call DSUDS(BY,COEF,CVEC,NFCC,NFCC,NFCC,KFLAG,MLSO,WORK,IWORK)
     if (KFLAG  /=  3) go to 100
     KFLAG = 1
     IFLAG = 1
  go to 90
  100 CONTINUE
  if (KFLAG  ==  4) IFLAG = 2
  call XSETF(NF)
  if (NFCC  ==  1) go to 180
     if (INHOMO  /=  3) go to 170
        if (IWORK(1)  <  NFCC) go to 140
           IFLAG = 3
           DO 110 K = 1, NFCC
              COEF(K) = 0.0D0
  110          CONTINUE
           COEF(NFCC) = 1.0D0
           NFCCM1 = NFCC - 1
           DO 130 K = 1, NFCCM1
              J = NFCC - K
              L = NFCC - J + 1
              GAM = DDOT(L,BY(J,J),NFCC,COEF(J),1)/(WORK(J)*BY(J,J))
              DO 120 I = J, NFCC
                 COEF(I) = COEF(I) + GAM*BY(J,I)
  120             CONTINUE
  130          CONTINUE
        go to 160
  140       CONTINUE
           DO 150 K = 1, NFCC
              KI = 4*NFCC + K
              COEF(K) = WORK(KI)
  150          CONTINUE
  160       CONTINUE
  170    CONTINUE
  go to 220
  180 CONTINUE
!
!            TESTING FOR EXISTENCE AND UNIQUENESS OF BOUNDARY-VALUE
!            PROBLEM SOLUTION IN A SCALAR CASE
!
     BN = 0.0D0
     UN = 0.0D0
     YPN = 0.0D0
     DO 190 K = 1, NCOMP
        UN = MAX(UN,ABS(YH(K,1)))
        YPN = MAX(YPN,ABS(YP(K)))
        BN = MAX(BN,ABS(B(1,K)))
  190    CONTINUE
     BBN = MAX(BN,ABS(BETA(1)))
     if (BYS  >  10.0D0*(RE*UN + AE)*BN) go to 200
        BRN = BBN/BN*BYS
        if (CONS  >=  0.1D0*BRN .AND. CONS  <=  10.0D0*BRN) &
           IFLAG = 1
        if (CONS  >  10.0D0*BRN) IFLAG = 2
        if (CONS  <=  RE*ABS(BETA(1)) + AE + (RE*YPN + AE)*BN) &
           IFLAG = 1
        if (INHOMO  ==  3) COEF(1) = 1.0D0
     go to 210
  200    CONTINUE
     if (INHOMO  /=  3) go to 210
        IFLAG = 3
        COEF(1) = 1.0D0
  210    CONTINUE
  220 CONTINUE
  return
end
