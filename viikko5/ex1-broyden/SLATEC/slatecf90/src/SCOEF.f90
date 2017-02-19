subroutine SCOEF (YH, YP, NCOMP, NROWB, NFC, NIC, B, BETA, COEF, &
     INHOMO, RE, AE, BY, CVEC, WORK, IWORK, IFLAG, NFCC)
!
!! SCOEF is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (SCOEF-S, DCOEF-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
! INPUT TO SCOEF
! **********************************************************************
!
!     YH = Matrix of homogeneous solutions.
!     YP = Vector containing particular solution.
!     NCOMP = Number of components per solution vector.
!     NROWB = First dimension of B in calling program.
!     NFC = Number of base solution vectors.
!     NFCC = 2*NFC for the special treatment of complex valued
!            equations. Otherwise, NFCC=NFC.
!     NIC = Number of specified initial conditions.
!     B = Boundary condition matrix at X = Xfinal.
!     BETA = Vector of nonhomogeneous boundary conditions at X = Xfinal.
!              1 - Nonzero particular solution
!     INHOMO = 2 - Zero particular solution
!              3 - Eigenvalue problem
!     RE = Relative error tolerance
!     AE = Absolute error tolerance
!     BY = Storage space for the matrix  B*YH
!     CVEC = Storage space for the vector  BETA-B*YP
!     WORK = Real array of internal storage. Dimension must be  >=
!            NFCC*(NFCC+4)
!     IWORK = Integer array of internal storage. Dimension must be  >=
!             3+NFCC
!
! **********************************************************************
! OUTPUT FROM SCOEF
! **********************************************************************
!
!     COEF = Array containing superposition constants.
!     IFLAG = Indicator of success from SUDS in solving the
!             boundary equations
!           = 0 Boundary equations are solved
!           = 1 Boundary equations appear to have many solutions
!           = 2 Boundary equations appear to be inconsistent
!           = 3 For this value of an eigenparameter, the boundary
!               equations have only the zero solution.
!
! **********************************************************************
!
!     Subroutine SCOEF solves for the superposition constants from the
!     linear equations defined by the boundary conditions at X = Xfinal.
!
!                          B*YP + B*YH*COEF = BETA
!
! **********************************************************************
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  SDOT, SUDS, XGETF, XSETF
!***COMMON BLOCKS    ML5MCO
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  SCOEF
!
  DIMENSION YH(NCOMP,*),YP(*),B(NROWB,*),BETA(*), &
            COEF(*),BY(NFCC,*),CVEC(*),WORK(*),IWORK(*)
!
  COMMON /ML5MCO/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
!
!     SET UP MATRIX  B*YH  AND VECTOR  BETA - B*YP
!
!***FIRST EXECUTABLE STATEMENT  SCOEF
  NCOMP2=NCOMP/2
  DO 7 K = 1,NFCC
  DO 1 J = 1,NFC
  L=J
  if (NFC  /=  NFCC) L=2*J-1
    1 BY(K,L) = SDOT(NCOMP,B(K,1),NROWB,YH(1,J),1)
  if (NFC  ==  NFCC) go to 3
  DO 2 J=1,NFC
  L=2*J
  BYKL=SDOT(NCOMP2,B(K,1),NROWB,YH(NCOMP2+1,J),1)
  BY(K,L)=SDOT(NCOMP2,B(K,NCOMP2+1),NROWB,YH(1,J),1) - BYKL
    2 CONTINUE
    3 go to (4,5,6), INHOMO
!     CASE 1
    4 CVEC(K) = BETA(K) - SDOT(NCOMP,B(K,1),NROWB,YP,1)
  go to 7
!     CASE 2
    5 CVEC(K) = BETA(K)
  go to 7
!     CASE 3
    6 CVEC(K) = 0.
    7 CONTINUE
  CONS=ABS(CVEC(1))
  BYS=ABS(BY(1,1))
!
! **********************************************************************
!     SOLVE LINEAR SYSTEM
!
  IFLAG=0
  MLSO=0
  if (INHOMO  ==  3) MLSO=1
  KFLAG = 0.5 * LOG10(EPS)
  call XGETF(NF)
  call XSETF(0)
   10 call SUDS(BY,COEF,CVEC,NFCC,NFCC,NFCC,KFLAG,MLSO,WORK,IWORK)
  if (KFLAG  /=  3) go to 13
  KFLAG=1
  IFLAG=1
  go to 10
   13 if (KFLAG  ==  4) IFLAG=2
  call XSETF(NF)
  if (NFCC  ==  1) go to 25
  if (INHOMO  /=  3) RETURN
  if (IWORK(1)  <  NFCC) go to 17
  IFLAG=3
  DO 14 K=1,NFCC
   14 COEF(K)=0.
  COEF(NFCC)=1.
  NFCCM1=NFCC-1
  DO 15 K=1,NFCCM1
  J=NFCC-K
  L=NFCC-J+1
  GAM=SDOT(L,BY(J,J),NFCC,COEF(J),1)/(WORK(J)*BY(J,J))
  DO 15 I=J,NFCC
   15 COEF(I)=COEF(I)+GAM*BY(J,I)
  return
   17 DO 20 K=1,NFCC
  KI=4*NFCC+K
   20 COEF(K)=WORK(KI)
  return
!
! **********************************************************************
!     TESTING FOR EXISTENCE AND UNIQUENESS OF BOUNDARY-VALUE PROBLEM
!     SOLUTION IN A SCALAR CASE
!
   25 BN = 0.
  UN = 0.
  YPN=0.
  DO 30 K = 1,NCOMP
  UN = MAX(UN,ABS(YH(K,1)))
  YPN=MAX(YPN,ABS(YP(K)))
   30 BN = MAX(BN,ABS(B(1,K)))
  BBN = MAX(BN,ABS(BETA(1)))
  if (BYS  >  10.*(RE*UN + AE)*BN)  go to 35
  BRN = BBN / BN * BYS
  if (CONS  >=  0.1*BRN  .AND.  CONS  <=  10.*BRN) IFLAG=1
  if (CONS  >  10.*BRN) IFLAG=2
  if (CONS   <=   RE*ABS(BETA(1))+AE + (RE*YPN+AE)*BN) IFLAG=1
  if (INHOMO  ==  3) COEF(1)=1.
  return
   35 if (INHOMO  /=  3) RETURN
  IFLAG=3
  COEF(1)=1.
  return
end
