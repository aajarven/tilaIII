subroutine CSPDI (AP, N, KPVT, DET, WORK, JOB)
!
!! CSPDI computes the determinant and inverse of a complex symmetric ...
!            matrix stored in packed form using the factors from CSPFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2C1, D3C1
!***TYPE      COMPLEX (SSPDI-S, DSPDI-D, CHPDI-C, CSPDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             PACKED, SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     CSPDI computes the determinant and inverse
!     of a complex symmetric matrix using the factors from CSPFA,
!     where the matrix is stored in packed form.
!
!     On Entry
!
!        AP      COMPLEX (N*(N+1)/2)
!                the output from CSPFA.
!
!        N       INTEGER
!                the order of the matrix A .
!
!        KVPT    INTEGER(N)
!                the pivot vector from CSPFA.
!
!        WORK    COMPLEX(N)
!                work vector.  Contents ignored.
!
!        JOB     INTEGER
!                JOB has the decimal expansion  AB  where
!                   if  B  /=  0, the inverse is computed,
!                   if  A  /=  0, the determinant is computed.
!
!                For example, JOB = 11  gives both.
!
!     On Return
!
!        Variables not requested by JOB are not used.
!
!        AP     contains the upper triangle of the inverse of
!               the original matrix, stored in packed form.
!               The columns of the upper triangle are stored
!               sequentially in a one-dimensional array.
!
!        DET    COMPLEX(2)
!               determinant of original matrix.
!               Determinant = DET(1) * 10.0**DET(2)
!               with 1.0  <=  ABS(DET(1))  <  10.0
!               or DET(1) = 0.0.
!
!     Error Condition
!
!        A division by zero will occur if the inverse is requested
!        and  CSPCO  has set RCOND  ==  0.0
!        or  CSPFA  has set  INFO  /=  0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CCOPY, CDOTU, CSWAP
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891107  Corrected category and modified routine equivalence
!           list.  (WRB)
!   891107  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CSPDI
  INTEGER N,JOB
  COMPLEX AP(*),WORK(*),DET(2)
  INTEGER KPVT(*)
!
  COMPLEX AK,AKKP1,AKP1,CDOTU,D,T,TEMP
  REAL TEN
  INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1
  INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP
  LOGICAL NOINV,NODET
  COMPLEX ZDUM
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
!
!***FIRST EXECUTABLE STATEMENT  CSPDI
  NOINV = MOD(JOB,10)  ==  0
  NODET = MOD(JOB,100)/10  ==  0
!
  if (NODET) go to 110
     DET(1) = (1.0E0,0.0E0)
     DET(2) = (0.0E0,0.0E0)
     TEN = 10.0E0
     T = (0.0E0,0.0E0)
     IK = 0
     DO 100 K = 1, N
        KK = IK + K
        D = AP(KK)
!
!           CHECK if 1 BY 1
!
        if (KPVT(K)  >  0) go to 30
!
!              2 BY 2 BLOCK
!              USE DET (D  T)  =  (D/T * C - T) * T
!                      (T  C)
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
!
           if (CABS1(T)  /=  0.0E0) go to 10
              IKP1 = IK + K
              KKP1 = IKP1 + K
              T = AP(KKP1)
              D = (D/T)*AP(KKP1+1) - T
           go to 20
   10          CONTINUE
              D = T
              T = (0.0E0,0.0E0)
   20          CONTINUE
   30       CONTINUE
!
        if (NODET) go to 90
           DET(1) = D*DET(1)
           if (CABS1(DET(1))  ==  0.0E0) go to 80
   40             if (CABS1(DET(1))  >=  1.0E0) go to 50
                 DET(1) = CMPLX(TEN,0.0E0)*DET(1)
                 DET(2) = DET(2) - (1.0E0,0.0E0)
              go to 40
   50             CONTINUE
   60             if (CABS1(DET(1))  <  TEN) go to 70
                 DET(1) = DET(1)/CMPLX(TEN,0.0E0)
                 DET(2) = DET(2) + (1.0E0,0.0E0)
              go to 60
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
        IK = IK + K
  100    CONTINUE
  110 CONTINUE
!
!     COMPUTE INVERSE(A)
!
  if (NOINV) go to 240
     K = 1
     IK = 0
  120    if (K  >  N) go to 230
        KM1 = K - 1
        KK = IK + K
        IKP1 = IK + K
        if (KPVT(K)  <  0) go to 150
!
!              1 BY 1
!
           AP(KK) = (1.0E0,0.0E0)/AP(KK)
           if (KM1  <  1) go to 140
              call CCOPY(KM1,AP(IK+1),1,WORK,1)
              IJ = 0
              DO 130 J = 1, KM1
                 JK = IK + J
                 AP(JK) = CDOTU(J,AP(IJ+1),1,WORK,1)
                 call CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                 IJ = IJ + J
  130             CONTINUE
              AP(KK) = AP(KK) + CDOTU(KM1,WORK,1,AP(IK+1),1)
  140          CONTINUE
           KSTEP = 1
        go to 190
  150       CONTINUE
!
!              2 BY 2
!
           KKP1 = IKP1 + K
           T = AP(KKP1)
           AK = AP(KK)/T
           AKP1 = AP(KKP1+1)/T
           AKKP1 = AP(KKP1)/T
           D = T*(AK*AKP1 - (1.0E0,0.0E0))
           AP(KK) = AKP1/D
           AP(KKP1+1) = AK/D
           AP(KKP1) = -AKKP1/D
           if (KM1  <  1) go to 180
              call CCOPY(KM1,AP(IKP1+1),1,WORK,1)
              IJ = 0
              DO 160 J = 1, KM1
                 JKP1 = IKP1 + J
                 AP(JKP1) = CDOTU(J,AP(IJ+1),1,WORK,1)
                 call CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)
                 IJ = IJ + J
  160             CONTINUE
              AP(KKP1+1) = AP(KKP1+1) &
                           + CDOTU(KM1,WORK,1,AP(IKP1+1),1)
              AP(KKP1) = AP(KKP1) &
                         + CDOTU(KM1,AP(IK+1),1,AP(IKP1+1),1)
              call CCOPY(KM1,AP(IK+1),1,WORK,1)
              IJ = 0
              DO 170 J = 1, KM1
                 JK = IK + J
                 AP(JK) = CDOTU(J,AP(IJ+1),1,WORK,1)
                 call CAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                 IJ = IJ + J
  170             CONTINUE
              AP(KK) = AP(KK) + CDOTU(KM1,WORK,1,AP(IK+1),1)
  180          CONTINUE
           KSTEP = 2
  190       CONTINUE
!
!           SWAP
!
        KS = ABS(KPVT(K))
        if (KS  ==  K) go to 220
           IKS = (KS*(KS - 1))/2
           call CSWAP(KS,AP(IKS+1),1,AP(IK+1),1)
           KSJ = IK + KS
           DO 200 JB = KS, K
              J = K + KS - JB
              JK = IK + J
              TEMP = AP(JK)
              AP(JK) = AP(KSJ)
              AP(KSJ) = TEMP
              KSJ = KSJ - (J - 1)
  200          CONTINUE
           if (KSTEP  ==  1) go to 210
              KSKP1 = IKP1 + KS
              TEMP = AP(KSKP1)
              AP(KSKP1) = AP(KKP1)
              AP(KKP1) = TEMP
  210          CONTINUE
  220       CONTINUE
        IK = IK + K
        if (KSTEP  ==  2) IK = IK + K + 1
        K = K + KSTEP
     go to 120
  230    CONTINUE
  240 CONTINUE
  return
end
