subroutine DSPDI (AP, N, KPVT, DET, INERT, WORK, JOB)
!
!! DSPDI computes the determinant, inertia, inverse of a real symmetric ...
!  matrix stored in packed form using the factors from DSPFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1A, D3B1A
!***TYPE      DOUBLE PRECISION (SSPDI-S, DSPDI-D, CHPDI-C, CSPDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             PACKED, SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     DSPDI computes the determinant, inertia and inverse
!     of a double precision symmetric matrix using the factors from
!     DSPFA, where the matrix is stored in packed form.
!
!     On Entry
!
!        AP      DOUBLE PRECISION (N*(N+1)/2)
!                the output from DSPFA.
!
!        N       INTEGER
!                the order of the matrix A.
!
!        KPVT    INTEGER(N)
!                the pivot vector from DSPFA.
!
!        WORK    DOUBLE PRECISION(N)
!                work vector.  Contents ignored.
!
!        JOB     INTEGER
!                JOB has the decimal expansion  ABC  where
!                   if  C  /=  0, the inverse is computed,
!                   if  B  /=  0, the determinant is computed,
!                   if  A  /=  0, the inertia is computed.
!
!                For example, JOB = 111  gives all three.
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
!        DET    DOUBLE PRECISION(2)
!               determinant of original matrix.
!               DETERMINANT = DET(1) * 10.0**DET(2)
!               with 1.0  <=  ABS(DET(1))  <  10.0
!               or DET(1) = 0.0.
!
!        INERT  INTEGER(3)
!               the inertia of the original matrix.
!               INERT(1)  =  number of positive eigenvalues.
!               INERT(2)  =  number of negative eigenvalues.
!               INERT(3)  =  number of zero eigenvalues.
!
!     Error Condition
!
!        A division by zero will occur if the inverse is requested
!        and  DSPCO  has set RCOND  ==  0.0
!        or  DSPFA  has set  INFO  /=  0 .
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DCOPY, DDOT, DSWAP
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891107  Modified routine equivalence list.  (WRB)
!   891107  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DSPDI
  INTEGER N,JOB
  DOUBLE PRECISION AP(*),WORK(*)
  DOUBLE PRECISION DET(2)
  INTEGER KPVT(*),INERT(3)
!
  DOUBLE PRECISION AKKP1,DDOT,TEMP
  DOUBLE PRECISION TEN,D,T,AK,AKP1
  INTEGER IJ,IK,IKP1,IKS,J,JB,JK,JKP1
  INTEGER K,KK,KKP1,KM1,KS,KSJ,KSKP1,KSTEP
  LOGICAL NOINV,NODET,NOERT
!***FIRST EXECUTABLE STATEMENT  DSPDI
  NOINV = MOD(JOB,10)  ==  0
  NODET = MOD(JOB,100)/10  ==  0
  NOERT = MOD(JOB,1000)/100  ==  0
!
  if (NODET .AND. NOERT) go to 140
     if (NOERT) go to 10
        INERT(1) = 0
        INERT(2) = 0
        INERT(3) = 0
   10    CONTINUE
     if (NODET) go to 20
        DET(1) = 1.0D0
        DET(2) = 0.0D0
        TEN = 10.0D0
   20    CONTINUE
     T = 0.0D0
     IK = 0
     DO 130 K = 1, N
        KK = IK + K
        D = AP(KK)
!
!           CHECK if 1 BY 1
!
        if (KPVT(K)  >  0) go to 50
!
!              2 BY 2 BLOCK
!              USE DET (D  S)  =  (D/T * C - T) * T  ,  T = ABS(S)
!                      (S  C)
!              TO AVOID UNDERFLOW/OVERFLOW TROUBLES.
!              TAKE TWO PASSES THROUGH SCALING.  USE  T  FOR FLAG.
!
           if (T  /=  0.0D0) go to 30
              IKP1 = IK + K
              KKP1 = IKP1 + K
              T = ABS(AP(KKP1))
              D = (D/T)*AP(KKP1+1) - T
           go to 40
   30          CONTINUE
              D = T
              T = 0.0D0
   40          CONTINUE
   50       CONTINUE
!
        if (NOERT) go to 60
           if (D  >  0.0D0) INERT(1) = INERT(1) + 1
           if (D  <  0.0D0) INERT(2) = INERT(2) + 1
           if (D  ==  0.0D0) INERT(3) = INERT(3) + 1
   60       CONTINUE
!
        if (NODET) go to 120
           DET(1) = D*DET(1)
           if (DET(1)  ==  0.0D0) go to 110
   70             if (ABS(DET(1))  >=  1.0D0) go to 80
                 DET(1) = TEN*DET(1)
                 DET(2) = DET(2) - 1.0D0
              go to 70
   80             CONTINUE
   90             if (ABS(DET(1))  <  TEN) go to 100
                 DET(1) = DET(1)/TEN
                 DET(2) = DET(2) + 1.0D0
              go to 90
  100             CONTINUE
  110          CONTINUE
  120       CONTINUE
        IK = IK + K
  130    CONTINUE
  140 CONTINUE
!
!     COMPUTE INVERSE(A)
!
  if (NOINV) go to 270
     K = 1
     IK = 0
  150    if (K  >  N) go to 260
        KM1 = K - 1
        KK = IK + K
        IKP1 = IK + K
        KKP1 = IKP1 + K
        if (KPVT(K)  <  0) go to 180
!
!              1 BY 1
!
           AP(KK) = 1.0D0/AP(KK)
           if (KM1  <  1) go to 170
              call DCOPY(KM1,AP(IK+1),1,WORK,1)
              IJ = 0
              DO 160 J = 1, KM1
                 JK = IK + J
                 AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
                 call DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                 IJ = IJ + J
  160             CONTINUE
              AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
  170          CONTINUE
           KSTEP = 1
        go to 220
  180       CONTINUE
!
!              2 BY 2
!
           T = ABS(AP(KKP1))
           AK = AP(KK)/T
           AKP1 = AP(KKP1+1)/T
           AKKP1 = AP(KKP1)/T
           D = T*(AK*AKP1 - 1.0D0)
           AP(KK) = AKP1/D
           AP(KKP1+1) = AK/D
           AP(KKP1) = -AKKP1/D
           if (KM1  <  1) go to 210
              call DCOPY(KM1,AP(IKP1+1),1,WORK,1)
              IJ = 0
              DO 190 J = 1, KM1
                 JKP1 = IKP1 + J
                 AP(JKP1) = DDOT(J,AP(IJ+1),1,WORK,1)
                 call DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IKP1+1),1)
                 IJ = IJ + J
  190             CONTINUE
              AP(KKP1+1) = AP(KKP1+1) &
                           + DDOT(KM1,WORK,1,AP(IKP1+1),1)
              AP(KKP1) = AP(KKP1) &
                         + DDOT(KM1,AP(IK+1),1,AP(IKP1+1),1)
              call DCOPY(KM1,AP(IK+1),1,WORK,1)
              IJ = 0
              DO 200 J = 1, KM1
                 JK = IK + J
                 AP(JK) = DDOT(J,AP(IJ+1),1,WORK,1)
                 call DAXPY(J-1,WORK(J),AP(IJ+1),1,AP(IK+1),1)
                 IJ = IJ + J
  200             CONTINUE
              AP(KK) = AP(KK) + DDOT(KM1,WORK,1,AP(IK+1),1)
  210          CONTINUE
           KSTEP = 2
  220       CONTINUE
!
!           SWAP
!
        KS = ABS(KPVT(K))
        if (KS  ==  K) go to 250
           IKS = (KS*(KS - 1))/2
           call DSWAP(KS,AP(IKS+1),1,AP(IK+1),1)
           KSJ = IK + KS
           DO 230 JB = KS, K
              J = K + KS - JB
              JK = IK + J
              TEMP = AP(JK)
              AP(JK) = AP(KSJ)
              AP(KSJ) = TEMP
              KSJ = KSJ - (J - 1)
  230          CONTINUE
           if (KSTEP  ==  1) go to 240
              KSKP1 = IKP1 + KS
              TEMP = AP(KSKP1)
              AP(KSKP1) = AP(KKP1)
              AP(KKP1) = TEMP
  240          CONTINUE
  250       CONTINUE
        IK = IK + K
        if (KSTEP  ==  2) IK = IK + K + 1
        K = K + KSTEP
     go to 150
  260    CONTINUE
  270 CONTINUE
  return
end
