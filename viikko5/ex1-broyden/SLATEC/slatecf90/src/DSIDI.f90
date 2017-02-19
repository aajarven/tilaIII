subroutine DSIDI (A, LDA, N, KPVT, DET, INERT, WORK, JOB)
!
!! DSIDI computes the determinant, inertia and inverse of a real symmetric ...
!  matrix using the factors from DSIFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1A, D3B1A
!***TYPE      DOUBLE PRECISION (SSIDI-S, DSIDI-D, CHIDI-C, CSIDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     DSIDI computes the determinant, inertia and inverse
!     of a double precision symmetric matrix using the factors from
!     DSIFA.
!
!     On Entry
!
!        A       DOUBLE PRECISION(LDA,N)
!                the output from DSIFA.
!
!        LDA     INTEGER
!                the leading dimension of the array A.
!
!        N       INTEGER
!                the order of the matrix A.
!
!        KPVT    INTEGER(N)
!                the pivot vector from DSIFA.
!
!        WORK    DOUBLE PRECISION(N)
!                work vector.  Contents destroyed.
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
!        A      contains the upper triangle of the inverse of
!               the original matrix.  The strict lower triangle
!               is never referenced.
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
!        A division by zero may occur if the inverse is requested
!        and  DSICO  has set RCOND  ==  0.0
!        or  DSIFA  has set  INFO  /=  0 .
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
!***END PROLOGUE  DSIDI
  INTEGER LDA,N,JOB
  DOUBLE PRECISION A(LDA,*),WORK(*)
  DOUBLE PRECISION DET(2)
  INTEGER KPVT(*),INERT(3)
!
  DOUBLE PRECISION AKKP1,DDOT,TEMP
  DOUBLE PRECISION TEN,D,T,AK,AKP1
  INTEGER J,JB,K,KM1,KS,KSTEP
  LOGICAL NOINV,NODET,NOERT
!***FIRST EXECUTABLE STATEMENT  DSIDI
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
     DO 130 K = 1, N
        D = A(K,K)
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
              T = ABS(A(K,K+1))
              D = (D/T)*A(K+1,K+1) - T
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
  130    CONTINUE
  140 CONTINUE
!
!     COMPUTE INVERSE(A)
!
  if (NOINV) go to 270
     K = 1
  150    if (K  >  N) go to 260
        KM1 = K - 1
        if (KPVT(K)  <  0) go to 180
!
!              1 BY 1
!
           A(K,K) = 1.0D0/A(K,K)
           if (KM1  <  1) go to 170
              call DCOPY(KM1,A(1,K),1,WORK,1)
              DO 160 J = 1, KM1
                 A(J,K) = DDOT(J,A(1,J),1,WORK,1)
                 call DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
              A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)
  170          CONTINUE
           KSTEP = 1
        go to 220
  180       CONTINUE
!
!              2 BY 2
!
           T = ABS(A(K,K+1))
           AK = A(K,K)/T
           AKP1 = A(K+1,K+1)/T
           AKKP1 = A(K,K+1)/T
           D = T*(AK*AKP1 - 1.0D0)
           A(K,K) = AKP1/D
           A(K+1,K+1) = AK/D
           A(K,K+1) = -AKKP1/D
           if (KM1  <  1) go to 210
              call DCOPY(KM1,A(1,K+1),1,WORK,1)
              DO 190 J = 1, KM1
                 A(J,K+1) = DDOT(J,A(1,J),1,WORK,1)
                 call DAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  190             CONTINUE
              A(K+1,K+1) = A(K+1,K+1) + DDOT(KM1,WORK,1,A(1,K+1),1)
              A(K,K+1) = A(K,K+1) + DDOT(KM1,A(1,K),1,A(1,K+1),1)
              call DCOPY(KM1,A(1,K),1,WORK,1)
              DO 200 J = 1, KM1
                 A(J,K) = DDOT(J,A(1,J),1,WORK,1)
                 call DAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  200             CONTINUE
              A(K,K) = A(K,K) + DDOT(KM1,WORK,1,A(1,K),1)
  210          CONTINUE
           KSTEP = 2
  220       CONTINUE
!
!           SWAP
!
        KS = ABS(KPVT(K))
        if (KS  ==  K) go to 250
           call DSWAP(KS,A(1,KS),1,A(1,K),1)
           DO 230 JB = KS, K
              J = K + KS - JB
              TEMP = A(J,K)
              A(J,K) = A(KS,J)
              A(KS,J) = TEMP
  230          CONTINUE
           if (KSTEP  ==  1) go to 240
              TEMP = A(KS,K+1)
              A(KS,K+1) = A(K,K+1)
              A(K,K+1) = TEMP
  240          CONTINUE
  250       CONTINUE
        K = K + KSTEP
     go to 150
  260    CONTINUE
  270 CONTINUE
  return
end
