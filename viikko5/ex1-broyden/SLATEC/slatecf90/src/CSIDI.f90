subroutine CSIDI (A, LDA, N, KPVT, DET, WORK, JOB)
!
!! CSIDI computes the determinant and inverse of a complex symmetric ...
!  matrix using the factors from CSIFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2C1, D3C1
!***TYPE      COMPLEX (SSIDI-S, DSIDI-D, CHIDI-C, CSIDI-C)
!***KEYWORDS  DETERMINANT, INVERSE, LINEAR ALGEBRA, LINPACK, MATRIX,
!             SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     CSIDI computes the determinant and inverse
!     of a complex symmetric matrix using the factors from CSIFA.
!
!     On Entry
!
!        A       COMPLEX(LDA,N)
!                the output from CSIFA.
!
!        LDA     INTEGER
!                the leading dimension of the array A .
!
!        N       INTEGER
!                the order of the matrix A .
!
!        KVPT    INTEGER(N)
!                the pivot vector from CSIFA.
!
!        WORK    COMPLEX(N)
!                work vector.  Contents destroyed.
!
!        JOB     INTEGER
!                JOB has the decimal expansion  AB  where
!                   If  B  /=  0, the inverse is computed,
!                   If  A  /=  0, the determinant is computed,
!
!                For example, JOB = 11  gives both.
!
!     On Return
!
!        Variables not requested by JOB are not used.
!
!        A      contains the upper triangle of the inverse of
!               the original matrix.  The strict lower triangle
!               is never referenced.
!
!        DET    COMPLEX(2)
!               determinant of original matrix.
!               Determinant = DET(1) * 10.0**DET(2)
!               with 1.0  <=  ABS(DET(1))  <  10.0
!               or DET(1) = 0.0.
!
!     Error Condition
!
!        A division by zero may occur if the inverse is requested
!        and  CSICO  has set RCOND  ==  0.0
!        or  CSIFA  has set  INFO  /=  0 .
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
!***END PROLOGUE  CSIDI
  INTEGER LDA,N,JOB
  COMPLEX A(LDA,*),DET(2),WORK(*)
  INTEGER KPVT(*)
!
  COMPLEX AK,AKP1,AKKP1,CDOTU,D,T,TEMP
  REAL TEN
  INTEGER J,JB,K,KM1,KS,KSTEP
  LOGICAL NOINV,NODET
  COMPLEX ZDUM
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
!
!***FIRST EXECUTABLE STATEMENT  CSIDI
  NOINV = MOD(JOB,10)  ==  0
  NODET = MOD(JOB,100)/10  ==  0
!
  if (NODET) go to 100
     DET(1) = (1.0E0,0.0E0)
     DET(2) = (0.0E0,0.0E0)
     TEN = 10.0E0
     T = (0.0E0,0.0E0)
     DO 90 K = 1, N
        D = A(K,K)
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
              T = A(K,K+1)
              D = (D/T)*A(K+1,K+1) - T
           go to 20
   10          CONTINUE
              D = T
              T = (0.0E0,0.0E0)
   20          CONTINUE
   30       CONTINUE
!
        DET(1) = D*DET(1)
        if (CABS1(DET(1))  ==  0.0E0) go to 80
   40          if (CABS1(DET(1))  >=  1.0E0) go to 50
              DET(1) = CMPLX(TEN,0.0E0)*DET(1)
              DET(2) = DET(2) - (1.0E0,0.0E0)
           go to 40
   50          CONTINUE
   60          if (CABS1(DET(1))  <  TEN) go to 70
              DET(1) = DET(1)/CMPLX(TEN,0.0E0)
              DET(2) = DET(2) + (1.0E0,0.0E0)
           go to 60
   70          CONTINUE
   80       CONTINUE
   90    CONTINUE
  100 CONTINUE
!
!     COMPUTE INVERSE(A)
!
  if (NOINV) go to 230
     K = 1
  110    if (K  >  N) go to 220
        KM1 = K - 1
        if (KPVT(K)  <  0) go to 140
!
!              1 BY 1
!
           A(K,K) = (1.0E0,0.0E0)/A(K,K)
           if (KM1  <  1) go to 130
              call CCOPY(KM1,A(1,K),1,WORK,1)
              DO 120 J = 1, KM1
                 A(J,K) = CDOTU(J,A(1,J),1,WORK,1)
                 call CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  120             CONTINUE
              A(K,K) = A(K,K) + CDOTU(KM1,WORK,1,A(1,K),1)
  130          CONTINUE
           KSTEP = 1
        go to 180
  140       CONTINUE
!
!              2 BY 2
!
           T = A(K,K+1)
           AK = A(K,K)/T
           AKP1 = A(K+1,K+1)/T
           AKKP1 = A(K,K+1)/T
           D = T*(AK*AKP1 - (1.0E0,0.0E0))
           A(K,K) = AKP1/D
           A(K+1,K+1) = AK/D
           A(K,K+1) = -AKKP1/D
           if (KM1  <  1) go to 170
              call CCOPY(KM1,A(1,K+1),1,WORK,1)
              DO 150 J = 1, KM1
                 A(J,K+1) = CDOTU(J,A(1,J),1,WORK,1)
                 call CAXPY(J-1,WORK(J),A(1,J),1,A(1,K+1),1)
  150             CONTINUE
              A(K+1,K+1) = A(K+1,K+1) &
                           + CDOTU(KM1,WORK,1,A(1,K+1),1)
              A(K,K+1) = A(K,K+1) + CDOTU(KM1,A(1,K),1,A(1,K+1),1)
              call CCOPY(KM1,A(1,K),1,WORK,1)
              DO 160 J = 1, KM1
                 A(J,K) = CDOTU(J,A(1,J),1,WORK,1)
                 call CAXPY(J-1,WORK(J),A(1,J),1,A(1,K),1)
  160             CONTINUE
              A(K,K) = A(K,K) + CDOTU(KM1,WORK,1,A(1,K),1)
  170          CONTINUE
           KSTEP = 2
  180       CONTINUE
!
!           SWAP
!
        KS = ABS(KPVT(K))
        if (KS  ==  K) go to 210
           call CSWAP(KS,A(1,KS),1,A(1,K),1)
           DO 190 JB = KS, K
              J = K + KS - JB
              TEMP = A(J,K)
              A(J,K) = A(KS,J)
              A(KS,J) = TEMP
  190          CONTINUE
           if (KSTEP  ==  1) go to 200
              TEMP = A(KS,K+1)
              A(KS,K+1) = A(K,K+1)
              A(K,K+1) = TEMP
  200          CONTINUE
  210       CONTINUE
        K = K + KSTEP
     go to 110
  220    CONTINUE
  230 CONTINUE
  return
end
