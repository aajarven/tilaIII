subroutine DSPSL (AP, N, KPVT, B)
!
!! DSPSL solves a real symmetric system using the factors obtained from DSPFA.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1A
!***TYPE      DOUBLE PRECISION (SSPSL-S, DSPSL-D, CHPSL-C, CSPSL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, PACKED, SOLVE, SYMMETRIC
!***AUTHOR  Bunch, J., (UCSD)
!***DESCRIPTION
!
!     DSISL solves the double precision symmetric system
!     A * X = B
!     using the factors computed by DSPFA.
!
!     On Entry
!
!        AP      DOUBLE PRECISION(N*(N+1)/2)
!                the output from DSPFA.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!        KPVT    INTEGER(N)
!                the pivot vector from DSPFA.
!
!        B       DOUBLE PRECISION(N)
!                the right hand side vector.
!
!     On Return
!
!        B       the solution vector  X .
!
!     Error Condition
!
!        A division by zero may occur if  DSPCO  has set RCOND  ==  0.0
!        or  DSPFA  has set INFO  /=  0  .
!
!     To compute  INVERSE(A) * C  where  C  is a matrix
!     with  P  columns
!           call DSPFA(AP,N,KPVT,INFO)
!           if (INFO  /=  0) go to ...
!           DO 10 J = 1, P
!              call DSPSL(AP,N,KPVT,C(1,J))
!        10 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
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
!***END PROLOGUE  DSPSL
  INTEGER N,KPVT(*)
  DOUBLE PRECISION AP(*),B(*)
!
  DOUBLE PRECISION AK,AKM1,BK,BKM1,DDOT,DENOM,TEMP
  INTEGER IK,IKM1,IKP1,K,KK,KM1K,KM1KM1,KP
!
!     LOOP BACKWARD APPLYING THE TRANSFORMATIONS AND
!     D INVERSE TO B.
!
!***FIRST EXECUTABLE STATEMENT  DSPSL
  K = N
  IK = (N*(N - 1))/2
   10 if (K  ==  0) go to 80
     KK = IK + K
     if (KPVT(K)  <  0) go to 40
!
!           1 X 1 PIVOT BLOCK.
!
        if (K  ==  1) go to 30
           KP = KPVT(K)
           if (KP  ==  K) go to 20
!
!                 INTERCHANGE.
!
              TEMP = B(K)
              B(K) = B(KP)
              B(KP) = TEMP
   20          CONTINUE
!
!              APPLY THE TRANSFORMATION.
!
           call DAXPY(K-1,B(K),AP(IK+1),1,B(1),1)
   30       CONTINUE
!
!           APPLY D INVERSE.
!
        B(K) = B(K)/AP(KK)
        K = K - 1
        IK = IK - K
     go to 70
   40    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
        IKM1 = IK - (K - 1)
        if (K  ==  2) go to 60
           KP = ABS(KPVT(K))
           if (KP  ==  K - 1) go to 50
!
!                 INTERCHANGE.
!
              TEMP = B(K-1)
              B(K-1) = B(KP)
              B(KP) = TEMP
   50          CONTINUE
!
!              APPLY THE TRANSFORMATION.
!
           call DAXPY(K-2,B(K),AP(IK+1),1,B(1),1)
           call DAXPY(K-2,B(K-1),AP(IKM1+1),1,B(1),1)
   60       CONTINUE
!
!           APPLY D INVERSE.
!
        KM1K = IK + K - 1
        KK = IK + K
        AK = AP(KK)/AP(KM1K)
        KM1KM1 = IKM1 + K - 1
        AKM1 = AP(KM1KM1)/AP(KM1K)
        BK = B(K)/AP(KM1K)
        BKM1 = B(K-1)/AP(KM1K)
        DENOM = AK*AKM1 - 1.0D0
        B(K) = (AKM1*BK - BKM1)/DENOM
        B(K-1) = (AK*BKM1 - BK)/DENOM
        K = K - 2
        IK = IK - (K + 1) - K
   70    CONTINUE
  go to 10
   80 CONTINUE
!
!     LOOP FORWARD APPLYING THE TRANSFORMATIONS.
!
  K = 1
  IK = 0
   90 if (K  >  N) go to 160
     if (KPVT(K)  <  0) go to 120
!
!           1 X 1 PIVOT BLOCK.
!
        if (K  ==  1) go to 110
!
!              APPLY THE TRANSFORMATION.
!
           B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
           KP = KPVT(K)
           if (KP  ==  K) go to 100
!
!                 INTERCHANGE.
!
              TEMP = B(K)
              B(K) = B(KP)
              B(KP) = TEMP
  100          CONTINUE
  110       CONTINUE
        IK = IK + K
        K = K + 1
     go to 150
  120    CONTINUE
!
!           2 X 2 PIVOT BLOCK.
!
        if (K  ==  1) go to 140
!
!              APPLY THE TRANSFORMATION.
!
           B(K) = B(K) + DDOT(K-1,AP(IK+1),1,B(1),1)
           IKP1 = IK + K
           B(K+1) = B(K+1) + DDOT(K-1,AP(IKP1+1),1,B(1),1)
           KP = ABS(KPVT(K))
           if (KP  ==  K) go to 130
!
!                 INTERCHANGE.
!
              TEMP = B(K)
              B(K) = B(KP)
              B(KP) = TEMP
  130          CONTINUE
  140       CONTINUE
        IK = IK + K + K + 1
        K = K + 2
  150    CONTINUE
  go to 90
  160 CONTINUE
  return
end
