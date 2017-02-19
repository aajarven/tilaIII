subroutine HSTCS1 (INTL, A, B, M, MBDCND, BDA, BDB, C, D, N, &
     NBDCND, BDC, BDD, ELMBDA, F, IDIMF, PERTRB, IERR1, AM, BM, CM, &
     AN, BN, CN, SNTH, RSQ, WRK)
!
!! HSTCS1 is subsidiary to HSTCSP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (HSTCS1-S)
!***AUTHOR  (UNKNOWN)
!***SEE ALSO  HSTCSP
!***ROUTINES CALLED  BLKTRI
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  HSTCS1
  DIMENSION       BDA(*)     ,BDB(*)     ,BDC(*)     ,BDD(*)     , &
                  F(IDIMF,*) ,AM(*)      ,BM(*)      ,CM(*)      , &
                  AN(*)      ,BN(*)      ,CN(*)      ,SNTH(*)    , &
                  RSQ(*)     ,WRK(*)
!***FIRST EXECUTABLE STATEMENT  HSTCS1
  DTH = (B-A)/M
  DTHSQ = DTH*DTH
  DO 101 I=1,M
     SNTH(I) = SIN(A+(I-0.5)*DTH)
  101 CONTINUE
  DR = (D-C)/N
  DO 102 J=1,N
     RSQ(J) = (C+(J-0.5)*DR)**2
  102 CONTINUE
!
!     MULTIPLY RIGHT SIDE BY R(J)**2
!
  DO 104 J=1,N
     X = RSQ(J)
     DO 103 I=1,M
        F(I,J) = X*F(I,J)
  103    CONTINUE
  104 CONTINUE
!
!      DEFINE COEFFICIENTS AM,BM,CM
!
  X = 1./(2.*COS(DTH/2.))
  DO 105 I=2,M
     AM(I) = (SNTH(I-1)+SNTH(I))*X
     CM(I-1) = AM(I)
  105 CONTINUE
  AM(1) = SIN(A)
  CM(M) = SIN(B)
  DO 106 I=1,M
     X = 1./SNTH(I)
     Y = X/DTHSQ
     AM(I) = AM(I)*Y
     CM(I) = CM(I)*Y
     BM(I) = ELMBDA*X*X-AM(I)-CM(I)
  106 CONTINUE
!
!     DEFINE COEFFICIENTS AN,BN,CN
!
  X = C/DR
  DO 107 J=1,N
     AN(J) = (X+J-1)**2
     CN(J) = (X+J)**2
     BN(J) = -(AN(J)+CN(J))
  107 CONTINUE
  ISW = 1
  NB = NBDCND
  if (C == 0. .AND. NB == 2) NB = 6
!
!     ENTER DATA ON THETA BOUNDARIES
!
  go to (108,108,110,110,112,112,108,110,112),MBDCND
  108 BM(1) = BM(1)-AM(1)
  X = 2.*AM(1)
  DO 109 J=1,N
     F(1,J) = F(1,J)-X*BDA(J)
  109 CONTINUE
  go to 112
  110 BM(1) = BM(1)+AM(1)
  X = DTH*AM(1)
  DO 111 J=1,N
     F(1,J) = F(1,J)+X*BDA(J)
  111 CONTINUE
  112 CONTINUE
  go to (113,115,115,113,113,115,117,117,117),MBDCND
  113 BM(M) = BM(M)-CM(M)
  X = 2.*CM(M)
  DO 114 J=1,N
     F(M,J) = F(M,J)-X*BDB(J)
  114 CONTINUE
  go to 117
  115 BM(M) = BM(M)+CM(M)
  X = DTH*CM(M)
  DO 116 J=1,N
     F(M,J) = F(M,J)-X*BDB(J)
  116 CONTINUE
  117 CONTINUE
!
!     ENTER DATA ON R BOUNDARIES
!
  go to (118,118,120,120,122,122),NB
  118 BN(1) = BN(1)-AN(1)
  X = 2.*AN(1)
  DO 119 I=1,M
     F(I,1) = F(I,1)-X*BDC(I)
  119 CONTINUE
  go to 122
  120 BN(1) = BN(1)+AN(1)
  X = DR*AN(1)
  DO 121 I=1,M
     F(I,1) = F(I,1)+X*BDC(I)
  121 CONTINUE
  122 CONTINUE
  go to (123,125,125,123,123,125),NB
  123 BN(N) = BN(N)-CN(N)
  X = 2.*CN(N)
  DO 124 I=1,M
     F(I,N) = F(I,N)-X*BDD(I)
  124 CONTINUE
  go to 127
  125 BN(N) = BN(N)+CN(N)
  X = DR*CN(N)
  DO 126 I=1,M
     F(I,N) = F(I,N)-X*BDD(I)
  126 CONTINUE
  127 CONTINUE
!
!     CHECK FOR SINGULAR PROBLEM.  if SINGULAR, PERTURB F.
!
  PERTRB = 0.
  go to (137,137,128,137,137,128,137,128,128),MBDCND
  128 go to (137,137,129,137,137,129),NB
  129 if (ELMBDA) 137,131,130
  130 IERR1 = 10
  go to 137
  131 CONTINUE
  ISW = 2
  DO 133 I=1,M
     X = 0.
     DO 132 J=1,N
        X = X+F(I,J)
  132    CONTINUE
     PERTRB = PERTRB+X*SNTH(I)
  133 CONTINUE
  X = 0.
  DO 134 J=1,N
     X = X+RSQ(J)
  134 CONTINUE
  PERTRB = 2.*(PERTRB*SIN(DTH/2.))/(X*(COS(A)-COS(B)))
  DO 136 J=1,N
     X = RSQ(J)*PERTRB
     DO 135 I=1,M
        F(I,J) = F(I,J)-X
  135    CONTINUE
  136 CONTINUE
  137 CONTINUE
  A2 = 0.
  DO 138 I=1,M
     A2 = A2+F(I,1)
  138 CONTINUE
  A2 = A2/RSQ(1)
!
!     INITIALIZE BLKTRI
!
  if (INTL  /=  0) go to 139
  call BLKTRI (0,1,N,AN,BN,CN,1,M,AM,BM,CM,IDIMF,F,IERR1,WRK)
  139 CONTINUE
!
!     call BLKTRI TO SOLVE SYSTEM OF EQUATIONS.
!
  call BLKTRI (1,1,N,AN,BN,CN,1,M,AM,BM,CM,IDIMF,F,IERR1,WRK)
  if (ISW /= 2 .OR. C /= 0. .OR. NBDCND /= 2) go to 143
  A1 = 0.
  A3 = 0.
  DO 140 I=1,M
     A1 = A1+SNTH(I)*F(I,1)
     A3 = A3+SNTH(I)
  140 CONTINUE
  A1 = A1+RSQ(1)*A2/2.
  if (MBDCND  ==  3) &
      A1 = A1+(SIN(B)*BDB(1)-SIN(A)*BDA(1))/(2.*(B-A))
  A1 = A1/A3
  A1 = BDC(1)-A1
  DO 142 I=1,M
     DO 141 J=1,N
        F(I,J) = F(I,J)+A1
  141    CONTINUE
  142 CONTINUE
  143 CONTINUE
  return
end
