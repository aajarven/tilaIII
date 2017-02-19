subroutine SSVDC (X, LDX, N, P, S, E, U, LDU, V, LDV, WORK, JOB, &
     INFO)
!
!! SSVDC performs the singular value decomposition of a rectangular matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D6
!***TYPE      SINGLE PRECISION (SSVDC-S, DSVDC-D, CSVDC-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX,
!             SINGULAR VALUE DECOMPOSITION
!***AUTHOR  Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     SSVDC is a subroutine to reduce a real NxP matrix X by orthogonal
!     transformations U and V to diagonal form.  The elements S(I) are
!     the singular values of X.  The columns of U are the corresponding
!     left singular vectors, and the columns of V the right singular
!     vectors.
!
!     On Entry
!
!         X         REAL(LDX,P), where LDX  >=  N.
!                   X contains the matrix whose singular value
!                   decomposition is to be computed.  X is
!                   destroyed by SSVDC.
!
!         LDX       INTEGER
!                   LDX is the leading dimension of the array X.
!
!         N         INTEGER
!                   N is the number of rows of the matrix X.
!
!         P         INTEGER
!                   P is the number of columns of the matrix X.
!
!         LDU       INTEGER
!                   LDU is the leading dimension of the array U.
!                   (See below).
!
!         LDV       INTEGER
!                   LDV is the leading dimension of the array V.
!                   (See below).
!
!         WORK      REAL(N)
!                   work is a scratch array.
!
!         JOB       INTEGER
!                   JOB controls the computation of the singular
!                   vectors.  It has the decimal expansion AB
!                   with the following meaning
!
!                        A  ==  0  Do not compute the left singular
!                                  vectors.
!                        A  ==  1  Return the N left singular vectors
!                                  in U.
!                        A  >=  2  Return the first MIN(N,P) singular
!                                  vectors in U.
!                        B  ==  0  Do not compute the right singular
!                                  vectors.
!                        B  ==  1  Return the right singular vectors
!                                  in V.
!
!     On Return
!
!         S         REAL(MM), where MM=MIN(N+1,P).
!                   The first MIN(N,P) entries of S contain the
!                   singular values of X arranged in descending
!                   order of magnitude.
!
!         E         REAL(P).
!                   E ordinarily contains zeros.  However, see the
!                   discussion of INFO for exceptions.
!
!         U         REAL(LDU,K), where LDU  >=  N.  If JOBA  ==  1, then
!                                   K  ==  N.  If JOBA  >=  2 , then
!                                   K  ==  MIN(N,P).
!                   U contains the matrix of right singular vectors.
!                   U is not referenced if JOBA  ==  0.  If N  <=  P
!                   or if JOBA  ==  2, then U may be identified with X
!                   in the subroutine call.
!
!         V         REAL(LDV,P), where LDV  >=  P.
!                   V contains the matrix of right singular vectors.
!                   V is not referenced if JOB  ==  0.  If P  <=  N,
!                   then V may be identified with X in the
!                   subroutine call.
!
!         INFO      INTEGER.
!                   the singular values (and their corresponding
!                   singular vectors) S(INFO+1),S(INFO+2),...,S(M)
!                   are correct (here M=MIN(N,P)).  Thus if
!                   INFO  ==  0, all the singular values and their
!                   vectors are correct.  In any event, the matrix
!                   B = TRANS(U)*X*V is the bidiagonal matrix
!                   with the elements of S on its diagonal and the
!                   elements of E on its super-diagonal (TRANS(U)
!                   is the transpose of U).  Thus the singular
!                   values of X and B are the same.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SAXPY, SDOT, SNRM2, SROT, SROTG, SSCAL, SSWAP
!***REVISION HISTORY  (YYMMDD)
!   790319  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SSVDC
  INTEGER LDX,N,P,LDU,LDV,JOB,INFO
  REAL X(LDX,*),S(*),E(*),U(LDU,*),V(LDV,*),WORK(*)
!
!
  INTEGER I,ITER,J,JOBU,K,KASE,KK,L,LL,LLS,LM1,LP1,LS,LU,M,MAXIT, &
          MM,MM1,MP1,NCT,NCTP1,NCU,NRT,NRTP1
  REAL SDOT,T
  REAL B,C,CS,EL,EMM1,F,G,SNRM2,SCALE,SHIFT,SL,SM,SN,SMM1,T1,TEST, &
       ZTEST
  LOGICAL WANTU,WANTV
!***FIRST EXECUTABLE STATEMENT  SSVDC
!
!     SET THE MAXIMUM NUMBER OF ITERATIONS.
!
  MAXIT = 30
!
!     DETERMINE WHAT IS TO BE COMPUTED.
!
  WANTU = .FALSE.
  WANTV = .FALSE.
  JOBU = MOD(JOB,100)/10
  NCU = N
  if (JOBU  >  1) NCU = MIN(N,P)
  if (JOBU  /=  0) WANTU = .TRUE.
  if (MOD(JOB,10)  /=  0) WANTV = .TRUE.
!
!     REDUCE X TO BIDIAGONAL FORM, STORING THE DIAGONAL ELEMENTS
!     IN S AND THE SUPER-DIAGONAL ELEMENTS IN E.
!
  INFO = 0
  NCT = MIN(N-1,P)
  NRT = MAX(0,MIN(P-2,N))
  LU = MAX(NCT,NRT)
  if (LU  <  1) go to 170
  DO 160 L = 1, LU
     LP1 = L + 1
     if (L  >  NCT) go to 20
!
!           COMPUTE THE TRANSFORMATION FOR THE L-TH COLUMN AND
!           PLACE THE L-TH DIAGONAL IN S(L).
!
        S(L) = SNRM2(N-L+1,X(L,L),1)
        if (S(L)  ==  0.0E0) go to 10
           if (X(L,L)  /=  0.0E0) S(L) = SIGN(S(L),X(L,L))
           call SSCAL(N-L+1,1.0E0/S(L),X(L,L),1)
           X(L,L) = 1.0E0 + X(L,L)
   10       CONTINUE
        S(L) = -S(L)
   20    CONTINUE
     if (P  <  LP1) go to 50
     DO 40 J = LP1, P
        if (L  >  NCT) go to 30
        if (S(L)  ==  0.0E0) go to 30
!
!              APPLY THE TRANSFORMATION.
!
           T = -SDOT(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
           call SAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
   30       CONTINUE
!
!           PLACE THE L-TH ROW OF X INTO  E FOR THE
!           SUBSEQUENT CALCULATION OF THE ROW TRANSFORMATION.
!
        E(J) = X(L,J)
   40    CONTINUE
   50    CONTINUE
     if (.NOT.WANTU .OR. L  >  NCT) go to 70
!
!           PLACE THE TRANSFORMATION IN U FOR SUBSEQUENT BACK
!           MULTIPLICATION.
!
        DO 60 I = L, N
           U(I,L) = X(I,L)
   60       CONTINUE
   70    CONTINUE
     if (L  >  NRT) go to 150
!
!           COMPUTE THE L-TH ROW TRANSFORMATION AND PLACE THE
!           L-TH SUPER-DIAGONAL IN E(L).
!
        E(L) = SNRM2(P-L,E(LP1),1)
        if (E(L)  ==  0.0E0) go to 80
           if (E(LP1)  /=  0.0E0) E(L) = SIGN(E(L),E(LP1))
           call SSCAL(P-L,1.0E0/E(L),E(LP1),1)
           E(LP1) = 1.0E0 + E(LP1)
   80       CONTINUE
        E(L) = -E(L)
        if (LP1  >  N .OR. E(L)  ==  0.0E0) go to 120
!
!              APPLY THE TRANSFORMATION.
!
           DO 90 I = LP1, N
              WORK(I) = 0.0E0
   90          CONTINUE
           DO 100 J = LP1, P
              call SAXPY(N-L,E(J),X(LP1,J),1,WORK(LP1),1)
  100          CONTINUE
           DO 110 J = LP1, P
              call SAXPY(N-L,-E(J)/E(LP1),WORK(LP1),1,X(LP1,J),1)
  110          CONTINUE
  120       CONTINUE
        if (.NOT.WANTV) go to 140
!
!              PLACE THE TRANSFORMATION IN V FOR SUBSEQUENT
!              BACK MULTIPLICATION.
!
           DO 130 I = LP1, P
              V(I,L) = E(I)
  130          CONTINUE
  140       CONTINUE
  150    CONTINUE
  160 CONTINUE
  170 CONTINUE
!
!     SET UP THE FINAL BIDIAGONAL MATRIX OR ORDER M.
!
  M = MIN(P,N+1)
  NCTP1 = NCT + 1
  NRTP1 = NRT + 1
  if (NCT  <  P) S(NCTP1) = X(NCTP1,NCTP1)
  if (N  <  M) S(M) = 0.0E0
  if (NRTP1  <  M) E(NRTP1) = X(NRTP1,M)
  E(M) = 0.0E0
!
!     if REQUIRED, GENERATE U.
!
  if (.NOT.WANTU) go to 300
     if (NCU  <  NCTP1) go to 200
     DO 190 J = NCTP1, NCU
        DO 180 I = 1, N
           U(I,J) = 0.0E0
  180       CONTINUE
        U(J,J) = 1.0E0
  190    CONTINUE
  200    CONTINUE
     if (NCT  <  1) go to 290
     DO 280 LL = 1, NCT
        L = NCT - LL + 1
        if (S(L)  ==  0.0E0) go to 250
           LP1 = L + 1
           if (NCU  <  LP1) go to 220
           DO 210 J = LP1, NCU
              T = -SDOT(N-L+1,U(L,L),1,U(L,J),1)/U(L,L)
              call SAXPY(N-L+1,T,U(L,L),1,U(L,J),1)
  210          CONTINUE
  220          CONTINUE
           call SSCAL(N-L+1,-1.0E0,U(L,L),1)
           U(L,L) = 1.0E0 + U(L,L)
           LM1 = L - 1
           if (LM1  <  1) go to 240
           DO 230 I = 1, LM1
              U(I,L) = 0.0E0
  230          CONTINUE
  240          CONTINUE
        go to 270
  250       CONTINUE
           DO 260 I = 1, N
              U(I,L) = 0.0E0
  260          CONTINUE
           U(L,L) = 1.0E0
  270       CONTINUE
  280    CONTINUE
  290    CONTINUE
  300 CONTINUE
!
!     if IT IS REQUIRED, GENERATE V.
!
  if (.NOT.WANTV) go to 350
     DO 340 LL = 1, P
        L = P - LL + 1
        LP1 = L + 1
        if (L  >  NRT) go to 320
        if (E(L)  ==  0.0E0) go to 320
           DO 310 J = LP1, P
              T = -SDOT(P-L,V(LP1,L),1,V(LP1,J),1)/V(LP1,L)
              call SAXPY(P-L,T,V(LP1,L),1,V(LP1,J),1)
  310          CONTINUE
  320       CONTINUE
        DO 330 I = 1, P
           V(I,L) = 0.0E0
  330       CONTINUE
        V(L,L) = 1.0E0
  340    CONTINUE
  350 CONTINUE
!
!     MAIN ITERATION LOOP FOR THE SINGULAR VALUES.
!
  MM = M
  ITER = 0
  360 CONTINUE
!
!        QUIT if ALL THE SINGULAR VALUES HAVE BEEN FOUND.
!
     if (M  ==  0) go to 620
!
!        if TOO MANY ITERATIONS HAVE BEEN PERFORMED, SET
!        FLAG AND RETURN.
!
     if (ITER  <  MAXIT) go to 370
        INFO = M
        go to 620
  370    CONTINUE
!
!        THIS SECTION OF THE PROGRAM INSPECTS FOR
!        NEGLIGIBLE ELEMENTS IN THE S AND E ARRAYS.  ON
!        COMPLETION THE VARIABLES KASE AND L ARE SET AS FOLLOWS.
!
!           KASE = 1     if S(M) AND E(L-1) ARE NEGLIGIBLE AND L < M
!           KASE = 2     if S(L) IS NEGLIGIBLE AND L < M
!           KASE = 3     if E(L-1) IS NEGLIGIBLE, L < M, AND
!                        S(L), ..., S(M) ARE NOT NEGLIGIBLE (QR STEP).
!           KASE = 4     if E(M-1) IS NEGLIGIBLE (CONVERGENCE).
!
     DO 390 LL = 1, M
        L = M - LL
        if (L  ==  0) go to 400
        TEST = ABS(S(L)) + ABS(S(L+1))
        ZTEST = TEST + ABS(E(L))
        if (ZTEST  /=  TEST) go to 380
           E(L) = 0.0E0
           go to 400
  380       CONTINUE
  390    CONTINUE
  400    CONTINUE
     if (L  /=  M - 1) go to 410
        KASE = 4
     go to 480
  410    CONTINUE
        LP1 = L + 1
        MP1 = M + 1
        DO 430 LLS = LP1, MP1
           LS = M - LLS + LP1
           if (LS  ==  L) go to 440
           TEST = 0.0E0
           if (LS  /=  M) TEST = TEST + ABS(E(LS))
           if (LS  /=  L + 1) TEST = TEST + ABS(E(LS-1))
           ZTEST = TEST + ABS(S(LS))
           if (ZTEST  /=  TEST) go to 420
              S(LS) = 0.0E0
              go to 440
  420          CONTINUE
  430       CONTINUE
  440       CONTINUE
        if (LS  /=  L) go to 450
           KASE = 3
        go to 470
  450       CONTINUE
        if (LS  /=  M) go to 460
           KASE = 1
        go to 470
  460       CONTINUE
           KASE = 2
           L = LS
  470       CONTINUE
  480    CONTINUE
     L = L + 1
!
!        PERFORM THE TASK INDICATED BY KASE.
!
     go to (490,520,540,570), KASE
!
!        DEFLATE NEGLIGIBLE S(M).
!
  490    CONTINUE
        MM1 = M - 1
        F = E(M-1)
        E(M-1) = 0.0E0
        DO 510 KK = L, MM1
           K = MM1 - KK + L
           T1 = S(K)
           call SROTG(T1,F,CS,SN)
           S(K) = T1
           if (K  ==  L) go to 500
              F = -SN*E(K-1)
              E(K-1) = CS*E(K-1)
  500          CONTINUE
           if (WANTV) call SROT(P,V(1,K),1,V(1,M),1,CS,SN)
  510       CONTINUE
     go to 610
!
!        SPLIT AT NEGLIGIBLE S(L).
!
  520    CONTINUE
        F = E(L-1)
        E(L-1) = 0.0E0
        DO 530 K = L, M
           T1 = S(K)
           call SROTG(T1,F,CS,SN)
           S(K) = T1
           F = -SN*E(K)
           E(K) = CS*E(K)
           if (WANTU) call SROT(N,U(1,K),1,U(1,L-1),1,CS,SN)
  530       CONTINUE
     go to 610
!
!        PERFORM ONE QR STEP.
!
  540    CONTINUE
!
!           CALCULATE THE SHIFT.
!
        SCALE = MAX(ABS(S(M)),ABS(S(M-1)),ABS(E(M-1)),ABS(S(L)), &
                      ABS(E(L)))
        SM = S(M)/SCALE
        SMM1 = S(M-1)/SCALE
        EMM1 = E(M-1)/SCALE
        SL = S(L)/SCALE
        EL = E(L)/SCALE
        B = ((SMM1 + SM)*(SMM1 - SM) + EMM1**2)/2.0E0
        C = (SM*EMM1)**2
        SHIFT = 0.0E0
        if (B  ==  0.0E0 .AND. C  ==  0.0E0) go to 550
           SHIFT = SQRT(B**2+C)
           if (B  <  0.0E0) SHIFT = -SHIFT
           SHIFT = C/(B + SHIFT)
  550       CONTINUE
        F = (SL + SM)*(SL - SM) - SHIFT
        G = SL*EL
!
!           CHASE ZEROS.
!
        MM1 = M - 1
        DO 560 K = L, MM1
           call SROTG(F,G,CS,SN)
           if (K  /=  L) E(K-1) = F
           F = CS*S(K) + SN*E(K)
           E(K) = CS*E(K) - SN*S(K)
           G = SN*S(K+1)
           S(K+1) = CS*S(K+1)
           if (WANTV) call SROT(P,V(1,K),1,V(1,K+1),1,CS,SN)
           call SROTG(F,G,CS,SN)
           S(K) = F
           F = CS*E(K) + SN*S(K+1)
           S(K+1) = -SN*E(K) + CS*S(K+1)
           G = SN*E(K+1)
           E(K+1) = CS*E(K+1)
           if (WANTU .AND. K  <  N) &
              call SROT(N,U(1,K),1,U(1,K+1),1,CS,SN)
  560       CONTINUE
        E(M-1) = F
        ITER = ITER + 1
     go to 610
!
!        CONVERGENCE.
!
  570    CONTINUE
!
!           MAKE THE SINGULAR VALUE  POSITIVE.
!
        if (S(L)  >=  0.0E0) go to 580
           S(L) = -S(L)
           if (WANTV) call SSCAL(P,-1.0E0,V(1,L),1)
  580       CONTINUE
!
!           ORDER THE SINGULAR VALUE.
!
  590       if (L  ==  MM) go to 600
           if (S(L)  >=  S(L+1)) go to 600
           T = S(L)
           S(L) = S(L+1)
           S(L+1) = T
           if (WANTV .AND. L  <  P) &
              call SSWAP(P,V(1,L),1,V(1,L+1),1)
           if (WANTU .AND. L  <  N) &
              call SSWAP(N,U(1,L),1,U(1,L+1),1)
           L = L + 1
        go to 590
  600       CONTINUE
        ITER = 0
        M = M - 1
  610    CONTINUE
  go to 360
  620 CONTINUE
  return
end
