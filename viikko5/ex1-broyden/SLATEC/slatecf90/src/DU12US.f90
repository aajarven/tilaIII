subroutine DU12US (A, MDA, M, N, B, MDB, NB, MODE, KRANK, RNORM, &
     H, W, IR, IC)
!
!! DU12US solves a QR factored system for DULSIA.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DULSIA
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (U12US-S, DU12US-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!        Given the Householder LQ factorization of A, this
!        subroutine solves the system AX=B. If the system
!        is of reduced rank, this routine returns a solution
!        according to the selected mode.
!
!       Note - If MODE /= 2, W is never accessed.
!
!***SEE ALSO  DULSIA
!***ROUTINES CALLED  DAXPY, DDOT, DNRM2, DSWAP
!***REVISION HISTORY  (YYMMDD)
!   810801  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DU12US
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  DOUBLE PRECISION DDOT,DNRM2
  DIMENSION A(MDA,*),B(MDB,*),RNORM(*),H(*),W(*)
  INTEGER IC(*),IR(*)
!***FIRST EXECUTABLE STATEMENT  DU12US
  K=KRANK
  KP1=K+1
!
!        RANK=0
!
  if ( K > 0) go to 410
  DO 404 JB=1,NB
  RNORM(JB)=DNRM2(M,B(1,JB),1)
  404 CONTINUE
  DO 406 JB=1,NB
  DO 406 I=1,N
  B(I,JB)=0.0D0
  406 CONTINUE
  return
!
!     REORDER B TO REFLECT ROW INTERCHANGES
!
  410 CONTINUE
  I=0
  412 I=I+1
  if ( I == M) go to 418
  J=IR(I)
  if ( J == I) go to 412
  if ( J < 0) go to 412
  IR(I)=-IR(I)
  DO 413 JB=1,NB
  RNORM(JB)=B(I,JB)
  413 CONTINUE
  IJ=I
  414 DO 415 JB=1,NB
  B(IJ,JB)=B(J,JB)
  415 CONTINUE
  IJ=J
  J=IR(IJ)
  IR(IJ)=-IR(IJ)
  if ( J /= I) go to 414
  DO 416 JB=1,NB
  B(IJ,JB)=RNORM(JB)
  416 CONTINUE
  go to 412
  418 CONTINUE
  DO 420 I=1,M
  IR(I)=ABS(IR(I))
  420 CONTINUE
!
!     if A IS OF REDUCED RANK AND MODE=2,
!     APPLY HOUSEHOLDER TRANSFORMATIONS TO B
!
  if ( MODE < 2 .OR. K == M) go to 440
  MMK=M-K
  DO 430 JB=1,NB
  DO 425 J=1,K
  I=KP1-J
  TT=-DDOT(MMK,A(KP1,I),1,B(KP1,JB),1)/W(I)
  TT=TT-B(I,JB)
  call DAXPY(MMK,TT,A(KP1,I),1,B(KP1,JB),1)
  B(I,JB)=B(I,JB)+TT*W(I)
  425 CONTINUE
  430 CONTINUE
!
!     FIND NORMS OF RESIDUAL VECTOR(S)..(BEFORE OVERWRITE B)
!
  440 DO 442 JB=1,NB
  RNORM(JB)=DNRM2((M-K),B(KP1,JB),1)
  442 CONTINUE
!
!     BACK SOLVE LOWER TRIANGULAR L
!
  DO 450 JB=1,NB
  DO 448 I=1,K
  B(I,JB)=B(I,JB)/A(I,I)
  if ( I == K) go to 450
  IP1=I+1
  call DAXPY(K-I,-B(I,JB),A(IP1,I),1,B(IP1,JB),1)
  448 CONTINUE
  450 CONTINUE
!
!
!      TRUNCATED SOLUTION
!
  if ( K == N) go to 462
  DO 460 JB=1,NB
  DO 460 I=KP1,N
  B(I,JB)=0.0D0
  460 CONTINUE
!
!     APPLY HOUSEHOLDER TRANSFORMATIONS TO B
!
  462 DO 470 I=1,K
  J=KP1-I
  TT=A(J,J)
  A(J,J)=H(J)
  DO 465 JB=1,NB
  BB=-DDOT(N-J+1,A(J,J),MDA,B(J,JB),1)/H(J)
  call DAXPY(N-J+1,BB,A(J,J),MDA,B(J,JB),1)
  465 CONTINUE
  A(J,J)=TT
  470 CONTINUE
!
!
!     REORDER B TO REFLECT COLUMN INTERCHANGES
!
  I=0
  482 I=I+1
  if ( I == N) go to 488
  J=IC(I)
  if ( J == I) go to 482
  if ( J < 0) go to 482
  IC(I)=-IC(I)
  484 call DSWAP(NB,B(J,1),MDB,B(I,1),MDB)
  IJ=IC(J)
  IC(J)=-IC(J)
  J=IJ
  if ( J == I) go to 482
  go to 484
  488 CONTINUE
  DO 490 I=1,N
  IC(I)=ABS(IC(I))
  490 CONTINUE
!
!        SOLUTION VECTORS ARE IN FIRST N ROWS OF B(,)
!
  return
end
