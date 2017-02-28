subroutine CQRDC (X, LDX, N, P, QRAUX, JPVT, WORK, JOB)
!
!! CQRDC uses Householder transformations to compute the QR factorization ...
!  of an N by P matrix.  Column pivoting is an option.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D5
!***TYPE      COMPLEX (SQRDC-S, DQRDC-D, CQRDC-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR,
!             QR DECOMPOSITION
!***AUTHOR  Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     CQRDC uses Householder transformations to compute the QR
!     factorization of an N by P matrix X.  Column pivoting
!     based on the 2-norms of the reduced columns may be
!     performed at the users option.
!
!     On Entry
!
!        X       COMPLEX(LDX,P), where LDX  >=  N.
!                X contains the matrix whose decomposition is to be
!                computed.
!
!        LDX     INTEGER.
!                LDX is the leading dimension of the array X.
!
!        N       INTEGER.
!                N is the number of rows of the matrix X.
!
!        P       INTEGER.
!                P is the number of columns of the matrix X.
!
!        JVPT    INTEGER(P).
!                JVPT contains integers that control the selection
!                of the pivot columns.  The K-th column X(K) of X
!                is placed in one of three classes according to the
!                value of JVPT(K).
!
!                   If JVPT(K)  >  0, then X(K) is an initial
!                                      column.
!
!                   If JVPT(K)  ==  0, then X(K) is a free column.
!
!                   If JVPT(K)  <  0, then X(K) is a final column.
!
!                Before the decomposition is computed, initial columns
!                are moved to the beginning of the array X and final
!                columns to the end.  Both initial and final columns
!                are frozen in place during the computation and only
!                free columns are moved.  At the K-th stage of the
!                reduction, if X(K) is occupied by a free column
!                it is interchanged with the free column of largest
!                reduced norm.  JVPT is not referenced if
!                JOB  ==  0.
!
!        WORK    COMPLEX(P).
!                WORK is a work array.  WORK is not referenced if
!                JOB  ==  0.
!
!        JOB     INTEGER.
!                JOB is an integer that initiates column pivoting.
!                If JOB  ==  0, no pivoting is done.
!                If JOB  /=  0, pivoting is done.
!
!     On Return
!
!        X       X contains in its upper triangle the upper
!                triangular matrix R of the QR factorization.
!                Below its diagonal X contains information from
!                which the unitary part of the decomposition
!                can be recovered.  Note that if pivoting has
!                been requested, the decomposition is not that
!                of the original matrix X but that of X
!                with its columns permuted as described by JVPT.
!
!        QRAUX   COMPLEX(P).
!                QRAUX contains further information required to recover
!                the unitary part of the decomposition.
!
!        JVPT    JVPT(K) contains the index of the column of the
!                original matrix that has been interchanged into
!                the K-th column, if pivoting was requested.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CDOTC, CSCAL, CSWAP, SCNRM2
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CQRDC
  INTEGER LDX,N,P,JOB
  INTEGER JPVT(*)
  COMPLEX X(LDX,*),QRAUX(*),WORK(*)
!
  INTEGER J,JP,L,LP1,LUP,MAXJ,PL,PU
  REAL MAXNRM,SCNRM2,TT
  COMPLEX CDOTC,NRMXL,T
  LOGICAL NEGJ,SWAPJ
  COMPLEX CSIGN,ZDUM,ZDUM1,ZDUM2
  REAL CABS1
  CSIGN(ZDUM1,ZDUM2) = ABS(ZDUM1)*(ZDUM2/ABS(ZDUM2))
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
!
!***FIRST EXECUTABLE STATEMENT  CQRDC
  PL = 1
  PU = 0
  if (JOB  ==  0) go to 60
!
!        PIVOTING HAS BEEN REQUESTED.  REARRANGE THE COLUMNS
!        ACCORDING TO JPVT.
!
     DO 20 J = 1, P
        SWAPJ = JPVT(J)  >  0
        NEGJ = JPVT(J)  <  0
        JPVT(J) = J
        if (NEGJ) JPVT(J) = -J
        if (.NOT.SWAPJ) go to 10
           if (J  /=  PL) call CSWAP(N,X(1,PL),1,X(1,J),1)
           JPVT(J) = JPVT(PL)
           JPVT(PL) = J
           PL = PL + 1
   10       CONTINUE
   20    CONTINUE
     PU = P
     DO 50 JJ = 1, P
        J = P - JJ + 1
        if (JPVT(J)  >=  0) go to 40
           JPVT(J) = -JPVT(J)
           if (J  ==  PU) go to 30
              call CSWAP(N,X(1,PU),1,X(1,J),1)
              JP = JPVT(PU)
              JPVT(PU) = JPVT(J)
              JPVT(J) = JP
   30          CONTINUE
           PU = PU - 1
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE
!
!     COMPUTE THE NORMS OF THE FREE COLUMNS.
!
  if (PU  <  PL) go to 80
  DO 70 J = PL, PU
     QRAUX(J) = CMPLX(SCNRM2(N,X(1,J),1),0.0E0)
     WORK(J) = QRAUX(J)
   70 CONTINUE
   80 CONTINUE
!
!     PERFORM THE HOUSEHOLDER REDUCTION OF X.
!
  LUP = MIN(N,P)
  DO 200 L = 1, LUP
     if (L  <  PL .OR. L  >=  PU) go to 120
!
!           LOCATE THE COLUMN OF LARGEST NORM AND BRING IT
!           INTO THE PIVOT POSITION.
!
        MAXNRM = 0.0E0
        MAXJ = L
        DO 100 J = L, PU
           if (REAL(QRAUX(J))  <=  MAXNRM) go to 90
              MAXNRM = REAL(QRAUX(J))
              MAXJ = J
   90          CONTINUE
  100       CONTINUE
        if (MAXJ  ==  L) go to 110
           call CSWAP(N,X(1,L),1,X(1,MAXJ),1)
           QRAUX(MAXJ) = QRAUX(L)
           WORK(MAXJ) = WORK(L)
           JP = JPVT(MAXJ)
           JPVT(MAXJ) = JPVT(L)
           JPVT(L) = JP
  110       CONTINUE
  120    CONTINUE
     QRAUX(L) = (0.0E0,0.0E0)
     if (L  ==  N) go to 190
!
!           COMPUTE THE HOUSEHOLDER TRANSFORMATION FOR COLUMN L.
!
        NRMXL = CMPLX(SCNRM2(N-L+1,X(L,L),1),0.0E0)
        if (CABS1(NRMXL)  ==  0.0E0) go to 180
           if (CABS1(X(L,L))  /=  0.0E0) &
              NRMXL = CSIGN(NRMXL,X(L,L))
           call CSCAL(N-L+1,(1.0E0,0.0E0)/NRMXL,X(L,L),1)
           X(L,L) = (1.0E0,0.0E0) + X(L,L)
!
!              APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS,
!              UPDATING THE NORMS.
!
           LP1 = L + 1
           if (P  <  LP1) go to 170
           DO 160 J = LP1, P
              T = -CDOTC(N-L+1,X(L,L),1,X(L,J),1)/X(L,L)
              call CAXPY(N-L+1,T,X(L,L),1,X(L,J),1)
              if (J  <  PL .OR. J  >  PU) go to 150
              if (CABS1(QRAUX(J))  ==  0.0E0) go to 150
                 TT = 1.0E0 - (ABS(X(L,J))/REAL(QRAUX(J)))**2
                 TT = MAX(TT,0.0E0)
                 T = CMPLX(TT,0.0E0)
                 TT = 1.0E0 &
                      + 0.05E0*TT*(REAL(QRAUX(J))/REAL(WORK(J)))**2
                 if (TT  ==  1.0E0) go to 130
                    QRAUX(J) = QRAUX(J)*SQRT(T)
                 go to 140
  130                CONTINUE
                    QRAUX(J) = CMPLX(SCNRM2(N-L,X(L+1,J),1),0.0E0)
                    WORK(J) = QRAUX(J)
  140                CONTINUE
  150             CONTINUE
  160          CONTINUE
  170          CONTINUE
!
!              SAVE THE TRANSFORMATION.
!
           QRAUX(L) = X(L,L)
           X(L,L) = -NRMXL
  180       CONTINUE
  190    CONTINUE
  200 CONTINUE
  return
end