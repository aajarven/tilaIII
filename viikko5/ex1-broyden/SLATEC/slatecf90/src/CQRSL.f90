subroutine CQRSL (X, LDX, N, K, QRAUX, Y, QY, QTY, B, RSD, XB, &
     JOB, INFO)
!
!! CQRSL applies the output of CQRDC to compute coordinate transformations, ...
!  projections, and least squares solutions.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D9, D2C1
!***TYPE      COMPLEX (SQRSL-S, DQRSL-D, CQRSL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, MATRIX, ORTHOGONAL TRIANGULAR,
!             SOLVE
!***AUTHOR  Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     CQRSL applies the output of CQRDC to compute coordinate
!     transformations, projections, and least squares solutions.
!     For K  <=  MIN(N,P), let XK be the matrix
!
!            XK = (X(JVPT(1)),X(JVPT(2)), ... ,X(JVPT(K)))
!
!     formed from columns JVPT(1), ... ,JVPT(K) of the original
!     N x P matrix X that was input to CQRDC (if no pivoting was
!     done, XK consists of the first K columns of X in their
!     original order).  CQRDC produces a factored unitary matrix Q
!     and an upper triangular matrix R such that
!
!              XK = Q * (R)
!                       (0)
!
!     This information is contained in coded form in the arrays
!     X and QRAUX.
!
!     On Entry
!
!        X      COMPLEX(LDX,P).
!               X contains the output of CQRDC.
!
!        LDX    INTEGER.
!               LDX is the leading dimension of the array X.
!
!        N      INTEGER.
!               N is the number of rows of the matrix XK.  It must
!               have the same value as N in CQRDC.
!
!        K      INTEGER.
!               K is the number of columns of the matrix XK.  K
!               must not be greater than (N,P), where P is the
!               same as in the calling sequence to CQRDC.
!
!        QRAUX  COMPLEX(P).
!               QRAUX contains the auxiliary output from CQRDC.
!
!        Y      COMPLEX(N)
!               Y contains an N-vector that is to be manipulated
!               by CQRSL.
!
!        JOB    INTEGER.
!               JOB specifies what is to be computed.  JOB has
!               the decimal expansion ABCDE, with the following
!               meaning.
!
!                    If A  /=  0, compute QY.
!                    If B,C,D, or E  /=  0, compute QTY.
!                    If C  /=  0, compute B.
!                    If D  /=  0, compute RSD .
!                    If E  /=  0, compute  XB.
!
!               Note that a request to compute B, RSD, or XB
!               automatically triggers the computation of QTY, for
!               which an array must be provided in the calling
!               sequence.
!
!     On Return
!
!        QY     COMPLEX(N).
!               QY contains Q*Y, if its computation has been
!               requested.
!
!        QTY    COMPLEX(N).
!               QTY contains CTRANS(Q)*Y, if its computation has
!               been requested.  Here CTRANS(Q) is the conjugate
!               transpose of the matrix Q.
!
!        B      COMPLEX(K)
!               B contains the solution of the least squares problem
!
!                    minimize NORM2(Y - XK*B),
!
!               if its computation has been requested.  (Note that
!               if pivoting was requested in CQRDC, the J-th
!               component of B will be associated with column JVPT(J)
!               of the original matrix X that was input into CQRDC.)
!
!        RSD    COMPLEX(N).
!               RSD contains the least squares residual Y - XK*B,
!               if its computation has been requested.  RSD is
!               also the orthogonal projection of Y onto the
!               orthogonal complement of the column space of XK.
!
!        XB     COMPLEX(N).
!               XB contains the least squares approximation XK*B,
!               if its computation has been requested.  XB is also
!               the orthogonal projection of Y onto the column space
!               of X.
!
!        INFO   INTEGER.
!               INFO is zero unless the computation of B has
!               been requested and R is exactly singular.  In
!               this case, INFO is the index of the first zero
!               diagonal element of R and B is left unaltered.
!
!     The parameters QY, QTY, B, RSD, and XB are not referenced
!     if their computation is not requested and in this case
!     can be replaced by dummy variables in the calling program.
!     To save storage, the user may in some cases use the same
!     array for different parameters in the calling sequence.  A
!     frequently occurring example is when one wishes to compute
!     any of B, RSD, or XB and does not need Y or QTY.  In this
!     case one may identify Y, QTY, and one of B, RSD, or XB, while
!     providing separate arrays for anything else that is to be
!     computed.  Thus the calling sequence
!
!          call CQRSL(X,LDX,N,K,QRAUX,Y,DUM,Y,B,Y,DUM,110,INFO)
!
!     will result in the computation of B and RSD, with RSD
!     overwriting Y.  More generally, each item in the following
!     list contains groups of permissible identifications for
!     a single calling sequence.
!
!          1. (Y,QTY,B) (RSD) (XB) (QY)
!
!          2. (Y,QTY,RSD) (B) (XB) (QY)
!
!          3. (Y,QTY,XB) (B) (RSD) (QY)
!
!          4. (Y,QY) (QTY,B) (RSD) (XB)
!
!          5. (Y,QY) (QTY,RSD) (B) (XB)
!
!          6. (Y,QY) (QTY,XB) (B) (RSD)
!
!     In any group the value returned in the array allocated to
!     the group corresponds to the last member of the group.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  CAXPY, CCOPY, CDOTC
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CQRSL
  INTEGER LDX,N,K,JOB,INFO
  COMPLEX X(LDX,*),QRAUX(*),Y(*),QY(*),QTY(*),B(*),RSD(*),XB(*)
!
  INTEGER I,J,JJ,JU,KP1
  COMPLEX CDOTC,T,TEMP
  LOGICAL CB,CQY,CQTY,CR,CXB
  COMPLEX ZDUM
  REAL CABS1
  CABS1(ZDUM) = ABS(REAL(ZDUM)) + ABS(AIMAG(ZDUM))
!***FIRST EXECUTABLE STATEMENT  CQRSL
!
!     SET INFO FLAG.
!
  INFO = 0
!
!     DETERMINE WHAT IS TO BE COMPUTED.
!
  CQY = JOB/10000  /=  0
  CQTY = MOD(JOB,10000)  /=  0
  CB = MOD(JOB,1000)/100  /=  0
  CR = MOD(JOB,100)/10  /=  0
  CXB = MOD(JOB,10)  /=  0
  JU = MIN(K,N-1)
!
!     SPECIAL ACTION WHEN N=1.
!
  if (JU  /=  0) go to 40
     if (CQY) QY(1) = Y(1)
     if (CQTY) QTY(1) = Y(1)
     if (CXB) XB(1) = Y(1)
     if (.NOT.CB) go to 30
        if (CABS1(X(1,1))  /=  0.0E0) go to 10
           INFO = 1
        go to 20
   10       CONTINUE
           B(1) = Y(1)/X(1,1)
   20       CONTINUE
   30    CONTINUE
     if (CR) RSD(1) = (0.0E0,0.0E0)
  go to 250
   40 CONTINUE
!
!        SET UP TO COMPUTE QY OR QTY.
!
     if (CQY) call CCOPY(N,Y,1,QY,1)
     if (CQTY) call CCOPY(N,Y,1,QTY,1)
     if (.NOT.CQY) go to 70
!
!           COMPUTE QY.
!
        DO 60 JJ = 1, JU
           J = JU - JJ + 1
           if (CABS1(QRAUX(J))  ==  0.0E0) go to 50
              TEMP = X(J,J)
              X(J,J) = QRAUX(J)
              T = -CDOTC(N-J+1,X(J,J),1,QY(J),1)/X(J,J)
              call CAXPY(N-J+1,T,X(J,J),1,QY(J),1)
              X(J,J) = TEMP
   50          CONTINUE
   60       CONTINUE
   70    CONTINUE
     if (.NOT.CQTY) go to 100
!
!           COMPUTE CTRANS(Q)*Y.
!
        DO 90 J = 1, JU
           if (CABS1(QRAUX(J))  ==  0.0E0) go to 80
              TEMP = X(J,J)
              X(J,J) = QRAUX(J)
              T = -CDOTC(N-J+1,X(J,J),1,QTY(J),1)/X(J,J)
              call CAXPY(N-J+1,T,X(J,J),1,QTY(J),1)
              X(J,J) = TEMP
   80          CONTINUE
   90       CONTINUE
  100    CONTINUE
!
!        SET UP TO COMPUTE B, RSD, OR XB.
!
     if (CB) call CCOPY(K,QTY,1,B,1)
     KP1 = K + 1
     if (CXB) call CCOPY(K,QTY,1,XB,1)
     if (CR .AND. K  <  N) call CCOPY(N-K,QTY(KP1),1,RSD(KP1),1)
     if (.NOT.CXB .OR. KP1  >  N) go to 120
        DO 110 I = KP1, N
           XB(I) = (0.0E0,0.0E0)
  110       CONTINUE
  120    CONTINUE
     if (.NOT.CR) go to 140
        DO 130 I = 1, K
           RSD(I) = (0.0E0,0.0E0)
  130       CONTINUE
  140    CONTINUE
     if (.NOT.CB) go to 190
!
!           COMPUTE B.
!
        DO 170 JJ = 1, K
           J = K - JJ + 1
           if (CABS1(X(J,J))  /=  0.0E0) go to 150
              INFO = J
              go to 180
  150          CONTINUE
           B(J) = B(J)/X(J,J)
           if (J  ==  1) go to 160
              T = -B(J)
              call CAXPY(J-1,T,X(1,J),1,B,1)
  160          CONTINUE
  170       CONTINUE
  180       CONTINUE
  190    CONTINUE
     if (.NOT.CR .AND. .NOT.CXB) go to 240
!
!           COMPUTE RSD OR XB AS REQUIRED.
!
        DO 230 JJ = 1, JU
           J = JU - JJ + 1
           if (CABS1(QRAUX(J))  ==  0.0E0) go to 220
              TEMP = X(J,J)
              X(J,J) = QRAUX(J)
              if (.NOT.CR) go to 200
                 T = -CDOTC(N-J+1,X(J,J),1,RSD(J),1)/X(J,J)
                 call CAXPY(N-J+1,T,X(J,J),1,RSD(J),1)
  200             CONTINUE
              if (.NOT.CXB) go to 210
                 T = -CDOTC(N-J+1,X(J,J),1,XB(J),1)/X(J,J)
                 call CAXPY(N-J+1,T,X(J,J),1,XB(J),1)
  210             CONTINUE
              X(J,J) = TEMP
  220          CONTINUE
  230       CONTINUE
  240    CONTINUE
  250 CONTINUE
  return
end
