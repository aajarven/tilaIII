subroutine SCHDC (A, LDA, P, WORK, JPVT, JOB, INFO)
!
!! SCHDC computes the Cholesky decomposition of a positive definite matrix.
!
!  A pivoting option allows the user to estimate the
!            condition number of a positive definite matrix or determine
!            the rank of a positive semidefinite matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B
!***TYPE      SINGLE PRECISION (SCHDC-S, DCHDC-D, CCHDC-C)
!***KEYWORDS  CHOLESKY DECOMPOSITION, LINEAR ALGEBRA, LINPACK, MATRIX,
!             POSITIVE DEFINITE
!***AUTHOR  Dongarra, J., (ANL)
!           Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     SCHDC computes the Cholesky decomposition of a positive definite
!     matrix.  A pivoting option allows the user to estimate the
!     condition of a positive definite matrix or determine the rank
!     of a positive semidefinite matrix.
!
!     On Entry
!
!         A      REAL(LDA,P).
!                A contains the matrix whose decomposition is to
!                be computed.  Only the upper half of A need be stored.
!                The lower part of the array A is not referenced.
!
!         LDA    INTEGER.
!                LDA is the leading dimension of the array A.
!
!         P      INTEGER.
!                P is the order of the matrix.
!
!         WORK   REAL.
!                WORK is a work array.
!
!         JPVT   INTEGER(P).
!                JPVT contains integers that control the selection
!                of the pivot elements, if pivoting has been requested.
!                Each diagonal element A(K,K)
!                is placed in one of three classes according to the
!                value of JPVT(K).
!
!                   If JPVT(K)  >  0, then X(K) is an initial
!                                      element.
!
!                   If JPVT(K)  ==  0, then X(K) is a free element.
!
!                   If JPVT(K)  <  0, then X(K) is a final element.
!
!                Before the decomposition is computed, initial elements
!                are moved by symmetric row and column interchanges to
!                the beginning of the array A and final
!                elements to the end.  Both initial and final elements
!                are frozen in place during the computation and only
!                free elements are moved.  At the K-th stage of the
!                reduction, if A(K,K) is occupied by a free element
!                it is interchanged with the largest free element
!                A(L,L) with L  >=  K.  JPVT is not referenced if
!                JOB  ==  0.
!
!        JOB     INTEGER.
!                JOB is an integer that initiates column pivoting.
!                If JOB  ==  0, no pivoting is done.
!                If JOB  /=  0, pivoting is done.
!
!     On Return
!
!         A      A contains in its upper half the Cholesky factor
!                of the matrix A as it has been permuted by pivoting.
!
!         JPVT   JPVT(J) contains the index of the diagonal element
!                of a that was moved into the J-th position,
!                provided pivoting was requested.
!
!         INFO   contains the index of the last positive diagonal
!                element of the Cholesky factor.
!
!     For positive definite matrices INFO = P is the normal return.
!     For pivoting with positive semidefinite matrices INFO will
!     in general be less than P.  However, INFO may be greater than
!     the rank of A, since rounding error can cause an otherwise zero
!     element to be positive.  Indefinite systems will always cause
!     INFO to be less than P.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SAXPY, SSWAP
!***REVISION HISTORY  (YYMMDD)
!   790319  DATE WRITTEN
!   890313  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SCHDC
  INTEGER LDA,P,JPVT(*),JOB,INFO
  REAL A(LDA,*),WORK(*)
!
  INTEGER PU,PL,PLP1,J,JP,JT,K,KB,KM1,KP1,L,MAXL
  REAL TEMP
  REAL MAXDIA
  LOGICAL SWAPK,NEGK
!***FIRST EXECUTABLE STATEMENT  SCHDC
  PL = 1
  PU = 0
  INFO = P
  if (JOB  ==  0) go to 160
!
!        PIVOTING HAS BEEN REQUESTED. REARRANGE THE
!        THE ELEMENTS ACCORDING TO JPVT.
!
     DO 70 K = 1, P
        SWAPK = JPVT(K)  >  0
        NEGK = JPVT(K)  <  0
        JPVT(K) = K
        if (NEGK) JPVT(K) = -JPVT(K)
        if (.NOT.SWAPK) go to 60
           if (K  ==  PL) go to 50
              call SSWAP(PL-1,A(1,K),1,A(1,PL),1)
              TEMP = A(K,K)
              A(K,K) = A(PL,PL)
              A(PL,PL) = TEMP
              PLP1 = PL + 1
              if (P  <  PLP1) go to 40
              DO 30 J = PLP1, P
                 if (J  >=  K) go to 10
                    TEMP = A(PL,J)
                    A(PL,J) = A(J,K)
                    A(J,K) = TEMP
                 go to 20
   10                CONTINUE
                 if (J  ==  K) go to 20
                    TEMP = A(K,J)
                    A(K,J) = A(PL,J)
                    A(PL,J) = TEMP
   20                CONTINUE
   30             CONTINUE
   40             CONTINUE
              JPVT(K) = JPVT(PL)
              JPVT(PL) = K
   50          CONTINUE
           PL = PL + 1
   60       CONTINUE
   70    CONTINUE
     PU = P
     if (P  <  PL) go to 150
     DO 140 KB = PL, P
        K = P - KB + PL
        if (JPVT(K)  >=  0) go to 130
           JPVT(K) = -JPVT(K)
           if (PU  ==  K) go to 120
              call SSWAP(K-1,A(1,K),1,A(1,PU),1)
              TEMP = A(K,K)
              A(K,K) = A(PU,PU)
              A(PU,PU) = TEMP
              KP1 = K + 1
              if (P  <  KP1) go to 110
              DO 100 J = KP1, P
                 if (J  >=  PU) go to 80
                    TEMP = A(K,J)
                    A(K,J) = A(J,PU)
                    A(J,PU) = TEMP
                 go to 90
   80                CONTINUE
                 if (J  ==  PU) go to 90
                    TEMP = A(K,J)
                    A(K,J) = A(PU,J)
                    A(PU,J) = TEMP
   90                CONTINUE
  100             CONTINUE
  110             CONTINUE
              JT = JPVT(K)
              JPVT(K) = JPVT(PU)
              JPVT(PU) = JT
  120          CONTINUE
           PU = PU - 1
  130       CONTINUE
  140    CONTINUE
  150    CONTINUE
  160 CONTINUE
  DO 270 K = 1, P
!
!        REDUCTION LOOP.
!
     MAXDIA = A(K,K)
     KP1 = K + 1
     MAXL = K
!
!        DETERMINE THE PIVOT ELEMENT.
!
     if (K  <  PL .OR. K  >=  PU) go to 190
        DO 180 L = KP1, PU
           if (A(L,L)  <=  MAXDIA) go to 170
              MAXDIA = A(L,L)
              MAXL = L
  170          CONTINUE
  180       CONTINUE
  190    CONTINUE
!
!        QUIT if THE PIVOT ELEMENT IS NOT POSITIVE.
!
     if (MAXDIA  >  0.0E0) go to 200
        INFO = K - 1
        go to 280
  200    CONTINUE
     if (K  ==  MAXL) go to 210
!
!           START THE PIVOTING AND UPDATE JPVT.
!
        KM1 = K - 1
        call SSWAP(KM1,A(1,K),1,A(1,MAXL),1)
        A(MAXL,MAXL) = A(K,K)
        A(K,K) = MAXDIA
        JP = JPVT(MAXL)
        JPVT(MAXL) = JPVT(K)
        JPVT(K) = JP
  210    CONTINUE
!
!        REDUCTION STEP. PIVOTING IS CONTAINED ACROSS THE ROWS.
!
     WORK(K) = SQRT(A(K,K))
     A(K,K) = WORK(K)
     if (P  <  KP1) go to 260
     DO 250 J = KP1, P
        if (K  ==  MAXL) go to 240
           if (J  >=  MAXL) go to 220
              TEMP = A(K,J)
              A(K,J) = A(J,MAXL)
              A(J,MAXL) = TEMP
           go to 230
  220          CONTINUE
           if (J  ==  MAXL) go to 230
              TEMP = A(K,J)
              A(K,J) = A(MAXL,J)
              A(MAXL,J) = TEMP
  230          CONTINUE
  240       CONTINUE
        A(K,J) = A(K,J)/WORK(K)
        WORK(J) = A(K,J)
        TEMP = -A(K,J)
        call SAXPY(J-K,TEMP,WORK(KP1),1,A(KP1,J),1)
  250    CONTINUE
  260    CONTINUE
  270 CONTINUE
  280 CONTINUE
  return
end
