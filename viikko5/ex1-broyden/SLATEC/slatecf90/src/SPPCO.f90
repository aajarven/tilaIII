subroutine SPPCO (AP, N, RCOND, Z, INFO)
!
!! SPPCO factors a symmetric positive definite matrix stored in packed form ...
!  and estimates the condition number of the matrix.
!
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2B1B
!***TYPE      SINGLE PRECISION (SPPCO-S, DPPCO-D, CPPCO-C)
!***KEYWORDS  CONDITION NUMBER, LINEAR ALGEBRA, LINPACK,
!             MATRIX FACTORIZATION, PACKED, POSITIVE DEFINITE
!***AUTHOR  Moler, C. B., (U. of New Mexico)
!***DESCRIPTION
!
!     SPPCO factors a real symmetric positive definite matrix
!     stored in packed form
!     and estimates the condition of the matrix.
!
!     If  RCOND  is not needed, SPPFA is slightly faster.
!     To solve  A*X = B , follow SPPCO by SPPSL.
!     To compute  INVERSE(A)*C , follow SPPCO by SPPSL.
!     To compute  DETERMINANT(A) , follow SPPCO by SPPDI.
!     To compute  INVERSE(A) , follow SPPCO by SPPDI.
!
!     On Entry
!
!        AP      REAL (N*(N+1)/2)
!                the packed form of a symmetric matrix  A .  The
!                columns of the upper triangle are stored sequentially
!                in a one-dimensional array of length  N*(N+1)/2 .
!                See comments below for details.
!
!        N       INTEGER
!                the order of the matrix  A .
!
!     On Return
!
!        AP      an upper triangular matrix  R , stored in packed
!                form, so that  A = TRANS(R)*R .
!                If  INFO  /=  0 , the factorization is not complete.
!
!        RCOND   REAL
!                an estimate of the reciprocal condition of  A .
!                For the system  A*X = B , relative perturbations
!                in  A  and  B  of size  EPSILON  may cause
!                relative perturbations in  X  of size  EPSILON/RCOND .
!                If  RCOND  is so small that the logical expression
!                           1.0 + RCOND  ==  1.0
!                is true, then  A  may be singular to working
!                precision.  In particular,  RCOND  is zero  if
!                exact singularity is detected or the estimate
!                underflows.  If INFO  /=  0 , RCOND is unchanged.
!
!        Z       REAL(N)
!                a work vector whose contents are usually unimportant.
!                If  A  is singular to working precision, then  Z  is
!                an approximate null vector in the sense that
!                NORM(A*Z) = RCOND*NORM(A)*NORM(Z) .
!                If  INFO  /=  0 , Z  is unchanged.
!
!        INFO    INTEGER
!                = 0  for normal return.
!                = K  signals an error condition.  The leading minor
!                     of order  K  is not positive definite.
!
!     Packed Storage
!
!          The following program segment will pack the upper
!          triangle of a symmetric matrix.
!
!                K = 0
!                DO 20 J = 1, N
!                   DO 10 I = 1, J
!                      K = K + 1
!                      AP(K) = A(I,J)
!             10    CONTINUE
!             20 CONTINUE
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  SASUM, SAXPY, SDOT, SPPFA, SSCAL
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  SPPCO
  INTEGER N,INFO
  REAL AP(*),Z(*)
  REAL RCOND
!
  REAL SDOT,EK,T,WK,WKM
  REAL ANORM,S,SASUM,SM,YNORM
  INTEGER I,IJ,J,JM1,J1,K,KB,KJ,KK,KP1
!
!     FIND NORM OF A
!
!***FIRST EXECUTABLE STATEMENT  SPPCO
  J1 = 1
  DO 30 J = 1, N
     Z(J) = SASUM(J,AP(J1),1)
     IJ = J1
     J1 = J1 + J
     JM1 = J - 1
     if (JM1  <  1) go to 20
     DO 10 I = 1, JM1
        Z(I) = Z(I) + ABS(AP(IJ))
        IJ = IJ + 1
   10    CONTINUE
   20    CONTINUE
   30 CONTINUE
  ANORM = 0.0E0
  DO 40 J = 1, N
     ANORM = MAX(ANORM,Z(J))
   40 CONTINUE
!
!     FACTOR
!
  call SPPFA(AP,N,INFO)
  if (INFO  /=  0) go to 180
!
!        RCOND = 1/(NORM(A)*(ESTIMATE OF NORM(INVERSE(A)))) .
!        ESTIMATE = NORM(Z)/NORM(Y) WHERE  A*Z = Y  AND  A*Y = E .
!        THE COMPONENTS OF  E  ARE CHOSEN TO CAUSE MAXIMUM LOCAL
!        GROWTH IN THE ELEMENTS OF W  WHERE  TRANS(R)*W = E .
!        THE VECTORS ARE FREQUENTLY RESCALED TO AVOID OVERFLOW.
!
!        SOLVE TRANS(R)*W = E
!
     EK = 1.0E0
     DO 50 J = 1, N
        Z(J) = 0.0E0
   50    CONTINUE
     KK = 0
     DO 110 K = 1, N
        KK = KK + K
        if (Z(K)  /=  0.0E0) EK = SIGN(EK,-Z(K))
        if (ABS(EK-Z(K))  <=  AP(KK)) go to 60
           S = AP(KK)/ABS(EK-Z(K))
           call SSCAL(N,S,Z,1)
           EK = S*EK
   60       CONTINUE
        WK = EK - Z(K)
        WKM = -EK - Z(K)
        S = ABS(WK)
        SM = ABS(WKM)
        WK = WK/AP(KK)
        WKM = WKM/AP(KK)
        KP1 = K + 1
        KJ = KK + K
        if (KP1  >  N) go to 100
           DO 70 J = KP1, N
              SM = SM + ABS(Z(J)+WKM*AP(KJ))
              Z(J) = Z(J) + WK*AP(KJ)
              S = S + ABS(Z(J))
              KJ = KJ + J
   70          CONTINUE
           if (S  >=  SM) go to 90
              T = WKM - WK
              WK = WKM
              KJ = KK + K
              DO 80 J = KP1, N
                 Z(J) = Z(J) + T*AP(KJ)
                 KJ = KJ + J
   80             CONTINUE
   90          CONTINUE
  100       CONTINUE
        Z(K) = WK
  110    CONTINUE
     S = 1.0E0/SASUM(N,Z,1)
     call SSCAL(N,S,Z,1)
!
!        SOLVE R*Y = W
!
     DO 130 KB = 1, N
        K = N + 1 - KB
        if (ABS(Z(K))  <=  AP(KK)) go to 120
           S = AP(KK)/ABS(Z(K))
           call SSCAL(N,S,Z,1)
  120       CONTINUE
        Z(K) = Z(K)/AP(KK)
        KK = KK - K
        T = -Z(K)
        call SAXPY(K-1,T,AP(KK+1),1,Z(1),1)
  130    CONTINUE
     S = 1.0E0/SASUM(N,Z,1)
     call SSCAL(N,S,Z,1)
!
     YNORM = 1.0E0
!
!        SOLVE TRANS(R)*V = Y
!
     DO 150 K = 1, N
        Z(K) = Z(K) - SDOT(K-1,AP(KK+1),1,Z(1),1)
        KK = KK + K
        if (ABS(Z(K))  <=  AP(KK)) go to 140
           S = AP(KK)/ABS(Z(K))
           call SSCAL(N,S,Z,1)
           YNORM = S*YNORM
  140       CONTINUE
        Z(K) = Z(K)/AP(KK)
  150    CONTINUE
     S = 1.0E0/SASUM(N,Z,1)
     call SSCAL(N,S,Z,1)
     YNORM = S*YNORM
!
!        SOLVE R*Z = V
!
     DO 170 KB = 1, N
        K = N + 1 - KB
        if (ABS(Z(K))  <=  AP(KK)) go to 160
           S = AP(KK)/ABS(Z(K))
           call SSCAL(N,S,Z,1)
           YNORM = S*YNORM
  160       CONTINUE
        Z(K) = Z(K)/AP(KK)
        KK = KK - K
        T = -Z(K)
        call SAXPY(K-1,T,AP(KK+1),1,Z(1),1)
  170    CONTINUE
!        MAKE ZNORM = 1.0
     S = 1.0E0/SASUM(N,Z,1)
     call SSCAL(N,S,Z,1)
     YNORM = S*YNORM
!
     if (ANORM  /=  0.0E0) RCOND = YNORM/ANORM
     if (ANORM  ==  0.0E0) RCOND = 0.0E0
  180 CONTINUE
  return
end
