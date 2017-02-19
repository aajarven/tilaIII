subroutine QZIT (NM, N, A, B, EPS1, MATZ, Z, IERR)
!
!! QZIT is the second step of the QZ algorithm for generalized eigenproblems.
!  Accepts an upper Hessenberg and an upper
!            triangular matrix and reduces the former to
!            quasi-triangular form while preserving the form of the
!            latter.  Usually preceded by QZHES and followed by QZVAL
!            and QZVEC.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C1B3
!***TYPE      SINGLE PRECISION (QZIT-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is the second step of the QZ algorithm
!     for solving generalized matrix eigenvalue problems,
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART,
!     as modified in technical note NASA TN D-7305(1973) by WARD.
!
!     This subroutine accepts a pair of REAL matrices, one of them
!     in upper Hessenberg form and the other in upper triangular form.
!     It reduces the Hessenberg matrix to quasi-triangular form using
!     orthogonal transformations while maintaining the triangular form
!     of the other matrix.  It is usually preceded by  QZHES  and
!     followed by  QZVAL  and, possibly,  QZVEC.
!
!     On Input
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, A, B, and Z, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrices A and B.  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        A contains a real upper Hessenberg matrix.  A is a two-
!          dimensional REAL array, dimensioned A(NM,N).
!
!        B contains a real upper triangular matrix.  B is a two-
!          dimensional REAL array, dimensioned B(NM,N).
!
!        EPS1 is a tolerance used to determine negligible elements.
!          EPS1 = 0.0 (or negative) may be input, in which case an
!          element will be neglected only if it is less than roundoff
!          error times the norm of its matrix.  If the input EPS1 is
!          positive, then an element will be considered negligible
!          if it is less than EPS1 times the norm of its matrix.  A
!          positive value of EPS1 may result in faster execution,
!          but less accurate results.  EPS1 is a REAL variable.
!
!        MATZ should be set to .TRUE. if the right hand transformations
!          are to be accumulated for later use in computing
!          eigenvectors, and to .FALSE. otherwise.  MATZ is a LOGICAL
!          variable.
!
!        Z contains, if MATZ has been set to .TRUE., the transformation
!          matrix produced in the reduction by  QZHES, if performed, or
!          else the identity matrix.  If MATZ has been set to .FALSE.,
!          Z is not referenced.  Z is a two-dimensional REAL array,
!          dimensioned Z(NM,N).
!
!     On Output
!
!        A has been reduced to quasi-triangular form.  The elements
!          below the first subdiagonal are still zero, and no two
!          consecutive subdiagonal elements are nonzero.
!
!        B is still in upper triangular form, although its elements
!          have been altered.  The location B(N,1) is used to store
!          EPS1 times the norm of B for later use by  QZVAL  and  QZVEC.
!
!        Z contains the product of the right hand transformations
!          (for both steps) if MATZ has been set to .TRUE.
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if neither A(J,J-1) nor A(J-1,J-2) has become
!                     zero after a total of 30*N iterations.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  QZIT
!
  INTEGER I,J,K,L,N,EN,K1,K2,LD,LL,L1,NA,NM,ISH,ITN,ITS,KM1,LM1
  INTEGER ENM2,IERR,LOR1,ENORN
  REAL A(NM,*),B(NM,*),Z(NM,*)
  REAL R,S,T,A1,A2,A3,EP,SH,U1,U2,U3,V1,V2,V3,ANI
  REAL A11,A12,A21,A22,A33,A34,A43,A44,BNI,B11
  REAL B12,B22,B33,B34,B44,EPSA,EPSB,EPS1,ANORM,BNORM
  LOGICAL MATZ,NOTLAS
!
!***FIRST EXECUTABLE STATEMENT  QZIT
  IERR = 0
!     .......... COMPUTE EPSA,EPSB ..........
  ANORM = 0.0E0
  BNORM = 0.0E0
!
  DO 30 I = 1, N
     ANI = 0.0E0
     if (I  /=  1) ANI = ABS(A(I,I-1))
     BNI = 0.0E0
!
     DO 20 J = I, N
        ANI = ANI + ABS(A(I,J))
        BNI = BNI + ABS(B(I,J))
   20    CONTINUE
!
     if (ANI  >  ANORM) ANORM = ANI
     if (BNI  >  BNORM) BNORM = BNI
   30 CONTINUE
!
  if (ANORM  ==  0.0E0) ANORM = 1.0E0
  if (BNORM  ==  0.0E0) BNORM = 1.0E0
  EP = EPS1
  if (EP  >  0.0E0) go to 50
!     .......... COMPUTE ROUNDOFF LEVEL if EPS1 IS ZERO ..........
  EP = 1.0E0
   40 EP = EP / 2.0E0
  if (1.0E0 + EP  >  1.0E0) go to 40
   50 EPSA = EP * ANORM
  EPSB = EP * BNORM
!     .......... REDUCE A TO QUASI-TRIANGULAR FORM, WHILE
!                KEEPING B TRIANGULAR ..........
  LOR1 = 1
  ENORN = N
  EN = N
  ITN = 30*N
!     .......... BEGIN QZ STEP ..........
   60 if (EN  <=  2) go to 1001
  if (.NOT. MATZ) ENORN = EN
  ITS = 0
  NA = EN - 1
  ENM2 = NA - 1
   70 ISH = 2
!     .......... CHECK FOR CONVERGENCE OR REDUCIBILITY.
!                FOR L=EN STEP -1 UNTIL 1 DO -- ..........
  DO 80 LL = 1, EN
     LM1 = EN - LL
     L = LM1 + 1
     if (L  ==  1) go to 95
     if (ABS(A(L,LM1))  <=  EPSA) go to 90
   80 CONTINUE
!
   90 A(L,LM1) = 0.0E0
  if (L  <  NA) go to 95
!     .......... 1-BY-1 OR 2-BY-2 BLOCK ISOLATED ..........
  EN = LM1
  go to 60
!     .......... CHECK FOR SMALL TOP OF B ..........
   95 LD = L
  100 L1 = L + 1
  B11 = B(L,L)
  if (ABS(B11)  >  EPSB) go to 120
  B(L,L) = 0.0E0
  S = ABS(A(L,L)) + ABS(A(L1,L))
  U1 = A(L,L) / S
  U2 = A(L1,L) / S
  R = SIGN(SQRT(U1*U1+U2*U2),U1)
  V1 = -(U1 + R) / R
  V2 = -U2 / R
  U2 = V2 / V1
!
  DO 110 J = L, ENORN
     T = A(L,J) + U2 * A(L1,J)
     A(L,J) = A(L,J) + T * V1
     A(L1,J) = A(L1,J) + T * V2
     T = B(L,J) + U2 * B(L1,J)
     B(L,J) = B(L,J) + T * V1
     B(L1,J) = B(L1,J) + T * V2
  110 CONTINUE
!
  if (L  /=  1) A(L,LM1) = -A(L,LM1)
  LM1 = L
  L = L1
  go to 90
  120 A11 = A(L,L) / B11
  A21 = A(L1,L) / B11
  if (ISH  ==  1) go to 140
!     .......... ITERATION STRATEGY ..........
  if (ITN  ==  0) go to 1000
  if (ITS  ==  10) go to 155
!     .......... DETERMINE TYPE OF SHIFT ..........
  B22 = B(L1,L1)
  if (ABS(B22)  <  EPSB) B22 = EPSB
  B33 = B(NA,NA)
  if (ABS(B33)  <  EPSB) B33 = EPSB
  B44 = B(EN,EN)
  if (ABS(B44)  <  EPSB) B44 = EPSB
  A33 = A(NA,NA) / B33
  A34 = A(NA,EN) / B44
  A43 = A(EN,NA) / B33
  A44 = A(EN,EN) / B44
  B34 = B(NA,EN) / B44
  T = 0.5E0 * (A43 * B34 - A33 - A44)
  R = T * T + A34 * A43 - A33 * A44
  if (R  <  0.0E0) go to 150
!     .......... DETERMINE SINGLE SHIFT ZEROTH COLUMN OF A ..........
  ISH = 1
  R = SQRT(R)
  SH = -T + R
  S = -T - R
  if (ABS(S-A44)  <  ABS(SH-A44)) SH = S
!     .......... LOOK FOR TWO CONSECUTIVE SMALL
!                SUB-DIAGONAL ELEMENTS OF A.
!                FOR L=EN-2 STEP -1 UNTIL LD DO -- ..........
  DO 130 LL = LD, ENM2
     L = ENM2 + LD - LL
     if (L  ==  LD) go to 140
     LM1 = L - 1
     L1 = L + 1
     T = A(L,L)
     if (ABS(B(L,L))  >  EPSB) T = T - SH * B(L,L)
     if (ABS(A(L,LM1))  <=  ABS(T/A(L1,L)) * EPSA) go to 100
  130 CONTINUE
!
  140 A1 = A11 - SH
  A2 = A21
  if (L  /=  LD) A(L,LM1) = -A(L,LM1)
  go to 160
!     .......... DETERMINE DOUBLE SHIFT ZEROTH COLUMN OF A ..........
  150 A12 = A(L,L1) / B22
  A22 = A(L1,L1) / B22
  B12 = B(L,L1) / B22
  A1 = ((A33 - A11) * (A44 - A11) - A34 * A43 + A43 * B34 * A11) &
       / A21 + A12 - A11 * B12
  A2 = (A22 - A11) - A21 * B12 - (A33 - A11) - (A44 - A11) &
       + A43 * B34
  A3 = A(L1+1,L1) / B22
  go to 160
!     .......... AD HOC SHIFT ..........
  155 A1 = 0.0E0
  A2 = 1.0E0
  A3 = 1.1605E0
  160 ITS = ITS + 1
  ITN = ITN - 1
  if (.NOT. MATZ) LOR1 = LD
!     .......... MAIN LOOP ..........
  DO 260 K = L, NA
     NOTLAS = K  /=  NA .AND. ISH  ==  2
     K1 = K + 1
     K2 = K + 2
     KM1 = MAX(K-1,L)
     LL = MIN(EN,K1+ISH)
     if (NOTLAS) go to 190
!     .......... ZERO A(K+1,K-1) ..........
     if (K  ==  L) go to 170
     A1 = A(K,KM1)
     A2 = A(K1,KM1)
  170    S = ABS(A1) + ABS(A2)
     if (S  ==  0.0E0) go to 70
     U1 = A1 / S
     U2 = A2 / S
     R = SIGN(SQRT(U1*U1+U2*U2),U1)
     V1 = -(U1 + R) / R
     V2 = -U2 / R
     U2 = V2 / V1
!
     DO 180 J = KM1, ENORN
        T = A(K,J) + U2 * A(K1,J)
        A(K,J) = A(K,J) + T * V1
        A(K1,J) = A(K1,J) + T * V2
        T = B(K,J) + U2 * B(K1,J)
        B(K,J) = B(K,J) + T * V1
        B(K1,J) = B(K1,J) + T * V2
  180    CONTINUE
!
     if (K  /=  L) A(K1,KM1) = 0.0E0
     go to 240
!     .......... ZERO A(K+1,K-1) AND A(K+2,K-1) ..........
  190    if (K  ==  L) go to 200
     A1 = A(K,KM1)
     A2 = A(K1,KM1)
     A3 = A(K2,KM1)
  200    S = ABS(A1) + ABS(A2) + ABS(A3)
     if (S  ==  0.0E0) go to 260
     U1 = A1 / S
     U2 = A2 / S
     U3 = A3 / S
     R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
     V1 = -(U1 + R) / R
     V2 = -U2 / R
     V3 = -U3 / R
     U2 = V2 / V1
     U3 = V3 / V1
!
     DO 210 J = KM1, ENORN
        T = A(K,J) + U2 * A(K1,J) + U3 * A(K2,J)
        A(K,J) = A(K,J) + T * V1
        A(K1,J) = A(K1,J) + T * V2
        A(K2,J) = A(K2,J) + T * V3
        T = B(K,J) + U2 * B(K1,J) + U3 * B(K2,J)
        B(K,J) = B(K,J) + T * V1
        B(K1,J) = B(K1,J) + T * V2
        B(K2,J) = B(K2,J) + T * V3
  210    CONTINUE
!
     if (K  ==  L) go to 220
     A(K1,KM1) = 0.0E0
     A(K2,KM1) = 0.0E0
!     .......... ZERO B(K+2,K+1) AND B(K+2,K) ..........
  220    S = ABS(B(K2,K2)) + ABS(B(K2,K1)) + ABS(B(K2,K))
     if (S  ==  0.0E0) go to 240
     U1 = B(K2,K2) / S
     U2 = B(K2,K1) / S
     U3 = B(K2,K) / S
     R = SIGN(SQRT(U1*U1+U2*U2+U3*U3),U1)
     V1 = -(U1 + R) / R
     V2 = -U2 / R
     V3 = -U3 / R
     U2 = V2 / V1
     U3 = V3 / V1
!
     DO 230 I = LOR1, LL
        T = A(I,K2) + U2 * A(I,K1) + U3 * A(I,K)
        A(I,K2) = A(I,K2) + T * V1
        A(I,K1) = A(I,K1) + T * V2
        A(I,K) = A(I,K) + T * V3
        T = B(I,K2) + U2 * B(I,K1) + U3 * B(I,K)
        B(I,K2) = B(I,K2) + T * V1
        B(I,K1) = B(I,K1) + T * V2
        B(I,K) = B(I,K) + T * V3
  230    CONTINUE
!
     B(K2,K) = 0.0E0
     B(K2,K1) = 0.0E0
     if (.NOT. MATZ) go to 240
!
     DO 235 I = 1, N
        T = Z(I,K2) + U2 * Z(I,K1) + U3 * Z(I,K)
        Z(I,K2) = Z(I,K2) + T * V1
        Z(I,K1) = Z(I,K1) + T * V2
        Z(I,K) = Z(I,K) + T * V3
  235    CONTINUE
!     .......... ZERO B(K+1,K) ..........
  240    S = ABS(B(K1,K1)) + ABS(B(K1,K))
     if (S  ==  0.0E0) go to 260
     U1 = B(K1,K1) / S
     U2 = B(K1,K) / S
     R = SIGN(SQRT(U1*U1+U2*U2),U1)
     V1 = -(U1 + R) / R
     V2 = -U2 / R
     U2 = V2 / V1
!
     DO 250 I = LOR1, LL
        T = A(I,K1) + U2 * A(I,K)
        A(I,K1) = A(I,K1) + T * V1
        A(I,K) = A(I,K) + T * V2
        T = B(I,K1) + U2 * B(I,K)
        B(I,K1) = B(I,K1) + T * V1
        B(I,K) = B(I,K) + T * V2
  250    CONTINUE
!
     B(K1,K) = 0.0E0
     if (.NOT. MATZ) go to 260
!
     DO 255 I = 1, N
        T = Z(I,K1) + U2 * Z(I,K)
        Z(I,K1) = Z(I,K1) + T * V1
        Z(I,K) = Z(I,K) + T * V2
  255    CONTINUE
!
  260 CONTINUE
!     .......... END QZ STEP ..........
  go to 70
!     .......... SET ERROR -- NEITHER BOTTOM SUBDIAGONAL ELEMENT
!                HAS BECOME NEGLIGIBLE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
!     .......... SAVE EPSB FOR USE BY QZVAL AND QZVEC ..........
 1001 if (N  >  1) B(N,1) = EPSB
  return
end
