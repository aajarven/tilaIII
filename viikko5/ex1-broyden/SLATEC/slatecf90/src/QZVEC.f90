subroutine QZVEC (NM, N, A, B, ALFR, ALFI, BETA, Z)
!
!! QZVEC is the fourth step of the QZ algorithm for generalized eigenproblems.
!  Accepts a matrix in
!            quasi-triangular form and another in upper triangular
!            and computes the eigenvectors of the triangular problem
!            and transforms them back to the original coordinates
!            Usually preceded by QZHES, QZIT, and QZVAL.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C3
!***TYPE      SINGLE PRECISION (QZVEC-S)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is the optional fourth step of the QZ algorithm
!     for solving generalized matrix eigenvalue problems,
!     SIAM J. NUMER. ANAL. 10, 241-256(1973) by MOLER and STEWART.
!
!     This subroutine accepts a pair of REAL matrices, one of them in
!     quasi-triangular form (in which each 2-by-2 block corresponds to
!     a pair of complex eigenvalues) and the other in upper triangular
!     form.  It computes the eigenvectors of the triangular problem and
!     transforms the results back to the original coordinate system.
!     It is usually preceded by  QZHES,  QZIT, and  QZVAL.
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
!        A contains a real upper quasi-triangular matrix.  A is a two-
!          dimensional REAL array, dimensioned A(NM,N).
!
!        B contains a real upper triangular matrix.  In addition,
!          location B(N,1) contains the tolerance quantity (EPSB)
!          computed and saved in  QZIT.  B is a two-dimensional REAL
!          array, dimensioned B(NM,N).
!
!        ALFR, ALFI, and BETA are one-dimensional REAL arrays with
!          components whose ratios ((ALFR+I*ALFI)/BETA) are the
!          generalized eigenvalues.  They are usually obtained from
!          QZVAL.  They are dimensioned ALFR(N), ALFI(N), and BETA(N).
!
!        Z contains the transformation matrix produced in the reductions
!          by  QZHES,  QZIT, and  QZVAL,  if performed.  If the
!          eigenvectors of the triangular problem are desired, Z must
!          contain the identity matrix.  Z is a two-dimensional REAL
!          array, dimensioned Z(NM,N).
!
!     On Output
!
!        A is unaltered.  Its subdiagonal elements provide information
!           about the storage of the complex eigenvectors.
!
!        B has been destroyed.
!
!        ALFR, ALFI, and BETA are unaltered.
!
!        Z contains the real and imaginary parts of the eigenvectors.
!          If ALFI(J)  ==  0.0, the J-th eigenvalue is real and
!            the J-th column of Z contains its eigenvector.
!          If ALFI(J)  /=  0.0, the J-th eigenvalue is complex.
!            If ALFI(J)  >  0.0, the eigenvalue is the first of
!              a complex pair and the J-th and (J+1)-th columns
!              of Z contain its eigenvector.
!            If ALFI(J)  <  0.0, the eigenvalue is the second of
!              a complex pair and the (J-1)-th and J-th columns
!              of Z contain the conjugate of its eigenvector.
!          Each eigenvector is normalized so that the modulus
!          of its largest component is 1.0 .
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
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  QZVEC
!
  INTEGER I,J,K,M,N,EN,II,JJ,NA,NM,NN,ISW,ENM2
  REAL A(NM,*),B(NM,*),ALFR(*),ALFI(*),BETA(*),Z(NM,*)
  REAL D,Q,R,S,T,W,X,Y,DI,DR,RA,RR,SA,TI,TR,T1,T2
  REAL W1,X1,ZZ,Z1,ALFM,ALMI,ALMR,BETM,EPSB
!
!***FIRST EXECUTABLE STATEMENT  QZVEC
  EPSB = B(N,1)
  ISW = 1
!     .......... FOR EN=N STEP -1 UNTIL 1 DO -- ..........
  DO 800 NN = 1, N
     EN = N + 1 - NN
     NA = EN - 1
     if (ISW  ==  2) go to 795
     if (ALFI(EN)  /=  0.0E0) go to 710
!     .......... REAL VECTOR ..........
     M = EN
     B(EN,EN) = 1.0E0
     if (NA  ==  0) go to 800
     ALFM = ALFR(M)
     BETM = BETA(M)
!     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
     DO 700 II = 1, NA
        I = EN - II
        W = BETM * A(I,I) - ALFM * B(I,I)
        R = 0.0E0
!
        DO 610 J = M, EN
  610       R = R + (BETM * A(I,J) - ALFM * B(I,J)) * B(J,EN)
!
        if (I  ==  1 .OR. ISW  ==  2) go to 630
        if (BETM * A(I,I-1)  ==  0.0E0) go to 630
        ZZ = W
        S = R
        go to 690
  630       M = I
        if (ISW  ==  2) go to 640
!     .......... REAL 1-BY-1 BLOCK ..........
        T = W
        if (W  ==  0.0E0) T = EPSB
        B(I,EN) = -R / T
        go to 700
!     .......... REAL 2-BY-2 BLOCK ..........
  640       X = BETM * A(I,I+1) - ALFM * B(I,I+1)
        Y = BETM * A(I+1,I)
        Q = W * ZZ - X * Y
        T = (X * S - ZZ * R) / Q
        B(I,EN) = T
        if (ABS(X)  <=  ABS(ZZ)) go to 650
        B(I+1,EN) = (-R - W * T) / X
        go to 690
  650       B(I+1,EN) = (-S - Y * T) / ZZ
  690       ISW = 3 - ISW
  700    CONTINUE
!     .......... END REAL VECTOR ..........
     go to 800
!     .......... COMPLEX VECTOR ..........
  710    M = NA
     ALMR = ALFR(M)
     ALMI = ALFI(M)
     BETM = BETA(M)
!     .......... LAST VECTOR COMPONENT CHOSEN IMAGINARY SO THAT
!                EIGENVECTOR MATRIX IS TRIANGULAR ..........
     Y = BETM * A(EN,NA)
     B(NA,NA) = -ALMI * B(EN,EN) / Y
     B(NA,EN) = (ALMR * B(EN,EN) - BETM * A(EN,EN)) / Y
     B(EN,NA) = 0.0E0
     B(EN,EN) = 1.0E0
     ENM2 = NA - 1
     if (ENM2  ==  0) go to 795
!     .......... FOR I=EN-2 STEP -1 UNTIL 1 DO -- ..........
     DO 790 II = 1, ENM2
        I = NA - II
        W = BETM * A(I,I) - ALMR * B(I,I)
        W1 = -ALMI * B(I,I)
        RA = 0.0E0
        SA = 0.0E0
!
        DO 760 J = M, EN
           X = BETM * A(I,J) - ALMR * B(I,J)
           X1 = -ALMI * B(I,J)
           RA = RA + X * B(J,NA) - X1 * B(J,EN)
           SA = SA + X * B(J,EN) + X1 * B(J,NA)
  760       CONTINUE
!
        if (I  ==  1 .OR. ISW  ==  2) go to 770
        if (BETM * A(I,I-1)  ==  0.0E0) go to 770
        ZZ = W
        Z1 = W1
        R = RA
        S = SA
        ISW = 2
        go to 790
  770       M = I
        if (ISW  ==  2) go to 780
!     .......... COMPLEX 1-BY-1 BLOCK ..........
        TR = -RA
        TI = -SA
  773       DR = W
        DI = W1
!     .......... COMPLEX DIVIDE (T1,T2) = (TR,TI) / (DR,DI) ..........
  775       if (ABS(DI)  >  ABS(DR)) go to 777
        RR = DI / DR
        D = DR + DI * RR
        T1 = (TR + TI * RR) / D
        T2 = (TI - TR * RR) / D
        go to (787,782), ISW
  777       RR = DR / DI
        D = DR * RR + DI
        T1 = (TR * RR + TI) / D
        T2 = (TI * RR - TR) / D
        go to (787,782), ISW
!     .......... COMPLEX 2-BY-2 BLOCK ..........
  780       X = BETM * A(I,I+1) - ALMR * B(I,I+1)
        X1 = -ALMI * B(I,I+1)
        Y = BETM * A(I+1,I)
        TR = Y * RA - W * R + W1 * S
        TI = Y * SA - W * S - W1 * R
        DR = W * ZZ - W1 * Z1 - X * Y
        DI = W * Z1 + W1 * ZZ - X1 * Y
        if (DR  ==  0.0E0 .AND. DI  ==  0.0E0) DR = EPSB
        go to 775
  782       B(I+1,NA) = T1
        B(I+1,EN) = T2
        ISW = 1
        if (ABS(Y)  >  ABS(W) + ABS(W1)) go to 785
        TR = -RA - X * B(I+1,NA) + X1 * B(I+1,EN)
        TI = -SA - X * B(I+1,EN) - X1 * B(I+1,NA)
        go to 773
  785       T1 = (-R - ZZ * B(I+1,NA) + Z1 * B(I+1,EN)) / Y
        T2 = (-S - ZZ * B(I+1,EN) - Z1 * B(I+1,NA)) / Y
  787       B(I,NA) = T1
        B(I,EN) = T2
  790    CONTINUE
!     .......... END COMPLEX VECTOR ..........
  795    ISW = 3 - ISW
  800 CONTINUE
!     .......... END BACK SUBSTITUTION.
!                TRANSFORM TO ORIGINAL COORDINATE SYSTEM.
!                FOR J=N STEP -1 UNTIL 1 DO -- ..........
  DO 880 JJ = 1, N
     J = N + 1 - JJ
!
     DO 880 I = 1, N
        ZZ = 0.0E0
!
        DO 860 K = 1, J
  860       ZZ = ZZ + Z(I,K) * B(K,J)
!
        Z(I,J) = ZZ
  880 CONTINUE
!     .......... NORMALIZE SO THAT MODULUS OF LARGEST
!                COMPONENT OF EACH VECTOR IS 1.
!                (ISW IS 1 INITIALLY FROM BEFORE) ..........
  DO 950 J = 1, N
     D = 0.0E0
     if (ISW  ==  2) go to 920
     if (ALFI(J)  /=  0.0E0) go to 945
!
     DO 890 I = 1, N
        if (ABS(Z(I,J))  >  D) D = ABS(Z(I,J))
  890    CONTINUE
!
     DO 900 I = 1, N
  900    Z(I,J) = Z(I,J) / D
!
     go to 950
!
  920    DO 930 I = 1, N
        R = ABS(Z(I,J-1)) + ABS(Z(I,J))
        if (R  /=  0.0E0) R = R * SQRT((Z(I,J-1)/R)**2 &
                                       +(Z(I,J)/R)**2)
        if (R  >  D) D = R
  930    CONTINUE
!
     DO 940 I = 1, N
        Z(I,J-1) = Z(I,J-1) / D
        Z(I,J) = Z(I,J) / D
  940    CONTINUE
!
  945    ISW = 3 - ISW
  950 CONTINUE
!
  return
end
