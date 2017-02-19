subroutine TEVLS (N, D, E2, IERR)
!
!! TEVLS finds eigenvalues of a symmetric tridiagonal matrix by rational QL.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to BLKTRI
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (TEVLS-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine finds the eigenvalues of a symmetric tridiagonal
!     matrix by the rational QL method.
!
!     On Input-
!
!        N is the order of the matrix,
!
!        D contains the diagonal elements of the input matrix,
!
!        E2 contains the subdiagonal elements of the input matrix
!           in its last N-1 positions.  E2(1) is arbitrary.
!
!      On Output-
!
!        D contains the eigenvalues in ascending order.  If an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...IERR-1, but may not be
!          the smallest eigenvalues,
!
!        E2 has been destroyed,
!
!        IERR is set to
!          ZERO       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after 30 iterations.
!
!***SEE ALSO  BLKTRI
!***REFERENCES  C. H. Reinsch, Eigenvalues of a real, symmetric, tri-
!                 diagonal matrix, Algorithm 464, Communications of the
!                 ACM 16, 11 (November 1973), pp. 689.
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    CBLKT
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!   920528  DESCRIPTION revised and REFERENCES section added.  (WRB)
!***END PROLOGUE  TEVLS
!
  INTEGER         I          ,J          ,L          ,M          , &
                  N          ,II         ,L1         ,MML        , &
                  IERR
  REAL            D(*)       ,E2(*)
  REAL            B          ,C          ,F          ,G          , &
                  H          ,P          ,R          ,S          , &
                  MACHEP
!
  COMMON /CBLKT/  NPP        ,K          ,MACHEP     ,CNV        , &
                  NM         ,NCMPLX     ,IK
!***FIRST EXECUTABLE STATEMENT  TEVLS
  IERR = 0
  if (N  ==  1) go to 115
!
  DO 101 I=2,N
     E2(I-1) = E2(I)*E2(I)
  101 CONTINUE
!
  F = 0.0
  B = 0.0
  E2(N) = 0.0
!
  DO 112 L=1,N
     J = 0
     H = MACHEP*(ABS(D(L))+SQRT(E2(L)))
     if (B  >  H) go to 102
     B = H
     C = B*B
!
!     ********** LOOK FOR SMALL SQUARED SUB-DIAGONAL ELEMENT **********
!
  102    DO 103 M=L,N
        if (E2(M)  <=  C) go to 104
!
!     ********** E2(N) IS ALWAYS ZERO, SO THERE IS NO EXIT
!                THROUGH THE BOTTOM OF THE LOOP **********
!
  103    CONTINUE
!
  104    if (M  ==  L) go to 108
  105    if (J  ==  30) go to 114
     J = J+1
!
!     ********** FORM SHIFT **********
!
     L1 = L+1
     S = SQRT(E2(L))
     G = D(L)
     P = (D(L1)-G)/(2.0*S)
     R = SQRT(P*P+1.0)
     D(L) = S/(P+SIGN(R,P))
     H = G-D(L)
!
     DO 106 I=L1,N
        D(I) = D(I)-H
  106    CONTINUE
!
     F = F+H
!
!     ********** RATIONAL QL TRANSFORMATION **********
!
     G = D(M)
     if (G  ==  0.0) G = B
     H = G
     S = 0.0
     MML = M-L
!
!     ********** FOR I=M-1 STEP -1 UNTIL L DO -- **********
!
     DO 107 II=1,MML
        I = M-II
        P = G*H
        R = P+E2(I)
        E2(I+1) = S*R
        S = E2(I)/R
        D(I+1) = H+S*(H+D(I))
        G = D(I)-E2(I)/G
        if (G  ==  0.0) G = B
        H = G*P/R
  107    CONTINUE
!
     E2(L) = S*G
     D(L) = H
!
!     ********** GUARD AGAINST UNDERFLOWED H **********
!
     if (H  ==  0.0) go to 108
     if (ABS(E2(L))  <=  ABS(C/H)) go to 108
     E2(L) = H*E2(L)
     if (E2(L)  /=  0.0) go to 105
  108    P = D(L)+F
!
!     ********** ORDER EIGENVALUES **********
!
     if (L  ==  1) go to 110
!
!     ********** FOR I=L STEP -1 UNTIL 2 DO -- **********
!
     DO 109 II=2,L
        I = L+2-II
        if (P  >=  D(I-1)) go to 111
        D(I) = D(I-1)
  109    CONTINUE
!
  110    I = 1
  111    D(I) = P
  112 CONTINUE
!
  if (ABS(D(N))  >=  ABS(D(1))) go to 115
  NHALF = N/2
  DO 113 I=1,NHALF
     NTOP = N-I
     DHOLD = D(I)
     D(I) = D(NTOP+1)
     D(NTOP+1) = DHOLD
  113 CONTINUE
  go to 115
!
!     ********** SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30 ITERATIONS **********
!
  114 IERR = L
  115 RETURN
!
!     ********** LAST CARD OF TQLRAT **********
!
end
