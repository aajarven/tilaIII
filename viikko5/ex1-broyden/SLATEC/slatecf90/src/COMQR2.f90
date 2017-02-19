subroutine COMQR2 (NM, N, LOW, IGH, ORTR, ORTI, HR, HI, WR, WI, &
     ZR, ZI, IERR)
!
!! COMQR2 computes the eigenvalues and eigenvectors of a complex upper ...
!  Hessenberg matrix.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C2B
!***TYPE      COMPLEX (HQR2-S, COMQR2-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of a unitary analogue of the
!     ALGOL procedure  COMLR2, NUM. MATH. 16, 181-204(1970) by Peters
!     and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!     The unitary analogue substitutes the QR algorithm of Francis
!     (COMP. JOUR. 4, 332-345(1962)) for the LR algorithm.
!
!     This subroutine finds the eigenvalues and eigenvectors
!     of a COMPLEX UPPER Hessenberg matrix by the QR
!     method.  The eigenvectors of a COMPLEX GENERAL matrix
!     can also be found if  CORTH  has been used to reduce
!     this general matrix to Hessenberg form.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, HR, HI, ZR, and ZI, as declared in the
!          calling program dimension statement.  NM is an INTEGER
!          variable.
!
!        N is the order of the matrix H=(HR,HI).  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  CBAL.  If  CBAL  has not been used,
!          set LOW=1 and IGH equal to the order of the matrix, N.
!
!        ORTR and ORTI contain information about the unitary trans-
!          formations used in the reduction by  CORTH, if performed.
!          Only elements LOW through IGH are used.  If the eigenvectors
!          of the Hessenberg matrix are desired, set ORTR(J) and
!          ORTI(J) to 0.0E0 for these elements.  ORTR and ORTI are
!          one-dimensional REAL arrays, dimensioned ORTR(IGH) and
!          ORTI(IGH).
!
!        HR and HI contain the real and imaginary parts, respectively,
!          of the complex upper Hessenberg matrix.  Their lower
!          triangles below the subdiagonal contain information about
!          the unitary transformations used in the reduction by  CORTH,
!          if performed.  If the eigenvectors of the Hessenberg matrix
!          are desired, these elements may be arbitrary.  HR and HI
!          are two-dimensional REAL arrays, dimensioned HR(NM,N) and
!          HI(NM,N).
!
!     On OUTPUT
!
!        ORTR, ORTI, and the upper Hessenberg portions of HR and HI
!          have been destroyed.
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues of the upper Hessenberg matrix.  If an
!          error exit is made, the eigenvalues should be correct for
!          indices IERR+1, IERR+2, ..., N.  WR and WI are one-
!          dimensional REAL arrays, dimensioned WR(N) and WI(N).
!
!        ZR and ZI contain the real and imaginary parts, respectively,
!          of the eigenvectors.  The eigenvectors are unnormalized.
!          If an error exit is made, none of the eigenvectors has been
!          found.  ZR and ZI are two-dimensional REAL arrays,
!          dimensioned ZR(NM,N) and ZI(NM,N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30*N iterations.
!                     The eigenvalues should be correct for indices
!                     IERR+1, IERR+2, ..., N, but no eigenvectors are
!                     computed.
!
!     Calls CSROOT for complex square root.
!     Calls PYTHAG(A,B) for sqrt(A**2 + B**2).
!     Calls CDIV for complex division.
!
!     Questions and comments should be directed to B. S. Garbow,
!     APPLIED MATHEMATICS DIVISION, ARGONNE NATIONAL LABORATORY
!     ------------------------------------------------------------------
!
!***REFERENCES  B. T. Smith, J. M. Boyle, J. J. Dongarra, B. S. Garbow,
!                 Y. Ikebe, V. C. Klema and C. B. Moler, Matrix Eigen-
!                 system Routines - EISPACK Guide, Springer-Verlag,
!                 1976.
!***ROUTINES CALLED  CDIV, CSROOT, PYTHAG
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  COMQR2
!
  INTEGER I,J,K,L,M,N,EN,II,JJ,LL,NM,NN,IGH,IP1
  INTEGER ITN,ITS,LOW,LP1,ENM1,IEND,IERR
  REAL HR(NM,*),HI(NM,*),WR(*),WI(*),ZR(NM,*),ZI(NM,*)
  REAL ORTR(*),ORTI(*)
  REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,S1,S2
  REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  COMQR2
  IERR = 0
!     .......... INITIALIZE EIGENVECTOR MATRIX ..........
  DO 100 I = 1, N
!
     DO 100 J = 1, N
        ZR(I,J) = 0.0E0
        ZI(I,J) = 0.0E0
        if (I  ==  J) ZR(I,J) = 1.0E0
  100 CONTINUE
!     .......... FORM THE MATRIX OF ACCUMULATED TRANSFORMATIONS
!                FROM THE INFORMATION LEFT BY CORTH ..........
  IEND = IGH - LOW - 1
  if (IEND) 180, 150, 105
!     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  105 DO 140 II = 1, IEND
     I = IGH - II
     if (ORTR(I)  ==  0.0E0 .AND. ORTI(I)  ==  0.0E0) go to 140
     if (HR(I,I-1)  ==  0.0E0 .AND. HI(I,I-1)  ==  0.0E0) go to 140
!     .......... NORM BELOW IS NEGATIVE OF H FORMED IN CORTH ..........
     NORM = HR(I,I-1) * ORTR(I) + HI(I,I-1) * ORTI(I)
     IP1 = I + 1
!
     DO 110 K = IP1, IGH
        ORTR(K) = HR(K,I-1)
        ORTI(K) = HI(K,I-1)
  110    CONTINUE
!
     DO 130 J = I, IGH
        SR = 0.0E0
        SI = 0.0E0
!
        DO 115 K = I, IGH
           SR = SR + ORTR(K) * ZR(K,J) + ORTI(K) * ZI(K,J)
           SI = SI + ORTR(K) * ZI(K,J) - ORTI(K) * ZR(K,J)
  115       CONTINUE
!
        SR = SR / NORM
        SI = SI / NORM
!
        DO 120 K = I, IGH
           ZR(K,J) = ZR(K,J) + SR * ORTR(K) - SI * ORTI(K)
           ZI(K,J) = ZI(K,J) + SR * ORTI(K) + SI * ORTR(K)
  120       CONTINUE
!
  130    CONTINUE
!
  140 CONTINUE
!     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  150 L = LOW + 1
!
  DO 170 I = L, IGH
     LL = MIN(I+1,IGH)
     if (HI(I,I-1)  ==  0.0E0) go to 170
     NORM = PYTHAG(HR(I,I-1),HI(I,I-1))
     YR = HR(I,I-1) / NORM
     YI = HI(I,I-1) / NORM
     HR(I,I-1) = NORM
     HI(I,I-1) = 0.0E0
!
     DO 155 J = I, N
        SI = YR * HI(I,J) - YI * HR(I,J)
        HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
        HI(I,J) = SI
  155    CONTINUE
!
     DO 160 J = 1, LL
        SI = YR * HI(J,I) + YI * HR(J,I)
        HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
        HI(J,I) = SI
  160    CONTINUE
!
     DO 165 J = LOW, IGH
        SI = YR * ZI(J,I) + YI * ZR(J,I)
        ZR(J,I) = YR * ZR(J,I) - YI * ZI(J,I)
        ZI(J,I) = SI
  165    CONTINUE
!
  170 CONTINUE
!     .......... STORE ROOTS ISOLATED BY CBAL ..........
  180 DO 200 I = 1, N
     if (I  >=  LOW .AND. I  <=  IGH) go to 200
     WR(I) = HR(I,I)
     WI(I) = HI(I,I)
  200 CONTINUE
!
  EN = IGH
  TR = 0.0E0
  TI = 0.0E0
  ITN = 30*N
!     .......... SEARCH FOR NEXT EIGENVALUE ..........
  220 if (EN  <  LOW) go to 680
  ITS = 0
  ENM1 = EN - 1
!     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!                FOR L=EN STEP -1 UNTIL LOW DO -- ..........
  240 DO 260 LL = LOW, EN
     L = EN + LOW - LL
     if (L  ==  LOW) go to 300
     S1 = ABS(HR(L-1,L-1)) + ABS(HI(L-1,L-1)) &
               + ABS(HR(L,L)) +ABS(HI(L,L))
     S2 = S1 + ABS(HR(L,L-1))
     if (S2  ==  S1) go to 300
  260 CONTINUE
!     .......... FORM SHIFT ..........
  300 if (L  ==  EN) go to 660
  if (ITN  ==  0) go to 1000
  if (ITS  ==  10 .OR. ITS  ==  20) go to 320
  SR = HR(EN,EN)
  SI = HI(EN,EN)
  XR = HR(ENM1,EN) * HR(EN,ENM1)
  XI = HI(ENM1,EN) * HR(EN,ENM1)
  if (XR  ==  0.0E0 .AND. XI  ==  0.0E0) go to 340
  YR = (HR(ENM1,ENM1) - SR) / 2.0E0
  YI = (HI(ENM1,ENM1) - SI) / 2.0E0
  call CSROOT(YR**2-YI**2+XR,2.0E0*YR*YI+XI,ZZR,ZZI)
  if (YR * ZZR + YI * ZZI  >=  0.0E0) go to 310
  ZZR = -ZZR
  ZZI = -ZZI
  310 call CDIV(XR,XI,YR+ZZR,YI+ZZI,XR,XI)
  SR = SR - XR
  SI = SI - XI
  go to 340
!     .......... FORM EXCEPTIONAL SHIFT ..........
  320 SR = ABS(HR(EN,ENM1)) + ABS(HR(ENM1,EN-2))
  SI = 0.0E0
!
  340 DO 360 I = LOW, EN
     HR(I,I) = HR(I,I) - SR
     HI(I,I) = HI(I,I) - SI
  360 CONTINUE
!
  TR = TR + SR
  TI = TI + SI
  ITS = ITS + 1
  ITN = ITN - 1
!     .......... REDUCE TO TRIANGLE (ROWS) ..........
  LP1 = L + 1
!
  DO 500 I = LP1, EN
     SR = HR(I,I-1)
     HR(I,I-1) = 0.0E0
     NORM = PYTHAG(PYTHAG(HR(I-1,I-1),HI(I-1,I-1)),SR)
     XR = HR(I-1,I-1) / NORM
     WR(I-1) = XR
     XI = HI(I-1,I-1) / NORM
     WI(I-1) = XI
     HR(I-1,I-1) = NORM
     HI(I-1,I-1) = 0.0E0
     HI(I,I-1) = SR / NORM
!
     DO 490 J = I, N
        YR = HR(I-1,J)
        YI = HI(I-1,J)
        ZZR = HR(I,J)
        ZZI = HI(I,J)
        HR(I-1,J) = XR * YR + XI * YI + HI(I,I-1) * ZZR
        HI(I-1,J) = XR * YI - XI * YR + HI(I,I-1) * ZZI
        HR(I,J) = XR * ZZR - XI * ZZI - HI(I,I-1) * YR
        HI(I,J) = XR * ZZI + XI * ZZR - HI(I,I-1) * YI
  490    CONTINUE
!
  500 CONTINUE
!
  SI = HI(EN,EN)
  if (SI  ==  0.0E0) go to 540
  NORM = PYTHAG(HR(EN,EN),SI)
  SR = HR(EN,EN) / NORM
  SI = SI / NORM
  HR(EN,EN) = NORM
  HI(EN,EN) = 0.0E0
  if (EN  ==  N) go to 540
  IP1 = EN + 1
!
  DO 520 J = IP1, N
     YR = HR(EN,J)
     YI = HI(EN,J)
     HR(EN,J) = SR * YR + SI * YI
     HI(EN,J) = SR * YI - SI * YR
  520 CONTINUE
!     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
     XR = WR(J-1)
     XI = WI(J-1)
!
     DO 580 I = 1, J
        YR = HR(I,J-1)
        YI = 0.0E0
        ZZR = HR(I,J)
        ZZI = HI(I,J)
        if (I  ==  J) go to 560
        YI = HI(I,J-1)
        HI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
  560       HR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
        HR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
        HI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  580    CONTINUE
!
     DO 590 I = LOW, IGH
        YR = ZR(I,J-1)
        YI = ZI(I,J-1)
        ZZR = ZR(I,J)
        ZZI = ZI(I,J)
        ZR(I,J-1) = XR * YR - XI * YI + HI(J,J-1) * ZZR
        ZI(I,J-1) = XR * YI + XI * YR + HI(J,J-1) * ZZI
        ZR(I,J) = XR * ZZR + XI * ZZI - HI(J,J-1) * YR
        ZI(I,J) = XR * ZZI - XI * ZZR - HI(J,J-1) * YI
  590    CONTINUE
!
  600 CONTINUE
!
  if (SI  ==  0.0E0) go to 240
!
  DO 630 I = 1, EN
     YR = HR(I,EN)
     YI = HI(I,EN)
     HR(I,EN) = SR * YR - SI * YI
     HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
!
  DO 640 I = LOW, IGH
     YR = ZR(I,EN)
     YI = ZI(I,EN)
     ZR(I,EN) = SR * YR - SI * YI
     ZI(I,EN) = SR * YI + SI * YR
  640 CONTINUE
!
  go to 240
!     .......... A ROOT FOUND ..........
  660 HR(EN,EN) = HR(EN,EN) + TR
  WR(EN) = HR(EN,EN)
  HI(EN,EN) = HI(EN,EN) + TI
  WI(EN) = HI(EN,EN)
  EN = ENM1
  go to 220
!     .......... ALL ROOTS FOUND.  BACKSUBSTITUTE TO FIND
!                VECTORS OF UPPER TRIANGULAR FORM ..........
  680 NORM = 0.0E0
!
  DO 720 I = 1, N
!
     DO 720 J = I, N
        NORM = NORM + ABS(HR(I,J)) + ABS(HI(I,J))
  720 CONTINUE
!
  if (N  ==  1 .OR. NORM  ==  0.0E0) go to 1001
!     .......... FOR EN=N STEP -1 UNTIL 2 DO -- ..........
  DO 800 NN = 2, N
     EN = N + 2 - NN
     XR = WR(EN)
     XI = WI(EN)
     ENM1 = EN - 1
!     .......... FOR I=EN-1 STEP -1 UNTIL 1 DO -- ..........
     DO 780 II = 1, ENM1
        I = EN - II
        ZZR = HR(I,EN)
        ZZI = HI(I,EN)
        if (I  ==  ENM1) go to 760
        IP1 = I + 1
!
        DO 740 J = IP1, ENM1
           ZZR = ZZR + HR(I,J) * HR(J,EN) - HI(I,J) * HI(J,EN)
           ZZI = ZZI + HR(I,J) * HI(J,EN) + HI(I,J) * HR(J,EN)
  740       CONTINUE
!
  760       YR = XR - WR(I)
        YI = XI - WI(I)
        if (YR  /=  0.0E0 .OR. YI  /=  0.0E0) go to 775
        YR = NORM
  770       YR = 0.5E0*YR
        if (NORM + YR  >  NORM) go to 770
        YR = 2.0E0*YR
  775       call CDIV(ZZR,ZZI,YR,YI,HR(I,EN),HI(I,EN))
  780    CONTINUE
!
  800 CONTINUE
!     .......... END BACKSUBSTITUTION ..........
  ENM1 = N - 1
!     .......... VECTORS OF ISOLATED ROOTS ..........
  DO  840 I = 1, ENM1
     if (I  >=  LOW .AND. I  <=  IGH) go to 840
     IP1 = I + 1
!
     DO 820 J = IP1, N
        ZR(I,J) = HR(I,J)
        ZI(I,J) = HI(I,J)
  820    CONTINUE
!
  840 CONTINUE
!     .......... MULTIPLY BY TRANSFORMATION MATRIX TO GIVE
!                VECTORS OF ORIGINAL FULL MATRIX.
!                FOR J=N STEP -1 UNTIL LOW+1 DO -- ..........
  DO 880 JJ = LOW, ENM1
     J = N + LOW - JJ
     M = MIN(J-1,IGH)
!
     DO 880 I = LOW, IGH
        ZZR = ZR(I,J)
        ZZI = ZI(I,J)
!
        DO 860 K = LOW, M
           ZZR = ZZR + ZR(I,K) * HR(K,J) - ZI(I,K) * HI(K,J)
           ZZI = ZZI + ZR(I,K) * HI(K,J) + ZI(I,K) * HR(K,J)
  860       CONTINUE
!
        ZR(I,J) = ZZR
        ZI(I,J) = ZZI
  880 CONTINUE
!
  go to 1001
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
end
