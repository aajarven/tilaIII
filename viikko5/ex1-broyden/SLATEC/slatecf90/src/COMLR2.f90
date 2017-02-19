subroutine COMLR2 (NM, N, LOW, IGH, INT, HR, HI, WR, WI, ZR, ZI, &
     IERR)
!
!! COMLR2 domputes the eigenvalues and eigenvectors of a complex upper ...
!  Hessenberg matrix using the modified LR method.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C2B
!***TYPE      COMPLEX (COMLR2-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK, LR METHOD
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure COMLR2,
!     NUM. MATH. 16, 181-204(1970) by Peters and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 372-395(1971).
!
!     This subroutine finds the eigenvalues and eigenvectors
!     of a COMPLEX UPPER Hessenberg matrix by the modified LR
!     method.  The eigenvectors of a COMPLEX GENERAL matrix
!     can also be found if  COMHES  has been used to reduce
!     this general matrix to Hessenberg form.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, HR, HI, ZR and ZI, as declared in the
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
!        INT contains information on the rows and columns
!          interchanged in the reduction by  COMHES, if performed.
!          Only elements LOW through IGH are used.  If you want the
!          eigenvectors of a complex general matrix, leave INT as it
!          came from  COMHES.  If the eigenvectors of the Hessenberg
!          matrix are desired, set INT(J)=J for these elements.  INT
!          is a one-dimensional INTEGER array, dimensioned INT(IGH).
!
!        HR and HI contain the real and imaginary parts, respectively,
!          of the complex upper Hessenberg matrix.  Their lower
!          triangles below the subdiagonal contain the multipliers
!          which were used in the reduction by  COMHES, if performed.
!          If the eigenvectors of a complex general matrix are
!          desired, leave these multipliers in the lower triangles.
!          If the eigenvectors of the Hessenberg matrix are desired,
!          these elements must be set to zero.  HR and HI are
!          two-dimensional REAL arrays, dimensioned HR(NM,N) and
!          HI(NM,N).
!
!     On OUTPUT
!
!        The upper Hessenberg portions of HR and HI have been
!          destroyed, but the location HR(1,1) contains the norm
!          of the triangularized matrix.
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
!***ROUTINES CALLED  CDIV, CSROOT
!***REVISION HISTORY  (YYMMDD)
!   760101  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  COMLR2
!
  INTEGER I,J,K,L,M,N,EN,II,JJ,LL,MM,NM,NN,IGH,IM1,IP1
  INTEGER ITN,ITS,LOW,MP1,ENM1,IEND,IERR
  REAL HR(NM,*),HI(NM,*),WR(*),WI(*),ZR(NM,*),ZI(NM,*)
  REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,S1,S2
  INTEGER INT(*)
!
!***FIRST EXECUTABLE STATEMENT  COMLR2
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
!                FROM THE INFORMATION LEFT BY COMHES ..........
  IEND = IGH - LOW - 1
  if (IEND  <=  0) go to 180
!     .......... FOR I=IGH-1 STEP -1 UNTIL LOW+1 DO -- ..........
  DO 160 II = 1, IEND
     I = IGH - II
     IP1 = I + 1
!
     DO 120 K = IP1, IGH
        ZR(K,I) = HR(K,I-1)
        ZI(K,I) = HI(K,I-1)
  120    CONTINUE
!
     J = INT(I)
     if (I  ==  J) go to 160
!
     DO 140 K = I, IGH
        ZR(I,K) = ZR(J,K)
        ZI(I,K) = ZI(J,K)
        ZR(J,K) = 0.0E0
        ZI(J,K) = 0.0E0
  140    CONTINUE
!
     ZR(J,I) = 1.0E0
  160 CONTINUE
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
               + ABS(HR(L,L)) + ABS(HI(L,L))
     S2 = S1 + ABS(HR(L,L-1)) + ABS(HI(L,L-1))
     if (S2  ==  S1) go to 300
  260 CONTINUE
!     .......... FORM SHIFT ..........
  300 if (L  ==  EN) go to 660
  if (ITN  ==  0) go to 1000
  if (ITS  ==  10 .OR. ITS  ==  20) go to 320
  SR = HR(EN,EN)
  SI = HI(EN,EN)
  XR = HR(ENM1,EN) * HR(EN,ENM1) - HI(ENM1,EN) * HI(EN,ENM1)
  XI = HR(ENM1,EN) * HI(EN,ENM1) + HI(ENM1,EN) * HR(EN,ENM1)
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
  SI = ABS(HI(EN,ENM1)) + ABS(HI(ENM1,EN-2))
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
!     .......... LOOK FOR TWO CONSECUTIVE SMALL
!                SUB-DIAGONAL ELEMENTS ..........
  XR = ABS(HR(ENM1,ENM1)) + ABS(HI(ENM1,ENM1))
  YR = ABS(HR(EN,ENM1)) + ABS(HI(EN,ENM1))
  ZZR = ABS(HR(EN,EN)) + ABS(HI(EN,EN))
!     .......... FOR M=EN-1 STEP -1 UNTIL L DO -- ..........
  DO 380 MM = L, ENM1
     M = ENM1 + L - MM
     if (M  ==  L) go to 420
     YI = YR
     YR = ABS(HR(M,M-1)) + ABS(HI(M,M-1))
     XI = ZZR
     ZZR = XR
     XR = ABS(HR(M-1,M-1)) + ABS(HI(M-1,M-1))
     S1 = ZZR / YI * (ZZR + XR + XI)
     S2 = S1 + YR
     if (S2  ==  S1) go to 420
  380 CONTINUE
!     .......... TRIANGULAR DECOMPOSITION H=L*R ..........
  420 MP1 = M + 1
!
  DO 520 I = MP1, EN
     IM1 = I - 1
     XR = HR(IM1,IM1)
     XI = HI(IM1,IM1)
     YR = HR(I,IM1)
     YI = HI(I,IM1)
     if (ABS(XR) + ABS(XI)  >=  ABS(YR) + ABS(YI)) go to 460
!     .......... INTERCHANGE ROWS OF HR AND HI ..........
     DO 440 J = IM1, N
        ZZR = HR(IM1,J)
        HR(IM1,J) = HR(I,J)
        HR(I,J) = ZZR
        ZZI = HI(IM1,J)
        HI(IM1,J) = HI(I,J)
        HI(I,J) = ZZI
  440    CONTINUE
!
     call CDIV(XR,XI,YR,YI,ZZR,ZZI)
     WR(I) = 1.0E0
     go to 480
  460    call CDIV(YR,YI,XR,XI,ZZR,ZZI)
     WR(I) = -1.0E0
  480    HR(I,IM1) = ZZR
     HI(I,IM1) = ZZI
!
     DO 500 J = I, N
        HR(I,J) = HR(I,J) - ZZR * HR(IM1,J) + ZZI * HI(IM1,J)
        HI(I,J) = HI(I,J) - ZZR * HI(IM1,J) - ZZI * HR(IM1,J)
  500    CONTINUE
!
  520 CONTINUE
!     .......... COMPOSITION R*L=H ..........
  DO 640 J = MP1, EN
     XR = HR(J,J-1)
     XI = HI(J,J-1)
     HR(J,J-1) = 0.0E0
     HI(J,J-1) = 0.0E0
!     .......... INTERCHANGE COLUMNS OF HR, HI, ZR, AND ZI,
!                if NECESSARY ..........
     if (WR(J)  <=  0.0E0) go to 580
!
     DO 540 I = 1, J
        ZZR = HR(I,J-1)
        HR(I,J-1) = HR(I,J)
        HR(I,J) = ZZR
        ZZI = HI(I,J-1)
        HI(I,J-1) = HI(I,J)
        HI(I,J) = ZZI
  540    CONTINUE
!
     DO 560 I = LOW, IGH
        ZZR = ZR(I,J-1)
        ZR(I,J-1) = ZR(I,J)
        ZR(I,J) = ZZR
        ZZI = ZI(I,J-1)
        ZI(I,J-1) = ZI(I,J)
        ZI(I,J) = ZZI
  560    CONTINUE
!
  580    DO 600 I = 1, J
        HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
        HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
!     .......... ACCUMULATE TRANSFORMATIONS ..........
     DO 620 I = LOW, IGH
        ZR(I,J-1) = ZR(I,J-1) + XR * ZR(I,J) - XI * ZI(I,J)
        ZI(I,J-1) = ZI(I,J-1) + XR * ZI(I,J) + XI * ZR(I,J)
  620    CONTINUE
!
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
  HR(1,1) = NORM
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
  DO 840 I = 1, ENM1
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
