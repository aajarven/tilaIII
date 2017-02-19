subroutine COMQR (NM, N, LOW, IGH, HR, HI, WR, WI, IERR)
!
!! COMQR computes the eigenvalues of complex upper Hessenberg matrix ...
!  using the QR method.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C2B
!***TYPE      COMPLEX (HQR-S, COMQR-C)
!***KEYWORDS  EIGENVALUES, EIGENVECTORS, EISPACK
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of a unitary analogue of the
!     ALGOL procedure  COMLR, NUM. MATH. 12, 369-376(1968) by Martin
!     and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
!     The unitary analogue substitutes the QR algorithm of Francis
!     (COMP. JOUR. 4, 332-345(1962)) for the LR algorithm.
!
!     This subroutine finds the eigenvalues of a COMPLEX
!     upper Hessenberg matrix by the QR method.
!
!     On INPUT
!
!        NM must be set to the row dimension of the two-dimensional
!          array parameters, HR and HI, as declared in the calling
!          program dimension statement.  NM is an INTEGER variable.
!
!        N is the order of the matrix H=(HR,HI).  N is an INTEGER
!          variable.  N must be less than or equal to NM.
!
!        LOW and IGH are two INTEGER variables determined by the
!          balancing subroutine  CBAL.  If  CBAL  has not been used,
!          set LOW=1 and IGH equal to the order of the matrix, N.
!
!        HR and HI contain the real and imaginary parts, respectively,
!          of the complex upper Hessenberg matrix.  Their lower
!          triangles below the subdiagonal contain information about
!          the unitary transformations used in the reduction by  CORTH,
!          if performed.  HR and HI are two-dimensional REAL arrays,
!          dimensioned HR(NM,N) and HI(NM,N).
!
!     On OUTPUT
!
!        The upper Hessenberg portions of HR and HI have been
!          destroyed.  Therefore, they must be saved before calling
!          COMQR  if subsequent calculation of eigenvectors is to
!          be performed.
!
!        WR and WI contain the real and imaginary parts, respectively,
!          of the eigenvalues of the upper Hessenberg matrix.  If an
!          error exit is made, the eigenvalues should be correct for
!          indices IERR+1, IERR+2, ..., N.  WR and WI are one-
!          dimensional REAL arrays, dimensioned WR(N) and WI(N).
!
!        IERR is an INTEGER flag set to
!          Zero       for normal return,
!          J          if the J-th eigenvalue has not been
!                     determined after a total of 30*N iterations.
!                     The eigenvalues should be correct for indices
!                     IERR+1, IERR+2, ..., N.
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
!***END PROLOGUE  COMQR
!
  INTEGER I,J,L,N,EN,LL,NM,IGH,ITN,ITS,LOW,LP1,ENM1,IERR
  REAL HR(NM,*),HI(NM,*),WR(*),WI(*)
  REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,NORM,S1,S2
  REAL PYTHAG
!
!***FIRST EXECUTABLE STATEMENT  COMQR
  IERR = 0
  if (LOW  ==  IGH) go to 180
!     .......... CREATE REAL SUBDIAGONAL ELEMENTS ..........
  L = LOW + 1
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
     DO 155 J = I, IGH
        SI = YR * HI(I,J) - YI * HR(I,J)
        HR(I,J) = YR * HR(I,J) + YI * HI(I,J)
        HI(I,J) = SI
  155    CONTINUE
!
     DO 160 J = LOW, LL
        SI = YR * HI(J,I) + YI * HR(J,I)
        HR(J,I) = YR * HR(J,I) - YI * HI(J,I)
        HI(J,I) = SI
  160    CONTINUE
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
  220 if (EN  <  LOW) go to 1001
  ITS = 0
  ENM1 = EN - 1
!     .......... LOOK FOR SINGLE SMALL SUB-DIAGONAL ELEMENT
!                FOR L=EN STEP -1 UNTIL LOWE0 -- ..........
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
     DO 490 J = I, EN
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
!     .......... INVERSE OPERATION (COLUMNS) ..........
  540 DO 600 J = LP1, EN
     XR = WR(J-1)
     XI = WI(J-1)
!
     DO 580 I = L, J
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
  600 CONTINUE
!
  if (SI  ==  0.0E0) go to 240
!
  DO 630 I = L, EN
     YR = HR(I,EN)
     YI = HI(I,EN)
     HR(I,EN) = SR * YR - SI * YI
     HI(I,EN) = SR * YI + SI * YR
  630 CONTINUE
!
  go to 240
!     .......... A ROOT FOUND ..........
  660 WR(EN) = HR(EN,EN) + TR
  WI(EN) = HI(EN,EN) + TI
  EN = ENM1
  go to 220
!     .......... SET ERROR -- NO CONVERGENCE TO AN
!                EIGENVALUE AFTER 30*N ITERATIONS ..........
 1000 IERR = EN
 1001 RETURN
end
