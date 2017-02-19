subroutine COMLR (NM, N, LOW, IGH, HR, HI, WR, WI, IERR)
!
!! COMLR computes the eigenvalues of a complex upper Hessenberg ...
!  matrix using the modified LR method.
!
!***LIBRARY   SLATEC (EISPACK)
!***CATEGORY  D4C2B
!***TYPE      COMPLEX (COMLR-C)
!***KEYWORDS  EIGENVALUES, EISPACK, LR METHOD
!***AUTHOR  Smith, B. T., et al.
!***DESCRIPTION
!
!     This subroutine is a translation of the ALGOL procedure COMLR,
!     NUM. MATH. 12, 369-376(1968) by Martin and Wilkinson.
!     HANDBOOK FOR AUTO. COMP., VOL.II-LINEAR ALGEBRA, 396-403(1971).
!
!     This subroutine finds the eigenvalues of a COMPLEX
!     UPPER Hessenberg matrix by the modified LR method.
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
!          triangles below the subdiagonal contain the multipliers
!          which were used in the reduction by  COMHES, if performed.
!          HR and HI are two-dimensional REAL arrays, dimensioned
!          HR(NM,N) and HI(NM,N).
!
!     On OUTPUT
!
!        The upper Hessenberg portions of HR and HI have been
!          destroyed.  Therefore, they must be saved before calling
!          COMLR  if subsequent calculation of eigenvectors is to
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
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  COMLR
!
  INTEGER I,J,L,M,N,EN,LL,MM,NM,IGH,IM1,ITN,ITS,LOW,MP1,ENM1,IERR
  REAL HR(NM,*),HI(NM,*),WR(*),WI(*)
  REAL SI,SR,TI,TR,XI,XR,YI,YR,ZZI,ZZR,S1,S2
!
!***FIRST EXECUTABLE STATEMENT  COMLR
  IERR = 0
!     .......... STORE ROOTS ISOLATED BY CBAL ..........
  DO 200 I = 1, N
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
     DO 440 J = IM1, EN
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
     DO 500 J = I, EN
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
!     .......... INTERCHANGE COLUMNS OF HR AND HI,
!                if NECESSARY ..........
     if (WR(J)  <=  0.0E0) go to 580
!
     DO 540 I = L, J
        ZZR = HR(I,J-1)
        HR(I,J-1) = HR(I,J)
        HR(I,J) = ZZR
        ZZI = HI(I,J-1)
        HI(I,J-1) = HI(I,J)
        HI(I,J) = ZZI
  540    CONTINUE
!
  580    DO 600 I = L, J
        HR(I,J-1) = HR(I,J-1) + XR * HR(I,J) - XI * HI(I,J)
        HI(I,J-1) = HI(I,J-1) + XR * HI(I,J) + XI * HR(I,J)
  600    CONTINUE
!
  640 CONTINUE
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
