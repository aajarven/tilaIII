subroutine CPRODP (ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, YY, M, &
     A, B, C, D, U, Y)
!
!! CPRODP is subsidiary to BLKTRI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CPRODP-S, CPROCP-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! PRODP applies a sequence of matrix operations to the vector X and
! stores the result in YY. (Periodic boundary conditions and COMPLEX
! case)
!
! BD,BM1,BM2     are arrays containing roots of certain B polynomials.
! ND,NM1,NM2     are the lengths of the arrays BD,BM1,BM2 respectively.
! AA             Array containing scalar multipliers of the vector X.
! NA             is the length of the array AA.
! X,YY      The matrix operations are applied to X and the result is YY.
! A,B,C          are arrays which contain the tridiagonal matrix.
! M              is the order of the matrix.
! D,U,Y          are working arrays.
! ISGN           determines whether or not a change in sign is made.
!
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CPRODP
!
  COMPLEX         Y          ,D          ,U          ,V          , &
                  DEN        ,BH         ,YM         ,AM         , &
                  Y1         ,Y2         ,YH         ,BD         , &
                  CRT
  DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       , &
                  Y(*)       ,D(*)       ,U(*)       ,BD(*)      , &
                  BM1(*)     ,BM2(*)     ,AA(*)      ,YY(*)
!***FIRST EXECUTABLE STATEMENT  CPRODP
  DO 101 J=1,M
     Y(J) = CMPLX(X(J),0.)
  101 CONTINUE
  MM = M-1
  MM2 = M-2
  ID = ND
  M1 = NM1
  M2 = NM2
  IA = NA
  102 IFLG = 0
  if (ID) 111,111,103
  103 CRT = BD(ID)
  ID = ID-1
  IFLG = 1
!
! BEGIN SOLUTION TO SYSTEM
!
  BH = B(M)-CRT
  YM = Y(M)
  DEN = B(1)-CRT
  D(1) = C(1)/DEN
  U(1) = A(1)/DEN
  Y(1) = Y(1)/DEN
  V = CMPLX(C(M),0.)
  if (MM2-2) 106,104,104
  104 DO 105 J=2,MM2
     DEN = B(J)-CRT-A(J)*D(J-1)
     D(J) = C(J)/DEN
     U(J) = -A(J)*U(J-1)/DEN
     Y(J) = (Y(J)-A(J)*Y(J-1))/DEN
     BH = BH-V*U(J-1)
     YM = YM-V*Y(J-1)
     V = -V*D(J-1)
  105 CONTINUE
  106 DEN = B(M-1)-CRT-A(M-1)*D(M-2)
  D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
  Y(M-1) = (Y(M-1)-A(M-1)*Y(M-2))/DEN
  AM = A(M)-V*D(M-2)
  BH = BH-V*U(M-2)
  YM = YM-V*Y(M-2)
  DEN = BH-AM*D(M-1)
  if (ABS(DEN)) 107,108,107
  107 Y(M) = (YM-AM*Y(M-1))/DEN
  go to 109
  108 Y(M) = (1.,0.)
  109 Y(M-1) = Y(M-1)-D(M-1)*Y(M)
  DO 110 J=2,MM
     K = M-J
     Y(K) = Y(K)-D(K)*Y(K+1)-U(K)*Y(M)
  110 CONTINUE
  111 if (M1) 112,112,114
  112 if (M2) 123,123,113
  113 RT = BM2(M2)
  M2 = M2-1
  go to 119
  114 if (M2) 115,115,116
  115 RT = BM1(M1)
  M1 = M1-1
  go to 119
  116 if (ABS(BM1(M1))-ABS(BM2(M2))) 118,118,117
  117 RT = BM1(M1)
  M1 = M1-1
  go to 119
  118 RT = BM2(M2)
  M2 = M2-1
!
! MATRIX MULTIPLICATION
!
  119 YH = Y(1)
  Y1 = (B(1)-RT)*Y(1)+C(1)*Y(2)+A(1)*Y(M)
  if (MM-2) 122,120,120
  120 DO 121 J=2,MM
     Y2 = A(J)*Y(J-1)+(B(J)-RT)*Y(J)+C(J)*Y(J+1)
     Y(J-1) = Y1
     Y1 = Y2
  121 CONTINUE
  122 Y(M) = A(M)*Y(M-1)+(B(M)-RT)*Y(M)+C(M)*YH
  Y(M-1) = Y1
  IFLG = 1
  go to 102
  123 if (IA) 126,126,124
  124 RT = AA(IA)
  IA = IA-1
  IFLG = 1
!
! SCALAR MULTIPLICATION
!
  DO 125 J=1,M
     Y(J) = RT*Y(J)
  125 CONTINUE
  126 if (IFLG) 127,127,102
  127 DO 128 J=1,M
     YY(J) = REAL(Y(J))
  128 CONTINUE
  return
end
