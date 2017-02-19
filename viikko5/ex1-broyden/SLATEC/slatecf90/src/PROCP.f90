subroutine PROCP (ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, Y, M, A, &
     B, C, D, U, W)
!
!! PROCP is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      COMPLEX (PRODP-C, PROCP-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! PROCP applies a sequence of matrix operations to the vector X and
! stores the result in Y (periodic boundary conditions).
!
! BD,BM1,BM2 are arrays containing roots of certain B polynomials.
! ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
! AA         Array containing scalar multipliers of the vector X.
! NA         is the length of the array AA.
! X,Y        The matrix operations are applied to X and the result is Y.
! A,B,C      are arrays which contain the tridiagonal matrix.
! M          is the order of the matrix.
! D,U,W      are working arrays.
! IS         determines whether or not a change in sign is made.
!
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PROCP
!
  DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       , &
                  Y(*)       ,D(*)       ,U(*)       ,BD(*)      , &
                  BM1(*)     ,BM2(*)     ,AA(*)      ,W(*)
  COMPLEX         X          ,Y          ,A          ,B          , &
                  C          ,D          ,U          ,W          , &
                  DEN        ,YM         ,V          ,BH         ,AM
!***FIRST EXECUTABLE STATEMENT  PROCP
  DO 101 J=1,M
     Y(J) = X(J)
     W(J) = Y(J)
  101 CONTINUE
  MM = M-1
  MM2 = M-2
  ID = ND
  IBR = 0
  M1 = NM1
  M2 = NM2
  IA = NA
  102 if (IA) 105,105,103
  103 RT = AA(IA)
  if (ND  ==  0) RT = -RT
  IA = IA-1
  DO 104 J=1,M
     Y(J) = RT*W(J)
  104 CONTINUE
  105 if (ID) 128,128,106
  106 RT = BD(ID)
  ID = ID-1
  if (ID  ==  0) IBR = 1
!
! BEGIN SOLUTION TO SYSTEM
!
  BH = B(M)-RT
  YM = Y(M)
  DEN = B(1)-RT
  D(1) = C(1)/DEN
  U(1) = A(1)/DEN
  W(1) = Y(1)/DEN
  V = C(M)
  if (MM2-2) 109,107,107
  107 DO 108 J=2,MM2
     DEN = B(J)-RT-A(J)*D(J-1)
     D(J) = C(J)/DEN
     U(J) = -A(J)*U(J-1)/DEN
     W(J) = (Y(J)-A(J)*W(J-1))/DEN
     BH = BH-V*U(J-1)
     YM = YM-V*W(J-1)
     V = -V*D(J-1)
  108 CONTINUE
  109 DEN = B(M-1)-RT-A(M-1)*D(M-2)
  D(M-1) = (C(M-1)-A(M-1)*U(M-2))/DEN
  W(M-1) = (Y(M-1)-A(M-1)*W(M-2))/DEN
  AM = A(M)-V*D(M-2)
  BH = BH-V*U(M-2)
  YM = YM-V*W(M-2)
  DEN = BH-AM*D(M-1)
  if (ABS(DEN)) 110,111,110
  110 W(M) = (YM-AM*W(M-1))/DEN
  go to 112
  111 W(M) = (1.,0.)
  112 W(M-1) = W(M-1)-D(M-1)*W(M)
  DO 113 J=2,MM
     K = M-J
     W(K) = W(K)-D(K)*W(K+1)-U(K)*W(M)
  113 CONTINUE
  if (NA) 116,116,102
  114 DO 115 J=1,M
     Y(J) = W(J)
  115 CONTINUE
  IBR = 1
  go to 102
  116 if (M1) 117,117,118
  117 if (M2) 114,114,123
  118 if (M2) 120,120,119
  119 if (ABS(BM1(M1))-ABS(BM2(M2))) 123,123,120
  120 if (IBR) 121,121,122
  121 if (ABS(BM1(M1)-BD(ID))-ABS(BM1(M1)-RT)) 114,122,122
  122 RT = RT-BM1(M1)
  M1 = M1-1
  go to 126
  123 if (IBR) 124,124,125
  124 if (ABS(BM2(M2)-BD(ID))-ABS(BM2(M2)-RT)) 114,125,125
  125 RT = RT-BM2(M2)
  M2 = M2-1
  126 DO 127 J=1,M
     Y(J) = Y(J)+RT*W(J)
  127 CONTINUE
  go to 102
  128 RETURN
end
