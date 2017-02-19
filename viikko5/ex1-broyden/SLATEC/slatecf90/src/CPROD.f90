subroutine CPROD (ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, YY, M, A, &
     B, C, D, W, Y)
!
!! CPROD is subsidiary to BLKTRI.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (CPROD-S, CPROC-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! PROD applies a sequence of matrix operations to the vector X and
! stores the result in YY.   (COMPLEX case)
! AA         array containing scalar multipliers of the vector X.
! ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
! BD,BM1,BM2 are arrays containing roots of certain B polynomials.
! NA         is the length of the array AA.
! X,YY      The matrix operations are applied to X and the result is YY.
! A,B,C      are arrays which contain the tridiagonal matrix.
! M          is the order of the matrix.
! D,W,Y      are working arrays.
! ISGN       determines whether or not a change in sign is made.
!
!***SEE ALSO  BLKTRI
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  CPROD
!
  COMPLEX         Y          ,D          ,W          ,BD         , &
                  CRT        ,DEN        ,Y1         ,Y2
  DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       , &
                  Y(*)       ,D(*)       ,W(*)       ,BD(*)      , &
                  BM1(*)     ,BM2(*)     ,AA(*)      ,YY(*)
!***FIRST EXECUTABLE STATEMENT  CPROD
  DO 101 J=1,M
     Y(J) = CMPLX(X(J),0.)
  101 CONTINUE
  MM = M-1
  ID = ND
  M1 = NM1
  M2 = NM2
  IA = NA
  102 IFLG = 0
  if (ID) 109,109,103
  103 CRT = BD(ID)
  ID = ID-1
!
! BEGIN SOLUTION TO SYSTEM
!
  D(M) = A(M)/(B(M)-CRT)
  W(M) = Y(M)/(B(M)-CRT)
  DO 104 J=2,MM
     K = M-J
     DEN = B(K+1)-CRT-C(K+1)*D(K+2)
     D(K+1) = A(K+1)/DEN
     W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
  104 CONTINUE
  DEN = B(1)-CRT-C(1)*D(2)
  if (ABS(DEN)) 105,106,105
  105 Y(1) = (Y(1)-C(1)*W(2))/DEN
  go to 107
  106 Y(1) = (1.,0.)
  107 DO 108 J=2,M
     Y(J) = W(J)-D(J)*Y(J-1)
  108 CONTINUE
  109 if (M1) 110,110,112
  110 if (M2) 121,121,111
  111 RT = BM2(M2)
  M2 = M2-1
  go to 117
  112 if (M2) 113,113,114
  113 RT = BM1(M1)
  M1 = M1-1
  go to 117
  114 if (ABS(BM1(M1))-ABS(BM2(M2))) 116,116,115
  115 RT = BM1(M1)
  M1 = M1-1
  go to 117
  116 RT = BM2(M2)
  M2 = M2-1
  117 Y1 = (B(1)-RT)*Y(1)+C(1)*Y(2)
  if (MM-2) 120,118,118
!
! MATRIX MULTIPLICATION
!
  118 DO 119 J=2,MM
     Y2 = A(J)*Y(J-1)+(B(J)-RT)*Y(J)+C(J)*Y(J+1)
     Y(J-1) = Y1
     Y1 = Y2
  119 CONTINUE
  120 Y(M) = A(M)*Y(M-1)+(B(M)-RT)*Y(M)
  Y(M-1) = Y1
  IFLG = 1
  go to 102
  121 if (IA) 124,124,122
  122 RT = AA(IA)
  IA = IA-1
  IFLG = 1
!
! SCALAR MULTIPLICATION
!
  DO 123 J=1,M
     Y(J) = RT*Y(J)
  123 CONTINUE
  124 if (IFLG) 125,125,102
  125 DO 126 J=1,M
     YY(J) = REAL(Y(J))
  126 CONTINUE
  return
end
