subroutine PROC (ND, BD, NM1, BM1, NM2, BM2, NA, AA, X, Y, M, A, &
     B, C, D, W, U)
!
!! PROC is subsidiary to CBLKTR.
!
!***LIBRARY   SLATEC
!***TYPE      COMPLEX (PROD-S, PROC-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
! PROC applies a sequence of matrix operations to the vector X and
!  stores the result in Y.
! BD,BM1,BM2 are arrays containing roots of certain B polynomials.
! ND,NM1,NM2 are the lengths of the arrays BD,BM1,BM2 respectively.
! AA         Array containing scalar multipliers of the vector X.
! NA         is the length of the array AA.
! X,Y        The matrix operations are applied to X and the result is Y.
! A,B,C      are arrays which contain the tridiagonal matrix.
! M          is the order of the matrix.
! D,W,U      are working arrays.
! IS         determines whether or not a change in sign is made.
!
!***SEE ALSO  CBLKTR
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PROC
!
  DIMENSION       A(*)       ,B(*)       ,C(*)       ,X(*)       , &
                  Y(*)       ,D(*)       ,W(*)       ,BD(*)      , &
                  BM1(*)     ,BM2(*)     ,AA(*)      ,U(*)
  COMPLEX         X          ,Y          ,A          ,B          , &
                  C          ,D          ,W          ,U          , &
                  DEN
!***FIRST EXECUTABLE STATEMENT  PROC
  DO 101 J=1,M
     W(J) = X(J)
     Y(J) = W(J)
  101 CONTINUE
  MM = M-1
  ID = ND
  IBR = 0
  M1 = NM1
  M2 = NM2
  IA = NA
  102 if (IA) 105,105,103
  103 RT = AA(IA)
  if (ND  ==  0) RT = -RT
  IA = IA-1
!
! SCALAR MULTIPLICATION
!
  DO 104 J=1,M
     Y(J) = RT*W(J)
  104 CONTINUE
  105 if (ID) 125,125,106
  106 RT = BD(ID)
  ID = ID-1
  if (ID  ==  0) IBR = 1
!
! BEGIN SOLUTION TO SYSTEM
!
  D(M) = A(M)/(B(M)-RT)
  W(M) = Y(M)/(B(M)-RT)
  DO 107 J=2,MM
     K = M-J
     DEN = B(K+1)-RT-C(K+1)*D(K+2)
     D(K+1) = A(K+1)/DEN
     W(K+1) = (Y(K+1)-C(K+1)*W(K+2))/DEN
  107 CONTINUE
  DEN = B(1)-RT-C(1)*D(2)
  W(1) = (1.,0.)
  if (ABS(DEN)) 108,109,108
  108 W(1) = (Y(1)-C(1)*W(2))/DEN
  109 DO 110 J=2,M
     W(J) = W(J)-D(J)*W(J-1)
  110 CONTINUE
  if (NA) 113,113,102
  111 DO 112 J=1,M
     Y(J) = W(J)
  112 CONTINUE
  IBR = 1
  go to 102
  113 if (M1) 114,114,115
  114 if (M2) 111,111,120
  115 if (M2) 117,117,116
  116 if (ABS(BM1(M1))-ABS(BM2(M2))) 120,120,117
  117 if (IBR) 118,118,119
  118 if (ABS(BM1(M1)-BD(ID))-ABS(BM1(M1)-RT)) 111,119,119
  119 RT = RT-BM1(M1)
  M1 = M1-1
  go to 123
  120 if (IBR) 121,121,122
  121 if (ABS(BM2(M2)-BD(ID))-ABS(BM2(M2)-RT)) 111,122,122
  122 RT = RT-BM2(M2)
  M2 = M2-1
  123 DO 124 J=1,M
     Y(J) = Y(J)+RT*W(J)
  124 CONTINUE
  go to 102
  125 RETURN
end
