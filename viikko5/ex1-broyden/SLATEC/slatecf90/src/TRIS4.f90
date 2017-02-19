subroutine TRIS4 (N, A, B, C, D, U, Z)
!
!! TRIS4 is subsidiary to SEPX4.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to SEPX4
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (TRIS4-S)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     This subroutine solves for a non-zero eigenvector corresponding
!     to the zero eigenvalue of the transpose of the rank
!     deficient ONE matrix with subdiagonal A, diagonal B, and
!     superdiagonal C , with A(1) in the (1,N) position, with
!     C(N) in the (N,1) position, AND all other elements zero.
!
!***SEE ALSO  SEPX4
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  TRIS4
!
  DIMENSION       A(*)       ,B(*)       ,C(*)       ,D(*)       , &
                  U(*)       ,Z(*)
!***FIRST EXECUTABLE STATEMENT  TRIS4
  BN = B(N)
  D(1) = A(2)/B(1)
  V = A(1)
  U(1) = C(N)/B(1)
  NM2 = N-2
  DO  10 J=2,NM2
     DEN = B(J)-C(J-1)*D(J-1)
     D(J) = A(J+1)/DEN
     U(J) = -C(J-1)*U(J-1)/DEN
     BN = BN-V*U(J-1)
     V = -V*D(J-1)
   10 CONTINUE
  DEN = B(N-1)-C(N-2)*D(N-2)
  D(N-1) = (A(N)-C(N-2)*U(N-2))/DEN
  AN = C(N-1)-V*D(N-2)
  BN = BN-V*U(N-2)
  DEN = BN-AN*D(N-1)
!
!     SET LAST COMPONENT EQUAL TO ONE
!
  Z(N) = 1.0
  Z(N-1) = -D(N-1)
  NM1 = N-1
  DO  20 J=2,NM1
     K = N-J
     Z(K) = -D(K)*Z(K+1)-U(K)*Z(N)
   20 CONTINUE
  return
end
