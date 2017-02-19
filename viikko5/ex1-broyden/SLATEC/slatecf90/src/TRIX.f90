subroutine TRIX (IDEGBR, IDEGCR, M, A, B, C, Y, TCOS, D, W)
!
!! TRIX is subsidiary to GENBUN.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (TRIX-S, CMPTRX-C)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!     Subroutine to solve a system of linear equations where the
!     coefficient matrix is a rational function in the matrix given by
!     TRIDIAGONAL  ( . . . , A(I), B(I), C(I), . . . ).
!
!***SEE ALSO  GENBUN
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  TRIX
!
  DIMENSION       A(*)       ,B(*)       ,C(*)       ,Y(*)       , &
                  TCOS(*)    ,D(*)       ,W(*)
  INTEGER KB, KC
!***FIRST EXECUTABLE STATEMENT  TRIX
  MM1 = M-1
  KB = IDEGBR+1
  KC = IDEGCR+1
  L = (IDEGBR+1)/(IDEGCR+1)
  LINT = 1
  DO 108 K=1,IDEGBR
     X = TCOS(K)
     if (K  /=  L) go to 102
     I = IDEGBR+LINT
     XX = X-TCOS(I)
     DO 101 I=1,M
        W(I) = Y(I)
        Y(I) = XX*Y(I)
  101    CONTINUE
  102    CONTINUE
     Z = 1./(B(1)-X)
     D(1) = C(1)*Z
     Y(1) = Y(1)*Z
     DO 103 I=2,MM1
        Z = 1./(B(I)-X-A(I)*D(I-1))
        D(I) = C(I)*Z
        Y(I) = (Y(I)-A(I)*Y(I-1))*Z
  103    CONTINUE
     Z = B(M)-X-A(M)*D(MM1)
     if (Z  /=  0.) go to 104
     Y(M) = 0.
     go to 105
  104    Y(M) = (Y(M)-A(M)*Y(MM1))/Z
  105    CONTINUE
     DO 106 IP=1,MM1
        I = M-IP
        Y(I) = Y(I)-D(I)*Y(I+1)
  106    CONTINUE
     if (K  /=  L) go to 108
     DO 107 I=1,M
        Y(I) = Y(I)+W(I)
  107    CONTINUE
     LINT = LINT+1
     L = (LINT*KB)/KC
  108 CONTINUE
  return
end
