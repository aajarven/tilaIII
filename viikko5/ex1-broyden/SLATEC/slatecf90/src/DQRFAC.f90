subroutine DQRFAC (M, N, A, LDA, PIVOT, IPVT, LIPVT, SIGMA, &
     ACNORM, WA)
!
!! DQRFAC computes the QR factorization of a rectangulr matrix.
!
!***SUBSIDIARY
!***PURPOSE  Subsidiary to DNLS1, DNLS1E, DNSQ and DNSQE
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (QRFAC-S, DQRFAC-D)
!***AUTHOR  (UNKNOWN)
!***DESCRIPTION
!
!   **** Double Precision version of QRFAC ****
!
!     This subroutine uses Householder transformations with column
!     pivoting (optional) to compute a QR factorization of the
!     M by N matrix A. That is, DQRFAC determines an orthogonal
!     matrix Q, a permutation matrix P, and an upper trapezoidal
!     matrix R with diagonal elements of nonincreasing magnitude,
!     such that A*P = Q*R. The Householder transformation for
!     column K, K = 1,2,...,MIN(M,N), is of the form
!
!                           T
!           I - (1/U(K))*U*U
!
!     where U has zeros in the first K-1 positions. The form of
!     this transformation and the method of pivoting first
!     appeared in the corresponding LINPACK subroutine.
!
!     The subroutine statement is
!
!       SUBROUTINE DQRFAC(M,N,A,LDA,PIVOT,IPVT,LIPVT,SIGMA,ACNORM,WA)
!
!     where
!
!       M is a positive integer input variable set to the number
!         of rows of A.
!
!       N is a positive integer input variable set to the number
!         of columns of A.
!
!       A is an M by N array. On input A contains the matrix for
!         which the QR factorization is to be computed. On output
!         the strict upper trapezoidal part of A contains the strict
!         upper trapezoidal part of R, and the lower trapezoidal
!         part of A contains a factored form of Q (the non-trivial
!         elements of the U vectors described above).
!
!       LDA is a positive integer input variable not less than M
!         which specifies the leading dimension of the array A.
!
!       PIVOT is a logical input variable. If pivot is set .TRUE.,
!         then column pivoting is enforced. If pivot is set .FALSE.,
!         then no column pivoting is done.
!
!       IPVT is an integer output array of length LIPVT. IPVT
!         defines the permutation matrix P such that A*P = Q*R.
!         Column J of P is column IPVT(J) of the identity matrix.
!         If pivot is .FALSE., IPVT is not referenced.
!
!       LIPVT is a positive integer input variable. If PIVOT is
!             .FALSE., then LIPVT may be as small as 1. If PIVOT is
!             .TRUE., then LIPVT must be at least N.
!
!       SIGMA is an output array of length N which contains the
!         diagonal elements of R.
!
!       ACNORM is an output array of length N which contains the
!         norms of the corresponding columns of the input matrix A.
!         If this information is not needed, then ACNORM can coincide
!         with SIGMA.
!
!       WA is a work array of length N. If pivot is .FALSE., then WA
!         can coincide with SIGMA.
!
!***SEE ALSO  DNLS1, DNLS1E, DNSQ, DNSQE
!***ROUTINES CALLED  D1MACH, DENORM
!***REVISION HISTORY  (YYMMDD)
!   800301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   900328  Added TYPE section.  (WRB)
!***END PROLOGUE  DQRFAC
  INTEGER M,N,LDA,LIPVT
  INTEGER IPVT(*)
  LOGICAL PIVOT
  SAVE ONE, P05, ZERO
  DOUBLE PRECISION A(LDA,*),SIGMA(*),ACNORM(*),WA(*)
  INTEGER I,J,JP1,K,KMAX,MINMN
  DOUBLE PRECISION AJNORM,EPSMCH,ONE,P05,SUM,TEMP,ZERO
  DOUBLE PRECISION D1MACH,DENORM
  DATA ONE,P05,ZERO /1.0D0,5.0D-2,0.0D0/
!***FIRST EXECUTABLE STATEMENT  DQRFAC
  EPSMCH = D1MACH(4)
!
!     COMPUTE THE INITIAL COLUMN NORMS AND INITIALIZE SEVERAL ARRAYS.
!
  DO 10 J = 1, N
     ACNORM(J) = DENORM(M,A(1,J))
     SIGMA(J) = ACNORM(J)
     WA(J) = SIGMA(J)
     if (PIVOT) IPVT(J) = J
   10    CONTINUE
!
!     REDUCE A TO R WITH HOUSEHOLDER TRANSFORMATIONS.
!
  MINMN = MIN(M,N)
  DO 110 J = 1, MINMN
     if (.NOT.PIVOT) go to 40
!
!        BRING THE COLUMN OF LARGEST NORM INTO THE PIVOT POSITION.
!
     KMAX = J
     DO 20 K = J, N
        if (SIGMA(K)  >  SIGMA(KMAX)) KMAX = K
   20       CONTINUE
     if (KMAX  ==  J) go to 40
     DO 30 I = 1, M
        TEMP = A(I,J)
        A(I,J) = A(I,KMAX)
        A(I,KMAX) = TEMP
   30       CONTINUE
     SIGMA(KMAX) = SIGMA(J)
     WA(KMAX) = WA(J)
     K = IPVT(J)
     IPVT(J) = IPVT(KMAX)
     IPVT(KMAX) = K
   40    CONTINUE
!
!        COMPUTE THE HOUSEHOLDER TRANSFORMATION TO REDUCE THE
!        J-TH COLUMN OF A TO A MULTIPLE OF THE J-TH UNIT VECTOR.
!
     AJNORM = DENORM(M-J+1,A(J,J))
     if (AJNORM  ==  ZERO) go to 100
     if (A(J,J)  <  ZERO) AJNORM = -AJNORM
     DO 50 I = J, M
        A(I,J) = A(I,J)/AJNORM
   50       CONTINUE
     A(J,J) = A(J,J) + ONE
!
!        APPLY THE TRANSFORMATION TO THE REMAINING COLUMNS
!        AND UPDATE THE NORMS.
!
     JP1 = J + 1
     if (N  <  JP1) go to 100
     DO 90 K = JP1, N
        SUM = ZERO
        DO 60 I = J, M
           SUM = SUM + A(I,J)*A(I,K)
   60          CONTINUE
        TEMP = SUM/A(J,J)
        DO 70 I = J, M
           A(I,K) = A(I,K) - TEMP*A(I,J)
   70          CONTINUE
        if (.NOT.PIVOT .OR. SIGMA(K)  ==  ZERO) go to 80
        TEMP = A(J,K)/SIGMA(K)
        SIGMA(K) = SIGMA(K)*SQRT(MAX(ZERO,ONE-TEMP**2))
        if (P05*(SIGMA(K)/WA(K))**2  >  EPSMCH) go to 80
        SIGMA(K) = DENORM(M-J,A(JP1,K))
        WA(K) = SIGMA(K)
   80       CONTINUE
   90       CONTINUE
  100    CONTINUE
     SIGMA(J) = -AJNORM
  110    CONTINUE
  return
!
!     LAST CARD OF SUBROUTINE DQRFAC.
!
end
