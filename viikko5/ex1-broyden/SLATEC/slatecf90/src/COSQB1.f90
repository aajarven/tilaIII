subroutine COSQB1 (N, X, W, XH)
!
!! COSQB1 computes the unnormalized inverse of COSQF1.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A3
!***TYPE      SINGLE PRECISION (COSQB1-S)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine COSQB1 computes the fast Fourier transform of quarter
!  wave data. That is, COSQB1 computes a sequence from its
!  representation in terms of a cosine series with odd wave numbers.
!  The transform is defined below at output parameter X.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  RFFTB
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  COSQB1
  DIMENSION X(*), W(*), XH(*)
!***FIRST EXECUTABLE STATEMENT  COSQB1
  NS2 = (N+1)/2
  NP2 = N+2
  DO 101 I=3,N,2
     XIM1 = X(I-1)+X(I)
     X(I) = X(I)-X(I-1)
     X(I-1) = XIM1
  101 CONTINUE
  X(1) = X(1)+X(1)
  MODN = MOD(N,2)
  if (MODN  ==  0) X(N) = X(N)+X(N)
  call RFFTB (N,X,XH)
  DO 102 K=2,NS2
     KC = NP2-K
     XH(K) = W(K-1)*X(KC)+W(KC-1)*X(K)
     XH(KC) = W(K-1)*X(K)-W(KC-1)*X(KC)
  102 CONTINUE
  if (MODN  ==  0) X(NS2+1) = W(NS2)*(X(NS2+1)+X(NS2+1))
  DO 103 K=2,NS2
     KC = NP2-K
     X(K) = XH(K)+XH(KC)
     X(KC) = XH(K)-XH(KC)
  103 CONTINUE
  X(1) = X(1)+X(1)
  return
end
