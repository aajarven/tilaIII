subroutine COSQF1 (N, X, W, XH)
!
!! COSQF1 computes the forward cosine transform with odd wave numbers.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A3
!***TYPE      SINGLE PRECISION (COSQF1-S)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine COSQF1 computes the fast Fourier transform of quarter
!  wave data. That is, COSQF1 computes the coefficients in a cosine
!  series representation with only odd wave numbers.  The transform
!  is defined below at Output Parameter X
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  RFFTF
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  COSQF1
  DIMENSION X(*), W(*), XH(*)
!***FIRST EXECUTABLE STATEMENT  COSQF1
  NS2 = (N+1)/2
  NP2 = N+2
  DO 101 K=2,NS2
     KC = NP2-K
     XH(K) = X(K)+X(KC)
     XH(KC) = X(K)-X(KC)
  101 CONTINUE
  MODN = MOD(N,2)
  if (MODN  ==  0) XH(NS2+1) = X(NS2+1)+X(NS2+1)
  DO 102 K=2,NS2
     KC = NP2-K
     X(K) = W(K-1)*XH(KC)+W(KC-1)*XH(K)
     X(KC) = W(K-1)*XH(K)-W(KC-1)*XH(KC)
  102 CONTINUE
  if (MODN  ==  0) X(NS2+1) = W(NS2)*XH(NS2+1)
  call RFFTF (N,X,XH)
  DO 103 I=3,N,2
     XIM1 = X(I-1)-X(I)
     X(I) = X(I-1)+X(I)
     X(I-1) = XIM1
  103 CONTINUE
  return
end
