subroutine COST (N, X, WSAVE)
!
!! COST computes the cosine transform of a real, even sequence.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A3
!***TYPE      SINGLE PRECISION (COST-S)
!***KEYWORDS  COSINE FOURIER TRANSFORM, FFTPACK
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!  Subroutine COST computes the discrete Fourier cosine transform
!  of an even sequence X(I).  The transform is defined below at output
!  parameter X.
!
!  COST is the unnormalized inverse of itself since a call of COST
!  followed by another call of COST will multiply the input sequence
!  X by 2*(N-1).  The transform is defined below at output parameter X.
!
!  The array WSAVE which is used by subroutine COST must be
!  initialized by calling subroutine COSTI(N,WSAVE).
!
!  Input Parameters
!
!  N       the length of the sequence X.  N must be greater than 1.
!          The method is most efficient when N-1 is a product of
!          small primes.
!
!  X       an array which contains the sequence to be transformed
!
!  WSAVE   a work array which must be dimensioned at least 3*N+15
!          in the program that calls COST.  The WSAVE array must be
!          initialized by calling subroutine COSTI(N,WSAVE), and a
!          different WSAVE array must be used for each different
!          value of N.  This initialization does not have to be
!          repeated so long as N remains unchanged.  Thus subsequent
!          transforms can be obtained faster than the first.
!
!  Output Parameters
!
!  X       For I=1,...,N
!
!             X(I) = X(1)+(-1)**(I-1)*X(N)
!
!               + the sum from K=2 to K=N-1
!
!                 2*X(K)*COS((K-1)*(I-1)*PI/(N-1))
!
!               A call of COST followed by another call of
!               COST will multiply the sequence X by 2*(N-1).
!               Hence COST is the unnormalized inverse
!               of itself.
!
!  WSAVE   contains initialization calculations which must not be
!          destroyed between calls of COST.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  RFFTF
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*)
!   861211  REVISION DATE from Version 3.2
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  COST
  DIMENSION X(*), WSAVE(*)
!***FIRST EXECUTABLE STATEMENT  COST
  NM1 = N-1
  NP1 = N+1
  NS2 = N/2
  if (N-2) 106,101,102
  101 X1H = X(1)+X(2)
  X(2) = X(1)-X(2)
  X(1) = X1H
  return
  102 if (N  >  3) go to 103
  X1P3 = X(1)+X(3)
  TX2 = X(2)+X(2)
  X(2) = X(1)-X(3)
  X(1) = X1P3+TX2
  X(3) = X1P3-TX2
  return
  103 C1 = X(1)-X(N)
  X(1) = X(1)+X(N)
  DO 104 K=2,NS2
     KC = NP1-K
     T1 = X(K)+X(KC)
     T2 = X(K)-X(KC)
     C1 = C1+WSAVE(KC)*T2
     T2 = WSAVE(K)*T2
     X(K) = T1-T2
     X(KC) = T1+T2
  104 CONTINUE
  MODN = MOD(N,2)
  if (MODN  /=  0) X(NS2+1) = X(NS2+1)+X(NS2+1)
  call RFFTF (NM1,X,WSAVE(N+1))
  XIM2 = X(2)
  X(2) = C1
  DO 105 I=4,N,2
     XI = X(I)
     X(I) = X(I-2)-X(I-1)
     X(I-1) = XIM2
     XIM2 = XI
  105 CONTINUE
  if (MODN  /=  0) X(N) = XIM2
  106 RETURN
end
