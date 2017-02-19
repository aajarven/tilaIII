subroutine RFFTF1 (N, C, CH, WA, IFAC)
!
!! RFFTF1 computes the forward transform of a real, periodic sequence.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A1
!***TYPE      SINGLE PRECISION (RFFTF1-S, CFFTF1-C)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!   Subroutine RFFTF1 computes the Fourier coefficients of a real
!   periodic sequence (Fourier analysis).  The transform is defined
!   below at output parameter C.
!
!   The arrays WA and IFAC which are used by subroutine RFFTB1 must be
!   initialized by calling subroutine RFFTI1.
!
!   Input Arguments
!
!   N       the length of the array R to be transformed.  The method
!           is most efficient when N is a product of small primes.
!           N may change so long as different work arrays are provided.
!
!   C       a real array of length N which contains the sequence
!           to be transformed.
!
!   CH      a real work array of length at least N.
!
!   WA      a real work array which must be dimensioned at least N.
!
!   IFAC    an integer work array which must be dimensioned at least 15.
!
!           The WA and IFAC arrays must be initialized by calling
!           subroutine RFFTI1, and different WA and IFAC arrays must be
!           used for each different value of N.  This initialization
!           does not have to be repeated so long as N remains unchanged.
!           Thus subsequent transforms can be obtained faster than the
!           first.  The same WA and IFAC arrays can be used by RFFTF1
!           and RFFTB1.
!
!   Output Argument
!
!   C       C(1) = the sum from I=1 to I=N of R(I)
!
!           If N is even set L = N/2; if N is odd set L = (N+1)/2
!
!             then for K = 2,...,L
!
!                C(2*K-2) = the sum from I = 1 to I = N of
!
!                     C(I)*COS((K-1)*(I-1)*2*PI/N)
!
!                C(2*K-1) = the sum from I = 1 to I = N of
!
!                    -C(I)*SIN((K-1)*(I-1)*2*PI/N)
!
!           If N is even
!
!                C(N) = the sum from I = 1 to I = N of
!
!                     (-1)**(I-1)*C(I)
!
!   Notes:  This transform is unnormalized since a call of RFFTF1
!           followed by a call of RFFTB1 will multiply the input
!           sequence by N.
!
!           WA and IFAC contain initialization calculations which must
!           not be destroyed between calls of subroutine RFFTF1 or
!           RFFTB1.
!
!***REFERENCES  P. N. Swarztrauber, Vectorizing the FFTs, in Parallel
!                 Computations (G. Rodrigue, ed.), Academic Press,
!                 1982, pp. 51-83.
!***ROUTINES CALLED  RADF2, RADF3, RADF4, RADF5, RADFG
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900131  Routine changed from subsidiary to user-callable.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RFFTF1
  DIMENSION CH(*), C(*), WA(*), IFAC(*)
!***FIRST EXECUTABLE STATEMENT  RFFTF1
  NF = IFAC(2)
  NA = 1
  L2 = N
  IW = N
  DO 111 K1=1,NF
     KH = NF-K1
     IP = IFAC(KH+3)
     L1 = L2/IP
     IDO = N/L2
     IDL1 = IDO*L1
     IW = IW-(IP-1)*IDO
     NA = 1-NA
     if (IP  /=  4) go to 102
     IX2 = IW+IDO
     IX3 = IX2+IDO
     if (NA  /=  0) go to 101
     call RADF4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
     go to 110
  101    call RADF4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
     go to 110
  102    if (IP  /=  2) go to 104
     if (NA  /=  0) go to 103
     call RADF2 (IDO,L1,C,CH,WA(IW))
     go to 110
  103    call RADF2 (IDO,L1,CH,C,WA(IW))
     go to 110
  104    if (IP  /=  3) go to 106
     IX2 = IW+IDO
     if (NA  /=  0) go to 105
     call RADF3 (IDO,L1,C,CH,WA(IW),WA(IX2))
     go to 110
  105    call RADF3 (IDO,L1,CH,C,WA(IW),WA(IX2))
     go to 110
  106    if (IP  /=  5) go to 108
     IX2 = IW+IDO
     IX3 = IX2+IDO
     IX4 = IX3+IDO
     if (NA  /=  0) go to 107
     call RADF5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     go to 110
  107    call RADF5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     go to 110
  108    if (IDO  ==  1) NA = 1-NA
     if (NA  /=  0) go to 109
     call RADFG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
     NA = 1
     go to 110
  109    call RADFG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
     NA = 0
  110    L2 = L1
  111 CONTINUE
  if (NA  ==  1) RETURN
  DO 112 I=1,N
     C(I) = CH(I)
  112 CONTINUE
  return
end
