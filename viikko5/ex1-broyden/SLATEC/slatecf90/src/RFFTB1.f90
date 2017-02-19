subroutine RFFTB1 (N, C, CH, WA, IFAC)
!
!! RFFTB1 computes the backward fast Fourier transform of a real array.
!
!***LIBRARY   SLATEC (FFTPACK)
!***CATEGORY  J1A1
!***TYPE      SINGLE PRECISION (RFFTB1-S, CFFTB1-C)
!***KEYWORDS  FFTPACK, FOURIER TRANSFORM
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***DESCRIPTION
!
!   Subroutine RFFTB1 computes the real periodic sequence from its
!   Fourier coefficients (Fourier synthesis).  The transform is defined
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
!   C       For N even and for I = 1,...,N
!
!                C(I) = C(1)+(-1)**(I-1)*C(N)
!
!                     plus the sum from K=2 to K=N/2 of
!
!                      2.*C(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
!
!                     -2.*C(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
!
!           For N odd and for I = 1,...,N
!
!                C(I) = C(1) plus the sum from K=2 to K=(N+1)/2 of
!
!                     2.*C(2*K-2)*COS((K-1)*(I-1)*2*PI/N)
!
!                    -2.*C(2*K-1)*SIN((K-1)*(I-1)*2*PI/N)
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
!***ROUTINES CALLED  RADB2, RADB3, RADB4, RADB5, RADBG
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900131  Routine changed from subsidiary to user-callable.  (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  RFFTB1
  DIMENSION CH(*), C(*), WA(*), IFAC(*)
!***FIRST EXECUTABLE STATEMENT  RFFTB1
  NF = IFAC(2)
  NA = 0
  L1 = 1
  IW = 1
  DO 116 K1=1,NF
     IP = IFAC(K1+2)
     L2 = IP*L1
     IDO = N/L2
     IDL1 = IDO*L1
     if (IP  /=  4) go to 103
     IX2 = IW+IDO
     IX3 = IX2+IDO
     if (NA  /=  0) go to 101
     call RADB4 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3))
     go to 102
  101    call RADB4 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3))
  102    NA = 1-NA
     go to 115
  103    if (IP  /=  2) go to 106
     if (NA  /=  0) go to 104
     call RADB2 (IDO,L1,C,CH,WA(IW))
     go to 105
  104    call RADB2 (IDO,L1,CH,C,WA(IW))
  105    NA = 1-NA
     go to 115
  106    if (IP  /=  3) go to 109
     IX2 = IW+IDO
     if (NA  /=  0) go to 107
     call RADB3 (IDO,L1,C,CH,WA(IW),WA(IX2))
     go to 108
  107    call RADB3 (IDO,L1,CH,C,WA(IW),WA(IX2))
  108    NA = 1-NA
     go to 115
  109    if (IP  /=  5) go to 112
     IX2 = IW+IDO
     IX3 = IX2+IDO
     IX4 = IX3+IDO
     if (NA  /=  0) go to 110
     call RADB5 (IDO,L1,C,CH,WA(IW),WA(IX2),WA(IX3),WA(IX4))
     go to 111
  110    call RADB5 (IDO,L1,CH,C,WA(IW),WA(IX2),WA(IX3),WA(IX4))
  111    NA = 1-NA
     go to 115
  112    if (NA  /=  0) go to 113
     call RADBG (IDO,IP,L1,IDL1,C,C,C,CH,CH,WA(IW))
     go to 114
  113    call RADBG (IDO,IP,L1,IDL1,CH,CH,CH,C,C,WA(IW))
  114    if (IDO  ==  1) NA = 1-NA
  115    L1 = L2
     IW = IW+(IP-1)*IDO
  116 CONTINUE
  if (NA  ==  0) RETURN
  DO 117 I=1,N
     C(I) = CH(I)
  117 CONTINUE
  return
end
