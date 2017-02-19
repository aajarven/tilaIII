subroutine RADB4 (IDO, L1, CC, CH, WA1, WA2, WA3)
!
!! RADB4 calculates the fast Fourier transform of subvectors of length four.
!
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (RADB4-S)
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing definition of variable SQRT2 by using
!               FORTRAN intrinsic function SQRT instead of a DATA
!               statement.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  RADB4
  DIMENSION CC(IDO,4,*), CH(IDO,L1,4), WA1(*), WA2(*), WA3(*)
!***FIRST EXECUTABLE STATEMENT  RADB4
  SQRT2 = SQRT(2.)
  DO 101 K=1,L1
     TR1 = CC(1,1,K)-CC(IDO,4,K)
     TR2 = CC(1,1,K)+CC(IDO,4,K)
     TR3 = CC(IDO,2,K)+CC(IDO,2,K)
     TR4 = CC(1,3,K)+CC(1,3,K)
     CH(1,K,1) = TR2+TR3
     CH(1,K,2) = TR1-TR4
     CH(1,K,3) = TR2-TR3
     CH(1,K,4) = TR1+TR4
  101 CONTINUE
  if (IDO-2) 107,105,102
  102 IDP2 = IDO+2
  if ( (IDO-1)/2 < L1) go to 108
  DO 104 K=1,L1
!DIR$ IVDEP
     DO 103 I=3,IDO,2
        IC = IDP2-I
        TI1 = CC(I,1,K)+CC(IC,4,K)
        TI2 = CC(I,1,K)-CC(IC,4,K)
        TI3 = CC(I,3,K)-CC(IC,2,K)
        TR4 = CC(I,3,K)+CC(IC,2,K)
        TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
        TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
        TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
        TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1-TR4
        CR4 = TR1+TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
        CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
        CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
        CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
        CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
        CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  103    CONTINUE
  104 CONTINUE
  go to 111
  108 DO 110 I=3,IDO,2
     IC = IDP2-I
!DIR$ IVDEP
     DO 109 K=1,L1
        TI1 = CC(I,1,K)+CC(IC,4,K)
        TI2 = CC(I,1,K)-CC(IC,4,K)
        TI3 = CC(I,3,K)-CC(IC,2,K)
        TR4 = CC(I,3,K)+CC(IC,2,K)
        TR1 = CC(I-1,1,K)-CC(IC-1,4,K)
        TR2 = CC(I-1,1,K)+CC(IC-1,4,K)
        TI4 = CC(I-1,3,K)-CC(IC-1,2,K)
        TR3 = CC(I-1,3,K)+CC(IC-1,2,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1-TR4
        CR4 = TR1+TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-2)*CR2-WA1(I-1)*CI2
        CH(I,K,2) = WA1(I-2)*CI2+WA1(I-1)*CR2
        CH(I-1,K,3) = WA2(I-2)*CR3-WA2(I-1)*CI3
        CH(I,K,3) = WA2(I-2)*CI3+WA2(I-1)*CR3
        CH(I-1,K,4) = WA3(I-2)*CR4-WA3(I-1)*CI4
        CH(I,K,4) = WA3(I-2)*CI4+WA3(I-1)*CR4
  109    CONTINUE
  110 CONTINUE
  111 if (MOD(IDO,2)  ==  1) RETURN
  105 DO 106 K=1,L1
     TI1 = CC(1,2,K)+CC(1,4,K)
     TI2 = CC(1,4,K)-CC(1,2,K)
     TR1 = CC(IDO,1,K)-CC(IDO,3,K)
     TR2 = CC(IDO,1,K)+CC(IDO,3,K)
     CH(IDO,K,1) = TR2+TR2
     CH(IDO,K,2) = SQRT2*(TR1-TI1)
     CH(IDO,K,3) = TI2+TI2
     CH(IDO,K,4) = -SQRT2*(TR1+TI1)
  106 CONTINUE
  107 RETURN
end
