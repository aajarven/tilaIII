subroutine PASSB4 (IDO, L1, CC, CH, WA1, WA2, WA3)
!
!! PASSB4 calculates the fast Fourier transform of subvectors of length four.
!
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (PASSB4-S)
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           changing dummy array size declarations (1) to (*).
!   881128  Modified by Dick Valent to meet prologue standards.
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PASSB4
  DIMENSION CC(IDO,4,*), CH(IDO,L1,4), WA1(*), WA2(*), WA3(*)
!***FIRST EXECUTABLE STATEMENT  PASSB4
  if (IDO  /=  2) go to 102
  DO 101 K=1,L1
     TI1 = CC(2,1,K)-CC(2,3,K)
     TI2 = CC(2,1,K)+CC(2,3,K)
     TR4 = CC(2,4,K)-CC(2,2,K)
     TI3 = CC(2,2,K)+CC(2,4,K)
     TR1 = CC(1,1,K)-CC(1,3,K)
     TR2 = CC(1,1,K)+CC(1,3,K)
     TI4 = CC(1,2,K)-CC(1,4,K)
     TR3 = CC(1,2,K)+CC(1,4,K)
     CH(1,K,1) = TR2+TR3
     CH(1,K,3) = TR2-TR3
     CH(2,K,1) = TI2+TI3
     CH(2,K,3) = TI2-TI3
     CH(1,K,2) = TR1+TR4
     CH(1,K,4) = TR1-TR4
     CH(2,K,2) = TI1+TI4
     CH(2,K,4) = TI1-TI4
  101 CONTINUE
  return
  102 if ( IDO/2 < L1) go to 105
  DO 104 K=1,L1
!DIR$ IVDEP
     DO 103 I=2,IDO,2
        TI1 = CC(I,1,K)-CC(I,3,K)
        TI2 = CC(I,1,K)+CC(I,3,K)
        TI3 = CC(I,2,K)+CC(I,4,K)
        TR4 = CC(I,4,K)-CC(I,2,K)
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)
        TI4 = CC(I-1,2,K)-CC(I-1,4,K)
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1+TR4
        CR4 = TR1-TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
        CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
        CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
        CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
        CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
        CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  103    CONTINUE
  104 CONTINUE
  return
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
     DO 106 K=1,L1
        TI1 = CC(I,1,K)-CC(I,3,K)
        TI2 = CC(I,1,K)+CC(I,3,K)
        TI3 = CC(I,2,K)+CC(I,4,K)
        TR4 = CC(I,4,K)-CC(I,2,K)
        TR1 = CC(I-1,1,K)-CC(I-1,3,K)
        TR2 = CC(I-1,1,K)+CC(I-1,3,K)
        TI4 = CC(I-1,2,K)-CC(I-1,4,K)
        TR3 = CC(I-1,2,K)+CC(I-1,4,K)
        CH(I-1,K,1) = TR2+TR3
        CR3 = TR2-TR3
        CH(I,K,1) = TI2+TI3
        CI3 = TI2-TI3
        CR2 = TR1+TR4
        CR4 = TR1-TR4
        CI2 = TI1+TI4
        CI4 = TI1-TI4
        CH(I-1,K,2) = WA1(I-1)*CR2-WA1(I)*CI2
        CH(I,K,2) = WA1(I-1)*CI2+WA1(I)*CR2
        CH(I-1,K,3) = WA2(I-1)*CR3-WA2(I)*CI3
        CH(I,K,3) = WA2(I-1)*CI3+WA2(I)*CR3
        CH(I-1,K,4) = WA3(I-1)*CR4-WA3(I)*CI4
        CH(I,K,4) = WA3(I-1)*CI4+WA3(I)*CR4
  106    CONTINUE
  107 CONTINUE
  return
end
