subroutine PASSB5 (IDO, L1, CC, CH, WA1, WA2, WA3, WA4)
!
!! PASSB5 calculates the fast Fourier transform of subvectors of length five.
!
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (PASSB5-S)
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing definition of variables PI, TI11, TI12,
!               TR11, TR12 by using FORTRAN intrinsic functions ATAN
!               and SIN instead of DATA statements.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  PASSB5
  DIMENSION CC(IDO,5,*), CH(IDO,L1,5), WA1(*), WA2(*), WA3(*), &
            WA4(*)
!***FIRST EXECUTABLE STATEMENT  PASSB5
  PI = 4.*ATAN(1.)
  TR11 = SIN(.1*PI)
  TI11 = SIN(.4*PI)
  TR12 = -SIN(.3*PI)
  TI12 = SIN(.2*PI)
  if (IDO  /=  2) go to 102
  DO 101 K=1,L1
     TI5 = CC(2,2,K)-CC(2,5,K)
     TI2 = CC(2,2,K)+CC(2,5,K)
     TI4 = CC(2,3,K)-CC(2,4,K)
     TI3 = CC(2,3,K)+CC(2,4,K)
     TR5 = CC(1,2,K)-CC(1,5,K)
     TR2 = CC(1,2,K)+CC(1,5,K)
     TR4 = CC(1,3,K)-CC(1,4,K)
     TR3 = CC(1,3,K)+CC(1,4,K)
     CH(1,K,1) = CC(1,1,K)+TR2+TR3
     CH(2,K,1) = CC(2,1,K)+TI2+TI3
     CR2 = CC(1,1,K)+TR11*TR2+TR12*TR3
     CI2 = CC(2,1,K)+TR11*TI2+TR12*TI3
     CR3 = CC(1,1,K)+TR12*TR2+TR11*TR3
     CI3 = CC(2,1,K)+TR12*TI2+TR11*TI3
     CR5 = TI11*TR5+TI12*TR4
     CI5 = TI11*TI5+TI12*TI4
     CR4 = TI12*TR5-TI11*TR4
     CI4 = TI12*TI5-TI11*TI4
     CH(1,K,2) = CR2-CI5
     CH(1,K,5) = CR2+CI5
     CH(2,K,2) = CI2+CR5
     CH(2,K,3) = CI3+CR4
     CH(1,K,3) = CR3-CI4
     CH(1,K,4) = CR3+CI4
     CH(2,K,4) = CI3-CR4
     CH(2,K,5) = CI2-CR5
  101 CONTINUE
  return
  102 if ( IDO/2 < L1) go to 105
  DO 104 K=1,L1
!DIR$ IVDEP
     DO 103 I=2,IDO,2
        TI5 = CC(I,2,K)-CC(I,5,K)
        TI2 = CC(I,2,K)+CC(I,5,K)
        TI4 = CC(I,3,K)-CC(I,4,K)
        TI3 = CC(I,3,K)+CC(I,4,K)
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
        CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
        CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
        CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
        CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  103    CONTINUE
  104 CONTINUE
  return
  105 DO 107 I=2,IDO,2
!DIR$ IVDEP
     DO 106 K=1,L1
        TI5 = CC(I,2,K)-CC(I,5,K)
        TI2 = CC(I,2,K)+CC(I,5,K)
        TI4 = CC(I,3,K)-CC(I,4,K)
        TI3 = CC(I,3,K)+CC(I,4,K)
        TR5 = CC(I-1,2,K)-CC(I-1,5,K)
        TR2 = CC(I-1,2,K)+CC(I-1,5,K)
        TR4 = CC(I-1,3,K)-CC(I-1,4,K)
        TR3 = CC(I-1,3,K)+CC(I-1,4,K)
        CH(I-1,K,1) = CC(I-1,1,K)+TR2+TR3
        CH(I,K,1) = CC(I,1,K)+TI2+TI3
        CR2 = CC(I-1,1,K)+TR11*TR2+TR12*TR3
        CI2 = CC(I,1,K)+TR11*TI2+TR12*TI3
        CR3 = CC(I-1,1,K)+TR12*TR2+TR11*TR3
        CI3 = CC(I,1,K)+TR12*TI2+TR11*TI3
        CR5 = TI11*TR5+TI12*TR4
        CI5 = TI11*TI5+TI12*TI4
        CR4 = TI12*TR5-TI11*TR4
        CI4 = TI12*TI5-TI11*TI4
        DR3 = CR3-CI4
        DR4 = CR3+CI4
        DI3 = CI3+CR4
        DI4 = CI3-CR4
        DR5 = CR2+CI5
        DR2 = CR2-CI5
        DI5 = CI2-CR5
        DI2 = CI2+CR5
        CH(I-1,K,2) = WA1(I-1)*DR2-WA1(I)*DI2
        CH(I,K,2) = WA1(I-1)*DI2+WA1(I)*DR2
        CH(I-1,K,3) = WA2(I-1)*DR3-WA2(I)*DI3
        CH(I,K,3) = WA2(I-1)*DI3+WA2(I)*DR3
        CH(I-1,K,4) = WA3(I-1)*DR4-WA3(I)*DI4
        CH(I,K,4) = WA3(I-1)*DI4+WA3(I)*DR4
        CH(I-1,K,5) = WA4(I-1)*DR5-WA4(I)*DI5
        CH(I,K,5) = WA4(I-1)*DI5+WA4(I)*DR5
  106    CONTINUE
  107 CONTINUE
  return
end