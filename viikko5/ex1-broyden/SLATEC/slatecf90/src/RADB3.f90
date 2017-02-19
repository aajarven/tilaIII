subroutine RADB3 (IDO, L1, CC, CH, WA1, WA2)
!
!! RADB3 calculates the fast Fourier transform of subvectors of length three.
!
!***LIBRARY   SLATEC (FFTPACK)
!***TYPE      SINGLE PRECISION (RADB3-S)
!***AUTHOR  Swarztrauber, P. N., (NCAR)
!***ROUTINES CALLED  (NONE)
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   830401  Modified to use SLATEC library source file format.
!   860115  Modified by Ron Boisvert to adhere to Fortran 77 by
!           (a) changing dummy array size declarations (1) to (*),
!           (b) changing definition of variable TAUI by using
!               FORTRAN intrinsic function SQRT instead of a DATA
!               statement.
!   881128  Modified by Dick Valent to meet prologue standards.
!   890831  Modified array declarations.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900402  Added TYPE section.  (WRB)
!***END PROLOGUE  RADB3
  DIMENSION CC(IDO,3,*), CH(IDO,L1,3), WA1(*), WA2(*)
!***FIRST EXECUTABLE STATEMENT  RADB3
  TAUR = -.5
  TAUI = .5*SQRT(3.)
  DO 101 K=1,L1
     TR2 = CC(IDO,2,K)+CC(IDO,2,K)
     CR2 = CC(1,1,K)+TAUR*TR2
     CH(1,K,1) = CC(1,1,K)+TR2
     CI3 = TAUI*(CC(1,3,K)+CC(1,3,K))
     CH(1,K,2) = CR2-CI3
     CH(1,K,3) = CR2+CI3
  101 CONTINUE
  if (IDO  ==  1) RETURN
  IDP2 = IDO+2
  if ( (IDO-1)/2 < L1) go to 104
  DO 103 K=1,L1
!DIR$ IVDEP
     DO 102 I=3,IDO,2
        IC = IDP2-I
        TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,3,K)-CC(IC,2,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
        CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
        CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
        CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
        CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  102    CONTINUE
  103 CONTINUE
  return
  104 DO 106 I=3,IDO,2
     IC = IDP2-I
!DIR$ IVDEP
     DO 105 K=1,L1
        TR2 = CC(I-1,3,K)+CC(IC-1,2,K)
        CR2 = CC(I-1,1,K)+TAUR*TR2
        CH(I-1,K,1) = CC(I-1,1,K)+TR2
        TI2 = CC(I,3,K)-CC(IC,2,K)
        CI2 = CC(I,1,K)+TAUR*TI2
        CH(I,K,1) = CC(I,1,K)+TI2
        CR3 = TAUI*(CC(I-1,3,K)-CC(IC-1,2,K))
        CI3 = TAUI*(CC(I,3,K)+CC(IC,2,K))
        DR2 = CR2-CI3
        DR3 = CR2+CI3
        DI2 = CI2+CR3
        DI3 = CI2-CR3
        CH(I-1,K,2) = WA1(I-2)*DR2-WA1(I-1)*DI2
        CH(I,K,2) = WA1(I-2)*DI2+WA1(I-1)*DR2
        CH(I-1,K,3) = WA2(I-2)*DR3-WA2(I-1)*DI3
        CH(I,K,3) = WA2(I-2)*DI3+WA2(I-1)*DR3
  105    CONTINUE
  106 CONTINUE
  return
end
