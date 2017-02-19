subroutine CPQR79 (NDEG, COEFF, ROOT, IERR, WORK)
!
!! CPQR79 finds the zeros of a polynomial with complex coefficients.
!
!***LIBRARY   SLATEC
!***CATEGORY  F1A1B
!***TYPE      COMPLEX (RPQR79-S, CPQR79-C)
!***KEYWORDS  COMPLEX POLYNOMIAL, POLYNOMIAL ROOTS, POLYNOMIAL ZEROS
!***AUTHOR  Vandevender, W. H., (SNLA)
!***DESCRIPTION
!
!   Abstract
!       This routine computes all zeros of a polynomial of degree NDEG
!       with complex coefficients by computing the eigenvalues of the
!       companion matrix.
!
!   Description of Parameters
!       The user must dimension all arrays appearing in the call list
!            COEFF(NDEG+1), ROOT(NDEG), WORK(2*NDEG*(NDEG+1))
!
!    --Input--
!      NDEG    degree of polynomial
!
!      COEFF   COMPLEX coefficients in descending order.  i.e.,
!              P(Z)= COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)
!
!      WORK    REAL work array of dimension at least 2*NDEG*(NDEG+1)
!
!   --Output--
!      ROOT    COMPLEX vector of roots
!
!      IERR    Output Error Code
!           - Normal Code
!          0  means the roots were computed.
!           - Abnormal Codes
!          1  more than 30 QR iterations on some eigenvalue of the
!             companion matrix
!          2  COEFF(1)=0.0
!          3  NDEG is invalid (less than or equal to 0)
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  COMQR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   791201  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   911010  Code reworked and simplified.  (RWC and WRB)
!***END PROLOGUE  CPQR79
  COMPLEX COEFF(*), ROOT(*), SCALE, C
  REAL WORK(*)
  INTEGER NDEG, IERR, K, KHR, KHI, KWR, KWI, KAD, KJ
!***FIRST EXECUTABLE STATEMENT  CPQR79
  IERR = 0
  if (ABS(COEFF(1))  ==  0.0) THEN
     IERR = 2
     call XERMSG ('SLATEC', 'CPQR79', &
        'LEADING COEFFICIENT IS ZERO.', 2, 1)
     return
  end if
!
  if (NDEG  <=  0) THEN
     IERR = 3
     call XERMSG ('SLATEC', 'CPQR79', 'DEGREE INVALID.', 3, 1)
     return
  end if
!
  if (NDEG  ==  1) THEN
     ROOT(1) = -COEFF(2)/COEFF(1)
     return
  end if
!
  SCALE = 1.0E0/COEFF(1)
  KHR = 1
  KHI = KHR+NDEG*NDEG
  KWR = KHI+KHI-KHR
  KWI = KWR+NDEG
!
  DO 10 K=1,KWR
     WORK(K) = 0.0E0
   10 CONTINUE
!
  DO 20 K=1,NDEG
     KAD = (K-1)*NDEG+1
     C = SCALE*COEFF(K+1)
     WORK(KAD) = -REAL(C)
     KJ = KHI+KAD-1
     WORK(KJ) = -AIMAG(C)
     if (K  /=  NDEG) WORK(KAD+K) = 1.0E0
   20 CONTINUE
!
  call COMQR (NDEG,NDEG,1,NDEG,WORK(KHR),WORK(KHI),WORK(KWR), &
     WORK(KWI),IERR)
!
  if (IERR  /=  0) THEN
     IERR = 1
     call XERMSG ('SLATEC', 'CPQR79', &
        'NO CONVERGENCE IN 30 QR ITERATIONS.', 1, 1)
     return
  end if
!
  DO 30 K=1,NDEG
     KM1 = K-1
     ROOT(K) = CMPLX(WORK(KWR+KM1),WORK(KWI+KM1))
   30 CONTINUE
  return
end
