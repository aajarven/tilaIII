subroutine RPQR79 (NDEG, COEFF, ROOT, IERR, WORK)
!
!! RPQR79 finds the zeros of a polynomial with real coefficients.
!
!***LIBRARY   SLATEC
!***CATEGORY  F1A1A
!***TYPE      SINGLE PRECISION (RPQR79-S, CPQR79-C)
!***KEYWORDS  COMPLEX POLYNOMIAL, POLYNOMIAL ROOTS, POLYNOMIAL ZEROS
!***AUTHOR  Vandevender, W. H., (SNLA)
!***DESCRIPTION
!
!   Abstract
!       This routine computes all zeros of a polynomial of degree NDEG
!       with real coefficients by computing the eigenvalues of the
!       companion matrix.
!
!   Description of Parameters
!       The user must dimension all arrays appearing in the call list
!            COEFF(NDEG+1), ROOT(NDEG), WORK(NDEG*(NDEG+2))
!
!    --Input--
!      NDEG    degree of polynomial
!
!      COEFF   REAL coefficients in descending order.  i.e.,
!              P(Z)= COEFF(1)*(Z**NDEG) + COEFF(NDEG)*Z + COEFF(NDEG+1)
!
!      WORK    REAL work array of dimension at least NDEG*(NDEG+2)
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
!***ROUTINES CALLED  HQR, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   800601  DATE WRITTEN
!   890505  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   911010  Code reworked and simplified.  (RWC and WRB)
!***END PROLOGUE  RPQR79
  REAL COEFF(*), WORK(*), SCALE
  COMPLEX ROOT(*)
  INTEGER NDEG, IERR, K, KH, KWR, KWI, KCOL
!***FIRST EXECUTABLE STATEMENT  RPQR79
  IERR = 0
  if (ABS(COEFF(1))  ==  0.0) THEN
     IERR = 2
     call XERMSG ('SLATEC', 'RPQR79', &
        'LEADING COEFFICIENT IS ZERO.', 2, 1)
     return
  end if
!
  if (NDEG  <=  0) THEN
     IERR = 3
     call XERMSG ('SLATEC', 'RPQR79', 'DEGREE INVALID.', 3, 1)
     return
  end if
!
  if (NDEG  ==  1) THEN
     ROOT(1) = CMPLX(-COEFF(2)/COEFF(1),0.0)
     return
  end if
!
  SCALE = 1.0E0/COEFF(1)
  KH = 1
  KWR = KH+NDEG*NDEG
  KWI = KWR+NDEG
  KWEND = KWI+NDEG-1
!
  DO 10 K=1,KWEND
     WORK(K) = 0.0E0
   10 CONTINUE
!
  DO 20 K=1,NDEG
     KCOL = (K-1)*NDEG+1
     WORK(KCOL) = -COEFF(K+1)*SCALE
     if (K  /=  NDEG) WORK(KCOL+K) = 1.0E0
   20 CONTINUE
!
  call HQR (NDEG,NDEG,1,NDEG,WORK(KH),WORK(KWR),WORK(KWI),IERR)
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
