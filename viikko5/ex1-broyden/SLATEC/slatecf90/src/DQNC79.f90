subroutine DQNC79 (FUN, A, B, ERR, ANS, IERR, K)
!
!! DQNC79 integrates a function using a 7-point adaptive Newton-Cotes rule.
!
!***LIBRARY   SLATEC
!***CATEGORY  H2A1A1
!***TYPE      DOUBLE PRECISION (QNC79-S, DQNC79-D)
!***KEYWORDS  ADAPTIVE QUADRATURE, INTEGRATION, NEWTON-COTES
!***AUTHOR  Kahaner, D. K., (NBS)
!           Jones, R. E., (SNLA)
!***DESCRIPTION
!
!     Abstract  *** a DOUBLE PRECISION routine ***
!       DQNC79 is a general purpose program for evaluation of
!       one dimensional integrals of user defined functions.
!       DQNC79 will pick its own points for evaluation of the
!       integrand and these will vary from problem to problem.
!       Thus, DQNC79 is not designed to integrate over data sets.
!       Moderately smooth integrands will be integrated efficiently
!       and reliably.  For problems with strong singularities,
!       oscillations etc., the user may wish to use more sophis-
!       ticated routines such as those in QUADPACK.  One measure
!       of the reliability of DQNC79 is the output parameter K,
!       giving the number of integrand evaluations that were needed.
!
!     Description of Arguments
!
!     --Input--* FUN, A, B, ERR are DOUBLE PRECISION *
!       FUN  - name of external function to be integrated.  This name
!              must be in an EXTERNAL statement in your calling
!              program.  You must write a Fortran function to evaluate
!              FUN.  This should be of the form
!                    DOUBLE PRECISION FUNCTION FUN (X)
!              C
!              C     X can vary from A to B
!              C     FUN(X) should be finite for all X on interval.
!              C
!                    FUN = ...
!                    return
!                    END
!       A    - lower limit of integration
!       B    - upper limit of integration (may be less than A)
!       ERR  - is a requested error tolerance.  Normally, pick a value
!              0  <  ERR  <  1.0D-8.
!
!     --Output--
!       ANS  - computed value of the integral.  Hopefully, ANS is
!              accurate to within ERR * integral of ABS(FUN(X)).
!       IERR - a status code
!            - Normal codes
!               1  ANS most likely meets requested error tolerance.
!              -1  A equals B, or A and B are too nearly equal to
!                  allow normal integration.  ANS is set to zero.
!            - Abnormal code
!               2  ANS probably does not meet requested error tolerance.
!       K    - the number of function evaluations actually used to do
!              the integration.  A value of K  >  1000 indicates a
!              difficult problem; other programs may be more efficient.
!              DQNC79 will gracefully give up if K exceeds 2000.
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  D1MACH, I1MACH, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   790601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890911  Removed unnecessary intrinsics.  (WRB)
!   890911  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920218  Code redone to parallel QNC79.  (WRB)
!   930120  Increase array size 80->99, and KMX 2000->5000 for SUN -r8
!           wordlength.  (RWC)
!***END PROLOGUE  DQNC79
!     .. Scalar Arguments ..
  DOUBLE PRECISION A, ANS, B, ERR
  INTEGER IERR, K
!     .. Function Arguments ..
  DOUBLE PRECISION FUN
  EXTERNAL FUN
!     .. Local Scalars ..
  DOUBLE PRECISION AE, AREA, BANK, BLOCAL, C, CE, EE, EF, EPS, Q13, &
                   Q7, Q7L, SQ2, TEST, TOL, VR, W1, W2, W3, W4
  INTEGER I, KML, KMX, L, LMN, LMX, NBITS, NIB, NLMN, NLMX
  LOGICAL FIRST
!     .. Local Arrays ..
  DOUBLE PRECISION AA(99), F(13), F1(99), F2(99), F3(99), F4(99), &
                   F5(99), F6(99), F7(99), HH(99), Q7R(99), VL(99)
  INTEGER LR(99)
!     .. External Functions ..
  DOUBLE PRECISION D1MACH
  INTEGER I1MACH
  EXTERNAL D1MACH, I1MACH
!     .. External Subroutines ..
  EXTERNAL XERMSG
!     .. Intrinsic Functions ..
  INTRINSIC ABS, LOG, MAX, MIN, SIGN, SQRT
!     .. Save statement ..
  SAVE NBITS, NLMX, FIRST, SQ2, W1, W2, W3, W4
!     .. Data statements ..
  DATA KML /7/, KMX /5000/, NLMN /2/
  DATA FIRST /.TRUE./
!***FIRST EXECUTABLE STATEMENT  DQNC79
  if (FIRST) THEN
    W1 = 41.0D0/140.0D0
    W2 = 216.0D0/140.0D0
    W3 = 27.0D0/140.0D0
    W4 = 272.0D0/140.0D0
    NBITS = D1MACH(5)*I1MACH(14)/0.30102000D0
    NLMX = MIN(99,(NBITS*4)/5)
    SQ2 = SQRT(2.0D0)
  end if
  FIRST = .FALSE.
  ANS = 0.0D0
  IERR = 1
  CE = 0.0D0
  if (A  ==  B) go to 260
  LMX = NLMX
  LMN = NLMN
  if (B  ==  0.0D0) go to 100
  if (SIGN(1.0D0,B)*A  <=  0.0D0) go to 100
  C = ABS(1.0D0-A/B)
  if (C  >  0.1D0) go to 100
  if (C  <=  0.0D0) go to 260
  NIB = 0.5D0 - LOG(C)/LOG(2.0D0)
  LMX = MIN(NLMX,NBITS-NIB-4)
  if (LMX  <  2) go to 260
  LMN = MIN(LMN,LMX)
  100 TOL = MAX(ABS(ERR),2.0D0**(5-NBITS))
  if (ERR  ==  0.0D0) TOL = SQRT(D1MACH(4))
  EPS = TOL
  HH(1) = (B-A)/12.0D0
  AA(1) = A
  LR(1) = 1
  DO 110 I = 1,11,2
    F(I) = FUN(A+(I-1)*HH(1))
  110 CONTINUE
  BLOCAL = B
  F(13) = FUN(BLOCAL)
  K = 7
  L = 1
  AREA = 0.0D0
  Q7 = 0.0D0
  EF = 256.0D0/255.0D0
  BANK = 0.0D0
!
!     Compute refined estimates, estimate the error, etc.
!
  120 DO 130 I = 2,12,2
    F(I) = FUN(AA(L)+(I-1)*HH(L))
  130 CONTINUE
  K = K + 6
!
!     Compute left and right half estimates
!
  Q7L = HH(L)*((W1*(F(1)+F(7))+W2*(F(2)+F(6)))+ &
        (W3*(F(3)+F(5))+W4*F(4)))
  Q7R(L) = HH(L)*((W1*(F(7)+F(13))+W2*(F(8)+F(12)))+ &
           (W3*(F(9)+F(11))+W4*F(10)))
!
!     Update estimate of integral of absolute value
!
  AREA = AREA + (ABS(Q7L)+ABS(Q7R(L))-ABS(Q7))
!
!     Do not bother to test convergence before minimum refinement level
!
  if (L  <  LMN) go to 180
!
!     Estimate the error in new value for whole interval, Q13
!
  Q13 = Q7L + Q7R(L)
  EE = ABS(Q7-Q13)*EF
!
!     Compute nominal allowed error
!
  AE = EPS*AREA
!
!     Borrow from bank account, but not too much
!
  TEST = MIN(AE+0.8D0*BANK,10.0D0*AE)
!
!     Don't ask for excessive accuracy
!
  TEST = MAX(TEST,TOL*ABS(Q13),0.00003D0*TOL*AREA)
!
!     Now, did this interval pass or not?
!
  if (EE-TEST) 150,150,170
!
!     Have hit maximum refinement level -- penalize the cumulative error
!
  140 CE = CE + (Q7-Q13)
  go to 160
!
!     On good intervals accumulate the theoretical estimate
!
  150 CE = CE + (Q7-Q13)/255.0D0
!
!     Update the bank account.  Don't go into debt.
!
  160 BANK = BANK + (AE-EE)
  if (BANK  <  0.0D0) BANK = 0.0D0
!
!     Did we just finish a left half or a right half?
!
  if (LR(L)) 190,190,210
!
!     Consider the left half of next deeper level
!
  170 if (K  >  KMX) LMX = MIN(KML,LMX)
  if (L  >=  LMX) go to 140
  180 L = L + 1
  EPS = EPS*0.5D0
  if (L  <=  17) EF = EF/SQ2
  HH(L) = HH(L-1)*0.5D0
  LR(L) = -1
  AA(L) = AA(L-1)
  Q7 = Q7L
  F1(L) = F(7)
  F2(L) = F(8)
  F3(L) = F(9)
  F4(L) = F(10)
  F5(L) = F(11)
  F6(L) = F(12)
  F7(L) = F(13)
  F(13) = F(7)
  F(11) = F(6)
  F(9) = F(5)
  F(7) = F(4)
  F(5) = F(3)
  F(3) = F(2)
  go to 120
!
!     Proceed to right half at this level
!
  190 VL(L) = Q13
  200 Q7 = Q7R(L-1)
  LR(L) = 1
  AA(L) = AA(L) + 12.0D0*HH(L)
  F(1) = F1(L)
  F(3) = F2(L)
  F(5) = F3(L)
  F(7) = F4(L)
  F(9) = F5(L)
  F(11) = F6(L)
  F(13) = F7(L)
  go to 120
!
!     Left and right halves are done, so go back up a level
!
  210 VR = Q13
  220 if (L  <=  1) go to 250
  if (L  <=  17) EF = EF*SQ2
  EPS = EPS*2.0D0
  L = L - 1
  if (LR(L)) 230,230,240
  230 VL(L) = VL(L+1) + VR
  go to 200
  240 VR = VL(L+1) + VR
  go to 220
!
!     Exit
!
  250 ANS = VR
  if (ABS(CE)  <=  2.0D0*TOL*AREA) go to 270
  IERR = 2
  call XERMSG ('SLATEC', 'DQNC79', &
     'ANS is probably insufficiently accurate.', 2, 1)
  go to 270
  260 IERR = -1
  call XERMSG ('SLATEC', 'DQNC79', &
     'A and B are too nearly equal to allow normal integration. $$' &
     // 'ANS is set to zero and IERR to -1.', -1, -1)
  270 RETURN
end
