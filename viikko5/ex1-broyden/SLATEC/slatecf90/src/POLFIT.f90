subroutine POLFIT (N, X, Y, W, MAXDEG, NDEG, EPS, R, IERR, A)
!
!! POLFIT fits discrete data in a least squares sense by polynomials ...
!            in one variable.
!
!***LIBRARY   SLATEC
!***CATEGORY  K1A1A2
!***TYPE      SINGLE PRECISION (POLFIT-S, DPOLFT-D)
!***KEYWORDS  CURVE FITTING, DATA FITTING, LEAST SQUARES, POLYNOMIAL FIT
!***AUTHOR  Shampine, L. F., (SNLA)
!           Davenport, S. M., (SNLA)
!           Huddleston, R. E., (SNLL)
!***DESCRIPTION
!
!     Abstract
!
!     Given a collection of points X(I) and a set of values Y(I) which
!     correspond to some function or measurement at each of the X(I),
!     subroutine  POLFIT  computes the weighted least-squares polynomial
!     fits of all degrees up to some degree either specified by the user
!     or determined by the routine.  The fits thus obtained are in
!     orthogonal polynomial form.  Subroutine  PVALUE  may then be
!     called to evaluate the fitted polynomials and any of their
!     derivatives at any point.  The subroutine  PCOEF  may be used to
!     express the polynomial fits as powers of (X-C) for any specified
!     point C.
!
!     The parameters for  POLFIT  are
!
!     Input --
!         N -      the number of data points.  The arrays X, Y and W
!                  must be dimensioned at least  N  (N  >=  1).
!         X -      array of values of the independent variable.  These
!                  values may appear in any order and need not all be
!                  distinct.
!         Y -      array of corresponding function values.
!         W -      array of positive values to be used as weights.  If
!                  W(1) is negative,  POLFIT  will set all the weights
!                  to 1.0, which means unweighted least squares error
!                  will be minimized.  To minimize relative error, the
!                  user should set the weights to:  W(I) = 1.0/Y(I)**2,
!                  I = 1,...,N .
!         MAXDEG - maximum degree to be allowed for polynomial fit.
!                  MAXDEG  may be any non-negative integer less than  N.
!                  Note -- MAXDEG  cannot be equal to  N-1  when a
!                  statistical test is to be used for degree selection,
!                  i.e., when input value of  EPS  is negative.
!         EPS -    specifies the criterion to be used in determining
!                  the degree of fit to be computed.
!                  (1)  If  EPS  is input negative,  POLFIT  chooses the
!                       degree based on a statistical F test of
!                       significance.  One of three possible
!                       significance levels will be used:  .01, .05 or
!                       .10.  If  EPS=-1.0 , the routine will
!                       automatically select one of these levels based
!                       on the number of data points and the maximum
!                       degree to be considered.  If  EPS  is input as
!                       -.01, -.05, or -.10, a significance level of
!                       .01, .05, or .10, respectively, will be used.
!                  (2)  If  EPS  is set to 0.,  POLFIT  computes the
!                       polynomials of degrees 0 through  MAXDEG .
!                  (3)  If  EPS  is input positive,  EPS  is the RMS
!                       error tolerance which must be satisfied by the
!                       fitted polynomial.  POLFIT  will increase the
!                       degree of fit until this criterion is met or
!                       until the maximum degree is reached.
!
!     Output --
!         NDEG -   degree of the highest degree fit computed.
!         EPS -    RMS error of the polynomial of degree  NDEG .
!         R -      vector of dimension at least NDEG containing values
!                  of the fit of degree  NDEG  at each of the  X(I) .
!                  Except when the statistical test is used, these
!                  values are more accurate than results from subroutine
!                  PVALUE  normally are.
!         IERR -   error flag with the following possible values.
!             1 -- indicates normal execution, i.e., either
!                  (1)  the input value of  EPS  was negative, and the
!                       computed polynomial fit of degree  NDEG
!                       satisfies the specified F test, or
!                  (2)  the input value of  EPS  was 0., and the fits of
!                       all degrees up to  MAXDEG  are complete, or
!                  (3)  the input value of  EPS  was positive, and the
!                       polynomial of degree  NDEG  satisfies the RMS
!                       error requirement.
!             2 -- invalid input parameter.  At least one of the input
!                  parameters has an illegal value and must be corrected
!                  before  POLFIT  can proceed.  Valid input results
!                  when the following restrictions are observed
!                       N  >=  1
!                       0  <=  MAXDEG  <=  N-1  for  EPS  >=  0.
!                       0  <=  MAXDEG  <=  N-2  for  EPS  <  0.
!                       W(1)=-1.0  or  W(I)  >  0., I=1,...,N .
!             3 -- cannot satisfy the RMS error requirement with a
!                  polynomial of degree no greater than  MAXDEG .  Best
!                  fit found is of degree  MAXDEG .
!             4 -- cannot satisfy the test for significance using
!                  current value of  MAXDEG .  Statistically, the
!                  best fit found is of order  NORD .  (In this case,
!                  NDEG will have one of the values:  MAXDEG-2,
!                  MAXDEG-1, or MAXDEG).  Using a higher value of
!                  MAXDEG  may result in passing the test.
!         A -      work and output array having at least 3N+3MAXDEG+3
!                  locations
!
!     Note - POLFIT  calculates all fits of degrees up to and including
!            NDEG .  Any or all of these fits can be evaluated or
!            expressed as powers of (X-C) using  PVALUE  and  PCOEF
!            after just one call to  POLFIT .
!
!***REFERENCES  L. F. Shampine, S. M. Davenport and R. E. Huddleston,
!                 Curve fitting by polynomials in one variable, Report
!                 SLA-74-0270, Sandia Laboratories, June 1974.
!***ROUTINES CALLED  PVALUE, XERMSG
!***REVISION HISTORY  (YYMMDD)
!   740601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
!   920501  Reformatted the REFERENCES section.  (WRB)
!   920527  Corrected erroneous statements in DESCRIPTION.  (WRB)
!***END PROLOGUE  POLFIT
  DOUBLE PRECISION TEMD1,TEMD2
  DIMENSION X(*), Y(*), W(*), R(*), A(*)
  DIMENSION CO(4,3)
  SAVE CO
  DATA  CO(1,1), CO(2,1), CO(3,1), CO(4,1), CO(1,2), CO(2,2), &
        CO(3,2), CO(4,2), CO(1,3), CO(2,3), CO(3,3), &
    CO(4,3)/-13.086850,-2.4648165,-3.3846535,-1.2973162, &
            -3.3381146,-1.7812271,-3.2578406,-1.6589279, &
            -1.6282703,-1.3152745,-3.2640179,-1.9829776/
!***FIRST EXECUTABLE STATEMENT  POLFIT
  M = ABS(N)
  if (M  ==  0) go to 30
  if (MAXDEG  <  0) go to 30
  A(1) = MAXDEG
  MOP1 = MAXDEG + 1
  if (M  <  MOP1) go to 30
  if (EPS  <  0.0  .AND.  M  ==  MOP1) go to 30
  XM = M
  ETST = EPS*EPS*XM
  if (W(1)  <  0.0) go to 2
  DO 1 I = 1,M
    if (W(I)  <=  0.0) go to 30
 1      CONTINUE
  go to 4
 2    DO 3 I = 1,M
 3      W(I) = 1.0
 4    if (EPS  >=  0.0) go to 8
!
! DETERMINE SIGNIFICANCE LEVEL INDEX TO BE USED IN STATISTICAL TEST FOR
! CHOOSING DEGREE OF POLYNOMIAL FIT
!
  if (EPS  >  (-.55)) go to 5
  IDEGF = M - MAXDEG - 1
  KSIG = 1
  if (IDEGF  <  10) KSIG = 2
  if (IDEGF  <  5) KSIG = 3
  go to 8
 5    KSIG = 1
  if (EPS  <  (-.03)) KSIG = 2
  if (EPS  <  (-.07)) KSIG = 3
!
! INITIALIZE INDEXES AND COEFFICIENTS FOR FITTING
!
 8    K1 = MAXDEG + 1
  K2 = K1 + MAXDEG
  K3 = K2 + MAXDEG + 2
  K4 = K3 + M
  K5 = K4 + M
  DO 9 I = 2,K4
 9      A(I) = 0.0
  W11 = 0.0
  if (N  <  0) go to 11
!
! UNCONSTRAINED CASE
!
  DO 10 I = 1,M
    K4PI = K4 + I
    A(K4PI) = 1.0
 10     W11 = W11 + W(I)
  go to 13
!
! CONSTRAINED CASE
!
 11   DO 12 I = 1,M
    K4PI = K4 + I
 12     W11 = W11 + W(I)*A(K4PI)**2
!
! COMPUTE FIT OF DEGREE ZERO
!
 13   TEMD1 = 0.0D0
  DO 14 I = 1,M
    K4PI = K4 + I
    TEMD1 = TEMD1 + DBLE(W(I))*DBLE(Y(I))*DBLE(A(K4PI))
 14     CONTINUE
  TEMD1 = TEMD1/DBLE(W11)
  A(K2+1) = TEMD1
  SIGJ = 0.0
  DO 15 I = 1,M
    K4PI = K4 + I
    K5PI = K5 + I
    TEMD2 = TEMD1*DBLE(A(K4PI))
    R(I) = TEMD2
    A(K5PI) = TEMD2 - DBLE(R(I))
 15     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
  J = 0
!
! SEE if POLYNOMIAL OF DEGREE 0 SATISFIES THE DEGREE SELECTION CRITERION
!
  if (EPS) 24,26,27
!
! INCREMENT DEGREE
!
 16   J = J + 1
  JP1 = J + 1
  K1PJ = K1 + J
  K2PJ = K2 + J
  SIGJM1 = SIGJ
!
! COMPUTE NEW B COEFFICIENT EXCEPT WHEN J = 1
!
  if (J  >  1) A(K1PJ) = W11/W1
!
! COMPUTE NEW A COEFFICIENT
!
  TEMD1 = 0.0D0
  DO 18 I = 1,M
    K4PI = K4 + I
    TEMD2 = A(K4PI)
    TEMD1 = TEMD1 + DBLE(X(I))*DBLE(W(I))*TEMD2*TEMD2
 18     CONTINUE
  A(JP1) = TEMD1/DBLE(W11)
!
! EVALUATE ORTHOGONAL POLYNOMIAL AT DATA POINTS
!
  W1 = W11
  W11 = 0.0
  DO 19 I = 1,M
    K3PI = K3 + I
    K4PI = K4 + I
    TEMP = A(K3PI)
    A(K3PI) = A(K4PI)
    A(K4PI) = (X(I)-A(JP1))*A(K3PI) - A(K1PJ)*TEMP
 19     W11 = W11 + W(I)*A(K4PI)**2
!
! GET NEW ORTHOGONAL POLYNOMIAL COEFFICIENT USING PARTIAL DOUBLE
! PRECISION
!
  TEMD1 = 0.0D0
  DO 20 I = 1,M
    K4PI = K4 + I
    K5PI = K5 + I
    TEMD2 = DBLE(W(I))*DBLE((Y(I)-R(I))-A(K5PI))*DBLE(A(K4PI))
 20     TEMD1 = TEMD1 + TEMD2
  TEMD1 = TEMD1/DBLE(W11)
  A(K2PJ+1) = TEMD1
!
! UPDATE POLYNOMIAL EVALUATIONS AT EACH OF THE DATA POINTS, AND
! ACCUMULATE SUM OF SQUARES OF ERRORS.  THE POLYNOMIAL EVALUATIONS ARE
! COMPUTED AND STORED IN EXTENDED PRECISION.  FOR THE I-TH DATA POINT,
! THE MOST SIGNIFICANT BITS ARE STORED IN  R(I) , AND THE LEAST
! SIGNIFICANT BITS ARE IN  A(K5PI) .
!
  SIGJ = 0.0
  DO 21 I = 1,M
    K4PI = K4 + I
    K5PI = K5 + I
    TEMD2 = DBLE(R(I)) + DBLE(A(K5PI)) + TEMD1*DBLE(A(K4PI))
    R(I) = TEMD2
    A(K5PI) = TEMD2 - DBLE(R(I))
 21     SIGJ = SIGJ + W(I)*((Y(I)-R(I)) - A(K5PI))**2
!
! SEE if DEGREE SELECTION CRITERION HAS BEEN SATISFIED OR IF DEGREE
! MAXDEG  HAS BEEN REACHED
!
  if (EPS) 23,26,27
!
! COMPUTE F STATISTICS  (INPUT EPS  <  0.)
!
 23   if (SIGJ  ==  0.0) go to 29
  DEGF = M - J - 1
  DEN = (CO(4,KSIG)*DEGF + 1.0)*DEGF
  FCRIT = (((CO(3,KSIG)*DEGF) + CO(2,KSIG))*DEGF + CO(1,KSIG))/DEN
  FCRIT = FCRIT*FCRIT
  F = (SIGJM1 - SIGJ)*DEGF/SIGJ
  if (F  <  FCRIT) go to 25
!
! POLYNOMIAL OF DEGREE J SATISFIES F TEST
!
 24   SIGPAS = SIGJ
  JPAS = J
  NFAIL = 0
  if (MAXDEG  ==  J) go to 32
  go to 16
!
! POLYNOMIAL OF DEGREE J FAILS F TEST.  if THERE HAVE BEEN THREE
! SUCCESSIVE FAILURES, A STATISTICALLY BEST DEGREE HAS BEEN FOUND.
!
 25   NFAIL = NFAIL + 1
  if (NFAIL  >=  3) go to 29
  if (MAXDEG  ==  J) go to 32
  go to 16
!
! RAISE THE DEGREE if DEGREE  MAXDEG  HAS NOT YET BEEN REACHED  (INPUT
! EPS = 0.)
!
 26   if (MAXDEG  ==  J) go to 28
  go to 16
!
! SEE if RMS ERROR CRITERION IS SATISFIED  (INPUT EPS  >  0.)
!
 27   if (SIGJ  <=  ETST) go to 28
  if (MAXDEG  ==  J) go to 31
  go to 16
!
! RETURNS
!
 28   IERR = 1
  NDEG = J
  SIG = SIGJ
  go to 33
 29   IERR = 1
  NDEG = JPAS
  SIG = SIGPAS
  go to 33
 30   IERR = 2
  call XERMSG ('SLATEC', 'POLFIT', 'INVALID INPUT PARAMETER.', 2, &
     1)
  go to 37
 31   IERR = 3
  NDEG = MAXDEG
  SIG = SIGJ
  go to 33
 32   IERR = 4
  NDEG = JPAS
  SIG = SIGPAS
!
 33   A(K3) = NDEG
!
! WHEN STATISTICAL TEST HAS BEEN USED, EVALUATE THE BEST POLYNOMIAL AT
! ALL THE DATA POINTS if  R  DOES NOT ALREADY CONTAIN THESE VALUES
!
  if ( EPS  >=  0.0  .OR.  NDEG  ==  MAXDEG) go to 36
  NDER = 0
  DO 35 I = 1,M
    call PVALUE (NDEG,NDER,X(I),R(I),YP,A)
 35     CONTINUE
 36   EPS = SQRT(SIG/XM)
 37   return
end
