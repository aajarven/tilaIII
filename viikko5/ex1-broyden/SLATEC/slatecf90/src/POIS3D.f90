subroutine POIS3D (LPEROD, L, C1, MPEROD, M, C2, NPEROD, N, A, B, &
     C, LDIMF, MDIMF, F, IERROR, W)
!
!! POIS3D solves a three-dimensional block tridiagonal linear system ...
!            which arises from a finite difference approximation to a
!            three-dimensional Poisson equation using the Fourier
!            transform package FFTPAK written by Paul Swarztrauber.
!
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B4B
!***TYPE      SINGLE PRECISION (POIS3D-S)
!***KEYWORDS  ELLIPTIC PDE, FISHPACK, HELMHOLTZ, POISSON
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine POIS3D solves the linear system of equations
!
!       C1*(X(I-1,J,K)-2.*X(I,J,K)+X(I+1,J,K))
!     + C2*(X(I,J-1,K)-2.*X(I,J,K)+X(I,J+1,K))
!     + A(K)*X(I,J,K-1)+B(K)*X(I,J,K)+C(K)*X(I,J,K+1) = F(I,J,K)
!
!     for  I=1,2,...,L , J=1,2,...,M , and K=1,2,...,N .
!
!     The indices K-1 and K+1 are evaluated modulo N, i.e.
!     X(I,J,0) = X(I,J,N) and X(I,J,N+1) = X(I,J,1). The unknowns
!     X(0,J,K), X(L+1,J,K), X(I,0,K), and X(I,M+1,K) are assumed to take
!     on certain prescribed values described below.
!
!    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!    * * * * * * * *    Parameter Description     * * * * * * * * * *
!
!
!            * * * * * *   On Input    * * * * * *
!
!     LPEROD   Indicates the values that X(0,J,K) and X(L+1,J,K) are
!              assumed to have.
!
!              = 0  If X(0,J,K) = X(L,J,K) and X(L+1,J,K) = X(1,J,K).
!              = 1  If X(0,J,K) = X(L+1,J,K) = 0.
!              = 2  If X(0,J,K) = 0  and X(L+1,J,K) = X(L-1,J,K).
!              = 3  If X(0,J,K) = X(2,J,K) and X(L+1,J,K) = X(L-1,J,K).
!              = 4  If X(0,J,K) = X(2,J,K) and X(L+1,J,K) = 0.
!
!     L        The number of unknowns in the I-direction. L must be at
!              least 3.
!
!     C1       The real constant that appears in the above equation.
!
!     MPEROD   Indicates the values that X(I,0,K) and X(I,M+1,K) are
!              assumed to have.
!
!              = 0  If X(I,0,K) = X(I,M,K) and X(I,M+1,K) = X(I,1,K).
!              = 1  If X(I,0,K) = X(I,M+1,K) = 0.
!              = 2  If X(I,0,K) = 0 and X(I,M+1,K) = X(I,M-1,K).
!              = 3  If X(I,0,K) = X(I,2,K) and X(I,M+1,K) = X(I,M-1,K).
!              = 4  If X(I,0,K) = X(I,2,K) and X(I,M+1,K) = 0.
!
!     M        The number of unknowns in the J-direction. M must be at
!              least 3.
!
!     C2       The real constant which appears in the above equation.
!
!     NPEROD   = 0  If A(1) and C(N) are not zero.
!              = 1  If A(1) = C(N) = 0.
!
!     N        The number of unknowns in the K-direction. N must be at
!              least 3.
!
!
!     A,B,C    One-dimensional arrays of length N that specify the
!              coefficients in the linear equations given above.
!
!              If NPEROD = 0 the array elements must not depend upon the
!              index K, but must be constant.  Specifically, the
!              subroutine checks the following condition
!
!                          A(K) = C(1)
!                          C(K) = C(1)
!                          B(K) = B(1)
!
!                  for K=1,2,...,N.
!
!     LDIMF    The row (or first) dimension of the three-dimensional
!              array F as it appears in the program calling POIS3D.
!              This parameter is used to specify the variable dimension
!              of F.  LDIMF must be at least L.
!
!     MDIMF    The column (or second) dimension of the three-dimensional
!              array F as it appears in the program calling POIS3D.
!              This parameter is used to specify the variable dimension
!              of F.  MDIMF must be at least M.
!
!     F        A three-dimensional array that specifies the values of
!              the right side of the linear system of equations given
!              above.  F must be dimensioned at least L x M x N.
!
!     W        A one-dimensional array that must be provided by the
!              user for work space.  The length of W must be at least
!              30 + L + M + 2*N + MAX(L,M,N) +
!              7*(INT((L+1)/2) + INT((M+1)/2)).
!
!
!            * * * * * *   On Output   * * * * * *
!
!     F        Contains the solution X.
!
!     IERROR   An error flag that indicates invalid input parameters.
!              Except for number zero, a solution is not attempted.
!              = 0  No error
!              = 1  If LPEROD  <  0 or  >  4
!              = 2  If L  <  3
!              = 3  If MPEROD  <  0 or  >  4
!              = 4  If M  <  3
!              = 5  If NPEROD  <  0 or  >  1
!              = 6  If N  <  3
!              = 7  If LDIMF  <  L
!              = 8  If MDIMF  <  M
!              = 9  If A(K)  /=  C(1) or C(K)  /=  C(1) or B(I)  /= B(1)
!                      for some K=1,2,...,N.
!              = 10 If NPEROD = 1 and A(1)  /=  0 or C(N)  /=  0
!
!              Since this is the only means of indicating a possibly
!              incorrect call to POIS3D, the user should test IERROR
!              after the call.
!
! *Long Description:
!
!    * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension of   A(N),B(N),C(N),F(LDIMF,MDIMF,N),
!     Arguments      W(see argument list)
!
!     Latest         December 1, 1978
!     Revision
!
!     Subprograms    POIS3D,POS3D1,TRIDQ,RFFTI,RFFTF,RFFTF1,RFFTB,
!     Required       RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,COSQF1
!                    COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,CFFTI1,
!                    CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,CFFTF,
!                    CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,PIMACH,
!
!     Special        NONE
!     Conditions
!
!     Common         NONE
!     Blocks
!
!     I/O            NONE
!
!     Precision      Single
!
!     Specialist     Roland Sweet
!
!     Language       FORTRAN
!
!     History        Written by Roland Sweet at NCAR in July 1977
!
!     Algorithm      This subroutine solves three-dimensional block
!                    tridiagonal linear systems arising from finite
!                    difference approximations to three-dimensional
!                    Poisson equations using the Fourier transform
!                    package FFTPAK written by Paul Swarztrauber.
!
!     Space          6561(decimal) = 14641(octal) locations on the
!     Required       NCAR Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine POIS3D is roughly proportional
!                    to L*M*N*(log2(L)+log2(M)+5), but also depends on
!                    input parameters LPEROD and MPEROD.  Some typical
!                    values are listed in the table below when NPEROD=0.
!                       To measure the accuracy of the algorithm a
!                    uniform random number generator was used to create
!                    a solution array X for the system given in the
!                    'PURPOSE' with
!
!                       A(K) = C(K) = -0.5*B(K) = 1,       K=1,2,...,N
!
!                    and, when NPEROD = 1
!
!                       A(1) = C(N) = 0
!                       A(N) = C(1) = 2.
!
!                    The solution X was substituted into the given sys-
!                    tem and, using double precision, a right side Y was
!                    computed.  Using this array Y subroutine POIS3D was
!                    called to produce an approximate solution Z.  Then
!                    the relative error, defined as
!
!                    E = MAX(ABS(Z(I,J,K)-X(I,J,K)))/MAX(ABS(X(I,J,K)))
!
!                    where the two maxima are taken over I=1,2,...,L,
!                    J=1,2,...,M and K=1,2,...,N, was computed.  The
!                    value of E is given in the table below for some
!                    typical values of L,M and N.
!
!
!                       L(=M=N)   LPEROD    MPEROD    T(MSECS)    E
!                       ------    ------    ------    --------  ------
!
!                         16        0         0         272     1.E-13
!                         15        1         1         287     4.E-13
!                         17        3         3         338     2.E-13
!                         32        0         0        1755     2.E-13
!                         31        1         1        1894     2.E-12
!                         33        3         3        2042     7.E-13
!
!
!     Portability    American National Standards Institute FORTRAN.
!                    The machine dependent constant PI is defined in
!                    function PIMACH.
!
!     Required       COS,SIN,ATAN
!     Resident
!     Routines
!
!     Reference      NONE
!
!    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  POS3D1
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  POIS3D
  DIMENSION       A(*)       ,B(*)       ,C(*)       , &
                  F(LDIMF,MDIMF,*)       ,W(*)       ,SAVE(6)
!***FIRST EXECUTABLE STATEMENT  POIS3D
  LP = LPEROD+1
  MP = MPEROD+1
  NP = NPEROD+1
!
!     CHECK FOR INVALID INPUT.
!
  IERROR = 0
  if (LP < 1 .OR. LP > 5) IERROR = 1
  if (L  <  3) IERROR = 2
  if (MP < 1 .OR. MP > 5) IERROR = 3
  if (M  <  3) IERROR = 4
  if (NP < 1 .OR. NP > 2) IERROR = 5
  if (N  <  3) IERROR = 6
  if (LDIMF  <  L) IERROR = 7
  if (MDIMF  <  M) IERROR = 8
  if (NP  /=  1) go to 103
  DO 101 K=1,N
     if (A(K)  /=  C(1)) go to 102
     if (C(K)  /=  C(1)) go to 102
     if (B(K)  /=  B(1)) go to 102
  101 CONTINUE
  go to 104
  102 IERROR = 9
  103 if (NPEROD == 1 .AND. (A(1) /= 0. .OR. C(N) /= 0.)) IERROR = 10
  104 if (IERROR  /=  0) go to 122
  IWYRT = L+1
  IWT = IWYRT+M
  IWD = IWT+MAX(L,M,N)+1
  IWBB = IWD+N
  IWX = IWBB+N
  IWY = IWX+7*((L+1)/2)+15
  go to (105,114),NP
!
!     REORDER UNKNOWNS WHEN NPEROD = 0.
!
  105 NH = (N+1)/2
  NHM1 = NH-1
  NODD = 1
  if (2*NH  ==  N) NODD = 2
  DO 111 I=1,L
     DO 110 J=1,M
        DO 106 K=1,NHM1
           NHPK = NH+K
           NHMK = NH-K
           W(K) = F(I,J,NHMK)-F(I,J,NHPK)
           W(NHPK) = F(I,J,NHMK)+F(I,J,NHPK)
  106       CONTINUE
        W(NH) = 2.*F(I,J,NH)
        go to (108,107),NODD
  107       W(N) = 2.*F(I,J,N)
  108       DO 109 K=1,N
           F(I,J,K) = W(K)
  109       CONTINUE
  110    CONTINUE
  111 CONTINUE
  SAVE(1) = C(NHM1)
  SAVE(2) = A(NH)
  SAVE(3) = C(NH)
  SAVE(4) = B(NHM1)
  SAVE(5) = B(N)
  SAVE(6) = A(N)
  C(NHM1) = 0.
  A(NH) = 0.
  C(NH) = 2.*C(NH)
  go to (112,113),NODD
  112 B(NHM1) = B(NHM1)-A(NH-1)
  B(N) = B(N)+A(N)
  go to 114
  113 A(N) = C(NH)
  114 CONTINUE
  call POS3D1 (LP,L,MP,M,N,A,B,C,LDIMF,MDIMF,F,W,W(IWYRT),W(IWT), &
               W(IWD),W(IWX),W(IWY),C1,C2,W(IWBB))
  go to (115,122),NP
  115 DO 121 I=1,L
     DO 120 J=1,M
        DO 116 K=1,NHM1
           NHMK = NH-K
           NHPK = NH+K
           W(NHMK) = .5*(F(I,J,NHPK)+F(I,J,K))
           W(NHPK) = .5*(F(I,J,NHPK)-F(I,J,K))
  116       CONTINUE
        W(NH) = .5*F(I,J,NH)
        go to (118,117),NODD
  117       W(N) = .5*F(I,J,N)
  118       DO 119 K=1,N
           F(I,J,K) = W(K)
  119       CONTINUE
  120    CONTINUE
  121 CONTINUE
  C(NHM1) = SAVE(1)
  A(NH) = SAVE(2)
  C(NH) = SAVE(3)
  B(NHM1) = SAVE(4)
  B(N) = SAVE(5)
  A(N) = SAVE(6)
  122 CONTINUE
  return
end
