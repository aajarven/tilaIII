subroutine HWSCSP (INTL, TS, TF, M, MBDCND, BDTS, BDTF, RS, RF, N, &
     NBDCND, BDRS, BDRF, ELMBDA, F, IDIMF, PERTRB, IERROR, W)
!
!! HWSCSP solves a finite difference approximation to the modified ...
!            Helmholtz equation in spherical coordinates assuming
!            axisymmetry  (no dependence on longitude).
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B1A1A
!***TYPE      SINGLE PRECISION (HWSCSP-S)
!***KEYWORDS  ELLIPTIC, FISHPACK, HELMHOLTZ, PDE, SPHERICAL
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine HWSCSP solves a finite difference approximation to the
!       modified Helmholtz equation in spherical coordinates assuming
!       axisymmetry  (no dependence on longitude)
!
!          (1/R**2)(d/dR)((R**2)(d/dR)U)
!
!             + (1/(R**2)SIN(THETA))(d/dTHETA)(SIN(THETA)(d/dTHETA)U)
!
!             + (LAMBDA/(RSIN(THETA))**2)U = F(THETA,R).
!
!     This two dimensional modified Helmholtz equation results from
!     the Fourier transform of the three dimensional Poisson equation
!
!     * * * * * * * * * *     On Input     * * * * * * * * * *
!
!     INTL
!       = 0  On initial entry to HWSCSP or if any of the arguments
!            RS, RF, N, NBDCND are changed from a previous call.
!       = 1  If RS, RF, N, NBDCND are all unchanged from previous call
!            to HWSCSP.
!
!       NOTE   A call with INTL=0 takes approximately 1.5 times as
!              much time as a call with INTL = 1.  Once a call with
!              INTL = 0 has been made then subsequent solutions
!              corresponding to different F, BDTS, BDTF, BDRS, BDRF can
!              be obtained faster with INTL = 1 since initialization is
!              not repeated.
!
!     TS,TF
!       The range of THETA (colatitude), i.e., TS  <=  THETA  <=  TF.
!       TS must be less than TF.  TS and TF are in radians.  A TS of
!       zero corresponds to the north pole and a TF of PI corresponds
!       to the south pole.
!
!     * * * * * * * * * * * * * * IMPORTANT * * * * * * * * * * * * * *
!
!     If TF is equal to PI then it must be computed using the statement
!     TF = PIMACH(DUM). This insures that TF in the users program is
!     equal to PI in this program which permits several tests of the
!     input parameters that otherwise would not be possible.
!
!     M
!       The number of panels into which the interval (TS,TF) is
!       subdivided.  Hence, there will be M+1 grid points in the
!       THETA-direction given by THETA(K) = (I-1)DTHETA+TS for
!       I = 1,2,...,M+1, where DTHETA = (TF-TS)/M is the panel width.
!
!     MBDCND
!       Indicates the type of boundary condition at THETA = TS and
!       THETA = TF.
!
!       = 1  If the solution is specified at THETA = TS and THETA = TF.
!       = 2  If the solution is specified at THETA = TS and the
!            derivative of the solution with respect to THETA is
!            specified at THETA = TF (see note 2 below).
!       = 3  If the derivative of the solution with respect to THETA is
!            specified at THETA = TS and THETA = TF (see notes 1,2
!            below).
!       = 4  If the derivative of the solution with respect to THETA is
!            specified at THETA = TS (see note 1 below) and the
!            solution is specified at THETA = TF.
!       = 5  If the solution is unspecified at THETA = TS = 0 and the
!            solution is specified at THETA = TF.
!       = 6  If the solution is unspecified at THETA = TS = 0 and the
!            derivative of the solution with respect to THETA is
!            specified at THETA = TF (see note 2 below).
!       = 7  If the solution is specified at THETA = TS and the
!            solution is unspecified at THETA = TF = PI.
!       = 8  If the derivative of the solution with respect to THETA is
!            specified at THETA = TS (see note 1 below) and the solution
!            is unspecified at THETA = TF = PI.
!       = 9  If the solution is unspecified at THETA = TS = 0 and
!            THETA = TF = PI.
!
!       NOTES:  1.  If TS = 0, do not use MBDCND = 3,4, or 8, but
!                   instead use MBDCND = 5,6, or 9  .
!               2.  If TF = PI, do not use MBDCND = 2,3, or 6, but
!                   instead use MBDCND = 7,8, or 9  .
!
!     BDTS
!       A one-dimensional array of length N+1 that specifies the values
!       of the derivative of the solution with respect to THETA at
!       THETA = TS.  When MBDCND = 3,4, or 8,
!
!            BDTS(J) = (d/dTHETA)U(TS,R(J)), J = 1,2,...,N+1  .
!
!       When MBDCND has any other value, BDTS is a dummy variable.
!
!     BDTF
!       A one-dimensional array of length N+1 that specifies the values
!       of the derivative of the solution with respect to THETA at
!       THETA = TF.  When MBDCND = 2,3, or 6,
!
!            BDTF(J) = (d/dTHETA)U(TF,R(J)), J = 1,2,...,N+1  .
!
!       When MBDCND has any other value, BDTF is a dummy variable.
!
!     RS,RF
!       The range of R, i.e., RS  <=  R  <  RF.  RS must be less than
!       RF.  RS must be non-negative.
!
!       N
!       The number of panels into which the interval (RS,RF) is
!       subdivided.  Hence, there will be N+1 grid points in the
!       R-direction given by R(J) = (J-1)DR+RS for J = 1,2,...,N+1,
!       where DR = (RF-RS)/N is the panel width.
!       N must be greater than 2
!
!     NBDCND
!       Indicates the type of boundary condition at R = RS and R = RF.
!
!       = 1  If the solution is specified at R = RS and R = RF.
!       = 2  If the solution is specified at R = RS and the derivative
!            of the solution with respect to R is specified at R = RF.
!       = 3  If the derivative of the solution with respect to R is
!            specified at R = RS and R = RF.
!       = 4  If the derivative of the solution with respect to R is
!            specified at RS and the solution is specified at R = RF.
!       = 5  If the solution is unspecified at R = RS = 0 (see note
!            below) and the solution is specified at R = RF.
!       = 6  If the solution is unspecified at R = RS = 0 (see note
!            below) and the derivative of the solution with respect to
!            R is specified at R = RF.
!
!       NOTE:  NBDCND = 5 or 6 cannot be used with
!              MBDCND = 1,2,4,5, or 7 (the former indicates that the
!                       solution is unspecified at R = 0, the latter
!                       indicates that the solution is specified).
!                       Use instead
!              NBDCND = 1 or 2  .
!
!     BDRS
!       A one-dimensional array of length M+1 that specifies the values
!       of the derivative of the solution with respect to R at R = RS.
!       When NBDCND = 3 or 4,
!
!            BDRS(I) = (d/dR)U(THETA(I),RS), I = 1,2,...,M+1  .
!
!       When NBDCND has any other value, BDRS is a dummy variable.
!
!     BDRF
!       A one-dimensional array of length M+1 that specifies the values
!       of the derivative of the solution with respect to R at R = RF.
!       When NBDCND = 2,3, or 6,
!
!            BDRF(I) = (d/dR)U(THETA(I),RF), I = 1,2,...,M+1  .
!
!       When NBDCND has any other value, BDRF is a dummy variable.
!
!     ELMBDA
!       The constant LAMBDA in the Helmholtz equation.  If
!       LAMBDA  >  0, a solution may not exist.  However, HWSCSP will
!       attempt to find a solution.  If NBDCND = 5 or 6 or
!       MBDCND = 5,6,7,8, or 9, ELMBDA must be zero.
!
!     F
!       A two-dimensional array that specifies the value of the right
!       side of the Helmholtz equation and boundary values (if any).
!       for I = 2,3,...,M and J = 2,3,...,N
!
!            F(I,J) = F(THETA(I),R(J)).
!
!       On the boundaries F is defined by
!
!            MBDCND   F(1,J)            F(M+1,J)
!            ------   ----------        ----------
!
!              1      U(TS,R(J))        U(TF,R(J))
!              2      U(TS,R(J))        F(TF,R(J))
!              3      F(TS,R(J))        F(TF,R(J))
!              4      F(TS,R(J))        U(TF,R(J))
!              5      F(0,R(J))         U(TF,R(J))   J = 1,2,...,N+1
!              6      F(0,R(J))         F(TF,R(J))
!              7      U(TS,R(J))        F(PI,R(J))
!              8      F(TS,R(J))        F(PI,R(J))
!              9      F(0,R(J))         F(PI,R(J))
!
!            NBDCND   F(I,1)            F(I,N+1)
!            ------   --------------    --------------
!
!              1      U(THETA(I),RS)    U(THETA(I),RF)
!              2      U(THETA(I),RS)    F(THETA(I),RF)
!              3      F(THETA(I),RS)    F(THETA(I),RF)
!              4      F(THETA(I),RS)    U(THETA(I),RF)   I = 1,2,...,M+1
!              5      F(TS,0)           U(THETA(I),RF)
!              6      F(TS,0)           F(THETA(I),RF)
!
!       F must be dimensioned at least (M+1)*(N+1).
!
!       NOTE
!
!       If the table calls for both the solution U and the right side F
!       at a corner then the solution must be specified.
!
!     IDIMF
!       The row (or first) dimension of the array F as it appears in the
!       program calling HWSCSP.  This parameter is used to specify the
!       variable dimension of F.  IDIMF must be at least M+1  .
!
!     W
!       A one-dimensional array that must be provided by the user for
!       work space. Its length can be computed from the formula below
!       which depends on the value of NBDCND.
!
!       If NBDCND=2,4 or 6 define NUNK=N
!       If NBDCND=1 or 5   define NUNK=N-1
!       If NBDCND=3        define NUNK=N+1
!
!       Now set K=INT(log2(NUNK))+1 and L=2**(K+1) then W must be
!       dimensioned at least (K-2)*L+K+5*(M+N)+MAX(2*N,6*M)+23
!
!       **IMPORTANT** For purposes of checking, the required length
!                     of W is computed by HWSCSP and stored in W(1)
!                     in floating point format.
!
!
!     * * * * * * * * * *     On Output     * * * * * * * * * *
!
!     F
!       Contains the solution U(I,J) of the finite difference
!       approximation for the grid point (THETA(I),R(J)),
!       I = 1,2,...,M+1,   J = 1,2,...,N+1  .
!
!     PERTRB
!       If a combination of periodic or derivative boundary conditions
!       is specified for a Poisson equation (LAMBDA = 0), a solution may
!       not exist.  PERTRB is a constant, calculated and subtracted from
!       F, which ensures that a solution exists.  HWSCSP then computes
!       this solution, which is a least squares solution to the original
!       approximation. This solution is not unique and is unnormalized.
!       The value of PERTRB should be small compared to the right side
!       F. Otherwise , a solution is obtained to an essentially
!       different problem. This comparison should always be made to
!       insure that a meaningful solution has been obtained.
!
!     IERROR
!       An error flag that indicates invalid input parameters.  Except
!       for numbers 0 and 10, a solution is not attempted.
!
!       = 1  TS < 0. or TF > PI
!       = 2  TS >= TF
!       = 3  M < 5
!       = 4  MBDCND < 1 or MBDCND > 9
!       = 5  RS < 0
!       = 6  RS >= RF
!       = 7  N < 5
!       = 8  NBDCND < 1 or NBDCND > 6
!       = 9  ELMBDA > 0
!       = 10 IDIMF < M+1
!       = 11 ELMBDA /= 0 and MBDCND >= 5
!       = 12 ELMBDA /= 0 and NBDCND equals 5 or 6
!       = 13 MBDCND equals 5,6 or 9 and TS /= 0
!       = 14 MBDCND >= 7 and TF /= PI
!       = 15 TS == 0 and MBDCND equals 3,4 or 8
!       = 16 TF == PI and MBDCND equals 2,3 or 6
!       = 17 NBDCND >= 5 and RS /= 0
!       = 18 NBDCND >= 5 and MBDCND equals 1,2,4,5 or 7
!
!       Since this is the only means of indicating a possibly incorrect
!       call to HWSCSP, the user should test IERROR after a call.
!
!     W
!       Contains intermediate values that must not be destroyed if
!       HWSCSP will be called again with INTL = 1.  W(1) contains the
!       number of locations which W must have.
!
! *Long Description:
!
!     * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension of   BDTS(N+1),BDTF(N+1),BDRS(M+1),BDRF(M+1),
!     Arguments      F(IDIMF,N+1),W(see argument list)
!
!     Latest         June 1979
!     Revision
!
!     Subprograms    HWSCSP,HWSCS1,BLKTRI,BLKTR1,PROD,PRODP,CPROD,CPRODP
!     Required       ,COMBP,PPADD,PSGF,BSRH,PPSGF,PPSPF,TEVLS,INDXA,
!                    ,INDXB,INDXC,R1MACH
!
!     Special
!     Conditions
!
!     Common         CBLKT
!     Blocks
!
!     I/O            NONE
!
!     Precision      Single
!
!     Specialist     Paul N Swarztrauber
!
!     Language       FORTRAN
!
!     History        Version 1 September 1973
!                    Version 2 April     1976
!                    Version 3 June      1979
!
!     Algorithm      The routine defines the finite difference
!                    equations, incorporates boundary data, and adjusts
!                    the right side of singular systems and then calls
!                    BLKTRI to solve the system.
!
!     Space
!     Required
!
!     Portability    American National Standards Institute FORTRAN.
!                    The machine accuracy is set using function R1MACH.
!
!     Required       NONE
!     Resident
!     Routines
!
!     Reference      Swarztrauber,P. and R. Sweet, 'Efficient FORTRAN
!                    Subprograms for The Solution Of Elliptic Equations'
!                    NCAR TN/IA-109, July, 1975, 138 pp.
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
!                 subprograms for the solution of elliptic equations,
!                 NCAR TN/IA-109, July 1975, 138 pp.
!***ROUTINES CALLED  HWSCS1, PIMACH
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  HWSCSP
!
  DIMENSION       F(IDIMF,*) ,BDTS(*)    ,BDTF(*)    ,BDRS(*)    , &
                  BDRF(*)    ,W(*)
!***FIRST EXECUTABLE STATEMENT  HWSCSP
  PI = PIMACH(DUM)
  IERROR = 0
  if (TS < 0. .OR. TF > PI) IERROR = 1
  if (TS  >=  TF) IERROR = 2
  if (M  <  5) IERROR = 3
  if (MBDCND < 1 .OR. MBDCND > 9) IERROR = 4
  if (RS  <  0.) IERROR = 5
  if (RS  >=  RF) IERROR = 6
  if (N  <  5) IERROR = 7
  if (NBDCND < 1 .OR. NBDCND > 6) IERROR = 8
  if (ELMBDA  >  0.) IERROR = 9
  if (IDIMF  <  M+1) IERROR = 10
  if (ELMBDA /= 0. .AND. MBDCND >= 5) IERROR = 11
  if (ELMBDA /= 0. .AND. (NBDCND == 5 .OR. NBDCND == 6)) IERROR = 12
  if ((MBDCND == 5 .OR. MBDCND == 6 .OR. MBDCND == 9) .AND. &
      TS /= 0.) IERROR = 13
  if (MBDCND >= 7 .AND. TF /= PI) IERROR = 14
  if (TS == 0. .AND. &
      (MBDCND == 4 .OR. MBDCND == 8 .OR. MBDCND == 3)) IERROR = 15
  if (TF == PI .AND. &
      (MBDCND == 2 .OR. MBDCND == 3 .OR. MBDCND == 6)) IERROR = 16
  if (NBDCND >= 5 .AND. RS /= 0.) IERROR = 17
  if (NBDCND >= 5 .AND. (MBDCND == 1 .OR. MBDCND == 2 .OR. &
                                      MBDCND == 5 .OR. MBDCND == 7)) &
      IERROR = 18
  if (IERROR /= 0 .AND. IERROR /= 9) RETURN
  NCK = N
  go to (101,103,102,103,101,103),NBDCND
  101 NCK = NCK-1
  go to 103
  102 NCK = NCK+1
  103 L = 2
  K = 1
  104 L = L+L
  K = K+1
  if (NCK-L) 105,105,104
  105 L = L+L
  NP1 = N+1
  MP1 = M+1
  I1 = (K-2)*L+K+MAX(2*N,6*M)+13
  I2 = I1+NP1
  I3 = I2+NP1
  I4 = I3+NP1
  I5 = I4+NP1
  I6 = I5+NP1
  I7 = I6+MP1
  I8 = I7+MP1
  I9 = I8+MP1
  I10 = I9+MP1
  W(1) = I10+M
  call HWSCS1 (INTL,TS,TF,M,MBDCND,BDTS,BDTF,RS,RF,N,NBDCND,BDRS, &
               BDRF,ELMBDA,F,IDIMF,PERTRB,W(2),W(I1),W(I2),W(I3), &
               W(I4),W(I5),W(I6),W(I7),W(I8),W(I9),W(I10))
  return
end
