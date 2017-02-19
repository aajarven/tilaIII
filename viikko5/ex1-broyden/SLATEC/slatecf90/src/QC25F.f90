subroutine QC25F (F, A, B, OMEGA, INTEGR, NRMOM, MAXP1, KSAVE, &
     RESULT, ABSERR, NEVAL, RESABS, RESASC, MOMCOM, CHEBMO)
!
!! QC25F computes the integral I=Integral of F(X) over (A,B) ...
!            Where W(X) = COS(OMEGA*X) Or (WX)=SIN(OMEGA*X)
!            and to compute J=Integral of ABS(F) over (A,B). For small
!            value of OMEGA or small intervals (A,B) 15-point GAUSS-
!            KRONROD Rule used. Otherwise generalized CLENSHAW-CURTIS is used.
!
!***LIBRARY   SLATEC (QUADPACK)
!***CATEGORY  H2A2A2
!***TYPE      SINGLE PRECISION (QC25F-S, DQC25F-D)
!***KEYWORDS  CLENSHAW-CURTIS METHOD, GAUSS-KRONROD RULES,
!             INTEGRATION RULES FOR FUNCTIONS WITH COS OR SIN FACTOR,
!             QUADPACK, QUADRATURE
!***AUTHOR  Piessens, Robert
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!           de Doncker, Elise
!             Applied Mathematics and Programming Division
!             K. U. Leuven
!***DESCRIPTION
!
!        Integration rules for functions with COS or SIN factor
!        Standard fortran subroutine
!        Real version
!
!        PARAMETERS
!         ON ENTRY
!           F      - Real
!                    Function subprogram defining the integrand
!                    function F(X). The actual name for F needs to
!                    be declared E X T E R N A L in the calling program.
!
!           A      - Real
!                    Lower limit of integration
!
!           B      - Real
!                    Upper limit of integration
!
!           OMEGA  - Real
!                    Parameter in the WEIGHT function
!
!           INTEGR - Integer
!                    Indicates which WEIGHT function is to be used
!                       INTEGR = 1   W(X) = COS(OMEGA*X)
!                       INTEGR = 2   W(X) = SIN(OMEGA*X)
!
!           NRMOM  - Integer
!                    The length of interval (A,B) is equal to the length
!                    of the original integration interval divided by
!                    2**NRMOM (we suppose that the routine is used in an
!                    adaptive integration process, otherwise set
!                    NRMOM = 0). NRMOM must be zero at the first call.
!
!           MAXP1  - Integer
!                    Gives an upper bound on the number of Chebyshev
!                    moments which can be stored, i.e. for the
!                    intervals of lengths ABS(BB-AA)*2**(-L),
!                    L = 0,1,2, ..., MAXP1-2.
!
!           KSAVE  - Integer
!                    Key which is one when the moments for the
!                    current interval have been computed
!
!         ON RETURN
!           RESULT - Real
!                    Approximation to the integral I
!
!           ABSERR - Real
!                    Estimate of the modulus of the absolute
!                    error, which should equal or exceed ABS(I-RESULT)
!
!           NEVAL  - Integer
!                    Number of integrand evaluations
!
!           RESABS - Real
!                    Approximation to the integral J
!
!           RESASC - Real
!                    Approximation to the integral of ABS(F-I/(B-A))
!
!         ON ENTRY AND RETURN
!           MOMCOM - Integer
!                    For each interval length we need to compute the
!                    Chebyshev moments. MOMCOM counts the number of
!                    intervals for which these moments have already been
!                    computed. If NRMOM < MOMCOM or KSAVE = 1, the
!                    Chebyshev moments for the interval (A,B) have
!                    already been computed and stored, otherwise we
!                    compute them and we increase MOMCOM.
!
!           CHEBMO - Real
!                    Array of dimension at least (MAXP1,25) containing
!                    the modified Chebyshev moments for the first MOMCOM
!                    MOMCOM interval lengths
!
!***REFERENCES  (NONE)
!***ROUTINES CALLED  QCHEB, QK15W, QWGTF, R1MACH, SGTSL
!***REVISION HISTORY  (YYMMDD)
!   810101  DATE WRITTEN
!   861211  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  QC25F
!
  REAL A,ABSERR,AC,AN,AN2,AS,ASAP,ASS,B,CENTR,CHEBMO, &
    CHEB12,CHEB24,CONC,CONS,COSPAR,D,QWGTF, &
    D1,R1MACH,D2,ESTC,ESTS,F,FVAL,HLGTH,OFLOW,OMEGA,PARINT,PAR2, &
    PAR22,P2,P3,P4,RESABS,RESASC,RESC12,RESC24,RESS12,RESS24, &
    RESULT,SINPAR,V,X
  INTEGER I,IERS,INTEGR,ISYM,J,K,KSAVE,M,MAXP1,MOMCOM,NEVAL, &
    NOEQU,NOEQ1,NRMOM
!
  DIMENSION CHEBMO(MAXP1,25),CHEB12(13),CHEB24(25),D(25),D1(25), &
    D2(25),FVAL(25),V(28),X(11)
!
  EXTERNAL F, QWGTF
!
!           THE VECTOR X CONTAINS THE VALUES COS(K*PI/24)
!           K = 1, ...,11, TO BE USED FOR THE CHEBYSHEV EXPANSION OF F
!
  SAVE X
  DATA X(1),X(2),X(3),X(4),X(5),X(6),X(7),X(8),X(9), &
    X(10),X(11)/ &
       0.9914448613738104E+00,     0.9659258262890683E+00, &
       0.9238795325112868E+00,     0.8660254037844386E+00, &
       0.7933533402912352E+00,     0.7071067811865475E+00, &
       0.6087614290087206E+00,     0.5000000000000000E+00, &
       0.3826834323650898E+00,     0.2588190451025208E+00, &
       0.1305261922200516E+00/
!
!           LIST OF MAJOR VARIABLES
!           -----------------------
!
!           CENTR  - MID POINT OF THE INTEGRATION INTERVAL
!           HLGTH  - HALF-LENGTH OF THE INTEGRATION INTERVAL
!           FVAL   - VALUE OF THE FUNCTION F AT THE POINTS
!                    (B-A)*0.5*COS(K*PI/12) + (B+A)*0.5,
!                    K = 0, ..., 24
!           CHEB12 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
!                    OF DEGREE 12, FOR THE FUNCTION F, IN THE
!                    INTERVAL (A,B)
!           CHEB24 - COEFFICIENTS OF THE CHEBYSHEV SERIES EXPANSION
!                    OF DEGREE 24, FOR THE FUNCTION F, IN THE
!                    INTERVAL (A,B)
!           RESC12 - APPROXIMATION TO THE INTEGRAL OF
!                    COS(0.5*(B-A)*OMEGA*X)*F(0.5*(B-A)*X+0.5*(B+A))
!                    OVER (-1,+1), USING THE CHEBYSHEV SERIES
!                    EXPANSION OF DEGREE 12
!           RESC24 - APPROXIMATION TO THE SAME INTEGRAL, USING THE
!                    CHEBYSHEV SERIES EXPANSION OF DEGREE 24
!           RESS12 - THE ANALOGUE OF RESC12 FOR THE SINE
!           RESS24 - THE ANALOGUE OF RESC24 FOR THE SINE
!
!
!           MACHINE DEPENDENT CONSTANT
!           --------------------------
!
!           OFLOW IS THE LARGEST POSITIVE MAGNITUDE.
!
!***FIRST EXECUTABLE STATEMENT  QC25F
  OFLOW = R1MACH(2)
!
  CENTR = 0.5E+00*(B+A)
  HLGTH = 0.5E+00*(B-A)
  PARINT = OMEGA*HLGTH
!
!           COMPUTE THE INTEGRAL USING THE 15-POINT GAUSS-KRONROD
!           FORMULA if THE VALUE OF THE PARAMETER IN THE INTEGRAND
!           IS SMALL.
!
  if ( ABS(PARINT) > 0.2E+01) go to 10
  call QK15W(F,QWGTF,OMEGA,P2,P3,P4,INTEGR,A,B,RESULT, &
    ABSERR,RESABS,RESASC)
  NEVAL = 15
  go to 170
!
!           COMPUTE THE INTEGRAL USING THE GENERALIZED CLENSHAW-
!           CURTIS METHOD.
!
   10 CONC = HLGTH*COS(CENTR*OMEGA)
  CONS = HLGTH*SIN(CENTR*OMEGA)
  RESASC = OFLOW
  NEVAL = 25
!
!           CHECK WHETHER THE CHEBYSHEV MOMENTS FOR THIS INTERVAL
!           HAVE ALREADY BEEN COMPUTED.
!
  if ( NRMOM < MOMCOM.OR.KSAVE == 1) go to 120
!
!           COMPUTE A NEW SET OF CHEBYSHEV MOMENTS.
!
  M = MOMCOM+1
  PAR2 = PARINT*PARINT
  PAR22 = PAR2+0.2E+01
  SINPAR = SIN(PARINT)
  COSPAR = COS(PARINT)
!
!           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO COSINE.
!
  V(1) = 0.2E+01*SINPAR/PARINT
  V(2) = (0.8E+01*COSPAR+(PAR2+PAR2-0.8E+01)*SINPAR/ &
    PARINT)/PAR2
  V(3) = (0.32E+02*(PAR2-0.12E+02)*COSPAR+(0.2E+01* &
    ((PAR2-0.80E+02)*PAR2+0.192E+03)*SINPAR)/ &
    PARINT)/(PAR2*PAR2)
  AC = 0.8E+01*COSPAR
  AS = 0.24E+02*PARINT*SINPAR
  if ( ABS(PARINT) > 0.24E+02) go to 30
!
!           COMPUTE THE CHEBYSHEV MOMENTS AS THE
!           SOLUTIONS OF A BOUNDARY VALUE PROBLEM WITH 1
!           INITIAL VALUE (V(3)) AND 1 END VALUE (COMPUTED
!           USING AN ASYMPTOTIC FORMULA).
!
  NOEQU = 25
  NOEQ1 = NOEQU-1
  AN = 0.6E+01
  DO 20 K = 1,NOEQ1
    AN2 = AN*AN
    D(K) = -0.2E+01*(AN2-0.4E+01)*(PAR22-AN2-AN2)
    D2(K) = (AN-0.1E+01)*(AN-0.2E+01)*PAR2
    D1(K+1) = (AN+0.3E+01)*(AN+0.4E+01)*PAR2
    V(K+3) = AS-(AN2-0.4E+01)*AC
    AN = AN+0.2E+01
   20 CONTINUE
  AN2 = AN*AN
  D(NOEQU) = -0.2E+01*(AN2-0.4E+01)*(PAR22-AN2-AN2)
  V(NOEQU+3) = AS-(AN2-0.4E+01)*AC
  V(4) = V(4)-0.56E+02*PAR2*V(3)
  ASS = PARINT*SINPAR
  ASAP = (((((0.210E+03*PAR2-0.1E+01)*COSPAR-(0.105E+03*PAR2 &
    -0.63E+02)*ASS)/AN2-(0.1E+01-0.15E+02*PAR2)*COSPAR &
    +0.15E+02*ASS)/AN2-COSPAR+0.3E+01*ASS)/AN2-COSPAR)/AN2
  V(NOEQU+3) = V(NOEQU+3)-0.2E+01*ASAP*PAR2*(AN-0.1E+01)* &
     (AN-0.2E+01)
!
!           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN
!           ELIMINATION WITH PARTIAL PIVOTING.
!
  call SGTSL(NOEQU,D1,D,D2,V(4),IERS)
  go to 50
!
!           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF FORWARD
!           RECURSION.
!
   30 AN = 0.4E+01
  DO 40 I = 4,13
    AN2 = AN*AN
    V(I) = ((AN2-0.4E+01)*(0.2E+01*(PAR22-AN2-AN2)*V(I-1)-AC) &
    +AS-PAR2*(AN+0.1E+01)*(AN+0.2E+01)*V(I-2))/ &
    (PAR2*(AN-0.1E+01)*(AN-0.2E+01))
    AN = AN+0.2E+01
   40 CONTINUE
   50 DO 60 J = 1,13
    CHEBMO(M,2*J-1) = V(J)
   60 CONTINUE
!
!           COMPUTE THE CHEBYSHEV MOMENTS WITH RESPECT TO SINE.
!
  V(1) = 0.2E+01*(SINPAR-PARINT*COSPAR)/PAR2
  V(2) = (0.18E+02-0.48E+02/PAR2)*SINPAR/PAR2 &
    +(-0.2E+01+0.48E+02/PAR2)*COSPAR/PARINT
  AC = -0.24E+02*PARINT*COSPAR
  AS = -0.8E+01*SINPAR
  if ( ABS(PARINT) > 0.24E+02) go to 80
!
!           COMPUTE THE CHEBYSHEV MOMENTS AS THE
!           SOLUTIONS OF A BOUNDARY VALUE PROBLEM WITH 1
!           INITIAL VALUE (V(2)) AND 1 END VALUE (COMPUTED
!           USING AN ASYMPTOTIC FORMULA).
!
  AN = 0.5E+01
  DO 70 K = 1,NOEQ1
    AN2 = AN*AN
    D(K) = -0.2E+01*(AN2-0.4E+01)*(PAR22-AN2-AN2)
    D2(K) = (AN-0.1E+01)*(AN-0.2E+01)*PAR2
    D1(K+1) = (AN+0.3E+01)*(AN+0.4E+01)*PAR2
    V(K+2) = AC+(AN2-0.4E+01)*AS
    AN = AN+0.2E+01
   70 CONTINUE
  AN2 = AN*AN
  D(NOEQU) = -0.2E+01*(AN2-0.4E+01)*(PAR22-AN2-AN2)
  V(NOEQU+2) = AC+(AN2-0.4E+01)*AS
  V(3) = V(3)-0.42E+02*PAR2*V(2)
  ASS = PARINT*COSPAR
  ASAP = (((((0.105E+03*PAR2-0.63E+02)*ASS+(0.210E+03*PAR2 &
    -0.1E+01)*SINPAR)/AN2+(0.15E+02*PAR2-0.1E+01)*SINPAR- &
    0.15E+02*ASS)/AN2-0.3E+01*ASS-SINPAR)/AN2-SINPAR)/AN2
  V(NOEQU+2) = V(NOEQU+2)-0.2E+01*ASAP*PAR2*(AN-0.1E+01) &
    *(AN-0.2E+01)
!
!           SOLVE THE TRIDIAGONAL SYSTEM BY MEANS OF GAUSSIAN
!           ELIMINATION WITH PARTIAL PIVOTING.
!
  call SGTSL(NOEQU,D1,D,D2,V(3),IERS)
  go to 100
!
!           COMPUTE THE CHEBYSHEV MOMENTS BY MEANS OF
!           FORWARD RECURSION.
!
   80 AN = 0.3E+01
  DO 90 I = 3,12
    AN2 = AN*AN
    V(I) = ((AN2-0.4E+01)*(0.2E+01*(PAR22-AN2-AN2)*V(I-1)+AS) &
    +AC-PAR2*(AN+0.1E+01)*(AN+0.2E+01)*V(I-2)) &
    /(PAR2*(AN-0.1E+01)*(AN-0.2E+01))
    AN = AN+0.2E+01
   90 CONTINUE
  100 DO 110 J = 1,12
    CHEBMO(M,2*J) = V(J)
  110 CONTINUE
  120 if (NRMOM < MOMCOM) M = NRMOM+1
   if (MOMCOM < MAXP1-1.AND.NRMOM >= MOMCOM) MOMCOM = MOMCOM+1
!
!           COMPUTE THE COEFFICIENTS OF THE CHEBYSHEV EXPANSIONS
!           OF DEGREES 12 AND 24 OF THE FUNCTION F.
!
  FVAL(1) = 0.5E+00*F(CENTR+HLGTH)
  FVAL(13) = F(CENTR)
  FVAL(25) = 0.5E+00*F(CENTR-HLGTH)
  DO 130 I = 2,12
    ISYM = 26-I
    FVAL(I) = F(HLGTH*X(I-1)+CENTR)
    FVAL(ISYM) = F(CENTR-HLGTH*X(I-1))
  130 CONTINUE
  call QCHEB(X,FVAL,CHEB12,CHEB24)
!
!           COMPUTE THE INTEGRAL AND ERROR ESTIMATES.
!
  RESC12 = CHEB12(13)*CHEBMO(M,13)
  RESS12 = 0.0E+00
  K = 11
  DO 140 J = 1,6
    RESC12 = RESC12+CHEB12(K)*CHEBMO(M,K)
    RESS12 = RESS12+CHEB12(K+1)*CHEBMO(M,K+1)
    K = K-2
  140 CONTINUE
  RESC24 = CHEB24(25)*CHEBMO(M,25)
  RESS24 = 0.0E+00
  RESABS = ABS(CHEB24(25))
  K = 23
  DO 150 J = 1,12
    RESC24 = RESC24+CHEB24(K)*CHEBMO(M,K)
    RESS24 = RESS24+CHEB24(K+1)*CHEBMO(M,K+1)
    RESABS = ABS(CHEB24(K))+ABS(CHEB24(K+1))
    K = K-2
  150 CONTINUE
  ESTC = ABS(RESC24-RESC12)
  ESTS = ABS(RESS24-RESS12)
  RESABS = RESABS*ABS(HLGTH)
  if ( INTEGR == 2) go to 160
  RESULT = CONC*RESC24-CONS*RESS24
  ABSERR = ABS(CONC*ESTC)+ABS(CONS*ESTS)
  go to 170
  160 RESULT = CONC*RESS24+CONS*RESC24
  ABSERR = ABS(CONC*ESTS)+ABS(CONS*ESTC)
  170 RETURN
end
