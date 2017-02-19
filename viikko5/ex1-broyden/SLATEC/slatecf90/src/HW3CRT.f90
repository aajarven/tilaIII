subroutine HW3CRT (XS, XF, L, LBDCND, BDXS, BDXF, YS, YF, M, &
     MBDCND, BDYS, BDYF, ZS, ZF, N, NBDCND, BDZS, BDZF, ELMBDA, &
     LDIMF, MDIMF, F, PERTRB, IERROR, W)
!
!! HW3CRT solves the standard seven-point finite difference ...
!            approximation to the Helmholtz equation in Cartesian
!            coordinates.
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B1A1A
!***TYPE      SINGLE PRECISION (HW3CRT-S)
!***KEYWORDS  CARTESIAN, ELLIPTIC, FISHPACK, HELMHOLTZ, PDE
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine HW3CRT solves the standard seven-point finite
!     difference approximation to the Helmholtz equation in Cartesian
!     coordinates:
!
!         (d/dX)(dU/dX) + (d/dY)(dU/dY) + (d/dZ)(dU/dZ)
!
!                    + LAMBDA*U = F(X,Y,Z) .
!
!    * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!
!    * * * * * * * *    Parameter Description     * * * * * * * * * *
!
!
!            * * * * * *   On Input    * * * * * *
!
!     XS,XF
!        The range of X, i.e. XS  <=  X  <=  XF .
!        XS must be less than XF.
!
!     L
!        The number of panels into which the interval (XS,XF) is
!        subdivided.  Hence, there will be L+1 grid points in the
!        X-direction given by X(I) = XS+(I-1)DX for I=1,2,...,L+1,
!        where DX = (XF-XS)/L is the panel width.  L must be at
!        least 5 .
!
!     LBDCND
!        Indicates the type of boundary conditions at X = XS and X = XF.
!
!        = 0  If the solution is periodic in X, i.e.
!             U(L+I,J,K) = U(I,J,K).
!        = 1  If the solution is specified at X = XS and X = XF.
!        = 2  If the solution is specified at X = XS and the derivative
!             of the solution with respect to X is specified at X = XF.
!        = 3  If the derivative of the solution with respect to X is
!             specified at X = XS and X = XF.
!        = 4  If the derivative of the solution with respect to X is
!             specified at X = XS and the solution is specified at X=XF.
!
!     BDXS
!        A two-dimensional array that specifies the values of the
!        derivative of the solution with respect to X at X = XS.
!        when LBDCND = 3 or 4,
!
!             BDXS(J,K) = (d/dX)U(XS,Y(J),Z(K)), J=1,2,...,M+1,
!                                                K=1,2,...,N+1.
!
!        When LBDCND has any other value, BDXS is a dummy variable.
!        BDXS must be dimensioned at least (M+1)*(N+1).
!
!     BDXF
!        A two-dimensional array that specifies the values of the
!        derivative of the solution with respect to X at X = XF.
!        When LBDCND = 2 or 3,
!
!             BDXF(J,K) = (d/dX)U(XF,Y(J),Z(K)), J=1,2,...,M+1,
!                                                K=1,2,...,N+1.
!
!        When LBDCND has any other value, BDXF is a dummy variable.
!        BDXF must be dimensioned at least (M+1)*(N+1).
!
!     YS,YF
!        The range of Y, i.e. YS  <=  Y  <=  YF.
!        YS must be less than YF.
!
!     M
!        The number of panels into which the interval (YS,YF) is
!        subdivided.  Hence, there will be M+1 grid points in the
!        Y-direction given by Y(J) = YS+(J-1)DY for J=1,2,...,M+1,
!        where DY = (YF-YS)/M is the panel width.  M must be at
!        least 5 .
!
!     MBDCND
!        Indicates the type of boundary conditions at Y = YS and Y = YF.
!
!        = 0  If the solution is periodic in Y, i.e.
!             U(I,M+J,K) = U(I,J,K).
!        = 1  If the solution is specified at Y = YS and Y = YF.
!        = 2  If the solution is specified at Y = YS and the derivative
!             of the solution with respect to Y is specified at Y = YF.
!        = 3  If the derivative of the solution with respect to Y is
!             specified at Y = YS and Y = YF.
!        = 4  If the derivative of the solution with respect to Y is
!             specified at Y = YS and the solution is specified at Y=YF.
!
!     BDYS
!        A two-dimensional array that specifies the values of the
!        derivative of the solution with respect to Y at Y = YS.
!        When MBDCND = 3 or 4,
!
!             BDYS(I,K) = (d/dY)U(X(I),YS,Z(K)), I=1,2,...,L+1,
!                                                K=1,2,...,N+1.
!
!        When MBDCND has any other value, BDYS is a dummy variable.
!        BDYS must be dimensioned at least (L+1)*(N+1).
!
!     BDYF
!        A two-dimensional array that specifies the values of the
!        derivative of the solution with respect to Y at Y = YF.
!        When MBDCND = 2 or 3,
!
!             BDYF(I,K) = (d/dY)U(X(I),YF,Z(K)), I=1,2,...,L+1,
!                                                K=1,2,...,N+1.
!
!        When MBDCND has any other value, BDYF is a dummy variable.
!        BDYF must be dimensioned at least (L+1)*(N+1).
!
!     ZS,ZF
!        The range of Z, i.e. ZS  <=  Z  <=  ZF.
!        ZS must be less than ZF.
!
!     N
!        The number of panels into which the interval (ZS,ZF) is
!        subdivided.  Hence, there will be N+1 grid points in the
!        Z-direction given by Z(K) = ZS+(K-1)DZ for K=1,2,...,N+1,
!        where DZ = (ZF-ZS)/N is the panel width.  N must be at least 5.
!
!     NBDCND
!        Indicates the type of boundary conditions at Z = ZS and Z = ZF.
!
!        = 0  If the solution is periodic in Z, i.e.
!             U(I,J,N+K) = U(I,J,K).
!        = 1  If the solution is specified at Z = ZS and Z = ZF.
!        = 2  If the solution is specified at Z = ZS and the derivative
!             of the solution with respect to Z is specified at Z = ZF.
!        = 3  If the derivative of the solution with respect to Z is
!             specified at Z = ZS and Z = ZF.
!        = 4  If the derivative of the solution with respect to Z is
!             specified at Z = ZS and the solution is specified at Z=ZF.
!
!     BDZS
!        A two-dimensional array that specifies the values of the
!        derivative of the solution with respect to Z at Z = ZS.
!        When NBDCND = 3 or 4,
!
!             BDZS(I,J) = (d/dZ)U(X(I),Y(J),ZS), I=1,2,...,L+1,
!                                                J=1,2,...,M+1.
!
!        When NBDCND has any other value, BDZS is a dummy variable.
!        BDZS must be dimensioned at least (L+1)*(M+1).
!
!     BDZF
!        A two-dimensional array that specifies the values of the
!        derivative of the solution with respect to Z at Z = ZF.
!        When NBDCND = 2 or 3,
!
!             BDZF(I,J) = (d/dZ)U(X(I),Y(J),ZF), I=1,2,...,L+1,
!                                                J=1,2,...,M+1.
!
!        When NBDCND has any other value, BDZF is a dummy variable.
!        BDZF must be dimensioned at least (L+1)*(M+1).
!
!     ELMBDA
!        The constant LAMBDA in the Helmholtz equation. If
!        LAMBDA  >  0, a solution may not exist.  However, HW3CRT will
!        attempt to find a solution.
!
!     F
!        A three-dimensional array that specifies the values of the
!        right side of the Helmholtz equation and boundary values (if
!        any).  For I=2,3,...,L, J=2,3,...,M, and K=2,3,...,N
!
!                   F(I,J,K) = F(X(I),Y(J),Z(K)).
!
!        On the boundaries F is defined by
!
!        LBDCND      F(1,J,K)         F(L+1,J,K)
!        ------   ---------------   ---------------
!
!          0      F(XS,Y(J),Z(K))   F(XS,Y(J),Z(K))
!          1      U(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
!          2      U(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   J=1,2,...,M+1
!          3      F(XS,Y(J),Z(K))   F(XF,Y(J),Z(K))   K=1,2,...,N+1
!          4      F(XS,Y(J),Z(K))   U(XF,Y(J),Z(K))
!
!        MBDCND      F(I,1,K)         F(I,M+1,K)
!        ------   ---------------   ---------------
!
!          0      F(X(I),YS,Z(K))   F(X(I),YS,Z(K))
!          1      U(X(I),YS,Z(K))   U(X(I),YF,Z(K))
!          2      U(X(I),YS,Z(K))   F(X(I),YF,Z(K))   I=1,2,...,L+1
!          3      F(X(I),YS,Z(K))   F(X(I),YF,Z(K))   K=1,2,...,N+1
!          4      F(X(I),YS,Z(K))   U(X(I),YF,Z(K))
!
!        NBDCND      F(I,J,1)         F(I,J,N+1)
!        ------   ---------------   ---------------
!
!          0      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZS)
!          1      U(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
!          2      U(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   I=1,2,...,L+1
!          3      F(X(I),Y(J),ZS)   F(X(I),Y(J),ZF)   J=1,2,...,M+1
!          4      F(X(I),Y(J),ZS)   U(X(I),Y(J),ZF)
!
!        F must be dimensioned at least (L+1)*(M+1)*(N+1).
!
!        NOTE:
!
!        If the table calls for both the solution U and the right side F
!        on a boundary, then the solution must be specified.
!
!     LDIMF
!        The row (or first) dimension of the arrays F,BDYS,BDYF,BDZS,
!        and BDZF as it appears in the program calling HW3CRT. this
!        parameter is used to specify the variable dimension of these
!        arrays.  LDIMF must be at least L+1.
!
!     MDIMF
!        The column (or second) dimension of the array F and the row (or
!        first) dimension of the arrays BDXS and BDXF as it appears in
!        the program calling HW3CRT.  This parameter is used to specify
!        the variable dimension of these arrays.
!        MDIMF must be at least M+1.
!
!     W
!        A one-dimensional array that must be provided by the user for
!        work space.  The length of W must be at least 30 + L + M + 5*N
!        + MAX(L,M,N) + 7*(INT((L+1)/2) + INT((M+1)/2))
!
!
!            * * * * * *   On Output   * * * * * *
!
!     F
!        Contains the solution U(I,J,K) of the finite difference
!        approximation for the grid point (X(I),Y(J),Z(K)) for
!        I=1,2,...,L+1, J=1,2,...,M+1, and K=1,2,...,N+1.
!
!     PERTRB
!        If a combination of periodic or derivative boundary conditions
!        is specified for a Poisson equation (LAMBDA = 0), a solution
!        may not exist.  PERTRB is a constant, calculated and subtracted
!        from F, which ensures that a solution exists.  PWSCRT then
!        computes this solution, which is a least squares solution to
!        the original approximation.  This solution is not unique and is
!        unnormalized.  The value of PERTRB should be small compared to
!        the right side F.  Otherwise, a solution is obtained to an
!        essentially different problem.  This comparison should always
!        be made to insure that a meaningful solution has been obtained.
!
!     IERROR
!        An error flag that indicates invalid input parameters.  Except
!        for numbers 0 and 12, a solution is not attempted.
!
!        =  0  No error
!        =  1  XS  >=  XF
!        =  2  L  <  5
!        =  3  LBDCND  <  0 .OR. LBDCND  >  4
!        =  4  YS  >=  YF
!        =  5  M  <  5
!        =  6  MBDCND  <  0 .OR. MBDCND  >  4
!        =  7  ZS  >=  ZF
!        =  8  N  <  5
!        =  9  NBDCND  <  0 .OR. NBDCND  >  4
!        = 10  LDIMF  <  L+1
!        = 11  MDIMF  <  M+1
!        = 12  LAMBDA  >  0
!
!        Since this is the only means of indicating a possibly incorrect
!        call to HW3CRT, the user should test IERROR after the call.
!
! *Long Description:
!
!    * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension of   BDXS(MDIMF,N+1),BDXF(MDIMF,N+1),BDYS(LDIMF,N+1),
!     Arguments      BDYF(LDIMF,N+1),BDZS(LDIMF,M+1),BDZF(LDIMF,M+1),
!                    F(LDIMF,MDIMF,N+1),W(see argument list)
!
!     Latest         December 1, 1978
!     Revision
!
!     Subprograms    HW3CRT,POIS3D,POS3D1,TRIDQ,RFFTI,RFFTF,RFFTF1,
!     Required       RFFTB,RFFTB1,COSTI,COST,SINTI,SINT,COSQI,COSQF,
!                    COSQF1,COSQB,COSQB1,SINQI,SINQF,SINQB,CFFTI,
!                    CFFTI1,CFFTB,CFFTB1,PASSB2,PASSB3,PASSB4,PASSB,
!                    CFFTF,CFFTF1,PASSF1,PASSF2,PASSF3,PASSF4,PASSF,
!                    PIMACH
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
!     Algorithm      This subroutine defines the finite difference
!                    equations, incorporates boundary data, and
!                    adjusts the right side of singular systems and
!                    then calls POIS3D to solve the system.
!
!     Space          7862(decimal) = 17300(octal) locations on the
!     Required       NCAR Control Data 7600
!
!     Timing and        The execution time T on the NCAR Control Data
!     Accuracy       7600 for subroutine HW3CRT is roughly proportional
!                    to L*M*N*(log2(L)+log2(M)+5), but also depends on
!                    input parameters LBDCND and MBDCND.  Some typical
!                    values are listed in the table below.
!                       The solution process employed results in a loss
!                    of no more than three significant digits for L,M
!                    and N as large as 32.  More detailed information
!                    about accuracy can be found in the documentation
!                    for subroutine POIS3D which is the routine that
!                    actually solves the finite difference equations.
!
!
!                       L(=M=N)     LBDCND(=MBDCND=NBDCND)      T(MSECS)
!                       -------     ----------------------      --------
!
!                         16                  0                    300
!                         16                  1                    302
!                         16                  3                    348
!                         32                  0                   1925
!                         32                  1                   1929
!                         32                  3                   2109
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
!***ROUTINES CALLED  POIS3D
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!***END PROLOGUE  HW3CRT
!
!
  DIMENSION       BDXS(MDIMF,*)          ,BDXF(MDIMF,*)          , &
                  BDYS(LDIMF,*)          ,BDYF(LDIMF,*)          , &
                  BDZS(LDIMF,*)          ,BDZF(LDIMF,*)          , &
                  F(LDIMF,MDIMF,*)       ,W(*)
!***FIRST EXECUTABLE STATEMENT  HW3CRT
  IERROR = 0
  if (XF  <=  XS) IERROR = 1
  if (L  <  5) IERROR = 2
  if (LBDCND < 0 .OR. LBDCND > 4) IERROR = 3
  if (YF  <=  YS) IERROR = 4
  if (M  <  5) IERROR = 5
  if (MBDCND < 0 .OR. MBDCND > 4) IERROR = 6
  if (ZF  <=  ZS) IERROR = 7
  if (N  <  5) IERROR = 8
  if (NBDCND < 0 .OR. NBDCND > 4) IERROR = 9
  if (LDIMF  <  L+1) IERROR = 10
  if (MDIMF  <  M+1) IERROR = 11
  if (IERROR  /=  0) go to 188
  DY = (YF-YS)/M
  TWBYDY = 2./DY
  C2 = 1./(DY**2)
  MSTART = 1
  MSTOP = M
  MP1 = M+1
  MP = MBDCND+1
  go to (104,101,101,102,102),MP
  101 MSTART = 2
  102 go to (104,104,103,103,104),MP
  103 MSTOP = MP1
  104 MUNK = MSTOP-MSTART+1
  DZ = (ZF-ZS)/N
  TWBYDZ = 2./DZ
  NP = NBDCND+1
  C3 = 1./(DZ**2)
  NP1 = N+1
  NSTART = 1
  NSTOP = N
  go to (108,105,105,106,106),NP
  105 NSTART = 2
  106 go to (108,108,107,107,108),NP
  107 NSTOP = NP1
  108 NUNK = NSTOP-NSTART+1
  LP1 = L+1
  DX = (XF-XS)/L
  C1 = 1./(DX**2)
  TWBYDX = 2./DX
  LP = LBDCND+1
  LSTART = 1
  LSTOP = L
!
!     ENTER BOUNDARY DATA FOR X-BOUNDARIES.
!
  go to (122,109,109,112,112),LP
  109 LSTART = 2
  DO 111 J=MSTART,MSTOP
     DO 110 K=NSTART,NSTOP
        F(2,J,K) = F(2,J,K)-C1*F(1,J,K)
  110    CONTINUE
  111 CONTINUE
  go to 115
  112 DO 114 J=MSTART,MSTOP
     DO 113 K=NSTART,NSTOP
        F(1,J,K) = F(1,J,K)+TWBYDX*BDXS(J,K)
  113    CONTINUE
  114 CONTINUE
  115 go to (122,116,119,119,116),LP
  116 DO 118 J=MSTART,MSTOP
     DO 117 K=NSTART,NSTOP
        F(L,J,K) = F(L,J,K)-C1*F(LP1,J,K)
  117    CONTINUE
  118 CONTINUE
  go to 122
  119 LSTOP = LP1
  DO 121 J=MSTART,MSTOP
     DO 120 K=NSTART,NSTOP
        F(LP1,J,K) = F(LP1,J,K)-TWBYDX*BDXF(J,K)
  120    CONTINUE
  121 CONTINUE
  122 LUNK = LSTOP-LSTART+1
!
!     ENTER BOUNDARY DATA FOR Y-BOUNDARIES.
!
  go to (136,123,123,126,126),MP
  123 DO 125 I=LSTART,LSTOP
     DO 124 K=NSTART,NSTOP
        F(I,2,K) = F(I,2,K)-C2*F(I,1,K)
  124    CONTINUE
  125 CONTINUE
  go to 129
  126 DO 128 I=LSTART,LSTOP
     DO 127 K=NSTART,NSTOP
        F(I,1,K) = F(I,1,K)+TWBYDY*BDYS(I,K)
  127    CONTINUE
  128 CONTINUE
  129 go to (136,130,133,133,130),MP
  130 DO 132 I=LSTART,LSTOP
     DO 131 K=NSTART,NSTOP
        F(I,M,K) = F(I,M,K)-C2*F(I,MP1,K)
  131    CONTINUE
  132 CONTINUE
  go to 136
  133 DO 135 I=LSTART,LSTOP
     DO 134 K=NSTART,NSTOP
        F(I,MP1,K) = F(I,MP1,K)-TWBYDY*BDYF(I,K)
  134    CONTINUE
  135 CONTINUE
  136 CONTINUE
!
!     ENTER BOUNDARY DATA FOR Z-BOUNDARIES.
!
  go to (150,137,137,140,140),NP
  137 DO 139 I=LSTART,LSTOP
     DO 138 J=MSTART,MSTOP
        F(I,J,2) = F(I,J,2)-C3*F(I,J,1)
  138    CONTINUE
  139 CONTINUE
  go to 143
  140 DO 142 I=LSTART,LSTOP
     DO 141 J=MSTART,MSTOP
        F(I,J,1) = F(I,J,1)+TWBYDZ*BDZS(I,J)
  141    CONTINUE
  142 CONTINUE
  143 go to (150,144,147,147,144),NP
  144 DO 146 I=LSTART,LSTOP
     DO 145 J=MSTART,MSTOP
        F(I,J,N) = F(I,J,N)-C3*F(I,J,NP1)
  145    CONTINUE
  146 CONTINUE
  go to 150
  147 DO 149 I=LSTART,LSTOP
     DO 148 J=MSTART,MSTOP
        F(I,J,NP1) = F(I,J,NP1)-TWBYDZ*BDZF(I,J)
  148    CONTINUE
  149 CONTINUE
!
!     DEFINE A,B,C COEFFICIENTS IN W-ARRAY.
!
  150 CONTINUE
  IWB = NUNK+1
  IWC = IWB+NUNK
  IWW = IWC+NUNK
  DO 151 K=1,NUNK
     I = IWC+K-1
     W(K) = C3
     W(I) = C3
     I = IWB+K-1
     W(I) = -2.*C3+ELMBDA
  151 CONTINUE
  go to (155,155,153,152,152),NP
  152 W(IWC) = 2.*C3
  153 go to (155,155,154,154,155),NP
  154 W(IWB-1) = 2.*C3
  155 CONTINUE
  PERTRB = 0.
!
!     FOR SINGULAR PROBLEMS ADJUST DATA TO INSURE A SOLUTION WILL EXIST.
!
  go to (156,172,172,156,172),LP
  156 go to (157,172,172,157,172),MP
  157 go to (158,172,172,158,172),NP
  158 if (ELMBDA) 172,160,159
  159 IERROR = 12
  go to 172
  160 CONTINUE
  MSTPM1 = MSTOP-1
  LSTPM1 = LSTOP-1
  NSTPM1 = NSTOP-1
  XLP = (2+LP)/3
  YLP = (2+MP)/3
  ZLP = (2+NP)/3
  S1 = 0.
  DO 164 K=2,NSTPM1
     DO 162 J=2,MSTPM1
        DO 161 I=2,LSTPM1
           S1 = S1+F(I,J,K)
  161       CONTINUE
        S1 = S1+(F(1,J,K)+F(LSTOP,J,K))/XLP
  162    CONTINUE
     S2 = 0.
     DO 163 I=2,LSTPM1
        S2 = S2+F(I,1,K)+F(I,MSTOP,K)
  163    CONTINUE
     S2 = (S2+(F(1,1,K)+F(1,MSTOP,K)+F(LSTOP,1,K)+F(LSTOP,MSTOP,K))/ &
                                                            XLP)/YLP
     S1 = S1+S2
  164 CONTINUE
  S = (F(1,1,1)+F(LSTOP,1,1)+F(1,1,NSTOP)+F(LSTOP,1,NSTOP)+ &
      F(1,MSTOP,1)+F(LSTOP,MSTOP,1)+F(1,MSTOP,NSTOP)+ &
                                     F(LSTOP,MSTOP,NSTOP))/(XLP*YLP)
  DO 166 J=2,MSTPM1
     DO 165 I=2,LSTPM1
        S = S+F(I,J,1)+F(I,J,NSTOP)
  165    CONTINUE
  166 CONTINUE
  S2 = 0.
  DO 167 I=2,LSTPM1
     S2 = S2+F(I,1,1)+F(I,1,NSTOP)+F(I,MSTOP,1)+F(I,MSTOP,NSTOP)
  167 CONTINUE
  S = S2/YLP+S
  S2 = 0.
  DO 168 J=2,MSTPM1
     S2 = S2+F(1,J,1)+F(1,J,NSTOP)+F(LSTOP,J,1)+F(LSTOP,J,NSTOP)
  168 CONTINUE
  S = S2/XLP+S
  PERTRB = (S/ZLP+S1)/((LUNK+1.-XLP)*(MUNK+1.-YLP)* &
                                                (NUNK+1.-ZLP))
  DO 171 I=1,LUNK
     DO 170 J=1,MUNK
        DO 169 K=1,NUNK
           F(I,J,K) = F(I,J,K)-PERTRB
  169       CONTINUE
  170    CONTINUE
  171 CONTINUE
  172 CONTINUE
  NPEROD = 0
  if (NBDCND  ==  0) go to 173
  NPEROD = 1
  W(1) = 0.
  W(IWW-1) = 0.
  173 CONTINUE
  call POIS3D (LBDCND,LUNK,C1,MBDCND,MUNK,C2,NPEROD,NUNK,W,W(IWB), &
               W(IWC),LDIMF,MDIMF,F(LSTART,MSTART,NSTART),IR,W(IWW))
!
!     FILL IN SIDES FOR PERIODIC BOUNDARY CONDITIONS.
!
  if (LP  /=  1) go to 180
  if (MP  /=  1) go to 175
  DO 174 K=NSTART,NSTOP
     F(1,MP1,K) = F(1,1,K)
  174 CONTINUE
  MSTOP = MP1
  175 if (NP  /=  1) go to 177
  DO 176 J=MSTART,MSTOP
     F(1,J,NP1) = F(1,J,1)
  176 CONTINUE
  NSTOP = NP1
  177 DO 179 J=MSTART,MSTOP
     DO 178 K=NSTART,NSTOP
        F(LP1,J,K) = F(1,J,K)
  178    CONTINUE
  179 CONTINUE
  180 CONTINUE
  if (MP  /=  1) go to 185
  if (NP  /=  1) go to 182
  DO 181 I=LSTART,LSTOP
     F(I,1,NP1) = F(I,1,1)
  181 CONTINUE
  NSTOP = NP1
  182 DO 184 I=LSTART,LSTOP
     DO 183 K=NSTART,NSTOP
        F(I,MP1,K) = F(I,1,K)
  183    CONTINUE
  184 CONTINUE
  185 CONTINUE
  if (NP  /=  1) go to 188
  DO 187 I=LSTART,LSTOP
     DO 186 J=MSTART,MSTOP
        F(I,J,NP1) = F(I,J,1)
  186    CONTINUE
  187 CONTINUE
  188 CONTINUE
  return
end
