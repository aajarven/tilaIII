subroutine CBLKTR (IFLG, NP, N, AN, BN, CN, MP, M, AM, BM, CM, &
     IDIMY, Y, IERROR, W)
!
!! CBLKTR solves a block tridiagonal system of linear equations ...
!            (usually resulting from the discretization of separable ...
!            two-dimensional elliptic equations).
!
!***LIBRARY   SLATEC (FISHPACK)
!***CATEGORY  I2B4B
!***TYPE      COMPLEX (BLKTRI-S, CBLKTR-C)
!***KEYWORDS  ELLIPTIC PDE, FISHPACK, TRIDIAGONAL LINEAR SYSTEM
!***AUTHOR  Adams, J., (NCAR)
!           Swarztrauber, P. N., (NCAR)
!           Sweet, R., (NCAR)
!***DESCRIPTION
!
!     Subroutine CBLKTR is a complex version of subroutine BLKTRI.
!     Both subroutines solve a system of linear equations of the form
!
!          AN(J)*X(I,J-1) + AM(I)*X(I-1,J) + (BN(J)+BM(I))*X(I,J)
!
!          + CN(J)*X(I,J+1) + CM(I)*X(I+1,J) = Y(I,J)
!
!               For I = 1,2,...,M  and  J = 1,2,...,N.
!
!     I+1 and I-1 are evaluated modulo M and J+1 and J-1 modulo N, i.e.,
!
!          X(I,0) = X(I,N),  X(I,N+1) = X(I,1),
!          X(0,J) = X(M,J),  X(M+1,J) = X(1,J).
!
!     These equations usually result from the discretization of
!     separable elliptic equations.  Boundary conditions may be
!     Dirichlet, Neumann, or periodic.
!
!
!     * * * * * * * * * *     On INPUT     * * * * * * * * * *
!
!     IFLG
!       = 0  Initialization only.  Certain quantities that depend on NP,
!            N, AN, BN, and CN are computed and stored in the work
!            array  W.
!       = 1  The quantities that were computed in the initialization are
!            used to obtain the solution X(I,J).
!
!       NOTE   A call with IFLG=0 takes approximately one half the time
!              time as a call with IFLG = 1.  However, the
!              initialization does not have to be repeated unless NP, N,
!              AN, BN, or CN change.
!
!     NP
!       = 0  If AN(1) and CN(N) are not zero, which corresponds to
!            periodic boundary conditions.
!       = 1  If AN(1) and CN(N) are zero.
!
!     N
!       The number of unknowns in the J-direction. N must be greater
!       than 4. The operation count is proportional to MNlog2(N), hence
!       N should be selected less than or equal to M.
!
!     AN,BN,CN
!       Real one-dimensional arrays of length N that specify the
!       coefficients in the linear equations given above.
!
!     MP
!       = 0  If AM(1) and CM(M) are not zero, which corresponds to
!            periodic boundary conditions.
!       = 1  If AM(1) = CM(M) = 0  .
!
!     M
!       The number of unknowns in the I-direction. M must be greater
!       than 4.
!
!     AM,BM,CM
!       Complex one-dimensional arrays of length M that specify the
!       coefficients in the linear equations given above.
!
!     IDIMY
!       The row (or first) dimension of the two-dimensional array Y as
!       it appears in the program calling BLKTRI.  This parameter is
!       used to specify the variable dimension of Y.  IDIMY must be at
!       least M.
!
!     Y
!       A complex two-dimensional array that specifies the values of
!       the right side of the linear system of equations given above.
!       Y must be dimensioned Y(IDIMY,N) with IDIMY  >=  M.
!
!     W
!       A one-dimensional array that must be provided by the user for
!       work space.
!             If NP=1 define K=INT(log2(N))+1 and set L=2**(K+1) then
!                     W must have dimension (K-2)*L+K+5+MAX(2N,12M)
!
!             If NP=0 define K=INT(log2(N-1))+1 and set L=2**(K+1) then
!                     W must have dimension (K-2)*L+K+5+2N+MAX(2N,12M)
!
!       **IMPORTANT** For purposes of checking, the required dimension
!                     of W is computed by BLKTRI and stored in W(1)
!                     in floating point format.
!
!     * * * * * * * * * *     On Output     * * * * * * * * * *
!
!     Y
!       Contains the solution X.
!
!     IERROR
!       An error flag that indicates invalid input parameters.  Except
!       for number zero, a solution is not attempted.
!
!       = 0  No error.
!       = 1  M is less than 5.
!       = 2  N is less than 5.
!       = 3  IDIMY is less than M.
!       = 4  BLKTRI failed while computing results that depend on the
!            coefficient arrays AN, BN, CN.  Check these arrays.
!       = 5  AN(J)*CN(J-1) is less than 0 for some J. Possible reasons
!            for this condition are
!            1. The arrays AN and CN are not correct.
!            2. Too large a grid spacing was used in the discretization
!               of the elliptic equation.
!            3. The linear equations resulted from a partial
!               differential equation which was not elliptic.
!
!     W
!       Contains intermediate values that must not be destroyed if
!       CBLKTR will be called again with IFLG=1.  W(1) contains the
!       number of locations required by W in floating point format.
!
! *Long Description:
!
!     * * * * * * *   Program Specifications    * * * * * * * * * * * *
!
!     Dimension of   AN(N),BN(N),CN(N),AM(M),BM(M),CM(M),Y(IDIMY,N)
!     Arguments      W(see argument list)
!
!     Latest         June 1979
!     Revision
!
!     Required       CBLKTR,CBLKT1,PROC,PROCP,CPROC,CPROCP,CCMPB,INXCA,
!     Subprograms    INXCB,INXCC,CPADD,PGSF,PPGSF,PPPSF,BCRH,TEVLC,
!                    R1MACH
!
!     Special        The algorithm may fail if ABS(BM(I)+BN(J)) is less
!     Conditions     than ABS(AM(I))+ABS(AN(J))+ABS(CM(I))+ABS(CN(J))
!                    for some I and J. The algorithm will also fail if
!                    AN(J)*CN(J-1) is less than zero for some J.
!                    See the description of the output parameter IERROR.
!
!     Common         CCBLK
!     Blocks
!
!     I/O            NONE
!
!     Precision      Single
!
!     Specialist     Paul Swarztrauber
!
!     Language       FORTRAN
!
!     History        CBLKTR is a complex version of BLKTRI (version 3)
!
!     Algorithm      Generalized Cyclic Reduction (see reference below)
!
!     Space
!     Required       CONTROL DATA 7600
!
!     Portability    American National Standards Institute FORTRAN.
!                    The machine accuracy is set using function R1MACH.
!
!     Required       NONE
!     Resident
!     Routines
!
!     References     Swarztrauber,P. and R. SWEET, 'Efficient Fortran
!                    Subprograms for the solution of elliptic equations'
!                    NCAR TN/IA-109, July, 1975, 138 PP.
!
!                    SWARZTRAUBER P. ,'A Direct Method for The Discrete
!                    Solution of Separable Elliptic Equations', SIAM
!                    J. Numer. Anal.,11(1974) PP. 1136-1150.
!
!     * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *
!
!***REFERENCES  P. N. Swarztrauber and R. Sweet, Efficient Fortran
!                 subprograms for the solution of elliptic equations,
!                 NCAR TN/IA-109, July 1975, 138 pp.
!               P. N. Swarztrauber, A direct method for the discrete
!                 solution of separable elliptic equations, SIAM Journal
!                 on Numerical Analysis 11, (1974), pp. 1136-1150.
!***ROUTINES CALLED  CBLKT1, CCMPB, CPROC, CPROCP, PROC, PROCP
!***COMMON BLOCKS    CCBLK
!***REVISION HISTORY  (YYMMDD)
!   801001  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890531  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  CBLKTR
!
  DIMENSION       AN(*)      ,BN(*)      ,CN(*)      ,AM(*)      , &
                  BM(*)      ,CM(*)      ,Y(IDIMY,*) ,W(*)
  EXTERNAL        PROC       ,PROCP      ,CPROC      ,CPROCP
  COMMON /CCBLK/  NPP        ,K          ,EPS        ,CNV        , &
                  NM         ,NCMPLX     ,IK
  COMPLEX         AM         ,BM         ,CM         ,Y
!***FIRST EXECUTABLE STATEMENT  CBLKTR
  NM = N
  M2 = M+M
  IERROR = 0
  if (M-5) 101,102,102
  101 IERROR = 1
  go to 119
  102 if (NM-3) 103,104,104
  103 IERROR = 2
  go to 119
  104 if (IDIMY-M) 105,106,106
  105 IERROR = 3
  go to 119
  106 NH = N
  NPP = NP
  if (NPP) 107,108,107
  107 NH = NH+1
  108 IK = 2
  K = 1
  109 IK = IK+IK
  K = K+1
  if (NH-IK) 110,110,109
  110 NL = IK
  IK = IK+IK
  NL = NL-1
  IWAH = (K-2)*IK+K+6
  if (NPP) 111,112,111
!
!     DIVIDE W INTO WORKING SUB ARRAYS
!
  111 IW1 = IWAH
  IWBH = IW1+NM
  W(1) = IW1-1+MAX(2*NM,12*M)
  go to 113
  112 IWBH = IWAH+NM+NM
  IW1 = IWBH
  W(1) = IW1-1+MAX(2*NM,12*M)
  NM = NM-1
!
! SUBROUTINE CCMPB COMPUTES THE ROOTS OF THE B POLYNOMIALS
!
  113 if (IERROR) 119,114,119
  114 IW2 = IW1+M2
  IW3 = IW2+M2
  IWD = IW3+M2
  IWW = IWD+M2
  IWU = IWW+M2
  if (IFLG) 116,115,116
  115 call CCMPB (NL,IERROR,AN,BN,CN,W(2),W(IWAH),W(IWBH))
  go to 119
  116 if (MP) 117,118,117
!
! SUBROUTINE CBLKT1 SOLVES THE LINEAR SYSTEM
!
  117 call CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2), &
               W(IW3),W(IWD),W(IWW),W(IWU),PROC,CPROC)
  go to 119
  118 call CBLKT1 (NL,AN,BN,CN,M,AM,BM,CM,IDIMY,Y,W(2),W(IW1),W(IW2), &
               W(IW3),W(IWD),W(IWW),W(IWU),PROCP,CPROCP)
  119 CONTINUE
  return
end
