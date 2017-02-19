subroutine DBVSUP (Y, NROWY, NCOMP, XPTS, NXPTS, A, NROWA, ALPHA, &
     NIC, B, NROWB, BETA, NFC, IGOFX, RE, AE, IFLAG, WORK, NDW, &
     IWORK, NDIW, NEQIVP)
!
!! DBVSUP solves a linear two-point boundary value problem using ...
!  superposition coupled with an orthonormalization procedure
!  and a variable-step integration scheme.
!
!***LIBRARY   SLATEC
!***CATEGORY  I1B1
!***TYPE      DOUBLE PRECISION (BVSUP-S, DBVSUP-D)
!***KEYWORDS  ORTHONORMALIZATION, SHOOTING,
!             TWO-POINT BOUNDARY VALUE PROBLEM
!***AUTHOR  Scott, M. R., (SNLA)
!           Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!
!     Subroutine DBVSUP solves a linear two-point boundary-value problem
!     of the form
!                        DY/DX = MATRIX(X,U)*Y(X) + G(X,U)
!                A*Y(XINITIAL) = ALPHA ,  B*Y(XFINAL) = BETA
!
!     coupled with the solution of the initial value problem
!
!                        DU/DX = F(X,U)
!                      U(XINITIAL) = ETA
!
! **********************************************************************
!     ABSTRACT
!        The method of solution uses superposition coupled with an
!     orthonormalization procedure and a variable-step integration
!     scheme.  Each time the superposition solutions start to
!     lose their numerical linear independence, the vectors are
!     reorthonormalized before integration proceeds.  The underlying
!     principle of the algorithm is then to piece together the
!     intermediate (orthogonalized) solutions, defined on the various
!     subintervals, to obtain the desired solutions.
!
! **********************************************************************
!     INPUT to DBVSUP
! **********************************************************************
!
!     NROWY = actual row dimension of Y in calling program.
!             NROWY must be  >=  NCOMP
!
!     NCOMP = number of components per solution vector.
!             NCOMP is equal to number of original differential
!             equations.  NCOMP = NIC + NFC.
!
!     XPTS = desired output points for solution. They must be monotonic.
!            XINITIAL = XPTS(1)
!            XFINAL = XPTS(NXPTS)
!
!     NXPTS = number of output points.
!
!     A(NROWA,NCOMP) = boundary condition matrix at XINITIAL
!                      must be contained in (NIC,NCOMP) sub-matrix.
!
!     NROWA = actual row dimension of A in calling program,
!             NROWA must be  >=  NIC.
!
!     ALPHA(NIC+NEQIVP) = boundary conditions at XINITIAL.
!                         If NEQIVP  >  0 (see below), the boundary
!                         conditions at XINITIAL for the initial value
!                         equations must be stored starting in
!                         position (NIC + 1) of ALPHA.
!                         Thus,  ALPHA(NIC+K) = ETA(K).
!
!     NIC = number of boundary conditions at XINITIAL.
!
!     B(NROWB,NCOMP) = boundary condition matrix at XFINAL.
!                      Must be contained in (NFC,NCOMP) sub-matrix.
!
!     NROWB = actual row dimension of B in calling program,
!             NROWB must be  >=  NFC.
!
!     BETA(NFC) = boundary conditions at XFINAL.
!
!     NFC = number of boundary conditions at XFINAL.
!
!     IGOFX =0 -- The inhomogeneous term G(X) is identically zero.
!           =1 -- The inhomogeneous term G(X) is not identically zero.
!                 (if IGOFX=1, then Subroutine DGVEC (or DUVEC) must be
!                  supplied).
!
!     RE = relative error tolerance used by the integrator.
!          (see one of the integrators)
!
!     AE = absolute error tolerance used by the integrator.
!          (see one of the integrators)
! **NOTE-  RE and AE should not both be zero.
!
!     IFLAG = a status parameter used principally for output.
!             However, for efficient solution of problems which
!             are originally defined as COMPLEX*16 valued (but
!             converted to double precision systems to use this code),
!             the user must set IFLAG=13 on input. See the comment
!             below for more information on solving such problems.
!
!     WORK(NDW) = floating point array used for internal storage.
!
!     NDW = actual dimension of work array allocated by user.
!           An estimate for NDW can be computed from the following
!            NDW = 130 + NCOMP**2 * (6 + NXPTS/2 + expected number of
!                                           orthonormalizations/8)
!           For the disk or tape storage mode,
!            NDW = 6 * NCOMP**2 + 10 * NCOMP + 130
!  However, when the ADAMS integrator is to be used, the estimates are
!            NDW = 130 + NCOMP**2 * (13 + NXPTS/2 + expected number of
!                                           orthonormalizations/8)
!    and     NDW = 13 * NCOMP**2 + 22 * NCOMP + 130   , respectively.
!
!     IWORK(NDIW) = integer array used for internal storage.
!
!     NDIW = actual dimension of IWORK array allocated by user.
!            An estimate for NDIW can be computed from the following
!            NDIW = 68 + NCOMP * (1 + expected number of
!                                            orthonormalizations)
! **NOTE --  the amount of storage required is problem dependent and may
!            be difficult to predict in advance.  Experience has shown
!            that for most problems 20 or fewer orthonormalizations
!            should suffice. If the problem cannot be completed with the
!            allotted storage, then a message will be printed which
!            estimates the amount of storage necessary. In any case, the
!            user can examine the IWORK array for the actual storage
!            requirements, as described in the output information below.
!
!     NEQIVP = number of auxiliary initial value equations being added
!              to the boundary value problem.
! **NOTE -- Occasionally the coefficients  matrix  and/or  G  may be
!           functions which depend on the independent variable  X  and
!           on  U, the solution of an auxiliary initial value problem.
!           In order to avoid the difficulties associated with
!           interpolation, the auxiliary equations may be solved
!           simultaneously with the given boundary value problem.
!           This initial value problem may be linear or nonlinear.
!                 See SAND77-1328 for an example.
!
!
!     The user must supply subroutines DFMAT, DGVEC, DUIVP and DUVEC,
!     when needed (they must be so named), to evaluate the derivatives
!     as follows
!
!        A. DFMAT must be supplied.
!
!              SUBROUTINE DFMAT(X,Y,YP)
!              X = independent variable (input to DFMAT)
!              Y = dependent variable vector (input to DFMAT)
!              YP = DY/DX = derivative vector (output from DFMAT)
!
!            Compute the derivatives for the homogeneous problem
!              YP(I) = DY(I)/DX = MATRIX(X) * Y(I)  , I = 1,...,NCOMP
!
!            When (NEQIVP  >  0) and  matrix  is dependent on  U  as
!            well as on  X, the following common statement must be
!            included in DFMAT
!                    COMMON /DMLIVP/ NOFST
!            for convenience, the  U  vector is stored at the bottom
!            of the  Y  array.  Thus, during any call to DFMAT,
!            U(I) is referenced by  Y(NOFST + I).
!
!
!            Subroutine DBVDER calls DFMAT NFC times to evaluate the
!            homogeneous equations and, if necessary, it calls DFMAT
!            once in evaluating the particular solution. since X remains
!            unchanged in this sequence of calls it is possible to
!            realize considerable computational savings for complicated
!            and expensive evaluations of the matrix entries. To do this
!            the user merely passes a variable, say XS, via common where
!            XS is defined in the main program to be any value except
!            the initial X. Then the non-constant elements of matrix(x)
!            appearing in the differential equations need only be
!            computed if X is unequal to XS, whereupon XS is reset to X.
!
!
!        B. If  NEQIVP  >  0 ,  DUIVP must also be supplied.
!
!              SUBROUTINE DUIVP(X,U,UP)
!              X = independent variable (input to DUIVP)
!              U = dependent variable vector (input to DUIVP)
!              UP = DU/DX = derivative vector (output from DUIVP)
!
!            Compute the derivatives for the auxiliary initial value eqs
!              UP(I) = DU(I)/DX, I = 1,...,NEQIVP.
!
!            Subroutine DBVDER calls DUIVP once to evaluate the
!            derivatives for the auxiliary initial value equations.
!
!
!        C. If  NEQIVP = 0  and  IGOFX = 1 ,  DGVEC must be supplied.
!
!              SUBROUTINE DGVEC(X,G)
!              X = independent variable (input to DGVEC)
!              G = vector of inhomogeneous terms G(X) (output from
!              DGVEC)
!
!            Compute the inhomogeneous terms G(X)
!                G(I) = G(X) values for I = 1,...,NCOMP.
!
!            Subroutine DBVDER calls DGVEC in evaluating the particular
!            solution provided G(X) is not identically zero. Thus, when
!            IGOFX=0, the user need not write a DGVEC subroutine. Also,
!            the user does not have to bother with the computational
!            savings scheme for DGVEC as this is automatically achieved
!            via the DBVDER subroutine.
!
!
!        D. If  NEQIVP  >  0  and  IGOFX = 1 ,  DUVEC must be supplied.
!
!             SUBROUTINE DUVEC(X,U,G)
!             X = independent variable (input to DUVEC)
!             U = dependent variable vector from the auxiliary initial
!                 value problem    (input to DUVEC)
!             G = array of inhomogeneous terms G(X,U)(output from DUVEC)
!
!            Compute the inhomogeneous terms G(X,U)
!                G(I) = G(X,U) values for I = 1,...,NCOMP.
!
!            Subroutine DBVDER calls DUVEC in evaluating the particular
!            solution provided G(X,U) is not identically zero.  Thus,
!            when IGOFX=0, the user need not write a DUVEC subroutine.
!
!
!
!     The following is optional input to DBVSUP to give user more
!     flexibility in use of code.  See SAND75-0198, SAND77-1328,
!     SAND77-1690, SAND78-0522, and SAND78-1501 for more information.
!
! ****CAUTION -- The user must zero out IWORK(1),...,IWORK(15)
!                prior to calling DBVSUP. These locations define
!                optional input and must be zero unless set to special
!                values by the user as described below.
!
!     IWORK(1) -- number of orthonormalization points.
!                 A value need be set only if IWORK(11) = 1
!
!     IWORK(9) -- integrator and orthonormalization parameter
!                 (default value is 1)
!                 1 = RUNGE-KUTTA-FEHLBERG code using GRAM-SCHMIDT test.
!                 2 = ADAMS code using GRAM-SCHMIDT test.
!
!     IWORK(11) -- orthonormalization points parameter
!                  (default value is 0)
!                  0 - orthonormalization points not pre-assigned.
!                  1 - orthonormalization points pre-assigned in
!                      the first IWORK(1) positions of work.
!
!     IWORK(12) -- storage parameter
!                  (default value is 0)
!                  0 - all storage in core.
!                  LUN - homogeneous and inhomogeneous solutions at
!                      output points and orthonormalization information
!                      are stored on disk.  The logical unit number to
!                      be used for disk I/O (NTAPE) is set to IWORK(12).
!
!     WORK(1),... -- pre-assigned orthonormalization points, stored
!                    monotonically, corresponding to the direction
!                    of integration.
!
!
!
!                 ******************************************************
!                 *** COMPLEX*16 VALUED PROBLEM ***
!                 ******************************************************
! **NOTE***
!       Suppose the original boundary value problem is NC equations
!     of the form
!                   DW/DX = MAT(X,U)*W(X) + H(X,U)
!                 R*W(XINITIAL)=GAMMA , S*W(XFINAL)=DELTA
!     where all variables are COMPLEX*16 valued. The DBVSUP code can be
!     used by converting to a double precision system of size 2*NC. To
!     solve the larger dimensioned problem efficiently, the user must
!     initialize IFLAG=13 on input and order the vector components
!     according to Y(1)=DOUBLE PRECISION(W(1)),...,Y(NC)=DOUBLE
!     PRECISION(W(NC)),Y(NC+1)=IMAG(W(1)),...., Y(2*NC)=IMAG(W(NC)).
!     Then define
!                        ...............................................
!                        . DOUBLE PRECISION(MAT)    -IMAG(MAT) .
!            MATRIX  =   .                         .
!                        . IMAG(MAT)     DOUBLE PRECISION(MAT) .
!                        ...............................................
!
!     The matrices A,B and vectors G,ALPHA,BETA must be defined
!     similarly. Further details can be found in SAND78-1501.
!
!
! **********************************************************************
!     OUTPUT from DBVSUP
! **********************************************************************
!
!     Y(NROWY,NXPTS) = solution at specified output points.
!
!     IFLAG Output Values
!            =-5 algorithm ,for obtaining starting vectors for the
!                special COMPLEX*16 problem structure, was unable to
!                obtain the initial vectors satisfying the necessary
!                independence criteria.
!            =-4 rank of boundary condition matrix A is less than NIC,
!                as determined by DLSSUD.
!            =-2 invalid input parameters.
!            =-1 insufficient number of storage locations allocated for
!                WORK or IWORK.
!
!            =0 indicates successful solution.
!
!            =1 a computed solution is returned but uniqueness of the
!               solution of the boundary-value problem is questionable.
!               For an eigenvalue problem, this should be treated as a
!               successful execution since this is the expected mode
!               of return.
!            =2 a computed solution is returned but the existence of the
!               solution to the boundary-value problem is questionable.
!            =3 a nontrivial solution approximation is returned although
!               the boundary condition matrix B*Y(XFINAL) is found to be
!               nonsingular (to the desired accuracy level) while the
!               right hand side vector is zero. To eliminate this type
!               of return, the accuracy of the eigenvalue parameter
!               must be improved.
!            ***NOTE-We attempt to diagnose the correct problem behavior
!               and report possible difficulties by the appropriate
!               error flag.  However, the user should probably resolve
!               the problem using smaller error tolerances and/or
!               perturbations in the boundary conditions or other
!               parameters. This will often reveal the correct
!               interpretation for the problem posed.
!
!            =13 maximum number of orthonormalizations attained before
!                reaching XFINAL.
!            =20-flag from integrator (DDERKF or DDEABM) values can
!                range from 21 to 25.
!            =30 solution vectors form a dependent set.
!
!     WORK(1),...,WORK(IWORK(1)) = orthonormalization points
!                                  determined by DBVPOR.
!
!     IWORK(1) = number of orthonormalizations performed by DBVPOR.
!
!     IWORK(2) = maximum number of orthonormalizations allowed as
!                calculated from storage allocated by user.
!
!     IWORK(3),IWORK(4),IWORK(5),IWORK(6)   give information about
!                actual storage requirements for WORK and IWORK
!                arrays.  In particular,
!                       required storage for  work array is
!        IWORK(3) + IWORK(4)*(expected number of orthonormalizations)
!
!                       required storage for IWORK array is
!        IWORK(5) + IWORK(6)*(expected number of orthonormalizations)
!
!     IWORK(8) = final value of exponent parameter used in tolerance
!                test for orthonormalization.
!
!     IWORK(16) = number of independent vectors returned from DMGSBV.
!                It is only of interest when IFLAG=30 is obtained.
!
!     IWORK(17) = numerically estimated rank of the boundary
!                 condition matrix defined from B*Y(XFINAL)
!
! **********************************************************************
!
!     Necessary machine constants are defined in the Function
!     Routine D1MACH. The user must make sure that the values
!     set in D1MACH are relevant to the computer being used.
!
! **********************************************************************
! **********************************************************************
!
!***REFERENCES  M. R. Scott and H. A. Watts, SUPORT - a computer code
!                 for two-point boundary-value problems via
!                 orthonormalization, SIAM Journal of Numerical
!                 Analysis 14, (1977), pp. 40-70.
!               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications
!                 of SUPORT, a linear boundary value problem solver
!                 Part I - pre-assigning orthonormalization points,
!                 auxiliary initial value problem, disk or tape storage,
!                 Report SAND77-1328, Sandia Laboratories, Albuquerque,
!                 New Mexico, 1977.
!               B. L. Darlow, M. R. Scott and H. A. Watts, Modifications
!                 of SUPORT, a linear boundary value problem solver
!                 Part II - inclusion of an Adams integrator, Report
!                 SAND77-1690, Sandia Laboratories, Albuquerque,
!                 New Mexico, 1977.
!               M. E. Lord and H. A. Watts, Modifications of SUPORT,
!                 a linear boundary value problem solver Part III -
!                 orthonormalization improvements, Report SAND78-0522,
!                 Sandia Laboratories, Albuquerque, New Mexico, 1978.
!               H. A. Watts, M. R. Scott and M. E. Lord, Computational
!                 solution of complex*16 valued boundary problems,
!                 Report SAND78-1501, Sandia Laboratories,
!                 Albuquerque, New Mexico, 1978.
!***ROUTINES CALLED  DEXBVP, DMACON, XERMSG
!***COMMON BLOCKS    DML15T, DML17B, DML18J, DML5MC, DML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   890921  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900510  Convert XERRWV calls to XERMSG calls, remove some extraneous
!           comments.  (RWC)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DBVSUP
! **********************************************************************
!
  INTEGER ICOCO, IFLAG, IGOFX, IGOFXD, INDPVT, INFO, INHOMO, INTEG, &
       IS, ISTKOP, IVP, IWORK(*), J, K, K1, K10, K11, K2, &
       K3, K4, K5, K6, K7, K8, K9, KKKCOE, KKKCOF, KKKG, KKKINT, &
       KKKS, KKKSTO, KKKSUD, KKKSVC, KKKU, KKKV, KKKWS, KKKYHP, &
       KKKZPW, KNSWOT, KOP, KPTS, L1, L2, LLLCOF, LLLINT, LLLIP, &
       LLLIWS, LLLSUD, LLLSVC, LOTJP, LPAR, MNSWOT, &
       MXNON, MXNONI, MXNONR, NCOMP, NCOMPD, NDEQ, NDISK, NDIW, &
       NDW, NEEDIW, NEEDW, NEQ, NEQIVD, NEQIVP, NFC, NFCC, &
       NFCD, NIC, NICD, NITEMP, NON, NOPG, NPS, NROWA, NROWB, &
       NROWY, NRTEMP, NSWOT, NTAPE, NTP, NUMORT, NXPTS, NXPTSD, &
       NXPTSM
  DOUBLE PRECISION A(NROWA,*), AE, AED, ALPHA(*), &
       B(NROWB,*), BETA(*), C, EPS, FOURU, PWCND, PX, RE, &
       RED, SQOVFL, SRU, TND, TOL, TWOU, URO, WORK(NDW), X, XBEG, &
       XEND, XOP, XOT, XPTS(*), XSAV, Y(NROWY,*)
  CHARACTER*8 XERN1, XERN2, XERN3, XERN4
!
!     ******************************************************************
!         THE COMMON BLOCK BELOW IS USED TO COMMUNICATE WITH SUBROUTINE
!         DBVDER.  THE USER SHOULD NOT ALTER OR USE THIS COMMON BLOCK IN
!         THE CALLING PROGRAM.
!
  COMMON /DML8SZ/ C,XSAV,IGOFXD,INHOMO,IVP,NCOMPD,NFCD
!
!     ******************************************************************
!         THESE COMMON BLOCKS AID IN REDUCING THE NUMBER OF SUBROUTINE
!         ARGUMENTS PREVALENT IN THIS MODULAR STRUCTURE
!
  COMMON /DML18J/ AED,RED,TOL,NXPTSD,NICD,NOPG,MXNON,NDISK,NTAPE, &
                  NEQ,INDPVT,INTEG,NPS,NTP,NEQIVD,NUMORT,NFCC, &
                  ICOCO
  COMMON /DML17B/ KKKZPW,NEEDW,NEEDIW,K1,K2,K3,K4,K5,K6,K7,K8,K9, &
                  K10,K11,L1,L2,KKKINT,LLLINT
!
!     ******************************************************************
!         THIS COMMON BLOCK IS USED IN SUBROUTINES DBVSUP,DBVPOR,DRKFAB,
!         DREORT, AND DSTWAY. IT CONTAINS INFORMATION NECESSARY
!         FOR THE ORTHONORMALIZATION TESTING PROCEDURE AND A BACKUP
!         RESTARTING CAPABILITY.
!
  COMMON /DML15T/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
!
!     ******************************************************************
!         THIS COMMON BLOCK CONTAINS THE MACHINE DEPENDENT PARAMETERS
!         USED BY THE CODE
!
  COMMON /DML5MC/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
!
!      *****************************************************************
!          SET UP MACHINE DEPENDENT CONSTANTS.
!
!***FIRST EXECUTABLE STATEMENT  DBVSUP
                    call DMACON
!
!                       ************************************************
!                           TEST FOR INVALID INPUT
!
                    if (NROWY  <  NCOMP) go to 80
                    if (NCOMP  /=  NIC + NFC) go to 80
                    if (NXPTS  <  2) go to 80
                    if (NIC  <=  0) go to 80
                    if (NROWA  <  NIC) go to 80
                    if (NFC  <=  0) go to 80
                    if (NROWB  <  NFC) go to 80
                    if (IGOFX  <  0 .OR. IGOFX  >  1) go to 80
                    if (RE  <  0.0D0) go to 80
                    if (AE  <  0.0D0) go to 80
                    if (RE  ==  0.0D0 .AND. AE  ==  0.0D0) go to 80
!                          BEGIN BLOCK PERMITTING ...EXITS TO 70
                          IS = 1
                          if (XPTS(NXPTS)  <  XPTS(1)) IS = 2
                          NXPTSM = NXPTS - 1
                          DO 30 K = 1, NXPTSM
                             if (IS  ==  2) go to 10
!                          .........EXIT
                                if (XPTS(K+1)  <=  XPTS(K)) go to 70
                             go to 20
   10                            CONTINUE
!                          .........EXIT
                                if (XPTS(K)  <=  XPTS(K+1)) go to 70
   20                            CONTINUE
   30                         CONTINUE
!
!                             ******************************************
!                                 CHECK FOR DISK STORAGE
!
                          KPTS = NXPTS
                          NDISK = 0
                          if (IWORK(12)  ==  0) go to 40
                             NTAPE = IWORK(12)
                             KPTS = 1
                             NDISK = 1
   40                         CONTINUE
!
!                             ******************************************
!                                 SET INTEG PARAMETER ACCORDING TO
!                                 CHOICE OF INTEGRATOR.
!
                          INTEG = 1
                          if (IWORK(9)  ==  2) INTEG = 2
!
!                             ******************************************
!                                 COMPUTE INHOMO
!
!                 ............EXIT
                          if (IGOFX  ==  1) go to 100
                          DO 50 J = 1, NIC
!                 ...............EXIT
                             if (ALPHA(J)  /=  0.0D0) go to 100
   50                         CONTINUE
                          DO 60 J = 1, NFC
!                    ............EXIT
                             if (BETA(J)  /=  0.0D0) go to 90
   60                         CONTINUE
                          INHOMO = 3
!              ...............EXIT
                          go to 110
   70                      CONTINUE
   80                   CONTINUE
                    IFLAG = -2
!     ..................EXIT
                    go to 220
   90                CONTINUE
                 INHOMO = 2
!              ......EXIT
                 go to 110
  100             CONTINUE
              INHOMO = 1
  110          CONTINUE
!
!              *********************************************************
!                  TO TAKE ADVANTAGE OF THE SPECIAL STRUCTURE WHEN
!                  SOLVING A COMPLEX*16 VALUED PROBLEM,WE INTRODUCE
!                  NFCC=NFC WHILE CHANGING THE INTERNAL VALUE OF NFC
!
           NFCC = NFC
           if (IFLAG  ==  13) NFC = NFC/2
!
!              *********************************************************
!                  DETERMINE NECESSARY STORAGE REQUIREMENTS
!
!              FOR BASIC ARRAYS IN DBVPOR
           KKKYHP = NCOMP*(NFC + 1) + NEQIVP
           KKKU = NCOMP*NFC*KPTS
           KKKV = NCOMP*KPTS
           KKKCOE = NFCC
           KKKS = NFC + 1
           KKKSTO = NCOMP*(NFC + 1) + NEQIVP + 1
           KKKG = NCOMP
!
!              FOR ORTHONORMALIZATION RELATED MATTERS
           NTP = (NFCC*(NFCC + 1))/2
           KKKZPW = 1 + NTP + NFCC
           LLLIP = NFCC
!
!              FOR ADDITIONAL REQUIRED WORK SPACE
!                (DLSSUD)
           KKKSUD = 4*NIC + (NROWA + 1)*NCOMP
           LLLSUD = NIC
!              (DVECS)
           KKKSVC = 1 + 4*NFCC + 2*NFCC**2
           LLLSVC = 2*NFCC
!
           NDEQ = NCOMP*NFC + NEQIVP
           if (INHOMO  ==  1) NDEQ = NDEQ + NCOMP
           go to (120,130), INTEG
!              (DDERKF)
  120          CONTINUE
              KKKINT = 33 + 7*NDEQ
              LLLINT = 34
           go to 140
!              (DDEABM)
  130          CONTINUE
              KKKINT = 130 + 21*NDEQ
              LLLINT = 51
  140          CONTINUE
!
!              (COEF)
           KKKCOF = 5*NFCC + NFCC**2
           LLLCOF = 3 + NFCC
!
           KKKWS = MAX(KKKSUD,KKKSVC,KKKINT,KKKCOF)
           LLLIWS = MAX(LLLSUD,LLLSVC,LLLINT,LLLCOF)
!
           NEEDW = KKKYHP + KKKU + KKKV + KKKCOE + KKKS + KKKSTO &
                   + KKKG + KKKZPW + KKKWS
           NEEDIW = 17 + LLLIP + LLLIWS
!              *********************************************************
!                  COMPUTE THE NUMBER OF POSSIBLE ORTHONORMALIZATIONS
!                  WITH THE ALLOTTED STORAGE
!
           IWORK(3) = NEEDW
           IWORK(4) = KKKZPW
           IWORK(5) = NEEDIW
           IWORK(6) = LLLIP
           NRTEMP = NDW - NEEDW
           NITEMP = NDIW - NEEDIW
!           ...EXIT
           if (NRTEMP  <  0) go to 180
!           ...EXIT
           if (NITEMP  <  0) go to 180
!
           if (NDISK  ==  0) go to 150
              NON = 0
              MXNON = NRTEMP
           go to 160
  150          CONTINUE
!
              MXNONR = NRTEMP/KKKZPW
              MXNONI = NITEMP/LLLIP
              MXNON = MIN(MXNONR,MXNONI)
              NON = MXNON
  160          CONTINUE
!
           IWORK(2) = MXNON
!
!              *********************************************************
!                  CHECK FOR PRE-ASSIGNED ORTHONORMALIZATION POINTS
!
           NOPG = 0
!        ......EXIT
           if (IWORK(11)  /=  1) go to 210
           if (MXNON  <  IWORK(1)) go to 170
              NOPG = 1
              MXNON = IWORK(1)
              WORK(MXNON+1) = 2.0D0*XPTS(NXPTS) - XPTS(1)
!        .........EXIT
              go to 210
  170          CONTINUE
  180       CONTINUE
!
        IFLAG = -1
  if (NDISK  /=  1) THEN
     WRITE (XERN1, '(I8)') NEEDW
     WRITE (XERN2, '(I8)') KKKZPW
     WRITE (XERN3, '(I8)') NEEDIW
     WRITE (XERN4, '(I8)') LLLIP
     call XERMSG ('SLATEC', 'DBVSUP', &
        'REQUIRED STORAGE FOR WORK ARRAY IS '  // XERN1 // ' + ' // &
        XERN2 // '*(EXPECTED NUMBER OF ORTHONORMALIZATIONS) $$'  // &
        'REQUIRED STORAGE FOR IWORK ARRAY IS ' // XERN3 // ' + ' // &
        XERN4 // '*(EXPECTED NUMBER OF ORTHONORMALIZATIONS)', 1, 0)
  ELSE
     WRITE (XERN1, '(I8)') NEEDW
     WRITE (XERN2, '(I8)') NEEDIW
     call XERMSG ('SLATEC', 'DBVSUP', &
        'REQUIRED STORAGE FOR WORK ARRAY IS '  // XERN1 // &
        ' + NUMBER OF ORTHONOMALIZATIONS. $$'  // &
        'REQUIRED STORAGE FOR IWORK ARRAY IS ' // XERN2, 1, 0)
  end if
  return
!
!        ***************************************************************
!            ALLOCATE STORAGE FROM WORK AND IWORK ARRAYS
!
!         (Z)
  210    K1 = 1 + (MXNON + 1)
!        (P)
     K2 = K1 + NTP*(NON + 1)
!        (W)
     K3 = K2 + NFCC*(NON + 1)
!        (YHP)
     K4 = K3 + KKKYHP
!        (U)
     K5 = K4 + KKKU
!        (V)
     K6 = K5 + KKKV
!        (COEF)
     K7 = K6 + KKKCOE
!        (S)
     K8 = K7 + KKKS
!        (STOWA)
     K9 = K8 + KKKSTO
!        (G)
     K10 = K9 + KKKG
     K11 = K10 + KKKWS
!                  REQUIRED ADDITIONAL DOUBLE PRECISION WORK SPACE
!                  STARTS AT WORK(K10) AND EXTENDS TO WORK(K11-1)
!
!           FIRST 17 LOCATIONS OF IWORK ARE USED FOR OPTIONAL
!           INPUT AND OUTPUT ITEMS
!        (IP)
     L1 = 18 + NFCC*(NON + 1)
     L2 = L1 + LLLIWS
!                   REQUIRED INTEGER WORK SPACE STARTS AT IWORK(L1)
!                   AND EXTENDS TO IWORK(L2-1)
!
!        ***************************************************************
!            SET INDICATOR FOR NORMALIZATION OF PARTICULAR SOLUTION
!
     NPS = 0
     if (IWORK(10)  ==  1) NPS = 1
!
!        ***************************************************************
!            SET PIVOTING PARAMETER
!
     INDPVT = 0
     if (IWORK(15)  ==  1) INDPVT = 1
!
!        ***************************************************************
!            SET OTHER COMMON BLOCK PARAMETERS
!
     NFCD = NFC
     NCOMPD = NCOMP
     IGOFXD = IGOFX
     NXPTSD = NXPTS
     NICD = NIC
     RED = RE
     AED = AE
     NEQIVD = NEQIVP
     MNSWOT = 20
     if (IWORK(13)  ==  -1) MNSWOT = MAX(1,IWORK(14))
     XBEG = XPTS(1)
     XEND = XPTS(NXPTS)
     XSAV = XEND
     ICOCO = 1
     if (INHOMO  ==  3 .AND. NOPG  ==  1) WORK(MXNON+1) = XEND
!
!        ***************************************************************
!
     call DEXBVP(Y,NROWY,XPTS,A,NROWA,ALPHA,B,NROWB,BETA,IFLAG,WORK, &
                 IWORK)
     NFC = NFCC
     IWORK(17) = IWORK(L1)
  220 CONTINUE
  return
end
