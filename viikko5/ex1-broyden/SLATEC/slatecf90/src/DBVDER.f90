subroutine DBVDER (X, Y, YP, G, IPAR)
!
!! DBVDER is subsidiary to DBVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (BVDER-S, DBVDER-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
! **********************************************************************
!     NFC = Number of base solution vectors
!
!     NCOMP = Number of components per solution vector
!
!              1 -- Nonzero particular solution
!     INHOMO =
!              2 or 3 -- Zero particular solution
!
!             0 -- Inhomogeneous vector term G(X) identically zero
!     IGOFX =
!             1 -- Inhomogeneous vector term G(X) not identically zero
!
!     G = Inhomogeneous vector term G(X)
!
!     XSAV = Previous value of X
!
!     C = Normalization factor for the particular solution
!
!           0   ( if  NEQIVP = 0 )
!     IVP =
!           Number of differential equations integrated due to
!           the original boundary value problem   ( if  NEQIVP  >  0 )
!
!     NOFST - For problems with auxiliary initial value equations,
!             NOFST communicates to the routine DFMAT how to access
!             the dependent variables corresponding to this initial
!             value problem.  For example, during any call to DFMAT,
!             the first dependent variable for the initial value
!             problem is in position  Y(NOFST + 1).
!             See example in SAND77-1328.
! **********************************************************************
!
!***SEE ALSO  DBVSUP
!***ROUTINES CALLED  (NONE)
!***COMMON BLOCKS    DML8SZ, DMLIVP
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910701  Corrected ROUTINES CALLED section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!   920618  Minor restructuring of code.  (RWC, WRB)
!***END PROLOGUE  DBVDER
  INTEGER IGOFX, INHOMO, IPAR, IVP, J, K, L, NA, NCOMP, NFC, NOFST
  DOUBLE PRECISION C, G(*), X, XSAV, Y(*), YP(*)
!
! **********************************************************************
!
  COMMON /DML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
!
! **********************************************************************
!     The COMMON block below is used to communicate with the user
!     supplied subroutine DFMAT.  The user should not alter this
!     COMMON block.
!
  COMMON /DMLIVP/ NOFST
! **********************************************************************
!
!***FIRST EXECUTABLE STATEMENT  DBVDER
  if (IVP  >  0) call DUIVP(X,Y(IVP+1),YP(IVP+1))
  NOFST = IVP
  NA = 1
  DO 10 K=1,NFC
     call DFMAT(X,Y(NA),YP(NA))
     NOFST = NOFST - NCOMP
     NA = NA + NCOMP
   10 CONTINUE
!
  if (INHOMO  /=  1) RETURN
  call DFMAT(X,Y(NA),YP(NA))
!
  if (IGOFX  ==  0) RETURN
  if (X  /=  XSAV) THEN
     if (IVP  ==  0) call DGVEC(X,G)
     if (IVP  >  0) call DUVEC(X,Y(IVP+1),G)
     XSAV = X
  end if
!
!     If the user has chosen not to normalize the particular
!     solution, then C is defined in DBVPOR to be 1.0
!
!     The following loop is just
!     call DAXPY (NCOMP, 1.0D0/C, G, 1, YP(NA), 1)
!
  DO 20 J=1,NCOMP
     L = NA + J - 1
     YP(L) = YP(L) + G(J)/C
   20 CONTINUE
  return
end
