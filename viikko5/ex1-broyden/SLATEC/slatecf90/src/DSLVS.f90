subroutine DSLVS (WM, IWM, X, TEM)
!
!! DSLVS is subsidiary to DDEBDF.
!
!***LIBRARY   SLATEC
!***TYPE      DOUBLE PRECISION (SLVS-S, DSLVS-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!   DSLVS solves the linear system in the iteration scheme for the
!   integrator package DDEBDF.
!
!***SEE ALSO  DDEBDF
!***ROUTINES CALLED  DGBSL, DGESL
!***COMMON BLOCKS    DDEBD1
!***REVISION HISTORY  (YYMMDD)
!   820301  DATE WRITTEN
!   890531  Changed all specific intrinsics to generic.  (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   910722  Updated AUTHOR section.  (ALS)
!   920422  Changed DIMENSION statement.  (WRB)
!***END PROLOGUE  DSLVS
!
  INTEGER I, IER, IOWND, IOWNS, IWM, JSTART, KFLAG, L, MAXORD, &
        MEBAND, METH, MITER, ML, MU, N, NFE, NJE, NQ, NQU, NST
  DOUBLE PRECISION DI, EL0, H, HL0, HMIN, HMXI, HU, PHL0, &
        R, ROWND, ROWNS, TEM, TN, UROUND, WM, X
  DIMENSION WM(*), IWM(*), X(*), TEM(*)
  COMMON /DDEBD1/ ROWND,ROWNS(210),EL0,H,HMIN,HMXI,HU,TN,UROUND, &
                  IOWND(14),IOWNS(6),IER,JSTART,KFLAG,L,METH,MITER, &
                  MAXORD,N,NQ,NST,NFE,NJE,NQU
!     ------------------------------------------------------------------
!      THIS ROUTINE MANAGES THE SOLUTION OF THE LINEAR SYSTEM ARISING
!      FROM A CHORD ITERATION.  IT IS CALLED BY DSTOD  if MITER  /=  0.
!      if MITER IS 1 OR 2, IT CALLS DGESL TO ACCOMPLISH THIS.
!      if MITER = 3 IT UPDATES THE COEFFICIENT H*EL0 IN THE DIAGONAL
!      MATRIX, AND THEN COMPUTES THE SOLUTION.
!      if MITER IS 4 OR 5, IT CALLS DGBSL.
!      COMMUNICATION WITH DSLVS USES THE FOLLOWING VARIABLES..
!      WM  = DOUBLE PRECISION WORK SPACE CONTAINING THE INVERSE DIAGONAL
!      MATRIX if MITER
!            IS 3 AND THE LU DECOMPOSITION OF THE MATRIX OTHERWISE.
!            STORAGE OF MATRIX ELEMENTS STARTS AT WM(3).
!            WM ALSO CONTAINS THE FOLLOWING MATRIX-RELATED DATA..
!            WM(1) = SQRT(UROUND) (NOT USED HERE),
!            WM(2) = HL0, THE PREVIOUS VALUE OF H*EL0, USED if MITER =
!            3.
!      IWM = INTEGER WORK SPACE CONTAINING PIVOT INFORMATION, STARTING
!            AT IWM(21), if MITER IS 1, 2, 4, OR 5.  IWM ALSO CONTAINS
!            THE BAND PARAMETERS ML = IWM(1) AND MU = IWM(2) if MITER IS
!            4 OR 5.
!      X   = THE RIGHT-HAND SIDE VECTOR ON INPUT, AND THE SOLUTION
!            VECTOR ON OUTPUT, OF LENGTH N.
!      TEM = VECTOR OF WORK SPACE OF LENGTH N, NOT USED IN THIS VERSION.
!      IER = OUTPUT FLAG (IN COMMON).  IER = 0 if NO TROUBLE OCCURRED.
!            IER = -1 if A SINGULAR MATRIX AROSE WITH MITER = 3.
!      THIS ROUTINE ALSO USES THE COMMON VARIABLES EL0, H, MITER, AND N.
!-----------------------------------------------------------------------
!     BEGIN BLOCK PERMITTING ...EXITS TO 80
!        BEGIN BLOCK PERMITTING ...EXITS TO 60
!***FIRST EXECUTABLE STATEMENT  DSLVS
        IER = 0
        go to (10,10,20,70,70), MITER
   10       CONTINUE
        call DGESL(WM(3),N,N,IWM(21),X,0)
!     ......EXIT
        go to 80
!
   20       CONTINUE
        PHL0 = WM(2)
        HL0 = H*EL0
        WM(2) = HL0
        if (HL0  ==  PHL0) go to 40
           R = HL0/PHL0
           DO 30 I = 1, N
              DI = 1.0D0 - R*(1.0D0 - 1.0D0/WM(I+2))
!        .........EXIT
              if (ABS(DI)  ==  0.0D0) go to 60
              WM(I+2) = 1.0D0/DI
   30          CONTINUE
   40       CONTINUE
        DO 50 I = 1, N
           X(I) = WM(I+2)*X(I)
   50       CONTINUE
!     ......EXIT
        go to 80
   60    CONTINUE
     IER = -1
!     ...EXIT
     go to 80
!
   70    CONTINUE
     ML = IWM(1)
     MU = IWM(2)
     MEBAND = 2*ML + MU + 1
     call DGBSL(WM(3),MEBAND,N,ML,MU,IWM(21),X,0)
   80 CONTINUE
  return
end
