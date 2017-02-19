subroutine EXBVP (Y, NROWY, XPTS, A, NROWA, ALPHA, B, NROWB, BETA, &
     IFLAG, WORK, IWORK)
!
!! EXBVP is subsidiary to BVSUP.
!
!***LIBRARY   SLATEC
!***TYPE      SINGLE PRECISION (EXBVP-S, DEXBVP-D)
!***AUTHOR  Watts, H. A., (SNLA)
!***DESCRIPTION
!
!  This subroutine is used to execute the basic technique for solving
!  the two-point boundary value problem
!
!***SEE ALSO  BVSUP
!***ROUTINES CALLED  BVPOR, XERMSG
!***COMMON BLOCKS    ML15TO, ML17BW, ML18JR, ML5MCO, ML8SZ
!***REVISION HISTORY  (YYMMDD)
!   750601  DATE WRITTEN
!   890921  Realigned order of variables in certain COMMON blocks.
!           (WRB)
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900328  Added TYPE section.  (WRB)
!   900510  Convert XERRWV calls to XERMSG calls.  (RWC)
!   910722  Updated AUTHOR section.  (ALS)
!***END PROLOGUE  EXBVP
!
  DIMENSION Y(NROWY,*),A(NROWA,*),ALPHA(*),B(NROWB,*),BETA(*), &
           WORK(*),IWORK(*),XPTS(*)
  CHARACTER*8 XERN1, XERN2
!
!     ****************************************************************
!
  COMMON /ML8SZ/ C,XSAV,IGOFX,INHOMO,IVP,NCOMP,NFC
  COMMON /ML18JR/ AE,RE,TOL,NXPTS,NIC,NOPG,MXNON,NDISK,NTAPE,NEQ, &
                  INDPVT,INTEG,NPS,NTP,NEQIVP,NUMORT,NFCC, &
                  ICOCO
  COMMON /ML15TO/ PX,PWCND,TND,X,XBEG,XEND,XOT,XOP,INFO(15),ISTKOP, &
                  KNSWOT,KOP,LOTJP,MNSWOT,NSWOT
  COMMON /ML17BW/ KKKZPW,NEEDW,NEEDIW,K1,K2,K3,K4,K5,K6,K7,K8,K9, &
                  K10,K11,L1,L2,KKKINT,LLLINT
!
  COMMON /ML5MCO/ URO,SRU,EPS,SQOVFL,TWOU,FOURU,LPAR
!
!***FIRST EXECUTABLE STATEMENT  EXBVP
  KOTC = 1
  IEXP = 0
  if (IWORK(7)  ==  -1) IEXP = IWORK(8)
!
!     COMPUTE ORTHONORMALIZATION TOLERANCES.
!
   10 TOL = 10.0**((-LPAR-IEXP)*2)
!
  IWORK(8) = IEXP
  MXNON = IWORK(2)
!
! **********************************************************************
! **********************************************************************
!
  call BVPOR(Y,NROWY,NCOMP,XPTS,NXPTS,A,NROWA,ALPHA,NIC,B, &
             NROWB,BETA,NFC,IFLAG,WORK(1),MXNON,WORK(K1),NTP, &
             IWORK(18),WORK(K2),IWORK(16),WORK(K3),WORK(K4), &
             WORK(K5),WORK(K6),WORK(K7),WORK(K8),WORK(K9), &
             WORK(K10),IWORK(L1),NFCC)
!
! **********************************************************************
! **********************************************************************
!     if MGSBV RETURNS WITH MESSAGE OF DEPENDENT VECTORS, WE REDUCE
!     ORTHONORMALIZATION TOLERANCE AND TRY AGAIN. THIS IS DONE
!     A MAXIMUM OF 2 TIMES.
!
  if (IFLAG  /=  30) go to 20
  if (KOTC  ==  3  .OR.  NOPG  ==  1) go to 30
  KOTC = KOTC + 1
  IEXP = IEXP - 2
  go to 10
!
! **********************************************************************
!     if BVPOR RETURNS MESSAGE THAT THE MAXIMUM NUMBER OF
!     ORTHONORMALIZATIONS HAS BEEN ATTAINED AND WE CANNOT CONTINUE, THEN
!     WE ESTIMATE THE NEW STORAGE REQUIREMENTS IN ORDER TO SOLVE PROBLEM
!
   20 if (IFLAG  /=  13) go to 30
  XL = ABS(XEND-XBEG)
  ZQUIT = ABS(X-XBEG)
  INC = 1.5 * XL/ZQUIT * (MXNON+1)
  if (NDISK  /=  1) THEN
     NSAFW = INC*KKKZPW + NEEDW
     NSAFIW = INC*NFCC + NEEDIW
  ELSE
     NSAFW = NEEDW + INC
     NSAFIW = NEEDIW
  end if
!
  WRITE (XERN1, '(I8)') NSAFW
  WRITE (XERN2, '(I8)') NSAFIW
  call XERMSG ('SLATEC', 'EXBVP', &
     'IN BVSUP, PREDICTED STORAGE ALLOCATION FOR WORK ARRAY IS ' // &
     XERN1 // ', PREDICTED STORAGE ALLOCATION FOR IWORK ARRAY IS ' &
     // XERN2, 1, 0)
!
   30 IWORK(1) = MXNON
  return
end
