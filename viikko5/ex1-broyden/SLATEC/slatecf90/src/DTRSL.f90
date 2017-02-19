subroutine DTRSL (T, LDT, N, B, JOB, INFO)
!
!! DTRSL solves a triangular system of linear equations.
!
!***PURPOSE  Solve a system of the form  T*X=B or TRANS(T)*X=B, where
!            T is a triangular matrix.
!***LIBRARY   SLATEC (LINPACK)
!***CATEGORY  D2A3
!***TYPE      DOUBLE PRECISION (STRSL-S, DTRSL-D, CTRSL-C)
!***KEYWORDS  LINEAR ALGEBRA, LINPACK, TRIANGULAR LINEAR SYSTEM,
!             TRIANGULAR MATRIX
!***AUTHOR  Stewart, G. W., (U. of Maryland)
!***DESCRIPTION
!
!     DTRSL solves systems of the form
!
!                   T * X = B
!     or
!                   TRANS(T) * X = B
!
!     where T is a triangular matrix of order N.  Here TRANS(T)
!     denotes the transpose of the matrix T.
!
!     On Entry
!
!         T         DOUBLE PRECISION(LDT,N)
!                   T contains the matrix of the system.  The zero
!                   elements of the matrix are not referenced, and
!                   the corresponding elements of the array can be
!                   used to store other information.
!
!         LDT       INTEGER
!                   LDT is the leading dimension of the array T.
!
!         N         INTEGER
!                   N is the order of the system.
!
!         B         DOUBLE PRECISION(N).
!                   B contains the right hand side of the system.
!
!         JOB       INTEGER
!                   JOB specifies what kind of system is to be solved.
!                   If JOB is
!
!                        00   solve T*X=B, T lower triangular,
!                        01   solve T*X=B, T upper triangular,
!                        10   solve TRANS(T)*X=B, T lower triangular,
!                        11   solve TRANS(T)*X=B, T upper triangular.
!
!     On Return
!
!         B         B contains the solution, if INFO  ==  0.
!                   Otherwise B is unaltered.
!
!         INFO      INTEGER
!                   INFO contains zero if the system is nonsingular.
!                   Otherwise INFO contains the index of
!                   the first zero diagonal element of T.
!
!***REFERENCES  J. J. Dongarra, J. R. Bunch, C. B. Moler, and G. W.
!                 Stewart, LINPACK Users' Guide, SIAM, 1979.
!***ROUTINES CALLED  DAXPY, DDOT
!***REVISION HISTORY  (YYMMDD)
!   780814  DATE WRITTEN
!   890831  Modified array declarations.  (WRB)
!   890831  REVISION DATE from Version 3.2
!   891214  Prologue converted to Version 4.0 format.  (BAB)
!   900326  Removed duplicate information from DESCRIPTION section.
!           (WRB)
!   920501  Reformatted the REFERENCES section.  (WRB)
!***END PROLOGUE  DTRSL
  INTEGER LDT,N,JOB,INFO
  DOUBLE PRECISION T(LDT,*),B(*)
!
!
  DOUBLE PRECISION DDOT,TEMP
  INTEGER CASE,J,JJ
!***FIRST EXECUTABLE STATEMENT  DTRSL
!
!        CHECK FOR ZERO DIAGONAL ELEMENTS.
!
     DO 10 INFO = 1, N
        if (T(INFO,INFO)  ==  0.0D0) go to 150
   10    CONTINUE
     INFO = 0
!
!        DETERMINE THE TASK AND go to IT.
!
     CASE = 1
     if (MOD(JOB,10)  /=  0) CASE = 2
     if (MOD(JOB,100)/10  /=  0) CASE = CASE + 2
     go to (20,50,80,110), CASE
!
!        SOLVE T*X=B FOR T LOWER TRIANGULAR
!
   20    CONTINUE
        B(1) = B(1)/T(1,1)
        if (N  <  2) go to 40
        DO 30 J = 2, N
           TEMP = -B(J-1)
           call DAXPY(N-J+1,TEMP,T(J,J-1),1,B(J),1)
           B(J) = B(J)/T(J,J)
   30       CONTINUE
   40       CONTINUE
     go to 140
!
!        SOLVE T*X=B FOR T UPPER TRIANGULAR.
!
   50    CONTINUE
        B(N) = B(N)/T(N,N)
        if (N  <  2) go to 70
        DO 60 JJ = 2, N
           J = N - JJ + 1
           TEMP = -B(J+1)
           call DAXPY(J,TEMP,T(1,J+1),1,B(1),1)
           B(J) = B(J)/T(J,J)
   60       CONTINUE
   70       CONTINUE
     go to 140
!
!        SOLVE TRANS(T)*X=B FOR T LOWER TRIANGULAR.
!
   80    CONTINUE
        B(N) = B(N)/T(N,N)
        if (N  <  2) go to 100
        DO 90 JJ = 2, N
           J = N - JJ + 1
           B(J) = B(J) - DDOT(JJ-1,T(J+1,J),1,B(J+1),1)
           B(J) = B(J)/T(J,J)
   90       CONTINUE
  100       CONTINUE
     go to 140
!
!        SOLVE TRANS(T)*X=B FOR T UPPER TRIANGULAR.
!
  110    CONTINUE
        B(1) = B(1)/T(1,1)
        if (N  <  2) go to 130
        DO 120 J = 2, N
           B(J) = B(J) - DDOT(J-1,T(1,J),1,B(1),1)
           B(J) = B(J)/T(J,J)
  120       CONTINUE
  130       CONTINUE
  140    CONTINUE
  150 CONTINUE
  return
end
