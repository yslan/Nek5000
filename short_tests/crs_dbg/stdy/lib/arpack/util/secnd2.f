      SUBROUTINE SECND2( T )
*
      REAL (8)   T
*
*  -- LAPACK auxiliary routine (preliminary version) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     July 26, 1991
*
*  Purpose
*  =======
*
*  SECOND returns the user time for a process in seconds.
*  This version gets the time from the system function ETIME.
*
*     .. Local Scalars ..
c     REAL(4)             T1
*     ..
*     .. Local Arrays ..
c     REAL(4)            TARRAY( 2 )
*     ..
*     .. External Functions ..
c     REAL(4)            ETIME
c     EXTERNAL           ETIME
*     ..
*     .. Executable Statements ..
*

c     T1 = ETIME( TARRAY )
c     T  = real(TARRAY( 1 ), 8)

      T = 0.0E0
      RETURN
*
*     End of SECND2
*
      END
