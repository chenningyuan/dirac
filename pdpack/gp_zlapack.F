      SUBROUTINE ZHEEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, RWORK,
     $                  INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
*  complex Hermitian matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,2*N-1).
*          For optimal efficiency, LWORK >= (NB+1)*N,
*          where NB is the blocksize for ZHETRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1, 3*N-2))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LOPT, LWKOPT, NB
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANHE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANHE
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTERF, XERBLA, ZHETRD, ZLASCL, ZSTEQR,
     $                   ZUNGTR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, 2*N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+1 )*N )
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHEEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 3
         IF( WANTZ )
     $      A( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = ZLANHE( 'M', UPLO, N, A, LDA, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL ZLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
*
*     Call ZHETRD to reduce Hermitian matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = 1
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL ZHETRD( UPLO, N, A, LDA, W, RWORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
      LOPT = N + WORK( INDWRK )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, first call
*     ZUNGTR to generate the unitary matrix, then call ZSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         CALL ZUNGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
         INDWRK = INDE + N
         CALL ZSTEQR( JOBZ, N, W, RWORK( INDE ), A, LDA,
     $                RWORK( INDWRK ), INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     Set WORK(1) to optimal complex workspace size.
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of ZHEEV
*
      END
      SUBROUTINE ZHETD2( UPLO, N, A, LDA, D, E, TAU, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAU( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHETD2 reduces a complex Hermitian matrix A to real symmetric
*  tridiagonal form T by a unitary similarity transformation:
*  Q' * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the unitary
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the unitary matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) COMPLEX*16 array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO, HALF
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   HALF = ( 0.5D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            UPPER
      INTEGER            I
      COMPLEX*16         ALPHA, TAUI
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZAXPY, ZHEMV, ZHER2, ZLARFG
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETD2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A
*
         A( N, N ) = DBLE( A( N, N ) )
         DO 10 I = N - 1, 1, -1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(1:i-1,i+1)
*
            ALPHA = A( I, I+1 )
            CALL ZLARFG( I, ALPHA, A( 1, I+1 ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(1:i,1:i)
*
               A( I, I+1 ) = ONE
*
*              Compute  x := tau * A * v  storing x in TAU(1:i)
*
               CALL ZHEMV( UPLO, I, TAUI, A, LDA, A( 1, I+1 ), 1, ZERO,
     $                     TAU, 1 )
*
*              Compute  w := x - 1/2 * tau * (x'*v) * v
*
               ALPHA = -HALF*TAUI*ZDOTC( I, TAU, 1, A( 1, I+1 ), 1 )
               CALL ZAXPY( I, ALPHA, A( 1, I+1 ), 1, TAU, 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL ZHER2( UPLO, I, -ONE, A( 1, I+1 ), 1, TAU, 1, A,
     $                     LDA )
*
            ELSE
               A( I, I ) = DBLE( A( I, I ) )
            END IF
            A( I, I+1 ) = E( I )
            D( I+1 ) = A( I+1, I+1 )
            TAU( I ) = TAUI
   10    CONTINUE
         D( 1 ) = A( 1, 1 )
      ELSE
*
*        Reduce the lower triangle of A
*
         A( 1, 1 ) = DBLE( A( 1, 1 ) )
         DO 20 I = 1, N - 1
*
*           Generate elementary reflector H(i) = I - tau * v * v'
*           to annihilate A(i+2:n,i)
*
            ALPHA = A( I+1, I )
            CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1, TAUI )
            E( I ) = ALPHA
*
            IF( TAUI.NE.ZERO ) THEN
*
*              Apply H(i) from both sides to A(i+1:n,i+1:n)
*
               A( I+1, I ) = ONE
*
*              Compute  x := tau * A * v  storing y in TAU(i:n-1)
*
               CALL ZHEMV( UPLO, N-I, TAUI, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, TAU( I ), 1 )
*
*              Compute  w := x - 1/2 * tau * (x'*v) * v
*
               ALPHA = -HALF*TAUI*ZDOTC( N-I, TAU( I ), 1, A( I+1, I ),
     $                 1 )
               CALL ZAXPY( N-I, ALPHA, A( I+1, I ), 1, TAU( I ), 1 )
*
*              Apply the transformation as a rank-2 update:
*                 A := A - v * w' - w * v'
*
               CALL ZHER2( UPLO, N-I, -ONE, A( I+1, I ), 1, TAU( I ), 1,
     $                     A( I+1, I+1 ), LDA )
*
            ELSE
               A( I+1, I+1 ) = DBLE( A( I+1, I+1 ) )
            END IF
            A( I+1, I ) = E( I )
            D( I ) = A( I, I )
            TAU( I ) = TAUI
   20    CONTINUE
         D( N ) = A( N, N )
      END IF
*
      RETURN
*
*     End of ZHETD2
*
      END
      SUBROUTINE ZHETRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHETRD reduces a complex Hermitian matrix A to real symmetric
*  tridiagonal form T by a unitary similarity transformation:
*  Q**H * A * Q = T.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit, if UPLO = 'U', the diagonal and first superdiagonal
*          of A are overwritten by the corresponding elements of the
*          tridiagonal matrix T, and the elements above the first
*          superdiagonal, with the array TAU, represent the unitary
*          matrix Q as a product of elementary reflectors; if UPLO
*          = 'L', the diagonal and first subdiagonal of A are over-
*          written by the corresponding elements of the tridiagonal
*          matrix T, and the elements below the first subdiagonal, with
*          the array TAU, represent the unitary matrix Q as a product
*          of elementary reflectors. See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*
*  TAU     (output) COMPLEX*16 array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= 1.
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n-1) . . . H(2) H(1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*  A(1:i-1,i+1), and tau in TAU(i).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(n-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*  and tau in TAU(i).
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  d   e   v2  v3  v4 )              (  d                  )
*    (      d   e   v3  v4 )              (  e   d              )
*    (          d   e   v4 )              (  v1  e   d          )
*    (              d   e  )              (  v1  v2  e   d      )
*    (                  d  )              (  v1  v2  v3  e   d  )
*
*  where d and e denote diagonal and off-diagonal elements of T, and vi
*  denotes an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D+0 )
      COMPLEX*16         CONE
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, IWS, J, KK, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZHER2K, ZHETD2, ZLATRD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.1 .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.
*
         NB = ILAENV( 1, 'ZHETRD', UPLO, N, -1, -1, -1 )
         LWKOPT = N*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHETRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NX = N
      IWS = 1
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
*
*        Determine when to cross over from blocked to unblocked code
*        (last block is always handled by unblocked code).
*
         NX = MAX( NB, ILAENV( 3, 'ZHETRD', UPLO, N, -1, -1, -1 ) )
         IF( NX.LT.N ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  determine the
*              minimum value of NB, and reduce NB or force use of
*              unblocked code by setting NX = N.
*
               NB = MAX( LWORK / LDWORK, 1 )
               NBMIN = ILAENV( 2, 'ZHETRD', UPLO, N, -1, -1, -1 )
               IF( NB.LT.NBMIN )
     $            NX = N
            END IF
         ELSE
            NX = N
         END IF
      ELSE
         NB = 1
      END IF
*
      IF( UPPER ) THEN
*
*        Reduce the upper triangle of A.
*        Columns 1:kk are handled by the unblocked method.
*
         KK = N - ( ( N-NX+NB-1 ) / NB )*NB
         DO 20 I = N - NB + 1, KK + 1, -NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL ZLATRD( UPLO, I+NB-1, NB, A, LDA, E, TAU, WORK,
     $                   LDWORK )
*
*           Update the unreduced submatrix A(1:i-1,1:i-1), using an
*           update of the form:  A := A - V*W' - W*V'
*
            CALL ZHER2K( UPLO, 'No transpose', I-1, NB, -CONE,
     $                   A( 1, I ), LDA, WORK, LDWORK, ONE, A, LDA )
*
*           Copy superdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 10 J = I, I + NB - 1
               A( J-1, J ) = E( J-1 )
               D( J ) = A( J, J )
   10       CONTINUE
   20    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL ZHETD2( UPLO, KK, A, LDA, D, E, TAU, IINFO )
      ELSE
*
*        Reduce the lower triangle of A
*
         DO 40 I = 1, N - NX, NB
*
*           Reduce columns i:i+nb-1 to tridiagonal form and form the
*           matrix W which is needed to update the unreduced part of
*           the matrix
*
            CALL ZLATRD( UPLO, N-I+1, NB, A( I, I ), LDA, E( I ),
     $                   TAU( I ), WORK, LDWORK )
*
*           Update the unreduced submatrix A(i+nb:n,i+nb:n), using
*           an update of the form:  A := A - V*W' - W*V'
*
            CALL ZHER2K( UPLO, 'No transpose', N-I-NB+1, NB, -CONE,
     $                   A( I+NB, I ), LDA, WORK( NB+1 ), LDWORK, ONE,
     $                   A( I+NB, I+NB ), LDA )
*
*           Copy subdiagonal elements back into A, and diagonal
*           elements into D
*
            DO 30 J = I, I + NB - 1
               A( J+1, J ) = E( J )
               D( J ) = A( J, J )
   30       CONTINUE
   40    CONTINUE
*
*        Use unblocked code to reduce the last or only block
*
         CALL ZHETD2( UPLO, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $                TAU( I ), IINFO )
      END IF
*
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of ZHETRD
*
      END
      SUBROUTINE ZLACGV( N, X, INCX )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLACGV conjugates a complex vector of length N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vector X.  N >= 0.
*
*  X       (input/output) COMPLEX*16 array, dimension
*                         (1+(N-1)*abs(INCX))
*          On entry, the vector of length N to be conjugated.
*          On exit, X is overwritten with conjg(X).
*
*  INCX    (input) INTEGER
*          The spacing between successive elements of X.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IOFF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. Executable Statements ..
*
      IF( INCX.EQ.1 ) THEN
         DO 10 I = 1, N
            X( I ) = DCONJG( X( I ) )
   10    CONTINUE
      ELSE
         IOFF = 1
         IF( INCX.LT.0 )
     $      IOFF = 1 - ( N-1 )*INCX
         DO 20 I = 1, N
            X( IOFF ) = DCONJG( X( IOFF ) )
            IOFF = IOFF + INCX
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZLACGV
*
      END
      DOUBLE COMPLEX   FUNCTION ZLADIV( X, Y )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      COMPLEX*16         X, Y
*     ..
*
*  Purpose
*  =======
*
*  ZLADIV := X / Y, where X and Y are complex.  The computation of X / Y
*  will not overflow on an intermediary step unless the results
*  overflows.
*
*  Arguments
*  =========
*
*  X       (input) COMPLEX*16
*  Y       (input) COMPLEX*16
*          The complex scalars X and Y.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   ZI, ZR
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DIMAG
*     ..
*     .. Executable Statements ..
*
      CALL DLADIV( DBLE( X ), DIMAG( X ), DBLE( Y ), DIMAG( Y ), ZR,
     $             ZI )
      ZLADIV = DCMPLX( ZR, ZI )
*
      RETURN
*
*     End of ZLADIV
*
      END
      DOUBLE PRECISION FUNCTION ZLANHE( NORM, UPLO, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLANHE  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  complex hermitian matrix A.
*
*  Description
*  ===========
*
*  ZLANHE returns the value
*
*     ZLANHE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in ZLANHE as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          hermitian matrix A is to be referenced.
*          = 'U':  Upper triangular part of A is referenced
*          = 'L':  Lower triangular part of A is referenced
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, ZLANHE is
*          set to zero.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The hermitian matrix A.  If UPLO = 'U', the leading n by n
*          upper triangular part of A contains the upper triangular part
*          of the matrix A, and the strictly lower triangular part of A
*          is not referenced.  If UPLO = 'L', the leading n by n lower
*          triangular part of A contains the lower triangular part of
*          the matrix A, and the strictly upper triangular part of A is
*          not referenced. Note that the imaginary parts of the diagonal
*          elements need not be set and are assumed to be zero.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
*          WORK is not referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = 1, J - 1
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10          CONTINUE
               VALUE = MAX( VALUE, ABS( DBLE( A( J, J ) ) ) )
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               VALUE = MAX( VALUE, ABS( DBLE( A( J, J ) ) ) )
               DO 30 I = J + 1, N
                  VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        Find normI(A) ( = norm1(A), since A is hermitian).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               DO 50 I = 1, J - 1
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( DBLE( A( J, J ) ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( DBLE( A( J, J ) ) )
               DO 90 I = J + 1, N
                  ABSA = ABS( A( I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 110 J = 2, N
               CALL ZLASSQ( J-1, A( 1, J ), 1, SCALE, SUM )
  110       CONTINUE
         ELSE
            DO 120 J = 1, N - 1
               CALL ZLASSQ( N-J, A( J+1, J ), 1, SCALE, SUM )
  120       CONTINUE
         END IF
         SUM = 2*SUM
         DO 130 I = 1, N
            IF( DBLE( A( I, I ) ).NE.ZERO ) THEN
               ABSA = ABS( DBLE( A( I, I ) ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      ZLANHE = VALUE
      RETURN
*
*     End of ZLANHE
*
      END
      SUBROUTINE ZLARF( SIDE, M, N, V, INCV, TAU, C, LDC, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE
      INTEGER            INCV, LDC, M, N
      COMPLEX*16         TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARF applies a complex elementary reflector H to a complex M-by-N
*  matrix C, from either the left or the right. H is represented in the
*  form
*
*        H = I - tau * v * v'
*
*  where tau is a complex scalar and v is a complex vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  To apply H' (the conjugate transpose of H), supply conjg(tau) instead
*  tau.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) COMPLEX*16 array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) COMPLEX*16
*          The value tau in the representation of H.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX*16 array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZGERC
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  H * C
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C' * v
*
            CALL ZGEMV( 'Conjugate transpose', M, N, ONE, C, LDC, V,
     $                  INCV, ZERO, WORK, 1 )
*
*           C := C - v * w'
*
            CALL ZGERC( M, N, -TAU, V, INCV, WORK, 1, C, LDC )
         END IF
      ELSE
*
*        Form  C * H
*
         IF( TAU.NE.ZERO ) THEN
*
*           w := C * v
*
            CALL ZGEMV( 'No transpose', M, N, ONE, C, LDC, V, INCV,
     $                  ZERO, WORK, 1 )
*
*           C := C - w * v'
*
            CALL ZGERC( M, N, -TAU, WORK, 1, V, INCV, C, LDC )
         END IF
      END IF
      RETURN
*
*     End of ZLARF
*
      END
      SUBROUTINE ZLARFB( SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV,
     $                   T, LDT, C, LDC, WORK, LDWORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, SIDE, STOREV, TRANS
      INTEGER            K, LDC, LDT, LDV, LDWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), T( LDT, * ), V( LDV, * ),
     $                   WORK( LDWORK, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARFB applies a complex block reflector H or its transpose H' to a
*  complex M-by-N matrix C, from either the left or the right.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H' from the Left
*          = 'R': apply H or H' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'C': apply H' (Conjugate transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  V       (input) COMPLEX*16 array, dimension
*                                (LDV,K) if STOREV = 'C'
*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*          if STOREV = 'R', LDV >= K.
*
*  T       (input) COMPLEX*16 array, dimension (LDT,K)
*          The triangular K-by-K matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      CHARACTER          TRANST
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZCOPY, ZGEMM, ZLACGV, ZTRMM
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( LSAME( TRANS, 'N' ) ) THEN
         TRANST = 'C'
      ELSE
         TRANST = 'N'
      END IF
*
      IF( LSAME( STOREV, 'C' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1 )    (first K rows)
*                     ( V2 )
*           where  V1  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C1'
*
               DO 10 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   10          CONTINUE
*
*              W := W * V1
*
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2
*
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N,
     $                        K, M-K, ONE, C( K+1, 1 ), LDC,
     $                        V( K+1, 1 ), LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2 * W'
*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose',
     $                        M-K, N, K, -ONE, V( K+1, 1 ), LDV, WORK,
     $                        LDWORK, ONE, C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 30 J = 1, K
                  DO 20 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
   20             CONTINUE
   30          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C1
*
               DO 40 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
   40          CONTINUE
*
*              W := W * V1
*
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2
*
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV,
     $                        ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2'
*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        N-K, K, -ONE, WORK, LDWORK, V( K+1, 1 ),
     $                        LDV, ONE, C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1'
*
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 60 J = 1, K
                  DO 50 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
   50             CONTINUE
   60          CONTINUE
            END IF
*
         ELSE
*
*           Let  V =  ( V1 )
*                     ( V2 )    (last K rows)
*           where  V2  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V  =  (C1'*V1 + C2'*V2)  (stored in WORK)
*
*              W := C2'
*
               DO 70 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
   70          CONTINUE
*
*              W := W * V2
*
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1
*
                  CALL ZGEMM( 'Conjugate transpose', 'No transpose', N,
     $                        K, M-K, ONE, C, LDC, V, LDV, ONE, WORK,
     $                        LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1 * W'
*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose',
     $                        M-K, N, K, -ONE, V, LDV, WORK, LDWORK,
     $                        ONE, C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK,
     $                     LDWORK )
*
*              C2 := C2 - W'
*
               DO 90 J = 1, K
                  DO 80 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -
     $                               DCONJG( WORK( I, J ) )
   80             CONTINUE
   90          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)
*
*              W := C2
*
               DO 100 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  100          CONTINUE
*
*              W := W * V2
*
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1
*
                  CALL ZGEMM( 'No transpose', 'No transpose', M, K, N-K,
     $                        ONE, C, LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V'
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1'
*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        N-K, K, -ONE, WORK, LDWORK, V, LDV, ONE,
     $                        C, LDC )
               END IF
*
*              W := W * V2'
*
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK,
     $                     LDWORK )
*
*              C2 := C2 - W
*
               DO 120 J = 1, K
                  DO 110 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  110             CONTINUE
  120          CONTINUE
            END IF
         END IF
*
      ELSE IF( LSAME( STOREV, 'R' ) ) THEN
*
         IF( LSAME( DIRECT, 'F' ) ) THEN
*
*           Let  V =  ( V1  V2 )    (V1: first K columns)
*           where  V1  is unit upper triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C1'
*
               DO 130 J = 1, K
                  CALL ZCOPY( N, C( J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  130          CONTINUE
*
*              W := W * V1'
*
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V, LDV, WORK, LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C2'*V2'
*
                  CALL ZGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', N, K, M-K, ONE,
     $                        C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE,
     $                        WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL ZTRMM( 'Right', 'Upper', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C2 := C2 - V2' * W'
*
                  CALL ZGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', M-K, N, K, -ONE,
     $                        V( 1, K+1 ), LDV, WORK, LDWORK, ONE,
     $                        C( K+1, 1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', N,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W'
*
               DO 150 J = 1, K
                  DO 140 I = 1, N
                     C( J, I ) = C( J, I ) - DCONJG( WORK( I, J ) )
  140             CONTINUE
  150          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C1
*
               DO 160 J = 1, K
                  CALL ZCOPY( M, C( 1, J ), 1, WORK( 1, J ), 1 )
  160          CONTINUE
*
*              W := W * V1'
*
               CALL ZTRMM( 'Right', 'Upper', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V, LDV, WORK, LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C2 * V2'
*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        K, N-K, ONE, C( 1, K+1 ), LDC,
     $                        V( 1, K+1 ), LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL ZTRMM( 'Right', 'Upper', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C2 := C2 - W * V2
*
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE,
     $                        C( 1, K+1 ), LDC )
               END IF
*
*              W := W * V1
*
               CALL ZTRMM( 'Right', 'Upper', 'No transpose', 'Unit', M,
     $                     K, ONE, V, LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 180 J = 1, K
                  DO 170 I = 1, M
                     C( I, J ) = C( I, J ) - WORK( I, J )
  170             CONTINUE
  180          CONTINUE
*
            END IF
*
         ELSE
*
*           Let  V =  ( V1  V2 )    (V2: last K columns)
*           where  V2  is unit lower triangular.
*
            IF( LSAME( SIDE, 'L' ) ) THEN
*
*              Form  H * C  or  H' * C  where  C = ( C1 )
*                                                  ( C2 )
*
*              W := C' * V'  =  (C1'*V1' + C2'*V2') (stored in WORK)
*
*              W := C2'
*
               DO 190 J = 1, K
                  CALL ZCOPY( N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 )
                  CALL ZLACGV( N, WORK( 1, J ), 1 )
  190          CONTINUE
*
*              W := W * V2'
*
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', N, K, ONE, V( 1, M-K+1 ), LDV, WORK,
     $                     LDWORK )
               IF( M.GT.K ) THEN
*
*                 W := W + C1'*V1'
*
                  CALL ZGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', N, K, M-K, ONE, C,
     $                        LDC, V, LDV, ONE, WORK, LDWORK )
               END IF
*
*              W := W * T'  or  W * T
*
               CALL ZTRMM( 'Right', 'Lower', TRANST, 'Non-unit', N, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - V' * W'
*
               IF( M.GT.K ) THEN
*
*                 C1 := C1 - V1' * W'
*
                  CALL ZGEMM( 'Conjugate transpose',
     $                        'Conjugate transpose', M-K, N, K, -ONE, V,
     $                        LDV, WORK, LDWORK, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', N,
     $                     K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK )
*
*              C2 := C2 - W'
*
               DO 210 J = 1, K
                  DO 200 I = 1, N
                     C( M-K+J, I ) = C( M-K+J, I ) -
     $                               DCONJG( WORK( I, J ) )
  200             CONTINUE
  210          CONTINUE
*
            ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*              Form  C * H  or  C * H'  where  C = ( C1  C2 )
*
*              W := C * V'  =  (C1*V1' + C2*V2')  (stored in WORK)
*
*              W := C2
*
               DO 220 J = 1, K
                  CALL ZCOPY( M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 )
  220          CONTINUE
*
*              W := W * V2'
*
               CALL ZTRMM( 'Right', 'Lower', 'Conjugate transpose',
     $                     'Unit', M, K, ONE, V( 1, N-K+1 ), LDV, WORK,
     $                     LDWORK )
               IF( N.GT.K ) THEN
*
*                 W := W + C1 * V1'
*
                  CALL ZGEMM( 'No transpose', 'Conjugate transpose', M,
     $                        K, N-K, ONE, C, LDC, V, LDV, ONE, WORK,
     $                        LDWORK )
               END IF
*
*              W := W * T  or  W * T'
*
               CALL ZTRMM( 'Right', 'Lower', TRANS, 'Non-unit', M, K,
     $                     ONE, T, LDT, WORK, LDWORK )
*
*              C := C - W * V
*
               IF( N.GT.K ) THEN
*
*                 C1 := C1 - W * V1
*
                  CALL ZGEMM( 'No transpose', 'No transpose', M, N-K, K,
     $                        -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC )
               END IF
*
*              W := W * V2
*
               CALL ZTRMM( 'Right', 'Lower', 'No transpose', 'Unit', M,
     $                     K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK )
*
*              C1 := C1 - W
*
               DO 240 J = 1, K
                  DO 230 I = 1, M
                     C( I, N-K+J ) = C( I, N-K+J ) - WORK( I, J )
  230             CONTINUE
  240          CONTINUE
*
            END IF
*
         END IF
      END IF
*
      RETURN
*
*     End of ZLARFB
*
      END
      SUBROUTINE ZLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX*16         ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARFG generates a complex elementary reflector H of order n, such
*  that
*
*        H' * ( alpha ) = ( beta ),   H' * H = I.
*             (   x   )   (   0  )
*
*  where alpha and beta are scalars, with beta real, and x is an
*  (n-1)-element complex vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a complex scalar and v is a complex (n-1)-element
*  vector. Note that H is not hermitian.
*
*  If the elements of x are all zero and alpha is real, then tau = 0
*  and H is taken to be the unit matrix.
*
*  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) COMPLEX*16
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) COMPLEX*16 array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) COMPLEX*16
*          The value tau.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            J, KNT
      DOUBLE PRECISION   ALPHI, ALPHR, BETA, RSAFMN, SAFMIN, XNORM
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY3, DZNRM2
      COMPLEX*16         ZLADIV
      EXTERNAL           DLAMCH, DLAPY3, DZNRM2, ZLADIV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, SIGN
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZDSCAL, ZSCAL
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 ) THEN
         TAU = ZERO
         RETURN
      END IF
*
      XNORM = DZNRM2( N-1, X, INCX )
      ALPHR = DBLE( ALPHA )
      ALPHI = DIMAG( ALPHA )
*
      IF( XNORM.EQ.ZERO .AND. ALPHI.EQ.ZERO ) THEN
*
*        H  =  I
*
         TAU = ZERO
      ELSE
*
*        general case
*
         BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
         SAFMIN = DLAMCH( 'S' ) / DLAMCH( 'E' )
         RSAFMN = ONE / SAFMIN
*
         IF( ABS( BETA ).LT.SAFMIN ) THEN
*
*           XNORM, BETA may be inaccurate; scale X and recompute them
*
            KNT = 0
   10       CONTINUE
            KNT = KNT + 1
            CALL ZDSCAL( N-1, RSAFMN, X, INCX )
            BETA = BETA*RSAFMN
            ALPHI = ALPHI*RSAFMN
            ALPHR = ALPHR*RSAFMN
            IF( ABS( BETA ).LT.SAFMIN )
     $         GO TO 10
*
*           New BETA is at most 1, at least SAFMIN
*
            XNORM = DZNRM2( N-1, X, INCX )
            ALPHA = DCMPLX( ALPHR, ALPHI )
            BETA = -SIGN( DLAPY3( ALPHR, ALPHI, XNORM ), ALPHR )
            TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
            CALL ZSCAL( N-1, ALPHA, X, INCX )
*
*           If ALPHA is subnormal, it may lose relative accuracy
*
            ALPHA = BETA
            DO 20 J = 1, KNT
               ALPHA = ALPHA*SAFMIN
   20       CONTINUE
         ELSE
            TAU = DCMPLX( ( BETA-ALPHR ) / BETA, -ALPHI / BETA )
            ALPHA = ZLADIV( DCMPLX( ONE ), ALPHA-BETA )
            CALL ZSCAL( N-1, ALPHA, X, INCX )
            ALPHA = BETA
         END IF
      END IF
*
      RETURN
*
*     End of ZLARFG
*
      END
      SUBROUTINE ZLARFT( DIRECT, STOREV, N, K, V, LDV, TAU, T, LDT )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, STOREV
      INTEGER            K, LDT, LDV, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         T( LDT, * ), TAU( * ), V( LDV, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARFT forms the triangular factor T of a complex block reflector H
*  of order n, which is defined as a product of k elementary reflectors.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the array V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the array V, and
*
*     H  =  I - V' * T * V
*
*  Arguments
*  =========
*
*  DIRECT  (input) CHARACTER*1
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise
*          = 'R': rowwise
*
*  N       (input) INTEGER
*          The order of the block reflector H. N >= 0.
*
*  K       (input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). K >= 1.
*
*  V       (input/output) COMPLEX*16 array, dimension
*                               (LDV,K) if STOREV = 'C'
*                               (LDV,N) if STOREV = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i).
*
*  T       (output) COMPLEX*16 array, dimension (LDT,K)
*          The k by k triangular factor T of the block reflector.
*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*          lower triangular. The rest of the array is not used.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  Further Details
*  ===============
*
*  The shape of the matrix V and the storage of the vectors which define
*  the H(i) is best illustrated by the following example with n = 5 and
*  k = 3. The elements equal to 1 are not stored; the corresponding
*  array elements are modified but restored on exit. The rest of the
*  array is not used.
*
*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*
*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*                   ( v1  1    )                     (     1 v2 v2 v2 )
*                   ( v1 v2  1 )                     (        1 v3 v3 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*
*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*                   (     1 v3 )
*                   (        1 )
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      COMPLEX*16         VII
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZLACGV, ZTRMV
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( DIRECT, 'F' ) ) THEN
         DO 20 I = 1, K
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 10 J = 1, I
                  T( J, I ) = ZERO
   10          CONTINUE
            ELSE
*
*              general case
*
               VII = V( I, I )
               V( I, I ) = ONE
               IF( LSAME( STOREV, 'C' ) ) THEN
*
*                 T(1:i-1,i) := - tau(i) * V(i:n,1:i-1)' * V(i:n,i)
*
                  CALL ZGEMV( 'Conjugate transpose', N-I+1, I-1,
     $                        -TAU( I ), V( I, 1 ), LDV, V( I, I ), 1,
     $                        ZERO, T( 1, I ), 1 )
               ELSE
*
*                 T(1:i-1,i) := - tau(i) * V(1:i-1,i:n) * V(i,i:n)'
*
                  IF( I.LT.N )
     $               CALL ZLACGV( N-I, V( I, I+1 ), LDV )
                  CALL ZGEMV( 'No transpose', I-1, N-I+1, -TAU( I ),
     $                        V( 1, I ), LDV, V( I, I ), LDV, ZERO,
     $                        T( 1, I ), 1 )
                  IF( I.LT.N )
     $               CALL ZLACGV( N-I, V( I, I+1 ), LDV )
               END IF
               V( I, I ) = VII
*
*              T(1:i-1,i) := T(1:i-1,1:i-1) * T(1:i-1,i)
*
               CALL ZTRMV( 'Upper', 'No transpose', 'Non-unit', I-1, T,
     $                     LDT, T( 1, I ), 1 )
               T( I, I ) = TAU( I )
            END IF
   20    CONTINUE
      ELSE
         DO 40 I = K, 1, -1
            IF( TAU( I ).EQ.ZERO ) THEN
*
*              H(i)  =  I
*
               DO 30 J = I, K
                  T( J, I ) = ZERO
   30          CONTINUE
            ELSE
*
*              general case
*
               IF( I.LT.K ) THEN
                  IF( LSAME( STOREV, 'C' ) ) THEN
                     VII = V( N-K+I, I )
                     V( N-K+I, I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(1:n-k+i,i+1:k)' * V(1:n-k+i,i)
*
                     CALL ZGEMV( 'Conjugate transpose', N-K+I, K-I,
     $                           -TAU( I ), V( 1, I+1 ), LDV, V( 1, I ),
     $                           1, ZERO, T( I+1, I ), 1 )
                     V( N-K+I, I ) = VII
                  ELSE
                     VII = V( I, N-K+I )
                     V( I, N-K+I ) = ONE
*
*                    T(i+1:k,i) :=
*                            - tau(i) * V(i+1:k,1:n-k+i) * V(i,1:n-k+i)'
*
                     CALL ZLACGV( N-K+I-1, V( I, 1 ), LDV )
                     CALL ZGEMV( 'No transpose', K-I, N-K+I, -TAU( I ),
     $                           V( I+1, 1 ), LDV, V( I, 1 ), LDV, ZERO,
     $                           T( I+1, I ), 1 )
                     CALL ZLACGV( N-K+I-1, V( I, 1 ), LDV )
                     V( I, N-K+I ) = VII
                  END IF
*
*                 T(i+1:k,i) := T(i+1:k,i+1:k) * T(i+1:k,i)
*
                  CALL ZTRMV( 'Lower', 'No transpose', 'Non-unit', K-I,
     $                        T( I+1, I+1 ), LDT, T( I+1, I ), 1 )
               END IF
               T( I, I ) = TAU( I )
            END IF
   40    CONTINUE
      END IF
      RETURN
*
*     End of ZLARFT
*
      END
      SUBROUTINE ZLASCL( TYPE, KL, KU, CFROM, CTO, M, N, A, LDA, INFO )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          TYPE
      INTEGER            INFO, KL, KU, LDA, M, N
      DOUBLE PRECISION   CFROM, CTO
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASCL multiplies the M by N complex matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) DOUBLE PRECISION
*  CTO     (input) DOUBLE PRECISION
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,M)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            DONE
      INTEGER            I, ITYPE, J, K1, K2, K3, K4
      DOUBLE PRECISION   BIGNUM, CFROM1, CFROMC, CTO1, CTOC, MUL, SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
*
      IF( LSAME( TYPE, 'G' ) ) THEN
         ITYPE = 0
      ELSE IF( LSAME( TYPE, 'L' ) ) THEN
         ITYPE = 1
      ELSE IF( LSAME( TYPE, 'U' ) ) THEN
         ITYPE = 2
      ELSE IF( LSAME( TYPE, 'H' ) ) THEN
         ITYPE = 3
      ELSE IF( LSAME( TYPE, 'B' ) ) THEN
         ITYPE = 4
      ELSE IF( LSAME( TYPE, 'Q' ) ) THEN
         ITYPE = 5
      ELSE IF( LSAME( TYPE, 'Z' ) ) THEN
         ITYPE = 6
      ELSE
         ITYPE = -1
      END IF
*
      IF( ITYPE.EQ.-1 ) THEN
         INFO = -1
      ELSE IF( CFROM.EQ.ZERO ) THEN
         INFO = -4
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( N.LT.0 .OR. ( ITYPE.EQ.4 .AND. N.NE.M ) .OR.
     $         ( ITYPE.EQ.5 .AND. N.NE.M ) ) THEN
         INFO = -7
      ELSE IF( ITYPE.LE.3 .AND. LDA.LT.MAX( 1, M ) ) THEN
         INFO = -9
      ELSE IF( ITYPE.GE.4 ) THEN
         IF( KL.LT.0 .OR. KL.GT.MAX( M-1, 0 ) ) THEN
            INFO = -2
         ELSE IF( KU.LT.0 .OR. KU.GT.MAX( N-1, 0 ) .OR.
     $            ( ( ITYPE.EQ.4 .OR. ITYPE.EQ.5 ) .AND. KL.NE.KU ) )
     $             THEN
            INFO = -3
         ELSE IF( ( ITYPE.EQ.4 .AND. LDA.LT.KL+1 ) .OR.
     $            ( ITYPE.EQ.5 .AND. LDA.LT.KU+1 ) .OR.
     $            ( ITYPE.EQ.6 .AND. LDA.LT.2*KL+KU+1 ) ) THEN
            INFO = -9
         END IF
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASCL', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 .OR. M.EQ.0 )
     $   RETURN
*
*     Get machine parameters
*
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
*
      CFROMC = CFROM
      CTOC = CTO
*
   10 CONTINUE
      CFROM1 = CFROMC*SMLNUM
      CTO1 = CTOC / BIGNUM
      IF( ABS( CFROM1 ).GT.ABS( CTOC ) .AND. CTOC.NE.ZERO ) THEN
         MUL = SMLNUM
         DONE = .FALSE.
         CFROMC = CFROM1
      ELSE IF( ABS( CTO1 ).GT.ABS( CFROMC ) ) THEN
         MUL = BIGNUM
         DONE = .FALSE.
         CTOC = CTO1
      ELSE
         MUL = CTOC / CFROMC
         DONE = .TRUE.
      END IF
*
      IF( ITYPE.EQ.0 ) THEN
*
*        Full matrix
*
         DO 30 J = 1, N
            DO 20 I = 1, M
               A( I, J ) = A( I, J )*MUL
   20       CONTINUE
   30    CONTINUE
*
      ELSE IF( ITYPE.EQ.1 ) THEN
*
*        Lower triangular matrix
*
         DO 50 J = 1, N
            DO 40 I = J, M
               A( I, J ) = A( I, J )*MUL
   40       CONTINUE
   50    CONTINUE
*
      ELSE IF( ITYPE.EQ.2 ) THEN
*
*        Upper triangular matrix
*
         DO 70 J = 1, N
            DO 60 I = 1, MIN( J, M )
               A( I, J ) = A( I, J )*MUL
   60       CONTINUE
   70    CONTINUE
*
      ELSE IF( ITYPE.EQ.3 ) THEN
*
*        Upper Hessenberg matrix
*
         DO 90 J = 1, N
            DO 80 I = 1, MIN( J+1, M )
               A( I, J ) = A( I, J )*MUL
   80       CONTINUE
   90    CONTINUE
*
      ELSE IF( ITYPE.EQ.4 ) THEN
*
*        Lower half of a symmetric band matrix
*
         K3 = KL + 1
         K4 = N + 1
         DO 110 J = 1, N
            DO 100 I = 1, MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  100       CONTINUE
  110    CONTINUE
*
      ELSE IF( ITYPE.EQ.5 ) THEN
*
*        Upper half of a symmetric band matrix
*
         K1 = KU + 2
         K3 = KU + 1
         DO 130 J = 1, N
            DO 120 I = MAX( K1-J, 1 ), K3
               A( I, J ) = A( I, J )*MUL
  120       CONTINUE
  130    CONTINUE
*
      ELSE IF( ITYPE.EQ.6 ) THEN
*
*        Band matrix
*
         K1 = KL + KU + 2
         K2 = KL + 1
         K3 = 2*KL + KU + 1
         K4 = KL + KU + 1 + M
         DO 150 J = 1, N
            DO 140 I = MAX( K1-J, K2 ), MIN( K3, K4-J )
               A( I, J ) = A( I, J )*MUL
  140       CONTINUE
  150    CONTINUE
*
      END IF
*
      IF( .NOT.DONE )
     $   GO TO 10
*
      RETURN
*
*     End of ZLASCL
*
      END
      SUBROUTINE ZLASET( UPLO, M, N, ALPHA, BETA, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, M, N
      COMPLEX*16         ALPHA, BETA
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASET initializes a 2-D array A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set. The lower triangle
*                      is unchanged.
*          = 'L':      Lower triangular part is set. The upper triangle
*                      is unchanged.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          On entry, M specifies the number of rows of A.
*
*  N       (input) INTEGER
*          On entry, N specifies the number of columns of A.
*
*  ALPHA   (input) COMPLEX*16
*          All the offdiagonal array elements are set to ALPHA.
*
*  BETA    (input) COMPLEX*16
*          All the diagonal array elements are set to BETA.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, A(i,j) = ALPHA, 1 <= i <= m, 1 <= j <= n, i.ne.j;
*                   A(i,i) = BETA , 1 <= i <= min(m,n)
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Set the diagonal to BETA and the strictly upper triangular
*        part of the array to ALPHA.
*
         DO 20 J = 2, N
            DO 10 I = 1, MIN( J-1, M )
               A( I, J ) = ALPHA
   10       CONTINUE
   20    CONTINUE
         DO 30 I = 1, MIN( N, M )
            A( I, I ) = BETA
   30    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
*
*        Set the diagonal to BETA and the strictly lower triangular
*        part of the array to ALPHA.
*
         DO 50 J = 1, MIN( M, N )
            DO 40 I = J + 1, M
               A( I, J ) = ALPHA
   40       CONTINUE
   50    CONTINUE
         DO 60 I = 1, MIN( N, M )
            A( I, I ) = BETA
   60    CONTINUE
*
      ELSE
*
*        Set the array to BETA on the diagonal and ALPHA on the
*        offdiagonal.
*
         DO 80 J = 1, N
            DO 70 I = 1, M
               A( I, J ) = ALPHA
   70       CONTINUE
   80    CONTINUE
         DO 90 I = 1, MIN( M, N )
            A( I, I ) = BETA
   90    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASET
*
      END
      SUBROUTINE ZLASR( SIDE, PIVOT, DIRECT, M, N, C, S, A, LDA )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          DIRECT, PIVOT, SIDE
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( * ), S( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASR   performs the transformation
*
*     A := P*A,   when SIDE = 'L' or 'l'  (  Left-hand side )
*
*     A := A*P',  when SIDE = 'R' or 'r'  ( Right-hand side )
*
*  where A is an m by n complex matrix and P is an orthogonal matrix,
*  consisting of a sequence of plane rotations determined by the
*  parameters PIVOT and DIRECT as follows ( z = m when SIDE = 'L' or 'l'
*  and z = n when SIDE = 'R' or 'r' ):
*
*  When  DIRECT = 'F' or 'f'  ( Forward sequence ) then
*
*     P = P( z - 1 )*...*P( 2 )*P( 1 ),
*
*  and when DIRECT = 'B' or 'b'  ( Backward sequence ) then
*
*     P = P( 1 )*P( 2 )*...*P( z - 1 ),
*
*  where  P( k ) is a plane rotation matrix for the following planes:
*
*     when  PIVOT = 'V' or 'v'  ( Variable pivot ),
*        the plane ( k, k + 1 )
*
*     when  PIVOT = 'T' or 't'  ( Top pivot ),
*        the plane ( 1, k + 1 )
*
*     when  PIVOT = 'B' or 'b'  ( Bottom pivot ),
*        the plane ( k, z )
*
*  c( k ) and s( k )  must contain the  cosine and sine that define the
*  matrix  P( k ).  The two by two plane rotation part of the matrix
*  P( k ), R( k ), is assumed to be of the form
*
*     R( k ) = (  c( k )  s( k ) ).
*              ( -s( k )  c( k ) )
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          Specifies whether the plane rotation matrix P is applied to
*          A on the left or the right.
*          = 'L':  Left, compute A := P*A
*          = 'R':  Right, compute A:= A*P'
*
*  DIRECT  (input) CHARACTER*1
*          Specifies whether P is a forward or backward sequence of
*          plane rotations.
*          = 'F':  Forward, P = P( z - 1 )*...*P( 2 )*P( 1 )
*          = 'B':  Backward, P = P( 1 )*P( 2 )*...*P( z - 1 )
*
*  PIVOT   (input) CHARACTER*1
*          Specifies the plane for which P(k) is a plane rotation
*          matrix.
*          = 'V':  Variable pivot, the plane (k,k+1)
*          = 'T':  Top pivot, the plane (1,k+1)
*          = 'B':  Bottom pivot, the plane (k,z)
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  If m <= 1, an immediate
*          return is effected.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  If n <= 1, an
*          immediate return is effected.
*
*  C, S    (input) DOUBLE PRECISION arrays, dimension
*                  (M-1) if SIDE = 'L'
*                  (N-1) if SIDE = 'R'
*          c(k) and s(k) contain the cosine and sine that define the
*          matrix P(k).  The two by two plane rotation part of the
*          matrix P(k), R(k), is assumed to be of the form
*          R( k ) = (  c( k )  s( k ) ).
*                   ( -s( k )  c( k ) )
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          The m by n matrix A.  On exit, A is overwritten by P*A if
*          SIDE = 'R' or by A*P' if SIDE = 'L'.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, INFO, J
      DOUBLE PRECISION   CTEMP, STEMP
      COMPLEX*16         TEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.( LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'R' ) ) ) THEN
         INFO = 1
      ELSE IF( .NOT.( LSAME( PIVOT, 'V' ) .OR. LSAME( PIVOT,
     $         'T' ) .OR. LSAME( PIVOT, 'B' ) ) ) THEN
         INFO = 2
      ELSE IF( .NOT.( LSAME( DIRECT, 'F' ) .OR. LSAME( DIRECT, 'B' ) ) )
     $          THEN
         INFO = 3
      ELSE IF( M.LT.0 ) THEN
         INFO = 4
      ELSE IF( N.LT.0 ) THEN
         INFO = 5
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = 9
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZLASR ', INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( ( M.EQ.0 ) .OR. ( N.EQ.0 ) )
     $   RETURN
      IF( LSAME( SIDE, 'L' ) ) THEN
*
*        Form  P * A
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 20 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 10 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   10                CONTINUE
                  END IF
   20          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 40 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 30 I = 1, N
                        TEMP = A( J+1, I )
                        A( J+1, I ) = CTEMP*TEMP - STEMP*A( J, I )
                        A( J, I ) = STEMP*TEMP + CTEMP*A( J, I )
   30                CONTINUE
                  END IF
   40          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 60 J = 2, M
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 50 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   50                CONTINUE
                  END IF
   60          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 80 J = M, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 70 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = CTEMP*TEMP - STEMP*A( 1, I )
                        A( 1, I ) = STEMP*TEMP + CTEMP*A( 1, I )
   70                CONTINUE
                  END IF
   80          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 100 J = 1, M - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 90 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
   90                CONTINUE
                  END IF
  100          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 120 J = M - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 110 I = 1, N
                        TEMP = A( J, I )
                        A( J, I ) = STEMP*A( M, I ) + CTEMP*TEMP
                        A( M, I ) = CTEMP*A( M, I ) - STEMP*TEMP
  110                CONTINUE
                  END IF
  120          CONTINUE
            END IF
         END IF
      ELSE IF( LSAME( SIDE, 'R' ) ) THEN
*
*        Form A * P'
*
         IF( LSAME( PIVOT, 'V' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 140 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 130 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  130                CONTINUE
                  END IF
  140          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 160 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 150 I = 1, M
                        TEMP = A( I, J+1 )
                        A( I, J+1 ) = CTEMP*TEMP - STEMP*A( I, J )
                        A( I, J ) = STEMP*TEMP + CTEMP*A( I, J )
  150                CONTINUE
                  END IF
  160          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'T' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 180 J = 2, N
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 170 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  170                CONTINUE
                  END IF
  180          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 200 J = N, 2, -1
                  CTEMP = C( J-1 )
                  STEMP = S( J-1 )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 190 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = CTEMP*TEMP - STEMP*A( I, 1 )
                        A( I, 1 ) = STEMP*TEMP + CTEMP*A( I, 1 )
  190                CONTINUE
                  END IF
  200          CONTINUE
            END IF
         ELSE IF( LSAME( PIVOT, 'B' ) ) THEN
            IF( LSAME( DIRECT, 'F' ) ) THEN
               DO 220 J = 1, N - 1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 210 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  210                CONTINUE
                  END IF
  220          CONTINUE
            ELSE IF( LSAME( DIRECT, 'B' ) ) THEN
               DO 240 J = N - 1, 1, -1
                  CTEMP = C( J )
                  STEMP = S( J )
                  IF( ( CTEMP.NE.ONE ) .OR. ( STEMP.NE.ZERO ) ) THEN
                     DO 230 I = 1, M
                        TEMP = A( I, J )
                        A( I, J ) = STEMP*A( I, N ) + CTEMP*TEMP
                        A( I, N ) = CTEMP*A( I, N ) - STEMP*TEMP
  230                CONTINUE
                  END IF
  240          CONTINUE
            END IF
         END IF
      END IF
*
      RETURN
*
*     End of ZLASR
*
      END
      SUBROUTINE ZLASSQ( N, X, INCX, SCALE, SUMSQ )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      DOUBLE PRECISION   SCALE, SUMSQ
*     ..
*     .. Array Arguments ..
      COMPLEX*16         X( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASSQ returns the values scl and ssq such that
*
*     ( scl**2 )*ssq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where x( i ) = abs( X( 1 + ( i - 1 )*INCX ) ). The value of sumsq is
*  assumed to be at least unity and the value of ssq will then satisfy
*
*     1.0 .le. ssq .le. ( sumsq + 2*n ).
*
*  scale is assumed to be non-negative and scl returns the value
*
*     scl = max( scale, abs( real( x( i ) ) ), abs( aimag( x( i ) ) ) ),
*            i
*
*  scale and sumsq must be supplied in SCALE and SUMSQ respectively.
*  SCALE and SUMSQ are overwritten by scl and ssq respectively.
*
*  The routine makes only one pass through the vector X.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) COMPLEX*16 array, dimension (N)
*          The vector x as described above.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with the value  scl .
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with the value  ssq .
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            IX
      DOUBLE PRECISION   TEMP1
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG
*     ..
*     .. Executable Statements ..
*
      IF( N.GT.0 ) THEN
         DO 10 IX = 1, 1 + ( N-1 )*INCX, INCX
            IF( DBLE( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DBLE( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
            IF( DIMAG( X( IX ) ).NE.ZERO ) THEN
               TEMP1 = ABS( DIMAG( X( IX ) ) )
               IF( SCALE.LT.TEMP1 ) THEN
                  SUMSQ = 1 + SUMSQ*( SCALE / TEMP1 )**2
                  SCALE = TEMP1
               ELSE
                  SUMSQ = SUMSQ + ( TEMP1 / SCALE )**2
               END IF
            END IF
   10    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASSQ
*
      END
      SUBROUTINE ZLATRD( UPLO, N, NB, A, LDA, E, TAU, W, LDW )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDW, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   E( * )
      COMPLEX*16         A( LDA, * ), TAU( * ), W( LDW, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLATRD reduces NB rows and columns of a complex Hermitian matrix A to
*  Hermitian tridiagonal form by a unitary similarity
*  transformation Q' * A * Q, and returns the matrices V and W which are
*  needed to apply the transformation to the unreduced part of A.
*
*  If UPLO = 'U', ZLATRD reduces the last NB rows and columns of a
*  matrix, of which the upper triangle is supplied;
*  if UPLO = 'L', ZLATRD reduces the first NB rows and columns of a
*  matrix, of which the lower triangle is supplied.
*
*  This is an auxiliary routine called by ZHETRD.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER
*          Specifies whether the upper or lower triangular part of the
*          Hermitian matrix A is stored:
*          = 'U': Upper triangular
*          = 'L': Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  NB      (input) INTEGER
*          The number of rows and columns to be reduced.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          n-by-n upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n-by-n lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*          On exit:
*          if UPLO = 'U', the last NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements above the diagonal
*            with the array TAU, represent the unitary matrix Q as a
*            product of elementary reflectors;
*          if UPLO = 'L', the first NB columns have been reduced to
*            tridiagonal form, with the diagonal elements overwriting
*            the diagonal elements of A; the elements below the diagonal
*            with the array TAU, represent the  unitary matrix Q as a
*            product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          If UPLO = 'U', E(n-nb:n-1) contains the superdiagonal
*          elements of the last NB columns of the reduced matrix;
*          if UPLO = 'L', E(1:nb) contains the subdiagonal elements of
*          the first NB columns of the reduced matrix.
*
*  TAU     (output) COMPLEX*16 array, dimension (N-1)
*          The scalar factors of the elementary reflectors, stored in
*          TAU(n-nb:n-1) if UPLO = 'U', and in TAU(1:nb) if UPLO = 'L'.
*          See Further Details.
*
*  W       (output) COMPLEX*16 array, dimension (LDW,NB)
*          The n-by-nb matrix W required to update the unreduced part
*          of A.
*
*  LDW     (input) INTEGER
*          The leading dimension of the array W. LDW >= max(1,N).
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(n) H(n-1) . . . H(n-nb+1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(i:n) = 0 and v(i-1) = 1; v(1:i-1) is stored on exit in A(1:i-1,i),
*  and tau in TAU(i-1).
*
*  If UPLO = 'L', the matrix Q is represented as a product of elementary
*  reflectors
*
*     Q = H(1) H(2) . . . H(nb).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i) = 0 and v(i+1) = 1; v(i+1:n) is stored on exit in A(i+1:n,i),
*  and tau in TAU(i).
*
*  The elements of the vectors v together form the n-by-nb matrix V
*  which is needed, with W, to apply the transformation to the unreduced
*  part of the matrix, using a Hermitian rank-2k update of the form:
*  A := A - V*W' - W*V'.
*
*  The contents of A on exit are illustrated by the following examples
*  with n = 5 and nb = 2:
*
*  if UPLO = 'U':                       if UPLO = 'L':
*
*    (  a   a   a   v4  v5 )              (  d                  )
*    (      a   a   v4  v5 )              (  1   d              )
*    (          a   1   v5 )              (  v1  1   a          )
*    (              d   1  )              (  v1  v2  a   a      )
*    (                  d  )              (  v1  v2  a   a   a  )
*
*  where d denotes a diagonal element of the reduced matrix, a denotes
*  an element of the original matrix that is unchanged, and vi denotes
*  an element of the vector defining H(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE, HALF
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ),
     $                   HALF = ( 0.5D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IW
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZAXPY, ZGEMV, ZHEMV, ZLACGV, ZLARFG, ZSCAL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      COMPLEX*16         ZDOTC
      EXTERNAL           LSAME, ZDOTC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
      IF( LSAME( UPLO, 'U' ) ) THEN
*
*        Reduce last NB columns of upper triangle
*
         DO 10 I = N, N - NB + 1, -1
            IW = I - N + NB
            IF( I.LT.N ) THEN
*
*              Update A(1:i,i)
*
               A( I, I ) = DBLE( A( I, I ) )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, A( 1, I+1 ),
     $                     LDA, W( I, IW+1 ), LDW, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, W( I, IW+1 ), LDW )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               CALL ZGEMV( 'No transpose', I, N-I, -ONE, W( 1, IW+1 ),
     $                     LDW, A( I, I+1 ), LDA, ONE, A( 1, I ), 1 )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I ) = DBLE( A( I, I ) )
            END IF
            IF( I.GT.1 ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(1:i-2,i)
*
               ALPHA = A( I-1, I )
               CALL ZLARFG( I-1, ALPHA, A( 1, I ), 1, TAU( I-1 ) )
               E( I-1 ) = ALPHA
               A( I-1, I ) = ONE
*
*              Compute W(1:i-1,i)
*
               CALL ZHEMV( 'Upper', I-1, ONE, A, LDA, A( 1, I ), 1,
     $                     ZERO, W( 1, IW ), 1 )
               IF( I.LT.N ) THEN
                  CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE,
     $                        W( 1, IW+1 ), LDW, A( 1, I ), 1, ZERO,
     $                        W( I+1, IW ), 1 )
                  CALL ZGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        A( 1, I+1 ), LDA, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
                  CALL ZGEMV( 'Conjugate transpose', I-1, N-I, ONE,
     $                        A( 1, I+1 ), LDA, A( 1, I ), 1, ZERO,
     $                        W( I+1, IW ), 1 )
                  CALL ZGEMV( 'No transpose', I-1, N-I, -ONE,
     $                        W( 1, IW+1 ), LDW, W( I+1, IW ), 1, ONE,
     $                        W( 1, IW ), 1 )
               END IF
               CALL ZSCAL( I-1, TAU( I-1 ), W( 1, IW ), 1 )
               ALPHA = -HALF*TAU( I-1 )*ZDOTC( I-1, W( 1, IW ), 1,
     $                 A( 1, I ), 1 )
               CALL ZAXPY( I-1, ALPHA, A( 1, I ), 1, W( 1, IW ), 1 )
            END IF
*
   10    CONTINUE
      ELSE
*
*        Reduce first NB columns of lower triangle
*
         DO 20 I = 1, NB
*
*           Update A(i:n,i)
*
            A( I, I ) = DBLE( A( I, I ) )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, W( I, 1 ), LDW, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, W( I, 1 ), LDW )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, W( I, 1 ),
     $                  LDW, A( I, 1 ), LDA, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            A( I, I ) = DBLE( A( I, I ) )
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:n,i)
*
               ALPHA = A( I+1, I )
               CALL ZLARFG( N-I, ALPHA, A( MIN( I+2, N ), I ), 1,
     $                      TAU( I ) )
               E( I ) = ALPHA
               A( I+1, I ) = ONE
*
*              Compute W(i+1:n,i)
*
               CALL ZHEMV( 'Lower', N-I, ONE, A( I+1, I+1 ), LDA,
     $                     A( I+1, I ), 1, ZERO, W( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE,
     $                     W( I+1, 1 ), LDW, A( I+1, I ), 1, ZERO,
     $                     W( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I, I-1, ONE,
     $                     A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO,
     $                     W( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, W( I+1, 1 ),
     $                     LDW, W( 1, I ), 1, ONE, W( I+1, I ), 1 )
               CALL ZSCAL( N-I, TAU( I ), W( I+1, I ), 1 )
               ALPHA = -HALF*TAU( I )*ZDOTC( N-I, W( I+1, I ), 1,
     $                 A( I+1, I ), 1 )
               CALL ZAXPY( N-I, ALPHA, A( I+1, I ), 1, W( I+1, I ), 1 )
            END IF
*
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLATRD
*
      END
      SUBROUTINE ZSTEQR( COMPZ, N, D, E, Z, LDZ, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          COMPZ
      INTEGER            INFO, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), WORK( * )
      COMPLEX*16         Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZSTEQR computes all eigenvalues and, optionally, eigenvectors of a
*  symmetric tridiagonal matrix using the implicit QL or QR method.
*  The eigenvectors of a full or band complex Hermitian matrix can also
*  be found if ZHETRD or ZHPTRD or ZHBTRD has been used to reduce this
*  matrix to tridiagonal form.
*
*  Arguments
*  =========
*
*  COMPZ   (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only.
*          = 'V':  Compute eigenvalues and eigenvectors of the original
*                  Hermitian matrix.  On entry, Z must contain the
*                  unitary matrix used to reduce the original matrix
*                  to tridiagonal form.
*          = 'I':  Compute eigenvalues and eigenvectors of the
*                  tridiagonal matrix.  Z is initialized to the identity
*                  matrix.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the diagonal elements of the tridiagonal matrix.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix.
*          On exit, E has been destroyed.
*
*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
*          On entry, if  COMPZ = 'V', then Z contains the unitary
*          matrix used in the reduction to tridiagonal form.
*          On exit, if INFO = 0, then if COMPZ = 'V', Z contains the
*          orthonormal eigenvectors of the original Hermitian matrix,
*          and if COMPZ = 'I', Z contains the orthonormal eigenvectors
*          of the symmetric tridiagonal matrix.
*          If COMPZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          eigenvectors are desired, then  LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,2*N-2))
*          If COMPZ = 'N', then WORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm has failed to find all the eigenvalues in
*                a total of 30*N iterations; if INFO = i, then i
*                elements of E have not converged to zero; on exit, D
*                and E contain the elements of a symmetric tridiagonal
*                matrix which is unitarily similar to the original
*                matrix.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, TWO, THREE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0, TWO = 2.0D0,
     $                   THREE = 3.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ),
     $                   CONE = ( 1.0D0, 0.0D0 ) )
      INTEGER            MAXIT
      PARAMETER          ( MAXIT = 30 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICOMPZ, II, ISCALE, J, JTOT, K, L, L1, LEND,
     $                   LENDM1, LENDP1, LENDSV, LM1, LSV, M, MM, MM1,
     $                   NM1, NMAXIT
      DOUBLE PRECISION   ANORM, B, C, EPS, EPS2, F, G, P, R, RT1, RT2,
     $                   S, SAFMAX, SAFMIN, SSFMAX, SSFMIN, TST
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, DLANST, DLAPY2
      EXTERNAL           LSAME, DLAMCH, DLANST, DLAPY2
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLAE2, DLAEV2, DLARTG, DLASCL, DLASRT, XERBLA,
     $                   ZLASET, ZLASR, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ICOMPZ = 0
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ICOMPZ = 2
      ELSE
         ICOMPZ = -1
      END IF
      IF( ICOMPZ.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( ( LDZ.LT.1 ) .OR. ( ICOMPZ.GT.0 .AND. LDZ.LT.MAX( 1,
     $         N ) ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZSTEQR', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( ICOMPZ.EQ.2 )
     $      Z( 1, 1 ) = CONE
         RETURN
      END IF
*
*     Determine the unit roundoff and over/underflow thresholds.
*
      EPS = DLAMCH( 'E' )
      EPS2 = EPS**2
      SAFMIN = DLAMCH( 'S' )
      SAFMAX = ONE / SAFMIN
      SSFMAX = SQRT( SAFMAX ) / THREE
      SSFMIN = SQRT( SAFMIN ) / EPS2
*
*     Compute the eigenvalues and eigenvectors of the tridiagonal
*     matrix.
*
      IF( ICOMPZ.EQ.2 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
      NMAXIT = N*MAXIT
      JTOT = 0
*
*     Determine where the matrix splits and choose QL or QR iteration
*     for each block, according to whether top or bottom diagonal
*     element is smaller.
*
      L1 = 1
      NM1 = N - 1
*
   10 CONTINUE
      IF( L1.GT.N )
     $   GO TO 160
      IF( L1.GT.1 )
     $   E( L1-1 ) = ZERO
      IF( L1.LE.NM1 ) THEN
         DO 20 M = L1, NM1
            TST = ABS( E( M ) )
            IF( TST.EQ.ZERO )
     $         GO TO 30
            IF( TST.LE.( SQRT( ABS( D( M ) ) )*SQRT( ABS( D( M+
     $          1 ) ) ) )*EPS ) THEN
               E( M ) = ZERO
               GO TO 30
            END IF
   20    CONTINUE
      END IF
      M = N
*
   30 CONTINUE
      L = L1
      LSV = L
      LEND = M
      LENDSV = LEND
      L1 = M + 1
      IF( LEND.EQ.L )
     $   GO TO 10
*
*     Scale submatrix in rows and columns L to LEND
*
      ANORM = DLANST( 'I', LEND-L+1, D( L ), E( L ) )
      ISCALE = 0
      IF( ANORM.EQ.ZERO )
     $   GO TO 10
      IF( ANORM.GT.SSFMAX ) THEN
         ISCALE = 1
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMAX, LEND-L, 1, E( L ), N,
     $                INFO )
      ELSE IF( ANORM.LT.SSFMIN ) THEN
         ISCALE = 2
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L+1, 1, D( L ), N,
     $                INFO )
         CALL DLASCL( 'G', 0, 0, ANORM, SSFMIN, LEND-L, 1, E( L ), N,
     $                INFO )
      END IF
*
*     Choose between QL and QR iteration
*
      IF( ABS( D( LEND ) ).LT.ABS( D( L ) ) ) THEN
         LEND = LSV
         L = LENDSV
      END IF
*
      IF( LEND.GT.L ) THEN
*
*        QL Iteration
*
*        Look for small subdiagonal element.
*
   40    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDM1 = LEND - 1
            DO 50 M = L, LENDM1
               TST = ABS( E( M ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M+1 ) )+
     $             SAFMIN )GO TO 60
   50       CONTINUE
         END IF
*
         M = LEND
*
   60    CONTINUE
         IF( M.LT.LEND )
     $      E( M ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 80
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L+1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L ), E( L ), D( L+1 ), RT1, RT2, C, S )
               WORK( L ) = C
               WORK( N-1+L ) = S
               CALL ZLASR( 'R', 'V', 'B', N, 2, WORK( L ),
     $                     WORK( N-1+L ), Z( 1, L ), LDZ )
            ELSE
               CALL DLAE2( D( L ), E( L ), D( L+1 ), RT1, RT2 )
            END IF
            D( L ) = RT1
            D( L+1 ) = RT2
            E( L ) = ZERO
            L = L + 2
            IF( L.LE.LEND )
     $         GO TO 40
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L+1 )-P ) / ( TWO*E( L ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         MM1 = M - 1
         DO 70 I = MM1, L, -1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M-1 )
     $         E( I+1 ) = R
            G = D( I+1 ) - P
            R = ( D( I )-G )*S + TWO*C*B
            P = S*R
            D( I+1 ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = -S
            END IF
*
   70    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = M - L + 1
            CALL ZLASR( 'R', 'V', 'B', N, MM, WORK( L ), WORK( N-1+L ),
     $                  Z( 1, L ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( L ) = G
         GO TO 40
*
*        Eigenvalue found.
*
   80    CONTINUE
         D( L ) = P
*
         L = L + 1
         IF( L.LE.LEND )
     $      GO TO 40
         GO TO 140
*
      ELSE
*
*        QR Iteration
*
*        Look for small superdiagonal element.
*
   90    CONTINUE
         IF( L.NE.LEND ) THEN
            LENDP1 = LEND + 1
            DO 100 M = L, LENDP1, -1
               TST = ABS( E( M-1 ) )**2
               IF( TST.LE.( EPS2*ABS( D( M ) ) )*ABS( D( M-1 ) )+
     $             SAFMIN )GO TO 110
  100       CONTINUE
         END IF
*
         M = LEND
*
  110    CONTINUE
         IF( M.GT.LEND )
     $      E( M-1 ) = ZERO
         P = D( L )
         IF( M.EQ.L )
     $      GO TO 130
*
*        If remaining matrix is 2-by-2, use DLAE2 or SLAEV2
*        to compute its eigensystem.
*
         IF( M.EQ.L-1 ) THEN
            IF( ICOMPZ.GT.0 ) THEN
               CALL DLAEV2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2, C, S )
               WORK( M ) = C
               WORK( N-1+M ) = S
               CALL ZLASR( 'R', 'V', 'F', N, 2, WORK( M ),
     $                     WORK( N-1+M ), Z( 1, L-1 ), LDZ )
            ELSE
               CALL DLAE2( D( L-1 ), E( L-1 ), D( L ), RT1, RT2 )
            END IF
            D( L-1 ) = RT1
            D( L ) = RT2
            E( L-1 ) = ZERO
            L = L - 2
            IF( L.GE.LEND )
     $         GO TO 90
            GO TO 140
         END IF
*
         IF( JTOT.EQ.NMAXIT )
     $      GO TO 140
         JTOT = JTOT + 1
*
*        Form shift.
*
         G = ( D( L-1 )-P ) / ( TWO*E( L-1 ) )
         R = DLAPY2( G, ONE )
         G = D( M ) - P + ( E( L-1 ) / ( G+SIGN( R, G ) ) )
*
         S = ONE
         C = ONE
         P = ZERO
*
*        Inner loop
*
         LM1 = L - 1
         DO 120 I = M, LM1
            F = S*E( I )
            B = C*E( I )
            CALL DLARTG( G, F, C, S, R )
            IF( I.NE.M )
     $         E( I-1 ) = R
            G = D( I ) - P
            R = ( D( I+1 )-G )*S + TWO*C*B
            P = S*R
            D( I ) = G + P
            G = C*R - B
*
*           If eigenvectors are desired, then save rotations.
*
            IF( ICOMPZ.GT.0 ) THEN
               WORK( I ) = C
               WORK( N-1+I ) = S
            END IF
*
  120    CONTINUE
*
*        If eigenvectors are desired, then apply saved rotations.
*
         IF( ICOMPZ.GT.0 ) THEN
            MM = L - M + 1
            CALL ZLASR( 'R', 'V', 'F', N, MM, WORK( M ), WORK( N-1+M ),
     $                  Z( 1, M ), LDZ )
         END IF
*
         D( L ) = D( L ) - P
         E( LM1 ) = G
         GO TO 90
*
*        Eigenvalue found.
*
  130    CONTINUE
         D( L ) = P
*
         L = L - 1
         IF( L.GE.LEND )
     $      GO TO 90
         GO TO 140
*
      END IF
*
*     Undo scaling if necessary
*
  140 CONTINUE
      IF( ISCALE.EQ.1 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMAX, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      ELSE IF( ISCALE.EQ.2 ) THEN
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV+1, 1,
     $                D( LSV ), N, INFO )
         CALL DLASCL( 'G', 0, 0, SSFMIN, ANORM, LENDSV-LSV, 1, E( LSV ),
     $                N, INFO )
      END IF
*
*     Check for no convergence to an eigenvalue after a total
*     of N*MAXIT iterations.
*
      IF( JTOT.EQ.NMAXIT ) THEN
         DO 150 I = 1, N - 1
            IF( E( I ).NE.ZERO )
     $         INFO = INFO + 1
  150    CONTINUE
         RETURN
      END IF
      GO TO 10
*
*     Order eigenvalues and eigenvectors.
*
  160 CONTINUE
      IF( ICOMPZ.EQ.0 ) THEN
*
*        Use Quick Sort
*
         CALL DLASRT( 'I', N, D, INFO )
*
      ELSE
*
*        Use Selection Sort to minimize swaps of eigenvectors
*
         DO 180 II = 2, N
            I = II - 1
            K = I
            P = D( I )
            DO 170 J = II, N
               IF( D( J ).LT.P ) THEN
                  K = J
                  P = D( J )
               END IF
  170       CONTINUE
            IF( K.NE.I ) THEN
               D( K ) = D( I )
               D( I ) = P
               CALL ZSWAP( N, Z( 1, I ), 1, Z( 1, K ), 1 )
            END IF
  180    CONTINUE
      END IF
      RETURN
*
*     End of ZSTEQR
*
      END
      SUBROUTINE ZUNG2L( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNG2L generates an m by n complex matrix Q with orthonormal columns,
*  which is defined as the last n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by ZGEQLF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the (n-k+i)-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by ZGEQLF in the last k columns of its array
*          argument A.
*          On exit, the m-by-n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQLF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, II, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2L', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns 1:n-k to columns of the unit matrix
*
      DO 20 J = 1, N - K
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( M-N+J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = 1, K
         II = N - K + I
*
*        Apply H(i) to A(1:m-k+i,1:n-k+i) from the left
*
         A( M-N+II, II ) = ONE
         CALL ZLARF( 'Left', M-N+II, II-1, A( 1, II ), 1, TAU( I ), A,
     $               LDA, WORK )
         CALL ZSCAL( M-N+II-1, -TAU( I ), A( 1, II ), 1 )
         A( M-N+II, II ) = ONE - TAU( I )
*
*        Set A(m-k+i+1:m,n-k+i) to zero
*
         DO 30 L = M - N + II + 1, M
            A( L, II ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNG2L
*
      END
      SUBROUTINE ZUNG2R( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNG2R generates an m by n complex matrix Q with orthonormal columns,
*  which is defined as the first n columns of a product of k elementary
*  reflectors of order m
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by ZGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by ZGEQRF in the first k columns of its array
*          argument A.
*          On exit, the m by n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQRF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNG2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 )
     $   RETURN
*
*     Initialise columns k+1:n to columns of the unit matrix
*
      DO 20 J = K + 1, N
         DO 10 L = 1, M
            A( L, J ) = ZERO
   10    CONTINUE
         A( J, J ) = ONE
   20 CONTINUE
*
      DO 40 I = K, 1, -1
*
*        Apply H(i) to A(i:m,i:n) from the left
*
         IF( I.LT.N ) THEN
            A( I, I ) = ONE
            CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1, TAU( I ),
     $                  A( I, I+1 ), LDA, WORK )
         END IF
         IF( I.LT.M )
     $      CALL ZSCAL( M-I, -TAU( I ), A( I+1, I ), 1 )
         A( I, I ) = ONE - TAU( I )
*
*        Set A(1:i-1,i) to zero
*
         DO 30 L = 1, I - 1
            A( L, I ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNG2R
*
      END
      SUBROUTINE ZUNGQL( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNGQL generates an M-by-N complex matrix Q with orthonormal columns,
*  which is defined as the last N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(k) . . . H(2) H(1)
*
*  as returned by ZGEQLF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the (n-k+i)-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by ZGEQLF in the last k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQLF.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KK, L, LDWORK, LWKOPT,
     $                   NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNG2L
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'ZUNGQL', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGQL', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'ZUNGQL', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGQL', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the first block.
*        The last kk columns are handled by the block method.
*
         KK = MIN( K, ( ( K-NX+NB-1 ) / NB )*NB )
*
*        Set A(m-kk+1:m,1:n-kk) to zero.
*
         DO 20 J = 1, N - KK
            DO 10 I = M - KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the first or only block.
*
      CALL ZUNG2L( M-KK, N-KK, K-KK, A, LDA, TAU, WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = K - KK + 1, K, NB
            IB = MIN( NB, K-I+1 )
            IF( N-K+I.GT.1 ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i+ib-1) . . . H(i+1) H(i)
*
               CALL ZLARFT( 'Backward', 'Columnwise', M-K+I+IB-1, IB,
     $                      A( 1, N-K+I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(1:m-k+i+ib-1,1:n-k+i-1) from the left
*
               CALL ZLARFB( 'Left', 'No transpose', 'Backward',
     $                      'Columnwise', M-K+I+IB-1, N-K+I-1, IB,
     $                      A( 1, N-K+I ), LDA, WORK, LDWORK, A, LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows 1:m-k+i+ib-1 of current block
*
            CALL ZUNG2L( M-K+I+IB-1, IB, IB, A( 1, N-K+I ), LDA,
     $                   TAU( I ), WORK, IINFO )
*
*           Set rows m-k+i+ib:m of current block to zero
*
            DO 40 J = N - K + I, N - K + I + IB - 1
               DO 30 L = M - K + I + IB, M
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZUNGQL
*
      END
      SUBROUTINE ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNGQR generates an M-by-N complex matrix Q with orthonormal columns,
*  which is defined as the first N columns of a product of K elementary
*  reflectors of order M
*
*        Q  =  H(1) H(2) . . . H(k)
*
*  as returned by ZGEQRF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by ZGEQRF in the first k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQRF.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNG2R
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'ZUNGQR', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, N )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 .OR. N.GT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.N ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'ZUNGQR', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGQR', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk columns are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(1:kk,kk+1:n) to zero.
*
         DO 20 J = KK + 1, N
            DO 10 I = 1, KK
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.N )
     $   CALL ZUNG2R( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL ZLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i:m,i+ib:n) from the left
*
               CALL ZLARFB( 'Left', 'No transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H to rows i:m of current block
*
            CALL ZUNG2R( M-I+1, IB, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set rows 1:i-1 of current block to zero
*
            DO 40 J = I, I + IB - 1
               DO 30 L = 1, I - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZUNGQR
*
      END
      SUBROUTINE ZUNGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNGTR generates a complex unitary matrix Q which is defined as the
*  product of n-1 elementary reflectors of order N, as returned by
*  ZHETRD:
*
*  if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
*
*  if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U': Upper triangle of A contains elementary reflectors
*                 from ZHETRD;
*          = 'L': Lower triangle of A contains elementary reflectors
*                 from ZHETRD.
*
*  N       (input) INTEGER
*          The order of the matrix Q. N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the vectors which define the elementary reflectors,
*          as returned by ZHETRD.
*          On exit, the N-by-N unitary matrix Q.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= N.
*
*  TAU     (input) COMPLEX*16 array, dimension (N-1)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZHETRD.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= N-1.
*          For optimum performance LWORK >= (N-1)*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, UPPER
      INTEGER            I, IINFO, J, LWKOPT, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNGQL, ZUNGQR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LQUERY = ( LWORK.EQ.-1 )
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( UPPER ) THEN
            NB = ILAENV( 1, 'ZUNGQL', ' ', N-1, N-1, N-1, -1 )
         ELSE
            NB = ILAENV( 1, 'ZUNGQR', ' ', N-1, N-1, N-1, -1 )
         END IF
         LWKOPT = MAX( 1, N-1 )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGTR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Q was determined by a call to ZHETRD with UPLO = 'U'
*
*        Shift the vectors which define the elementary reflectors one
*        column to the left, and set the last row and column of Q to
*        those of the unit matrix
*
         DO 20 J = 1, N - 1
            DO 10 I = 1, J - 1
               A( I, J ) = A( I, J+1 )
   10       CONTINUE
            A( N, J ) = ZERO
   20    CONTINUE
         DO 30 I = 1, N - 1
            A( I, N ) = ZERO
   30    CONTINUE
         A( N, N ) = ONE
*
*        Generate Q(1:n-1,1:n-1)
*
         CALL ZUNGQL( N-1, N-1, N-1, A, LDA, TAU, WORK, LWORK, IINFO )
*
      ELSE
*
*        Q was determined by a call to ZHETRD with UPLO = 'L'.
*
*        Shift the vectors which define the elementary reflectors one
*        column to the right, and set the first row and column of Q to
*        those of the unit matrix
*
         DO 50 J = N, 2, -1
            A( 1, J ) = ZERO
            DO 40 I = J + 1, N
               A( I, J ) = A( I, J-1 )
   40       CONTINUE
   50    CONTINUE
         A( 1, 1 ) = ONE
         DO 60 I = 2, N
            A( I, 1 ) = ZERO
   60    CONTINUE
         IF( N.GT.1 ) THEN
*
*           Generate Q(2:n,2:n)
*
            CALL ZUNGQR( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                   LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of ZUNGTR
*
      END
      SUBROUTINE DLADIV( A, B, C, D, P, Q )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   A, B, C, D, P, Q
*     ..
*
*  Purpose
*  =======
*
*  DLADIV performs complex division in  real arithmetic
*
*                        a + i*b
*             p + i*q = ---------
*                        c + i*d
*
*  The algorithm is due to Robert L. Smith and can be found
*  in D. Knuth, The art of Computer Programming, Vol.2, p.195
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*  B       (input) DOUBLE PRECISION
*  C       (input) DOUBLE PRECISION
*  D       (input) DOUBLE PRECISION
*          The scalars a, b, c, and d in the above expression.
*
*  P       (output) DOUBLE PRECISION
*  Q       (output) DOUBLE PRECISION
*          The scalars p and q in the above expression.
*
*  =====================================================================
*
*     .. Local Scalars ..
      DOUBLE PRECISION   E, F
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS
*     ..
*     .. Executable Statements ..
*
      IF( ABS( D ).LT.ABS( C ) ) THEN
         E = D / C
         F = C + D*E
         P = ( A+B*E ) / F
         Q = ( B-A*E ) / F
      ELSE
         E = C / D
         F = D + C*E
         P = ( B+A*E ) / F
         Q = ( -A+B*E ) / F
      END IF
*
      RETURN
*
*     End of DLADIV
*
      END
      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), W( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DSYEV computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,3*N-1).
*          For optimal efficiency, LWORK >= (NB+2)*N,
*          where NB is the blocksize for DSYTRD returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, LQUERY, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDTAU, INDWRK, ISCALE,
     $                   LLWORK, LOPT, LWKOPT, NB
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DLANSY
      EXTERNAL           LSAME, ILAENV, DLAMCH, DLANSY
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASCL, DORGTR, DSCAL, DSTEQR, DSTERF, DSYTRD,
     $                   XERBLA
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
      LQUERY = ( LWORK.EQ.-1 )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, 3*N-1 ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
*
      IF( INFO.EQ.0 ) THEN
         NB = ILAENV( 1, 'DSYTRD', UPLO, N, -1, -1, -1 )
         LWKOPT = MAX( 1, ( NB+2 )*N )
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DSYEV ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( N.EQ.1 ) THEN
         W( 1 ) = A( 1, 1 )
         WORK( 1 ) = 3
         IF( WANTZ )
     $      A( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = DLANSY( 'M', UPLO, N, A, LDA, WORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 )
     $   CALL DLASCL( UPLO, 0, 0, ONE, SIGMA, N, N, A, LDA, INFO )
*
*     Call DSYTRD to reduce symmetric matrix to tridiagonal form.
*
      INDE = 1
      INDTAU = INDE + N
      INDWRK = INDTAU + N
      LLWORK = LWORK - INDWRK + 1
      CALL DSYTRD( UPLO, N, A, LDA, W, WORK( INDE ), WORK( INDTAU ),
     $             WORK( INDWRK ), LLWORK, IINFO )
      LOPT = 2*N + WORK( INDWRK )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, first call
*     DORGTR to generate the orthogonal matrix, then call DSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, WORK( INDE ), INFO )
      ELSE
         CALL DORGTR( UPLO, N, A, LDA, WORK( INDTAU ), WORK( INDWRK ),
     $                LLWORK, IINFO )
         CALL DSTEQR( JOBZ, N, W, WORK( INDE ), A, LDA, WORK( INDTAU ),
     $                INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
*     Set WORK(1) to optimal workspace size.
*
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of DSYEV
*
      END

      SUBROUTINE ZGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHA, BETA,
     $                  VL, LDVL, VR, LDVR, WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*     8-15-00:  Improve consistency of WS calculations (eca)
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGGEV computes for a pair of N-by-N complex nonsymmetric matrices
*  (A,B), the generalized eigenvalues, and optionally, the left and/or
*  right generalized eigenvectors.
*
*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
*  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
*  singular. It is usually represented as the pair (alpha,beta), as
*  there is a reasonable interpretation for beta=0, and even for both
*  being zero.
*
*  The right generalized eigenvector v(j) corresponding to the
*  generalized eigenvalue lambda(j) of (A,B) satisfies
*
*               A * v(j) = lambda(j) * B * v(j).
*
*  The left generalized eigenvector u(j) corresponding to the
*  generalized eigenvalues lambda(j) of (A,B) satisfies
*
*               u(j)**H * A = lambda(j) * u(j)**H * B
*
*  where u(j)**H is the conjugate-transpose of u(j).
*
*  Arguments
*  =========
*
*  JOBVL   (input) CHARACTER*1
*          = 'N':  do not compute the left generalized eigenvectors;
*          = 'V':  compute the left generalized eigenvectors.
*
*  JOBVR   (input) CHARACTER*1
*          = 'N':  do not compute the right generalized eigenvectors;
*          = 'V':  compute the right generalized eigenvectors.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, VL, and VR.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the matrix A in the pair (A,B).
*          On exit, A has been overwritten.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  LDA >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
*          On entry, the matrix B in the pair (A,B).
*          On exit, B has been overwritten.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  LDB >= max(1,N).
*
*  ALPHA   (output) COMPLEX*16 array, dimension (N)
*  BETA    (output) COMPLEX*16 array, dimension (N)
*          On exit, ALPHA(j)/BETA(j), j=1,...,N, will be the
*          generalized eigenvalues.
*
*          Note: the quotients ALPHA(j)/BETA(j) may easily over- or
*          underflow, and BETA(j) may even be zero.  Thus, the user
*          should avoid naively computing the ratio alpha/beta.
*          However, ALPHA will be always less than and usually
*          comparable with norm(A) in magnitude, and BETA always less
*          than and usually comparable with norm(B).
*
*  VL      (output) COMPLEX*16 array, dimension (LDVL,N)
*          If JOBVL = 'V', the left generalized eigenvectors u(j) are
*          stored one after another in the columns of VL, in the same
*          order as their eigenvalues.
*          Each eigenvector will be scaled so the largest component
*          will have abs(real part) + abs(imag. part) = 1.
*          Not referenced if JOBVL = 'N'.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the matrix VL. LDVL >= 1, and
*          if JOBVL = 'V', LDVL >= N.
*
*  VR      (output) COMPLEX*16 array, dimension (LDVR,N)
*          If JOBVR = 'V', the right generalized eigenvectors v(j) are
*          stored one after another in the columns of VR, in the same
*          order as their eigenvalues.
*          Each eigenvector will be scaled so the largest component
*          will have abs(real part) + abs(imag. part) = 1.
*          Not referenced if JOBVR = 'N'.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the matrix VR. LDVR >= 1, and
*          if JOBVR = 'V', LDVR >= N.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,2*N).
*          For good performance, LWORK must generally be larger.
*
*          If LWORK = -1, a workspace query is assumed.  The optimal
*          size for the WORK array is calculated and stored in WORK(1),
*          and no other work except argument checking is performed.
*
*  RWORK   (workspace/output) DOUBLE PRECISION array, dimension (8*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          =1,...,N:
*                The QZ iteration failed.  No eigenvectors have been
*                calculated, but ALPHA(j) and BETA(j) should be
*                correct for j=INFO+1,...,N.
*          > N:  =N+1: other then QZ iteration failed in DHGEQZ,
*                =N+2: error return from DTGEVC.
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            LQUERV
      PARAMETER          ( LQUERV = -1 )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ),
     $                   CONE = ( 1.0D0, 0.0D0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILASCL, ILBSCL, ILV, ILVL, ILVR
      CHARACTER          CHTEMP
      INTEGER            ICOLS, IERR, IHI, IJOBVL, IJOBVR, ILEFT, ILO,
     $                   IN, IRIGHT, IROWS, IRWRK, ITAU, IWRK, JC, JR,
     $                   LWKMIN, LWKOPT
      DOUBLE PRECISION   ANRM, ANRMTO, BIGNUM, BNRM, BNRMTO, EPS,
     $                   SMLNUM, TEMP
      COMPLEX*16         X
*     ..
*     .. Local Arrays ..
      LOGICAL            LDUMMA( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZGEQRF, ZGGBAK, ZGGBAL, ZGGHRD,
     $                   ZHGEQZ, ZLACPY, ZLASCL, ZLASET, ZTGEVC, ZUNGQR,
     $                   ZUNMQR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, MAX, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode the input arguments
*
      IF( LSAME( JOBVL, 'N' ) ) THEN
         IJOBVL = 1
         ILVL = .FALSE.
      ELSE IF( LSAME( JOBVL, 'V' ) ) THEN
         IJOBVL = 2
         ILVL = .TRUE.
      ELSE
         IJOBVL = -1
         ILVL = .FALSE.
      END IF
*
      IF( LSAME( JOBVR, 'N' ) ) THEN
         IJOBVR = 1
         ILVR = .FALSE.
      ELSE IF( LSAME( JOBVR, 'V' ) ) THEN
         IJOBVR = 2
         ILVR = .TRUE.
      ELSE
         IJOBVR = -1
         ILVR = .FALSE.
      END IF
      ILV = ILVL .OR. ILVR
*
*     Test the input arguments
*
      INFO = 0
      IF( IJOBVL.LE.0 ) THEN
         INFO = -1
      ELSE IF( IJOBVR.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDVL.LT.1 .OR. ( ILVL .AND. LDVL.LT.N ) ) THEN
         INFO = -11
      ELSE IF( LDVR.LT.1 .OR. ( ILVR .AND. LDVR.LT.N ) ) THEN
         INFO = -13
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       NB refers to the optimal block size for the immediately
*       following subroutine, as returned by ILAENV. The workspace is
*       computed assuming ILO = 1 and IHI = N, the worst case.)
*
      LWKMIN = 1
      IF( INFO.EQ.0 ) THEN
         LWKOPT = N + N*ILAENV( 1, 'ZGEQRF', ' ', N, 1, N, 0 )
         LWKMIN = MAX( 1, 2*N )
         WORK( 1 ) = LWKOPT
         IF( LWORK.LT.LWKMIN .AND. LWORK.NE.LQUERV )
     $      INFO = -15
      END IF
*
*     Quick returns
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGEV ', -INFO )
         RETURN
      END IF
      IF( LWORK.EQ.LQUERV )
     $   RETURN
      IF( N.EQ.0 )
     $   RETURN
*
*     Get machine constants
*
      EPS = DLAMCH( 'E' )*DLAMCH( 'B' )
      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      CALL DLABAD( SMLNUM, BIGNUM )
      SMLNUM = SQRT( SMLNUM ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', N, N, A, LDA, RWORK )
      ILASCL = .FALSE.
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ANRMTO = SMLNUM
         ILASCL = .TRUE.
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ANRMTO = BIGNUM
         ILASCL = .TRUE.
      END IF
      IF( ILASCL )
     $   CALL ZLASCL( 'G', 0, 0, ANRM, ANRMTO, N, N, A, LDA, IERR )
*
*     Scale B if max element outside range [SMLNUM,BIGNUM]
*
      BNRM = ZLANGE( 'M', N, N, B, LDB, RWORK )
      ILBSCL = .FALSE.
      IF( BNRM.GT.ZERO .AND. BNRM.LT.SMLNUM ) THEN
         BNRMTO = SMLNUM
         ILBSCL = .TRUE.
      ELSE IF( BNRM.GT.BIGNUM ) THEN
         BNRMTO = BIGNUM
         ILBSCL = .TRUE.
      END IF
      IF( ILBSCL )
     $   CALL ZLASCL( 'G', 0, 0, BNRM, BNRMTO, N, N, B, LDB, IERR )
*
*     Permute the matrices A, B to isolate eigenvalues if possible
*     (Real Workspace: need 6*N)
*
      ILEFT = 1
      IRIGHT = N + 1
      IRWRK = IRIGHT + N
      CALL ZGGBAL( 'P', N, A, LDA, B, LDB, ILO, IHI, RWORK( ILEFT ),
     $             RWORK( IRIGHT ), RWORK( IRWRK ), IERR )
*
*     Reduce B to triangular form (QR decomposition of B)
*     (Complex Workspace: need N, prefer N*NB)
*
      IROWS = IHI + 1 - ILO
      IF( ILV ) THEN
         ICOLS = N + 1 - ILO
      ELSE
         ICOLS = IROWS
      END IF
      ITAU = 1
      IWRK = ITAU + IROWS
      CALL ZGEQRF( IROWS, ICOLS, B( ILO, ILO ), LDB, WORK( ITAU ),
     $             WORK( IWRK ), LWORK+1-IWRK, IERR )
*
*     Apply the orthogonal transformation to matrix A
*     (Complex Workspace: need N, prefer N*NB)
*
      CALL ZUNMQR( 'L', 'C', IROWS, ICOLS, IROWS, B( ILO, ILO ), LDB,
     $             WORK( ITAU ), A( ILO, ILO ), LDA, WORK( IWRK ),
     $             LWORK+1-IWRK, IERR )
*
*     Initialize VL
*     (Complex Workspace: need N, prefer N*NB)
*
      IF( ILVL ) THEN
         CALL ZLASET( 'Full', N, N, CZERO, CONE, VL, LDVL )
         CALL ZLACPY( 'L', IROWS-1, IROWS-1, B( ILO+1, ILO ), LDB,
     $                VL( ILO+1, ILO ), LDVL )
         CALL ZUNGQR( IROWS, IROWS, IROWS, VL( ILO, ILO ), LDVL,
     $                WORK( ITAU ), WORK( IWRK ), LWORK+1-IWRK, IERR )
      END IF
*
*     Initialize VR
*
      IF( ILVR )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, VR, LDVR )
*
*     Reduce to generalized Hessenberg form
*
      IF( ILV ) THEN
*
*        Eigenvectors requested -- work on whole matrix.
*
         CALL ZGGHRD( JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB, VL,
     $                LDVL, VR, LDVR, IERR )
      ELSE
         CALL ZGGHRD( 'N', 'N', IROWS, 1, IROWS, A( ILO, ILO ), LDA,
     $                B( ILO, ILO ), LDB, VL, LDVL, VR, LDVR, IERR )
      END IF
*
*     Perform QZ algorithm (Compute eigenvalues, and optionally, the
*     Schur form and Schur vectors)
*     (Complex Workspace: need N)
*     (Real Workspace: need N)
*
      IWRK = ITAU
      IF( ILV ) THEN
         CHTEMP = 'S'
      ELSE
         CHTEMP = 'E'
      END IF
      CALL ZHGEQZ( CHTEMP, JOBVL, JOBVR, N, ILO, IHI, A, LDA, B, LDB,
     $             ALPHA, BETA, VL, LDVL, VR, LDVR, WORK( IWRK ),
     $             LWORK+1-IWRK, RWORK( IRWRK ), IERR )
      IF( IERR.NE.0 ) THEN
         IF( IERR.GT.0 .AND. IERR.LE.N ) THEN
            INFO = IERR
         ELSE IF( IERR.GT.N .AND. IERR.LE.2*N ) THEN
            INFO = IERR - N
         ELSE
            INFO = N + 1
         END IF
         GO TO 70
      END IF
*
*     Compute Eigenvectors
*     (Real Workspace: need 2*N)
*     (Complex Workspace: need 2*N)
*
      IF( ILV ) THEN
         IF( ILVL ) THEN
            IF( ILVR ) THEN
               CHTEMP = 'B'
            ELSE
               CHTEMP = 'L'
            END IF
         ELSE
            CHTEMP = 'R'
         END IF
*
         CALL ZTGEVC( CHTEMP, 'B', LDUMMA, N, A, LDA, B, LDB, VL, LDVL,
     $                VR, LDVR, N, IN, WORK( IWRK ), RWORK( IRWRK ),
     $                IERR )
         IF( IERR.NE.0 ) THEN
            INFO = N + 2
            GO TO 70
         END IF
*
*        Undo balancing on VL and VR and normalization
*        (Workspace: none needed)
*
         IF( ILVL ) THEN
            CALL ZGGBAK( 'P', 'L', N, ILO, IHI, RWORK( ILEFT ),
     $                   RWORK( IRIGHT ), N, VL, LDVL, IERR )
            DO 30 JC = 1, N
               TEMP = ZERO
               DO 10 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VL( JR, JC ) ) )
   10          CONTINUE
               IF( TEMP.LT.SMLNUM )
     $            GO TO 30
               TEMP = ONE / TEMP
               DO 20 JR = 1, N
                  VL( JR, JC ) = VL( JR, JC )*TEMP
   20          CONTINUE
   30       CONTINUE
         END IF
         IF( ILVR ) THEN
            CALL ZGGBAK( 'P', 'R', N, ILO, IHI, RWORK( ILEFT ),
     $                   RWORK( IRIGHT ), N, VR, LDVR, IERR )
            DO 60 JC = 1, N
               TEMP = ZERO
               DO 40 JR = 1, N
                  TEMP = MAX( TEMP, ABS1( VR( JR, JC ) ) )
   40          CONTINUE
               IF( TEMP.LT.SMLNUM )
     $            GO TO 60
               TEMP = ONE / TEMP
               DO 50 JR = 1, N
                  VR( JR, JC ) = VR( JR, JC )*TEMP
   50          CONTINUE
   60       CONTINUE
         END IF
      END IF
*
*     Undo scaling if necessary
*
      IF( ILASCL )
     $   CALL ZLASCL( 'G', 0, 0, ANRMTO, ANRM, N, 1, ALPHA, N, IERR )
*
      IF( ILBSCL )
     $   CALL ZLASCL( 'G', 0, 0, BNRMTO, BNRM, N, 1, BETA, N, IERR )
*
   70 CONTINUE
      WORK( 1 ) = LWKOPT
*
      RETURN
*
*     End of ZGGEV
*
      END

      DOUBLE PRECISION FUNCTION ZLANGE( NORM, M, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLANGE  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  complex matrix A.
*
*  Description
*  ===========
*
*  ZLANGE returns the value
*
*     ZLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in ZLANGE as described
*          above.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.  When M = 0,
*          ZLANGE is set to zero.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.  When N = 0,
*          ZLANGE is set to zero.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( MIN( M, N ).EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, M
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, M
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 I = 1, M
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, M
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, M
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL ZLASSQ( M, A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      ZLANGE = VALUE
      RETURN
*
*     End of ZLANGE
*
      END

      SUBROUTINE ZGGBAL( JOB, N, A, LDA, B, LDB, ILO, IHI, LSCALE,
     $                   RSCALE, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     April 12, 2001
*
*     .. Scalar Arguments ..
      CHARACTER          JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   LSCALE( * ), RSCALE( * ), WORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGGBAL balances a pair of general complex matrices (A,B).  This
*  involves, first, permuting A and B by similarity transformations to
*  isolate eigenvalues in the first 1 to ILO$-$1 and last IHI+1 to N
*  elements on the diagonal; and second, applying a diagonal similarity
*  transformation to rows and columns ILO to IHI to make the rows
*  and columns as close in norm as possible. Both steps are optional.
*
*  Balancing may reduce the 1-norm of the matrices, and improve the
*  accuracy of the computed eigenvalues and/or eigenvectors in the
*  generalized eigenvalue problem A*x = lambda*B*x.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies the operations to be performed on A and B:
*          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0
*                  and RSCALE(I) = 1.0 for i=1,...,N;
*          = 'P':  permute only;
*          = 'S':  scale only;
*          = 'B':  both permute and scale.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the input matrix A.
*          On exit, A is overwritten by the balanced matrix.
*          If JOB = 'N', A is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB,N)
*          On entry, the input matrix B.
*          On exit, B is overwritten by the balanced matrix.
*          If JOB = 'N', B is not referenced.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  ILO     (output) INTEGER
*  IHI     (output) INTEGER
*          ILO and IHI are set to integers such that on exit
*          A(i,j) = 0 and B(i,j) = 0 if i > j and
*          j = 1,...,ILO-1 or i = IHI+1,...,N.
*          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
*
*  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and scaling factors applied
*          to the left side of A and B.  If P(j) is the index of the
*          row interchanged with row j, and D(j) is the scaling factor
*          applied to row j, then
*            LSCALE(j) = P(j)    for J = 1,...,ILO-1
*                      = D(j)    for J = ILO,...,IHI
*                      = P(j)    for J = IHI+1,...,N.
*          The order in which the interchanges are made is N to IHI+1,
*          then 1 to ILO-1.
*
*  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and scaling factors applied
*          to the right side of A and B.  If P(j) is the index of the
*          column interchanged with column j, and D(j) is the scaling
*          factor applied to column j, then
*            RSCALE(j) = P(j)    for J = 1,...,ILO-1
*                      = D(j)    for J = ILO,...,IHI
*                      = P(j)    for J = IHI+1,...,N.
*          The order in which the interchanges are made is N to IHI+1,
*          then 1 to ILO-1.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  See R.C. WARD, Balancing the generalized eigenvalue problem,
*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, HALF, ONE
      PARAMETER          ( ZERO = 0.0D+0, HALF = 0.5D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   THREE, SCLFAC
      PARAMETER          ( THREE = 3.0D+0, SCLFAC = 1.0D+1 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, ICAB, IFLOW, IP1, IR, IRAB, IT, J, JC, JP1,
     $                   K, KOUNT, L, LCAB, LM1, LRAB, LSFMAX, LSFMIN,
     $                   M, NR, NRP2
      DOUBLE PRECISION   ALPHA, BASL, BETA, CAB, CMAX, COEF, COEF2,
     $                   COEF5, COR, EW, EWC, GAMMA, PGAMMA, RAB, SFMAX,
     $                   SFMIN, SUM, T, TA, TB, TC
      COMPLEX*16         CDUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            IZAMAX
      DOUBLE PRECISION   DDOT, DLAMCH
      EXTERNAL           LSAME, IZAMAX, DDOT, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DAXPY, DSCAL, XERBLA, ZDSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DIMAG, INT, LOG10, MAX, MIN, SIGN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CDUM ) = ABS( DBLE( CDUM ) ) + ABS( DIMAG( CDUM ) )
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGBAL', -INFO )
         RETURN
      END IF
*
      K = 1
      L = N
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( LSAME( JOB, 'N' ) ) THEN
         ILO = 1
         IHI = N
         DO 10 I = 1, N
            LSCALE( I ) = ONE
            RSCALE( I ) = ONE
   10    CONTINUE
         RETURN
      END IF
*
      IF( K.EQ.L ) THEN
         ILO = 1
         IHI = 1
         LSCALE( 1 ) = ONE
         RSCALE( 1 ) = ONE
         RETURN
      END IF
*
      IF( LSAME( JOB, 'S' ) )
     $   GO TO 190
*
      GO TO 30
*
*     Permute the matrices A and B to isolate the eigenvalues.
*
*     Find row with one nonzero in columns 1 through L
*
   20 CONTINUE
      L = LM1
      IF( L.NE.1 )
     $   GO TO 30
*
      RSCALE( 1 ) = ONE
      LSCALE( 1 ) = ONE
      GO TO 190
*
   30 CONTINUE
      LM1 = L - 1
      DO 80 I = L, 1, -1
         DO 40 J = 1, LM1
            JP1 = J + 1
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO )
     $         GO TO 50
   40    CONTINUE
         J = L
         GO TO 70
*
   50    CONTINUE
         DO 60 J = JP1, L
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO )
     $         GO TO 80
   60    CONTINUE
         J = JP1 - 1
*
   70    CONTINUE
         M = L
         IFLOW = 1
         GO TO 160
   80 CONTINUE
      GO TO 100
*
*     Find column with one nonzero in rows K through N
*
   90 CONTINUE
      K = K + 1
*
  100 CONTINUE
      DO 150 J = K, L
         DO 110 I = K, LM1
            IP1 = I + 1
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO )
     $         GO TO 120
  110    CONTINUE
         I = L
         GO TO 140
  120    CONTINUE
         DO 130 I = IP1, L
            IF( A( I, J ).NE.CZERO .OR. B( I, J ).NE.CZERO )
     $         GO TO 150
  130    CONTINUE
         I = IP1 - 1
  140    CONTINUE
         M = K
         IFLOW = 2
         GO TO 160
  150 CONTINUE
      GO TO 190
*
*     Permute rows M and I
*
  160 CONTINUE
      LSCALE( M ) = DBLE( I )
      IF( I.EQ.M )
     $   GO TO 170
      CALL ZSWAP( N-K+1, A( I, K ), LDA, A( M, K ), LDA )
      CALL ZSWAP( N-K+1, B( I, K ), LDB, B( M, K ), LDB )
*
*     Permute columns M and J
*
  170 CONTINUE
      RSCALE( M ) = DBLE( J )
      IF( J.EQ.M )
     $   GO TO 180
      CALL ZSWAP( L, A( 1, J ), 1, A( 1, M ), 1 )
      CALL ZSWAP( L, B( 1, J ), 1, B( 1, M ), 1 )
*
  180 CONTINUE
      GO TO ( 20, 90 )IFLOW
*
  190 CONTINUE
      ILO = K
      IHI = L
*
      IF( ILO.EQ.IHI )
     $   RETURN
*
      IF( LSAME( JOB, 'P' ) )
     $   RETURN
*
*     Balance the submatrix in rows ILO to IHI.
*
      NR = IHI - ILO + 1
      DO 200 I = ILO, IHI
         RSCALE( I ) = ZERO
         LSCALE( I ) = ZERO
*
         WORK( I ) = ZERO
         WORK( I+N ) = ZERO
         WORK( I+2*N ) = ZERO
         WORK( I+3*N ) = ZERO
         WORK( I+4*N ) = ZERO
         WORK( I+5*N ) = ZERO
  200 CONTINUE
*
*     Compute right side vector in resulting linear equations
*
      BASL = LOG10( SCLFAC )
      DO 240 I = ILO, IHI
         DO 230 J = ILO, IHI
            IF( A( I, J ).EQ.CZERO ) THEN
               TA = ZERO
               GO TO 210
            END IF
            TA = LOG10( CABS1( A( I, J ) ) ) / BASL
*
  210       CONTINUE
            IF( B( I, J ).EQ.CZERO ) THEN
               TB = ZERO
               GO TO 220
            END IF
            TB = LOG10( CABS1( B( I, J ) ) ) / BASL
*
  220       CONTINUE
            WORK( I+4*N ) = WORK( I+4*N ) - TA - TB
            WORK( J+5*N ) = WORK( J+5*N ) - TA - TB
  230    CONTINUE
  240 CONTINUE
*
      COEF = ONE / DBLE( 2*NR )
      COEF2 = COEF*COEF
      COEF5 = HALF*COEF2
      NRP2 = NR + 2
      BETA = ZERO
      IT = 1
*
*     Start generalized conjugate gradient iteration
*
  250 CONTINUE
*
      GAMMA = DDOT( NR, WORK( ILO+4*N ), 1, WORK( ILO+4*N ), 1 ) +
     $        DDOT( NR, WORK( ILO+5*N ), 1, WORK( ILO+5*N ), 1 )
*
      EW = ZERO
      EWC = ZERO
      DO 260 I = ILO, IHI
         EW = EW + WORK( I+4*N )
         EWC = EWC + WORK( I+5*N )
  260 CONTINUE
*
      GAMMA = COEF*GAMMA - COEF2*( EW**2+EWC**2 ) - COEF5*( EW-EWC )**2
      IF( GAMMA.EQ.ZERO )
     $   GO TO 350
      IF( IT.NE.1 )
     $   BETA = GAMMA / PGAMMA
      T = COEF5*( EWC-THREE*EW )
      TC = COEF5*( EW-THREE*EWC )
*
      CALL DSCAL( NR, BETA, WORK( ILO ), 1 )
      CALL DSCAL( NR, BETA, WORK( ILO+N ), 1 )
*
      CALL DAXPY( NR, COEF, WORK( ILO+4*N ), 1, WORK( ILO+N ), 1 )
      CALL DAXPY( NR, COEF, WORK( ILO+5*N ), 1, WORK( ILO ), 1 )
*
      DO 270 I = ILO, IHI
         WORK( I ) = WORK( I ) + TC
         WORK( I+N ) = WORK( I+N ) + T
  270 CONTINUE
*
*     Apply matrix to vector
*
      DO 300 I = ILO, IHI
         KOUNT = 0
         SUM = ZERO
         DO 290 J = ILO, IHI
            IF( A( I, J ).EQ.CZERO )
     $         GO TO 280
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  280       CONTINUE
            IF( B( I, J ).EQ.CZERO )
     $         GO TO 290
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( J )
  290    CONTINUE
         WORK( I+2*N ) = DBLE( KOUNT )*WORK( I+N ) + SUM
  300 CONTINUE
*
      DO 330 J = ILO, IHI
         KOUNT = 0
         SUM = ZERO
         DO 320 I = ILO, IHI
            IF( A( I, J ).EQ.CZERO )
     $         GO TO 310
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  310       CONTINUE
            IF( B( I, J ).EQ.CZERO )
     $         GO TO 320
            KOUNT = KOUNT + 1
            SUM = SUM + WORK( I+N )
  320    CONTINUE
         WORK( J+3*N ) = DBLE( KOUNT )*WORK( J ) + SUM
  330 CONTINUE
*
      SUM = DDOT( NR, WORK( ILO+N ), 1, WORK( ILO+2*N ), 1 ) +
     $      DDOT( NR, WORK( ILO ), 1, WORK( ILO+3*N ), 1 )
      ALPHA = GAMMA / SUM
*
*     Determine correction to current iteration
*
      CMAX = ZERO
      DO 340 I = ILO, IHI
         COR = ALPHA*WORK( I+N )
         IF( ABS( COR ).GT.CMAX )
     $      CMAX = ABS( COR )
         LSCALE( I ) = LSCALE( I ) + COR
         COR = ALPHA*WORK( I )
         IF( ABS( COR ).GT.CMAX )
     $      CMAX = ABS( COR )
         RSCALE( I ) = RSCALE( I ) + COR
  340 CONTINUE
      IF( CMAX.LT.HALF )
     $   GO TO 350
*
      CALL DAXPY( NR, -ALPHA, WORK( ILO+2*N ), 1, WORK( ILO+4*N ), 1 )
      CALL DAXPY( NR, -ALPHA, WORK( ILO+3*N ), 1, WORK( ILO+5*N ), 1 )
*
      PGAMMA = GAMMA
      IT = IT + 1
      IF( IT.LE.NRP2 )
     $   GO TO 250
*
*     End generalized conjugate gradient iteration
*
  350 CONTINUE
      SFMIN = DLAMCH( 'S' )
      SFMAX = ONE / SFMIN
      LSFMIN = INT( LOG10( SFMIN ) / BASL+ONE )
      LSFMAX = INT( LOG10( SFMAX ) / BASL )
      DO 360 I = ILO, IHI
         IRAB = IZAMAX( N-ILO+1, A( I, ILO ), LDA )
         RAB = ABS( A( I, IRAB+ILO-1 ) )
         IRAB = IZAMAX( N-ILO+1, B( I, ILO ), LDB )
         RAB = MAX( RAB, ABS( B( I, IRAB+ILO-1 ) ) )
         LRAB = INT( LOG10( RAB+SFMIN ) / BASL+ONE )
         IR = LSCALE( I ) + SIGN( HALF, LSCALE( I ) )
         IR = MIN( MAX( IR, LSFMIN ), LSFMAX, LSFMAX-LRAB )
         LSCALE( I ) = SCLFAC**IR
         ICAB = IZAMAX( IHI, A( 1, I ), 1 )
         CAB = ABS( A( ICAB, I ) )
         ICAB = IZAMAX( IHI, B( 1, I ), 1 )
         CAB = MAX( CAB, ABS( B( ICAB, I ) ) )
         LCAB = INT( LOG10( CAB+SFMIN ) / BASL+ONE )
         JC = RSCALE( I ) + SIGN( HALF, RSCALE( I ) )
         JC = MIN( MAX( JC, LSFMIN ), LSFMAX, LSFMAX-LCAB )
         RSCALE( I ) = SCLFAC**JC
  360 CONTINUE
*
*     Row scaling of matrices A and B
*
      DO 370 I = ILO, IHI
         CALL ZDSCAL( N-ILO+1, LSCALE( I ), A( I, ILO ), LDA )
         CALL ZDSCAL( N-ILO+1, LSCALE( I ), B( I, ILO ), LDB )
  370 CONTINUE
*
*     Column scaling of matrices A and B
*
      DO 380 J = ILO, IHI
         CALL ZDSCAL( IHI, RSCALE( J ), A( 1, J ), 1 )
         CALL ZDSCAL( IHI, RSCALE( J ), B( 1, J ), 1 )
  380 CONTINUE
*
      RETURN
*
*     End of ZGGBAL
*
      END


      SUBROUTINE ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNMQR overwrites the general complex M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'C':      Q**H * C       C * Q**H
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by ZGEQRF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**H from the Left;
*          = 'R': apply Q or Q**H from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'C':  Conjugate transpose, apply Q**H.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          ZGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQRF.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK,
     $                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. Local Arrays ..
      COMPLEX*16         T( LDT, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNM2R
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.  NB may be at most NBMAX, where NBMAX
*        is used to define the local array T.
*
         NB = MIN( NBMAX, ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M, N, K,
     $        -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMQR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZUNMQR', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        Use blocked code
*
         IF( ( LEFT .AND. .NOT.NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL ZLARFT( 'Forward', 'Columnwise', NQ-I+1, IB, A( I, I ),
     $                   LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H or H' is applied to C(i:m,1:n)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H or H' is applied to C(1:m,i:n)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           Apply H or H'
*
            CALL ZLARFB( SIDE, TRANS, 'Forward', 'Columnwise', MI, NI,
     $                   IB, A( I, I ), LDA, T, LDT, C( IC, JC ), LDC,
     $                   WORK, LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of ZUNMQR
*
      END

      SUBROUTINE ZGGHRD( COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB, Q,
     $                   LDQ, Z, LDZ, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), Q( LDQ, * ),
     $                   Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGGHRD reduces a pair of complex matrices (A,B) to generalized upper
*  Hessenberg form using unitary transformations, where A is a
*  general matrix and B is upper triangular:  Q' * A * Z = H and
*  Q' * B * Z = T, where H is upper Hessenberg, T is upper triangular,
*  and Q and Z are unitary, and ' means conjugate transpose.
*
*  The unitary matrices Q and Z are determined as products of Givens
*  rotations.  They may either be formed explicitly, or they may be
*  postmultiplied into input matrices Q1 and Z1, so that
*
*       Q1 * A * Z1' = (Q1*Q) * H * (Z1*Z)'
*       Q1 * B * Z1' = (Q1*Q) * T * (Z1*Z)'
*
*  Arguments
*  =========
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': do not compute Q;
*          = 'I': Q is initialized to the unit matrix, and the
*                 unitary matrix Q is returned;
*          = 'V': Q must contain a unitary matrix Q1 on entry,
*                 and the product Q1*Q is returned.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': do not compute Q;
*          = 'I': Q is initialized to the unit matrix, and the
*                 unitary matrix Q is returned;
*          = 'V': Q must contain a unitary matrix Q1 on entry,
*                 and the product Q1*Q is returned.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that A is already upper triangular in rows and
*          columns 1:ILO-1 and IHI+1:N.  ILO and IHI are normally set
*          by a previous call to ZGGBAL; otherwise they should be set
*          to 1 and N respectively.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the N-by-N general matrix to be reduced.
*          On exit, the upper triangle and the first subdiagonal of A
*          are overwritten with the upper Hessenberg matrix H, and the
*          rest is set to zero.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
*          On entry, the N-by-N upper triangular matrix B.
*          On exit, the upper triangular matrix T = Q' B Z.  The
*          elements below the diagonal are set to zero.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  Q       (input/output) COMPLEX*16 array, dimension (LDQ, N)
*          If COMPQ='N':  Q is not referenced.
*          If COMPQ='I':  on entry, Q need not be set, and on exit it
*                         contains the unitary matrix Q, where Q'
*                         is the product of the Givens transformations
*                         which are applied to A and B on the left.
*          If COMPQ='V':  on entry, Q must contain a unitary matrix
*                         Q1, and on exit this is overwritten by Q1*Q.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.
*
*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
*          If COMPZ='N':  Z is not referenced.
*          If COMPZ='I':  on entry, Z need not be set, and on exit it
*                         contains the unitary matrix Z, which is
*                         the product of the Givens transformations
*                         which are applied to A and B on the right.
*          If COMPZ='V':  on entry, Z must contain a unitary matrix
*                         Z1, and on exit this is overwritten by Z1*Z.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.
*          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  This routine reduces A to Hessenberg and B to triangular form by
*  an unblocked reduction, as described in _Matrix_Computations_,
*  by Golub and van Loan (Johns Hopkins Press).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CONE, CZERO
      PARAMETER          ( CONE = ( 1.0D+0, 0.0D+0 ),
     $                   CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILQ, ILZ
      INTEGER            ICOMPQ, ICOMPZ, JCOL, JROW
      DOUBLE PRECISION   C
      COMPLEX*16         CTEMP, S
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARTG, ZLASET, ZROT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Decode COMPQ
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
*     Decode COMPZ
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Test the input parameters.
*
      INFO = 0
      IF( ICOMPQ.LE.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPZ.LE.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -7
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -9
      ELSE IF( ( ILQ .AND. LDQ.LT.N ) .OR. LDQ.LT.1 ) THEN
         INFO = -11
      ELSE IF( ( ILZ .AND. LDZ.LT.N ) .OR. LDZ.LT.1 ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGHRD', -INFO )
         RETURN
      END IF
*
*     Initialize Q and Z if desired.
*
      IF( ICOMPQ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
*     Quick return if possible
*
      IF( N.LE.1 )
     $   RETURN
*
*     Zero out lower triangle of B
*
      DO 20 JCOL = 1, N - 1
         DO 10 JROW = JCOL + 1, N
            B( JROW, JCOL ) = CZERO
   10    CONTINUE
   20 CONTINUE
*
*     Reduce A and B
*
      DO 40 JCOL = ILO, IHI - 2
*
         DO 30 JROW = IHI, JCOL + 2, -1
*
*           Step 1: rotate rows JROW-1, JROW to kill A(JROW,JCOL)
*
            CTEMP = A( JROW-1, JCOL )
            CALL ZLARTG( CTEMP, A( JROW, JCOL ), C, S,
     $                   A( JROW-1, JCOL ) )
            A( JROW, JCOL ) = CZERO
            CALL ZROT( N-JCOL, A( JROW-1, JCOL+1 ), LDA,
     $                 A( JROW, JCOL+1 ), LDA, C, S )
            CALL ZROT( N+2-JROW, B( JROW-1, JROW-1 ), LDB,
     $                 B( JROW, JROW-1 ), LDB, C, S )
            IF( ILQ )
     $         CALL ZROT( N, Q( 1, JROW-1 ), 1, Q( 1, JROW ), 1, C,
     $                    DCONJG( S ) )
*
*           Step 2: rotate columns JROW, JROW-1 to kill B(JROW,JROW-1)
*
            CTEMP = B( JROW, JROW )
            CALL ZLARTG( CTEMP, B( JROW, JROW-1 ), C, S,
     $                   B( JROW, JROW ) )
            B( JROW, JROW-1 ) = CZERO
            CALL ZROT( IHI, A( 1, JROW ), 1, A( 1, JROW-1 ), 1, C, S )
            CALL ZROT( JROW-1, B( 1, JROW ), 1, B( 1, JROW-1 ), 1, C,
     $                 S )
            IF( ILZ )
     $         CALL ZROT( N, Z( 1, JROW ), 1, Z( 1, JROW-1 ), 1, C, S )
   30    CONTINUE
   40 CONTINUE
*
      RETURN
*
*     End of ZGGHRD
*
      END

      SUBROUTINE ZROT( N, CX, INCX, CY, INCY, C, S )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   C
      COMPLEX*16         S
*     ..
*     .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
*     ..
*
*  Purpose
*  =======
*
*  ZROT   applies a plane rotation, where the cos (C) is real and the
*  sin (S) is complex, and the vectors CX and CY are complex.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements in the vectors CX and CY.
*
*  CX      (input/output) COMPLEX*16 array, dimension (N)
*          On input, the vector X.
*          On output, CX is overwritten with C*X + S*Y.
*
*  INCX    (input) INTEGER
*          The increment between successive values of CY.  INCX <> 0.
*
*  CY      (input/output) COMPLEX*16 array, dimension (N)
*          On input, the vector Y.
*          On output, CY is overwritten with -CONJG(S)*X + C*Y.
*
*  INCY    (input) INTEGER
*          The increment between successive values of CY.  INCX <> 0.
*
*  C       (input) DOUBLE PRECISION
*  S       (input) COMPLEX*16
*          C and S define a rotation
*             [  C          S  ]
*             [ -conjg(S)   C  ]
*          where C*C + S*CONJG(S) = 1.0.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX*16         STEMP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 )
     $   GO TO 20
*
*     Code for unequal increments or equal increments not equal to 1
*
      IX = 1
      IY = 1
      IF( INCX.LT.0 )
     $   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 )
     $   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         STEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - DCONJG( S )*CX( IX )
         CX( IX ) = STEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
*
*     Code for both increments equal to 1
*
   20 CONTINUE
      DO 30 I = 1, N
         STEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - DCONJG( S )*CX( I )
         CX( I ) = STEMP
   30 CONTINUE
      RETURN
      END

      SUBROUTINE ZUNM2R( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNM2R overwrites the general complex m-by-n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'C', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'C',
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by ZGEQRF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'C': apply Q' (Conjugate transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          ZGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGEQRF.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the m-by-n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX*16 array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX*16         AII, TAUI
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, NQ ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNM2R', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. .NOT.NOTRAN .OR. .NOT.LEFT .AND. NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) or H(i)' is applied to C(i:m,1:n)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(i) or H(i)' is applied to C(1:m,i:n)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        Apply H(i) or H(i)'
*
         IF( NOTRAN ) THEN
            TAUI = TAU( I )
         ELSE
            TAUI = DCONJG( TAU( I ) )
         END IF
         AII = A( I, I )
         A( I, I ) = ONE
         CALL ZLARF( SIDE, MI, NI, A( I, I ), 1, TAUI, C( IC, JC ), LDC,
     $               WORK )
         A( I, I ) = AII
   10 CONTINUE
      RETURN
*
*     End of ZUNM2R
*
      END

      SUBROUTINE ZGEQRF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEQRF computes a QR factorization of a complex M-by-N matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of min(m,n) elementary reflectors (see Further
*          Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEQR2, ZLARFB, ZLARFT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQRF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = N
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'ZGEQRF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = N
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZGEQRF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the QR factorization of the current block
*           A(i:m,i:i+ib-1)
*
            CALL ZGEQR2( M-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.N ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL ZLARFT( 'Forward', 'Columnwise', M-I+1, IB,
     $                      A( I, I ), LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H' to A(i:m,i+ib:n) from the left
*
               CALL ZLARFB( 'Left', 'Conjugate transpose', 'Forward',
     $                      'Columnwise', M-I+1, N-I-IB+1, IB,
     $                      A( I, I ), LDA, WORK, LDWORK, A( I, I+IB ),
     $                      LDA, WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL ZGEQR2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZGEQRF
*
      END

      SUBROUTINE ZGEQR2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEQR2 computes a QR factorization of a complex m by n matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(m,n) by n upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARF, ZLARFG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGEQR2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
         CALL ZLARFG( M-I+1, A( I, I ), A( MIN( I+1, M ), I ), 1,
     $                TAU( I ) )
         IF( I.LT.N ) THEN
*
*           Apply H(i)' to A(i:m,i+1:n) from the left
*
            ALPHA = A( I, I )
            A( I, I ) = ONE
            CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1,
     $                  DCONJG( TAU( I ) ), A( I, I+1 ), LDA, WORK )
            A( I, I ) = ALPHA
         END IF
   10 CONTINUE
      RETURN
*
*     End of ZGEQR2
*
      END

      SUBROUTINE ZHGEQZ( JOB, COMPQ, COMPZ, N, ILO, IHI, A, LDA, B, LDB,
     $                   ALPHA, BETA, Q, LDQ, Z, LDZ, WORK, LWORK,
     $                   RWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          COMPQ, COMPZ, JOB
      INTEGER            IHI, ILO, INFO, LDA, LDB, LDQ, LDZ, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), ALPHA( * ), B( LDB, * ),
     $                   BETA( * ), Q( LDQ, * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZHGEQZ implements a single-shift version of the QZ
*  method for finding the generalized eigenvalues w(i)=ALPHA(i)/BETA(i)
*  of the equation
*
*       det( A - w(i) B ) = 0
*
*  If JOB='S', then the pair (A,B) is simultaneously
*  reduced to Schur form (i.e., A and B are both upper triangular) by
*  applying one unitary tranformation (usually called Q) on the left and
*  another (usually called Z) on the right.  The diagonal elements of
*  A are then ALPHA(1),...,ALPHA(N), and of B are BETA(1),...,BETA(N).
*
*  If JOB='S' and COMPQ and COMPZ are 'V' or 'I', then the unitary
*  transformations used to reduce (A,B) are accumulated into the arrays
*  Q and Z s.t.:
*
*       Q(in) A(in) Z(in)* = Q(out) A(out) Z(out)*
*       Q(in) B(in) Z(in)* = Q(out) B(out) Z(out)*
*
*  Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*       Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*       pp. 241--256.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          = 'E': compute only ALPHA and BETA.  A and B will not
*                 necessarily be put into generalized Schur form.
*          = 'S': put A and B into generalized Schur form, as well
*                 as computing ALPHA and BETA.
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': do not modify Q.
*          = 'V': multiply the array Q on the right by the conjugate
*                 transpose of the unitary tranformation that is
*                 applied to the left side of A and B to reduce them
*                 to Schur form.
*          = 'I': like COMPQ='V', except that Q will be initialized to
*                 the identity first.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': do not modify Z.
*          = 'V': multiply the array Z on the right by the unitary
*                 tranformation that is applied to the right side of
*                 A and B to reduce them to Schur form.
*          = 'I': like COMPZ='V', except that Z will be initialized to
*                 the identity first.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, Q, and Z.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          It is assumed that A is already upper triangular in rows and
*          columns 1:ILO-1 and IHI+1:N.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA, N)
*          On entry, the N-by-N upper Hessenberg matrix A.  Elements
*          below the subdiagonal must be zero.
*          If JOB='S', then on exit A and B will have been
*             simultaneously reduced to upper triangular form.
*          If JOB='E', then on exit A will have been destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max( 1, N ).
*
*  B       (input/output) COMPLEX*16 array, dimension (LDB, N)
*          On entry, the N-by-N upper triangular matrix B.  Elements
*          below the diagonal must be zero.
*          If JOB='S', then on exit A and B will have been
*             simultaneously reduced to upper triangular form.
*          If JOB='E', then on exit B will have been destroyed.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max( 1, N ).
*
*  ALPHA   (output) COMPLEX*16 array, dimension (N)
*          The diagonal elements of A when the pair (A,B) has been
*          reduced to Schur form.  ALPHA(i)/BETA(i) i=1,...,N
*          are the generalized eigenvalues.
*
*  BETA    (output) COMPLEX*16 array, dimension (N)
*          The diagonal elements of B when the pair (A,B) has been
*          reduced to Schur form.  ALPHA(i)/BETA(i) i=1,...,N
*          are the generalized eigenvalues.  A and B are normalized
*          so that BETA(1),...,BETA(N) are non-negative real numbers.
*
*  Q       (input/output) COMPLEX*16 array, dimension (LDQ, N)
*          If COMPQ='N', then Q will not be referenced.
*          If COMPQ='V' or 'I', then the conjugate transpose of the
*             unitary transformations which are applied to A and B on
*             the left will be applied to the array Q on the right.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= 1.
*          If COMPQ='V' or 'I', then LDQ >= N.
*
*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
*          If COMPZ='N', then Z will not be referenced.
*          If COMPZ='V' or 'I', then the unitary transformations which
*             are applied to A and B on the right will be applied to the
*             array Z on the right.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If COMPZ='V' or 'I', then LDZ >= N.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1,...,N: the QZ iteration did not converge.  (A,B) is not
*                     in Schur form, but ALPHA(i) and BETA(i),
*                     i=INFO+1,...,N should be correct.
*          = N+1,...,2*N: the shift calculation failed.  (A,B) is not
*                     in Schur form, but ALPHA(i) and BETA(i),
*                     i=INFO-N+1,...,N should be correct.
*          > 2*N:     various "impossible" errors.
*
*  Further Details
*  ===============
*
*  We assume that complex ABS works as long as its value is less than
*  overflow.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      DOUBLE PRECISION   HALF
      PARAMETER          ( HALF = 0.5D+0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ILAZR2, ILAZRO, ILQ, ILSCHR, ILZ, LQUERY
      INTEGER            ICOMPQ, ICOMPZ, IFIRST, IFRSTM, IITER, ILAST,
     $                   ILASTM, IN, ISCHUR, ISTART, J, JC, JCH, JITER,
     $                   JR, MAXIT
      DOUBLE PRECISION   ABSB, ANORM, ASCALE, ATOL, BNORM, BSCALE, BTOL,
     $                   C, SAFMIN, TEMP, TEMP2, TEMPR, ULP
      COMPLEX*16         ABI22, AD11, AD12, AD21, AD22, CTEMP, CTEMP2,
     $                   CTEMP3, ESHIFT, RTDISC, S, SHIFT, SIGNBC, T,
     $                   U12, X
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHS
      EXTERNAL           LSAME, DLAMCH, ZLANHS
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARTG, ZLASET, ZROT, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN,
     $                   SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode JOB, COMPQ, COMPZ
*
      IF( LSAME( JOB, 'E' ) ) THEN
         ILSCHR = .FALSE.
         ISCHUR = 1
      ELSE IF( LSAME( JOB, 'S' ) ) THEN
         ILSCHR = .TRUE.
         ISCHUR = 2
      ELSE
         ISCHUR = 0
      END IF
*
      IF( LSAME( COMPQ, 'N' ) ) THEN
         ILQ = .FALSE.
         ICOMPQ = 1
      ELSE IF( LSAME( COMPQ, 'V' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 2
      ELSE IF( LSAME( COMPQ, 'I' ) ) THEN
         ILQ = .TRUE.
         ICOMPQ = 3
      ELSE
         ICOMPQ = 0
      END IF
*
      IF( LSAME( COMPZ, 'N' ) ) THEN
         ILZ = .FALSE.
         ICOMPZ = 1
      ELSE IF( LSAME( COMPZ, 'V' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 2
      ELSE IF( LSAME( COMPZ, 'I' ) ) THEN
         ILZ = .TRUE.
         ICOMPZ = 3
      ELSE
         ICOMPZ = 0
      END IF
*
*     Check Argument Values
*
      INFO = 0
      WORK( 1 ) = MAX( 1, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( ISCHUR.EQ.0 ) THEN
         INFO = -1
      ELSE IF( ICOMPQ.EQ.0 ) THEN
         INFO = -2
      ELSE IF( ICOMPZ.EQ.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -5
      ELSE IF( IHI.GT.N .OR. IHI.LT.ILO-1 ) THEN
         INFO = -6
      ELSE IF( LDA.LT.N ) THEN
         INFO = -8
      ELSE IF( LDB.LT.N ) THEN
         INFO = -10
      ELSE IF( LDQ.LT.1 .OR. ( ILQ .AND. LDQ.LT.N ) ) THEN
         INFO = -14
      ELSE IF( LDZ.LT.1 .OR. ( ILZ .AND. LDZ.LT.N ) ) THEN
         INFO = -16
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -18
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHGEQZ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
*     WORK( 1 ) = CMPLX( 1 )
      IF( N.LE.0 ) THEN
         WORK( 1 ) = DCMPLX( 1 )
         RETURN
      END IF
*
*     Initialize Q and Z
*
      IF( ICOMPQ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
      IF( ICOMPZ.EQ.3 )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Z, LDZ )
*
*     Machine Constants
*
      IN = IHI + 1 - ILO
      SAFMIN = DLAMCH( 'S' )
      ULP = DLAMCH( 'E' )*DLAMCH( 'B' )
      ANORM = ZLANHS( 'F', IN, A( ILO, ILO ), LDA, RWORK )
      BNORM = ZLANHS( 'F', IN, B( ILO, ILO ), LDB, RWORK )
      ATOL = MAX( SAFMIN, ULP*ANORM )
      BTOL = MAX( SAFMIN, ULP*BNORM )
      ASCALE = ONE / MAX( SAFMIN, ANORM )
      BSCALE = ONE / MAX( SAFMIN, BNORM )
*
*
*     Set Eigenvalues IHI+1:N
*
      DO 10 J = IHI + 1, N
         ABSB = ABS( B( J, J ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( B( J, J ) / ABSB )
            B( J, J ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( J-1, SIGNBC, B( 1, J ), 1 )
               CALL ZSCAL( J, SIGNBC, A( 1, J ), 1 )
            ELSE
               A( J, J ) = A( J, J )*SIGNBC
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, J ), 1 )
         ELSE
            B( J, J ) = CZERO
         END IF
         ALPHA( J ) = A( J, J )
         BETA( J ) = B( J, J )
   10 CONTINUE
*
*     If IHI < ILO, skip QZ steps
*
      IF( IHI.LT.ILO )
     $   GO TO 190
*
*     MAIN QZ ITERATION LOOP
*
*     Initialize dynamic indices
*
*     Eigenvalues ILAST+1:N have been found.
*        Column operations modify rows IFRSTM:whatever
*        Row operations modify columns whatever:ILASTM
*
*     If only eigenvalues are being computed, then
*        IFRSTM is the row of the last splitting row above row ILAST;
*        this is always at least ILO.
*     IITER counts iterations since the last eigenvalue was found,
*        to tell when to use an extraordinary shift.
*     MAXIT is the maximum number of QZ sweeps allowed.
*
      ILAST = IHI
      IF( ILSCHR ) THEN
         IFRSTM = 1
         ILASTM = N
      ELSE
         IFRSTM = ILO
         ILASTM = IHI
      END IF
      IITER = 0
      ESHIFT = CZERO
      MAXIT = 30*( IHI-ILO+1 )
*
      DO 170 JITER = 1, MAXIT
*
*        Check for too many iterations.
*
         IF( JITER.GT.MAXIT )
     $      GO TO 180
*
*        Split the matrix if possible.
*
*        Two tests:
*           1: A(j,j-1)=0  or  j=ILO
*           2: B(j,j)=0
*
*        Special case: j=ILAST
*
         IF( ILAST.EQ.ILO ) THEN
            GO TO 60
         ELSE
            IF( ABS1( A( ILAST, ILAST-1 ) ).LE.ATOL ) THEN
               A( ILAST, ILAST-1 ) = CZERO
               GO TO 60
            END IF
         END IF
*
         IF( ABS( B( ILAST, ILAST ) ).LE.BTOL ) THEN
            B( ILAST, ILAST ) = CZERO
            GO TO 50
         END IF
*
*        General case: j<ILAST
*
         DO 40 J = ILAST - 1, ILO, -1
*
*           Test 1: for A(j,j-1)=0 or j=ILO
*
            IF( J.EQ.ILO ) THEN
               ILAZRO = .TRUE.
            ELSE
               IF( ABS1( A( J, J-1 ) ).LE.ATOL ) THEN
                  A( J, J-1 ) = CZERO
                  ILAZRO = .TRUE.
               ELSE
                  ILAZRO = .FALSE.
               END IF
            END IF
*
*           Test 2: for B(j,j)=0
*
            IF( ABS( B( J, J ) ).LT.BTOL ) THEN
               B( J, J ) = CZERO
*
*              Test 1a: Check for 2 consecutive small subdiagonals in A
*
               ILAZR2 = .FALSE.
               IF( .NOT.ILAZRO ) THEN
                  IF( ABS1( A( J, J-1 ) )*( ASCALE*ABS1( A( J+1,
     $                J ) ) ).LE.ABS1( A( J, J ) )*( ASCALE*ATOL ) )
     $                ILAZR2 = .TRUE.
               END IF
*
*              If both tests pass (1 & 2), i.e., the leading diagonal
*              element of B in the block is zero, split a 1x1 block off
*              at the top. (I.e., at the J-th row/column) The leading
*              diagonal element of the remainder can also be zero, so
*              this may have to be done repeatedly.
*
               IF( ILAZRO .OR. ILAZR2 ) THEN
                  DO 20 JCH = J, ILAST - 1
                     CTEMP = A( JCH, JCH )
                     CALL ZLARTG( CTEMP, A( JCH+1, JCH ), C, S,
     $                            A( JCH, JCH ) )
                     A( JCH+1, JCH ) = CZERO
                     CALL ZROT( ILASTM-JCH, A( JCH, JCH+1 ), LDA,
     $                          A( JCH+1, JCH+1 ), LDA, C, S )
                     CALL ZROT( ILASTM-JCH, B( JCH, JCH+1 ), LDB,
     $                          B( JCH+1, JCH+1 ), LDB, C, S )
                     IF( ILQ )
     $                  CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, DCONJG( S ) )
                     IF( ILAZR2 )
     $                  A( JCH, JCH-1 ) = A( JCH, JCH-1 )*C
                     ILAZR2 = .FALSE.
                     IF( ABS1( B( JCH+1, JCH+1 ) ).GE.BTOL ) THEN
                        IF( JCH+1.GE.ILAST ) THEN
                           GO TO 60
                        ELSE
                           IFIRST = JCH + 1
                           GO TO 70
                        END IF
                     END IF
                     B( JCH+1, JCH+1 ) = CZERO
   20             CONTINUE
                  GO TO 50
               ELSE
*
*                 Only test 2 passed -- chase the zero to B(ILAST,ILAST)
*                 Then process as in the case B(ILAST,ILAST)=0
*
                  DO 30 JCH = J, ILAST - 1
                     CTEMP = B( JCH, JCH+1 )
                     CALL ZLARTG( CTEMP, B( JCH+1, JCH+1 ), C, S,
     $                            B( JCH, JCH+1 ) )
                     B( JCH+1, JCH+1 ) = CZERO
                     IF( JCH.LT.ILASTM-1 )
     $                  CALL ZROT( ILASTM-JCH-1, B( JCH, JCH+2 ), LDB,
     $                             B( JCH+1, JCH+2 ), LDB, C, S )
                     CALL ZROT( ILASTM-JCH+2, A( JCH, JCH-1 ), LDA,
     $                          A( JCH+1, JCH-1 ), LDA, C, S )
                     IF( ILQ )
     $                  CALL ZROT( N, Q( 1, JCH ), 1, Q( 1, JCH+1 ), 1,
     $                             C, DCONJG( S ) )
                     CTEMP = A( JCH+1, JCH )
                     CALL ZLARTG( CTEMP, A( JCH+1, JCH-1 ), C, S,
     $                            A( JCH+1, JCH ) )
                     A( JCH+1, JCH-1 ) = CZERO
                     CALL ZROT( JCH+1-IFRSTM, A( IFRSTM, JCH ), 1,
     $                          A( IFRSTM, JCH-1 ), 1, C, S )
                     CALL ZROT( JCH-IFRSTM, B( IFRSTM, JCH ), 1,
     $                          B( IFRSTM, JCH-1 ), 1, C, S )
                     IF( ILZ )
     $                  CALL ZROT( N, Z( 1, JCH ), 1, Z( 1, JCH-1 ), 1,
     $                             C, S )
   30             CONTINUE
                  GO TO 50
               END IF
            ELSE IF( ILAZRO ) THEN
*
*              Only test 1 passed -- work on J:ILAST
*
               IFIRST = J
               GO TO 70
            END IF
*
*           Neither test passed -- try next J
*
   40    CONTINUE
*
*        (Drop-through is "impossible")
*
         INFO = 2*N + 1
         GO TO 210
*
*        B(ILAST,ILAST)=0 -- clear A(ILAST,ILAST-1) to split off a
*        1x1 block.
*
   50    CONTINUE
         CTEMP = A( ILAST, ILAST )
         CALL ZLARTG( CTEMP, A( ILAST, ILAST-1 ), C, S,
     $                A( ILAST, ILAST ) )
         A( ILAST, ILAST-1 ) = CZERO
         CALL ZROT( ILAST-IFRSTM, A( IFRSTM, ILAST ), 1,
     $              A( IFRSTM, ILAST-1 ), 1, C, S )
         CALL ZROT( ILAST-IFRSTM, B( IFRSTM, ILAST ), 1,
     $              B( IFRSTM, ILAST-1 ), 1, C, S )
         IF( ILZ )
     $      CALL ZROT( N, Z( 1, ILAST ), 1, Z( 1, ILAST-1 ), 1, C, S )
*
*        A(ILAST,ILAST-1)=0 -- Standardize B, set ALPHA and BETA
*
   60    CONTINUE
         ABSB = ABS( B( ILAST, ILAST ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( B( ILAST, ILAST ) / ABSB )
            B( ILAST, ILAST ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( ILAST-IFRSTM, SIGNBC, B( IFRSTM, ILAST ), 1 )
               CALL ZSCAL( ILAST+1-IFRSTM, SIGNBC, A( IFRSTM, ILAST ),
     $                     1 )
            ELSE
               A( ILAST, ILAST ) = A( ILAST, ILAST )*SIGNBC
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, ILAST ), 1 )
         ELSE
            B( ILAST, ILAST ) = CZERO
         END IF
         ALPHA( ILAST ) = A( ILAST, ILAST )
         BETA( ILAST ) = B( ILAST, ILAST )
*
*        Go to next block -- exit if finished.
*
         ILAST = ILAST - 1
         IF( ILAST.LT.ILO )
     $      GO TO 190
*
*        Reset counters
*
         IITER = 0
         ESHIFT = CZERO
         IF( .NOT.ILSCHR ) THEN
            ILASTM = ILAST
            IF( IFRSTM.GT.ILAST )
     $         IFRSTM = ILO
         END IF
         GO TO 160
*
*        QZ step
*
*        This iteration only involves rows/columns IFIRST:ILAST.  We
*        assume IFIRST < ILAST, and that the diagonal of B is non-zero.
*
   70    CONTINUE
         IITER = IITER + 1
         IF( .NOT.ILSCHR ) THEN
            IFRSTM = IFIRST
         END IF
*
*        Compute the Shift.
*
*        At this point, IFIRST < ILAST, and the diagonal elements of
*        B(IFIRST:ILAST,IFIRST,ILAST) are larger than BTOL (in
*        magnitude)
*
         IF( ( IITER / 10 )*10.NE.IITER ) THEN
*
*           The Wilkinson shift (AEP p.512), i.e., the eigenvalue of
*           the bottom-right 2x2 block of A inv(B) which is nearest to
*           the bottom-right element.
*
*           We factor B as U*D, where U has unit diagonals, and
*           compute (A*inv(D))*inv(U).
*
            U12 = ( BSCALE*B( ILAST-1, ILAST ) ) /
     $            ( BSCALE*B( ILAST, ILAST ) )
            AD11 = ( ASCALE*A( ILAST-1, ILAST-1 ) ) /
     $             ( BSCALE*B( ILAST-1, ILAST-1 ) )
            AD21 = ( ASCALE*A( ILAST, ILAST-1 ) ) /
     $             ( BSCALE*B( ILAST-1, ILAST-1 ) )
            AD12 = ( ASCALE*A( ILAST-1, ILAST ) ) /
     $             ( BSCALE*B( ILAST, ILAST ) )
            AD22 = ( ASCALE*A( ILAST, ILAST ) ) /
     $             ( BSCALE*B( ILAST, ILAST ) )
            ABI22 = AD22 - U12*AD21
*
            T = HALF*( AD11+ABI22 )
            RTDISC = SQRT( T**2+AD12*AD21-AD11*AD22 )
            TEMP = DBLE( T-ABI22 )*DBLE( RTDISC ) +
     $             DIMAG( T-ABI22 )*DIMAG( RTDISC )
            IF( TEMP.LE.ZERO ) THEN
               SHIFT = T + RTDISC
            ELSE
               SHIFT = T - RTDISC
            END IF
         ELSE
*
*           Exceptional shift.  Chosen for no particularly good reason.
*
            ESHIFT = ESHIFT + DCONJG( ( ASCALE*A( ILAST-1, ILAST ) ) /
     $               ( BSCALE*B( ILAST-1, ILAST-1 ) ) )
            SHIFT = ESHIFT
         END IF
*
*        Now check for two consecutive small subdiagonals.
*
         DO 80 J = ILAST - 1, IFIRST + 1, -1
            ISTART = J
            CTEMP = ASCALE*A( J, J ) - SHIFT*( BSCALE*B( J, J ) )
            TEMP = ABS1( CTEMP )
            TEMP2 = ASCALE*ABS1( A( J+1, J ) )
            TEMPR = MAX( TEMP, TEMP2 )
            IF( TEMPR.LT.ONE .AND. TEMPR.NE.ZERO ) THEN
               TEMP = TEMP / TEMPR
               TEMP2 = TEMP2 / TEMPR
            END IF
            IF( ABS1( A( J, J-1 ) )*TEMP2.LE.TEMP*ATOL )
     $         GO TO 90
   80    CONTINUE
*
         ISTART = IFIRST
         CTEMP = ASCALE*A( IFIRST, IFIRST ) -
     $           SHIFT*( BSCALE*B( IFIRST, IFIRST ) )
   90    CONTINUE
*
*        Do an implicit-shift QZ sweep.
*
*        Initial Q
*
         CTEMP2 = ASCALE*A( ISTART+1, ISTART )
         CALL ZLARTG( CTEMP, CTEMP2, C, S, CTEMP3 )
*
*        Sweep
*
         DO 150 J = ISTART, ILAST - 1
            IF( J.GT.ISTART ) THEN
               CTEMP = A( J, J-1 )
               CALL ZLARTG( CTEMP, A( J+1, J-1 ), C, S, A( J, J-1 ) )
               A( J+1, J-1 ) = CZERO
            END IF
*
            DO 100 JC = J, ILASTM
               CTEMP = C*A( J, JC ) + S*A( J+1, JC )
               A( J+1, JC ) = -DCONJG( S )*A( J, JC ) + C*A( J+1, JC )
               A( J, JC ) = CTEMP
               CTEMP2 = C*B( J, JC ) + S*B( J+1, JC )
               B( J+1, JC ) = -DCONJG( S )*B( J, JC ) + C*B( J+1, JC )
               B( J, JC ) = CTEMP2
  100       CONTINUE
            IF( ILQ ) THEN
               DO 110 JR = 1, N
                  CTEMP = C*Q( JR, J ) + DCONJG( S )*Q( JR, J+1 )
                  Q( JR, J+1 ) = -S*Q( JR, J ) + C*Q( JR, J+1 )
                  Q( JR, J ) = CTEMP
  110          CONTINUE
            END IF
*
            CTEMP = B( J+1, J+1 )
            CALL ZLARTG( CTEMP, B( J+1, J ), C, S, B( J+1, J+1 ) )
            B( J+1, J ) = CZERO
*
            DO 120 JR = IFRSTM, MIN( J+2, ILAST )
               CTEMP = C*A( JR, J+1 ) + S*A( JR, J )
               A( JR, J ) = -DCONJG( S )*A( JR, J+1 ) + C*A( JR, J )
               A( JR, J+1 ) = CTEMP
  120       CONTINUE
            DO 130 JR = IFRSTM, J
               CTEMP = C*B( JR, J+1 ) + S*B( JR, J )
               B( JR, J ) = -DCONJG( S )*B( JR, J+1 ) + C*B( JR, J )
               B( JR, J+1 ) = CTEMP
  130       CONTINUE
            IF( ILZ ) THEN
               DO 140 JR = 1, N
                  CTEMP = C*Z( JR, J+1 ) + S*Z( JR, J )
                  Z( JR, J ) = -DCONJG( S )*Z( JR, J+1 ) + C*Z( JR, J )
                  Z( JR, J+1 ) = CTEMP
  140          CONTINUE
            END IF
  150    CONTINUE
*
  160    CONTINUE
*
  170 CONTINUE
*
*     Drop-through = non-convergence
*
  180 CONTINUE
      INFO = ILAST
      GO TO 210
*
*     Successful completion of all QZ steps
*
  190 CONTINUE
*
*     Set Eigenvalues 1:ILO-1
*
      DO 200 J = 1, ILO - 1
         ABSB = ABS( B( J, J ) )
         IF( ABSB.GT.SAFMIN ) THEN
            SIGNBC = DCONJG( B( J, J ) / ABSB )
            B( J, J ) = ABSB
            IF( ILSCHR ) THEN
               CALL ZSCAL( J-1, SIGNBC, B( 1, J ), 1 )
               CALL ZSCAL( J, SIGNBC, A( 1, J ), 1 )
            ELSE
               A( J, J ) = A( J, J )*SIGNBC
            END IF
            IF( ILZ )
     $         CALL ZSCAL( N, SIGNBC, Z( 1, J ), 1 )
         ELSE
            B( J, J ) = CZERO
         END IF
         ALPHA( J ) = A( J, J )
         BETA( J ) = B( J, J )
  200 CONTINUE
*
*     Normal Termination
*
      INFO = 0
*
*     Exit (other than argument error) -- return optimal workspace size
*
  210 CONTINUE
      WORK( 1 ) = DCMPLX( N )
      RETURN
*
*     End of ZHGEQZ
*
      END

      SUBROUTINE ZLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   CS
      COMPLEX*16         F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  ZLARTG generates a plane rotation so that
*
*     [  CS  SN  ]     [ F ]     [ R ]
*     [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a faster version of the BLAS1 routine ZROTG, except for
*  the following differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0, then CS=0 and SN is chosen so that R is real.
*
*  Arguments
*  =========
*
*  F       (input) COMPLEX*16
*          The first component of vector to be rotated.
*
*  G       (input) COMPLEX*16
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) COMPLEX*16
*          The sine of the rotation.
*
*  R       (output) COMPLEX*16
*          The nonzero component of the rotated vector.
*
*  Further Details
*  ======= =======
*
*  3-5-96 - Modified with a new algorithm by W. Kahan and J. Demmel
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   TWO, ONE, ZERO
      PARAMETER          ( TWO = 2.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            FIRST
      INTEGER            COUNT, I
      DOUBLE PRECISION   D, DI, DR, EPS, F2, F2S, G2, G2S, SAFMIN,
     $                   SAFMN2, SAFMX2, SCALE
      COMPLEX*16         FF, FS, GS
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, INT, LOG,
     $                   MAX, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1, ABSSQ
*     ..
*     .. Save statement ..
      SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
      DATA               FIRST / .TRUE. /
*     ..
*     .. Statement Function definitions ..
      ABS1( FF ) = MAX( ABS( DBLE( FF ) ), ABS( DIMAG( FF ) ) )
      ABSSQ( FF ) = DBLE( FF )**2 + DIMAG( FF )**2
*     ..
*     .. Executable Statements ..
*
      IF( FIRST ) THEN
         FIRST = .FALSE.
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
      END IF
      SCALE = MAX( ABS1( F ), ABS1( G ) )
      FS = F
      GS = G
      COUNT = 0
      IF( SCALE.GE.SAFMX2 ) THEN
   10    CONTINUE
         COUNT = COUNT + 1
         FS = FS*SAFMN2
         GS = GS*SAFMN2
         SCALE = SCALE*SAFMN2
         IF( SCALE.GE.SAFMX2 )
     $      GO TO 10
      ELSE IF( SCALE.LE.SAFMN2 ) THEN
         IF( G.EQ.CZERO ) THEN
            CS = ONE
            SN = CZERO
            R = F
            RETURN
         END IF
   20    CONTINUE
         COUNT = COUNT - 1
         FS = FS*SAFMX2
         GS = GS*SAFMX2
         SCALE = SCALE*SAFMX2
         IF( SCALE.LE.SAFMN2 )
     $      GO TO 20
      END IF
      F2 = ABSSQ( FS )
      G2 = ABSSQ( GS )
      IF( F2.LE.MAX( G2, ONE )*SAFMIN ) THEN
*
*        This is a rare case: F is very small.
*
         IF( F.EQ.CZERO ) THEN
            CS = ZERO
            R = DLAPY2( DBLE( G ), DIMAG( G ) )
*           Do complex/real division explicitly with two real divisions
            D = DLAPY2( DBLE( GS ), DIMAG( GS ) )
            SN = DCMPLX( DBLE( GS ) / D, -DIMAG( GS ) / D )
            RETURN
         END IF
         F2S = DLAPY2( DBLE( FS ), DIMAG( FS ) )
*        G2 and G2S are accurate
*        G2 is at least SAFMIN, and G2S is at least SAFMN2
         G2S = SQRT( G2 )
*        Error in CS from underflow in F2S is at most
*        UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS
*        If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
*        and so CS .lt. sqrt(SAFMIN)
*        If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
*        and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
*        Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
         CS = F2S / G2S
*        Make sure abs(FF) = 1
*        Do complex/real division explicitly with 2 real divisions
         IF( ABS1( F ).GT.ONE ) THEN
            D = DLAPY2( DBLE( F ), DIMAG( F ) )
            FF = DCMPLX( DBLE( F ) / D, DIMAG( F ) / D )
         ELSE
            DR = SAFMX2*DBLE( F )
            DI = SAFMX2*DIMAG( F )
            D = DLAPY2( DR, DI )
            FF = DCMPLX( DR / D, DI / D )
         END IF
         SN = FF*DCMPLX( DBLE( GS ) / G2S, -DIMAG( GS ) / G2S )
         R = CS*F + SN*G
      ELSE
*
*        This is the most common case.
*        Neither F2 nor F2/G2 are less than SAFMIN
*        F2S cannot overflow, and it is accurate
*
         F2S = SQRT( ONE+G2 / F2 )
*        Do the F2S(real)*FS(complex) multiply with two real multiplies
         R = DCMPLX( F2S*DBLE( FS ), F2S*DIMAG( FS ) )
         CS = ONE / F2S
         D = F2 + G2
*        Do complex/real division explicitly with two real divisions
         SN = DCMPLX( DBLE( R ) / D, DIMAG( R ) / D )
         SN = SN*DCONJG( GS )
         IF( COUNT.NE.0 ) THEN
            IF( COUNT.GT.0 ) THEN
               DO 30 I = 1, COUNT
                  R = R*SAFMX2
   30          CONTINUE
            ELSE
               DO 40 I = 1, -COUNT
                  R = R*SAFMN2
   40          CONTINUE
            END IF
         END IF
      END IF
      RETURN
*
*     End of ZLARTG
*
      END

      DOUBLE PRECISION FUNCTION ZLANHS( NORM, N, A, LDA, WORK )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          NORM
      INTEGER            LDA, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLANHS  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  Hessenberg matrix A.
*
*  Description
*  ===========
*
*  ZLANHS returns the value
*
*     ZLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in ZLANHS as described
*          above.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, ZLANHS is
*          set to zero.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The n by n upper Hessenberg matrix A; the part of A below the
*          first sub-diagonal is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
*          referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J
      DOUBLE PRECISION   SCALE, SUM, VALUE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         DO 20 J = 1, N
            DO 10 I = 1, MIN( N, J+1 )
               VALUE = MAX( VALUE, ABS( A( I, J ) ) )
   10       CONTINUE
   20    CONTINUE
      ELSE IF( ( LSAME( NORM, 'O' ) ) .OR. ( NORM.EQ.'1' ) ) THEN
*
*        Find norm1(A).
*
         VALUE = ZERO
         DO 40 J = 1, N
            SUM = ZERO
            DO 30 I = 1, MIN( N, J+1 )
               SUM = SUM + ABS( A( I, J ) )
   30       CONTINUE
            VALUE = MAX( VALUE, SUM )
   40    CONTINUE
      ELSE IF( LSAME( NORM, 'I' ) ) THEN
*
*        Find normI(A).
*
         DO 50 I = 1, N
            WORK( I ) = ZERO
   50    CONTINUE
         DO 70 J = 1, N
            DO 60 I = 1, MIN( N, J+1 )
               WORK( I ) = WORK( I ) + ABS( A( I, J ) )
   60       CONTINUE
   70    CONTINUE
         VALUE = ZERO
         DO 80 I = 1, N
            VALUE = MAX( VALUE, WORK( I ) )
   80    CONTINUE
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         DO 90 J = 1, N
            CALL ZLASSQ( MIN( N, J+1 ), A( 1, J ), 1, SCALE, SUM )
   90    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      ZLANHS = VALUE
      RETURN
*
*     End of ZLANHS
*
      END

      SUBROUTINE ZTGEVC( SIDE, HOWMNY, SELECT, N, A, LDA, B, LDB, VL,
     $                   LDVL, VR, LDVR, MM, M, WORK, RWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX*16         A( LDA, * ), B( LDB, * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
*     ..
*
*
*  Purpose
*  =======
*
*  ZTGEVC computes some or all of the right and/or left generalized
*  eigenvectors of a pair of complex upper triangular matrices (A,B).
*
*  The right generalized eigenvector x and the left generalized
*  eigenvector y of (A,B) corresponding to a generalized eigenvalue
*  w are defined by:
*
*          (A - wB) * x = 0  and  y**H * (A - wB) = 0
*
*  where y**H denotes the conjugate tranpose of y.
*
*  If an eigenvalue w is determined by zero diagonal elements of both A
*  and B, a unit vector is returned as the corresponding eigenvector.
*
*  If all eigenvectors are requested, the routine may either return
*  the matrices X and/or Y of right or left eigenvectors of (A,B), or
*  the products Z*X and/or Q*Y, where Z and Q are input unitary
*  matrices.  If (A,B) was obtained from the generalized Schur
*  factorization of an original pair of matrices
*     (A0,B0) = (Q*A*Z**H,Q*B*Z**H),
*  then Z*X and Q*Y are the matrices of right or left eigenvectors of
*  A.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'R': compute right eigenvectors only;
*          = 'L': compute left eigenvectors only;
*          = 'B': compute both right and left eigenvectors.
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A': compute all right and/or left eigenvectors;
*          = 'B': compute all right and/or left eigenvectors, and
*                 backtransform them using the input matrices supplied
*                 in VR and/or VL;
*          = 'S': compute selected right and/or left eigenvectors,
*                 specified by the logical array SELECT.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          If HOWMNY='S', SELECT specifies the eigenvectors to be
*          computed.
*          If HOWMNY='A' or 'B', SELECT is not referenced.
*          To select the eigenvector corresponding to the j-th
*          eigenvalue, SELECT(j) must be set to .TRUE..
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The upper triangular matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of array A.  LDA >= max(1,N).
*
*  B       (input) COMPLEX*16 array, dimension (LDB,N)
*          The upper triangular matrix B.  B must have real diagonal
*          elements.
*
*  LDB     (input) INTEGER
*          The leading dimension of array B.  LDB >= max(1,N).
*
*  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM)
*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*          contain an N-by-N matrix Q (usually the unitary matrix Q
*          of left Schur vectors returned by ZHGEQZ).
*          On exit, if SIDE = 'L' or 'B', VL contains:
*          if HOWMNY = 'A', the matrix Y of left eigenvectors of (A,B);
*          if HOWMNY = 'B', the matrix Q*Y;
*          if HOWMNY = 'S', the left eigenvectors of (A,B) specified by
*                      SELECT, stored consecutively in the columns of
*                      VL, in the same order as their eigenvalues.
*          If SIDE = 'R', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of array VL.
*          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.
*
*  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)
*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*          contain an N-by-N matrix Q (usually the unitary matrix Z
*          of right Schur vectors returned by ZHGEQZ).
*          On exit, if SIDE = 'R' or 'B', VR contains:
*          if HOWMNY = 'A', the matrix X of right eigenvectors of (A,B);
*          if HOWMNY = 'B', the matrix Z*X;
*          if HOWMNY = 'S', the right eigenvectors of (A,B) specified by
*                      SELECT, stored consecutively in the columns of
*                      VR, in the same order as their eigenvalues.
*          If SIDE = 'L', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.
*          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.
*
*  MM      (input) INTEGER
*          The number of columns in the arrays VL and/or VR. MM >= M.
*
*  M       (output) INTEGER
*          The number of columns in the arrays VL and/or VR actually
*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
*          is set to N.  Each selected eigenvector occupies one column.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (2*N)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (2*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            COMPL, COMPR, ILALL, ILBACK, ILBBAD, ILCOMP,
     $                   LSA, LSB
      INTEGER            I, IBEG, IEIG, IEND, IHWMNY, IM, ISIDE, ISRC,
     $                   J, JE, JR
      DOUBLE PRECISION   ACOEFA, ACOEFF, ANORM, ASCALE, BCOEFA, BIG,
     $                   BIGNUM, BNORM, BSCALE, DMIN, SAFMIN, SBETA,
     $                   SCALE, SMALL, TEMP, ULP, XMAX
      COMPLEX*16         BCOEFF, CA, CB, D, SALPHA, SUM, SUMA, SUMB, X
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      COMPLEX*16         ZLADIV
      EXTERNAL           LSAME, DLAMCH, ZLADIV
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLABAD, XERBLA, ZGEMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, MAX, MIN
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1
*     ..
*     .. Statement Function definitions ..
      ABS1( X ) = ABS( DBLE( X ) ) + ABS( DIMAG( X ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and Test the input parameters
*
      IF( LSAME( HOWMNY, 'A' ) ) THEN
         IHWMNY = 1
         ILALL = .TRUE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( HOWMNY, 'S' ) ) THEN
         IHWMNY = 2
         ILALL = .FALSE.
         ILBACK = .FALSE.
      ELSE IF( LSAME( HOWMNY, 'B' ) .OR. LSAME( HOWMNY, 'T' ) ) THEN
         IHWMNY = 3
         ILALL = .TRUE.
         ILBACK = .TRUE.
      ELSE
         IHWMNY = -1
      END IF
*
      IF( LSAME( SIDE, 'R' ) ) THEN
         ISIDE = 1
         COMPL = .FALSE.
         COMPR = .TRUE.
      ELSE IF( LSAME( SIDE, 'L' ) ) THEN
         ISIDE = 2
         COMPL = .TRUE.
         COMPR = .FALSE.
      ELSE IF( LSAME( SIDE, 'B' ) ) THEN
         ISIDE = 3
         COMPL = .TRUE.
         COMPR = .TRUE.
      ELSE
         ISIDE = -1
      END IF
*
      INFO = 0
      IF( ISIDE.LT.0 ) THEN
         INFO = -1
      ELSE IF( IHWMNY.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDB.LT.MAX( 1, N ) ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGEVC', -INFO )
         RETURN
      END IF
*
*     Count the number of eigenvectors
*
      IF( .NOT.ILALL ) THEN
         IM = 0
         DO 10 J = 1, N
            IF( SELECT( J ) )
     $         IM = IM + 1
   10    CONTINUE
      ELSE
         IM = N
      END IF
*
*     Check diagonal of B
*
      ILBBAD = .FALSE.
      DO 20 J = 1, N
         IF( DIMAG( B( J, J ) ).NE.ZERO )
     $      ILBBAD = .TRUE.
   20 CONTINUE
*
      IF( ILBBAD ) THEN
         INFO = -7
      ELSE IF( COMPL .AND. LDVL.LT.N .OR. LDVL.LT.1 ) THEN
         INFO = -10
      ELSE IF( COMPR .AND. LDVR.LT.N .OR. LDVR.LT.1 ) THEN
         INFO = -12
      ELSE IF( MM.LT.IM ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTGEVC', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      M = IM
      IF( N.EQ.0 )
     $   RETURN
*
*     Machine Constants
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      BIG = ONE / SAFMIN
      CALL DLABAD( SAFMIN, BIG )
      ULP = DLAMCH( 'Epsilon' )*DLAMCH( 'Base' )
      SMALL = SAFMIN*N / ULP
      BIG = ONE / SMALL
      BIGNUM = ONE / ( SAFMIN*N )
*
*     Compute the 1-norm of each column of the strictly upper triangular
*     part of A and B to check for possible overflow in the triangular
*     solver.
*
      ANORM = ABS1( A( 1, 1 ) )
      BNORM = ABS1( B( 1, 1 ) )
      RWORK( 1 ) = ZERO
      RWORK( N+1 ) = ZERO
      DO 40 J = 2, N
         RWORK( J ) = ZERO
         RWORK( N+J ) = ZERO
         DO 30 I = 1, J - 1
            RWORK( J ) = RWORK( J ) + ABS1( A( I, J ) )
            RWORK( N+J ) = RWORK( N+J ) + ABS1( B( I, J ) )
   30    CONTINUE
         ANORM = MAX( ANORM, RWORK( J )+ABS1( A( J, J ) ) )
         BNORM = MAX( BNORM, RWORK( N+J )+ABS1( B( J, J ) ) )
   40 CONTINUE
*
      ASCALE = ONE / MAX( ANORM, SAFMIN )
      BSCALE = ONE / MAX( BNORM, SAFMIN )
*
*     Left eigenvectors
*
      IF( COMPL ) THEN
         IEIG = 0
*
*        Main loop over eigenvalues
*
         DO 140 JE = 1, N
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( ILCOMP ) THEN
               IEIG = IEIG + 1
*
               IF( ABS1( A( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( DBLE( B( JE, JE ) ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- return unit eigenvector
*
                  DO 50 JR = 1, N
                     VL( JR, IEIG ) = CZERO
   50             CONTINUE
                  VL( IEIG, IEIG ) = CONE
                  GO TO 140
               END IF
*
*              Non-singular eigenvalue:
*              Compute coefficients  a  and  b  in
*                   H
*                 y  ( a A - b B ) = 0
*
               TEMP = ONE / MAX( ABS1( A( JE, JE ) )*ASCALE,
     $                ABS( DBLE( B( JE, JE ) ) )*BSCALE, SAFMIN )
               SALPHA = ( TEMP*A( JE, JE ) )*ASCALE
               SBETA = ( TEMP*DBLE( B( JE, JE ) ) )*BSCALE
               ACOEFF = SBETA*ASCALE
               BCOEFF = SALPHA*BSCALE
*
*              Scale to avoid underflow
*
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEFF ).LT.SMALL
               LSB = ABS1( SALPHA ).GE.SAFMIN .AND. ABS1( BCOEFF ).LT.
     $               SMALL
*
               SCALE = ONE
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEFF ),
     $                    ABS1( BCOEFF ) ) ) )
                  IF( LSA ) THEN
                     ACOEFF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEFF = SCALE*ACOEFF
                  END IF
                  IF( LSB ) THEN
                     BCOEFF = BSCALE*( SCALE*SALPHA )
                  ELSE
                     BCOEFF = SCALE*BCOEFF
                  END IF
               END IF
*
               ACOEFA = ABS( ACOEFF )
               BCOEFA = ABS1( BCOEFF )
               XMAX = ONE
               DO 60 JR = 1, N
                  WORK( JR ) = CZERO
   60          CONTINUE
               WORK( JE ) = CONE
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*                                              H
*              Triangular solve of  (a A - b B)  y = 0
*
*                                      H
*              (rowwise in  (a A - b B) , or columnwise in a A - b B)
*
               DO 100 J = JE + 1, N
*
*                 Compute
*                       j-1
*                 SUM = sum  conjg( a*A(k,j) - b*B(k,j) )*x(k)
*                       k=je
*                 (Scale if necessary)
*
                  TEMP = ONE / XMAX
                  IF( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ).GT.BIGNUM*
     $                TEMP ) THEN
                     DO 70 JR = JE, J - 1
                        WORK( JR ) = TEMP*WORK( JR )
   70                CONTINUE
                     XMAX = ONE
                  END IF
                  SUMA = CZERO
                  SUMB = CZERO
*
                  DO 80 JR = JE, J - 1
                     SUMA = SUMA + DCONJG( A( JR, J ) )*WORK( JR )
                     SUMB = SUMB + DCONJG( B( JR, J ) )*WORK( JR )
   80             CONTINUE
                  SUM = ACOEFF*SUMA - DCONJG( BCOEFF )*SUMB
*
*                 Form x(j) = - SUM / conjg( a*A(j,j) - b*B(j,j) )
*
*                 with scaling and perturbation of the denominator
*
                  D = DCONJG( ACOEFF*A( J, J )-BCOEFF*B( J, J ) )
                  IF( ABS1( D ).LE.DMIN )
     $               D = DCMPLX( DMIN )
*
                  IF( ABS1( D ).LT.ONE ) THEN
                     IF( ABS1( SUM ).GE.BIGNUM*ABS1( D ) ) THEN
                        TEMP = ONE / ABS1( SUM )
                        DO 90 JR = JE, J - 1
                           WORK( JR ) = TEMP*WORK( JR )
   90                   CONTINUE
                        XMAX = TEMP*XMAX
                        SUM = TEMP*SUM
                     END IF
                  END IF
                  WORK( J ) = ZLADIV( -SUM, D )
                  XMAX = MAX( XMAX, ABS1( WORK( J ) ) )
  100          CONTINUE
*
*              Back transform eigenvector if HOWMNY='B'.
*
               IF( ILBACK ) THEN
                  CALL ZGEMV( 'N', N, N+1-JE, CONE, VL( 1, JE ), LDVL,
     $                        WORK( JE ), 1, CZERO, WORK( N+1 ), 1 )
                  ISRC = 2
                  IBEG = 1
               ELSE
                  ISRC = 1
                  IBEG = JE
               END IF
*
*              Copy and scale eigenvector into column of VL
*
               XMAX = ZERO
               DO 110 JR = IBEG, N
                  XMAX = MAX( XMAX, ABS1( WORK( ( ISRC-1 )*N+JR ) ) )
  110          CONTINUE
*
               IF( XMAX.GT.SAFMIN ) THEN
                  TEMP = ONE / XMAX
                  DO 120 JR = IBEG, N
                     VL( JR, IEIG ) = TEMP*WORK( ( ISRC-1 )*N+JR )
  120             CONTINUE
               ELSE
                  IBEG = N + 1
               END IF
*
               DO 130 JR = 1, IBEG - 1
                  VL( JR, IEIG ) = CZERO
  130          CONTINUE
*
            END IF
  140    CONTINUE
      END IF
*
*     Right eigenvectors
*
      IF( COMPR ) THEN
         IEIG = IM + 1
*
*        Main loop over eigenvalues
*
         DO 250 JE = N, 1, -1
            IF( ILALL ) THEN
               ILCOMP = .TRUE.
            ELSE
               ILCOMP = SELECT( JE )
            END IF
            IF( ILCOMP ) THEN
               IEIG = IEIG - 1
*
               IF( ABS1( A( JE, JE ) ).LE.SAFMIN .AND.
     $             ABS( DBLE( B( JE, JE ) ) ).LE.SAFMIN ) THEN
*
*                 Singular matrix pencil -- return unit eigenvector
*
                  DO 150 JR = 1, N
                     VR( JR, IEIG ) = CZERO
  150             CONTINUE
                  VR( IEIG, IEIG ) = CONE
                  GO TO 250
               END IF
*
*              Non-singular eigenvalue:
*              Compute coefficients  a  and  b  in
*
*              ( a A - b B ) x  = 0
*
               TEMP = ONE / MAX( ABS1( A( JE, JE ) )*ASCALE,
     $                ABS( DBLE( B( JE, JE ) ) )*BSCALE, SAFMIN )
               SALPHA = ( TEMP*A( JE, JE ) )*ASCALE
               SBETA = ( TEMP*DBLE( B( JE, JE ) ) )*BSCALE
               ACOEFF = SBETA*ASCALE
               BCOEFF = SALPHA*BSCALE
*
*              Scale to avoid underflow
*
               LSA = ABS( SBETA ).GE.SAFMIN .AND. ABS( ACOEFF ).LT.SMALL
               LSB = ABS1( SALPHA ).GE.SAFMIN .AND. ABS1( BCOEFF ).LT.
     $               SMALL
*
               SCALE = ONE
               IF( LSA )
     $            SCALE = ( SMALL / ABS( SBETA ) )*MIN( ANORM, BIG )
               IF( LSB )
     $            SCALE = MAX( SCALE, ( SMALL / ABS1( SALPHA ) )*
     $                    MIN( BNORM, BIG ) )
               IF( LSA .OR. LSB ) THEN
                  SCALE = MIN( SCALE, ONE /
     $                    ( SAFMIN*MAX( ONE, ABS( ACOEFF ),
     $                    ABS1( BCOEFF ) ) ) )
                  IF( LSA ) THEN
                     ACOEFF = ASCALE*( SCALE*SBETA )
                  ELSE
                     ACOEFF = SCALE*ACOEFF
                  END IF
                  IF( LSB ) THEN
                     BCOEFF = BSCALE*( SCALE*SALPHA )
                  ELSE
                     BCOEFF = SCALE*BCOEFF
                  END IF
               END IF
*
               ACOEFA = ABS( ACOEFF )
               BCOEFA = ABS1( BCOEFF )
               XMAX = ONE
               DO 160 JR = 1, N
                  WORK( JR ) = CZERO
  160          CONTINUE
               WORK( JE ) = CONE
               DMIN = MAX( ULP*ACOEFA*ANORM, ULP*BCOEFA*BNORM, SAFMIN )
*
*              Triangular solve of  (a A - b B) x = 0  (columnwise)
*
*              WORK(1:j-1) contains sums w,
*              WORK(j+1:JE) contains x
*
               DO 170 JR = 1, JE - 1
                  WORK( JR ) = ACOEFF*A( JR, JE ) - BCOEFF*B( JR, JE )
  170          CONTINUE
               WORK( JE ) = CONE
*
               DO 210 J = JE - 1, 1, -1
*
*                 Form x(j) := - w(j) / d
*                 with scaling and perturbation of the denominator
*
                  D = ACOEFF*A( J, J ) - BCOEFF*B( J, J )
                  IF( ABS1( D ).LE.DMIN )
     $               D = DCMPLX( DMIN )
*
                  IF( ABS1( D ).LT.ONE ) THEN
                     IF( ABS1( WORK( J ) ).GE.BIGNUM*ABS1( D ) ) THEN
                        TEMP = ONE / ABS1( WORK( J ) )
                        DO 180 JR = 1, JE
                           WORK( JR ) = TEMP*WORK( JR )
  180                   CONTINUE
                     END IF
                  END IF
*
                  WORK( J ) = ZLADIV( -WORK( J ), D )
*
                  IF( J.GT.1 ) THEN
*
*                    w = w + x(j)*(a A(*,j) - b B(*,j) ) with scaling
*
                     IF( ABS1( WORK( J ) ).GT.ONE ) THEN
                        TEMP = ONE / ABS1( WORK( J ) )
                        IF( ACOEFA*RWORK( J )+BCOEFA*RWORK( N+J ).GE.
     $                      BIGNUM*TEMP ) THEN
                           DO 190 JR = 1, JE
                              WORK( JR ) = TEMP*WORK( JR )
  190                      CONTINUE
                        END IF
                     END IF
*
                     CA = ACOEFF*WORK( J )
                     CB = BCOEFF*WORK( J )
                     DO 200 JR = 1, J - 1
                        WORK( JR ) = WORK( JR ) + CA*A( JR, J ) -
     $                               CB*B( JR, J )
  200                CONTINUE
                  END IF
  210          CONTINUE
*
*              Back transform eigenvector if HOWMNY='B'.
*
               IF( ILBACK ) THEN
                  CALL ZGEMV( 'N', N, JE, CONE, VR, LDVR, WORK, 1,
     $                        CZERO, WORK( N+1 ), 1 )
                  ISRC = 2
                  IEND = N
               ELSE
                  ISRC = 1
                  IEND = JE
               END IF
*
*              Copy and scale eigenvector into column of VR
*
               XMAX = ZERO
               DO 220 JR = 1, IEND
                  XMAX = MAX( XMAX, ABS1( WORK( ( ISRC-1 )*N+JR ) ) )
  220          CONTINUE
*
               IF( XMAX.GT.SAFMIN ) THEN
                  TEMP = ONE / XMAX
                  DO 230 JR = 1, IEND
                     VR( JR, IEIG ) = TEMP*WORK( ( ISRC-1 )*N+JR )
  230             CONTINUE
               ELSE
                  IEND = 0
               END IF
*
               DO 240 JR = IEND + 1, N
                  VR( JR, IEIG ) = CZERO
  240          CONTINUE
*
            END IF
  250    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTGEVC
*
      END

      SUBROUTINE ZGGBAK( JOB, SIDE, N, ILO, IHI, LSCALE, RSCALE, M, V,
     $                   LDV, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          JOB, SIDE
      INTEGER            IHI, ILO, INFO, LDV, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   LSCALE( * ), RSCALE( * )
      COMPLEX*16         V( LDV, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGGBAK forms the right or left eigenvectors of a complex generalized
*  eigenvalue problem A*x = lambda*B*x, by backward transformation on
*  the computed eigenvectors of the balanced pair of matrices output by
*  ZGGBAL.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies the type of backward transformation required:
*          = 'N':  do nothing, return immediately;
*          = 'P':  do backward transformation for permutation only;
*          = 'S':  do backward transformation for scaling only;
*          = 'B':  do backward transformations for both permutation and
*                  scaling.
*          JOB must be the same as the argument JOB supplied to ZGGBAL.
*
*  SIDE    (input) CHARACTER*1
*          = 'R':  V contains right eigenvectors;
*          = 'L':  V contains left eigenvectors.
*
*  N       (input) INTEGER
*          The number of rows of the matrix V.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          The integers ILO and IHI determined by ZGGBAL.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  LSCALE  (input) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and/or scaling factors applied
*          to the left side of A and B, as returned by ZGGBAL.
*
*  RSCALE  (input) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and/or scaling factors applied
*          to the right side of A and B, as returned by ZGGBAL.
*
*  M       (input) INTEGER
*          The number of columns of the matrix V.  M >= 0.
*
*  V       (input/output) COMPLEX*16 array, dimension (LDV,M)
*          On entry, the matrix of right or left eigenvectors to be
*          transformed, as returned by ZTGEVC.
*          On exit, V is overwritten by the transformed eigenvectors.
*
*  LDV     (input) INTEGER
*          The leading dimension of the matrix V. LDV >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  See R.C. Ward, Balancing the generalized eigenvalue problem,
*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV
      INTEGER            I, K
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZDSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      RIGHTV = LSAME( SIDE, 'R' )
      LEFTV = LSAME( SIDE, 'L' )
*
      INFO = 0
      IF( .NOT.LSAME( JOB, 'N' ) .AND. .NOT.LSAME( JOB, 'P' ) .AND.
     $    .NOT.LSAME( JOB, 'S' ) .AND. .NOT.LSAME( JOB, 'B' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( ILO.LT.1 ) THEN
         INFO = -4
      ELSE IF( IHI.LT.ILO .OR. IHI.GT.MAX( 1, N ) ) THEN
         INFO = -5
      ELSE IF( M.LT.0 ) THEN
         INFO = -6
      ELSE IF( LDV.LT.MAX( 1, N ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGBAK', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
      IF( M.EQ.0 )
     $   RETURN
      IF( LSAME( JOB, 'N' ) )
     $   RETURN
*
      IF( ILO.EQ.IHI )
     $   GO TO 30
*
*     Backward balance
*
      IF( LSAME( JOB, 'S' ) .OR. LSAME( JOB, 'B' ) ) THEN
*
*        Backward transformation on right eigenvectors
*
         IF( RIGHTV ) THEN
            DO 10 I = ILO, IHI
               CALL ZDSCAL( M, RSCALE( I ), V( I, 1 ), LDV )
   10       CONTINUE
         END IF
*
*        Backward transformation on left eigenvectors
*
         IF( LEFTV ) THEN
            DO 20 I = ILO, IHI
               CALL ZDSCAL( M, LSCALE( I ), V( I, 1 ), LDV )
   20       CONTINUE
         END IF
      END IF
*
*     Backward permutation
*
   30 CONTINUE
      IF( LSAME( JOB, 'P' ) .OR. LSAME( JOB, 'B' ) ) THEN
*
*        Backward permutation on right eigenvectors
*
         IF( RIGHTV ) THEN
            IF( ILO.EQ.1 )
     $         GO TO 50
            DO 40 I = ILO - 1, 1, -1
               K = RSCALE( I )
               IF( K.EQ.I )
     $            GO TO 40
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   40       CONTINUE
*
   50       CONTINUE
            IF( IHI.EQ.N )
     $         GO TO 70
            DO 60 I = IHI + 1, N
               K = RSCALE( I )
               IF( K.EQ.I )
     $            GO TO 60
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   60       CONTINUE
         END IF
*
*        Backward permutation on left eigenvectors
*
   70    CONTINUE
         IF( LEFTV ) THEN
            IF( ILO.EQ.1 )
     $         GO TO 90
            DO 80 I = ILO - 1, 1, -1
               K = LSCALE( I )
               IF( K.EQ.I )
     $            GO TO 80
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
   80       CONTINUE
*
   90       CONTINUE
            IF( IHI.EQ.N )
     $         GO TO 110
            DO 100 I = IHI + 1, N
               K = LSCALE( I )
               IF( K.EQ.I )
     $            GO TO 100
               CALL ZSWAP( M, V( I, 1 ), LDV, V( K, 1 ), LDV )
  100       CONTINUE
         END IF
      END IF
*
  110 CONTINUE
*
      RETURN
*
*     End of ZGGBAK
*
      END

      SUBROUTINE ZLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLACPY copies all or part of a two-dimensional matrix A to another
*  matrix B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the matrix A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) COMPLEX*16 array, dimension (LDA,N)
*          The m by n matrix A.  If UPLO = 'U', only the upper trapezium
*          is accessed; if UPLO = 'L', only the lower trapezium is
*          accessed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (output) COMPLEX*16 array, dimension (LDB,N)
*          On exit, B = A in the locations specified by UPLO.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, J
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
      IF( LSAME( UPLO, 'U' ) ) THEN
         DO 20 J = 1, N
            DO 10 I = 1, MIN( J, M )
               B( I, J ) = A( I, J )
   10       CONTINUE
   20    CONTINUE
*
      ELSE IF( LSAME( UPLO, 'L' ) ) THEN
         DO 40 J = 1, N
            DO 30 I = J, M
               B( I, J ) = A( I, J )
   30       CONTINUE
   40    CONTINUE
*
      ELSE
         DO 60 J = 1, N
            DO 50 I = 1, M
               B( I, J ) = A( I, J )
   50       CONTINUE
   60    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLACPY
*
      END

      SUBROUTINE ZGESVD( JOBU, JOBVT, M, N, A, LDA, S, U, LDU, VT, LDVT,
     $                   WORK, LWORK, RWORK, INFO )
*
*  -- LAPACK driver routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          JOBU, JOBVT
      INTEGER            INFO, LDA, LDU, LDVT, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), S( * )
      COMPLEX*16         A( LDA, * ), U( LDU, * ), VT( LDVT, * ),
     $                   WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGESVD computes the singular value decomposition (SVD) of a complex
*  M-by-N matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written
*
*       A = U * SIGMA * conjugate-transpose(V)
*
*  where SIGMA is an M-by-N matrix which is zero except for its
*  min(m,n) diagonal elements, U is an M-by-M unitary matrix, and
*  V is an N-by-N unitary matrix.  The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and
*  are returned in descending order.  The first min(m,n) columns of
*  U and V are the left and right singular vectors of A.
*
*  Note that the routine returns V**H, not V.
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          Specifies options for computing all or part of the matrix U:
*          = 'A':  all M columns of U are returned in array U:
*          = 'S':  the first min(m,n) columns of U (the left singular
*                  vectors) are returned in the array U;
*          = 'O':  the first min(m,n) columns of U (the left singular
*                  vectors) are overwritten on the array A;
*          = 'N':  no columns of U (no left singular vectors) are
*                  computed.
*
*  JOBVT   (input) CHARACTER*1
*          Specifies options for computing all or part of the matrix
*          V**H:
*          = 'A':  all N rows of V**H are returned in the array VT;
*          = 'S':  the first min(m,n) rows of V**H (the right singular
*                  vectors) are returned in the array VT;
*          = 'O':  the first min(m,n) rows of V**H (the right singular
*                  vectors) are overwritten on the array A;
*          = 'N':  no rows of V**H (no right singular vectors) are
*                  computed.
*
*          JOBVT and JOBU cannot both be 'O'.
*
*  M       (input) INTEGER
*          The number of rows of the input matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the input matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit,
*          if JOBU = 'O',  A is overwritten with the first min(m,n)
*                          columns of U (the left singular vectors,
*                          stored columnwise);
*          if JOBVT = 'O', A is overwritten with the first min(m,n)
*                          rows of V**H (the right singular vectors,
*                          stored rowwise);
*          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
*                          are destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The singular values of A, sorted so that S(i) >= S(i+1).
*
*  U       (output) COMPLEX*16 array, dimension (LDU,UCOL)
*          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
*          If JOBU = 'A', U contains the M-by-M unitary matrix U;
*          if JOBU = 'S', U contains the first min(m,n) columns of U
*          (the left singular vectors, stored columnwise);
*          if JOBU = 'N' or 'O', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= 1; if
*          JOBU = 'S' or 'A', LDU >= M.
*
*  VT      (output) COMPLEX*16 array, dimension (LDVT,N)
*          If JOBVT = 'A', VT contains the N-by-N unitary matrix
*          V**H;
*          if JOBVT = 'S', VT contains the first min(m,n) rows of
*          V**H (the right singular vectors, stored rowwise);
*          if JOBVT = 'N' or 'O', VT is not referenced.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT.  LDVT >= 1; if
*          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= 1.
*          LWORK >=  2*MIN(M,N)+MAX(M,N).
*          For good performance, LWORK should generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (5*min(M,N))
*          On exit, if INFO > 0, RWORK(1:MIN(M,N)-1) contains the
*          unconverged superdiagonal elements of an upper bidiagonal
*          matrix B whose diagonal is in S (not necessarily sorted).
*          B satisfies A = U * B * VT, so it has the same singular
*          values as A, and singular vectors related by U and VT.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if ZBDSQR did not converge, INFO specifies how many
*                superdiagonals of an intermediate bidiagonal form B
*                did not converge to zero. See the description of RWORK
*                above for details.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D0, 0.0D0 ),
     $                   CONE = ( 1.0D0, 0.0D0 ) )
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, WNTUA, WNTUAS, WNTUN, WNTUO, WNTUS,
     $                   WNTVA, WNTVAS, WNTVN, WNTVO, WNTVS
      INTEGER            BLK, CHUNK, I, IE, IERR, IR, IRWORK, ISCL,
     $                   ITAU, ITAUP, ITAUQ, IU, IWORK, LDWRKR, LDWRKU,
     $                   MAXWRK, MINMN, MINWRK, MNTHR, NCU, NCVT, NRU,
     $                   NRVT, WRKBL
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, SMLNUM
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   DUM( 1 )
      COMPLEX*16         CDUM( 1 )
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLASCL, XERBLA, ZBDSQR, ZGEBRD, ZGELQF, ZGEMM,
     $                   ZGEQRF, ZLACPY, ZLASCL, ZLASET, ZUNGBR, ZUNGLQ,
     $                   ZUNGQR, ZUNMBR
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, ZLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, ZLANGE
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      MINMN = MIN( M, N )
      MNTHR = ILAENV( 6, 'ZGESVD', JOBU // JOBVT, M, N, 0, 0 )
      WNTUA = LSAME( JOBU, 'A' )
      WNTUS = LSAME( JOBU, 'S' )
      WNTUAS = WNTUA .OR. WNTUS
      WNTUO = LSAME( JOBU, 'O' )
      WNTUN = LSAME( JOBU, 'N' )
      WNTVA = LSAME( JOBVT, 'A' )
      WNTVS = LSAME( JOBVT, 'S' )
      WNTVAS = WNTVA .OR. WNTVS
      WNTVO = LSAME( JOBVT, 'O' )
      WNTVN = LSAME( JOBVT, 'N' )
      MINWRK = 1
      LQUERY = ( LWORK.EQ.-1 )
*
      IF( .NOT.( WNTUA .OR. WNTUS .OR. WNTUO .OR. WNTUN ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( WNTVA .OR. WNTVS .OR. WNTVO .OR. WNTVN ) .OR.
     $         ( WNTVO .AND. WNTUO ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LDU.LT.1 .OR. ( WNTUAS .AND. LDU.LT.M ) ) THEN
         INFO = -9
      ELSE IF( LDVT.LT.1 .OR. ( WNTVA .AND. LDVT.LT.N ) .OR.
     $         ( WNTVS .AND. LDVT.LT.MINMN ) ) THEN
         INFO = -11
      END IF
*
*     Compute workspace
*      (Note: Comments in the code beginning "Workspace:" describe the
*       minimal amount of workspace needed at that point in the code,
*       as well as the preferred amount for good performance.
*       CWorkspace refers to complex workspace, and RWorkspace to
*       real workspace. NB refers to the optimal block size for the
*       immediately following subroutine, as returned by ILAENV.)
*
      IF( INFO.EQ.0 .AND. ( LWORK.GE.1 .OR. LQUERY ) .AND. M.GT.0 .AND.
     $    N.GT.0 ) THEN
         IF( M.GE.N ) THEN
*
*           Space needed for ZBDSQR is BDSPAC = 5*N
*
            IF( M.GE.MNTHR ) THEN
               IF( WNTUN ) THEN
*
*                 Path 1 (M much larger than N, JOBU='N')
*
                  MAXWRK = N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1,
     $                     -1 )
                  MAXWRK = MAX( MAXWRK, 2*N+2*N*
     $                     ILAENV( 1, 'ZGEBRD', ' ', N, N, -1, -1 ) )
                  IF( WNTVO .OR. WNTVAS )
     $               MAXWRK = MAX( MAXWRK, 2*N+( N-1 )*
     $                        ILAENV( 1, 'ZUNGBR', 'P', N, N, N, -1 ) )
                  MINWRK = 3*N
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTUO .AND. WNTVN ) THEN
*
*                 Path 2 (M much larger than N, JOBU='O', JOBVT='N')
*
                  WRKBL = N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, N+N*ILAENV( 1, 'ZUNGQR', ' ', M,
     $                    N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+2*N*
     $                    ILAENV( 1, 'ZGEBRD', ' ', N, N, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+N*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', N, N, N, -1 ) )
                  MAXWRK = MAX( N*N+WRKBL, N*N+M*N )
                  MINWRK = 2*N + M
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTUO .AND. WNTVAS ) THEN
*
*                 Path 3 (M much larger than N, JOBU='O', JOBVT='S' or
*                 'A')
*
                  WRKBL = N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, N+N*ILAENV( 1, 'ZUNGQR', ' ', M,
     $                    N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+2*N*
     $                    ILAENV( 1, 'ZGEBRD', ' ', N, N, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+N*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', N, N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+( N-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', N, N, N, -1 ) )
                  MAXWRK = MAX( N*N+WRKBL, N*N+M*N )
                  MINWRK = 2*N + M
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTUS .AND. WNTVN ) THEN
*
*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
*
                  WRKBL = N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, N+N*ILAENV( 1, 'ZUNGQR', ' ', M,
     $                    N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+2*N*
     $                    ILAENV( 1, 'ZGEBRD', ' ', N, N, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+N*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', N, N, N, -1 ) )
                  MAXWRK = N*N + WRKBL
                  MINWRK = 2*N + M
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTUS .AND. WNTVO ) THEN
*
*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
*
                  WRKBL = N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, N+N*ILAENV( 1, 'ZUNGQR', ' ', M,
     $                    N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+2*N*
     $                    ILAENV( 1, 'ZGEBRD', ' ', N, N, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+N*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', N, N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+( N-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', N, N, N, -1 ) )
                  MAXWRK = 2*N*N + WRKBL
                  MINWRK = 2*N + M
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTUS .AND. WNTVAS ) THEN
*
*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S' or
*                 'A')
*
                  WRKBL = N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, N+N*ILAENV( 1, 'ZUNGQR', ' ', M,
     $                    N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+2*N*
     $                    ILAENV( 1, 'ZGEBRD', ' ', N, N, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+N*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', N, N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+( N-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', N, N, N, -1 ) )
                  MAXWRK = N*N + WRKBL
                  MINWRK = 2*N + M
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTUA .AND. WNTVN ) THEN
*
*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
*
                  WRKBL = N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, N+M*ILAENV( 1, 'ZUNGQR', ' ', M,
     $                    M, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+2*N*
     $                    ILAENV( 1, 'ZGEBRD', ' ', N, N, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+N*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', N, N, N, -1 ) )
                  MAXWRK = N*N + WRKBL
                  MINWRK = 2*N + M
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTUA .AND. WNTVO ) THEN
*
*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
*
                  WRKBL = N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, N+M*ILAENV( 1, 'ZUNGQR', ' ', M,
     $                    M, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+2*N*
     $                    ILAENV( 1, 'ZGEBRD', ' ', N, N, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+N*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', N, N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+( N-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', N, N, N, -1 ) )
                  MAXWRK = 2*N*N + WRKBL
                  MINWRK = 2*N + M
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTUA .AND. WNTVAS ) THEN
*
*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S' or
*                 'A')
*
                  WRKBL = N + N*ILAENV( 1, 'ZGEQRF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, N+M*ILAENV( 1, 'ZUNGQR', ' ', M,
     $                    M, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+2*N*
     $                    ILAENV( 1, 'ZGEBRD', ' ', N, N, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+N*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', N, N, N, -1 ) )
                  WRKBL = MAX( WRKBL, 2*N+( N-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', N, N, N, -1 ) )
                  MAXWRK = N*N + WRKBL
                  MINWRK = 2*N + M
                  MAXWRK = MAX( MINWRK, MAXWRK )
               END IF
            ELSE
*
*              Path 10 (M at least N, but not much larger)
*
               MAXWRK = 2*N + ( M+N )*ILAENV( 1, 'ZGEBRD', ' ', M, N,
     $                  -1, -1 )
               IF( WNTUS .OR. WNTUO )
     $            MAXWRK = MAX( MAXWRK, 2*N+N*
     $                     ILAENV( 1, 'ZUNGBR', 'Q', M, N, N, -1 ) )
               IF( WNTUA )
     $            MAXWRK = MAX( MAXWRK, 2*N+M*
     $                     ILAENV( 1, 'ZUNGBR', 'Q', M, M, N, -1 ) )
               IF( .NOT.WNTVN )
     $            MAXWRK = MAX( MAXWRK, 2*N+( N-1 )*
     $                     ILAENV( 1, 'ZUNGBR', 'P', N, N, N, -1 ) )
               MINWRK = 2*N + M
               MAXWRK = MAX( MINWRK, MAXWRK )
            END IF
         ELSE
*
*           Space needed for ZBDSQR is BDSPAC = 5*M
*
            IF( N.GE.MNTHR ) THEN
               IF( WNTVN ) THEN
*
*                 Path 1t(N much larger than M, JOBVT='N')
*
                  MAXWRK = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1,
     $                     -1 )
                  MAXWRK = MAX( MAXWRK, 2*M+2*M*
     $                     ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )
                  IF( WNTUO .OR. WNTUAS )
     $               MAXWRK = MAX( MAXWRK, 2*M+M*
     $                        ILAENV( 1, 'ZUNGBR', 'Q', M, M, M, -1 ) )
                  MINWRK = 3*M
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTVO .AND. WNTUN ) THEN
*
*                 Path 2t(N much larger than M, JOBU='N', JOBVT='O')
*
                  WRKBL = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, M+M*ILAENV( 1, 'ZUNGLQ', ' ', M,
     $                    N, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+2*M*
     $                    ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+( M-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', M, M, M, -1 ) )
                  MAXWRK = MAX( M*M+WRKBL, M*M+M*N )
                  MINWRK = 2*M + N
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTVO .AND. WNTUAS ) THEN
*
*                 Path 3t(N much larger than M, JOBU='S' or 'A',
*                 JOBVT='O')
*
                  WRKBL = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, M+M*ILAENV( 1, 'ZUNGLQ', ' ', M,
     $                    N, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+2*M*
     $                    ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+( M-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', M, M, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+M*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', M, M, M, -1 ) )
                  MAXWRK = MAX( M*M+WRKBL, M*M+M*N )
                  MINWRK = 2*M + N
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTVS .AND. WNTUN ) THEN
*
*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
*
                  WRKBL = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, M+M*ILAENV( 1, 'ZUNGLQ', ' ', M,
     $                    N, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+2*M*
     $                    ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+( M-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', M, M, M, -1 ) )
                  MAXWRK = M*M + WRKBL
                  MINWRK = 2*M + N
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTVS .AND. WNTUO ) THEN
*
*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
*
                  WRKBL = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, M+M*ILAENV( 1, 'ZUNGLQ', ' ', M,
     $                    N, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+2*M*
     $                    ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+( M-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', M, M, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+M*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', M, M, M, -1 ) )
                  MAXWRK = 2*M*M + WRKBL
                  MINWRK = 2*M + N
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTVS .AND. WNTUAS ) THEN
*
*                 Path 6t(N much larger than M, JOBU='S' or 'A',
*                 JOBVT='S')
*
                  WRKBL = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, M+M*ILAENV( 1, 'ZUNGLQ', ' ', M,
     $                    N, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+2*M*
     $                    ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+( M-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', M, M, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+M*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', M, M, M, -1 ) )
                  MAXWRK = M*M + WRKBL
                  MINWRK = 2*M + N
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTVA .AND. WNTUN ) THEN
*
*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
*
                  WRKBL = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, M+N*ILAENV( 1, 'ZUNGLQ', ' ', N,
     $                    N, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+2*M*
     $                    ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+( M-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', M, M, M, -1 ) )
                  MAXWRK = M*M + WRKBL
                  MINWRK = 2*M + N
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTVA .AND. WNTUO ) THEN
*
*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
*
                  WRKBL = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, M+N*ILAENV( 1, 'ZUNGLQ', ' ', N,
     $                    N, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+2*M*
     $                    ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+( M-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', M, M, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+M*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', M, M, M, -1 ) )
                  MAXWRK = 2*M*M + WRKBL
                  MINWRK = 2*M + N
                  MAXWRK = MAX( MINWRK, MAXWRK )
               ELSE IF( WNTVA .AND. WNTUAS ) THEN
*
*                 Path 9t(N much larger than M, JOBU='S' or 'A',
*                 JOBVT='A')
*
                  WRKBL = M + M*ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
                  WRKBL = MAX( WRKBL, M+N*ILAENV( 1, 'ZUNGLQ', ' ', N,
     $                    N, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+2*M*
     $                    ILAENV( 1, 'ZGEBRD', ' ', M, M, -1, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+( M-1 )*
     $                    ILAENV( 1, 'ZUNGBR', 'P', M, M, M, -1 ) )
                  WRKBL = MAX( WRKBL, 2*M+M*
     $                    ILAENV( 1, 'ZUNGBR', 'Q', M, M, M, -1 ) )
                  MAXWRK = M*M + WRKBL
                  MINWRK = 2*M + N
                  MAXWRK = MAX( MINWRK, MAXWRK )
               END IF
            ELSE
*
*              Path 10t(N greater than M, but not much larger)
*
               MAXWRK = 2*M + ( M+N )*ILAENV( 1, 'ZGEBRD', ' ', M, N,
     $                  -1, -1 )
               IF( WNTVS .OR. WNTVO )
     $            MAXWRK = MAX( MAXWRK, 2*M+M*
     $                     ILAENV( 1, 'ZUNGBR', 'P', M, N, M, -1 ) )
               IF( WNTVA )
     $            MAXWRK = MAX( MAXWRK, 2*M+N*
     $                     ILAENV( 1, 'ZUNGBR', 'P', N, N, M, -1 ) )
               IF( .NOT.WNTUN )
     $            MAXWRK = MAX( MAXWRK, 2*M+( M-1 )*
     $                     ILAENV( 1, 'ZUNGBR', 'Q', M, M, M, -1 ) )
               MINWRK = 2*M + N
               MAXWRK = MAX( MINWRK, MAXWRK )
            END IF
         END IF
         WORK( 1 ) = MAXWRK
      END IF
*
      IF( LWORK.LT.MINWRK .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGESVD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         IF( LWORK.GE.1 )
     $      WORK( 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants
*
      EPS = DLAMCH( 'P' )
      SMLNUM = SQRT( DLAMCH( 'S' ) ) / EPS
      BIGNUM = ONE / SMLNUM
*
*     Scale A if max element outside range [SMLNUM,BIGNUM]
*
      ANRM = ZLANGE( 'M', M, N, A, LDA, DUM )
      ISCL = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.SMLNUM ) THEN
         ISCL = 1
         CALL ZLASCL( 'G', 0, 0, ANRM, SMLNUM, M, N, A, LDA, IERR )
      ELSE IF( ANRM.GT.BIGNUM ) THEN
         ISCL = 1
         CALL ZLASCL( 'G', 0, 0, ANRM, BIGNUM, M, N, A, LDA, IERR )
      END IF
*
      IF( M.GE.N ) THEN
*
*        A has at least as many rows as columns. If A has sufficiently
*        more rows than columns, first reduce using the QR
*        decomposition (if sufficient workspace available)
*
         IF( M.GE.MNTHR ) THEN
*
            IF( WNTUN ) THEN
*
*              Path 1 (M much larger than N, JOBU='N')
*              No left singular vectors to be computed
*
               ITAU = 1
               IWORK = ITAU + N
*
*              Compute A=Q*R
*              (CWorkspace: need 2*N, prefer N+N*NB)
*              (RWorkspace: need 0)
*
               CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ),
     $                      LWORK-IWORK+1, IERR )
*
*              Zero out below R
*
               CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO, A( 2, 1 ),
     $                      LDA )
               IE = 1
               ITAUQ = 1
               ITAUP = ITAUQ + N
               IWORK = ITAUP + N
*
*              Bidiagonalize R in A
*              (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
*              (RWorkspace: need N)
*
               CALL ZGEBRD( N, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ),
     $                      WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1,
     $                      IERR )
               NCVT = 0
               IF( WNTVO .OR. WNTVAS ) THEN
*
*                 If right singular vectors desired, generate P'.
*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  NCVT = N
               END IF
               IRWORK = IE + N
*
*              Perform bidiagonal QR iteration, computing right
*              singular vectors of A in A if desired
*              (CWorkspace: 0)
*              (RWorkspace: need BDSPAC)
*
               CALL ZBDSQR( 'U', N, NCVT, 0, 0, S, RWORK( IE ), A, LDA,
     $                      CDUM, 1, CDUM, 1, RWORK( IRWORK ), INFO )
*
*              If right singular vectors desired in VT, copy them there
*
               IF( WNTVAS )
     $            CALL ZLACPY( 'F', N, N, A, LDA, VT, LDVT )
*
            ELSE IF( WNTUO .AND. WNTVN ) THEN
*
*              Path 2 (M much larger than N, JOBU='O', JOBVT='N')
*              N left singular vectors to be overwritten on A and
*              no right singular vectors to be computed
*
               IF( LWORK.GE.N*N+3*N ) THEN
*
*                 Sufficient workspace for a fast algorithm
*
                  IR = 1
                  IF( LWORK.GE.MAX( WRKBL, LDA*N )+LDA*N ) THEN
*
*                    WORK(IU) is LDA by N, WORK(IR) is LDA by N
*
                     LDWRKU = LDA
                     LDWRKR = LDA
                  ELSE IF( LWORK.GE.MAX( WRKBL, LDA*N )+N*N ) THEN
*
*                    WORK(IU) is LDA by N, WORK(IR) is N by N
*
                     LDWRKU = LDA
                     LDWRKR = N
                  ELSE
*
*                    WORK(IU) is LDWRKU by N, WORK(IR) is N by N
*
                     LDWRKU = ( LWORK-N*N ) / N
                     LDWRKR = N
                  END IF
                  ITAU = IR + LDWRKR*N
                  IWORK = ITAU + N
*
*                 Compute A=Q*R
*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                 (RWorkspace: 0)
*
                  CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Copy R to WORK(IR) and zero out below it
*
                  CALL ZLACPY( 'U', N, N, A, LDA, WORK( IR ), LDWRKR )
                  CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                         WORK( IR+1 ), LDWRKR )
*
*                 Generate Q in A
*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGQR( M, N, N, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = 1
                  ITAUQ = ITAU
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N
*
*                 Bidiagonalize R in WORK(IR)
*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
*                 (RWorkspace: need N)
*
                  CALL ZGEBRD( N, N, WORK( IR ), LDWRKR, S, RWORK( IE ),
     $                         WORK( ITAUQ ), WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Generate left vectors bidiagonalizing R
*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
*                 (RWorkspace: need 0)
*
                  CALL ZUNGBR( 'Q', N, N, N, WORK( IR ), LDWRKR,
     $                         WORK( ITAUQ ), WORK( IWORK ),
     $                         LWORK-IWORK+1, IERR )
                  IRWORK = IE + N
*
*                 Perform bidiagonal QR iteration, computing left
*                 singular vectors of R in WORK(IR)
*                 (CWorkspace: need N*N)
*                 (RWorkspace: need BDSPAC)
*
                  CALL ZBDSQR( 'U', N, 0, N, 0, S, RWORK( IE ), CDUM, 1,
     $                         WORK( IR ), LDWRKR, CDUM, 1,
     $                         RWORK( IRWORK ), INFO )
                  IU = ITAUQ
*
*                 Multiply Q in A by left singular vectors of R in
*                 WORK(IR), storing result in WORK(IU) and copying to A
*                 (CWorkspace: need N*N+N, prefer N*N+M*N)
*                 (RWorkspace: 0)
*
                  DO 10 I = 1, M, LDWRKU
                     CHUNK = MIN( M-I+1, LDWRKU )
                     CALL ZGEMM( 'N', 'N', CHUNK, N, N, CONE, A( I, 1 ),
     $                           LDA, WORK( IR ), LDWRKR, CZERO,
     $                           WORK( IU ), LDWRKU )
                     CALL ZLACPY( 'F', CHUNK, N, WORK( IU ), LDWRKU,
     $                            A( I, 1 ), LDA )
   10             CONTINUE
*
               ELSE
*
*                 Insufficient workspace for a fast algorithm
*
                  IE = 1
                  ITAUQ = 1
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N
*
*                 Bidiagonalize A
*                 (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
*                 (RWorkspace: N)
*
                  CALL ZGEBRD( M, N, A, LDA, S, RWORK( IE ),
     $                         WORK( ITAUQ ), WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Generate left vectors bidiagonalizing A
*                 (CWorkspace: need 3*N, prefer 2*N+N*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'Q', M, N, N, A, LDA, WORK( ITAUQ ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IRWORK = IE + N
*
*                 Perform bidiagonal QR iteration, computing left
*                 singular vectors of A in A
*                 (CWorkspace: need 0)
*                 (RWorkspace: need BDSPAC)
*
                  CALL ZBDSQR( 'U', N, 0, M, 0, S, RWORK( IE ), CDUM, 1,
     $                         A, LDA, CDUM, 1, RWORK( IRWORK ), INFO )
*
               END IF
*
            ELSE IF( WNTUO .AND. WNTVAS ) THEN
*
*              Path 3 (M much larger than N, JOBU='O', JOBVT='S' or 'A')
*              N left singular vectors to be overwritten on A and
*              N right singular vectors to be computed in VT
*
               IF( LWORK.GE.N*N+3*N ) THEN
*
*                 Sufficient workspace for a fast algorithm
*
                  IR = 1
                  IF( LWORK.GE.MAX( WRKBL, LDA*N )+LDA*N ) THEN
*
*                    WORK(IU) is LDA by N and WORK(IR) is LDA by N
*
                     LDWRKU = LDA
                     LDWRKR = LDA
                  ELSE IF( LWORK.GE.MAX( WRKBL, LDA*N )+N*N ) THEN
*
*                    WORK(IU) is LDA by N and WORK(IR) is N by N
*
                     LDWRKU = LDA
                     LDWRKR = N
                  ELSE
*
*                    WORK(IU) is LDWRKU by N and WORK(IR) is N by N
*
                     LDWRKU = ( LWORK-N*N ) / N
                     LDWRKR = N
                  END IF
                  ITAU = IR + LDWRKR*N
                  IWORK = ITAU + N
*
*                 Compute A=Q*R
*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                 (RWorkspace: 0)
*
                  CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Copy R to VT, zeroing out below it
*
                  CALL ZLACPY( 'U', N, N, A, LDA, VT, LDVT )
                  CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO, VT( 2, 1 ),
     $                         LDVT )
*
*                 Generate Q in A
*                 (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGQR( M, N, N, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = 1
                  ITAUQ = ITAU
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N
*
*                 Bidiagonalize R in VT, copying result to WORK(IR)
*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
*                 (RWorkspace: need N)
*
                  CALL ZGEBRD( N, N, VT, LDVT, S, RWORK( IE ),
     $                         WORK( ITAUQ ), WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  CALL ZLACPY( 'L', N, N, VT, LDVT, WORK( IR ), LDWRKR )
*
*                 Generate left vectors bidiagonalizing R in WORK(IR)
*                 (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'Q', N, N, N, WORK( IR ), LDWRKR,
     $                         WORK( ITAUQ ), WORK( IWORK ),
     $                         LWORK-IWORK+1, IERR )
*
*                 Generate right vectors bidiagonalizing R in VT
*                 (CWorkspace: need N*N+3*N-1, prefer N*N+2*N+(N-1)*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IRWORK = IE + N
*
*                 Perform bidiagonal QR iteration, computing left
*                 singular vectors of R in WORK(IR) and computing right
*                 singular vectors of R in VT
*                 (CWorkspace: need N*N)
*                 (RWorkspace: need BDSPAC)
*
                  CALL ZBDSQR( 'U', N, N, N, 0, S, RWORK( IE ), VT,
     $                         LDVT, WORK( IR ), LDWRKR, CDUM, 1,
     $                         RWORK( IRWORK ), INFO )
                  IU = ITAUQ
*
*                 Multiply Q in A by left singular vectors of R in
*                 WORK(IR), storing result in WORK(IU) and copying to A
*                 (CWorkspace: need N*N+N, prefer N*N+M*N)
*                 (RWorkspace: 0)
*
                  DO 20 I = 1, M, LDWRKU
                     CHUNK = MIN( M-I+1, LDWRKU )
                     CALL ZGEMM( 'N', 'N', CHUNK, N, N, CONE, A( I, 1 ),
     $                           LDA, WORK( IR ), LDWRKR, CZERO,
     $                           WORK( IU ), LDWRKU )
                     CALL ZLACPY( 'F', CHUNK, N, WORK( IU ), LDWRKU,
     $                            A( I, 1 ), LDA )
   20             CONTINUE
*
               ELSE
*
*                 Insufficient workspace for a fast algorithm
*
                  ITAU = 1
                  IWORK = ITAU + N
*
*                 Compute A=Q*R
*                 (CWorkspace: need 2*N, prefer N+N*NB)
*                 (RWorkspace: 0)
*
                  CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Copy R to VT, zeroing out below it
*
                  CALL ZLACPY( 'U', N, N, A, LDA, VT, LDVT )
                  CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO, VT( 2, 1 ),
     $                         LDVT )
*
*                 Generate Q in A
*                 (CWorkspace: need 2*N, prefer N+N*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGQR( M, N, N, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = 1
                  ITAUQ = ITAU
                  ITAUP = ITAUQ + N
                  IWORK = ITAUP + N
*
*                 Bidiagonalize R in VT
*                 (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
*                 (RWorkspace: N)
*
                  CALL ZGEBRD( N, N, VT, LDVT, S, RWORK( IE ),
     $                         WORK( ITAUQ ), WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Multiply Q in A by left vectors bidiagonalizing R
*                 (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNMBR( 'Q', 'R', 'N', M, N, N, VT, LDVT,
     $                         WORK( ITAUQ ), A, LDA, WORK( IWORK ),
     $                         LWORK-IWORK+1, IERR )
*
*                 Generate right vectors bidiagonalizing R in VT
*                 (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IRWORK = IE + N
*
*                 Perform bidiagonal QR iteration, computing left
*                 singular vectors of A in A and computing right
*                 singular vectors of A in VT
*                 (CWorkspace: 0)
*                 (RWorkspace: need BDSPAC)
*
                  CALL ZBDSQR( 'U', N, N, M, 0, S, RWORK( IE ), VT,
     $                         LDVT, A, LDA, CDUM, 1, RWORK( IRWORK ),
     $                         INFO )
*
               END IF
*
            ELSE IF( WNTUS ) THEN
*
               IF( WNTVN ) THEN
*
*                 Path 4 (M much larger than N, JOBU='S', JOBVT='N')
*                 N left singular vectors to be computed in U and
*                 no right singular vectors to be computed
*
                  IF( LWORK.GE.N*N+3*N ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IR = 1
                     IF( LWORK.GE.WRKBL+LDA*N ) THEN
*
*                       WORK(IR) is LDA by N
*
                        LDWRKR = LDA
                     ELSE
*
*                       WORK(IR) is N by N
*
                        LDWRKR = N
                     END IF
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N
*
*                    Compute A=Q*R
*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy R to WORK(IR), zeroing out below it
*
                     CALL ZLACPY( 'U', N, N, A, LDA, WORK( IR ),
     $                            LDWRKR )
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            WORK( IR+1 ), LDWRKR )
*
*                    Generate Q in A
*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, N, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Bidiagonalize R in WORK(IR)
*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, WORK( IR ), LDWRKR, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate left vectors bidiagonalizing R in WORK(IR)
*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', N, N, N, WORK( IR ), LDWRKR,
     $                            WORK( ITAUQ ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of R in WORK(IR)
*                    (CWorkspace: need N*N)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, 0, N, 0, S, RWORK( IE ), CDUM,
     $                            1, WORK( IR ), LDWRKR, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
*                    Multiply Q in A by left singular vectors of R in
*                    WORK(IR), storing result in U
*                    (CWorkspace: need N*N)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, N, CONE, A, LDA,
     $                           WORK( IR ), LDWRKR, CZERO, U, LDU )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + N
*
*                    Compute A=Q*R, copying result to U
*                    (CWorkspace: need 2*N, prefer N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
*
*                    Generate Q in U
*                    (CWorkspace: need 2*N, prefer N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, N, N, U, LDU, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Zero out below R in A
*
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            A( 2, 1 ), LDA )
*
*                    Bidiagonalize R in A
*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, A, LDA, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply Q in U by left vectors bidiagonalizing R
*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'Q', 'R', 'N', M, N, N, A, LDA,
     $                            WORK( ITAUQ ), U, LDU, WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in U
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, 0, M, 0, S, RWORK( IE ), CDUM,
     $                            1, U, LDU, CDUM, 1, RWORK( IRWORK ),
     $                            INFO )
*
                  END IF
*
               ELSE IF( WNTVO ) THEN
*
*                 Path 5 (M much larger than N, JOBU='S', JOBVT='O')
*                 N left singular vectors to be computed in U and
*                 N right singular vectors to be overwritten on A
*
                  IF( LWORK.GE.2*N*N+3*N ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IU = 1
                     IF( LWORK.GE.WRKBL+2*LDA*N ) THEN
*
*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N
*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = LDA
                     ELSE IF( LWORK.GE.WRKBL+( LDA+N )*N ) THEN
*
*                       WORK(IU) is LDA by N and WORK(IR) is N by N
*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     ELSE
*
*                       WORK(IU) is N by N and WORK(IR) is N by N
*
                        LDWRKU = N
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     END IF
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N
*
*                    Compute A=Q*R
*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy R to WORK(IU), zeroing out below it
*
                     CALL ZLACPY( 'U', N, N, A, LDA, WORK( IU ),
     $                            LDWRKU )
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            WORK( IU+1 ), LDWRKU )
*
*                    Generate Q in A
*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, N, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Bidiagonalize R in WORK(IU), copying result to
*                    WORK(IR)
*                    (CWorkspace: need   2*N*N+3*N,
*                                 prefer 2*N*N+2*N+2*N*NB)
*                    (RWorkspace: need   N)
*
                     CALL ZGEBRD( N, N, WORK( IU ), LDWRKU, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', N, N, WORK( IU ), LDWRKU,
     $                            WORK( IR ), LDWRKR )
*
*                    Generate left bidiagonalizing vectors in WORK(IU)
*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', N, N, N, WORK( IU ), LDWRKU,
     $                            WORK( ITAUQ ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right bidiagonalizing vectors in WORK(IR)
*                    (CWorkspace: need   2*N*N+3*N-1,
*                                 prefer 2*N*N+2*N+(N-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', N, N, N, WORK( IR ), LDWRKR,
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of R in WORK(IU) and computing
*                    right singular vectors of R in WORK(IR)
*                    (CWorkspace: need 2*N*N)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, N, N, 0, S, RWORK( IE ),
     $                            WORK( IR ), LDWRKR, WORK( IU ),
     $                            LDWRKU, CDUM, 1, RWORK( IRWORK ),
     $                            INFO )
*
*                    Multiply Q in A by left singular vectors of R in
*                    WORK(IU), storing result in U
*                    (CWorkspace: need N*N)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, N, CONE, A, LDA,
     $                           WORK( IU ), LDWRKU, CZERO, U, LDU )
*
*                    Copy right singular vectors of R to A
*                    (CWorkspace: need N*N)
*                    (RWorkspace: 0)
*
                     CALL ZLACPY( 'F', N, N, WORK( IR ), LDWRKR, A,
     $                            LDA )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + N
*
*                    Compute A=Q*R, copying result to U
*                    (CWorkspace: need 2*N, prefer N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
*
*                    Generate Q in U
*                    (CWorkspace: need 2*N, prefer N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, N, N, U, LDU, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Zero out below R in A
*
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            A( 2, 1 ), LDA )
*
*                    Bidiagonalize R in A
*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, A, LDA, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply Q in U by left vectors bidiagonalizing R
*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'Q', 'R', 'N', M, N, N, A, LDA,
     $                            WORK( ITAUQ ), U, LDU, WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right vectors bidiagonalizing R in A
*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in U and computing right
*                    singular vectors of A in A
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, N, M, 0, S, RWORK( IE ), A,
     $                            LDA, U, LDU, CDUM, 1, RWORK( IRWORK ),
     $                            INFO )
*
                  END IF
*
               ELSE IF( WNTVAS ) THEN
*
*                 Path 6 (M much larger than N, JOBU='S', JOBVT='S'
*                         or 'A')
*                 N left singular vectors to be computed in U and
*                 N right singular vectors to be computed in VT
*
                  IF( LWORK.GE.N*N+3*N ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IU = 1
                     IF( LWORK.GE.WRKBL+LDA*N ) THEN
*
*                       WORK(IU) is LDA by N
*
                        LDWRKU = LDA
                     ELSE
*
*                       WORK(IU) is N by N
*
                        LDWRKU = N
                     END IF
                     ITAU = IU + LDWRKU*N
                     IWORK = ITAU + N
*
*                    Compute A=Q*R
*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy R to WORK(IU), zeroing out below it
*
                     CALL ZLACPY( 'U', N, N, A, LDA, WORK( IU ),
     $                            LDWRKU )
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            WORK( IU+1 ), LDWRKU )
*
*                    Generate Q in A
*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, N, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Bidiagonalize R in WORK(IU), copying result to VT
*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, WORK( IU ), LDWRKU, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', N, N, WORK( IU ), LDWRKU, VT,
     $                            LDVT )
*
*                    Generate left bidiagonalizing vectors in WORK(IU)
*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', N, N, N, WORK( IU ), LDWRKU,
     $                            WORK( ITAUQ ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right bidiagonalizing vectors in VT
*                    (CWorkspace: need   N*N+3*N-1,
*                                 prefer N*N+2*N+(N-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of R in WORK(IU) and computing
*                    right singular vectors of R in VT
*                    (CWorkspace: need N*N)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, N, N, 0, S, RWORK( IE ), VT,
     $                            LDVT, WORK( IU ), LDWRKU, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
*                    Multiply Q in A by left singular vectors of R in
*                    WORK(IU), storing result in U
*                    (CWorkspace: need N*N)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, N, CONE, A, LDA,
     $                           WORK( IU ), LDWRKU, CZERO, U, LDU )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + N
*
*                    Compute A=Q*R, copying result to U
*                    (CWorkspace: need 2*N, prefer N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
*
*                    Generate Q in U
*                    (CWorkspace: need 2*N, prefer N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, N, N, U, LDU, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy R to VT, zeroing out below it
*
                     CALL ZLACPY( 'U', N, N, A, LDA, VT, LDVT )
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            VT( 2, 1 ), LDVT )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Bidiagonalize R in VT
*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, VT, LDVT, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply Q in U by left bidiagonalizing vectors
*                    in VT
*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'Q', 'R', 'N', M, N, N, VT, LDVT,
     $                            WORK( ITAUQ ), U, LDU, WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right bidiagonalizing vectors in VT
*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in U and computing right
*                    singular vectors of A in VT
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, N, M, 0, S, RWORK( IE ), VT,
     $                            LDVT, U, LDU, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
                  END IF
*
               END IF
*
            ELSE IF( WNTUA ) THEN
*
               IF( WNTVN ) THEN
*
*                 Path 7 (M much larger than N, JOBU='A', JOBVT='N')
*                 M left singular vectors to be computed in U and
*                 no right singular vectors to be computed
*
                  IF( LWORK.GE.N*N+MAX( N+M, 3*N ) ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IR = 1
                     IF( LWORK.GE.WRKBL+LDA*N ) THEN
*
*                       WORK(IR) is LDA by N
*
                        LDWRKR = LDA
                     ELSE
*
*                       WORK(IR) is N by N
*
                        LDWRKR = N
                     END IF
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N
*
*                    Compute A=Q*R, copying result to U
*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
*
*                    Copy R to WORK(IR), zeroing out below it
*
                     CALL ZLACPY( 'U', N, N, A, LDA, WORK( IR ),
     $                            LDWRKR )
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            WORK( IR+1 ), LDWRKR )
*
*                    Generate Q in U
*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, M, N, U, LDU, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Bidiagonalize R in WORK(IR)
*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, WORK( IR ), LDWRKR, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate left bidiagonalizing vectors in WORK(IR)
*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', N, N, N, WORK( IR ), LDWRKR,
     $                            WORK( ITAUQ ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of R in WORK(IR)
*                    (CWorkspace: need N*N)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, 0, N, 0, S, RWORK( IE ), CDUM,
     $                            1, WORK( IR ), LDWRKR, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
*                    Multiply Q in U by left singular vectors of R in
*                    WORK(IR), storing result in A
*                    (CWorkspace: need N*N)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, N, CONE, U, LDU,
     $                           WORK( IR ), LDWRKR, CZERO, A, LDA )
*
*                    Copy left singular vectors of A from A to U
*
                     CALL ZLACPY( 'F', M, N, A, LDA, U, LDU )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + N
*
*                    Compute A=Q*R, copying result to U
*                    (CWorkspace: need 2*N, prefer N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
*
*                    Generate Q in U
*                    (CWorkspace: need N+M, prefer N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, M, N, U, LDU, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Zero out below R in A
*
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            A( 2, 1 ), LDA )
*
*                    Bidiagonalize R in A
*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, A, LDA, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply Q in U by left bidiagonalizing vectors
*                    in A
*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'Q', 'R', 'N', M, N, N, A, LDA,
     $                            WORK( ITAUQ ), U, LDU, WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in U
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, 0, M, 0, S, RWORK( IE ), CDUM,
     $                            1, U, LDU, CDUM, 1, RWORK( IRWORK ),
     $                            INFO )
*
                  END IF
*
               ELSE IF( WNTVO ) THEN
*
*                 Path 8 (M much larger than N, JOBU='A', JOBVT='O')
*                 M left singular vectors to be computed in U and
*                 N right singular vectors to be overwritten on A
*
                  IF( LWORK.GE.2*N*N+MAX( N+M, 3*N ) ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IU = 1
                     IF( LWORK.GE.WRKBL+2*LDA*N ) THEN
*
*                       WORK(IU) is LDA by N and WORK(IR) is LDA by N
*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = LDA
                     ELSE IF( LWORK.GE.WRKBL+( LDA+N )*N ) THEN
*
*                       WORK(IU) is LDA by N and WORK(IR) is N by N
*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     ELSE
*
*                       WORK(IU) is N by N and WORK(IR) is N by N
*
                        LDWRKU = N
                        IR = IU + LDWRKU*N
                        LDWRKR = N
                     END IF
                     ITAU = IR + LDWRKR*N
                     IWORK = ITAU + N
*
*                    Compute A=Q*R, copying result to U
*                    (CWorkspace: need 2*N*N+2*N, prefer 2*N*N+N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
*
*                    Generate Q in U
*                    (CWorkspace: need 2*N*N+N+M, prefer 2*N*N+N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, M, N, U, LDU, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy R to WORK(IU), zeroing out below it
*
                     CALL ZLACPY( 'U', N, N, A, LDA, WORK( IU ),
     $                            LDWRKU )
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            WORK( IU+1 ), LDWRKU )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Bidiagonalize R in WORK(IU), copying result to
*                    WORK(IR)
*                    (CWorkspace: need   2*N*N+3*N,
*                                 prefer 2*N*N+2*N+2*N*NB)
*                    (RWorkspace: need   N)
*
                     CALL ZGEBRD( N, N, WORK( IU ), LDWRKU, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', N, N, WORK( IU ), LDWRKU,
     $                            WORK( IR ), LDWRKR )
*
*                    Generate left bidiagonalizing vectors in WORK(IU)
*                    (CWorkspace: need 2*N*N+3*N, prefer 2*N*N+2*N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', N, N, N, WORK( IU ), LDWRKU,
     $                            WORK( ITAUQ ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right bidiagonalizing vectors in WORK(IR)
*                    (CWorkspace: need   2*N*N+3*N-1,
*                                 prefer 2*N*N+2*N+(N-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', N, N, N, WORK( IR ), LDWRKR,
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of R in WORK(IU) and computing
*                    right singular vectors of R in WORK(IR)
*                    (CWorkspace: need 2*N*N)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, N, N, 0, S, RWORK( IE ),
     $                            WORK( IR ), LDWRKR, WORK( IU ),
     $                            LDWRKU, CDUM, 1, RWORK( IRWORK ),
     $                            INFO )
*
*                    Multiply Q in U by left singular vectors of R in
*                    WORK(IU), storing result in A
*                    (CWorkspace: need N*N)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, N, CONE, U, LDU,
     $                           WORK( IU ), LDWRKU, CZERO, A, LDA )
*
*                    Copy left singular vectors of A from A to U
*
                     CALL ZLACPY( 'F', M, N, A, LDA, U, LDU )
*
*                    Copy right singular vectors of R from WORK(IR) to A
*
                     CALL ZLACPY( 'F', N, N, WORK( IR ), LDWRKR, A,
     $                            LDA )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + N
*
*                    Compute A=Q*R, copying result to U
*                    (CWorkspace: need 2*N, prefer N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
*
*                    Generate Q in U
*                    (CWorkspace: need N+M, prefer N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, M, N, U, LDU, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Zero out below R in A
*
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            A( 2, 1 ), LDA )
*
*                    Bidiagonalize R in A
*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, A, LDA, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply Q in U by left bidiagonalizing vectors
*                    in A
*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'Q', 'R', 'N', M, N, N, A, LDA,
     $                            WORK( ITAUQ ), U, LDU, WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right bidiagonalizing vectors in A
*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in U and computing right
*                    singular vectors of A in A
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, N, M, 0, S, RWORK( IE ), A,
     $                            LDA, U, LDU, CDUM, 1, RWORK( IRWORK ),
     $                            INFO )
*
                  END IF
*
               ELSE IF( WNTVAS ) THEN
*
*                 Path 9 (M much larger than N, JOBU='A', JOBVT='S'
*                         or 'A')
*                 M left singular vectors to be computed in U and
*                 N right singular vectors to be computed in VT
*
                  IF( LWORK.GE.N*N+MAX( N+M, 3*N ) ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IU = 1
                     IF( LWORK.GE.WRKBL+LDA*N ) THEN
*
*                       WORK(IU) is LDA by N
*
                        LDWRKU = LDA
                     ELSE
*
*                       WORK(IU) is N by N
*
                        LDWRKU = N
                     END IF
                     ITAU = IU + LDWRKU*N
                     IWORK = ITAU + N
*
*                    Compute A=Q*R, copying result to U
*                    (CWorkspace: need N*N+2*N, prefer N*N+N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
*
*                    Generate Q in U
*                    (CWorkspace: need N*N+N+M, prefer N*N+N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, M, N, U, LDU, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy R to WORK(IU), zeroing out below it
*
                     CALL ZLACPY( 'U', N, N, A, LDA, WORK( IU ),
     $                            LDWRKU )
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            WORK( IU+1 ), LDWRKU )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Bidiagonalize R in WORK(IU), copying result to VT
*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, WORK( IU ), LDWRKU, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', N, N, WORK( IU ), LDWRKU, VT,
     $                            LDVT )
*
*                    Generate left bidiagonalizing vectors in WORK(IU)
*                    (CWorkspace: need N*N+3*N, prefer N*N+2*N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', N, N, N, WORK( IU ), LDWRKU,
     $                            WORK( ITAUQ ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right bidiagonalizing vectors in VT
*                    (CWorkspace: need   N*N+3*N-1,
*                                 prefer N*N+2*N+(N-1)*NB)
*                    (RWorkspace: need   0)
*
                     CALL ZUNGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of R in WORK(IU) and computing
*                    right singular vectors of R in VT
*                    (CWorkspace: need N*N)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, N, N, 0, S, RWORK( IE ), VT,
     $                            LDVT, WORK( IU ), LDWRKU, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
*                    Multiply Q in U by left singular vectors of R in
*                    WORK(IU), storing result in A
*                    (CWorkspace: need N*N)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, N, CONE, U, LDU,
     $                           WORK( IU ), LDWRKU, CZERO, A, LDA )
*
*                    Copy left singular vectors of A from A to U
*
                     CALL ZLACPY( 'F', M, N, A, LDA, U, LDU )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + N
*
*                    Compute A=Q*R, copying result to U
*                    (CWorkspace: need 2*N, prefer N+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGEQRF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
*
*                    Generate Q in U
*                    (CWorkspace: need N+M, prefer N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGQR( M, M, N, U, LDU, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy R from A to VT, zeroing out below it
*
                     CALL ZLACPY( 'U', N, N, A, LDA, VT, LDVT )
                     CALL ZLASET( 'L', N-1, N-1, CZERO, CZERO,
     $                            VT( 2, 1 ), LDVT )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + N
                     IWORK = ITAUP + N
*
*                    Bidiagonalize R in VT
*                    (CWorkspace: need 3*N, prefer 2*N+2*N*NB)
*                    (RWorkspace: need N)
*
                     CALL ZGEBRD( N, N, VT, LDVT, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply Q in U by left bidiagonalizing vectors
*                    in VT
*                    (CWorkspace: need 2*N+M, prefer 2*N+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'Q', 'R', 'N', M, N, N, VT, LDVT,
     $                            WORK( ITAUQ ), U, LDU, WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right bidiagonalizing vectors in VT
*                    (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + N
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in U and computing right
*                    singular vectors of A in VT
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', N, N, M, 0, S, RWORK( IE ), VT,
     $                            LDVT, U, LDU, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
                  END IF
*
               END IF
*
            END IF
*
         ELSE
*
*           M .LT. MNTHR
*
*           Path 10 (M at least N, but not much larger)
*           Reduce to bidiagonal form without QR decomposition
*
            IE = 1
            ITAUQ = 1
            ITAUP = ITAUQ + N
            IWORK = ITAUP + N
*
*           Bidiagonalize A
*           (CWorkspace: need 2*N+M, prefer 2*N+(M+N)*NB)
*           (RWorkspace: need N)
*
            CALL ZGEBRD( M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ),
     $                   WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1,
     $                   IERR )
            IF( WNTUAS ) THEN
*
*              If left singular vectors desired in U, copy result to U
*              and generate left bidiagonalizing vectors in U
*              (CWorkspace: need 2*N+NCU, prefer 2*N+NCU*NB)
*              (RWorkspace: 0)
*
               CALL ZLACPY( 'L', M, N, A, LDA, U, LDU )
               IF( WNTUS )
     $            NCU = N
               IF( WNTUA )
     $            NCU = M
               CALL ZUNGBR( 'Q', M, NCU, N, U, LDU, WORK( ITAUQ ),
     $                      WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTVAS ) THEN
*
*              If right singular vectors desired in VT, copy result to
*              VT and generate right bidiagonalizing vectors in VT
*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*              (RWorkspace: 0)
*
               CALL ZLACPY( 'U', N, N, A, LDA, VT, LDVT )
               CALL ZUNGBR( 'P', N, N, N, VT, LDVT, WORK( ITAUP ),
     $                      WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTUO ) THEN
*
*              If left singular vectors desired in A, generate left
*              bidiagonalizing vectors in A
*              (CWorkspace: need 3*N, prefer 2*N+N*NB)
*              (RWorkspace: 0)
*
               CALL ZUNGBR( 'Q', M, N, N, A, LDA, WORK( ITAUQ ),
     $                      WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTVO ) THEN
*
*              If right singular vectors desired in A, generate right
*              bidiagonalizing vectors in A
*              (CWorkspace: need 3*N-1, prefer 2*N+(N-1)*NB)
*              (RWorkspace: 0)
*
               CALL ZUNGBR( 'P', N, N, N, A, LDA, WORK( ITAUP ),
     $                      WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IRWORK = IE + N
            IF( WNTUAS .OR. WNTUO )
     $         NRU = M
            IF( WNTUN )
     $         NRU = 0
            IF( WNTVAS .OR. WNTVO )
     $         NCVT = N
            IF( WNTVN )
     $         NCVT = 0
            IF( ( .NOT.WNTUO ) .AND. ( .NOT.WNTVO ) ) THEN
*
*              Perform bidiagonal QR iteration, if desired, computing
*              left singular vectors in U and computing right singular
*              vectors in VT
*              (CWorkspace: 0)
*              (RWorkspace: need BDSPAC)
*
               CALL ZBDSQR( 'U', N, NCVT, NRU, 0, S, RWORK( IE ), VT,
     $                      LDVT, U, LDU, CDUM, 1, RWORK( IRWORK ),
     $                      INFO )
            ELSE IF( ( .NOT.WNTUO ) .AND. WNTVO ) THEN
*
*              Perform bidiagonal QR iteration, if desired, computing
*              left singular vectors in U and computing right singular
*              vectors in A
*              (CWorkspace: 0)
*              (RWorkspace: need BDSPAC)
*
               CALL ZBDSQR( 'U', N, NCVT, NRU, 0, S, RWORK( IE ), A,
     $                      LDA, U, LDU, CDUM, 1, RWORK( IRWORK ),
     $                      INFO )
            ELSE
*
*              Perform bidiagonal QR iteration, if desired, computing
*              left singular vectors in A and computing right singular
*              vectors in VT
*              (CWorkspace: 0)
*              (RWorkspace: need BDSPAC)
*
               CALL ZBDSQR( 'U', N, NCVT, NRU, 0, S, RWORK( IE ), VT,
     $                      LDVT, A, LDA, CDUM, 1, RWORK( IRWORK ),
     $                      INFO )
            END IF
*
         END IF
*
      ELSE
*
*        A has more columns than rows. If A has sufficiently more
*        columns than rows, first reduce using the LQ decomposition (if
*        sufficient workspace available)
*
         IF( N.GE.MNTHR ) THEN
*
            IF( WNTVN ) THEN
*
*              Path 1t(N much larger than M, JOBVT='N')
*              No right singular vectors to be computed
*
               ITAU = 1
               IWORK = ITAU + M
*
*              Compute A=L*Q
*              (CWorkspace: need 2*M, prefer M+M*NB)
*              (RWorkspace: 0)
*
               CALL ZGELQF( M, N, A, LDA, WORK( ITAU ), WORK( IWORK ),
     $                      LWORK-IWORK+1, IERR )
*
*              Zero out above L
*
               CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO, A( 1, 2 ),
     $                      LDA )
               IE = 1
               ITAUQ = 1
               ITAUP = ITAUQ + M
               IWORK = ITAUP + M
*
*              Bidiagonalize L in A
*              (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
*              (RWorkspace: need M)
*
               CALL ZGEBRD( M, M, A, LDA, S, RWORK( IE ), WORK( ITAUQ ),
     $                      WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1,
     $                      IERR )
               IF( WNTUO .OR. WNTUAS ) THEN
*
*                 If left singular vectors desired, generate Q
*                 (CWorkspace: need 3*M, prefer 2*M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'Q', M, M, M, A, LDA, WORK( ITAUQ ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
               END IF
               IRWORK = IE + M
               NRU = 0
               IF( WNTUO .OR. WNTUAS )
     $            NRU = M
*
*              Perform bidiagonal QR iteration, computing left singular
*              vectors of A in A if desired
*              (CWorkspace: 0)
*              (RWorkspace: need BDSPAC)
*
               CALL ZBDSQR( 'U', M, 0, NRU, 0, S, RWORK( IE ), CDUM, 1,
     $                      A, LDA, CDUM, 1, RWORK( IRWORK ), INFO )
*
*              If left singular vectors desired in U, copy them there
*
               IF( WNTUAS )
     $            CALL ZLACPY( 'F', M, M, A, LDA, U, LDU )
*
            ELSE IF( WNTVO .AND. WNTUN ) THEN
*
*              Path 2t(N much larger than M, JOBU='N', JOBVT='O')
*              M right singular vectors to be overwritten on A and
*              no left singular vectors to be computed
*
               IF( LWORK.GE.M*M+3*M ) THEN
*
*                 Sufficient workspace for a fast algorithm
*
                  IR = 1
                  IF( LWORK.GE.MAX( WRKBL, LDA*N )+LDA*M ) THEN
*
*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
*
                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = LDA
                  ELSE IF( LWORK.GE.MAX( WRKBL, LDA*N )+M*M ) THEN
*
*                    WORK(IU) is LDA by N and WORK(IR) is M by M
*
                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = M
                  ELSE
*
*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M
*
                     LDWRKU = M
                     CHUNK = ( LWORK-M*M ) / M
                     LDWRKR = M
                  END IF
                  ITAU = IR + LDWRKR*M
                  IWORK = ITAU + M
*
*                 Compute A=L*Q
*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Copy L to WORK(IR) and zero out above it
*
                  CALL ZLACPY( 'L', M, M, A, LDA, WORK( IR ), LDWRKR )
                  CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                         WORK( IR+LDWRKR ), LDWRKR )
*
*                 Generate Q in A
*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGLQ( M, N, M, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = 1
                  ITAUQ = ITAU
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M
*
*                 Bidiagonalize L in WORK(IR)
*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
*                 (RWorkspace: need M)
*
                  CALL ZGEBRD( M, M, WORK( IR ), LDWRKR, S, RWORK( IE ),
     $                         WORK( ITAUQ ), WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Generate right vectors bidiagonalizing L
*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'P', M, M, M, WORK( IR ), LDWRKR,
     $                         WORK( ITAUP ), WORK( IWORK ),
     $                         LWORK-IWORK+1, IERR )
                  IRWORK = IE + M
*
*                 Perform bidiagonal QR iteration, computing right
*                 singular vectors of L in WORK(IR)
*                 (CWorkspace: need M*M)
*                 (RWorkspace: need BDSPAC)
*
                  CALL ZBDSQR( 'U', M, M, 0, 0, S, RWORK( IE ),
     $                         WORK( IR ), LDWRKR, CDUM, 1, CDUM, 1,
     $                         RWORK( IRWORK ), INFO )
                  IU = ITAUQ
*
*                 Multiply right singular vectors of L in WORK(IR) by Q
*                 in A, storing result in WORK(IU) and copying to A
*                 (CWorkspace: need M*M+M, prefer M*M+M*N)
*                 (RWorkspace: 0)
*
                  DO 30 I = 1, N, CHUNK
                     BLK = MIN( N-I+1, CHUNK )
                     CALL ZGEMM( 'N', 'N', M, BLK, M, CONE, WORK( IR ),
     $                           LDWRKR, A( 1, I ), LDA, CZERO,
     $                           WORK( IU ), LDWRKU )
                     CALL ZLACPY( 'F', M, BLK, WORK( IU ), LDWRKU,
     $                            A( 1, I ), LDA )
   30             CONTINUE
*
               ELSE
*
*                 Insufficient workspace for a fast algorithm
*
                  IE = 1
                  ITAUQ = 1
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M
*
*                 Bidiagonalize A
*                 (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
*                 (RWorkspace: need M)
*
                  CALL ZGEBRD( M, N, A, LDA, S, RWORK( IE ),
     $                         WORK( ITAUQ ), WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Generate right vectors bidiagonalizing A
*                 (CWorkspace: need 3*M, prefer 2*M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IRWORK = IE + M
*
*                 Perform bidiagonal QR iteration, computing right
*                 singular vectors of A in A
*                 (CWorkspace: 0)
*                 (RWorkspace: need BDSPAC)
*
                  CALL ZBDSQR( 'L', M, N, 0, 0, S, RWORK( IE ), A, LDA,
     $                         CDUM, 1, CDUM, 1, RWORK( IRWORK ), INFO )
*
               END IF
*
            ELSE IF( WNTVO .AND. WNTUAS ) THEN
*
*              Path 3t(N much larger than M, JOBU='S' or 'A', JOBVT='O')
*              M right singular vectors to be overwritten on A and
*              M left singular vectors to be computed in U
*
               IF( LWORK.GE.M*M+3*M ) THEN
*
*                 Sufficient workspace for a fast algorithm
*
                  IR = 1
                  IF( LWORK.GE.MAX( WRKBL, LDA*N )+LDA*M ) THEN
*
*                    WORK(IU) is LDA by N and WORK(IR) is LDA by M
*
                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = LDA
                  ELSE IF( LWORK.GE.MAX( WRKBL, LDA*N )+M*M ) THEN
*
*                    WORK(IU) is LDA by N and WORK(IR) is M by M
*
                     LDWRKU = LDA
                     CHUNK = N
                     LDWRKR = M
                  ELSE
*
*                    WORK(IU) is M by CHUNK and WORK(IR) is M by M
*
                     LDWRKU = M
                     CHUNK = ( LWORK-M*M ) / M
                     LDWRKR = M
                  END IF
                  ITAU = IR + LDWRKR*M
                  IWORK = ITAU + M
*
*                 Compute A=L*Q
*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Copy L to U, zeroing about above it
*
                  CALL ZLACPY( 'L', M, M, A, LDA, U, LDU )
                  CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO, U( 1, 2 ),
     $                         LDU )
*
*                 Generate Q in A
*                 (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGLQ( M, N, M, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = 1
                  ITAUQ = ITAU
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M
*
*                 Bidiagonalize L in U, copying result to WORK(IR)
*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
*                 (RWorkspace: need M)
*
                  CALL ZGEBRD( M, M, U, LDU, S, RWORK( IE ),
     $                         WORK( ITAUQ ), WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  CALL ZLACPY( 'U', M, M, U, LDU, WORK( IR ), LDWRKR )
*
*                 Generate right vectors bidiagonalizing L in WORK(IR)
*                 (CWorkspace: need M*M+3*M-1, prefer M*M+2*M+(M-1)*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'P', M, M, M, WORK( IR ), LDWRKR,
     $                         WORK( ITAUP ), WORK( IWORK ),
     $                         LWORK-IWORK+1, IERR )
*
*                 Generate left vectors bidiagonalizing L in U
*                 (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IRWORK = IE + M
*
*                 Perform bidiagonal QR iteration, computing left
*                 singular vectors of L in U, and computing right
*                 singular vectors of L in WORK(IR)
*                 (CWorkspace: need M*M)
*                 (RWorkspace: need BDSPAC)
*
                  CALL ZBDSQR( 'U', M, M, M, 0, S, RWORK( IE ),
     $                         WORK( IR ), LDWRKR, U, LDU, CDUM, 1,
     $                         RWORK( IRWORK ), INFO )
                  IU = ITAUQ
*
*                 Multiply right singular vectors of L in WORK(IR) by Q
*                 in A, storing result in WORK(IU) and copying to A
*                 (CWorkspace: need M*M+M, prefer M*M+M*N))
*                 (RWorkspace: 0)
*
                  DO 40 I = 1, N, CHUNK
                     BLK = MIN( N-I+1, CHUNK )
                     CALL ZGEMM( 'N', 'N', M, BLK, M, CONE, WORK( IR ),
     $                           LDWRKR, A( 1, I ), LDA, CZERO,
     $                           WORK( IU ), LDWRKU )
                     CALL ZLACPY( 'F', M, BLK, WORK( IU ), LDWRKU,
     $                            A( 1, I ), LDA )
   40             CONTINUE
*
               ELSE
*
*                 Insufficient workspace for a fast algorithm
*
                  ITAU = 1
                  IWORK = ITAU + M
*
*                 Compute A=L*Q
*                 (CWorkspace: need 2*M, prefer M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Copy L to U, zeroing out above it
*
                  CALL ZLACPY( 'L', M, M, A, LDA, U, LDU )
                  CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO, U( 1, 2 ),
     $                         LDU )
*
*                 Generate Q in A
*                 (CWorkspace: need 2*M, prefer M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGLQ( M, N, M, A, LDA, WORK( ITAU ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IE = 1
                  ITAUQ = ITAU
                  ITAUP = ITAUQ + M
                  IWORK = ITAUP + M
*
*                 Bidiagonalize L in U
*                 (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
*                 (RWorkspace: need M)
*
                  CALL ZGEBRD( M, M, U, LDU, S, RWORK( IE ),
     $                         WORK( ITAUQ ), WORK( ITAUP ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                 Multiply right vectors bidiagonalizing L by Q in A
*                 (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNMBR( 'P', 'L', 'C', M, N, M, U, LDU,
     $                         WORK( ITAUP ), A, LDA, WORK( IWORK ),
     $                         LWORK-IWORK+1, IERR )
*
*                 Generate left vectors bidiagonalizing L in U
*                 (CWorkspace: need 3*M, prefer 2*M+M*NB)
*                 (RWorkspace: 0)
*
                  CALL ZUNGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ),
     $                         WORK( IWORK ), LWORK-IWORK+1, IERR )
                  IRWORK = IE + M
*
*                 Perform bidiagonal QR iteration, computing left
*                 singular vectors of A in U and computing right
*                 singular vectors of A in A
*                 (CWorkspace: 0)
*                 (RWorkspace: need BDSPAC)
*
                  CALL ZBDSQR( 'U', M, N, M, 0, S, RWORK( IE ), A, LDA,
     $                         U, LDU, CDUM, 1, RWORK( IRWORK ), INFO )
*
               END IF
*
            ELSE IF( WNTVS ) THEN
*
               IF( WNTUN ) THEN
*
*                 Path 4t(N much larger than M, JOBU='N', JOBVT='S')
*                 M right singular vectors to be computed in VT and
*                 no left singular vectors to be computed
*
                  IF( LWORK.GE.M*M+3*M ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IR = 1
                     IF( LWORK.GE.WRKBL+LDA*M ) THEN
*
*                       WORK(IR) is LDA by M
*
                        LDWRKR = LDA
                     ELSE
*
*                       WORK(IR) is M by M
*
                        LDWRKR = M
                     END IF
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M
*
*                    Compute A=L*Q
*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy L to WORK(IR), zeroing out above it
*
                     CALL ZLACPY( 'L', M, M, A, LDA, WORK( IR ),
     $                            LDWRKR )
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            WORK( IR+LDWRKR ), LDWRKR )
*
*                    Generate Q in A
*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( M, N, M, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Bidiagonalize L in WORK(IR)
*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, WORK( IR ), LDWRKR, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right vectors bidiagonalizing L in
*                    WORK(IR)
*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', M, M, M, WORK( IR ), LDWRKR,
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing right
*                    singular vectors of L in WORK(IR)
*                    (CWorkspace: need M*M)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, M, 0, 0, S, RWORK( IE ),
     $                            WORK( IR ), LDWRKR, CDUM, 1, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
*                    Multiply right singular vectors of L in WORK(IR) by
*                    Q in A, storing result in VT
*                    (CWorkspace: need M*M)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, M, CONE, WORK( IR ),
     $                           LDWRKR, A, LDA, CZERO, VT, LDVT )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + M
*
*                    Compute A=L*Q
*                    (CWorkspace: need 2*M, prefer M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy result to VT
*
                     CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
*
*                    Generate Q in VT
*                    (CWorkspace: need 2*M, prefer M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( M, N, M, VT, LDVT, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Zero out above L in A
*
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            A( 1, 2 ), LDA )
*
*                    Bidiagonalize L in A
*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, A, LDA, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply right vectors bidiagonalizing L by Q in VT
*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'P', 'L', 'C', M, N, M, A, LDA,
     $                            WORK( ITAUP ), VT, LDVT,
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing right
*                    singular vectors of A in VT
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, N, 0, 0, S, RWORK( IE ), VT,
     $                            LDVT, CDUM, 1, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
                  END IF
*
               ELSE IF( WNTUO ) THEN
*
*                 Path 5t(N much larger than M, JOBU='O', JOBVT='S')
*                 M right singular vectors to be computed in VT and
*                 M left singular vectors to be overwritten on A
*
                  IF( LWORK.GE.2*M*M+3*M ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IU = 1
                     IF( LWORK.GE.WRKBL+2*LDA*M ) THEN
*
*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = LDA
                     ELSE IF( LWORK.GE.WRKBL+( LDA+M )*M ) THEN
*
*                       WORK(IU) is LDA by M and WORK(IR) is M by M
*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     ELSE
*
*                       WORK(IU) is M by M and WORK(IR) is M by M
*
                        LDWRKU = M
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     END IF
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M
*
*                    Compute A=L*Q
*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy L to WORK(IU), zeroing out below it
*
                     CALL ZLACPY( 'L', M, M, A, LDA, WORK( IU ),
     $                            LDWRKU )
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            WORK( IU+LDWRKU ), LDWRKU )
*
*                    Generate Q in A
*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( M, N, M, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Bidiagonalize L in WORK(IU), copying result to
*                    WORK(IR)
*                    (CWorkspace: need   2*M*M+3*M,
*                                 prefer 2*M*M+2*M+2*M*NB)
*                    (RWorkspace: need   M)
*
                     CALL ZGEBRD( M, M, WORK( IU ), LDWRKU, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, M, WORK( IU ), LDWRKU,
     $                            WORK( IR ), LDWRKR )
*
*                    Generate right bidiagonalizing vectors in WORK(IU)
*                    (CWorkspace: need   2*M*M+3*M-1,
*                                 prefer 2*M*M+2*M+(M-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', M, M, M, WORK( IU ), LDWRKU,
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate left bidiagonalizing vectors in WORK(IR)
*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', M, M, M, WORK( IR ), LDWRKR,
     $                            WORK( ITAUQ ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of L in WORK(IR) and computing
*                    right singular vectors of L in WORK(IU)
*                    (CWorkspace: need 2*M*M)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, M, M, 0, S, RWORK( IE ),
     $                            WORK( IU ), LDWRKU, WORK( IR ),
     $                            LDWRKR, CDUM, 1, RWORK( IRWORK ),
     $                            INFO )
*
*                    Multiply right singular vectors of L in WORK(IU) by
*                    Q in A, storing result in VT
*                    (CWorkspace: need M*M)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, M, CONE, WORK( IU ),
     $                           LDWRKU, A, LDA, CZERO, VT, LDVT )
*
*                    Copy left singular vectors of L to A
*                    (CWorkspace: need M*M)
*                    (RWorkspace: 0)
*
                     CALL ZLACPY( 'F', M, M, WORK( IR ), LDWRKR, A,
     $                            LDA )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + M
*
*                    Compute A=L*Q, copying result to VT
*                    (CWorkspace: need 2*M, prefer M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
*
*                    Generate Q in VT
*                    (CWorkspace: need 2*M, prefer M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( M, N, M, VT, LDVT, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Zero out above L in A
*
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            A( 1, 2 ), LDA )
*
*                    Bidiagonalize L in A
*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, A, LDA, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply right vectors bidiagonalizing L by Q in VT
*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'P', 'L', 'C', M, N, M, A, LDA,
     $                            WORK( ITAUP ), VT, LDVT,
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Generate left bidiagonalizing vectors of L in A
*                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', M, M, M, A, LDA, WORK( ITAUQ ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in A and computing right
*                    singular vectors of A in VT
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, N, M, 0, S, RWORK( IE ), VT,
     $                            LDVT, A, LDA, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
                  END IF
*
               ELSE IF( WNTUAS ) THEN
*
*                 Path 6t(N much larger than M, JOBU='S' or 'A',
*                         JOBVT='S')
*                 M right singular vectors to be computed in VT and
*                 M left singular vectors to be computed in U
*
                  IF( LWORK.GE.M*M+3*M ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IU = 1
                     IF( LWORK.GE.WRKBL+LDA*M ) THEN
*
*                       WORK(IU) is LDA by N
*
                        LDWRKU = LDA
                     ELSE
*
*                       WORK(IU) is LDA by M
*
                        LDWRKU = M
                     END IF
                     ITAU = IU + LDWRKU*M
                     IWORK = ITAU + M
*
*                    Compute A=L*Q
*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy L to WORK(IU), zeroing out above it
*
                     CALL ZLACPY( 'L', M, M, A, LDA, WORK( IU ),
     $                            LDWRKU )
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            WORK( IU+LDWRKU ), LDWRKU )
*
*                    Generate Q in A
*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( M, N, M, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Bidiagonalize L in WORK(IU), copying result to U
*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, WORK( IU ), LDWRKU, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, M, WORK( IU ), LDWRKU, U,
     $                            LDU )
*
*                    Generate right bidiagonalizing vectors in WORK(IU)
*                    (CWorkspace: need   M*M+3*M-1,
*                                 prefer M*M+2*M+(M-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', M, M, M, WORK( IU ), LDWRKU,
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate left bidiagonalizing vectors in U
*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of L in U and computing right
*                    singular vectors of L in WORK(IU)
*                    (CWorkspace: need M*M)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, M, M, 0, S, RWORK( IE ),
     $                            WORK( IU ), LDWRKU, U, LDU, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
*                    Multiply right singular vectors of L in WORK(IU) by
*                    Q in A, storing result in VT
*                    (CWorkspace: need M*M)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, M, CONE, WORK( IU ),
     $                           LDWRKU, A, LDA, CZERO, VT, LDVT )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + M
*
*                    Compute A=L*Q, copying result to VT
*                    (CWorkspace: need 2*M, prefer M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
*
*                    Generate Q in VT
*                    (CWorkspace: need 2*M, prefer M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( M, N, M, VT, LDVT, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy L to U, zeroing out above it
*
                     CALL ZLACPY( 'L', M, M, A, LDA, U, LDU )
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            U( 1, 2 ), LDU )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Bidiagonalize L in U
*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, U, LDU, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply right bidiagonalizing vectors in U by Q
*                    in VT
*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'P', 'L', 'C', M, N, M, U, LDU,
     $                            WORK( ITAUP ), VT, LDVT,
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Generate left bidiagonalizing vectors in U
*                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in U and computing right
*                    singular vectors of A in VT
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, N, M, 0, S, RWORK( IE ), VT,
     $                            LDVT, U, LDU, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
                  END IF
*
               END IF
*
            ELSE IF( WNTVA ) THEN
*
               IF( WNTUN ) THEN
*
*                 Path 7t(N much larger than M, JOBU='N', JOBVT='A')
*                 N right singular vectors to be computed in VT and
*                 no left singular vectors to be computed
*
                  IF( LWORK.GE.M*M+MAX( N+M, 3*M ) ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IR = 1
                     IF( LWORK.GE.WRKBL+LDA*M ) THEN
*
*                       WORK(IR) is LDA by M
*
                        LDWRKR = LDA
                     ELSE
*
*                       WORK(IR) is M by M
*
                        LDWRKR = M
                     END IF
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M
*
*                    Compute A=L*Q, copying result to VT
*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
*
*                    Copy L to WORK(IR), zeroing out above it
*
                     CALL ZLACPY( 'L', M, M, A, LDA, WORK( IR ),
     $                            LDWRKR )
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            WORK( IR+LDWRKR ), LDWRKR )
*
*                    Generate Q in VT
*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( N, N, M, VT, LDVT, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Bidiagonalize L in WORK(IR)
*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, WORK( IR ), LDWRKR, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate right bidiagonalizing vectors in WORK(IR)
*                    (CWorkspace: need   M*M+3*M-1,
*                                 prefer M*M+2*M+(M-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', M, M, M, WORK( IR ), LDWRKR,
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing right
*                    singular vectors of L in WORK(IR)
*                    (CWorkspace: need M*M)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, M, 0, 0, S, RWORK( IE ),
     $                            WORK( IR ), LDWRKR, CDUM, 1, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
*                    Multiply right singular vectors of L in WORK(IR) by
*                    Q in VT, storing result in A
*                    (CWorkspace: need M*M)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, M, CONE, WORK( IR ),
     $                           LDWRKR, VT, LDVT, CZERO, A, LDA )
*
*                    Copy right singular vectors of A from A to VT
*
                     CALL ZLACPY( 'F', M, N, A, LDA, VT, LDVT )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + M
*
*                    Compute A=L*Q, copying result to VT
*                    (CWorkspace: need 2*M, prefer M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
*
*                    Generate Q in VT
*                    (CWorkspace: need M+N, prefer M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( N, N, M, VT, LDVT, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Zero out above L in A
*
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            A( 1, 2 ), LDA )
*
*                    Bidiagonalize L in A
*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, A, LDA, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply right bidiagonalizing vectors in A by Q
*                    in VT
*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'P', 'L', 'C', M, N, M, A, LDA,
     $                            WORK( ITAUP ), VT, LDVT,
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing right
*                    singular vectors of A in VT
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, N, 0, 0, S, RWORK( IE ), VT,
     $                            LDVT, CDUM, 1, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
                  END IF
*
               ELSE IF( WNTUO ) THEN
*
*                 Path 8t(N much larger than M, JOBU='O', JOBVT='A')
*                 N right singular vectors to be computed in VT and
*                 M left singular vectors to be overwritten on A
*
                  IF( LWORK.GE.2*M*M+MAX( N+M, 3*M ) ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IU = 1
                     IF( LWORK.GE.WRKBL+2*LDA*M ) THEN
*
*                       WORK(IU) is LDA by M and WORK(IR) is LDA by M
*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = LDA
                     ELSE IF( LWORK.GE.WRKBL+( LDA+M )*M ) THEN
*
*                       WORK(IU) is LDA by M and WORK(IR) is M by M
*
                        LDWRKU = LDA
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     ELSE
*
*                       WORK(IU) is M by M and WORK(IR) is M by M
*
                        LDWRKU = M
                        IR = IU + LDWRKU*M
                        LDWRKR = M
                     END IF
                     ITAU = IR + LDWRKR*M
                     IWORK = ITAU + M
*
*                    Compute A=L*Q, copying result to VT
*                    (CWorkspace: need 2*M*M+2*M, prefer 2*M*M+M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
*
*                    Generate Q in VT
*                    (CWorkspace: need 2*M*M+M+N, prefer 2*M*M+M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( N, N, M, VT, LDVT, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy L to WORK(IU), zeroing out above it
*
                     CALL ZLACPY( 'L', M, M, A, LDA, WORK( IU ),
     $                            LDWRKU )
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            WORK( IU+LDWRKU ), LDWRKU )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Bidiagonalize L in WORK(IU), copying result to
*                    WORK(IR)
*                    (CWorkspace: need   2*M*M+3*M,
*                                 prefer 2*M*M+2*M+2*M*NB)
*                    (RWorkspace: need   M)
*
                     CALL ZGEBRD( M, M, WORK( IU ), LDWRKU, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, M, WORK( IU ), LDWRKU,
     $                            WORK( IR ), LDWRKR )
*
*                    Generate right bidiagonalizing vectors in WORK(IU)
*                    (CWorkspace: need   2*M*M+3*M-1,
*                                 prefer 2*M*M+2*M+(M-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', M, M, M, WORK( IU ), LDWRKU,
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate left bidiagonalizing vectors in WORK(IR)
*                    (CWorkspace: need 2*M*M+3*M, prefer 2*M*M+2*M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', M, M, M, WORK( IR ), LDWRKR,
     $                            WORK( ITAUQ ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of L in WORK(IR) and computing
*                    right singular vectors of L in WORK(IU)
*                    (CWorkspace: need 2*M*M)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, M, M, 0, S, RWORK( IE ),
     $                            WORK( IU ), LDWRKU, WORK( IR ),
     $                            LDWRKR, CDUM, 1, RWORK( IRWORK ),
     $                            INFO )
*
*                    Multiply right singular vectors of L in WORK(IU) by
*                    Q in VT, storing result in A
*                    (CWorkspace: need M*M)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, M, CONE, WORK( IU ),
     $                           LDWRKU, VT, LDVT, CZERO, A, LDA )
*
*                    Copy right singular vectors of A from A to VT
*
                     CALL ZLACPY( 'F', M, N, A, LDA, VT, LDVT )
*
*                    Copy left singular vectors of A from WORK(IR) to A
*
                     CALL ZLACPY( 'F', M, M, WORK( IR ), LDWRKR, A,
     $                            LDA )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + M
*
*                    Compute A=L*Q, copying result to VT
*                    (CWorkspace: need 2*M, prefer M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
*
*                    Generate Q in VT
*                    (CWorkspace: need M+N, prefer M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( N, N, M, VT, LDVT, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Zero out above L in A
*
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            A( 1, 2 ), LDA )
*
*                    Bidiagonalize L in A
*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, A, LDA, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply right bidiagonalizing vectors in A by Q
*                    in VT
*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'P', 'L', 'C', M, N, M, A, LDA,
     $                            WORK( ITAUP ), VT, LDVT,
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Generate left bidiagonalizing vectors in A
*                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', M, M, M, A, LDA, WORK( ITAUQ ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in A and computing right
*                    singular vectors of A in VT
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, N, M, 0, S, RWORK( IE ), VT,
     $                            LDVT, A, LDA, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
                  END IF
*
               ELSE IF( WNTUAS ) THEN
*
*                 Path 9t(N much larger than M, JOBU='S' or 'A',
*                         JOBVT='A')
*                 N right singular vectors to be computed in VT and
*                 M left singular vectors to be computed in U
*
                  IF( LWORK.GE.M*M+MAX( N+M, 3*M ) ) THEN
*
*                    Sufficient workspace for a fast algorithm
*
                     IU = 1
                     IF( LWORK.GE.WRKBL+LDA*M ) THEN
*
*                       WORK(IU) is LDA by M
*
                        LDWRKU = LDA
                     ELSE
*
*                       WORK(IU) is M by M
*
                        LDWRKU = M
                     END IF
                     ITAU = IU + LDWRKU*M
                     IWORK = ITAU + M
*
*                    Compute A=L*Q, copying result to VT
*                    (CWorkspace: need M*M+2*M, prefer M*M+M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
*
*                    Generate Q in VT
*                    (CWorkspace: need M*M+M+N, prefer M*M+M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( N, N, M, VT, LDVT, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy L to WORK(IU), zeroing out above it
*
                     CALL ZLACPY( 'L', M, M, A, LDA, WORK( IU ),
     $                            LDWRKU )
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            WORK( IU+LDWRKU ), LDWRKU )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Bidiagonalize L in WORK(IU), copying result to U
*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, WORK( IU ), LDWRKU, S,
     $                            RWORK( IE ), WORK( ITAUQ ),
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'L', M, M, WORK( IU ), LDWRKU, U,
     $                            LDU )
*
*                    Generate right bidiagonalizing vectors in WORK(IU)
*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+(M-1)*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'P', M, M, M, WORK( IU ), LDWRKU,
     $                            WORK( ITAUP ), WORK( IWORK ),
     $                            LWORK-IWORK+1, IERR )
*
*                    Generate left bidiagonalizing vectors in U
*                    (CWorkspace: need M*M+3*M, prefer M*M+2*M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of L in U and computing right
*                    singular vectors of L in WORK(IU)
*                    (CWorkspace: need M*M)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, M, M, 0, S, RWORK( IE ),
     $                            WORK( IU ), LDWRKU, U, LDU, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
*                    Multiply right singular vectors of L in WORK(IU) by
*                    Q in VT, storing result in A
*                    (CWorkspace: need M*M)
*                    (RWorkspace: 0)
*
                     CALL ZGEMM( 'N', 'N', M, N, M, CONE, WORK( IU ),
     $                           LDWRKU, VT, LDVT, CZERO, A, LDA )
*
*                    Copy right singular vectors of A from A to VT
*
                     CALL ZLACPY( 'F', M, N, A, LDA, VT, LDVT )
*
                  ELSE
*
*                    Insufficient workspace for a fast algorithm
*
                     ITAU = 1
                     IWORK = ITAU + M
*
*                    Compute A=L*Q, copying result to VT
*                    (CWorkspace: need 2*M, prefer M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZGELQF( M, N, A, LDA, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
*
*                    Generate Q in VT
*                    (CWorkspace: need M+N, prefer M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGLQ( N, N, M, VT, LDVT, WORK( ITAU ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Copy L to U, zeroing out above it
*
                     CALL ZLACPY( 'L', M, M, A, LDA, U, LDU )
                     CALL ZLASET( 'U', M-1, M-1, CZERO, CZERO,
     $                            U( 1, 2 ), LDU )
                     IE = 1
                     ITAUQ = ITAU
                     ITAUP = ITAUQ + M
                     IWORK = ITAUP + M
*
*                    Bidiagonalize L in U
*                    (CWorkspace: need 3*M, prefer 2*M+2*M*NB)
*                    (RWorkspace: need M)
*
                     CALL ZGEBRD( M, M, U, LDU, S, RWORK( IE ),
     $                            WORK( ITAUQ ), WORK( ITAUP ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Multiply right bidiagonalizing vectors in U by Q
*                    in VT
*                    (CWorkspace: need 2*M+N, prefer 2*M+N*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNMBR( 'P', 'L', 'C', M, N, M, U, LDU,
     $                            WORK( ITAUP ), VT, LDVT,
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
*
*                    Generate left bidiagonalizing vectors in U
*                    (CWorkspace: need 3*M, prefer 2*M+M*NB)
*                    (RWorkspace: 0)
*
                     CALL ZUNGBR( 'Q', M, M, M, U, LDU, WORK( ITAUQ ),
     $                            WORK( IWORK ), LWORK-IWORK+1, IERR )
                     IRWORK = IE + M
*
*                    Perform bidiagonal QR iteration, computing left
*                    singular vectors of A in U and computing right
*                    singular vectors of A in VT
*                    (CWorkspace: 0)
*                    (RWorkspace: need BDSPAC)
*
                     CALL ZBDSQR( 'U', M, N, M, 0, S, RWORK( IE ), VT,
     $                            LDVT, U, LDU, CDUM, 1,
     $                            RWORK( IRWORK ), INFO )
*
                  END IF
*
               END IF
*
            END IF
*
         ELSE
*
*           N .LT. MNTHR
*
*           Path 10t(N greater than M, but not much larger)
*           Reduce to bidiagonal form without LQ decomposition
*
            IE = 1
            ITAUQ = 1
            ITAUP = ITAUQ + M
            IWORK = ITAUP + M
*
*           Bidiagonalize A
*           (CWorkspace: need 2*M+N, prefer 2*M+(M+N)*NB)
*           (RWorkspace: M)
*
            CALL ZGEBRD( M, N, A, LDA, S, RWORK( IE ), WORK( ITAUQ ),
     $                   WORK( ITAUP ), WORK( IWORK ), LWORK-IWORK+1,
     $                   IERR )
            IF( WNTUAS ) THEN
*
*              If left singular vectors desired in U, copy result to U
*              and generate left bidiagonalizing vectors in U
*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
*              (RWorkspace: 0)
*
               CALL ZLACPY( 'L', M, M, A, LDA, U, LDU )
               CALL ZUNGBR( 'Q', M, M, N, U, LDU, WORK( ITAUQ ),
     $                      WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTVAS ) THEN
*
*              If right singular vectors desired in VT, copy result to
*              VT and generate right bidiagonalizing vectors in VT
*              (CWorkspace: need 2*M+NRVT, prefer 2*M+NRVT*NB)
*              (RWorkspace: 0)
*
               CALL ZLACPY( 'U', M, N, A, LDA, VT, LDVT )
               IF( WNTVA )
     $            NRVT = N
               IF( WNTVS )
     $            NRVT = M
               CALL ZUNGBR( 'P', NRVT, N, M, VT, LDVT, WORK( ITAUP ),
     $                      WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTUO ) THEN
*
*              If left singular vectors desired in A, generate left
*              bidiagonalizing vectors in A
*              (CWorkspace: need 3*M-1, prefer 2*M+(M-1)*NB)
*              (RWorkspace: 0)
*
               CALL ZUNGBR( 'Q', M, M, N, A, LDA, WORK( ITAUQ ),
     $                      WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IF( WNTVO ) THEN
*
*              If right singular vectors desired in A, generate right
*              bidiagonalizing vectors in A
*              (CWorkspace: need 3*M, prefer 2*M+M*NB)
*              (RWorkspace: 0)
*
               CALL ZUNGBR( 'P', M, N, M, A, LDA, WORK( ITAUP ),
     $                      WORK( IWORK ), LWORK-IWORK+1, IERR )
            END IF
            IRWORK = IE + M
            IF( WNTUAS .OR. WNTUO )
     $         NRU = M
            IF( WNTUN )
     $         NRU = 0
            IF( WNTVAS .OR. WNTVO )
     $         NCVT = N
            IF( WNTVN )
     $         NCVT = 0
            IF( ( .NOT.WNTUO ) .AND. ( .NOT.WNTVO ) ) THEN
*
*              Perform bidiagonal QR iteration, if desired, computing
*              left singular vectors in U and computing right singular
*              vectors in VT
*              (CWorkspace: 0)
*              (RWorkspace: need BDSPAC)
*
               CALL ZBDSQR( 'L', M, NCVT, NRU, 0, S, RWORK( IE ), VT,
     $                      LDVT, U, LDU, CDUM, 1, RWORK( IRWORK ),
     $                      INFO )
            ELSE IF( ( .NOT.WNTUO ) .AND. WNTVO ) THEN
*
*              Perform bidiagonal QR iteration, if desired, computing
*              left singular vectors in U and computing right singular
*              vectors in A
*              (CWorkspace: 0)
*              (RWorkspace: need BDSPAC)
*
               CALL ZBDSQR( 'L', M, NCVT, NRU, 0, S, RWORK( IE ), A,
     $                      LDA, U, LDU, CDUM, 1, RWORK( IRWORK ),
     $                      INFO )
            ELSE
*
*              Perform bidiagonal QR iteration, if desired, computing
*              left singular vectors in A and computing right singular
*              vectors in VT
*              (CWorkspace: 0)
*              (RWorkspace: need BDSPAC)
*
               CALL ZBDSQR( 'L', M, NCVT, NRU, 0, S, RWORK( IE ), VT,
     $                      LDVT, A, LDA, CDUM, 1, RWORK( IRWORK ),
     $                      INFO )
            END IF
*
         END IF
*
      END IF
*
*     Undo scaling if necessary
*
      IF( ISCL.EQ.1 ) THEN
         IF( ANRM.GT.BIGNUM )
     $      CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN, 1, S, MINMN,
     $                   IERR )
         IF( INFO.NE.0 .AND. ANRM.GT.BIGNUM )
     $      CALL DLASCL( 'G', 0, 0, BIGNUM, ANRM, MINMN-1, 1,
     $                   RWORK( IE ), MINMN, IERR )
         IF( ANRM.LT.SMLNUM )
     $      CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN, 1, S, MINMN,
     $                   IERR )
         IF( INFO.NE.0 .AND. ANRM.LT.SMLNUM )
     $      CALL DLASCL( 'G', 0, 0, SMLNUM, ANRM, MINMN-1, 1,
     $                   RWORK( IE ), MINMN, IERR )
      END IF
*
*     Return optimal workspace in WORK(1)
*
      WORK( 1 ) = MAXWRK
*
      RETURN
*
*     End of ZGESVD
*
      END
 

C
C     The following is from linpack, via www.netlib.org
C
      subroutine zgefa(a,lda,n,ipvt,info)
      integer lda,n,ipvt(1),info
      complex*16 a(lda,1)
c
c     zgefa factors a complex*16 matrix by gaussian elimination.
c
c     zgefa is usually called by zgeco, but it can be called
c     directly with a saving in time if  rcond  is not needed.
c     (time for zgeco) = (1 + 9/n)*(time for zgefa) .
c
c     on entry
c
c        a       complex*16(lda, n)
c                the matrix to be factored.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c     on return
c
c        a       an upper triangular matrix and the multipliers
c                which were used to obtain it.
c                the factorization can be written  a = l*u  where
c                l  is a product of permutation and unit lower
c                triangular matrices and  u  is upper triangular.
c
c        ipvt    integer(n)
c                an integer vector of pivot indices.
c
c        info    integer
c                = 0  normal value.
c                = k  if  u(k,k) .eq. 0.0 .  this is not an error
c                     condition for this subroutine, but it does
c                     indicate that zgesl or zgedi will divide by zero
c                     if called.  use  rcond  in zgeco for a reliable
c                     indication of singularity.
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,izamax
c     fortran dabs
c
c     internal variables
c
      complex*16 t
      integer izamax,j,k,kp1,l,nm1
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     gaussian elimination with partial pivoting
c
      info = 0
      nm1 = n - 1
      if (nm1 .lt. 1) go to 70
      do 60 k = 1, nm1
         kp1 = k + 1
c
c        find l = pivot index
c
         l = izamax(n-k+1,a(k,k),1) + k - 1
         ipvt(k) = l
c
c        zero pivot implies this column already triangularized
c
         if (cabs1(a(l,k)) .eq. 0.0d0) go to 40
c
c           interchange if necessary
c
            if (l .eq. k) go to 10
               t = a(l,k)
               a(l,k) = a(k,k)
               a(k,k) = t
   10       continue
c
c           compute multipliers
c
            t = -(1.0d0,0.0d0)/a(k,k)
            call zscal(n-k,t,a(k+1,k),1)
c
c           row elimination with column indexing
c
            do 30 j = kp1, n
               t = a(l,j)
               if (l .eq. k) go to 20
                  a(l,j) = a(k,j)
                  a(k,j) = t
   20          continue
               call zaxpy(n-k,t,a(k+1,k),1,a(k+1,j),1)
   30       continue
         go to 50
   40    continue
            info = k
   50    continue
   60 continue
   70 continue
      ipvt(n) = n
      if (cabs1(a(n,n)) .eq. 0.0d0) info = n
      return
      end


      subroutine zgedi(a,lda,n,ipvt,det,work,job)
      integer lda,n,ipvt(1),job
      complex*16 a(lda,1),det(2),work(1)
c
c     zgedi computes the determinant and inverse of a matrix
c     using the factors computed by zgeco or zgefa.
c
c     on entry
c
c        a       complex*16(lda, n)
c                the output from zgeco or zgefa.
c
c        lda     integer
c                the leading dimension of the array  a .
c
c        n       integer
c                the order of the matrix  a .
c
c        ipvt    integer(n)
c                the pivot vector from zgeco or zgefa.
c
c        work    complex*16(n)
c                work vector.  contents destroyed.
c
c        job     integer
c                = 11   both determinant and inverse.
c                = 01   inverse only.
c                = 10   determinant only.
c
c     on return
c
c        a       inverse of original matrix if requested.
c                otherwise unchanged.
c
c        det     complex*16(2)
c                determinant of original matrix if requested.
c                otherwise not referenced.
c                determinant = det(1) * 10.0**det(2)
c                with  1.0 .le. cabs1(det(1)) .lt. 10.0
c                or  det(1) .eq. 0.0 .
c
c     error condition
c
c        a division by zero will occur if the input factor contains
c        a zero on the diagonal and the inverse is requested.
c        it will not occur if the subroutines are called correctly
c        and if zgeco has set rcond .gt. 0.0 or zgefa has set
c        info .eq. 0 .
c
c     linpack. this version dated 08/14/78 .
c     cleve moler, university of new mexico, argonne national lab.
c
c     subroutines and functions
c
c     blas zaxpy,zscal,zswap
c     fortran dabs,dcmplx,mod
c
c     internal variables
c
      complex*16 t
      double precision ten
      integer i,j,k,kb,kp1,l,nm1
c
      complex*16 zdum
      double precision cabs1
      double precision dreal,dimag
      complex*16 zdumr,zdumi
      dreal(zdumr) = zdumr
      dimag(zdumi) = (0.0d0,-1.0d0)*zdumi
      cabs1(zdum) = dabs(dreal(zdum)) + dabs(dimag(zdum))
c
c     compute determinant
c
      if (job/10 .eq. 0) go to 70
         det(1) = (1.0d0,0.0d0)
         det(2) = (0.0d0,0.0d0)
         ten = 10.0d0
         do 50 i = 1, n
            if (ipvt(i) .ne. i) det(1) = -det(1)
            det(1) = a(i,i)*det(1)
c        ...exit
            if (cabs1(det(1)) .eq. 0.0d0) go to 60
   10       if (cabs1(det(1)) .ge. 1.0d0) go to 20
               det(1) = dcmplx(ten,0.0d0)*det(1)
               det(2) = det(2) - (1.0d0,0.0d0)
            go to 10
   20       continue
   30       if (cabs1(det(1)) .lt. ten) go to 40
               det(1) = det(1)/dcmplx(ten,0.0d0)
               det(2) = det(2) + (1.0d0,0.0d0)
            go to 30
   40       continue
   50    continue
   60    continue
   70 continue
c
c     compute inverse(u)
c
      if (mod(job,10) .eq. 0) go to 150
         do 100 k = 1, n
            a(k,k) = (1.0d0,0.0d0)/a(k,k)
            t = -a(k,k)
            call zscal(k-1,t,a(1,k),1)
            kp1 = k + 1
            if (n .lt. kp1) go to 90
            do 80 j = kp1, n
               t = a(k,j)
               a(k,j) = (0.0d0,0.0d0)
               call zaxpy(k,t,a(1,k),1,a(1,j),1)
   80       continue
   90       continue
  100    continue
c
c        form inverse(u)*inverse(l)
c
         nm1 = n - 1
         if (nm1 .lt. 1) go to 140
         do 130 kb = 1, nm1
            k = n - kb
            kp1 = k + 1
            do 110 i = kp1, n
               work(i) = a(i,k)
               a(i,k) = (0.0d0,0.0d0)
  110       continue
            do 120 j = kp1, n
               t = work(j)
               call zaxpy(n,t,a(1,j),1,a(1,k),1)
  120       continue
            l = ipvt(k)
            if (l .ne. k) call zswap(n,a(1,k),1,a(1,l),1)
  130    continue
  140    continue
  150 continue
      return
      end

      SUBROUTINE ZBDSQR( UPLO, N, NCVT, NRU, NCC, D, E, VT, LDVT, U,
     $                   LDU, C, LDC, RWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            INFO, LDC, LDU, LDVT, N, NCC, NCVT, NRU
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * ), RWORK( * )
      COMPLEX*16         C( LDC, * ), U( LDU, * ), VT( LDVT, * )
*     ..
*
*  Purpose
*  =======
*
*  ZBDSQR computes the singular value decomposition (SVD) of a real
*  N-by-N (upper or lower) bidiagonal matrix B:  B = Q * S * P' (P'
*  denotes the transpose of P), where S is a diagonal matrix with
*  non-negative diagonal elements (the singular values of B), and Q
*  and P are orthogonal matrices.
*
*  The routine computes S, and optionally computes U * Q, P' * VT,
*  or Q' * C, for given complex input matrices U, VT, and C.
*
*  See "Computing  Small Singular Values of Bidiagonal Matrices With
*  Guaranteed High Relative Accuracy," by J. Demmel and W. Kahan,
*  LAPACK Working Note #3 (or SIAM J. Sci. Statist. Comput. vol. 11,
*  no. 5, pp. 873-912, Sept 1990) and
*  "Accurate singular values and differential qd algorithms," by
*  B. Parlett and V. Fernando, Technical Report CPAM-554, Mathematics
*  Department, University of California at Berkeley, July 1992
*  for a detailed description of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  B is upper bidiagonal;
*          = 'L':  B is lower bidiagonal.
*
*  N       (input) INTEGER
*          The order of the matrix B.  N >= 0.
*
*  NCVT    (input) INTEGER
*          The number of columns of the matrix VT. NCVT >= 0.
*
*  NRU     (input) INTEGER
*          The number of rows of the matrix U. NRU >= 0.
*
*  NCC     (input) INTEGER
*          The number of columns of the matrix C. NCC >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the bidiagonal matrix B.
*          On exit, if INFO=0, the singular values of B in decreasing
*          order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the elements of E contain the
*          offdiagonal elements of of the bidiagonal matrix whose SVD
*          is desired. On normal exit (INFO = 0), E is destroyed.
*          If the algorithm does not converge (INFO > 0), D and E
*          will contain the diagonal and superdiagonal elements of a
*          bidiagonal matrix orthogonally equivalent to the one given
*          as input. E(N) is used for workspace.
*
*  VT      (input/output) COMPLEX*16 array, dimension (LDVT, NCVT)
*          On entry, an N-by-NCVT matrix VT.
*          On exit, VT is overwritten by P' * VT.
*          VT is not referenced if NCVT = 0.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT.
*          LDVT >= max(1,N) if NCVT > 0; LDVT >= 1 if NCVT = 0.
*
*  U       (input/output) COMPLEX*16 array, dimension (LDU, N)
*          On entry, an NRU-by-N matrix U.
*          On exit, U is overwritten by U * Q.
*          U is not referenced if NRU = 0.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= max(1,NRU).
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC, NCC)
*          On entry, an N-by-NCC matrix C.
*          On exit, C is overwritten by Q' * C.
*          C is not referenced if NCC = 0.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C.
*          LDC >= max(1,N) if NCC > 0; LDC >=1 if NCC = 0.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (4*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  If INFO = -i, the i-th argument had an illegal value
*          > 0:  the algorithm did not converge; D and E contain the
*                elements of a bidiagonal matrix which is orthogonally
*                similar to the input matrix B;  if INFO = i, i
*                elements of E have not converged to zero.
*
*  Internal Parameters
*  ===================
*
*  TOLMUL  DOUBLE PRECISION, default = max(10,min(100,EPS**(-1/8)))
*          TOLMUL controls the convergence criterion of the QR loop.
*          If it is positive, TOLMUL*EPS is the desired relative
*             precision in the computed singular values.
*          If it is negative, abs(TOLMUL*EPS*sigma_max) is the
*             desired absolute accuracy in the computed singular
*             values (corresponds to relative accuracy
*             abs(TOLMUL*EPS) in the largest singular value.
*          abs(TOLMUL) should be between 1 and 1/EPS, and preferably
*             between 10 (for fast convergence) and .1/EPS
*             (for there to be some accuracy in the results).
*          Default is to lose at either one eighth or 2 of the
*             available decimal digits in each computed singular value
*             (whichever is smaller).
*
*  MAXITR  INTEGER, default = 6
*          MAXITR controls the maximum number of passes of the
*          algorithm through its inner loop. The algorithms stops
*          (and so fails to converge) if the number of passes
*          through the inner loop exceeds MAXITR*N**2.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D0 )
      DOUBLE PRECISION   ONE
      PARAMETER          ( ONE = 1.0D0 )
      DOUBLE PRECISION   NEGONE
      PARAMETER          ( NEGONE = -1.0D0 )
      DOUBLE PRECISION   HNDRTH
      PARAMETER          ( HNDRTH = 0.01D0 )
      DOUBLE PRECISION   TEN
      PARAMETER          ( TEN = 10.0D0 )
      DOUBLE PRECISION   HNDRD
      PARAMETER          ( HNDRD = 100.0D0 )
      DOUBLE PRECISION   MEIGTH
      PARAMETER          ( MEIGTH = -0.125D0 )
      INTEGER            MAXITR
      PARAMETER          ( MAXITR = 6 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, ROTATE
      INTEGER            I, IDIR, ISUB, ITER, J, LL, LLL, M, MAXIT, NM1,
     $                   NM12, NM13, OLDLL, OLDM
      DOUBLE PRECISION   ABSE, ABSS, COSL, COSR, CS, EPS, F, G, H, MU,
     $                   OLDCS, OLDSN, R, SHIFT, SIGMN, SIGMX, SINL,
     $                   SINR, SLL, SMAX, SMIN, SMINL, SMINLO, SMINOA,
     $                   SN, THRESH, TOL, TOLMUL, UNFL
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           LSAME, DLAMCH
*     ..
*     .. External Subroutines ..
      EXTERNAL           DLARTG, DLAS2, DLASQ1, DLASV2, XERBLA, ZDROT,
     $                   ZDSCAL, ZLASR, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SIGN, SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      LOWER = LSAME( UPLO, 'L' )
      IF( .NOT.LSAME( UPLO, 'U' ) .AND. .NOT.LOWER ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( NCVT.LT.0 ) THEN
         INFO = -3
      ELSE IF( NRU.LT.0 ) THEN
         INFO = -4
      ELSE IF( NCC.LT.0 ) THEN
         INFO = -5
      ELSE IF( ( NCVT.EQ.0 .AND. LDVT.LT.1 ) .OR.
     $         ( NCVT.GT.0 .AND. LDVT.LT.MAX( 1, N ) ) ) THEN
         INFO = -9
      ELSE IF( LDU.LT.MAX( 1, NRU ) ) THEN
         INFO = -11
      ELSE IF( ( NCC.EQ.0 .AND. LDC.LT.1 ) .OR.
     $         ( NCC.GT.0 .AND. LDC.LT.MAX( 1, N ) ) ) THEN
         INFO = -13
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZBDSQR', -INFO )
         RETURN
      END IF
      IF( N.EQ.0 )
     $   RETURN
      IF( N.EQ.1 )
     $   GO TO 160
*
*     ROTATE is true if any singular vectors desired, false otherwise
*
      ROTATE = ( NCVT.GT.0 ) .OR. ( NRU.GT.0 ) .OR. ( NCC.GT.0 )
*
*     If no singular vectors desired, use qd algorithm
*
      IF( .NOT.ROTATE ) THEN
         CALL DLASQ1( N, D, E, RWORK, INFO )
         RETURN
      END IF
*
      NM1 = N - 1
      NM12 = NM1 + NM1
      NM13 = NM12 + NM1
      IDIR = 0
*
*     Get machine constants
*
      EPS = DLAMCH( 'Epsilon' )
      UNFL = DLAMCH( 'Safe minimum' )
*
*     If matrix lower bidiagonal, rotate to be upper bidiagonal
*     by applying Givens rotations on the left
*
      IF( LOWER ) THEN
         DO 10 I = 1, N - 1
            CALL DLARTG( D( I ), E( I ), CS, SN, R )
            D( I ) = R
            E( I ) = SN*D( I+1 )
            D( I+1 ) = CS*D( I+1 )
            RWORK( I ) = CS
            RWORK( NM1+I ) = SN
   10    CONTINUE
*
*        Update singular vectors if desired
*
         IF( NRU.GT.0 )
     $      CALL ZLASR( 'R', 'V', 'F', NRU, N, RWORK( 1 ), RWORK( N ),
     $                  U, LDU )
         IF( NCC.GT.0 )
     $      CALL ZLASR( 'L', 'V', 'F', N, NCC, RWORK( 1 ), RWORK( N ),
     $                  C, LDC )
      END IF
*
*     Compute singular values to relative accuracy TOL
*     (By setting TOL to be negative, algorithm will compute
*     singular values to absolute accuracy ABS(TOL)*norm(input matrix))
*
      TOLMUL = MAX( TEN, MIN( HNDRD, EPS**MEIGTH ) )
      TOL = TOLMUL*EPS
*
*     Compute approximate maximum, minimum singular values
*
      SMAX = ZERO
      DO 20 I = 1, N
         SMAX = MAX( SMAX, ABS( D( I ) ) )
   20 CONTINUE
      DO 30 I = 1, N - 1
         SMAX = MAX( SMAX, ABS( E( I ) ) )
   30 CONTINUE
      SMINL = ZERO
      IF( TOL.GE.ZERO ) THEN
*
*        Relative accuracy desired
*
         SMINOA = ABS( D( 1 ) )
         IF( SMINOA.EQ.ZERO )
     $      GO TO 50
         MU = SMINOA
         DO 40 I = 2, N
            MU = ABS( D( I ) )*( MU / ( MU+ABS( E( I-1 ) ) ) )
            SMINOA = MIN( SMINOA, MU )
            IF( SMINOA.EQ.ZERO )
     $         GO TO 50
   40    CONTINUE
   50    CONTINUE
         SMINOA = SMINOA / SQRT( DBLE( N ) )
         THRESH = MAX( TOL*SMINOA, MAXITR*N*N*UNFL )
      ELSE
*
*        Absolute accuracy desired
*
         THRESH = MAX( ABS( TOL )*SMAX, MAXITR*N*N*UNFL )
      END IF
*
*     Prepare for main iteration loop for the singular values
*     (MAXIT is the maximum number of passes through the inner
*     loop permitted before nonconvergence signalled.)
*
      MAXIT = MAXITR*N*N
      ITER = 0
      OLDLL = -1
      OLDM = -1
*
*     M points to last element of unconverged part of matrix
*
      M = N
*
*     Begin main iteration loop
*
   60 CONTINUE
*
*     Check for convergence or exceeding iteration count
*
      IF( M.LE.1 )
     $   GO TO 160
      IF( ITER.GT.MAXIT )
     $   GO TO 200
*
*     Find diagonal block of matrix to work on
*
      IF( TOL.LT.ZERO .AND. ABS( D( M ) ).LE.THRESH )
     $   D( M ) = ZERO
      SMAX = ABS( D( M ) )
      SMIN = SMAX
      DO 70 LLL = 1, M - 1
         LL = M - LLL
         ABSS = ABS( D( LL ) )
         ABSE = ABS( E( LL ) )
         IF( TOL.LT.ZERO .AND. ABSS.LE.THRESH )
     $      D( LL ) = ZERO
         IF( ABSE.LE.THRESH )
     $      GO TO 80
         SMIN = MIN( SMIN, ABSS )
         SMAX = MAX( SMAX, ABSS, ABSE )
   70 CONTINUE
      LL = 0
      GO TO 90
   80 CONTINUE
      E( LL ) = ZERO
*
*     Matrix splits since E(LL) = 0
*
      IF( LL.EQ.M-1 ) THEN
*
*        Convergence of bottom singular value, return to top of loop
*
         M = M - 1
         GO TO 60
      END IF
   90 CONTINUE
      LL = LL + 1
*
*     E(LL) through E(M-1) are nonzero, E(LL-1) is zero
*
      IF( LL.EQ.M-1 ) THEN
*
*        2 by 2 block, handle separately
*
         CALL DLASV2( D( M-1 ), E( M-1 ), D( M ), SIGMN, SIGMX, SINR,
     $                COSR, SINL, COSL )
         D( M-1 ) = SIGMX
         E( M-1 ) = ZERO
         D( M ) = SIGMN
*
*        Compute singular vectors, if desired
*
         IF( NCVT.GT.0 )
     $      CALL ZDROT( NCVT, VT( M-1, 1 ), LDVT, VT( M, 1 ), LDVT,
     $                  COSR, SINR )
         IF( NRU.GT.0 )
     $      CALL ZDROT( NRU, U( 1, M-1 ), 1, U( 1, M ), 1, COSL, SINL )
         IF( NCC.GT.0 )
     $      CALL ZDROT( NCC, C( M-1, 1 ), LDC, C( M, 1 ), LDC, COSL,
     $                  SINL )
         M = M - 2
         GO TO 60
      END IF
*
*     If working on new submatrix, choose shift direction
*     (from larger end diagonal element towards smaller)
*
      IF( LL.GT.OLDM .OR. M.LT.OLDLL ) THEN
         IF( ABS( D( LL ) ).GE.ABS( D( M ) ) ) THEN
*
*           Chase bulge from top (big end) to bottom (small end)
*
            IDIR = 1
         ELSE
*
*           Chase bulge from bottom (big end) to top (small end)
*
            IDIR = 2
         END IF
      END IF
*
*     Apply convergence tests
*
      IF( IDIR.EQ.1 ) THEN
*
*        Run convergence test in forward direction
*        First apply standard test to bottom of matrix
*
         IF( ABS( E( M-1 ) ).LE.ABS( TOL )*ABS( D( M ) ) .OR.
     $       ( TOL.LT.ZERO .AND. ABS( E( M-1 ) ).LE.THRESH ) ) THEN
            E( M-1 ) = ZERO
            GO TO 60
         END IF
*
         IF( TOL.GE.ZERO ) THEN
*
*           If relative accuracy desired,
*           apply convergence criterion forward
*
            MU = ABS( D( LL ) )
            SMINL = MU
            DO 100 LLL = LL, M - 1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 60
               END IF
               SMINLO = SMINL
               MU = ABS( D( LLL+1 ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
  100       CONTINUE
         END IF
*
      ELSE
*
*        Run convergence test in backward direction
*        First apply standard test to top of matrix
*
         IF( ABS( E( LL ) ).LE.ABS( TOL )*ABS( D( LL ) ) .OR.
     $       ( TOL.LT.ZERO .AND. ABS( E( LL ) ).LE.THRESH ) ) THEN
            E( LL ) = ZERO
            GO TO 60
         END IF
*
         IF( TOL.GE.ZERO ) THEN
*
*           If relative accuracy desired,
*           apply convergence criterion backward
*
            MU = ABS( D( M ) )
            SMINL = MU
            DO 110 LLL = M - 1, LL, -1
               IF( ABS( E( LLL ) ).LE.TOL*MU ) THEN
                  E( LLL ) = ZERO
                  GO TO 60
               END IF
               SMINLO = SMINL
               MU = ABS( D( LLL ) )*( MU / ( MU+ABS( E( LLL ) ) ) )
               SMINL = MIN( SMINL, MU )
  110       CONTINUE
         END IF
      END IF
      OLDLL = LL
      OLDM = M
*
*     Compute shift.  First, test if shifting would ruin relative
*     accuracy, and if so set the shift to zero.
*
      IF( TOL.GE.ZERO .AND. N*TOL*( SMINL / SMAX ).LE.
     $    MAX( EPS, HNDRTH*TOL ) ) THEN
*
*        Use a zero shift to avoid loss of relative accuracy
*
         SHIFT = ZERO
      ELSE
*
*        Compute the shift from 2-by-2 block at end of matrix
*
         IF( IDIR.EQ.1 ) THEN
            SLL = ABS( D( LL ) )
            CALL DLAS2( D( M-1 ), E( M-1 ), D( M ), SHIFT, R )
         ELSE
            SLL = ABS( D( M ) )
            CALL DLAS2( D( LL ), E( LL ), D( LL+1 ), SHIFT, R )
         END IF
*
*        Test if shift negligible, and if so set to zero
*
         IF( SLL.GT.ZERO ) THEN
            IF( ( SHIFT / SLL )**2.LT.EPS )
     $         SHIFT = ZERO
         END IF
      END IF
*
*     Increment iteration count
*
      ITER = ITER + M - LL
*
*     If SHIFT = 0, do simplified QR iteration
*
      IF( SHIFT.EQ.ZERO ) THEN
         IF( IDIR.EQ.1 ) THEN
*
*           Chase bulge from top to bottom
*           Save cosines and sines for later singular vector updates
*
            CS = ONE
            OLDCS = ONE
            DO 120 I = LL, M - 1
               CALL DLARTG( D( I )*CS, E( I ), CS, SN, R )
               IF( I.GT.LL )
     $            E( I-1 ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I+1 )*SN, OLDCS, OLDSN, D( I ) )
               RWORK( I-LL+1 ) = CS
               RWORK( I-LL+1+NM1 ) = SN
               RWORK( I-LL+1+NM12 ) = OLDCS
               RWORK( I-LL+1+NM13 ) = OLDSN
  120       CONTINUE
            H = D( M )*CS
            D( M ) = H*OLDCS
            E( M-1 ) = H*OLDSN
*
*           Update singular vectors
*
            IF( NCVT.GT.0 )
     $         CALL ZLASR( 'L', 'V', 'F', M-LL+1, NCVT, RWORK( 1 ),
     $                     RWORK( N ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL ZLASR( 'R', 'V', 'F', NRU, M-LL+1, RWORK( NM12+1 ),
     $                     RWORK( NM13+1 ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL ZLASR( 'L', 'V', 'F', M-LL+1, NCC, RWORK( NM12+1 ),
     $                     RWORK( NM13+1 ), C( LL, 1 ), LDC )
*
*           Test convergence
*
            IF( ABS( E( M-1 ) ).LE.THRESH )
     $         E( M-1 ) = ZERO
*
         ELSE
*
*           Chase bulge from bottom to top
*           Save cosines and sines for later singular vector updates
*
            CS = ONE
            OLDCS = ONE
            DO 130 I = M, LL + 1, -1
               CALL DLARTG( D( I )*CS, E( I-1 ), CS, SN, R )
               IF( I.LT.M )
     $            E( I ) = OLDSN*R
               CALL DLARTG( OLDCS*R, D( I-1 )*SN, OLDCS, OLDSN, D( I ) )
               RWORK( I-LL ) = CS
               RWORK( I-LL+NM1 ) = -SN
               RWORK( I-LL+NM12 ) = OLDCS
               RWORK( I-LL+NM13 ) = -OLDSN
  130       CONTINUE
            H = D( LL )*CS
            D( LL ) = H*OLDCS
            E( LL ) = H*OLDSN
*
*           Update singular vectors
*
            IF( NCVT.GT.0 )
     $         CALL ZLASR( 'L', 'V', 'B', M-LL+1, NCVT, RWORK( NM12+1 ),
     $                     RWORK( NM13+1 ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL ZLASR( 'R', 'V', 'B', NRU, M-LL+1, RWORK( 1 ),
     $                     RWORK( N ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL ZLASR( 'L', 'V', 'B', M-LL+1, NCC, RWORK( 1 ),
     $                     RWORK( N ), C( LL, 1 ), LDC )
*
*           Test convergence
*
            IF( ABS( E( LL ) ).LE.THRESH )
     $         E( LL ) = ZERO
         END IF
      ELSE
*
*        Use nonzero shift
*
         IF( IDIR.EQ.1 ) THEN
*
*           Chase bulge from top to bottom
*           Save cosines and sines for later singular vector updates
*
            F = ( ABS( D( LL ) )-SHIFT )*
     $          ( SIGN( ONE, D( LL ) )+SHIFT / D( LL ) )
            G = E( LL )
            DO 140 I = LL, M - 1
               CALL DLARTG( F, G, COSR, SINR, R )
               IF( I.GT.LL )
     $            E( I-1 ) = R
               F = COSR*D( I ) + SINR*E( I )
               E( I ) = COSR*E( I ) - SINR*D( I )
               G = SINR*D( I+1 )
               D( I+1 ) = COSR*D( I+1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I ) + SINL*D( I+1 )
               D( I+1 ) = COSL*D( I+1 ) - SINL*E( I )
               IF( I.LT.M-1 ) THEN
                  G = SINL*E( I+1 )
                  E( I+1 ) = COSL*E( I+1 )
               END IF
               RWORK( I-LL+1 ) = COSR
               RWORK( I-LL+1+NM1 ) = SINR
               RWORK( I-LL+1+NM12 ) = COSL
               RWORK( I-LL+1+NM13 ) = SINL
  140       CONTINUE
            E( M-1 ) = F
*
*           Update singular vectors
*
            IF( NCVT.GT.0 )
     $         CALL ZLASR( 'L', 'V', 'F', M-LL+1, NCVT, RWORK( 1 ),
     $                     RWORK( N ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL ZLASR( 'R', 'V', 'F', NRU, M-LL+1, RWORK( NM12+1 ),
     $                     RWORK( NM13+1 ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL ZLASR( 'L', 'V', 'F', M-LL+1, NCC, RWORK( NM12+1 ),
     $                     RWORK( NM13+1 ), C( LL, 1 ), LDC )
*
*           Test convergence
*
            IF( ABS( E( M-1 ) ).LE.THRESH )
     $         E( M-1 ) = ZERO
*
         ELSE
*
*           Chase bulge from bottom to top
*           Save cosines and sines for later singular vector updates
*
            F = ( ABS( D( M ) )-SHIFT )*( SIGN( ONE, D( M ) )+SHIFT /
     $          D( M ) )
            G = E( M-1 )
            DO 150 I = M, LL + 1, -1
               CALL DLARTG( F, G, COSR, SINR, R )
               IF( I.LT.M )
     $            E( I ) = R
               F = COSR*D( I ) + SINR*E( I-1 )
               E( I-1 ) = COSR*E( I-1 ) - SINR*D( I )
               G = SINR*D( I-1 )
               D( I-1 ) = COSR*D( I-1 )
               CALL DLARTG( F, G, COSL, SINL, R )
               D( I ) = R
               F = COSL*E( I-1 ) + SINL*D( I-1 )
               D( I-1 ) = COSL*D( I-1 ) - SINL*E( I-1 )
               IF( I.GT.LL+1 ) THEN
                  G = SINL*E( I-2 )
                  E( I-2 ) = COSL*E( I-2 )
               END IF
               RWORK( I-LL ) = COSR
               RWORK( I-LL+NM1 ) = -SINR
               RWORK( I-LL+NM12 ) = COSL
               RWORK( I-LL+NM13 ) = -SINL
  150       CONTINUE
            E( LL ) = F
*
*           Test convergence
*
            IF( ABS( E( LL ) ).LE.THRESH )
     $         E( LL ) = ZERO
*
*           Update singular vectors if desired
*
            IF( NCVT.GT.0 )
     $         CALL ZLASR( 'L', 'V', 'B', M-LL+1, NCVT, RWORK( NM12+1 ),
     $                     RWORK( NM13+1 ), VT( LL, 1 ), LDVT )
            IF( NRU.GT.0 )
     $         CALL ZLASR( 'R', 'V', 'B', NRU, M-LL+1, RWORK( 1 ),
     $                     RWORK( N ), U( 1, LL ), LDU )
            IF( NCC.GT.0 )
     $         CALL ZLASR( 'L', 'V', 'B', M-LL+1, NCC, RWORK( 1 ),
     $                     RWORK( N ), C( LL, 1 ), LDC )
         END IF
      END IF
*
*     QR iteration finished, go back and check convergence
*
      GO TO 60
*
*     All singular values converged, so make them positive
*
  160 CONTINUE
      DO 170 I = 1, N
         IF( D( I ).LT.ZERO ) THEN
            D( I ) = -D( I )
*
*           Change sign of singular vectors, if desired
*
            IF( NCVT.GT.0 )
     $         CALL ZDSCAL( NCVT, NEGONE, VT( I, 1 ), LDVT )
         END IF
  170 CONTINUE
*
*     Sort the singular values into decreasing order (insertion sort on
*     singular values, but only one transposition per singular vector)
*
      DO 190 I = 1, N - 1
*
*        Scan for smallest D(I)
*
         ISUB = 1
         SMIN = D( 1 )
         DO 180 J = 2, N + 1 - I
            IF( D( J ).LE.SMIN ) THEN
               ISUB = J
               SMIN = D( J )
            END IF
  180    CONTINUE
         IF( ISUB.NE.N+1-I ) THEN
*
*           Swap singular values and vectors
*
            D( ISUB ) = D( N+1-I )
            D( N+1-I ) = SMIN
            IF( NCVT.GT.0 )
     $         CALL ZSWAP( NCVT, VT( ISUB, 1 ), LDVT, VT( N+1-I, 1 ),
     $                     LDVT )
            IF( NRU.GT.0 )
     $         CALL ZSWAP( NRU, U( 1, ISUB ), 1, U( 1, N+1-I ), 1 )
            IF( NCC.GT.0 )
     $         CALL ZSWAP( NCC, C( ISUB, 1 ), LDC, C( N+1-I, 1 ), LDC )
         END IF
  190 CONTINUE
      GO TO 220
*
*     Maximum number of iterations exceeded, failure to converge
*
  200 CONTINUE
      INFO = 0
      DO 210 I = 1, N - 1
         IF( E( I ).NE.ZERO )
     $      INFO = INFO + 1
  210 CONTINUE
  220 CONTINUE
      RETURN
*
*     End of ZBDSQR
*
      END

      SUBROUTINE ZGEBD2( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEBD2 reduces a complex general m by n matrix A to upper or lower
*  real bidiagonal form B by a unitary transformation: Q' * A * P = B.
*
*  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows in the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns in the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n general matrix to be reduced.
*          On exit,
*          if m >= n, the diagonal and the first superdiagonal are
*            overwritten with the upper bidiagonal matrix B; the
*            elements below the diagonal, with the array TAUQ, represent
*            the unitary matrix Q as a product of elementary
*            reflectors, and the elements above the first superdiagonal,
*            with the array TAUP, represent the unitary matrix P as
*            a product of elementary reflectors;
*          if m < n, the diagonal and the first subdiagonal are
*            overwritten with the lower bidiagonal matrix B; the
*            elements below the first subdiagonal, with the array TAUQ,
*            represent the unitary matrix Q as a product of
*            elementary reflectors, and the elements above the diagonal,
*            with the array TAUP, represent the unitary matrix P as
*            a product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The diagonal elements of the bidiagonal matrix B:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
*          The off-diagonal elements of the bidiagonal matrix B:
*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
*
*  TAUQ    (output) COMPLEX*16 array dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the unitary matrix Q. See Further Details.
*
*  TAUP    (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the unitary matrix P. See Further Details.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (max(M,N))
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*  If m >= n,
*
*     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are complex scalars, and v and u are complex
*  vectors; v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in
*  A(i+1:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in
*  A(i,i+2:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  If m < n,
*
*     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are complex scalars, v and u are complex vectors;
*  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
*  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
*  tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  The contents of A on exit are illustrated by the following examples:
*
*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*
*    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
*    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
*    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
*    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
*    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
*    (  v1  v2  v3  v4  v5 )
*
*  where d and e denote diagonal and off-diagonal elements of B, vi
*  denotes an element of the vector defining H(i), and ui an element of
*  the vector defining G(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLARF, ZLARFG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'ZGEBD2', -INFO )
         RETURN
      END IF
*
      IF( M.GE.N ) THEN
*
*        Reduce to upper bidiagonal form
*
         DO 10 I = 1, N
*
*           Generate elementary reflector H(i) to annihilate A(i+1:m,i)
*
            ALPHA = A( I, I )
            CALL ZLARFG( M-I+1, ALPHA, A( MIN( I+1, M ), I ), 1,
     $                   TAUQ( I ) )
            D( I ) = ALPHA
            A( I, I ) = ONE
*
*           Apply H(i)' to A(i:m,i+1:n) from the left
*
            CALL ZLARF( 'Left', M-I+1, N-I, A( I, I ), 1,
     $                  DCONJG( TAUQ( I ) ), A( I, I+1 ), LDA, WORK )
            A( I, I ) = D( I )
*
            IF( I.LT.N ) THEN
*
*              Generate elementary reflector G(i) to annihilate
*              A(i,i+2:n)
*
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               ALPHA = A( I, I+1 )
               CALL ZLARFG( N-I, ALPHA, A( I, MIN( I+2, N ) ), LDA,
     $                      TAUP( I ) )
               E( I ) = ALPHA
               A( I, I+1 ) = ONE
*
*              Apply G(i) to A(i+1:m,i+1:n) from the right
*
               CALL ZLARF( 'Right', M-I, N-I, A( I, I+1 ), LDA,
     $                     TAUP( I ), A( I+1, I+1 ), LDA, WORK )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               A( I, I+1 ) = E( I )
            ELSE
               TAUP( I ) = ZERO
            END IF
   10    CONTINUE
      ELSE
*
*        Reduce to lower bidiagonal form
*
         DO 20 I = 1, M
*
*           Generate elementary reflector G(i) to annihilate A(i,i+1:n)
*
            CALL ZLACGV( N-I+1, A( I, I ), LDA )
            ALPHA = A( I, I )
            CALL ZLARFG( N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA,
     $                   TAUP( I ) )
            D( I ) = ALPHA
            A( I, I ) = ONE
*
*           Apply G(i) to A(i+1:m,i:n) from the right
*
            CALL ZLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAUP( I ),
     $                  A( MIN( I+1, M ), I ), LDA, WORK )
            CALL ZLACGV( N-I+1, A( I, I ), LDA )
            A( I, I ) = D( I )
*
            IF( I.LT.M ) THEN
*
*              Generate elementary reflector H(i) to annihilate
*              A(i+2:m,i)
*
               ALPHA = A( I+1, I )
               CALL ZLARFG( M-I, ALPHA, A( MIN( I+2, M ), I ), 1,
     $                      TAUQ( I ) )
               E( I ) = ALPHA
               A( I+1, I ) = ONE
*
*              Apply H(i)' to A(i+1:m,i+1:n) from the left
*
               CALL ZLARF( 'Left', M-I, N-I, A( I+1, I ), 1,
     $                     DCONJG( TAUQ( I ) ), A( I+1, I+1 ), LDA,
     $                     WORK )
               A( I+1, I ) = E( I )
            ELSE
               TAUQ( I ) = ZERO
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZGEBD2
*
      END
      SUBROUTINE ZGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK,
     $                   INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGEBRD reduces a general complex M-by-N matrix A to upper or lower
*  bidiagonal form B by a unitary transformation: Q**H * A * P = B.
*
*  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows in the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns in the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N general matrix to be reduced.
*          On exit,
*          if m >= n, the diagonal and the first superdiagonal are
*            overwritten with the upper bidiagonal matrix B; the
*            elements below the diagonal, with the array TAUQ, represent
*            the unitary matrix Q as a product of elementary
*            reflectors, and the elements above the first superdiagonal,
*            with the array TAUP, represent the unitary matrix P as
*            a product of elementary reflectors;
*          if m < n, the diagonal and the first subdiagonal are
*            overwritten with the lower bidiagonal matrix B; the
*            elements below the first subdiagonal, with the array TAUQ,
*            represent the unitary matrix Q as a product of
*            elementary reflectors, and the elements above the diagonal,
*            with the array TAUP, represent the unitary matrix P as
*            a product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The diagonal elements of the bidiagonal matrix B:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
*          The off-diagonal elements of the bidiagonal matrix B:
*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
*
*  TAUQ    (output) COMPLEX*16 array dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the unitary matrix Q. See Further Details.
*
*  TAUP    (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the unitary matrix P. See Further Details.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,M,N).
*          For optimum performance LWORK >= (M+N)*NB, where NB
*          is the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*  If m >= n,
*
*     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are complex scalars, and v and u are complex
*  vectors; v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in
*  A(i+1:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in
*  A(i,i+2:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  If m < n,
*
*     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are complex scalars, and v and u are complex
*  vectors; v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in
*  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in
*  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  The contents of A on exit are illustrated by the following examples:
*
*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*
*    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
*    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
*    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
*    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
*    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
*    (  v1  v2  v3  v4  v5 )
*
*  where d and e denote diagonal and off-diagonal elements of B, vi
*  denotes an element of the vector defining H(i), and ui an element of
*  the vector defining G(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IINFO, J, LDWRKX, LDWRKY, LWKOPT, MINMN, NB,
     $                   NBMIN, NX
      DOUBLE PRECISION   WS
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEBD2, ZGEMM, ZLABRD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INFO = 0
      NB = MAX( 1, ILAENV( 1, 'ZGEBRD', ' ', M, N, -1, -1 ) )
      LWKOPT = ( M+N )*NB
      WORK( 1 ) = DBLE( LWKOPT )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'ZGEBRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      MINMN = MIN( M, N )
      IF( MINMN.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      WS = MAX( M, N )
      LDWRKX = M
      LDWRKY = N
*
      IF( NB.GT.1 .AND. NB.LT.MINMN ) THEN
*
*        Set the crossover point NX.
*
         NX = MAX( NB, ILAENV( 3, 'ZGEBRD', ' ', M, N, -1, -1 ) )
*
*        Determine when to switch from blocked to unblocked code.
*
         IF( NX.LT.MINMN ) THEN
            WS = ( M+N )*NB
            IF( LWORK.LT.WS ) THEN
*
*              Not enough work space for the optimal NB, consider using
*              a smaller block size.
*
               NBMIN = ILAENV( 2, 'ZGEBRD', ' ', M, N, -1, -1 )
               IF( LWORK.GE.( M+N )*NBMIN ) THEN
                  NB = LWORK / ( M+N )
               ELSE
                  NB = 1
                  NX = MINMN
               END IF
            END IF
         END IF
      ELSE
         NX = MINMN
      END IF
*
      DO 30 I = 1, MINMN - NX, NB
*
*        Reduce rows and columns i:i+ib-1 to bidiagonal form and return
*        the matrices X and Y which are needed to update the unreduced
*        part of the matrix
*
         CALL ZLABRD( M-I+1, N-I+1, NB, A( I, I ), LDA, D( I ), E( I ),
     $                TAUQ( I ), TAUP( I ), WORK, LDWRKX,
     $                WORK( LDWRKX*NB+1 ), LDWRKY )
*
*        Update the trailing submatrix A(i+ib:m,i+ib:n), using
*        an update of the form  A := A - V*Y' - X*U'
*
         CALL ZGEMM( 'No transpose', 'Conjugate transpose', M-I-NB+1,
     $               N-I-NB+1, NB, -ONE, A( I+NB, I ), LDA,
     $               WORK( LDWRKX*NB+NB+1 ), LDWRKY, ONE,
     $               A( I+NB, I+NB ), LDA )
         CALL ZGEMM( 'No transpose', 'No transpose', M-I-NB+1, N-I-NB+1,
     $               NB, -ONE, WORK( NB+1 ), LDWRKX, A( I, I+NB ), LDA,
     $               ONE, A( I+NB, I+NB ), LDA )
*
*        Copy diagonal and off-diagonal elements of B back into A
*
         IF( M.GE.N ) THEN
            DO 10 J = I, I + NB - 1
               A( J, J ) = D( J )
               A( J, J+1 ) = E( J )
   10       CONTINUE
         ELSE
            DO 20 J = I, I + NB - 1
               A( J, J ) = D( J )
               A( J+1, J ) = E( J )
   20       CONTINUE
         END IF
   30 CONTINUE
*
*     Use unblocked code to reduce the remainder of the matrix
*
      CALL ZGEBD2( M-I+1, N-I+1, A( I, I ), LDA, D( I ), E( I ),
     $             TAUQ( I ), TAUP( I ), WORK, IINFO )
      WORK( 1 ) = WS
      RETURN
*
*     End of ZGEBRD
*
      END
      SUBROUTINE ZGELQ2( M, N, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGELQ2 computes an LQ factorization of a complex m by n matrix A:
*  A = L * Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and below the diagonal of the array
*          contain the m by min(m,n) lower trapezoidal matrix L (L is
*          lower triangular if m <= n); the elements above the diagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (M)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(k)' . . . H(2)' H(1)', where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
*  A(i,i+1:n), and tau in TAU(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, K
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLARF, ZLARFG
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGELQ2', -INFO )
         RETURN
      END IF
*
      K = MIN( M, N )
*
      DO 10 I = 1, K
*
*        Generate elementary reflector H(i) to annihilate A(i,i+1:n)
*
         CALL ZLACGV( N-I+1, A( I, I ), LDA )
         ALPHA = A( I, I )
         CALL ZLARFG( N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA,
     $                TAU( I ) )
         IF( I.LT.M ) THEN
*
*           Apply H(i) to A(i+1:m,i:n) from the right
*
            A( I, I ) = ONE
            CALL ZLARF( 'Right', M-I, N-I+1, A( I, I ), LDA, TAU( I ),
     $                  A( I+1, I ), LDA, WORK )
         END IF
         A( I, I ) = ALPHA
         CALL ZLACGV( N-I+1, A( I, I ), LDA )
   10 CONTINUE
      RETURN
*
*     End of ZGELQ2
*
      END
      SUBROUTINE ZGELQF( M, N, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGELQF computes an LQ factorization of a complex M-by-N matrix A:
*  A = L * Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and below the diagonal of the array
*          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
*          lower triangular if m <= n); the elements above the diagonal,
*          with the array TAU, represent the unitary matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) COMPLEX*16 array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,M).
*          For optimum performance LWORK >= M*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(k)' . . . H(2)' H(1)', where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a complex scalar, and v is a complex vector with
*  v(1:i-1) = 0 and v(i) = 1; conjg(v(i+1:n)) is stored on exit in
*  A(i,i+1:n), and tau in TAU(i).
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, K, LDWORK, LWKOPT, NB,
     $                   NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGELQ2, ZLARFB, ZLARFT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'ZGELQF', ' ', M, N, -1, -1 )
      LWKOPT = M*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -7
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGELQF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      K = MIN( M, N )
      IF( K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'ZGELQF', ' ', M, N, -1, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZGELQF', ' ', M, N, -1,
     $                 -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code initially
*
         DO 10 I = 1, K - NX, NB
            IB = MIN( K-I+1, NB )
*
*           Compute the LQ factorization of the current block
*           A(i:i+ib-1,i:n)
*
            CALL ZGELQ2( IB, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
            IF( I+IB.LE.M ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL ZLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ),
     $                      LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H to A(i+ib:m,i:n) from the right
*
               CALL ZLARFB( 'Right', 'No transpose', 'Forward',
     $                      'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ),
     $                      LDA, WORK, LDWORK, A( I+IB, I ), LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
   10    CONTINUE
      ELSE
         I = 1
      END IF
*
*     Use unblocked code to factor the last or only block.
*
      IF( I.LE.K )
     $   CALL ZGELQ2( M-I+1, N-I+1, A( I, I ), LDA, TAU( I ), WORK,
     $                IINFO )
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZGELQF
*
      END

      SUBROUTINE ZLABRD( M, N, NB, A, LDA, D, E, TAUQ, TAUP, X, LDX, Y,
     $                   LDY )
*
*  -- LAPACK auxiliary routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            LDA, LDX, LDY, M, N, NB
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         A( LDA, * ), TAUP( * ), TAUQ( * ), X( LDX, * ),
     $                   Y( LDY, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLABRD reduces the first NB rows and columns of a complex general
*  m by n matrix A to upper or lower real bidiagonal form by a unitary
*  transformation Q' * A * P, and returns the matrices X and Y which
*  are needed to apply the transformation to the unreduced part of A.
*
*  If m >= n, A is reduced to upper bidiagonal form; if m < n, to lower
*  bidiagonal form.
*
*  This is an auxiliary routine called by ZGEBRD
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows in the matrix A.
*
*  N       (input) INTEGER
*          The number of columns in the matrix A.
*
*  NB      (input) INTEGER
*          The number of leading rows and columns of A to be reduced.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n general matrix to be reduced.
*          On exit, the first NB rows and columns of the matrix are
*          overwritten; the rest of the array is unchanged.
*          If m >= n, elements on and below the diagonal in the first NB
*            columns, with the array TAUQ, represent the unitary
*            matrix Q as a product of elementary reflectors; and
*            elements above the diagonal in the first NB rows, with the
*            array TAUP, represent the unitary matrix P as a product
*            of elementary reflectors.
*          If m < n, elements below the diagonal in the first NB
*            columns, with the array TAUQ, represent the unitary
*            matrix Q as a product of elementary reflectors, and
*            elements on and above the diagonal in the first NB rows,
*            with the array TAUP, represent the unitary matrix P as
*            a product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  D       (output) DOUBLE PRECISION array, dimension (NB)
*          The diagonal elements of the first NB rows and columns of
*          the reduced matrix.  D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (NB)
*          The off-diagonal elements of the first NB rows and columns of
*          the reduced matrix.
*
*  TAUQ    (output) COMPLEX*16 array dimension (NB)
*          The scalar factors of the elementary reflectors which
*          represent the unitary matrix Q. See Further Details.
*
*  TAUP    (output) COMPLEX*16 array, dimension (NB)
*          The scalar factors of the elementary reflectors which
*          represent the unitary matrix P. See Further Details.
*
*  X       (output) COMPLEX*16 array, dimension (LDX,NB)
*          The m-by-nb matrix X required to update the unreduced part
*          of A.
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X. LDX >= max(1,M).
*
*  Y       (output) COMPLEX*16 array, dimension (LDY,NB)
*          The n-by-nb matrix Y required to update the unreduced part
*          of A.
*
*  LDY     (output) INTEGER
*          The leading dimension of the array Y. LDY >= max(1,N).
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*     Q = H(1) H(2) . . . H(nb)  and  P = G(1) G(2) . . . G(nb)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are complex scalars, and v and u are complex
*  vectors.
*
*  If m >= n, v(1:i-1) = 0, v(i) = 1, and v(i:m) is stored on exit in
*  A(i:m,i); u(1:i) = 0, u(i+1) = 1, and u(i+1:n) is stored on exit in
*  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  If m < n, v(1:i) = 0, v(i+1) = 1, and v(i+1:m) is stored on exit in
*  A(i+2:m,i); u(1:i-1) = 0, u(i) = 1, and u(i:n) is stored on exit in
*  A(i,i+1:n); tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  The elements of the vectors v and u together form the m-by-nb matrix
*  V and the nb-by-n matrix U' which are needed, with X and Y, to apply
*  the transformation to the unreduced part of the matrix, using a block
*  update of the form:  A := A - V*Y' - X*U'.
*
*  The contents of A on exit are illustrated by the following examples
*  with nb = 2:
*
*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*
*    (  1   1   u1  u1  u1 )           (  1   u1  u1  u1  u1  u1 )
*    (  v1  1   1   u2  u2 )           (  1   1   u2  u2  u2  u2 )
*    (  v1  v2  a   a   a  )           (  v1  1   a   a   a   a  )
*    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
*    (  v1  v2  a   a   a  )           (  v1  v2  a   a   a   a  )
*    (  v1  v2  a   a   a  )
*
*  where a denotes an element of the original matrix which is unchanged,
*  vi denotes an element of the vector defining H(i), and ui an element
*  of the vector defining G(i).
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I
      COMPLEX*16         ALPHA
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZGEMV, ZLACGV, ZLARFG, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MIN
*     ..
*     .. Executable Statements ..
*
*     Quick return if possible
*
      IF( M.LE.0 .OR. N.LE.0 )
     $   RETURN
*
      IF( M.GE.N ) THEN
*
*        Reduce to upper bidiagonal form
*
         DO 10 I = 1, NB
*
*           Update A(i:m,i)
*
            CALL ZLACGV( I-1, Y( I, 1 ), LDY )
            CALL ZGEMV( 'No transpose', M-I+1, I-1, -ONE, A( I, 1 ),
     $                  LDA, Y( I, 1 ), LDY, ONE, A( I, I ), 1 )
            CALL ZLACGV( I-1, Y( I, 1 ), LDY )
            CALL ZGEMV( 'No transpose', M-I+1, I-1, -ONE, X( I, 1 ),
     $                  LDX, A( 1, I ), 1, ONE, A( I, I ), 1 )
*
*           Generate reflection Q(i) to annihilate A(i+1:m,i)
*
            ALPHA = A( I, I )
            CALL ZLARFG( M-I+1, ALPHA, A( MIN( I+1, M ), I ), 1,
     $                   TAUQ( I ) )
            D( I ) = ALPHA
            IF( I.LT.N ) THEN
               A( I, I ) = ONE
*
*              Compute Y(i+1:n,i)
*
               CALL ZGEMV( 'Conjugate transpose', M-I+1, N-I, ONE,
     $                     A( I, I+1 ), LDA, A( I, I ), 1, ZERO,
     $                     Y( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', M-I+1, I-1, ONE,
     $                     A( I, 1 ), LDA, A( I, I ), 1, ZERO,
     $                     Y( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ),
     $                     LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', M-I+1, I-1, ONE,
     $                     X( I, 1 ), LDX, A( I, I ), 1, ZERO,
     $                     Y( 1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', I-1, N-I, -ONE,
     $                     A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE,
     $                     Y( I+1, I ), 1 )
               CALL ZSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )
*
*              Update A(i,i+1:n)
*
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
               CALL ZLACGV( I, A( I, 1 ), LDA )
               CALL ZGEMV( 'No transpose', N-I, I, -ONE, Y( I+1, 1 ),
     $                     LDY, A( I, 1 ), LDA, ONE, A( I, I+1 ), LDA )
               CALL ZLACGV( I, A( I, 1 ), LDA )
               CALL ZLACGV( I-1, X( I, 1 ), LDX )
               CALL ZGEMV( 'Conjugate transpose', I-1, N-I, -ONE,
     $                     A( 1, I+1 ), LDA, X( I, 1 ), LDX, ONE,
     $                     A( I, I+1 ), LDA )
               CALL ZLACGV( I-1, X( I, 1 ), LDX )
*
*              Generate reflection P(i) to annihilate A(i,i+2:n)
*
               ALPHA = A( I, I+1 )
               CALL ZLARFG( N-I, ALPHA, A( I, MIN( I+2, N ) ), LDA,
     $                      TAUP( I ) )
               E( I ) = ALPHA
               A( I, I+1 ) = ONE
*
*              Compute X(i+1:m,i)
*
               CALL ZGEMV( 'No transpose', M-I, N-I, ONE, A( I+1, I+1 ),
     $                     LDA, A( I, I+1 ), LDA, ZERO, X( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I, I, ONE,
     $                     Y( I+1, 1 ), LDY, A( I, I+1 ), LDA, ZERO,
     $                     X( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', M-I, I, -ONE, A( I+1, 1 ),
     $                     LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL ZGEMV( 'No transpose', I-1, N-I, ONE, A( 1, I+1 ),
     $                     LDA, A( I, I+1 ), LDA, ZERO, X( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ),
     $                     LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL ZSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )
               CALL ZLACGV( N-I, A( I, I+1 ), LDA )
            END IF
   10    CONTINUE
      ELSE
*
*        Reduce to lower bidiagonal form
*
         DO 20 I = 1, NB
*
*           Update A(i,i:n)
*
            CALL ZLACGV( N-I+1, A( I, I ), LDA )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            CALL ZGEMV( 'No transpose', N-I+1, I-1, -ONE, Y( I, 1 ),
     $                  LDY, A( I, 1 ), LDA, ONE, A( I, I ), LDA )
            CALL ZLACGV( I-1, A( I, 1 ), LDA )
            CALL ZLACGV( I-1, X( I, 1 ), LDX )
            CALL ZGEMV( 'Conjugate transpose', I-1, N-I+1, -ONE,
     $                  A( 1, I ), LDA, X( I, 1 ), LDX, ONE, A( I, I ),
     $                  LDA )
            CALL ZLACGV( I-1, X( I, 1 ), LDX )
*
*           Generate reflection P(i) to annihilate A(i,i+1:n)
*
            ALPHA = A( I, I )
            CALL ZLARFG( N-I+1, ALPHA, A( I, MIN( I+1, N ) ), LDA,
     $                   TAUP( I ) )
            D( I ) = ALPHA
            IF( I.LT.M ) THEN
               A( I, I ) = ONE
*
*              Compute X(i+1:m,i)
*
               CALL ZGEMV( 'No transpose', M-I, N-I+1, ONE, A( I+1, I ),
     $                     LDA, A( I, I ), LDA, ZERO, X( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', N-I+1, I-1, ONE,
     $                     Y( I, 1 ), LDY, A( I, I ), LDA, ZERO,
     $                     X( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL ZGEMV( 'No transpose', I-1, N-I+1, ONE, A( 1, I ),
     $                     LDA, A( I, I ), LDA, ZERO, X( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', M-I, I-1, -ONE, X( I+1, 1 ),
     $                     LDX, X( 1, I ), 1, ONE, X( I+1, I ), 1 )
               CALL ZSCAL( M-I, TAUP( I ), X( I+1, I ), 1 )
               CALL ZLACGV( N-I+1, A( I, I ), LDA )
*
*              Update A(i+1:m,i)
*
               CALL ZLACGV( I-1, Y( I, 1 ), LDY )
               CALL ZGEMV( 'No transpose', M-I, I-1, -ONE, A( I+1, 1 ),
     $                     LDA, Y( I, 1 ), LDY, ONE, A( I+1, I ), 1 )
               CALL ZLACGV( I-1, Y( I, 1 ), LDY )
               CALL ZGEMV( 'No transpose', M-I, I, -ONE, X( I+1, 1 ),
     $                     LDX, A( 1, I ), 1, ONE, A( I+1, I ), 1 )
*
*              Generate reflection Q(i) to annihilate A(i+2:m,i)
*
               ALPHA = A( I+1, I )
               CALL ZLARFG( M-I, ALPHA, A( MIN( I+2, M ), I ), 1,
     $                      TAUQ( I ) )
               E( I ) = ALPHA
               A( I+1, I ) = ONE
*
*              Compute Y(i+1:n,i)
*
               CALL ZGEMV( 'Conjugate transpose', M-I, N-I, ONE,
     $                     A( I+1, I+1 ), LDA, A( I+1, I ), 1, ZERO,
     $                     Y( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', M-I, I-1, ONE,
     $                     A( I+1, 1 ), LDA, A( I+1, I ), 1, ZERO,
     $                     Y( 1, I ), 1 )
               CALL ZGEMV( 'No transpose', N-I, I-1, -ONE, Y( I+1, 1 ),
     $                     LDY, Y( 1, I ), 1, ONE, Y( I+1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', M-I, I, ONE,
     $                     X( I+1, 1 ), LDX, A( I+1, I ), 1, ZERO,
     $                     Y( 1, I ), 1 )
               CALL ZGEMV( 'Conjugate transpose', I, N-I, -ONE,
     $                     A( 1, I+1 ), LDA, Y( 1, I ), 1, ONE,
     $                     Y( I+1, I ), 1 )
               CALL ZSCAL( N-I, TAUQ( I ), Y( I+1, I ), 1 )
            ELSE
               CALL ZLACGV( N-I+1, A( I, I ), LDA )
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZLABRD
*
      END







      SUBROUTINE ZUNGBR( VECT, M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          VECT
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNGBR generates one of the complex unitary matrices Q or P**H
*  determined by ZGEBRD when reducing a complex matrix A to bidiagonal
*  form: A = Q * B * P**H.  Q and P**H are defined as products of
*  elementary reflectors H(i) or G(i) respectively.
*
*  If VECT = 'Q', A is assumed to have been an M-by-K matrix, and Q
*  is of order M:
*  if m >= k, Q = H(1) H(2) . . . H(k) and ZUNGBR returns the first n
*  columns of Q, where m >= n >= k;
*  if m < k, Q = H(1) H(2) . . . H(m-1) and ZUNGBR returns Q as an
*  M-by-M matrix.
*
*  If VECT = 'P', A is assumed to have been a K-by-N matrix, and P**H
*  is of order N:
*  if k < n, P**H = G(k) . . . G(2) G(1) and ZUNGBR returns the first m
*  rows of P**H, where n >= m >= k;
*  if k >= n, P**H = G(n-1) . . . G(2) G(1) and ZUNGBR returns P**H as
*  an N-by-N matrix.
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          Specifies whether the matrix Q or the matrix P**H is
*          required, as defined in the transformation applied by ZGEBRD:
*          = 'Q':  generate Q;
*          = 'P':  generate P**H.
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q or P**H to be returned.
*          M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q or P**H to be returned.
*          N >= 0.
*          If VECT = 'Q', M >= N >= min(M,K);
*          if VECT = 'P', N >= M >= min(N,K).
*
*  K       (input) INTEGER
*          If VECT = 'Q', the number of columns in the original M-by-K
*          matrix reduced by ZGEBRD.
*          If VECT = 'P', the number of rows in the original K-by-N
*          matrix reduced by ZGEBRD.
*          K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the vectors which define the elementary reflectors,
*          as returned by ZGEBRD.
*          On exit, the M-by-N matrix Q or P**H.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= M.
*
*  TAU     (input) COMPLEX*16 array, dimension
*                                (min(M,K)) if VECT = 'Q'
*                                (min(N,K)) if VECT = 'P'
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i) or G(i), which determines Q or P**H, as
*          returned by ZGEBRD in its array argument TAUQ or TAUP.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,min(M,N)).
*          For optimum performance LWORK >= min(M,N)*NB, where NB
*          is the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY, WANTQ
      INTEGER            I, IINFO, J, LWKOPT, MN, NB
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNGLQ, ZUNGQR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      WANTQ = LSAME( VECT, 'Q' )
      MN = MIN( M, N )
      LQUERY = ( LWORK.EQ.-1 )
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
         INFO = -1
      ELSE IF( M.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 .OR. ( WANTQ .AND. ( N.GT.M .OR. N.LT.MIN( M,
     $         K ) ) ) .OR. ( .NOT.WANTQ .AND. ( M.GT.N .OR. M.LT.
     $         MIN( N, K ) ) ) ) THEN
         INFO = -3
      ELSE IF( K.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -6
      ELSE IF( LWORK.LT.MAX( 1, MN ) .AND. .NOT.LQUERY ) THEN
         INFO = -9
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( WANTQ ) THEN
            NB = ILAENV( 1, 'ZUNGQR', ' ', M, N, K, -1 )
         ELSE
            NB = ILAENV( 1, 'ZUNGLQ', ' ', M, N, K, -1 )
         END IF
         LWKOPT = MAX( 1, MN )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGBR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      IF( WANTQ ) THEN
*
*        Form Q, determined by a call to ZGEBRD to reduce an m-by-k
*        matrix
*
         IF( M.GE.K ) THEN
*
*           If m >= k, assume m >= n >= k
*
            CALL ZUNGQR( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
*
         ELSE
*
*           If m < k, assume m = n
*
*           Shift the vectors which define the elementary reflectors one
*           column to the right, and set the first row and column of Q
*           to those of the unit matrix
*
            DO 20 J = M, 2, -1
               A( 1, J ) = ZERO
               DO 10 I = J + 1, M
                  A( I, J ) = A( I, J-1 )
   10          CONTINUE
   20       CONTINUE
            A( 1, 1 ) = ONE
            DO 30 I = 2, M
               A( I, 1 ) = ZERO
   30       CONTINUE
            IF( M.GT.1 ) THEN
*
*              Form Q(2:m,2:m)
*
               CALL ZUNGQR( M-1, M-1, M-1, A( 2, 2 ), LDA, TAU, WORK,
     $                      LWORK, IINFO )
            END IF
         END IF
      ELSE
*
*        Form P', determined by a call to ZGEBRD to reduce a k-by-n
*        matrix
*
         IF( K.LT.N ) THEN
*
*           If k < n, assume k <= m <= n
*
            CALL ZUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, IINFO )
*
         ELSE
*
*           If k >= n, assume m = n
*
*           Shift the vectors which define the elementary reflectors one
*           row downward, and set the first row and column of P' to
*           those of the unit matrix
*
            A( 1, 1 ) = ONE
            DO 40 I = 2, N
               A( I, 1 ) = ZERO
   40       CONTINUE
            DO 60 J = 2, N
               DO 50 I = J - 1, 2, -1
                  A( I, J ) = A( I-1, J )
   50          CONTINUE
               A( 1, J ) = ZERO
   60       CONTINUE
            IF( N.GT.1 ) THEN
*
*              Form P'(2:n,2:n)
*
               CALL ZUNGLQ( N-1, N-1, N-1, A( 2, 2 ), LDA, TAU, WORK,
     $                      LWORK, IINFO )
            END IF
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of ZUNGBR
*
      END
      SUBROUTINE ZUNGL2( M, N, K, A, LDA, TAU, WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNGL2 generates an m-by-n complex matrix Q with orthonormal rows,
*  which is defined as the first m rows of a product of k elementary
*  reflectors of order n
*
*        Q  =  H(k)' . . . H(2)' H(1)'
*
*  as returned by ZGELQF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. N >= M.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. M >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the i-th row must contain the vector which defines
*          the elementary reflector H(i), for i = 1,2,...,k, as returned
*          by ZGELQF in the first k rows of its array argument A.
*          On exit, the m by n matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGELQF.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (M)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLARF, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGL2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 )
     $   RETURN
*
      IF( K.LT.M ) THEN
*
*        Initialise rows k+1:m to rows of the unit matrix
*
         DO 20 J = 1, N
            DO 10 L = K + 1, M
               A( L, J ) = ZERO
   10       CONTINUE
            IF( J.GT.K .AND. J.LE.M )
     $         A( J, J ) = ONE
   20    CONTINUE
      END IF
*
      DO 40 I = K, 1, -1
*
*        Apply H(i)' to A(i:m,i:n) from the right
*
         IF( I.LT.N ) THEN
            CALL ZLACGV( N-I, A( I, I+1 ), LDA )
            IF( I.LT.M ) THEN
               A( I, I ) = ONE
               CALL ZLARF( 'Right', M-I, N-I+1, A( I, I ), LDA,
     $                     DCONJG( TAU( I ) ), A( I+1, I ), LDA, WORK )
            END IF
            CALL ZSCAL( N-I, -TAU( I ), A( I, I+1 ), LDA )
            CALL ZLACGV( N-I, A( I, I+1 ), LDA )
         END IF
         A( I, I ) = ONE - DCONJG( TAU( I ) )
*
*        Set A(i,1:i-1) to zero
*
         DO 30 L = 1, I - 1
            A( I, L ) = ZERO
   30    CONTINUE
   40 CONTINUE
      RETURN
*
*     End of ZUNGL2
*
      END
      SUBROUTINE ZUNGLQ( M, N, K, A, LDA, TAU, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      INTEGER            INFO, K, LDA, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNGLQ generates an M-by-N complex matrix Q with orthonormal rows,
*  which is defined as the first M rows of a product of K elementary
*  reflectors of order N
*
*        Q  =  H(k)' . . . H(2)' H(1)'
*
*  as returned by ZGELQF.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix Q. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix Q. N >= M.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. M >= K >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the i-th row must contain the vector which defines
*          the elementary reflector H(i), for i = 1,2,...,k, as returned
*          by ZGELQF in the first k rows of its array argument A.
*          On exit, the M-by-N matrix Q.
*
*  LDA     (input) INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGELQF.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,M).
*          For optimum performance LWORK >= M*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit;
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IB, IINFO, IWS, J, KI, KK, L, LDWORK,
     $                   LWKOPT, NB, NBMIN, NX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNGL2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      NB = ILAENV( 1, 'ZUNGLQ', ' ', M, N, K, -1 )
      LWKOPT = MAX( 1, M )*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.M ) THEN
         INFO = -2
      ELSE IF( K.LT.0 .OR. K.GT.M ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LWORK.LT.MAX( 1, M ) .AND. .NOT.LQUERY ) THEN
         INFO = -8
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNGLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.LE.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      NX = 0
      IWS = M
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
*
*        Determine when to cross over from blocked to unblocked code.
*
         NX = MAX( 0, ILAENV( 3, 'ZUNGLQ', ' ', M, N, K, -1 ) )
         IF( NX.LT.K ) THEN
*
*           Determine if workspace is large enough for blocked code.
*
            LDWORK = M
            IWS = LDWORK*NB
            IF( LWORK.LT.IWS ) THEN
*
*              Not enough workspace to use optimal NB:  reduce NB and
*              determine the minimum value of NB.
*
               NB = LWORK / LDWORK
               NBMIN = MAX( 2, ILAENV( 2, 'ZUNGLQ', ' ', M, N, K, -1 ) )
            END IF
         END IF
      END IF
*
      IF( NB.GE.NBMIN .AND. NB.LT.K .AND. NX.LT.K ) THEN
*
*        Use blocked code after the last block.
*        The first kk rows are handled by the block method.
*
         KI = ( ( K-NX-1 ) / NB )*NB
         KK = MIN( K, KI+NB )
*
*        Set A(kk+1:m,1:kk) to zero.
*
         DO 20 J = 1, KK
            DO 10 I = KK + 1, M
               A( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      ELSE
         KK = 0
      END IF
*
*     Use unblocked code for the last or only block.
*
      IF( KK.LT.M )
     $   CALL ZUNGL2( M-KK, N-KK, K-KK, A( KK+1, KK+1 ), LDA,
     $                TAU( KK+1 ), WORK, IINFO )
*
      IF( KK.GT.0 ) THEN
*
*        Use blocked code
*
         DO 50 I = KI + 1, 1, -NB
            IB = MIN( NB, K-I+1 )
            IF( I+IB.LE.M ) THEN
*
*              Form the triangular factor of the block reflector
*              H = H(i) H(i+1) . . . H(i+ib-1)
*
               CALL ZLARFT( 'Forward', 'Rowwise', N-I+1, IB, A( I, I ),
     $                      LDA, TAU( I ), WORK, LDWORK )
*
*              Apply H' to A(i+ib:m,i:n) from the right
*
               CALL ZLARFB( 'Right', 'Conjugate transpose', 'Forward',
     $                      'Rowwise', M-I-IB+1, N-I+1, IB, A( I, I ),
     $                      LDA, WORK, LDWORK, A( I+IB, I ), LDA,
     $                      WORK( IB+1 ), LDWORK )
            END IF
*
*           Apply H' to columns i:n of current block
*
            CALL ZUNGL2( IB, N-I+1, IB, A( I, I ), LDA, TAU( I ), WORK,
     $                   IINFO )
*
*           Set columns 1:i-1 of current block to zero
*
            DO 40 J = 1, I - 1
               DO 30 L = I, I + IB - 1
                  A( L, J ) = ZERO
   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
      END IF
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZUNGLQ
*
      END

      SUBROUTINE ZUNMBR( VECT, SIDE, TRANS, M, N, K, A, LDA, TAU, C,
     $                   LDC, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS, VECT
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  If VECT = 'Q', ZUNMBR overwrites the general complex M-by-N matrix C
*  with
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'C':      Q**H * C       C * Q**H
*
*  If VECT = 'P', ZUNMBR overwrites the general complex M-by-N matrix C
*  with
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      P * C          C * P
*  TRANS = 'C':      P**H * C       C * P**H
*
*  Here Q and P**H are the unitary matrices determined by ZGEBRD when
*  reducing a complex matrix A to bidiagonal form: A = Q * B * P**H. Q
*  and P**H are defined as products of elementary reflectors H(i) and
*  G(i) respectively.
*
*  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
*  order of the unitary matrix Q or P**H that is applied.
*
*  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
*  if nq >= k, Q = H(1) H(2) . . . H(k);
*  if nq < k, Q = H(1) H(2) . . . H(nq-1).
*
*  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
*  if k < nq, P = G(1) G(2) . . . G(k);
*  if k >= nq, P = G(1) G(2) . . . G(nq-1).
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'Q': apply Q or Q**H;
*          = 'P': apply P or P**H.
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q, Q**H, P or P**H from the Left;
*          = 'R': apply Q, Q**H, P or P**H from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q or P;
*          = 'C':  Conjugate transpose, apply Q**H or P**H.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          If VECT = 'Q', the number of columns in the original
*          matrix reduced by ZGEBRD.
*          If VECT = 'P', the number of rows in the original
*          matrix reduced by ZGEBRD.
*          K >= 0.
*
*  A       (input) COMPLEX*16 array, dimension
*                                (LDA,min(nq,K)) if VECT = 'Q'
*                                (LDA,nq)        if VECT = 'P'
*          The vectors which define the elementary reflectors H(i) and
*          G(i), whose products determine the matrices Q and P, as
*          returned by ZGEBRD.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If VECT = 'Q', LDA >= max(1,nq);
*          if VECT = 'P', LDA >= max(1,min(nq,K)).
*
*  TAU     (input) COMPLEX*16 array, dimension (min(nq,K))
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i) or G(i) which determines Q or P, as returned
*          by ZGEBRD in the array argument TAUQ or TAUP.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q
*          or P*C or P**H*C or C*P or C*P**H.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Local Scalars ..
      LOGICAL            APPLYQ, LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I1, I2, IINFO, LWKOPT, MI, NB, NI, NQ, NW
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZUNMLQ, ZUNMQR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      APPLYQ = LSAME( VECT, 'Q' )
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q or P and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.APPLYQ .AND. .NOT.LSAME( VECT, 'P' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -2
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -3
      ELSE IF( M.LT.0 ) THEN
         INFO = -4
      ELSE IF( N.LT.0 ) THEN
         INFO = -5
      ELSE IF( K.LT.0 ) THEN
         INFO = -6
      ELSE IF( ( APPLYQ .AND. LDA.LT.MAX( 1, NQ ) ) .OR.
     $         ( .NOT.APPLYQ .AND. LDA.LT.MAX( 1, MIN( NQ, K ) ) ) )
     $          THEN
         INFO = -8
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -13
      END IF
*
      IF( INFO.EQ.0 ) THEN
         IF( APPLYQ ) THEN
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M-1, N, M-1,
     $              -1 )
            ELSE
               NB = ILAENV( 1, 'ZUNMQR', SIDE // TRANS, M, N-1, N-1,
     $              -1 )
            END IF
         ELSE
            IF( LEFT ) THEN
               NB = ILAENV( 1, 'ZUNMLQ', SIDE // TRANS, M-1, N, M-1,
     $              -1 )
            ELSE
               NB = ILAENV( 1, 'ZUNMLQ', SIDE // TRANS, M, N-1, N-1,
     $              -1 )
            END IF
         END IF
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMBR', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
      END IF
*
*     Quick return if possible
*
      WORK( 1 ) = 1
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
      IF( APPLYQ ) THEN
*
*        Apply Q
*
         IF( NQ.GE.K ) THEN
*
*           Q was determined by a call to ZGEBRD with nq >= k
*
            CALL ZUNMQR( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN
*
*           Q was determined by a call to ZGEBRD with nq < k
*
            IF( LEFT ) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL ZUNMQR( SIDE, TRANS, MI, NI, NQ-1, A( 2, 1 ), LDA, TAU,
     $                   C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         END IF
      ELSE
*
*        Apply P
*
         IF( NOTRAN ) THEN
            TRANST = 'C'
         ELSE
            TRANST = 'N'
         END IF
         IF( NQ.GT.K ) THEN
*
*           P was determined by a call to ZGEBRD with nq > k
*
            CALL ZUNMLQ( SIDE, TRANST, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, IINFO )
         ELSE IF( NQ.GT.1 ) THEN
*
*           P was determined by a call to ZGEBRD with nq <= k
*
            IF( LEFT ) THEN
               MI = M - 1
               NI = N
               I1 = 2
               I2 = 1
            ELSE
               MI = M
               NI = N - 1
               I1 = 1
               I2 = 2
            END IF
            CALL ZUNMLQ( SIDE, TRANST, MI, NI, NQ-1, A( 1, 2 ), LDA,
     $                   TAU, C( I1, I2 ), LDC, WORK, LWORK, IINFO )
         END IF
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of ZUNMBR
*
      END
      SUBROUTINE ZUNML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNML2 overwrites the general complex m-by-n matrix C with
*
*        Q * C  if SIDE = 'L' and TRANS = 'N', or
*
*        Q'* C  if SIDE = 'L' and TRANS = 'C', or
*
*        C * Q  if SIDE = 'R' and TRANS = 'N', or
*
*        C * Q' if SIDE = 'R' and TRANS = 'C',
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(k)' . . . H(2)' H(1)'
*
*  as returned by ZGELQF. Q is of order m if SIDE = 'L' and of order n
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q' from the Left
*          = 'R': apply Q or Q' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply Q  (No transpose)
*          = 'C': apply Q' (Conjugate transpose)
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) COMPLEX*16 array, dimension
*                               (LDA,M) if SIDE = 'L',
*                               (LDA,N) if SIDE = 'R'
*          The i-th row must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          ZGELQF in the first k rows of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,K).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGELQF.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the m-by-n matrix C.
*          On exit, C is overwritten by Q*C or Q'*C or C*Q' or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) COMPLEX*16 array, dimension
*                                   (N) if SIDE = 'L',
*                                   (M) if SIDE = 'R'
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, NOTRAN
      INTEGER            I, I1, I2, I3, IC, JC, MI, NI, NQ
      COMPLEX*16         AII, TAUI
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLARF
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG, MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
*
*     NQ is the order of Q
*
      IF( LEFT ) THEN
         NQ = M
      ELSE
         NQ = N
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNML2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 )
     $   RETURN
*
      IF( ( LEFT .AND. NOTRAN .OR. .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
         I1 = 1
         I2 = K
         I3 = 1
      ELSE
         I1 = K
         I2 = 1
         I3 = -1
      END IF
*
      IF( LEFT ) THEN
         NI = N
         JC = 1
      ELSE
         MI = M
         IC = 1
      END IF
*
      DO 10 I = I1, I2, I3
         IF( LEFT ) THEN
*
*           H(i) or H(i)' is applied to C(i:m,1:n)
*
            MI = M - I + 1
            IC = I
         ELSE
*
*           H(i) or H(i)' is applied to C(1:m,i:n)
*
            NI = N - I + 1
            JC = I
         END IF
*
*        Apply H(i) or H(i)'
*
         IF( NOTRAN ) THEN
            TAUI = DCONJG( TAU( I ) )
         ELSE
            TAUI = TAU( I )
         END IF
         IF( I.LT.NQ )
     $      CALL ZLACGV( NQ-I, A( I, I+1 ), LDA )
         AII = A( I, I )
         A( I, I ) = ONE
         CALL ZLARF( SIDE, MI, NI, A( I, I ), LDA, TAUI, C( IC, JC ),
     $               LDC, WORK )
         A( I, I ) = AII
         IF( I.LT.NQ )
     $      CALL ZLACGV( NQ-I, A( I, I+1 ), LDA )
   10 CONTINUE
      RETURN
*
*     End of ZUNML2
*
      END
      SUBROUTINE ZUNMLQ( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC,
     $                   WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     June 30, 1999
*
*     .. Scalar Arguments ..
      CHARACTER          SIDE, TRANS
      INTEGER            INFO, K, LDA, LDC, LWORK, M, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), C( LDC, * ), TAU( * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZUNMLQ overwrites the general complex M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'C':      Q**H * C       C * Q**H
*
*  where Q is a complex unitary matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(k)' . . . H(2)' H(1)'
*
*  as returned by ZGELQF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**H from the Left;
*          = 'R': apply Q or Q**H from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'C':  Conjugate transpose, apply Q**H.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) COMPLEX*16 array, dimension
*                               (LDA,M) if SIDE = 'L',
*                               (LDA,N) if SIDE = 'R'
*          The i-th row must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          ZGELQF in the first k rows of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,K).
*
*  TAU     (input) COMPLEX*16 array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by ZGELQF.
*
*  C       (input/output) COMPLEX*16 array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**H*C or C*Q**H or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      INTEGER            NBMAX, LDT
      PARAMETER          ( NBMAX = 64, LDT = NBMAX+1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFT, LQUERY, NOTRAN
      CHARACTER          TRANST
      INTEGER            I, I1, I2, I3, IB, IC, IINFO, IWS, JC, LDWORK,
     $                   LWKOPT, MI, NB, NBMIN, NI, NQ, NW
*     ..
*     .. Local Arrays ..
      COMPLEX*16         T( LDT, NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLARFB, ZLARFT, ZUNML2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input arguments
*
      INFO = 0
      LEFT = LSAME( SIDE, 'L' )
      NOTRAN = LSAME( TRANS, 'N' )
      LQUERY = ( LWORK.EQ.-1 )
*
*     NQ is the order of Q and NW is the minimum dimension of WORK
*
      IF( LEFT ) THEN
         NQ = M
         NW = N
      ELSE
         NQ = N
         NW = M
      END IF
      IF( .NOT.LEFT .AND. .NOT.LSAME( SIDE, 'R' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOTRAN .AND. .NOT.LSAME( TRANS, 'C' ) ) THEN
         INFO = -2
      ELSE IF( M.LT.0 ) THEN
         INFO = -3
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( K.LT.0 .OR. K.GT.NQ ) THEN
         INFO = -5
      ELSE IF( LDA.LT.MAX( 1, K ) ) THEN
         INFO = -7
      ELSE IF( LDC.LT.MAX( 1, M ) ) THEN
         INFO = -10
      ELSE IF( LWORK.LT.MAX( 1, NW ) .AND. .NOT.LQUERY ) THEN
         INFO = -12
      END IF
*
      IF( INFO.EQ.0 ) THEN
*
*        Determine the block size.  NB may be at most NBMAX, where NBMAX
*        is used to define the local array T.
*
         NB = MIN( NBMAX, ILAENV( 1, 'ZUNMLQ', SIDE // TRANS, M, N, K,
     $        -1 ) )
         LWKOPT = MAX( 1, NW )*NB
         WORK( 1 ) = LWKOPT
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZUNMLQ', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 .OR. K.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF
*
      NBMIN = 2
      LDWORK = NW
      IF( NB.GT.1 .AND. NB.LT.K ) THEN
         IWS = NW*NB
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZUNMLQ', SIDE // TRANS, M, N, K,
     $              -1 ) )
         END IF
      ELSE
         IWS = NW
      END IF
*
      IF( NB.LT.NBMIN .OR. NB.GE.K ) THEN
*
*        Use unblocked code
*
         CALL ZUNML2( SIDE, TRANS, M, N, K, A, LDA, TAU, C, LDC, WORK,
     $                IINFO )
      ELSE
*
*        Use blocked code
*
         IF( ( LEFT .AND. NOTRAN ) .OR.
     $       ( .NOT.LEFT .AND. .NOT.NOTRAN ) ) THEN
            I1 = 1
            I2 = K
            I3 = NB
         ELSE
            I1 = ( ( K-1 ) / NB )*NB + 1
            I2 = 1
            I3 = -NB
         END IF
*
         IF( LEFT ) THEN
            NI = N
            JC = 1
         ELSE
            MI = M
            IC = 1
         END IF
*
         IF( NOTRAN ) THEN
            TRANST = 'C'
         ELSE
            TRANST = 'N'
         END IF
*
         DO 10 I = I1, I2, I3
            IB = MIN( NB, K-I+1 )
*
*           Form the triangular factor of the block reflector
*           H = H(i) H(i+1) . . . H(i+ib-1)
*
            CALL ZLARFT( 'Forward', 'Rowwise', NQ-I+1, IB, A( I, I ),
     $                   LDA, TAU( I ), T, LDT )
            IF( LEFT ) THEN
*
*              H or H' is applied to C(i:m,1:n)
*
               MI = M - I + 1
               IC = I
            ELSE
*
*              H or H' is applied to C(1:m,i:n)
*
               NI = N - I + 1
               JC = I
            END IF
*
*           Apply H or H'
*
            CALL ZLARFB( SIDE, TRANST, 'Forward', 'Rowwise', MI, NI, IB,
     $                   A( I, I ), LDA, T, LDT, C( IC, JC ), LDC, WORK,
     $                   LDWORK )
   10    CONTINUE
      END IF
      WORK( 1 ) = LWKOPT
      RETURN
*
*     End of ZUNMLQ
*
      END

      SUBROUTINE ZDROT( N, CX, INCX, CY, INCY, C, S )
*
*     .. Scalar Arguments ..
      INTEGER            INCX, INCY, N
      DOUBLE PRECISION   C, S
*     ..
*     .. Array Arguments ..
      COMPLEX*16         CX( * ), CY( * )
*     ..
*
*  Purpose
*  =======
*
*  Applies a plane rotation, where the cos and sin (c and s) are real
*  and the vectors cx and cy are complex.
*  jack dongarra, linpack, 3/11/78.
*
*  Arguments
*  ==========
*
*  N        (input) INTEGER
*           On entry, N specifies the order of the vectors cx and cy.
*           N must be at least zero.
*           Unchanged on exit.
*
*  CX       (input) COMPLEX*16 array, dimension at least
*           ( 1 + ( N - 1 )*abs( INCX ) ).
*           Before entry, the incremented array CX must contain the n
*           element vector cx. On exit, CX is overwritten by the updated
*           vector cx.
*
*  INCX     (input) INTEGER
*           On entry, INCX specifies the increment for the elements of
*           CX. INCX must not be zero.
*           Unchanged on exit.
*
*  CY       (input) COMPLEX*16 array, dimension at least
*           ( 1 + ( N - 1 )*abs( INCY ) ).
*           Before entry, the incremented array CY must contain the n
*           element vector cy. On exit, CY is overwritten by the updated
*           vector cy.
*
*  INCY     (input) INTEGER
*           On entry, INCY specifies the increment for the elements of
*           CY. INCY must not be zero.
*           Unchanged on exit.
*
*  C        (input) DOUBLE PRECISION
*           On entry, C specifies the cosine, cos.
*           Unchanged on exit.
*
*  S        (input) DOUBLE PRECISION
*           On entry, S specifies the sine, sin.
*           Unchanged on exit.
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IX, IY
      COMPLEX*16         CTEMP
*     ..
*     .. Executable Statements ..
*
      IF( N.LE.0 )
     $   RETURN
      IF( INCX.EQ.1 .AND. INCY.EQ.1 )
     $   GO TO 20
*
*        code for unequal increments or equal increments not equal
*          to 1
*
      IX = 1
      IY = 1
      IF( INCX.LT.0 )
     $   IX = ( -N+1 )*INCX + 1
      IF( INCY.LT.0 )
     $   IY = ( -N+1 )*INCY + 1
      DO 10 I = 1, N
         CTEMP = C*CX( IX ) + S*CY( IY )
         CY( IY ) = C*CY( IY ) - S*CX( IX )
         CX( IX ) = CTEMP
         IX = IX + INCX
         IY = IY + INCY
   10 CONTINUE
      RETURN
*
*        code for both increments equal to 1
*
   20 CONTINUE
      DO 30 I = 1, N
         CTEMP = C*CX( I ) + S*CY( I )
         CY( I ) = C*CY( I ) - S*CX( I )
         CX( I ) = CTEMP
   30 CONTINUE
      RETURN
      END
      SUBROUTINE ZGETRF( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRF computes an LU factorization of a general M-by-N matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the M-by-N matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      INTEGER            I, IINFO, J, JB, NB
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGETF2, ZLASWP, ZTRSM
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRF', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZGETRF', ' ', M, N, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.MIN( M, N ) ) THEN
*
*        Use unblocked code.
*
         CALL ZGETF2( M, N, A, LDA, IPIV, INFO )
      ELSE
*
*        Use blocked code.
*
         DO 20 J = 1, MIN( M, N ), NB
            JB = MIN( MIN( M, N )-J+1, NB )
*
*           Factor diagonal and subdiagonal blocks and test for exact
*           singularity.
*
            CALL ZGETF2( M-J+1, JB, A( J, J ), LDA, IPIV( J ), IINFO )
*
*           Adjust INFO and the pivot indices.
*
            IF( INFO.EQ.0 .AND. IINFO.GT.0 )
     $         INFO = IINFO + J - 1
            DO 10 I = J, MIN( M, J+JB-1 )
               IPIV( I ) = J - 1 + IPIV( I )
   10       CONTINUE
*
*           Apply interchanges to columns 1:J-1.
*
            CALL ZLASWP( J-1, A, LDA, J, J+JB-1, IPIV, 1 )
*
            IF( J+JB.LE.N ) THEN
*
*              Apply interchanges to columns J+JB:N.
*
               CALL ZLASWP( N-J-JB+1, A( 1, J+JB ), LDA, J, J+JB-1,
     $                      IPIV, 1 )
*
*              Compute block row of U.
*
               CALL ZTRSM( 'Left', 'Lower', 'No transpose', 'Unit', JB,
     $                     N-J-JB+1, ONE, A( J, J ), LDA, A( J, J+JB ),
     $                     LDA )
               IF( J+JB.LE.M ) THEN
*
*                 Update trailing submatrix.
*
                  CALL ZGEMM( 'No transpose', 'No transpose', M-J-JB+1,
     $                        N-J-JB+1, JB, -ONE, A( J+JB, J ), LDA,
     $                        A( J, J+JB ), LDA, ONE, A( J+JB, J+JB ),
     $                        LDA )
               END IF
            END IF
   20    CONTINUE
      END IF
      RETURN
*
*     End of ZGETRF
*
      END
      SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, LWORK, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETRI computes the inverse of a matrix using the LU factorization
*  computed by ZGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by ZGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from ZGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ZERO, ONE
      PARAMETER          ( ZERO = ( 0.0D+0, 0.0D+0 ),
     $                   ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            LQUERY
      INTEGER            I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB,
     $                   NBMIN, NN
*     ..
*     .. External Functions ..
      INTEGER            ILAENV
      EXTERNAL           ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGEMM, ZGEMV, ZSWAP, ZTRSM, ZTRTRI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      NB = ILAENV( 1, 'ZGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = N*NB
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( N.LT.0 ) THEN
         INFO = -1
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -3
      ELSE IF( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -6
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETRI', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Form inv(U).  If INFO > 0 from ZTRTRI, then U is singular,
*     and the inverse is not computed.
*
      CALL ZTRTRI( 'Upper', 'Non-unit', N, A, LDA, INFO )
      IF( INFO.GT.0 )
     $   RETURN
*
      NBMIN = 2
      LDWORK = N
      IF( NB.GT.1 .AND. NB.LT.N ) THEN
         IWS = MAX( LDWORK*NB, 1 )
         IF( LWORK.LT.IWS ) THEN
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZGETRI', ' ', N, -1, -1, -1 ) )
         END IF
      ELSE
         IWS = N
      END IF
*
*     Solve the equation inv(A)*L = inv(U) for inv(A).
*
      IF( NB.LT.NBMIN .OR. NB.GE.N ) THEN
*
*        Use unblocked code.
*
         DO 20 J = N, 1, -1
*
*           Copy current column of L to WORK and replace with zeros.
*
            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE
*
*           Compute current column of inv(A).
*
            IF( J.LT.N )
     $         CALL ZGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ),
     $                     LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      ELSE
*
*        Use blocked code.
*
         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )
*
*           Copy current block column of L to WORK and replace with
*           zeros.
*
            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE
*
*           Compute current block column of inv(A).
*
            IF( J+JB.LE.N )
     $         CALL ZGEMM( 'No transpose', 'No transpose', N, JB,
     $                     N-J-JB+1, -ONE, A( 1, J+JB ), LDA,
     $                     WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            CALL ZTRSM( 'Right', 'Lower', 'No transpose', 'Unit', N, JB,
     $                  ONE, WORK( J ), LDWORK, A( 1, J ), LDA )
   50    CONTINUE
      END IF
*
*     Apply column interchanges.
*
      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J )
     $      CALL ZSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE
*
      WORK( 1 ) = IWS
      RETURN
*
*     End of ZGETRI
*
      END
      SUBROUTINE ZHBEV( JOBZ, UPLO, N, KD, AB, LDAB, W, Z, LDZ, WORK,
     $                  RWORK, INFO )
*
*  -- LAPACK driver routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, KD, LDAB, LDZ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   RWORK( * ), W( * )
      COMPLEX*16         AB( LDAB, * ), WORK( * ), Z( LDZ, * )
*     ..
*
*  Purpose
*  =======
*
*  ZHBEV computes all the eigenvalues and, optionally, eigenvectors of
*  a complex Hermitian band matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) COMPLEX*16 array, dimension (LDAB, N)
*          On entry, the upper or lower triangle of the Hermitian band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, AB is overwritten by values generated during the
*          reduction to tridiagonal form.  If UPLO = 'U', the first
*          superdiagonal and the diagonal of the tridiagonal matrix T
*          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
*          the diagonal and first subdiagonal of T are returned in the
*          first two rows of AB.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD + 1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1,3*N-2))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE
      PARAMETER          ( ZERO = 0.0D0, ONE = 1.0D0 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LOWER, WANTZ
      INTEGER            IINFO, IMAX, INDE, INDRWK, ISCALE
      DOUBLE PRECISION   ANRM, BIGNUM, EPS, RMAX, RMIN, SAFMIN, SIGMA,
     $                   SMLNUM
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      DOUBLE PRECISION   DLAMCH, ZLANHB
      EXTERNAL           LSAME, DLAMCH, ZLANHB
*     ..
*     .. External Subroutines ..
      EXTERNAL           DSCAL, DSTERF, XERBLA, ZHBTRD, ZLASCL, ZSTEQR
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          SQRT
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      WANTZ = LSAME( JOBZ, 'V' )
      LOWER = LSAME( UPLO, 'L' )
*
      INFO = 0
      IF( .NOT.( WANTZ .OR. LSAME( JOBZ, 'N' ) ) ) THEN
         INFO = -1
      ELSE IF( .NOT.( LOWER .OR. LSAME( UPLO, 'U' ) ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KD.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD+1 ) THEN
         INFO = -6
      ELSE IF( LDZ.LT.1 .OR. ( WANTZ .AND. LDZ.LT.N ) ) THEN
         INFO = -9
      END IF
*
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHBEV ', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
      IF( N.EQ.1 ) THEN
         IF( LOWER ) THEN
            W( 1 ) = AB( 1, 1 )
         ELSE
            W( 1 ) = AB( KD+1, 1 )
         END IF
         IF( WANTZ )
     $      Z( 1, 1 ) = ONE
         RETURN
      END IF
*
*     Get machine constants.
*
      SAFMIN = DLAMCH( 'Safe minimum' )
      EPS = DLAMCH( 'Precision' )
      SMLNUM = SAFMIN / EPS
      BIGNUM = ONE / SMLNUM
      RMIN = SQRT( SMLNUM )
      RMAX = SQRT( BIGNUM )
*
*     Scale matrix to allowable range, if necessary.
*
      ANRM = ZLANHB( 'M', UPLO, N, KD, AB, LDAB, RWORK )
      ISCALE = 0
      IF( ANRM.GT.ZERO .AND. ANRM.LT.RMIN ) THEN
         ISCALE = 1
         SIGMA = RMIN / ANRM
      ELSE IF( ANRM.GT.RMAX ) THEN
         ISCALE = 1
         SIGMA = RMAX / ANRM
      END IF
      IF( ISCALE.EQ.1 ) THEN
         IF( LOWER ) THEN
            CALL ZLASCL( 'B', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         ELSE
            CALL ZLASCL( 'Q', KD, KD, ONE, SIGMA, N, N, AB, LDAB, INFO )
         END IF
      END IF
*
*     Call ZHBTRD to reduce Hermitian band matrix to tridiagonal form.
*
      INDE = 1
      CALL ZHBTRD( JOBZ, UPLO, N, KD, AB, LDAB, W, RWORK( INDE ), Z,
     $             LDZ, WORK, IINFO )
*
*     For eigenvalues only, call DSTERF.  For eigenvectors, call ZSTEQR.
*
      IF( .NOT.WANTZ ) THEN
         CALL DSTERF( N, W, RWORK( INDE ), INFO )
      ELSE
         INDRWK = INDE + N
         CALL ZSTEQR( JOBZ, N, W, RWORK( INDE ), Z, LDZ,
     $                RWORK( INDRWK ), INFO )
      END IF
*
*     If matrix was scaled, then rescale eigenvalues appropriately.
*
      IF( ISCALE.EQ.1 ) THEN
         IF( INFO.EQ.0 ) THEN
            IMAX = N
         ELSE
            IMAX = INFO - 1
         END IF
         CALL DSCAL( IMAX, ONE / SIGMA, W, 1 )
      END IF
*
      RETURN
*
*     End of ZHBEV
*
      END
      SUBROUTINE ZGETF2( M, N, A, LDA, IPIV, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INFO, LDA, M, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZGETF2 computes an LU factorization of a general m-by-n matrix A
*  using partial pivoting with row interchanges.
*
*  The factorization has the form
*     A = P * L * U
*  where P is a permutation matrix, L is lower triangular with unit
*  diagonal elements (lower trapezoidal if m > n), and U is upper
*  triangular (upper trapezoidal if m < n).
*
*  This is the right-looking Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the m by n matrix to be factored.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, U(k,k) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      DOUBLE PRECISION   SFMIN
      INTEGER            I, J, JP
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH
      INTEGER            IZAMAX
      EXTERNAL           DLAMCH, IZAMAX
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZGERU, ZSCAL, ZSWAP
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGETF2', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( M.EQ.0 .OR. N.EQ.0 )
     $   RETURN
*
*     Compute machine safe minimum
*
      SFMIN = DLAMCH('S') 
*
      DO 10 J = 1, MIN( M, N )
*
*        Find pivot and test for singularity.
*
         JP = J - 1 + IZAMAX( M-J+1, A( J, J ), 1 )
         IPIV( J ) = JP
         IF( A( JP, J ).NE.ZERO ) THEN
*
*           Apply the interchange to columns 1:N.
*
            IF( JP.NE.J )
     $         CALL ZSWAP( N, A( J, 1 ), LDA, A( JP, 1 ), LDA )
*
*           Compute elements J+1:M of J-th column.
*
            IF( J.LT.M ) THEN
               IF( ABS(A( J, J )) .GE. SFMIN ) THEN
                  CALL ZSCAL( M-J, ONE / A( J, J ), A( J+1, J ), 1 )
               ELSE
                  DO 20 I = 1, M-J
                     A( J+I, J ) = A( J+I, J ) / A( J, J )
   20             CONTINUE
               END IF
            END IF
*
         ELSE IF( INFO.EQ.0 ) THEN
*
            INFO = J
         END IF
*
         IF( J.LT.MIN( M, N ) ) THEN
*
*           Update trailing submatrix.
*
            CALL ZGERU( M-J, N-J, -ONE, A( J+1, J ), 1, A( J, J+1 ),
     $                  LDA, A( J+1, J+1 ), LDA )
         END IF
   10 CONTINUE
      RETURN
*
*     End of ZGETF2
*
      END
      SUBROUTINE ZHBTRD( VECT, UPLO, N, KD, AB, LDAB, D, E, Q, LDQ,
     $                   WORK, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO, VECT
      INTEGER            INFO, KD, LDAB, LDQ, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   D( * ), E( * )
      COMPLEX*16         AB( LDAB, * ), Q( LDQ, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  ZHBTRD reduces a complex Hermitian band matrix A to real symmetric
*  tridiagonal form T by a unitary similarity transformation:
*  Q**H * A * Q = T.
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'N':  do not form Q;
*          = 'V':  form Q;
*          = 'U':  update a matrix X, by forming X*Q.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) COMPLEX*16 array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the Hermitian band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*          On exit, the diagonal elements of AB are overwritten by the
*          diagonal elements of the tridiagonal matrix T; if KD > 0, the
*          elements on the first superdiagonal (if UPLO = 'U') or the
*          first subdiagonal (if UPLO = 'L') are overwritten by the
*          off-diagonal elements of T; the rest of AB is overwritten by
*          values generated during the reduction.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T.
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
*
*  Q       (input/output) COMPLEX*16 array, dimension (LDQ,N)
*          On entry, if VECT = 'U', then Q must contain an N-by-N
*          matrix X; if VECT = 'N' or 'V', then Q need not be set.
*
*          On exit:
*          if VECT = 'V', Q contains the N-by-N unitary matrix Q;
*          if VECT = 'U', Q contains the product X*Q;
*          if VECT = 'N', the array Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  Modified by Linda Kaufman, Bell Labs.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      COMPLEX*16         CZERO, CONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                   CONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            INITQ, UPPER, WANTQ
      INTEGER            I, I2, IBL, INCA, INCX, IQAEND, IQB, IQEND, J,
     $                   J1, J1END, J1INC, J2, JEND, JIN, JINC, K, KD1,
     $                   KDM1, KDN, L, LAST, LEND, NQ, NR, NRT
      DOUBLE PRECISION   ABST
      COMPLEX*16         T, TEMP
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZLACGV, ZLAR2V, ZLARGV, ZLARTG, ZLARTV,
     $                   ZLASET, ZROT, ZSCAL
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCONJG, MAX, MIN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters
*
      INITQ = LSAME( VECT, 'V' )
      WANTQ = INITQ .OR. LSAME( VECT, 'U' )
      UPPER = LSAME( UPLO, 'U' )
      KD1 = KD + 1
      KDM1 = KD - 1
      INCX = LDAB - 1
      IQEND = 1
*
      INFO = 0
      IF( .NOT.WANTQ .AND. .NOT.LSAME( VECT, 'N' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( KD.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDAB.LT.KD1 ) THEN
         INFO = -6
      ELSE IF( LDQ.LT.MAX( 1, N ) .AND. WANTQ ) THEN
         INFO = -10
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZHBTRD', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Initialize Q to the unit matrix, if needed
*
      IF( INITQ )
     $   CALL ZLASET( 'Full', N, N, CZERO, CONE, Q, LDQ )
*
*     Wherever possible, plane rotations are generated and applied in
*     vector operations of length NR over the index set J1:J2:KD1.
*
*     The real cosines and complex sines of the plane rotations are
*     stored in the arrays D and WORK.
*
      INCA = KD1*LDAB
      KDN = MIN( N-1, KD )
      IF( UPPER ) THEN
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to complex Hermitian tridiagonal form, working with
*           the upper triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            AB( KD1, 1 ) = DBLE( AB( KD1, 1 ) )
            DO 90 I = 1, N - 2
*
*              Reduce i-th row of matrix to tridiagonal form
*
               DO 80 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
                     CALL ZLARGV( NR, AB( 1, J1-1 ), INCA, WORK( J1 ),
     $                            KD1, D( J1 ), KD1 )
*
*                    apply rotations from the right
*
*
*                    Dependent on the the number of diagonals either
*                    ZLARTV or ZROT is used
*
                     IF( NR.GE.2*KD-1 ) THEN
                        DO 10 L = 1, KD - 1
                           CALL ZLARTV( NR, AB( L+1, J1-1 ), INCA,
     $                                  AB( L, J1 ), INCA, D( J1 ),
     $                                  WORK( J1 ), KD1 )
   10                   CONTINUE
*
                     ELSE
                        JEND = J1 + ( NR-1 )*KD1
                        DO 20 JINC = J1, JEND, KD1
                           CALL ZROT( KDM1, AB( 2, JINC-1 ), 1,
     $                                AB( 1, JINC ), 1, D( JINC ),
     $                                WORK( JINC ) )
   20                   CONTINUE
                     END IF
                  END IF
*
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i,i+k-1)
*                       within the band
*
                        CALL ZLARTG( AB( KD-K+3, I+K-2 ),
     $                               AB( KD-K+2, I+K-1 ), D( I+K-1 ),
     $                               WORK( I+K-1 ), TEMP )
                        AB( KD-K+3, I+K-2 ) = TEMP
*
*                       apply rotation from the right
*
                        CALL ZROT( K-3, AB( KD-K+4, I+K-2 ), 1,
     $                             AB( KD-K+3, I+K-1 ), 1, D( I+K-1 ),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
                  IF( NR.GT.0 )
     $               CALL ZLAR2V( NR, AB( KD1, J1-1 ), AB( KD1, J1 ),
     $                            AB( KD, J1 ), INCA, D( J1 ),
     $                            WORK( J1 ), KD1 )
*
*                 apply plane rotations from the left
*
                  IF( NR.GT.0 ) THEN
                     CALL ZLACGV( NR, WORK( J1 ), KD1 )
                     IF( 2*KD-1.LT.NR ) THEN
*
*                    Dependent on the the number of diagonals either
*                    ZLARTV or ZROT is used
*
                        DO 30 L = 1, KD - 1
                           IF( J2+L.GT.N ) THEN
                              NRT = NR - 1
                           ELSE
                              NRT = NR
                           END IF
                           IF( NRT.GT.0 )
     $                        CALL ZLARTV( NRT, AB( KD-L, J1+L ), INCA,
     $                                     AB( KD-L+1, J1+L ), INCA,
     $                                     D( J1 ), WORK( J1 ), KD1 )
   30                   CONTINUE
                     ELSE
                        J1END = J1 + KD1*( NR-2 )
                        IF( J1END.GE.J1 ) THEN
                           DO 40 JIN = J1, J1END, KD1
                              CALL ZROT( KD-1, AB( KD-1, JIN+1 ), INCX,
     $                                   AB( KD, JIN+1 ), INCX,
     $                                   D( JIN ), WORK( JIN ) )
   40                      CONTINUE
                        END IF
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 )
     $                     CALL ZROT( LEND, AB( KD-1, LAST+1 ), INCX,
     $                                AB( KD, LAST+1 ), INCX, D( LAST ),
     $                                WORK( LAST ) )
                     END IF
                  END IF
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
                     IF( INITQ ) THEN
*
*                 take advantage of the fact that Q was
*                 initially the Identity matrix
*
                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 )
     $                     IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 50 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
                           CALL ZROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
     $                                1, D( J ), DCONJG( WORK( J ) ) )
   50                   CONTINUE
                     ELSE
*
                        DO 60 J = J1, J2, KD1
                           CALL ZROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                                D( J ), DCONJG( WORK( J ) ) )
   60                   CONTINUE
                     END IF
*
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 70 J = J1, J2, KD1
*
*                    create nonzero element a(j-1,j+kd) outside the band
*                    and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( 1, J+KD )
                     AB( 1, J+KD ) = D( J )*AB( 1, J+KD )
   70             CONTINUE
   80          CONTINUE
   90       CONTINUE
         END IF
*
         IF( KD.GT.0 ) THEN
*
*           make off-diagonal elements real and copy them to E
*
            DO 100 I = 1, N - 1
               T = AB( KD, I+1 )
               ABST = ABS( T )
               AB( KD, I+1 ) = ABST
               E( I ) = ABST
               IF( ABST.NE.ZERO ) THEN
                  T = T / ABST
               ELSE
                  T = CONE
               END IF
               IF( I.LT.N-1 )
     $            AB( KD, I+2 ) = AB( KD, I+2 )*T
               IF( WANTQ ) THEN
                  CALL ZSCAL( N, DCONJG( T ), Q( 1, I+1 ), 1 )
               END IF
  100       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 110 I = 1, N - 1
               E( I ) = ZERO
  110       CONTINUE
         END IF
*
*        copy diagonal elements to D
*
         DO 120 I = 1, N
            D( I ) = AB( KD1, I )
  120    CONTINUE
*
      ELSE
*
         IF( KD.GT.1 ) THEN
*
*           Reduce to complex Hermitian tridiagonal form, working with
*           the lower triangle
*
            NR = 0
            J1 = KDN + 2
            J2 = 1
*
            AB( 1, 1 ) = DBLE( AB( 1, 1 ) )
            DO 210 I = 1, N - 2
*
*              Reduce i-th column of matrix to tridiagonal form
*
               DO 200 K = KDN + 1, 2, -1
                  J1 = J1 + KDN
                  J2 = J2 + KDN
*
                  IF( NR.GT.0 ) THEN
*
*                    generate plane rotations to annihilate nonzero
*                    elements which have been created outside the band
*
                     CALL ZLARGV( NR, AB( KD1, J1-KD1 ), INCA,
     $                            WORK( J1 ), KD1, D( J1 ), KD1 )
*
*                    apply plane rotations from one side
*
*
*                    Dependent on the the number of diagonals either
*                    ZLARTV or ZROT is used
*
                     IF( NR.GT.2*KD-1 ) THEN
                        DO 130 L = 1, KD - 1
                           CALL ZLARTV( NR, AB( KD1-L, J1-KD1+L ), INCA,
     $                                  AB( KD1-L+1, J1-KD1+L ), INCA,
     $                                  D( J1 ), WORK( J1 ), KD1 )
  130                   CONTINUE
                     ELSE
                        JEND = J1 + KD1*( NR-1 )
                        DO 140 JINC = J1, JEND, KD1
                           CALL ZROT( KDM1, AB( KD, JINC-KD ), INCX,
     $                                AB( KD1, JINC-KD ), INCX,
     $                                D( JINC ), WORK( JINC ) )
  140                   CONTINUE
                     END IF
*
                  END IF
*
                  IF( K.GT.2 ) THEN
                     IF( K.LE.N-I+1 ) THEN
*
*                       generate plane rotation to annihilate a(i+k-1,i)
*                       within the band
*
                        CALL ZLARTG( AB( K-1, I ), AB( K, I ),
     $                               D( I+K-1 ), WORK( I+K-1 ), TEMP )
                        AB( K-1, I ) = TEMP
*
*                       apply rotation from the left
*
                        CALL ZROT( K-3, AB( K-2, I+1 ), LDAB-1,
     $                             AB( K-1, I+1 ), LDAB-1, D( I+K-1 ),
     $                             WORK( I+K-1 ) )
                     END IF
                     NR = NR + 1
                     J1 = J1 - KDN - 1
                  END IF
*
*                 apply plane rotations from both sides to diagonal
*                 blocks
*
                  IF( NR.GT.0 )
     $               CALL ZLAR2V( NR, AB( 1, J1-1 ), AB( 1, J1 ),
     $                            AB( 2, J1-1 ), INCA, D( J1 ),
     $                            WORK( J1 ), KD1 )
*
*                 apply plane rotations from the right
*
*
*                    Dependent on the the number of diagonals either
*                    ZLARTV or ZROT is used
*
                  IF( NR.GT.0 ) THEN
                     CALL ZLACGV( NR, WORK( J1 ), KD1 )
                     IF( NR.GT.2*KD-1 ) THEN
                        DO 150 L = 1, KD - 1
                           IF( J2+L.GT.N ) THEN
                              NRT = NR - 1
                           ELSE
                              NRT = NR
                           END IF
                           IF( NRT.GT.0 )
     $                        CALL ZLARTV( NRT, AB( L+2, J1-1 ), INCA,
     $                                     AB( L+1, J1 ), INCA, D( J1 ),
     $                                     WORK( J1 ), KD1 )
  150                   CONTINUE
                     ELSE
                        J1END = J1 + KD1*( NR-2 )
                        IF( J1END.GE.J1 ) THEN
                           DO 160 J1INC = J1, J1END, KD1
                              CALL ZROT( KDM1, AB( 3, J1INC-1 ), 1,
     $                                   AB( 2, J1INC ), 1, D( J1INC ),
     $                                   WORK( J1INC ) )
  160                      CONTINUE
                        END IF
                        LEND = MIN( KDM1, N-J2 )
                        LAST = J1END + KD1
                        IF( LEND.GT.0 )
     $                     CALL ZROT( LEND, AB( 3, LAST-1 ), 1,
     $                                AB( 2, LAST ), 1, D( LAST ),
     $                                WORK( LAST ) )
                     END IF
                  END IF
*
*
*
                  IF( WANTQ ) THEN
*
*                    accumulate product of plane rotations in Q
*
                     IF( INITQ ) THEN
*
*                 take advantage of the fact that Q was
*                 initially the Identity matrix
*
                        IQEND = MAX( IQEND, J2 )
                        I2 = MAX( 0, K-3 )
                        IQAEND = 1 + I*KD
                        IF( K.EQ.2 )
     $                     IQAEND = IQAEND + KD
                        IQAEND = MIN( IQAEND, IQEND )
                        DO 170 J = J1, J2, KD1
                           IBL = I - I2 / KDM1
                           I2 = I2 + 1
                           IQB = MAX( 1, J-IBL )
                           NQ = 1 + IQAEND - IQB
                           IQAEND = MIN( IQAEND+KD, IQEND )
                           CALL ZROT( NQ, Q( IQB, J-1 ), 1, Q( IQB, J ),
     $                                1, D( J ), WORK( J ) )
  170                   CONTINUE
                     ELSE
*
                        DO 180 J = J1, J2, KD1
                           CALL ZROT( N, Q( 1, J-1 ), 1, Q( 1, J ), 1,
     $                                D( J ), WORK( J ) )
  180                   CONTINUE
                     END IF
                  END IF
*
                  IF( J2+KDN.GT.N ) THEN
*
*                    adjust J2 to keep within the bounds of the matrix
*
                     NR = NR - 1
                     J2 = J2 - KDN - 1
                  END IF
*
                  DO 190 J = J1, J2, KD1
*
*                    create nonzero element a(j+kd,j-1) outside the
*                    band and store it in WORK
*
                     WORK( J+KD ) = WORK( J )*AB( KD1, J )
                     AB( KD1, J ) = D( J )*AB( KD1, J )
  190             CONTINUE
  200          CONTINUE
  210       CONTINUE
         END IF
*
         IF( KD.GT.0 ) THEN
*
*           make off-diagonal elements real and copy them to E
*
            DO 220 I = 1, N - 1
               T = AB( 2, I )
               ABST = ABS( T )
               AB( 2, I ) = ABST
               E( I ) = ABST
               IF( ABST.NE.ZERO ) THEN
                  T = T / ABST
               ELSE
                  T = CONE
               END IF
               IF( I.LT.N-1 )
     $            AB( 2, I+1 ) = AB( 2, I+1 )*T
               IF( WANTQ ) THEN
                  CALL ZSCAL( N, T, Q( 1, I+1 ), 1 )
               END IF
  220       CONTINUE
         ELSE
*
*           set E to zero if original matrix was diagonal
*
            DO 230 I = 1, N - 1
               E( I ) = ZERO
  230       CONTINUE
         END IF
*
*        copy diagonal elements to D
*
         DO 240 I = 1, N
            D( I ) = AB( 1, I )
  240    CONTINUE
      END IF
*
      RETURN
*
*     End of ZHBTRD
*
      END
      DOUBLE PRECISION FUNCTION ZLANHB( NORM, UPLO, N, K, AB, LDAB,
     $                 WORK )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          NORM, UPLO
      INTEGER            K, LDAB, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   WORK( * )
      COMPLEX*16         AB( LDAB, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLANHB  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the element of  largest absolute value  of an
*  n by n hermitian band matrix A,  with k super-diagonals.
*
*  Description
*  ===========
*
*  ZLANHB returns the value
*
*     ZLANHB = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a consistent matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in ZLANHB as described
*          above.
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the upper or lower triangular part of the
*          band matrix A is supplied.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, ZLANHB is
*          set to zero.
*
*  K       (input) INTEGER
*          The number of super-diagonals or sub-diagonals of the
*          band matrix A.  K >= 0.
*
*  AB      (input) COMPLEX*16 array, dimension (LDAB,N)
*          The upper or lower triangle of the hermitian band matrix A,
*          stored in the first K+1 rows of AB.  The j-th column of A is
*          stored in the j-th column of the array AB as follows:
*          if UPLO = 'U', AB(k+1+i-j,j) = A(i,j) for max(1,j-k)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)   = A(i,j) for j<=i<=min(n,j+k).
*          Note that the imaginary parts of the diagonal elements need
*          not be set and are assumed to be zero.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= K+1.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (MAX(1,LWORK)),
*          where LWORK >= N when NORM = 'I' or '1' or 'O'; otherwise,
*          WORK is not referenced.
*
* =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ONE, ZERO
      PARAMETER          ( ONE = 1.0D+0, ZERO = 0.0D+0 )
*     ..
*     .. Local Scalars ..
      INTEGER            I, J, L
      DOUBLE PRECISION   ABSA, SCALE, SUM, VALUE
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           ZLASSQ
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, MAX, MIN, SQRT
*     ..
*     .. Executable Statements ..
*
      IF( N.EQ.0 ) THEN
         VALUE = ZERO
      ELSE IF( LSAME( NORM, 'M' ) ) THEN
*
*        Find max(abs(A(i,j))).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 20 J = 1, N
               DO 10 I = MAX( K+2-J, 1 ), K
                  VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
   10          CONTINUE
               VALUE = MAX( VALUE, ABS( DBLE( AB( K+1, J ) ) ) )
   20       CONTINUE
         ELSE
            DO 40 J = 1, N
               VALUE = MAX( VALUE, ABS( DBLE( AB( 1, J ) ) ) )
               DO 30 I = 2, MIN( N+1-J, K+1 )
                  VALUE = MAX( VALUE, ABS( AB( I, J ) ) )
   30          CONTINUE
   40       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'I' ) ) .OR. ( LSAME( NORM, 'O' ) ) .OR.
     $         ( NORM.EQ.'1' ) ) THEN
*
*        Find normI(A) ( = norm1(A), since A is hermitian).
*
         VALUE = ZERO
         IF( LSAME( UPLO, 'U' ) ) THEN
            DO 60 J = 1, N
               SUM = ZERO
               L = K + 1 - J
               DO 50 I = MAX( 1, J-K ), J - 1
                  ABSA = ABS( AB( L+I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   50          CONTINUE
               WORK( J ) = SUM + ABS( DBLE( AB( K+1, J ) ) )
   60       CONTINUE
            DO 70 I = 1, N
               VALUE = MAX( VALUE, WORK( I ) )
   70       CONTINUE
         ELSE
            DO 80 I = 1, N
               WORK( I ) = ZERO
   80       CONTINUE
            DO 100 J = 1, N
               SUM = WORK( J ) + ABS( DBLE( AB( 1, J ) ) )
               L = 1 - J
               DO 90 I = J + 1, MIN( N, J+K )
                  ABSA = ABS( AB( L+I, J ) )
                  SUM = SUM + ABSA
                  WORK( I ) = WORK( I ) + ABSA
   90          CONTINUE
               VALUE = MAX( VALUE, SUM )
  100       CONTINUE
         END IF
      ELSE IF( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) THEN
*
*        Find normF(A).
*
         SCALE = ZERO
         SUM = ONE
         IF( K.GT.0 ) THEN
            IF( LSAME( UPLO, 'U' ) ) THEN
               DO 110 J = 2, N
                  CALL ZLASSQ( MIN( J-1, K ), AB( MAX( K+2-J, 1 ), J ),
     $                         1, SCALE, SUM )
  110          CONTINUE
               L = K + 1
            ELSE
               DO 120 J = 1, N - 1
                  CALL ZLASSQ( MIN( N-J, K ), AB( 2, J ), 1, SCALE,
     $                         SUM )
  120          CONTINUE
               L = 1
            END IF
            SUM = 2*SUM
         ELSE
            L = 1
         END IF
         DO 130 J = 1, N
            IF( DBLE( AB( L, J ) ).NE.ZERO ) THEN
               ABSA = ABS( DBLE( AB( L, J ) ) )
               IF( SCALE.LT.ABSA ) THEN
                  SUM = ONE + SUM*( SCALE / ABSA )**2
                  SCALE = ABSA
               ELSE
                  SUM = SUM + ( ABSA / SCALE )**2
               END IF
            END IF
  130    CONTINUE
         VALUE = SCALE*SQRT( SUM )
      END IF
*
      ZLANHB = VALUE
      RETURN
*
*     End of ZLANHB
*
      END
      SUBROUTINE ZLASWP( N, A, LDA, K1, K2, IPIV, INCX )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCX, K1, K2, LDA, N
*     ..
*     .. Array Arguments ..
      INTEGER            IPIV( * )
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*
*  Further Details
*  ===============
*
*  Modified by
*   R. C. Whaley, Computer Science Dept., Univ. of Tenn., Knoxville, USA
*
* =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, I1, I2, INC, IP, IX, IX0, J, K, N32
      COMPLEX*16         TEMP
*     ..
*     .. Executable Statements ..
*
*     Interchange row I with row IPIV(I) for each of rows K1 through K2.
*
      IF( INCX.GT.0 ) THEN
         IX0 = K1
         I1 = K1
         I2 = K2
         INC = 1
      ELSE IF( INCX.LT.0 ) THEN
         IX0 = 1 + ( 1-K2 )*INCX
         I1 = K2
         I2 = K1
         INC = -1
      ELSE
         RETURN
      END IF
*
      N32 = ( N / 32 )*32
      IF( N32.NE.0 ) THEN
         DO 30 J = 1, N32, 32
            IX = IX0
            DO 20 I = I1, I2, INC
               IP = IPIV( IX )
               IF( IP.NE.I ) THEN
                  DO 10 K = J, J + 31
                     TEMP = A( I, K )
                     A( I, K ) = A( IP, K )
                     A( IP, K ) = TEMP
   10             CONTINUE
               END IF
               IX = IX + INCX
   20       CONTINUE
   30    CONTINUE
      END IF
      IF( N32.NE.N ) THEN
         N32 = N32 + 1
         IX = IX0
         DO 50 I = I1, I2, INC
            IP = IPIV( IX )
            IF( IP.NE.I ) THEN
               DO 40 K = N32, N
                  TEMP = A( I, K )
                  A( I, K ) = A( IP, K )
                  A( IP, K ) = TEMP
   40          CONTINUE
            END IF
            IX = IX + INCX
   50    CONTINUE
      END IF
*
      RETURN
*
*     End of ZLASWP
*
      END
      SUBROUTINE ZTRTRI( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTRTRI computes the inverse of a complex upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ),
     $                   ZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J, JB, NB, NN
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      EXTERNAL           LSAME, ILAENV
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZTRMM, ZTRSM, ZTRTI2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX, MIN
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRTRI', -INFO )
         RETURN
      END IF
*
*     Quick return if possible
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Check for singularity if non-unit.
*
      IF( NOUNIT ) THEN
         DO 10 INFO = 1, N
            IF( A( INFO, INFO ).EQ.ZERO )
     $         RETURN
   10    CONTINUE
         INFO = 0
      END IF
*
*     Determine the block size for this environment.
*
      NB = ILAENV( 1, 'ZTRTRI', UPLO // DIAG, N, -1, -1, -1 )
      IF( NB.LE.1 .OR. NB.GE.N ) THEN
*
*        Use unblocked code
*
         CALL ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
      ELSE
*
*        Use blocked code
*
         IF( UPPER ) THEN
*
*           Compute inverse of upper triangular matrix
*
            DO 20 J = 1, N, NB
               JB = MIN( NB, N-J+1 )
*
*              Compute rows 1:j-1 of current block column
*
               CALL ZTRMM( 'Left', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, ONE, A, LDA, A( 1, J ), LDA )
               CALL ZTRSM( 'Right', 'Upper', 'No transpose', DIAG, J-1,
     $                     JB, -ONE, A( J, J ), LDA, A( 1, J ), LDA )
*
*              Compute inverse of current diagonal block
*
               CALL ZTRTI2( 'Upper', DIAG, JB, A( J, J ), LDA, INFO )
   20       CONTINUE
         ELSE
*
*           Compute inverse of lower triangular matrix
*
            NN = ( ( N-1 ) / NB )*NB + 1
            DO 30 J = NN, 1, -NB
               JB = MIN( NB, N-J+1 )
               IF( J+JB.LE.N ) THEN
*
*                 Compute rows j+jb:n of current block column
*
                  CALL ZTRMM( 'Left', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, ONE, A( J+JB, J+JB ), LDA,
     $                        A( J+JB, J ), LDA )
                  CALL ZTRSM( 'Right', 'Lower', 'No transpose', DIAG,
     $                        N-J-JB+1, JB, -ONE, A( J, J ), LDA,
     $                        A( J+JB, J ), LDA )
               END IF
*
*              Compute inverse of current diagonal block
*
               CALL ZTRTI2( 'Lower', DIAG, JB, A( J, J ), LDA, INFO )
   30       CONTINUE
         END IF
      END IF
*
      RETURN
*
*     End of ZTRTRI
*
      END
      SUBROUTINE ZLAR2V( N, X, Y, Z, INCX, C, S, INCC )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCC, INCX, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( * )
      COMPLEX*16         S( * ), X( * ), Y( * ), Z( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLAR2V applies a vector of complex plane rotations with real cosines
*  from both sides to a sequence of 2-by-2 complex Hermitian matrices,
*  defined by the elements of the vectors x, y and z. For i = 1,2,...,n
*
*     (       x(i)  z(i) ) :=
*     ( conjg(z(i)) y(i) )
*
*       (  c(i) conjg(s(i)) ) (       x(i)  z(i) ) ( c(i) -conjg(s(i)) )
*       ( -s(i)       c(i)  ) ( conjg(z(i)) y(i) ) ( s(i)        c(i)  )
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of plane rotations to be applied.
*
*  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
*          The vector x; the elements of x are assumed to be real.
*
*  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
*          The vector y; the elements of y are assumed to be real.
*
*  Z       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
*          The vector z.
*
*  INCX    (input) INTEGER
*          The increment between elements of X, Y and Z. INCX > 0.
*
*  C       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
*          The cosines of the plane rotations.
*
*  S       (input) COMPLEX*16 array, dimension (1+(N-1)*INCC)
*          The sines of the plane rotations.
*
*  INCC    (input) INTEGER
*          The increment between elements of C and S. INCC > 0.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IX
      DOUBLE PRECISION   CI, SII, SIR, T1I, T1R, T5, T6, XI, YI, ZII,
     $                   ZIR
      COMPLEX*16         SI, T2, T3, T4, ZI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE, DCMPLX, DCONJG, DIMAG
*     ..
*     .. Executable Statements ..
*
      IX = 1
      IC = 1
      DO 10 I = 1, N
         XI = DBLE( X( IX ) )
         YI = DBLE( Y( IX ) )
         ZI = Z( IX )
         ZIR = DBLE( ZI )
         ZII = DIMAG( ZI )
         CI = C( IC )
         SI = S( IC )
         SIR = DBLE( SI )
         SII = DIMAG( SI )
         T1R = SIR*ZIR - SII*ZII
         T1I = SIR*ZII + SII*ZIR
         T2 = CI*ZI
         T3 = T2 - DCONJG( SI )*XI
         T4 = DCONJG( T2 ) + SI*YI
         T5 = CI*XI + T1R
         T6 = CI*YI - T1R
         X( IX ) = CI*T5 + ( SIR*DBLE( T4 )+SII*DIMAG( T4 ) )
         Y( IX ) = CI*T6 - ( SIR*DBLE( T3 )-SII*DIMAG( T3 ) )
         Z( IX ) = CI*T3 + DCONJG( SI )*DCMPLX( T6, T1I )
         IX = IX + INCX
         IC = IC + INCC
   10 CONTINUE
      RETURN
*
*     End of ZLAR2V
*
      END
      SUBROUTINE ZLARGV( N, X, INCX, Y, INCY, C, INCC )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCC, INCX, INCY, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( * )
      COMPLEX*16         X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARGV generates a vector of complex plane rotations with real
*  cosines, determined by elements of the complex vectors x and y.
*  For i = 1,2,...,n
*
*     (        c(i)   s(i) ) ( x(i) ) = ( r(i) )
*     ( -conjg(s(i))  c(i) ) ( y(i) ) = (   0  )
*
*     where c(i)**2 + ABS(s(i))**2 = 1
*
*  The following conventions are used (these are the same as in ZLARTG,
*  but differ from the BLAS1 routine ZROTG):
*     If y(i)=0, then c(i)=1 and s(i)=0.
*     If x(i)=0, then c(i)=0 and s(i) is chosen so that r(i) is real.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of plane rotations to be generated.
*
*  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
*          On entry, the vector x.
*          On exit, x(i) is overwritten by r(i), for i = 1,...,n.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCY)
*          On entry, the vector y.
*          On exit, the sines of the plane rotations.
*
*  INCY    (input) INTEGER
*          The increment between elements of Y. INCY > 0.
*
*  C       (output) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
*          The cosines of the plane rotations.
*
*  INCC    (input) INTEGER
*          The increment between elements of C. INCC > 0.
*
*  Further Details
*  ======= =======
*
*  6-6-96 - Modified with a new algorithm by W. Kahan and J. Demmel
*
*  This version has a few statements commented out for thread safety
*  (machine parameters are computed on each entry). 10 feb 03, SJH.
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   TWO, ONE, ZERO
      PARAMETER          ( TWO = 2.0D+0, ONE = 1.0D+0, ZERO = 0.0D+0 )
      COMPLEX*16         CZERO
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
*     LOGICAL            FIRST

      INTEGER            COUNT, I, IC, IX, IY, J
      DOUBLE PRECISION   CS, D, DI, DR, EPS, F2, F2S, G2, G2S, SAFMIN,
     $                   SAFMN2, SAFMX2, SCALE
      COMPLEX*16         F, FF, FS, G, GS, R, SN
*     ..
*     .. External Functions ..
      DOUBLE PRECISION   DLAMCH, DLAPY2
      EXTERNAL           DLAMCH, DLAPY2
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DCONJG, DIMAG, INT, LOG,
     $                   MAX, SQRT
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   ABS1, ABSSQ
*     ..
*     .. Save statement ..
*     SAVE               FIRST, SAFMX2, SAFMIN, SAFMN2
*     ..
*     .. Data statements ..
*     DATA               FIRST / .TRUE. /
*     ..
*     .. Statement Function definitions ..
      ABS1( FF ) = MAX( ABS( DBLE( FF ) ), ABS( DIMAG( FF ) ) )
      ABSSQ( FF ) = DBLE( FF )**2 + DIMAG( FF )**2
*     ..
*     .. Executable Statements ..
*
*     IF( FIRST ) THEN
*        FIRST = .FALSE.
         SAFMIN = DLAMCH( 'S' )
         EPS = DLAMCH( 'E' )
         SAFMN2 = DLAMCH( 'B' )**INT( LOG( SAFMIN / EPS ) /
     $            LOG( DLAMCH( 'B' ) ) / TWO )
         SAFMX2 = ONE / SAFMN2
*     END IF
      IX = 1
      IY = 1
      IC = 1
      DO 60 I = 1, N
         F = X( IX )
         G = Y( IY )
*
*        Use identical algorithm as in ZLARTG
*
         SCALE = MAX( ABS1( F ), ABS1( G ) )
         FS = F
         GS = G
         COUNT = 0
         IF( SCALE.GE.SAFMX2 ) THEN
   10       CONTINUE
            COUNT = COUNT + 1
            FS = FS*SAFMN2
            GS = GS*SAFMN2
            SCALE = SCALE*SAFMN2
            IF( SCALE.GE.SAFMX2 )
     $         GO TO 10
         ELSE IF( SCALE.LE.SAFMN2 ) THEN
            IF( G.EQ.CZERO ) THEN
               CS = ONE
               SN = CZERO
               R = F
               GO TO 50
            END IF
   20       CONTINUE
            COUNT = COUNT - 1
            FS = FS*SAFMX2
            GS = GS*SAFMX2
            SCALE = SCALE*SAFMX2
            IF( SCALE.LE.SAFMN2 )
     $         GO TO 20
         END IF
         F2 = ABSSQ( FS )
         G2 = ABSSQ( GS )
         IF( F2.LE.MAX( G2, ONE )*SAFMIN ) THEN
*
*           This is a rare case: F is very small.
*
            IF( F.EQ.CZERO ) THEN
               CS = ZERO
               R = DLAPY2( DBLE( G ), DIMAG( G ) )
*              Do complex/real division explicitly with two real
*              divisions
               D = DLAPY2( DBLE( GS ), DIMAG( GS ) )
               SN = DCMPLX( DBLE( GS ) / D, -DIMAG( GS ) / D )
               GO TO 50
            END IF
            F2S = DLAPY2( DBLE( FS ), DIMAG( FS ) )
*           G2 and G2S are accurate
*           G2 is at least SAFMIN, and G2S is at least SAFMN2
            G2S = SQRT( G2 )
*           Error in CS from underflow in F2S is at most
*           UNFL / SAFMN2 .lt. sqrt(UNFL*EPS) .lt. EPS
*           If MAX(G2,ONE)=G2, then F2 .lt. G2*SAFMIN,
*           and so CS .lt. sqrt(SAFMIN)
*           If MAX(G2,ONE)=ONE, then F2 .lt. SAFMIN
*           and so CS .lt. sqrt(SAFMIN)/SAFMN2 = sqrt(EPS)
*           Therefore, CS = F2S/G2S / sqrt( 1 + (F2S/G2S)**2 ) = F2S/G2S
            CS = F2S / G2S
*           Make sure abs(FF) = 1
*           Do complex/real division explicitly with 2 real divisions
            IF( ABS1( F ).GT.ONE ) THEN
               D = DLAPY2( DBLE( F ), DIMAG( F ) )
               FF = DCMPLX( DBLE( F ) / D, DIMAG( F ) / D )
            ELSE
               DR = SAFMX2*DBLE( F )
               DI = SAFMX2*DIMAG( F )
               D = DLAPY2( DR, DI )
               FF = DCMPLX( DR / D, DI / D )
            END IF
            SN = FF*DCMPLX( DBLE( GS ) / G2S, -DIMAG( GS ) / G2S )
            R = CS*F + SN*G
         ELSE
*
*           This is the most common case.
*           Neither F2 nor F2/G2 are less than SAFMIN
*           F2S cannot overflow, and it is accurate
*
            F2S = SQRT( ONE+G2 / F2 )
*           Do the F2S(real)*FS(complex) multiply with two real
*           multiplies
            R = DCMPLX( F2S*DBLE( FS ), F2S*DIMAG( FS ) )
            CS = ONE / F2S
            D = F2 + G2
*           Do complex/real division explicitly with two real divisions
            SN = DCMPLX( DBLE( R ) / D, DIMAG( R ) / D )
            SN = SN*DCONJG( GS )
            IF( COUNT.NE.0 ) THEN
               IF( COUNT.GT.0 ) THEN
                  DO 30 J = 1, COUNT
                     R = R*SAFMX2
   30             CONTINUE
               ELSE
                  DO 40 J = 1, -COUNT
                     R = R*SAFMN2
   40             CONTINUE
               END IF
            END IF
         END IF
   50    CONTINUE
         C( IC ) = CS
         Y( IY ) = SN
         X( IX ) = R
         IC = IC + INCC
         IY = IY + INCY
         IX = IX + INCX
   60 CONTINUE
      RETURN
*
*     End of ZLARGV
*
      END
      SUBROUTINE ZLARTV( N, X, INCX, Y, INCY, C, S, INCC )
*
*  -- LAPACK auxiliary routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      INTEGER            INCC, INCX, INCY, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   C( * )
      COMPLEX*16         S( * ), X( * ), Y( * )
*     ..
*
*  Purpose
*  =======
*
*  ZLARTV applies a vector of complex plane rotations with real cosines
*  to elements of the complex vectors x and y. For i = 1,2,...,n
*
*     ( x(i) ) := (        c(i)   s(i) ) ( x(i) )
*     ( y(i) )    ( -conjg(s(i))  c(i) ) ( y(i) )
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of plane rotations to be applied.
*
*  X       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCX)
*          The vector x.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  Y       (input/output) COMPLEX*16 array, dimension (1+(N-1)*INCY)
*          The vector y.
*
*  INCY    (input) INTEGER
*          The increment between elements of Y. INCY > 0.
*
*  C       (input) DOUBLE PRECISION array, dimension (1+(N-1)*INCC)
*          The cosines of the plane rotations.
*
*  S       (input) COMPLEX*16 array, dimension (1+(N-1)*INCC)
*          The sines of the plane rotations.
*
*  INCC    (input) INTEGER
*          The increment between elements of C and S. INCC > 0.
*
*  =====================================================================
*
*     .. Local Scalars ..
      INTEGER            I, IC, IX, IY
      COMPLEX*16         XI, YI
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCONJG
*     ..
*     .. Executable Statements ..
*
      IX = 1
      IY = 1
      IC = 1
      DO 10 I = 1, N
         XI = X( IX )
         YI = Y( IY )
         X( IX ) = C( IC )*XI + S( IC )*YI
         Y( IY ) = C( IC )*YI - DCONJG( S( IC ) )*XI
         IX = IX + INCX
         IY = IY + INCY
         IC = IC + INCC
   10 CONTINUE
      RETURN
*
*     End of ZLARTV
*
      END
      SUBROUTINE ZTRTI2( UPLO, DIAG, N, A, LDA, INFO )
*
*  -- LAPACK routine (version 3.1) --
*     Univ. of Tennessee, Univ. of California Berkeley and NAG Ltd..
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          DIAG, UPLO
      INTEGER            INFO, LDA, N
*     ..
*     .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
*     ..
*
*  Purpose
*  =======
*
*  ZTRTI2 computes the inverse of a complex upper or lower triangular
*  matrix.
*
*  This is the Level 2 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies whether the matrix A is upper or lower triangular.
*          = 'U':  Upper triangular
*          = 'L':  Lower triangular
*
*  DIAG    (input) CHARACTER*1
*          Specifies whether or not the matrix A is unit triangular.
*          = 'N':  Non-unit triangular
*          = 'U':  Unit triangular
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading n by n upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading n by n lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*
*  =====================================================================
*
*     .. Parameters ..
      COMPLEX*16         ONE
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ) )
*     ..
*     .. Local Scalars ..
      LOGICAL            NOUNIT, UPPER
      INTEGER            J
      COMPLEX*16         AJJ
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      EXTERNAL           LSAME
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, ZSCAL, ZTRMV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          MAX
*     ..
*     .. Executable Statements ..
*
*     Test the input parameters.
*
      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      NOUNIT = LSAME( DIAG, 'N' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( .NOT.NOUNIT .AND. .NOT.LSAME( DIAG, 'U' ) ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -5
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZTRTI2', -INFO )
         RETURN
      END IF
*
      IF( UPPER ) THEN
*
*        Compute inverse of upper triangular matrix.
*
         DO 10 J = 1, N
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
*
*           Compute elements 1:j-1 of j-th column.
*
            CALL ZTRMV( 'Upper', 'No transpose', DIAG, J-1, A, LDA,
     $                  A( 1, J ), 1 )
            CALL ZSCAL( J-1, AJJ, A( 1, J ), 1 )
   10    CONTINUE
      ELSE
*
*        Compute inverse of lower triangular matrix.
*
         DO 20 J = N, 1, -1
            IF( NOUNIT ) THEN
               A( J, J ) = ONE / A( J, J )
               AJJ = -A( J, J )
            ELSE
               AJJ = -ONE
            END IF
            IF( J.LT.N ) THEN
*
*              Compute elements j+1:n of j-th column.
*
               CALL ZTRMV( 'Lower', 'No transpose', DIAG, N-J,
     $                     A( J+1, J+1 ), LDA, A( J+1, J ), 1 )
               CALL ZSCAL( N-J, AJJ, A( J+1, J ), 1 )
            END IF
   20    CONTINUE
      END IF
*
      RETURN
*
*     End of ZTRTI2
*
      END

