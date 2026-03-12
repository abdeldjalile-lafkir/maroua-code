      SUBROUTINE CG(N, B, X,   ITER, RESID, MATVEC,
     $                    PSOLVE, INFO )

c     This subroutine solves the linear system AX = B
c     with A a square N*N symmetric matrix, while
c     B is the right-hand-side. The subroutine uses
c     the Conjugate Gradient method.
c     Here N denotes the size of the system
c     B is a vector of dimension N (the RHS)
c     X is a vector of dimension N (the unkown)
c     ITER : an integer. 
c           As input: the maximal number of iterations (for example 1000)
c           As output:the effective number of iterations
c     RESID : a real
c             As input : the stopping residue (for example 1.0D-6)
c             Convergence holds if the relative difference between
c             two consecutive terms is smaller then RESID
c             As output : the effective residue
c     MATVEC : a subroutine (external)   
c              MATVEC has four parameters
c              (alpha, X, beta, Y)
c              where alpha and beta are reals
c              and X and Y are vectors 
c              MATVEC changes Y and Y = alpha*A*X + beta*Y
c     (this is the most important subroutine because it
c      represents the matrix A)
c     PSOLVE : is a preconditonning subroutine which
c              has two parameters X and Y
c              X and Y are vectors
c              At output Y = M*X
c               with M an approximation of the inverse of A
c               If no precontionning, M = I
c     INFO    : an integer which gives the information about convergence

      INTEGER           N, LDW, ITER, INFO
      double precision  RESID
      double precision  X( * ), B( * ), WORK(N+100, 7)
      EXTERNAL          MATVEC, PSOLVE
      double precision              ZERO, ONE
      PARAMETER       ( ZERO = 0.0D0, ONE = 1.0D0)
      INTEGER           MAXIT, R, Z, P, Q
      double precision  TOL, ALPHA, BETA,RHO,RHO1,BNRM2,SDOT,SNRM2 
      EXTERNAL          SAXPY, SCOPY, SDOT, SNRM2
  
      LDW = N + 100
      INFO = 0
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT. MAX(1, N ) ) THEN 
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE. 0 ) RETURN
      MAXIT = ITER
      TOL = RESID
      R = 1
      Z = 2
      P = 3
      Q = 4
      CALL SCOPY( N, B, 1, WORK(1,R), 1)
      IF ( SNRM2( N, X, 1 ).NE.ZERO ) THEN
        CALL MATVEC( -ONE, X, ONE, WORK(1,R) )	
        IF ( SNRM2(N, WORK(1,R), 1 ).LE.TOL ) GO TO 30
      ENDIF
      BNRM2 = SNRM2 ( N, B, 1)
      IF ( BNRM2 .EQ. ZERO ) BNRM2 = ONE
      ITER = 0
10    CONTINUE
      ITER = ITER + 1
      CALL PSOLVE( WORK(1,Z) , WORK(1,R) )
      RHO = SDOT(N, WORK(1,R), 1, WORK(1,Z), 1 )
      IF ( ITER.GT.1 ) THEN
         BETA =  RHO / RHO1 
         CALL SAXPY( N, BETA, WORK(1,P), 1,WORK(1,Z),1 )
         CALL SCOPY( N, WORK(1,Z),1,WORK(1,P), 1 ) 
      ELSE
         CALL SCOPY (N, WORK(1,Z),1,WORK(1,P), 1 ) 
      ENDIF 
      CALL MATVEC( ONE, WORK(1,P), ZERO,WORK(1,Q) )
      ALPHA = RHO / SDOT( N, WORK(1,P), 1, WORK(1,Q), 1 )
      CALL SAXPY( N, ALPHA, WORK(1,P), 1, X, 1 ) 
      CALL SAXPY( N, -ALPHA, WORK(1,Q), 1, WORK(1, R), 1 )
      RESID = SNRM2( N, WORK(1,R), 1 ) / BNRM2
       if ((iter - (iter/400)*400).eq.0) 
     &             WRITE(*, *) 'ITERATION :', ITER ,  ', RESIDU:', RESID
      IF ( RESID.LE. TOL ) GO TO 30
      IF ( ITER.EQ.MAXIT ) GO TO 20
      RHO1 = RHO 
      GO TO 10
    
20    CONTINUE
      INFO = ITER
      RETURN

30    CONTINUE
      RETURN
      END




 
c---------------------------------------------------------------------

      SUBROUTINE BICGSTAB(N, B, X,  ITER, RESID, MATVEC,
     $                    PSOLVE, INFO )
      INTEGER           N, LDW, ITER, INFO
      double precision              RESID
      double precision              X( * ), B( * ), WORK(N+100, 7)
      EXTERNAL          MATVEC, PSOLVE
      double precision  ZERO, ONE
      PARAMETER       ( ZERO = 0.0D0,ONE = 1.0D0)
      INTEGER           R, RTLD, P, PHAT, V, S, SHAT, T, MAXIT
      double precision  TOL, ALPHA, BETA, RHO, RHO1, BNRM2, OMEGA,
     $                  RHOTOL, OMEGATOL, SDOT, SNRM2 
      EXTERNAL          SAXPY, SCOPY, SDOT, SNRM2, SSCAL
      INTRINSIC         ABS, MAX


      LDW = N + 100
      INFO = 0
      IF ( N.LT.0 ) THEN
         INFO = -1
      ELSE IF ( LDW.LT. MAX(1, N ) ) THEN
         INFO = -2
      ELSE IF ( ITER.LE.0 ) THEN
         INFO = -3
      ENDIF
      IF ( INFO.NE. 0 ) RETURN
      MAXIT = ITER
      TOL = RESID
      R = 1
      RTLD = 2 
      P = 3
      V =  4
      T = 5
      PHAT = 6
      SHAT = 7
      S = 1
      RHOTOL = 1.D-13
      OMEGATOL = 1.D-13
      CALL SCOPY( N, B, 1, WORK(1,R), 1)
      IF ( SNRM2( N, X, 1 ).NE.ZERO ) THEN 
          CALL MATVEC(-ONE, X, ONE, WORK(1,R))
         IF ( SNRM2(N, WORK(1,R), 1 ).LE.TOL ) GO TO 30
      ENDIF
      CALL SCOPY (N, WORK(1,R), 1, WORK(1,RTLD), 1 )
      BNRM2 = SNRM2 ( N, B, 1)
      IF ( BNRM2 .EQ. ZERO ) BNRM2 = ONE
      ITER = 0
10    CONTINUE
      ITER = ITER + 1
      RHO = SDOT(N, WORK(1,RTLD), 1, WORK(1,R), 1 )
      IF ( ABS( RHO ).LT.RHOTOL ) GO TO 25
      IF ( ITER.GT.1 ) THEN
         BETA = ( RHO / RHO1 ) * ( ALPHA / OMEGA )
         CALL SAXPY( N, -OMEGA, WORK(1,V), 1,WORK(1,P),1 )
         CALL SSCAL ( N, BETA, WORK(1,P),1)
         CALL SAXPY(N, ONE, WORK(1,R), 1, WORK(1,P), 1 )
         ELSE
            CALL SCOPY (N, WORK(1,R),1,WORK(1,P), 1 ) 
          ENDIF 
      CALL PSOLVE( WORK(1,PHAT) , WORK(1,P) )
      CALL MATVEC( ONE, WORK(1,PHAT), ZERO,WORK(1,V) )
      ALPHA = RHO / SDOT( N, WORK(1,RTLD), 1, WORK(1,V), 1 )
      CALL SAXPY( N, -ALPHA, WORK(1,V), 1, WORK(1,R), 1 ) 
      CALL SCOPY( N, WORK(1,R), 1, WORK(1,S), 1 )
      IF ( SNRM2 ( N, WORK(1,S), 1 ) .LE.TOL ) THEN
         CALL SAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
         RESID = SNRM2( N, WORK(1,S), 1 ) / BNRM2
         GO TO 30
      ELSE
         CALL PSOLVE( WORK(1,SHAT), WORK(1,S) )
         CALL MATVEC( ONE, WORK(1,SHAT), ZERO, WORK(1,T) )
         OMEGA = SDOT( N, WORK(1,T), 1, WORK(1,S), 1 ) /
     $           SDOT( N, WORK(1,T), 1, WORK(1,T), 1 ) 
         CALL SAXPY( N, ALPHA, WORK(1,PHAT), 1, X, 1 )
         CALL SAXPY( N, OMEGA, WORK(1,SHAT), 1, X, 1 )
         CALL SAXPY( N, -OMEGA, WORK(1,T), 1, WORK(1,R), 1 )
         RESID = SNRM2( N, WORK(1,R), 1 ) / BNRM2
c         WRITE(*, *) 'ITERATION :',ITER, ', RESIDU:', RESID
         IF ( RESID.LE. TOL ) GO TO 30
         IF ( ITER.EQ.MAXIT ) GO TO 20
      ENDIF
      IF (ABS( OMEGA ) .LT.OMEGATOL ) THEN
      GO TO 25
      ELSE
         RHO1 = RHO 
         GO TO 10
      ENDIF
20    CONTINUE
      INFO = ITER
      RETURN
25    CONTINUE
      IF ( ABS( RHO ) .LT.RHOTOL ) THEN
         INFO = -10
      ELSE IF ( ABS( OMEGA ) .LT. OMEGATOL ) THEN
         INFO = -11
      ENDIF
      RETURN
	
30    CONTINUE
      RETURN
      END

c-------------------------------------------------------------
      subroutine scopy(n,sx,incx,sy,incy)
      implicit none
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
      double precision sx(1),sy(1)
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c     code for both increments equal to 1
c
c
c     clean-up loop
c
   20 m = mod(n,7)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sx(i)
   30 continue
      if( n .lt. 7 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,7
       sy(i) = sx(i)
       sy(i + 1) = sx(i + 1)
       sy(i + 2) = sx(i + 2)
       sy(i + 3) = sx(i + 3)
       sy(i + 4) = sx(i + 4)
       sy(i + 5) = sx(i + 5)
       sy(i + 6) = sx(i + 6)
   50 continue
      return
      end
c-------------------------------------------------------------
      double precision function sdot(n,sx,incx,sy,incy)
      implicit none
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      double precision sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      stemp = 0.0d0
      sdot = 0.0d0
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = stemp + sx(ix)*sy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      sdot = stemp
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)*sy(i)
   30 continue
      if( n .lt. 5 ) go to 60
   40 mp1 = m + 1
      do 55 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     &  sx(i+2)*sy(i+2) + sx(i+3)*sy(i+3) + sx(i+4)*sy(i+4)
   55 continue
   60 sdot = stemp
      return
      end
c-------------------------------------------------------------
      double precision function snrm2 ( n, sx, incx)
	implicit none
      integer i, incx, ix, j, n, next
      double precision   sx(1), cutlo,cuthi, hitest,sum,xmax,zero, one
      data   zero, one /0.0D0, 1.0D0/
c
c     euclidean norm of the n-vector stored in sx() with storage
c     increment incx .
c     if    n .le. 0 return with result = 0.
c     if n .ge. 1 then incx must be .ge. 1
c
c           c.l.lawson, 1978 jan 08
c     modified to correct problem with negative increment, 8/21/90.
c
c     four phase method     using two built-in constants that are
c     hopefully applicable to all machines.
c         cutlo = maximum of  sqrt(u/eps)  over all known machines.
c         cuthi = minimum of  sqrt(v)      over all known machines.
c     where
c         eps = smallest no. such that eps + 1. .gt. 1.
c         u   = smallest positive no.   (underflow limit)
c         v   = largest  no.            (overflow  limit)
c
c     brief outline of algorithm..
c
c     phase 1    scans zero components.
c     move to phase 2 when a component is nonzero and .le. cutlo
c     move to phase 3 when a component is .gt. cutlo
c     move to phase 4 when a component is .ge. cuthi/m
c     where m = n for x() double precision and m = 2*n for complex.
c
c     values for cutlo and cuthi..
c     from the environmental parameters listed in the imsl converter
c     document the limiting values are as follows..
c     cutlo, s.p.   u/eps = 2**(-102) for  honeywell.  close seconds are
c                   univac and dec at 2**(-103)
c                   thus cutlo = 2**(-51) = 4.44089e-16
c     cuthi, s.p.   v = 2**127 for univac, honeywell, and dec.
c                   thus cuthi = 2**(63.5) = 1.30438e19
c     cutlo, d.p.   u/eps = 2**(-67) for honeywell and dec.
c                   thus cutlo = 2**(-33.5) = 8.23181d-11
c     cuthi, d.p.   same as s.p.  cuthi = 1.30438d19
c     data cutlo, cuthi / 8.232d-11,  1.304d19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
c     data cutlo, cuthi / 4.441e-16,  1.304e19 /
c...  from Ed Anderson (for Cray)
c     data cutlo / 0300315520236314774737b /
c     data cuthi / 0500004000000000000000b /
c...  from Ed Anderson (for Sun4)
      data cutlo /   0.44408921D-15 /
      data cuthi /   0.18446743D+20 /
c
      if(n .gt. 0) go to 10
         snrm2  = zero
         go to 300
c
   10 assign 30 to next
      sum = zero
      i = 1
      if(incx.lt.0) i = (-n+1)*incx + 1
      ix = 1
c                                       begin main loop
   20    go to next, (30, 50, 70, 110)
   30 if( abs(sx(i)) .gt. cutlo) go to 85
      assign 50 to next
      xmax = zero
c
c                        phase 1.  sum is zero
c
   50 if( sx(i) .eq. zero) go to 200
      if( abs(sx(i)) .gt. cutlo) go to 85
c
c                                prepare for phase 2.
      assign 70 to next
      go to 105
c
c                                prepare for phase 4.
c
  100 continue
      assign 110 to next
      sum = (sum / sx(i)) / sx(i)
  105 xmax = abs(sx(i))
      go to 115
c
c                   phase 2.  sum is small.
c                             scale to avoid destructive underflow.
c
   70 if( abs(sx(i)) .gt. cutlo ) go to 75
c
c                     common code for phases 2 and 4.
c                     in phase 4 sum is large.  scale to avoid overflow.
c
  110 if( abs(sx(i)) .le. xmax ) go to 115
         sum = one + sum * (xmax / sx(i))**2
         xmax = abs(sx(i))
         go to 200
c
  115 sum = sum + (sx(i)/xmax)**2
      go to 200
c
c
c                  prepare for phase 3.
c
   75 sum = (sum * xmax) * xmax
c
c
c     for real or d.p. set hitest = cuthi/n
c     for complex      set hitest = cuthi/(2*n)
c
   85 hitest = cuthi/float( n )
c
c                   phase 3.  sum is mid-range.  no scaling.
c
      do 95 j = ix, n
         if(abs(sx(i)) .ge. hitest) go to 100
         sum = sum + sx(i)**2
         i = i + incx
   95 continue
      snrm2 = sqrt( sum )
      go to 300
c
  200 continue
      ix = ix + 1
      i = i + incx
      if( ix .le. n ) go to 20
c
c              end of main loop.
c
c              compute square root and adjust for scaling.
c
      snrm2 = xmax * sqrt(sum)
  300 continue
      return
      end

c-------------------------------------------------------------
      subroutine sscal(n,sa,sx,incx)
      implicit none
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c     modified to correct problem with negative increment, 8/21/90.
c
      double precision sa,sx(1)
      integer i,incx,ix,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1 
      if(incx.lt.0)ix = (-n+1)*incx + 1 
      do 10 i = 1,n 
        sx(ix) = sa*sx(ix)
        ix = ix + incx 
   10 continue
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,5)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sx(i) = sa*sx(i)
   30 continue
      if( n .lt. 5 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,5
        sx(i) = sa*sx(i)
        sx(i + 1) = sa*sx(i + 1)
        sx(i + 2) = sa*sx(i + 2)
        sx(i + 3) = sa*sx(i + 3)
        sx(i + 4) = sa*sx(i + 4)
   50 continue
      return
      end

c-------------------------------------------------------------
      subroutine saxpy(n,sa,sx,incx,sy,incy)
c
c     constant times a vector plus a vector.
c     uses unrolled loop for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
      implicit none
      double precision sx(1),sy(1),sa
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (sa .eq. 0.0) return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c        code for unequal increments or equal increments
c          not equal to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        sy(iy) = sy(iy) + sa*sx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,4)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        sy(i) = sy(i) + sa*sx(i)
   30 continue
      if( n .lt. 4 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,4
        sy(i) = sy(i) + sa*sx(i)
        sy(i + 1) = sy(i + 1) + sa*sx(i + 1)
        sy(i + 2) = sy(i + 2) + sa*sx(i + 2)
        sy(i + 3) = sy(i + 3) + sa*sx(i + 3)
   50 continue
      return
      end

c-------------------------------------------------------------





