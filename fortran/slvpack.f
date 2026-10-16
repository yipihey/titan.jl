      subroutine bglsdc(n, l, m, a, ia, iflag)
c
c     description
c
c     given a banded matrix a with l lower subdiagonals and m super
c     diagonals, this routine factors the matrix into a = lu.  these
c     factors can then be used by routine bglssl to solve the system
c     a*x = lu*x = y.  the upper and lower bandwidths must both be
c     less than or equal to 64.  no pivoting is performed.
c
c     the m super diagonals of the matrix a are stored in the first
c     m rows of the array a.  the main diagonal of the matrix a is
c     stored in the m+1 row of the array a.  the l subdiagonals of the
c     matrix a are stored in rows m+2 through m+l+1 of the array a.
c     the columns of the matrix a correspond to the columns of array a.
c     for example, if n=4, l=2, and m=1, the matrix elements a(i,j)
c     will be stored in array a as follows:
c
c          0       a(1,2)  a(2,3)  a(3,4)
c          a(1,1)  a(2,2)  a(3,3)  a(4,4)
c          a(2,1)  a(3,2)  a(4,3)  0
c          a(3,1)  a(4,2)  0       0
c
c     on entry
c        n = the order of the system
c        l = the number of lower diagonals  l .le. 64
c        m = the number of upper diagonals  m .le. 64
c        a = a 2d array that contains the matrix to be factored
c            the diagonals are stored in the rows of a as
c            described above
c       ia = the first dimension of a      ia .gt. m + l
c
c     on return
c         a = lu, where l is in place of the lower diagonals and u
c             is in rows above l
c     iflag = 0 if ok
c           = 1 if (n .le. 0 .or. m + l .ge. ia)
c           = 2 if a is singular
c
      implicit real*8  (a - h, o - z)
      implicit integer (i - n)
c
      dimension a(ia,n)
c
      iflag = 0
c
      if (n .le. 0 .or. m + l .ge. ia) then
				       iflag = 1
				       return
      end if
c
      mp1 = m + 1
      mp2 = m + 2
c
      if (n .gt. 1) then
c
           do 1 j = 2, n
           jm1 = j - 1
           if ( a(mp1, jm1) .eq. 0.0D0 ) then
	                                 iflag = 2
                                         return
           end if
           a(mp1, jm1) = 1.0D0 / a(mp1, jm1)
           jp = max(1, j -  m )
           lp = min(l, n - jm1)
           call sscal (    lp,     a(mp1, jm1),     a(mp2     , jm1), 1)
           call bglsrd(n,l,jp,jm1, a(mp2, jp ), ia, a(mp1-j+jp, j  )   )
    1      continue
           return
      end if
c
      if (n .eq. 1) then
c
           if ( a(mp1,   n) .eq. 0.0D0 ) then
	                                 iflag = 2
                                         return
           end if
           a(mp1, n) = 1.0D0 / a(mp1, n)
           return
      end if
c
      end
      subroutine bglsrd(n, jm1, m, k, a, ia, b)
c
      return
      end
      subroutine bglssl(n, l, m, a, ia, y, iflag)
c
c     purpose: to solve the general banded linear system
c              of equations a*x=b, where a has been previously
c              factored into a=lu.
c
c     description:
c
c     this routine solves the general banded linear system of
c     equations a*x=b, using the lu factorization of a done by
c     the routine bglsdc.  the upper and lower bandwidths must
c     each be less than or equal to 64.  no pivoting is performed.
c     therefore, there are dangers if the matrix does not have
c     sufficient diagonal dominance.
c
c     on entry:
c        n = the order of the system
c        l = the number of lower diagonals  l .le. 64
c        m = the number of upper diagonals  m .le. 64
c        a = lu, where l and u are the factorization from
c            routine bglsdc
c       ia = the first dimension of a      ia .gt. m + l
c        y = the right hand side of the equations
c
c     on return:
c         y = the solution to the system of equations
c     iflag = 0 if ok
c           = 1 if (n .le. 0 .or. m + l .ge. ia)
c
      implicit real*8  (a - h, o - z)
      implicit integer (i - n)
c
      dimension a(ia, n), y(n)
c
      iflag = 0
c
      if (n .le. 0 .or. m + l .ge. ia) then
				       iflag = 1
                                       return
      end if
c
      mp1 = m + 1
      mp2 = m + 2
c
      if (n .gt. 1) then
                      call bglsrd(n,  l, 1, n-1, a(mp2, 1),  ia, y(1))
                      call bglsrd(n, -m, 1, n-1, a(m  , n), -ia, y(n))
		      return
      end if
c
      if (n .eq. 1) then
                      y(1) = y(1) * a(mp1, 1)
                      return
      end if
c
      end
      subroutine sgefad(a, lda, n, info)
c
c***********************************************************************
c     factors a real matrix by gaussian elimination with no pivoting
c***********************************************************************
c
      include 'titan.imp'
      dimension a(lda,1)
c
c-----------------------------------------------------------------------
c on entry:
c        a      real(lda, n): matrix to be factored.
c        lda    integer: leading dimension of array a.
c        n      integer: order of matrix a.
c
c on return:
c        a      factorization a = l*u, where l is product of permutation
c               and unit lower triangular matrices and u is upper 
c               triangular.
c        info   integer:
c                = 0  normal value.
c                = k  if u(k,k) .eq. 0.0.
c                not an error condition for this subroutine, but 
c                indicates that sgesl will divide by zero if called. 
c-----------------------------------------------------------------------
c
      info = 0
      nm1  = n - 1
      if (nm1 .lt. 1) then
                      info = nm1
                      return
      end if
c
      do 60 k = 1, nm1
      kp1 = k + 1
c
c     zero pivot implies this column already triangularized
c
      if (a(k, k) .eq. 0.0D0) then
                              info = k
                              return
      end if
c
c     compute multipliers
c
      t = - 1.0D0 / a(k, k)
      call sscal( n-k, t, a(k+1, k), 1 )
c
c     row elimination with column indexing
c
      do 30 j = kp1, n
      t = a(k, j)
      call saxpy( n-k, t, a(k+1, k), 1, a(k+1, j), 1 )
   30 continue
   60 continue
c
      return
      end
      subroutine sgefaj(a, ia, n, info)
c
c***********************************************************************
c     sgefaj computes the l-u decomposition of a general real matrix   c
c     using gaussian elimination. no pivoting. adapted from code by    c
c     t. l. jordan                                                     c
c***********************************************************************
c
      include 'titan.imp'
      dimension a(ia, n)
c
      info = 0
      if (n .eq. 1) then
                    info = 1
                    return
      end if
c
c     perform reduction on each column of a
c
      do 103 j = 2, n
c
c     scale the preceding column
c
            jm1 = j - 1
      if (a(jm1, jm1) .ne. 0.0E0) then
c
                   call sscal(n - jm1, 1.0E0/a(jm1, jm1), a(j, jm1), 1)
c
c                  reduce the jth column of a
c
                   call trslbl(n, jm1, a(1, 1), ia, 1, a(1, j))
      end if
  103 continue
c
      info = jm1
c
c     check for singular matrix
c
      if (a(n,n) .ne. 0.0E0) return
      info = n
c
      return
      end
      subroutine sgesld(a, lda, n, b, job)
c
c***********************************************************************
c     solves the real system a*x = b or trans(a)*x = b using the factors
c     computed by sgefa. no pivoting.
c***********************************************************************
c
      include 'titan.imp'
      dimension a(lda, 1), b(1)
c
c----------------------------------------------------------------------
c     on entry:
c        a       real(lda, n), output from sgefa
c        lda     integer, leading dimension of array a
c        n       integer, the order of the matrix a
c        b       real(n), right hand side vector 
c        job     integer
c                = 0         solve a*x = b
c                = nonzero   solve trans(a)*x = b
c
c     on return:
c        b       solution vector  x 
c
c     error condition
c        division by zero will occur if input factor contains a
c        zero on the diagonal 
c----------------------------------------------------------------------
c
      nm1 = n - 1
      if (job .eq. 0) then
c
c        job = zero , solve  a * x = b. first solve  l*y = b
c
         if (nm1 .ge. 1) then
            do 20 k = 1, nm1
               do 15 ili = 1, n-k
                  b(k+ili) = b(k+ili) + b(k) * a(k+ili, k)
   15          continue
   20       continue
         end if
c
c        now solve  u*x = y
c
         do 40 kb = 1, n
            k = n + 1 - kb
            b(k) = b(k) / a(k, k)
            do 35 ili = 1, k-1
               b(ili) = b(ili) - b(k) * a(ili, k)
   35       continue
   40    continue
      return
      end if
c
c        job = nonzero, solve trans(a) * x = b
c        first solve trans(u)*y = b
c
      if (job .ne. 0) then
         do 60 k = 1, n
            b(k) = b(k) / a(k,k)
cdir$ ivdep
            do 60 ili = k+1, n
               b(ili) = b(ili) - a(k, ili) * b(k)
   60       continue
c
c        now solve trans(l)*x = y
c
         if (nm1 .ge. 1) then
            do 80 kb = 1, nm1
               k = n - kb
               b(k) = b(k) + sdot( n-k, a(k+1, k), 1, b(k+1), 1)
   80       continue
         end if
c
      return
      end if
c
      end
      subroutine sgeslj( a, ia, n, y )
c
c***********************************************************************
c     sgeslj computes the solution of a general real system which has  c
c     been decomposed by sgefaj. no pivoting. adapted from code by     c
c     t. l. jordan                                                     c
c***********************************************************************
c
      include 'titan.imp'
      dimension a(ia, 1), y(1)
c
c==================================================================
c     solution x=y/a
c==================================================================
c
      if (n .eq. 1) then
                    y(1) = y(1) / a(1,1)
                    return
      endif
c
c==================================================================
c     solve ax=y
c==================================================================
c
c     perform the forward substitution
c
      call trslbl ( n, n - 1, a(1, 1),  ia,  1, y(1) )
c
c     perform the back substitution
c
      call trslbl ( n, n    , a(n, n), -ia, -1, y(n) )
c
      return
      end
      subroutine sgeslmj(a, ia, n, b, ib, lot)
c
c***********************************************************************
c     solve a system of linear equations with many right hand sides.   c
c     no pivoting. adapted from code by t.l. jordan. modified to use   c
c     standard linpack sgemv calling sequence                          c
c***********************************************************************
c
c     description:
c     solves real system of equations a*x=b, where x and b are matrices.
c     lu  factorization of routine sgefaj assumed. no pivoting performed
c     internally.
c
c     on entry:
c        a = output from sgefaj
c       ia = leading dimension of array a
c        n = order of matrix a
c        b = matrix of right hand sides
c       ib = leading dimension of the array b
c      lot = number of right hand sides
c
c     on return
c        b = solution matrix x
c
c***********************************************************************
c
      include 'titan.imp'
      dimension a(ia,1), b(ib,1)
c
c=======================================================================
c
      if (lot .lt. 20 .and. n .gt. 4*lot) then
c
c     for small number of long rhs's use sgeslj
c
	    do 10 j = 1, lot
	    call sgeslj(a, ia, n, b(1, j))
   10       continue
      else
c
c     for large number of rhs's use sgemv to solve l*y = b
c
            do 20 k = 1, n - 1
            call sgemv('t',   k, lot, -1.0D0, b(  1, 1), ib, a(k+1, 1),
     .                                         ia, 1.0D0, b(k+1, 1), ib)
   20       continue
c
c     then solve u*x = y
c
            do 30 k = n, 1, -1
            call sgemv('t', n-k, lot, -1.0D0, b(k+1, 1), ib, a(k, k+1),
     .                                         ia, 1.0D0, b(k  , 1), ib)
            call sscal(lot, 1.0D0 / a(k, k), b(k, 1), ib)
   30       continue
      end if
c
      return
      end
