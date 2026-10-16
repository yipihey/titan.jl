c***********************************************************************
c          cray vector merges
c***********************************************************************
      real*8 function cvmgp(a, b, c)
      implicit real*8 (a-h, o-z)
                       cvmgp = a
      if(c .lt. 0.0d0) cvmgp = b
      return
      end
      real*8 function cvmgm(a, b, c)
      implicit real*8 (a-h, o-z)
                       cvmgm = a
      if(c .ge. 0.0d0) cvmgm = b
      return
      end
      real*8 function cvmgz(a, b, c)
      implicit real*8 (a-h, o-z)
                       cvmgz = a
      if(c .ne. 0.0d0) cvmgz = b
      return
      end
      real*8 function cvmgn(a, b, c)
      implicit real*8 (a-h, o-z)
                       cvmgn = a
      if(c .eq. 0.0d0) cvmgn = b
      return
      end
      real*8 function cvmgt(a, b, c)
      implicit real*8 (a-h, o-z)
      logical c
              cvmgt = b
      if( c ) cvmgt = a
      return
      end
      subroutine trslbl(n, jm1, a, ia, m, b)
c
      return
      end
      real*8 function ssum(n,sx,incx)
      implicit real*8 (a-h, o-z)
c
c     takes the sum of the values.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
c
      real*8 sx(1),stemp
      integer i,incx,m,mp1,n,nincx
c
      ssum = 0.0d0
      stemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        stemp = stemp + sx(i)
   10 continue
      ssum = stemp
      return
c
c        code for increment equal to 1
c
c
c        clean-up loop
c
   20 m = mod(n,6)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = stemp + sx(i)
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        stemp = stemp + sx(i) + sx(i + 1) + sx(i + 2)
     *  + sx(i + 3) + sx(i + 4) + sx(i + 5)
   50 continue
   60 ssum = stemp
      return
      end
      integer function ismax (n,sx,incx)
c
      implicit real*8 (a-h, o-z)
c
c***********************************************************************
c     finds the index of element having maximum value.
c***********************************************************************
c
      real*8 sx(1),smax
      integer i,incx,ix,n
c
      ismax = 0
      if (n .lt. 1) return
      ismax = 1
      if (n .eq. 1) return
c
      if (incx .ne. 1) then
c
c          code for increment not equal to 1
c
           ix = 1
           smax = sx(1)
           ix = ix + incx
           do 10 i = 2, n
           if ( sx(ix) .le. smax ) go to 5
           ismax = i
           smax = sx(ix)
    5      ix = ix + incx
   10      continue
           return
      endif
c
      if (incx .eq. 1) then
c
c          code for increment equal to 1
c
           smax = sx(1)
           do 30 i = 2, n
           if ( sx(i) .le. smax ) go to 30
           ismax = i
           smax = sx(i)
   30      continue
           return
      endif
c
      end
      integer function ismin (n,sx,incx)
c
      implicit real*8 (a-h, o-z)
c
c***********************************************************************
c     finds the index of element having minimum value.
c***********************************************************************
c
      real*8 sx(1),smin
      integer i,incx,ix,n
c
      ismin = 0
      if (n .lt. 1) return
      ismin = 1
      if (n .eq. 1) return
c
      if (incx .ne. 1) then
c
c          code for increment not equal to 1
c
           ix = 1
           smin = sx(1)
           ix = ix + incx
           do 10 i = 2, n
           if ( sx(ix) .ge. smin ) go to 5
           ismin = i
           smin = sx(ix)
    5      ix = ix + incx
   10      continue
           return
      endif
c
      if (incx .eq. 1) then
c
c          code for increment equal to 1
c
           smin = sx(1)
           do 30 i = 2, n
           if ( sx(i) .ge. smin ) go to 30
           ismin = i
           smin = sx(i)
   30      continue
           return
      endif
c
      end
