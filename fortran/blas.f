c***********************************************************************
c                blas1
c************************************************************************
      real*8 function snrm2( n, cx, incx)
      implicit real*8 (a-h, o-z)
c
        integer n,incx
        real*8 cx(*)
c
c  find the euclidean norm of a vector of reals. do it in a way
c  to help avoid overflow and underflow problems.
c
c  input:
c    n          number of elements.
c    incx       increment between elements.
c    cx         real vector.
c
c  output:
c    snrm2      the euclidean norm.
c
c------------------------------------------------------------------------
        integer i
        real*8 scal,maxv,sum
c
c  externals.
c
        integer isamax
c
        if(n.le.0)then
          snrm2 = 0
        else
          i = isamax(n,cx,incx)
          i = (i-1)*incx + 1
          maxv = abs( cx(i) )
          scal = 1./maxv
          sum = 0
          do 10 i=1,n*incx,incx
            sum = sum + ( scal*cx(i) )**2
   10     continue
          snrm2 = maxv * sqrt( sum )
        endif
        end
c***********************************************************************
      integer function isamax(n,sx,incx)
      implicit real*8 (a-h, o-z)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
c-
      real*8 sx(1),smax
      integer i,incx,ix,n
c
      isamax = 0
      if( n .lt. 1 ) return
      isamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = abs(sx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(abs(sx(ix)).le.smax) go to 5
         isamax = i
         smax = abs(sx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = abs(sx(1))
      do 30 i = 2,n
         if(abs(sx(i)).le.smax) go to 30
         isamax = i
         smax = abs(sx(i))
   30 continue
      return
      end
c***********************************************************************
      real*8 function sasum(n,sx,incx)
      implicit real*8 (a-h, o-z)
c
c     takes the sum of the absolute values.
c     uses unrolled loops for increment equal to one.
c     jack dongarra, linpack, 3/11/78.
c
c-
      real*8 sx(1),stemp
      integer i,incx,m,mp1,n,nincx
c
      sasum = 0.0d0
      stemp = 0.0d0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        stemp = stemp + abs(sx(i))
   10 continue
      sasum = stemp
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
        stemp = stemp + abs(sx(i))
   30 continue
      if( n .lt. 6 ) go to 60
   40 mp1 = m + 1
      do 50 i = mp1,n,6
        stemp = stemp + abs(sx(i)) + abs(sx(i + 1)) + abs(sx(i + 2))
     *  + abs(sx(i + 3)) + abs(sx(i + 4)) + abs(sx(i + 5))
   50 continue
   60 sasum = stemp
      return
      end
c***********************************************************************
      subroutine saxpy(n,sa,sx,incx,sy,incy)
      implicit real*8 (a-h, o-z)
c
c     constant times a vector plus a vector.
c     uses unrolled loop for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c
c-
      real*8 sx(1),sy(1),sa
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if (sa .eq. 0.0d0) return
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
c***********************************************************************
      subroutine  scopy(n,sx,incx,sy,incy)
      implicit real*8 (a-h, o-z)
c
c     copies a vector, x, to a vector, y.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c-
c
      real*8 sx(1),sy(1)
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
c        code for both increments equal to 1
c
c
c        clean-up loop
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
c***********************************************************************
      real*8 function sdot(n,sx,incx,sy,incy)
      implicit real*8 (a-h, o-z)
c
c     forms the dot product of two vectors.
c     uses unrolled loops for increments equal to one.
c     jack dongarra, linpack, 3/11/78.
c-
c
      real*8 sx(1),sy(1),stemp
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
      do 50 i = mp1,n,5
        stemp = stemp + sx(i)*sy(i) + sx(i + 1)*sy(i + 1) +
     *   sx(i + 2)*sy(i + 2) + sx(i + 3)*sy(i + 3) + sx(i + 4)*sy(i + 4)
   50 continue
   60 sdot = stemp
      return
      end
c***********************************************************************
      real*8 function smach(job)
      implicit real*8 (a-h, o-z)
      integer job
c
c     smach computes machine parameters of floating point
c     arithmetic for use in testing only.  not required by
c     linpack proper.
c
c     if trouble with automatic computation of these quantities,
c     they can be set by direct assignment statements.
c     assume the computer has
c
c        b = base of arithmetic
c        t = number of base  b  digits
c        l = smallest possible exponent
c        u = largest possible exponent
c
c     then
c
c        eps = b**(1-t)
c        tiny = 100.0*b**(-l+t)
c        huge = 0.01*b**(u-t)
c
c     job is 1, 2 or 3 for epsilon, tiny and huge, respectively.
c
c
      real*8 eps,tiny,huge,s
c
      eps = 1.0d0
   10 eps = eps/2.0d0
      s = 1.0d0 + eps
      if (s .gt. 1.0d0) go to 10
      eps = 2.0d0*eps
c
      s = 1.0d0
   20 tiny = s
      s = s/16.0d0
      if (s*100.0d0 .ne. 0.0d0) go to 20
      tiny = (tiny/eps)*100.0d0
      huge = 1.0d0/tiny
c
      if (job .eq. 1) smach = eps
      if (job .eq. 2) smach = tiny
      if (job .eq. 3) smach = huge
      return
      end
c***********************************************************************
      subroutine  srot (n,sx,incx,sy,incy,c,s)
      implicit real*8 (a-h, o-z)
c
c     applies a plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
c-
      real*8 sx(1),sy(1),stemp,c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = c*sx(ix) + s*sy(iy)
        sy(iy) = c*sy(iy) - s*sx(ix)
        sx(ix) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        stemp = c*sx(i) + s*sy(i)
        sy(i) = c*sy(i) - s*sx(i)
        sx(i) = stemp
   30 continue
      return
      end
c***********************************************************************
      subroutine  sscal(n,sa,sx,incx)
      implicit real*8 (a-h, o-z)
c
c     scales a vector by a constant.
c     uses unrolled loops for increment equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
c-
      real*8 sa,sx(1)
      integer i,incx,m,mp1,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        sx(i) = sa*sx(i)
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
c***********************************************************************
      subroutine  sswap (n,sx,incx,sy,incy)
      implicit real*8 (a-h, o-z)
c
c     interchanges two vectors.
c     uses unrolled loops for increments equal to 1.
c     jack dongarra, linpack, 3/11/78.
c
c-
      real*8 sx(1),sy(1),stemp
      integer i,incx,incy,ix,iy,m,mp1,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        stemp = sx(ix)
        sx(ix) = sy(iy)
        sy(iy) = stemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
c
c       clean-up loop
c
   20 m = mod(n,3)
      if( m .eq. 0 ) go to 40
      do 30 i = 1,m
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
   30 continue
      if( n .lt. 3 ) return
   40 mp1 = m + 1
      do 50 i = mp1,n,3
        stemp = sx(i)
        sx(i) = sy(i)
        sy(i) = stemp
        stemp = sx(i + 1)
        sx(i + 1) = sy(i + 1)
        sy(i + 1) = stemp
        stemp = sx(i + 2)
        sx(i + 2) = sy(i + 2)
        sy(i + 2) = stemp
   50 continue
      return
      end
c***********************************************************************
      subroutine srotg(sa,sb,c,s)
      implicit real*8 (a-h, o-z)
c
c     construct givens plane rotation.
c     jack dongarra, linpack, 3/11/78.
c
c-
      real*8 sa,sb,c,s,roe,scale,r,z
c
      roe = sb
      if( abs(sa) .gt. abs(sb) ) roe = sa
      scale = abs(sa) + abs(sb)
      if( scale .ne. 0.0d0 ) go to 10
         c = 1.0d0
         s = 0.0d0
         r = 0.0d0
         go to 20
   10 r = scale*sqrt((sa/scale)**2 + (sb/scale)**2)
      r = sign(1.0d0,roe)*r
      c = sa/r
      s = sb/r
   20 z = 1.0d0
      if( abs(sa) .gt. abs(sb) ) z = s
      if( abs(sb) .ge. abs(sa) .and. c .ne. 0.0d0 ) z = 1.0d0/c
      sa = r
      sb = z
      return
      end
c***********************************************************************
      real function scnrm2(n,cx,incx)
c
        integer n,incx
        complex cx(*)
c
c  find the euclidean norm of a complex vector.
c
c------------------------------------------------------------------------
        integer i
        real sum,maxv,scal
        complex temp
c
c  externals.
c
        integer icamax
c
        if(n.le.0)then
          scnrm2 = 0
        else
          i = icamax(n, cx, incx)
          i = (i-1)*incx + 1
          maxv = abs(cx(i))
          scal = 1/maxv
          sum = 0
          do 10 i=1,n*incx,incx
            temp = scal * cx(i)
            sum = sum + real(temp)**2 + aimag(temp)**2
   10     continue
          scnrm2 = maxv * sqrt(sum)
        endif
        end
c***********************************************************************
      subroutine caxpy(n,ca,cx,incx,cy,incy)
c
c     constant times a vector plus a vector.
c     jack dongarra, linpack, 3/11/78.
c
c-
      complex cx(1),cy(1),ca
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if (abs(real(ca)) + abs(aimag(ca)) .eq. 0.0 ) return
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
        cy(iy) = cy(iy) + ca*cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        cy(i) = cy(i) + ca*cx(i)
   30 continue
      return
      end
c***********************************************************************
      subroutine  ccopy(n,cx,incx,cy,incy)
c
c     copies a vector, x, to a vector, y.
c     jack dongarra, linpack, 3/11/78.
c
c-
      complex cx(1),cy(1)
      integer i,incx,incy,ix,iy,n
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
        cy(iy) = cx(ix)
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        cy(i) = cx(i)
   30 continue
      return
      end
c***********************************************************************
      complex function cdotc(n,cx,incx,cy,incy)
c
c     forms the dot product of two vectors, conjugating the first
c     vector.
c     jack dongarra, linpack,  3/11/78.
c
c-
      complex cx(1),cy(1),ctemp
      integer i,incx,incy,ix,iy,n
c
      ctemp = (0.0,0.0)
      cdotc = (0.0,0.0)
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
        ctemp = ctemp + conjg(cx(ix))*cy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      cdotc = ctemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ctemp = ctemp + conjg(cx(i))*cy(i)
   30 continue
      cdotc = ctemp
      return
      end
c***********************************************************************
      complex function cdotu(n,cx,incx,cy,incy)
c
c     forms the dot product of two vectors.
c     jack dongarra, linpack, 3/11/78.
c
c-
      complex cx(1),cy(1),ctemp
      integer i,incx,incy,ix,iy,n
c
      ctemp = (0.0,0.0)
      cdotu = (0.0,0.0)
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
        ctemp = ctemp + cx(ix)*cy(iy)
        ix = ix + incx
        iy = iy + incy
   10 continue
      cdotu = ctemp
      return
c
c        code for both increments equal to 1
c
   20 do 30 i = 1,n
        ctemp = ctemp + cx(i)*cy(i)
   30 continue
      cdotu = ctemp
      return
      end
c***********************************************************************
      real function cmach(job)
      integer job
c
c     smach computes machine parameters of floating point
c     arithmetic for use in testing only.  not required by
c     linpack proper.
c
c     if trouble with automatic computation of these quantities,
c     they can be set by direct assignment statements.
c     assume the computer has
c
c        b = base of arithmetic
c        t = number of base  b  digits
c        l = smallest possible exponent
c        u = largest possible exponent
c
c     then
c
c        eps = b**(1-t)
c        tiny = 100.0*b**(-l+t)
c        huge = 0.01*b**(u-t)
c
c     cmach same as smach except if complex division
c     is done by
c
c        1/(x+i*y) = (x-i*y)/(x**2+y**2)
c
c     then
c
c        tiny = sqrt(tiny)
c        huge = sqrt(huge)
c
c
c     job is 1, 2 or 3 for epsilon, tiny and huge, respectively.
c
c
      real eps,tiny,huge,s
c
c-
      eps = 1.0
   10 eps = eps/2.0
      s = 1.0 + eps
      if (s .gt. 1.0) go to 10
      eps = 2.0*eps
      cmach =eps
      if( job .eq. 1) return
c
      s = 1.0
   20 tiny = s
      s = s/16.0
      if (s*1.0 .ne. 0.0) go to 20
      tiny = (tiny/eps)*100.
      s = real((1.0,0.0)/cmplx(tiny,0.0))
      if (s .ne. 1.0/tiny) tiny = sqrt(tiny)
      huge = 1.0/tiny
      if (job .eq. 1) cmach = eps
      if (job .eq. 2) cmach = tiny
      if (job .eq. 3) cmach = huge
      return
      end
c***********************************************************************
      subroutine crotg(ca,cb,c,s)
      complex ca,cb,s
c-
      real c
      real norm,scale
      complex alpha
      if (cabs(ca) .ne. 0.) go to 10
         c = 0.
         s = (1.,0.)
         ca = cb
         go to 20
   10 continue
         scale = cabs(ca) + cabs(cb)
         norm = scale * sqrt((cabs(ca/scale))**2 + (cabs(cb/scale))**2)
         alpha = ca /cabs(ca)
         c = cabs(ca) / norm
         s = alpha * conjg(cb) / norm
         ca = alpha * norm
   20 continue
      return
      end
c***********************************************************************
      subroutine  cscal(n,ca,cx,incx)
c
c     scales a vector by a constant.
c     jack dongarra, linpack,  3/11/78.
c
c-
      complex ca,cx(1)
      integer i,incx,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = ca*cx(i)
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        cx(i) = ca*cx(i)
   30 continue
      return
      end
c***********************************************************************
      subroutine  csrot (n,cx,incx,cy,incy,c,s)
c
c     applies a plane rotation, where the cos and sin (c and s) are real
c     and the vectors cx and cy are complex.
c     jack dongarra, linpack, 3/11/78.
c
c-
      complex cx(1),cy(1),ctemp
      real c,s
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = c*cx(ix) + s*cy(iy)
        cy(iy) = c*cy(iy) - s*cx(ix)
        cx(ix) = ctemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
c
   20 do 30 i = 1,n
        ctemp = c*cx(i) + s*cy(i)
        cy(i) = c*cy(i) - s*cx(i)
        cx(i) = ctemp
   30 continue
      return
      end
c***********************************************************************
      subroutine  csscal(n,sa,cx,incx)
c
c     scales a complex vector by a real constant.
c     jack dongarra, linpack, 3/11/78.
c
c-
      complex cx(1)
      real sa
      integer i,incx,n,nincx
c
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i)))
   10 continue
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        cx(i) = cmplx(sa*real(cx(i)),sa*aimag(cx(i)))
   30 continue
      return
      end
c***********************************************************************
      subroutine  cswap (n,cx,incx,cy,incy)
c
c     interchanges two vectors.
c     jack dongarra, linpack, 3/11/78.
c
c-
      complex cx(1),cy(1),ctemp
      integer i,incx,incy,ix,iy,n
c
      if(n.le.0)return
      if(incx.eq.1.and.incy.eq.1)go to 20
c
c       code for unequal increments or equal increments not equal
c         to 1
c
      ix = 1
      iy = 1
      if(incx.lt.0)ix = (-n+1)*incx + 1
      if(incy.lt.0)iy = (-n+1)*incy + 1
      do 10 i = 1,n
        ctemp = cx(ix)
        cx(ix) = cy(iy)
        cy(iy) = ctemp
        ix = ix + incx
        iy = iy + incy
   10 continue
      return
c
c       code for both increments equal to 1
   20 do 30 i = 1,n
        ctemp = cx(i)
        cx(i) = cy(i)
        cy(i) = ctemp
   30 continue
      return
      end
c***********************************************************************
      integer function icamax(n,cx,incx)
c
c     finds the index of element having max. absolute value.
c     jack dongarra, linpack, 3/11/78.
c
c-
      complex cx(1)
      real smax
      integer i,incx,ix,n
      complex zdum
      real cabs1
      cabs1(zdum) = abs(real(zdum)) + abs(aimag(zdum))
c
      icamax = 0
      if( n .lt. 1 ) return
      icamax = 1
      if(n.eq.1)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      ix = 1
      smax = cabs1(cx(1))
      ix = ix + incx
      do 10 i = 2,n
         if(cabs1(cx(ix)).le.smax) go to 5
         icamax = i
         smax = cabs1(cx(ix))
    5    ix = ix + incx
   10 continue
      return
c
c        code for increment equal to 1
c
   20 smax = cabs1(cx(1))
      do 30 i = 2,n
         if(cabs1(cx(i)).le.smax) go to 30
         icamax = i
         smax = cabs1(cx(i))
   30 continue
      return
      end
c***********************************************************************
      real function scasum(n,cx,incx)
c
c     takes the sum of the absolute values of a complex vector and
c     returns a single precision result.
c     jack dongarra, linpack, 3/11/78.
c
c-
      complex cx(1)
      real stemp
      integer i,incx,n,nincx
c
      scasum = 0.0e0
      stemp = 0.0e0
      if(n.le.0)return
      if(incx.eq.1)go to 20
c
c        code for increment not equal to 1
c
      nincx = n*incx
      do 10 i = 1,nincx,incx
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
   10 continue
      scasum = stemp
      return
c
c        code for increment equal to 1
c
   20 do 30 i = 1,n
        stemp = stemp + abs(real(cx(i))) + abs(aimag(cx(i)))
   30 continue
      scasum = stemp
      return
      end
c***********************************************************************
c                             blas2
c***********************************************************************
c
c     see:
c
c        dongarra j. j., du croz j. j., hammarling s.  and hanson r. j..
c        an  extended  set of fortran  basic linear algebra subprograms.
c
c        technical  memoranda  nos. 41 (revision 3) and 81,  mathematics
c        and  computer science  division,  argonne  national laboratory,
c        9700 south cass avenue, argonne, illinois 60439, us.
c
c        or
c
c        nag  technical reports tr3/87 and tr4/87,  numerical algorithms
c        group  ltd.,  nag  central  office,  256  banbury  road, oxford
c        ox2 7de, uk,  and  numerical algorithms group inc.,  1101  31st
c        street,  suite 100,  downers grove,  illinois 60515-1263,  usa.
c
c***********************************************************************
c
      subroutine sgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx,
     $                   beta, y, incy )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha, beta
      integer            incx, incy, kl, ku, lda, m, n
      character*1        trans
c     .. array arguments ..
      real*8             a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  sgbmv  performs one of the matrix-vector operations
c
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
c
c  parameters
c  ==========
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
c
c              trans = 't' or 't'   y := alpha*a'*x + beta*y.
c
c              trans = 'c' or 'c'   y := alpha*a'*x + beta*y.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  kl     - integer.
c           on entry, kl specifies the number of sub-diagonals of the
c           matrix a. kl must satisfy  0 .le. kl.
c           unchanged on exit.
c
c  ku     - integer.
c           on entry, ku specifies the number of super-diagonals of the
c           matrix a. ku must satisfy  0 .le. ku.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry, the leading ( kl + ku + 1 ) by n part of the
c           array a must contain the matrix of coefficients, supplied
c           column by column, with the leading diagonal of the matrix in
c           row ( ku + 1 ) of the array, the first super-diagonal
c           starting at position 2 in row ku, the first sub-diagonal
c           starting at position 1 in row ( ku + 2 ), and so on.
c           elements in the array a that do not correspond to elements
c           in the band matrix (such as the top left ku by ku triangle)
c           are not referenced.
c           the following program segment will transfer a band matrix
c           from conventional full matrix storage to band storage:
c
c                 do 20, j = 1, n
c                    k = ku + 1 - j
c                    do 10, i = max( 1, j - ku ), min( m, j + kl )
c                       a( k + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( kl + ku + 1 ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
c           before entry, the incremented array y must contain the
c           vector y. on exit, y is overwritten by the updated vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, iy, j, jx, jy, k, kup1, kx, ky,
     $                   lenx, leny
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max, min
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( kl.lt.0 )then
         info = 4
      else if( ku.lt.0 )then
         info = 5
      else if( lda.lt.( kl + ku + 1 ) )then
         info = 8
      else if( incx.eq.0 )then
         info = 10
      else if( incy.eq.0 )then
         info = 13
      end if
      if( info.ne.0 )then
         call xerbla( 'sgbmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set  lenx  and  leny, the lengths of the vectors x and y, and set
c     up the start points in  x  and  y.
c
      if( lsame( trans, 'n' ) )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the band part of a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      kup1 = ku + 1
      if( lsame( trans, 'n' ) )then
c
c        form  y := alpha*a*x + y.
c
         jx = kx
         if( incy.eq.1 )then
            do 60, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  k    = kup1 - j
                  do 50, i = max( 1, j - ku ), min( m, j + kl )
                     y( i ) = y( i ) + temp*a( k + i, j )
   50             continue
               end if
               jx = jx + incx
   60       continue
         else
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy   = ky
                  k    = kup1 - j
                  do 70, i = max( 1, j - ku ), min( m, j + kl )
                     y( iy ) = y( iy ) + temp*a( k + i, j )
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
               if( j.gt.ku )
     $            ky = ky + incy
   80       continue
         end if
      else
c
c        form  y := alpha*a'*x + y.
c
         jy = ky
         if( incx.eq.1 )then
            do 100, j = 1, n
               temp = zero
               k    = kup1 - j
               do 90, i = max( 1, j - ku ), min( m, j + kl )
                  temp = temp + a( k + i, j )*x( i )
   90          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  100       continue
         else
            do 120, j = 1, n
               temp = zero
               ix   = kx
               k    = kup1 - j
               do 110, i = max( 1, j - ku ), min( m, j + kl )
                  temp = temp + a( k + i, j )*x( ix )
                  ix   = ix   + incx
  110          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
               if( j.gt.ku )
     $            kx = kx + incx
  120       continue
         end if
      end if
c
      return
c
c     end of sgbmv .
c
      end
c
c***********************************************************************
c
      subroutine sgemv ( trans, m, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha, beta
      integer            incx, incy, lda, m, n
      character*1        trans
c     .. array arguments ..
      real*8             a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  sgemv  performs one of the matrix-vector operations
c
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n matrix.
c
c  parameters
c  ==========
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
c
c              trans = 't' or 't'   y := alpha*a'*x + beta*y.
c
c              trans = 'c' or 'c'   y := alpha*a'*x + beta*y.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
c           before entry with beta non-zero, the incremented array y
c           must contain the vector y. on exit, y is overwritten by the
c           updated vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( lda.lt.max( 1, m ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'sgemv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set  lenx  and  leny, the lengths of the vectors x and y, and set
c     up the start points in  x  and  y.
c
      if( lsame( trans, 'n' ) )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( lsame( trans, 'n' ) )then
c
c        form  y := alpha*a*x + y.
c
         jx = kx
         if( incy.eq.1 )then
            do 60, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  do 50, i = 1, m
                     y( i ) = y( i ) + temp*a( i, j )
   50             continue
               end if
               jx = jx + incx
   60       continue
         else
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy   = ky
                  do 70, i = 1, m
                     y( iy ) = y( iy ) + temp*a( i, j )
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
   80       continue
         end if
      else
c
c        form  y := alpha*a'*x + y.
c
         jy = ky
         if( incx.eq.1 )then
            do 100, j = 1, n
               temp = zero
               do 90, i = 1, m
                  temp = temp + a( i, j )*x( i )
   90          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  100       continue
         else
            do 120, j = 1, n
               temp = zero
               ix   = kx
               do 110, i = 1, m
                  temp = temp + a( i, j )*x( ix )
                  ix   = ix   + incx
  110          continue
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  120       continue
         end if
      end if
c
      return
c
c     end of sgemv .
c
      end
c
c***********************************************************************
c
      subroutine sger  ( m, n, alpha, x, incx, y, incy, a, lda )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha
      integer            incx, incy, lda, m, n
c     .. array arguments ..
      real*8             a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  sger   performs the rank 1 operation
c
c     a := alpha*x*y' + a,
c
c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and a is an m by n matrix.
c
c  parameters
c  ==========
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( m - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the m
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y.
c           unchanged on exit.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients. on exit, a is
c           overwritten by the updated matrix.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, j, jy, kx
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( m.lt.0 )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      else if( lda.lt.max( 1, m ) )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'sger  ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( incy.gt.0 )then
         jy = 1
      else
         jy = 1 - ( n - 1 )*incy
      end if
      if( incx.eq.1 )then
         do 20, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               do 10, i = 1, m
                  a( i, j ) = a( i, j ) + x( i )*temp
   10          continue
            end if
            jy = jy + incy
   20    continue
      else
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( m - 1 )*incx
         end if
         do 40, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               ix   = kx
               do 30, i = 1, m
                  a( i, j ) = a( i, j ) + x( ix )*temp
                  ix        = ix        + incx
   30          continue
            end if
            jy = jy + incy
   40    continue
      end if
c
      return
c
c     end of sger  .
c
      end
c
c***********************************************************************
c
      subroutine ssbmv ( uplo, n, k, alpha, a, lda, x, incx,
     $                   beta, y, incy )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha, beta
      integer            incx, incy, k, lda, n
      character*1        uplo
c     .. array arguments ..
      real*8             a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  ssbmv  performs the matrix-vector  operation
c
c     y := alpha*a*x + beta*y,
c
c  where alpha and beta are scalars, x and y are n element vectors and
c  a is an n by n symmetric band matrix, with k super-diagonals.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the band matrix a is being supplied as
c           follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  being supplied.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  being supplied.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry, k specifies the number of super-diagonals of the
c           matrix a. k must satisfy  0 .le. k.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
c           by n part of the array a must contain the upper triangular
c           band part of the symmetric matrix, supplied column by
c           column, with the leading diagonal of the matrix in row
c           ( k + 1 ) of the array, the first super-diagonal starting at
c           position 2 in row k, and so on. the top left k by k triangle
c           of the array a is not referenced.
c           the following program segment will transfer the upper
c           triangular part of a symmetric band matrix from conventional
c           full matrix storage to band storage:
c
c                 do 20, j = 1, n
c                    m = k + 1 - j
c                    do 10, i = max( 1, j - k ), j
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           before entry with uplo = 'l' or 'l', the leading ( k + 1 )
c           by n part of the array a must contain the lower triangular
c           band part of the symmetric matrix, supplied column by
c           column, with the leading diagonal of the matrix in row 1 of
c           the array, the first sub-diagonal starting at position 1 in
c           row 2, and so on. the bottom right k by k triangle of the
c           array a is not referenced.
c           the following program segment will transfer the lower
c           triangular part of a symmetric band matrix from conventional
c           full matrix storage to band storage:
c
c                 do 20, j = 1, n
c                    m = 1 - j
c                    do 10, i = j, min( n, j + k )
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( k + 1 ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the
c           vector y. on exit, y is overwritten by the updated vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max, min
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( k.lt.0 )then
         info = 3
      else if( lda.lt.( k + 1 ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'ssbmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set up the start points in  x  and  y.
c
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( n - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( n - 1 )*incy
      end if
c
c     start the operations. in this version the elements of the array a
c     are accessed sequentially with one pass through a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, n
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, n
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, n
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, n
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( lsame( uplo, 'u' ) )then
c
c        form  y  when upper triangle of a is stored.
c
         kplus1 = k + 1
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               l     = kplus1 - j
               do 50, i = max( 1, j - k ), j - 1
                  y( i ) = y( i ) + temp1*a( l + i, j )
                  temp2  = temp2  + a( l + i, j )*x( i )
   50          continue
               y( j ) = y( j ) + temp1*a( kplus1, j ) + alpha*temp2
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               l     = kplus1 - j
               do 70, i = max( 1, j - k ), j - 1
                  y( iy ) = y( iy ) + temp1*a( l + i, j )
                  temp2   = temp2   + a( l + i, j )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*a( kplus1, j ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
               if( j.gt.k )then
                  kx = kx + incx
                  ky = ky + incy
               end if
   80       continue
         end if
      else
c
c        form  y  when lower triangle of a is stored.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j )       + temp1*a( 1, j )
               l      = 1            - j
               do 90, i = j + 1, min( n, j + k )
                  y( i ) = y( i ) + temp1*a( l + i, j )
                  temp2  = temp2  + a( l + i, j )*x( i )
   90          continue
               y( j ) = y( j ) + alpha*temp2
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy )       + temp1*a( 1, j )
               l       = 1             - j
               ix      = jx
               iy      = jy
               do 110, i = j + 1, min( n, j + k )
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*a( l + i, j )
                  temp2   = temp2   + a( l + i, j )*x( ix )
  110          continue
               y( jy ) = y( jy ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
  120       continue
         end if
      end if
c
      return
c
c     end of ssbmv .
c
      end
c
c***********************************************************************
c
      subroutine sspmv ( uplo, n, alpha, ap, x, incx, beta, y, incy )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha, beta
      integer            incx, incy, n
      character*1        uplo
c     .. array arguments ..
      real*8             ap( * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  sspmv  performs the matrix-vector operation
c
c     y := alpha*a*x + beta*y,
c
c  where alpha and beta are scalars, x and y are n element vectors and
c  a is an n by n symmetric matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the matrix a is supplied in the packed
c           array ap as follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  supplied in ap.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  supplied in ap.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  ap     - real             array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with uplo = 'u' or 'u', the array ap must
c           contain the upper triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
c           and a( 2, 2 ) respectively, and so on.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
c           and a( 3, 1 ) respectively, and so on.
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y. on exit, y is overwritten by the updated
c           vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 6
      else if( incy.eq.0 )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'sspmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set up the start points in  x  and  y.
c
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( n - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( n - 1 )*incy
      end if
c
c     start the operations. in this version the elements of the array ap
c     are accessed sequentially with one pass through ap.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, n
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, n
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, n
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, n
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      kk = 1
      if( lsame( uplo, 'u' ) )then
c
c        form  y  when ap contains the upper triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               k     = kk
               do 50, i = 1, j - 1
                  y( i ) = y( i ) + temp1*ap( k )
                  temp2  = temp2  + ap( k )*x( i )
                  k      = k      + 1
   50          continue
               y( j ) = y( j ) + temp1*ap( kk + j - 1 ) + alpha*temp2
               kk     = kk     + j
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               do 70, k = kk, kk + j - 2
                  y( iy ) = y( iy ) + temp1*ap( k )
                  temp2   = temp2   + ap( k )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*ap( kk + j - 1 ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
               kk      = kk      + j
   80       continue
         end if
      else
c
c        form  y  when ap contains the lower triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j )       + temp1*ap( kk )
               k      = kk           + 1
               do 90, i = j + 1, n
                  y( i ) = y( i ) + temp1*ap( k )
                  temp2  = temp2  + ap( k )*x( i )
                  k      = k      + 1
   90          continue
               y( j ) = y( j ) + alpha*temp2
               kk     = kk     + ( n - j + 1 )
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy )       + temp1*ap( kk )
               ix      = jx
               iy      = jy
               do 110, k = kk + 1, kk + n - j
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*ap( k )
                  temp2   = temp2   + ap( k )*x( ix )
  110          continue
               y( jy ) = y( jy ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
               kk      = kk      + ( n - j + 1 )
  120       continue
         end if
      end if
c
      return
c
c     end of sspmv .
c
      end
c
c***********************************************************************
c
      subroutine sspr  ( uplo, n, alpha, x, incx, ap )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha
      integer            incx, n
      character*1        uplo
c     .. array arguments ..
      real*8             ap( * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  sspr    performs the symmetric rank 1 operation
c
c     a := alpha*x*x' + a,
c
c  where alpha is a real scalar, x is an n element vector and a is an
c  n by n symmetric matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the matrix a is supplied in the packed
c           array ap as follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  supplied in ap.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  supplied in ap.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  ap     - real             array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
c           and a( 2, 2 ) respectively, and so on. on exit, the array
c           ap is overwritten by the upper triangular part of the
c           updated matrix.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
c           and a( 3, 1 ) respectively, and so on. on exit, the array
c           ap is overwritten by the lower triangular part of the
c           updated matrix.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, j, jx, k, kk, kx
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      end if
      if( info.ne.0 )then
         call xerbla( 'sspr  ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     set the start point in x if the increment is not unity.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of the array ap
c     are accessed sequentially with one pass through ap.
c
      kk = 1
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when upper triangle is stored in ap.
c
         if( incx.eq.1 )then
            do 20, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*x( j )
                  k    = kk
                  do 10, i = 1, j
                     ap( k ) = ap( k ) + x( i )*temp
                     k       = k       + 1
   10             continue
               end if
               kk = kk + j
   20       continue
         else
            jx = kx
            do 40, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  ix   = kx
                  do 30, k = kk, kk + j - 1
                     ap( k ) = ap( k ) + x( ix )*temp
                     ix      = ix      + incx
   30             continue
               end if
               jx = jx + incx
               kk = kk + j
   40       continue
         end if
      else
c
c        form  a  when lower triangle is stored in ap.
c
         if( incx.eq.1 )then
            do 60, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*x( j )
                  k    = kk
                  do 50, i = j, n
                     ap( k ) = ap( k ) + x( i )*temp
                     k       = k       + 1
   50             continue
               end if
               kk = kk + n - j + 1
   60       continue
         else
            jx = kx
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  ix   = jx
                  do 70, k = kk, kk + n - j
                     ap( k ) = ap( k ) + x( ix )*temp
                     ix      = ix      + incx
   70             continue
               end if
               jx = jx + incx
               kk = kk + n - j + 1
   80       continue
         end if
      end if
c
      return
c
c     end of sspr  .
c
      end
c
c***********************************************************************
c
      subroutine sspr2 ( uplo, n, alpha, x, incx, y, incy, ap )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha
      integer            incx, incy, n
      character*1        uplo
c     .. array arguments ..
      real*8             ap( * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  sspr2  performs the symmetric rank 2 operation
c
c     a := alpha*x*y' + alpha*y*x' + a,
c
c  where alpha is a scalar, x and y are n element vectors and a is an
c  n by n symmetric matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the matrix a is supplied in the packed
c           array ap as follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  supplied in ap.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  supplied in ap.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y.
c           unchanged on exit.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c  ap     - real             array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
c           and a( 2, 2 ) respectively, and so on. on exit, the array
c           ap is overwritten by the upper triangular part of the
c           updated matrix.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular part of the symmetric matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
c           and a( 3, 1 ) respectively, and so on. on exit, the array
c           ap is overwritten by the lower triangular part of the
c           updated matrix.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'sspr2 ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     set up the start points in x and y if the increments are not both
c     unity.
c
      if( ( incx.ne.1 ).or.( incy.ne.1 ) )then
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            ky = 1
         else
            ky = 1 - ( n - 1 )*incy
         end if
         jx = kx
         jy = ky
      end if
c
c     start the operations. in this version the elements of the array ap
c     are accessed sequentially with one pass through ap.
c
      kk = 1
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when upper triangle is stored in ap.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 20, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1 = alpha*y( j )
                  temp2 = alpha*x( j )
                  k     = kk
                  do 10, i = 1, j
                     ap( k ) = ap( k ) + x( i )*temp1 + y( i )*temp2
                     k       = k       + 1
   10             continue
               end if
               kk = kk + j
   20       continue
         else
            do 40, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1 = alpha*y( jy )
                  temp2 = alpha*x( jx )
                  ix    = kx
                  iy    = ky
                  do 30, k = kk, kk + j - 1
                     ap( k ) = ap( k ) + x( ix )*temp1 + y( iy )*temp2
                     ix      = ix      + incx
                     iy      = iy      + incy
   30             continue
               end if
               jx = jx + incx
               jy = jy + incy
               kk = kk + j
   40       continue
         end if
      else
c
c        form  a  when lower triangle is stored in ap.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1 = alpha*y( j )
                  temp2 = alpha*x( j )
                  k     = kk
                  do 50, i = j, n
                     ap( k ) = ap( k ) + x( i )*temp1 + y( i )*temp2
                     k       = k       + 1
   50             continue
               end if
               kk = kk + n - j + 1
   60       continue
         else
            do 80, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1 = alpha*y( jy )
                  temp2 = alpha*x( jx )
                  ix    = jx
                  iy    = jy
                  do 70, k = kk, kk + n - j
                     ap( k ) = ap( k ) + x( ix )*temp1 + y( iy )*temp2
                     ix      = ix      + incx
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
               jy = jy + incy
               kk = kk + n - j + 1
   80       continue
         end if
      end if
c
      return
c
c     end of sspr2 .
c
      end
c
c***********************************************************************
c
      subroutine ssymv ( uplo, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha, beta
      integer            incx, incy, lda, n
      character*1        uplo
c     .. array arguments ..
      real*8             a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  ssymv  performs the matrix-vector  operation
c
c     y := alpha*a*x + beta*y,
c
c  where alpha and beta are scalars, x and y are n element vectors and
c  a is an n by n symmetric matrix.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the array a is to be referenced as
c           follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of a
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of a
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular part of the symmetric matrix and the strictly
c           lower triangular part of a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular part of the symmetric matrix and the strictly
c           upper triangular part of a is not referenced.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y. on exit, y is overwritten by the updated
c           vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( lda.lt.max( 1, n ) )then
         info = 5
      else if( incx.eq.0 )then
         info = 7
      else if( incy.eq.0 )then
         info = 10
      end if
      if( info.ne.0 )then
         call xerbla( 'ssymv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set up the start points in  x  and  y.
c
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( n - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( n - 1 )*incy
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the triangular part
c     of a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, n
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, n
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, n
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, n
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( lsame( uplo, 'u' ) )then
c
c        form  y  when a is stored in upper triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               do 50, i = 1, j - 1
                  y( i ) = y( i ) + temp1*a( i, j )
                  temp2  = temp2  + a( i, j )*x( i )
   50          continue
               y( j ) = y( j ) + temp1*a( j, j ) + alpha*temp2
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               do 70, i = 1, j - 1
                  y( iy ) = y( iy ) + temp1*a( i, j )
                  temp2   = temp2   + a( i, j )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*a( j, j ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
   80       continue
         end if
      else
c
c        form  y  when a is stored in lower triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j )       + temp1*a( j, j )
               do 90, i = j + 1, n
                  y( i ) = y( i ) + temp1*a( i, j )
                  temp2  = temp2  + a( i, j )*x( i )
   90          continue
               y( j ) = y( j ) + alpha*temp2
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy )       + temp1*a( j, j )
               ix      = jx
               iy      = jy
               do 110, i = j + 1, n
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*a( i, j )
                  temp2   = temp2   + a( i, j )*x( ix )
  110          continue
               y( jy ) = y( jy ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
  120       continue
         end if
      end if
c
      return
c
c     end of ssymv .
c
      end
c
c***********************************************************************
c
      subroutine ssyr  ( uplo, n, alpha, x, incx, a, lda )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha
      integer            incx, lda, n
      character*1        uplo
c     .. array arguments ..
      real*8             a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  ssyr   performs the symmetric rank 1 operation
c
c     a := alpha*x*x' + a,
c
c  where alpha is a real scalar, x is an n element vector and a is an
c  n by n symmetric matrix.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the array a is to be referenced as
c           follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of a
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of a
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular part of the symmetric matrix and the strictly
c           lower triangular part of a is not referenced. on exit, the
c           upper triangular part of the array a is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular part of the symmetric matrix and the strictly
c           upper triangular part of a is not referenced. on exit, the
c           lower triangular part of the array a is overwritten by the
c           lower triangular part of the updated matrix.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, j, jx, kx
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( lda.lt.max( 1, n ) )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'ssyr  ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     set the start point in x if the increment is not unity.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the triangular part
c     of a.
c
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when a is stored in upper triangle.
c
         if( incx.eq.1 )then
            do 20, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*x( j )
                  do 10, i = 1, j
                     a( i, j ) = a( i, j ) + x( i )*temp
   10             continue
               end if
   20       continue
         else
            jx = kx
            do 40, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  ix   = kx
                  do 30, i = 1, j
                     a( i, j ) = a( i, j ) + x( ix )*temp
                     ix        = ix        + incx
   30             continue
               end if
               jx = jx + incx
   40       continue
         end if
      else
c
c        form  a  when a is stored in lower triangle.
c
         if( incx.eq.1 )then
            do 60, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*x( j )
                  do 50, i = j, n
                     a( i, j ) = a( i, j ) + x( i )*temp
   50             continue
               end if
   60       continue
         else
            jx = kx
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  ix   = jx
                  do 70, i = j, n
                     a( i, j ) = a( i, j ) + x( ix )*temp
                     ix        = ix        + incx
   70             continue
               end if
               jx = jx + incx
   80       continue
         end if
      end if
c
      return
c
c     end of ssyr  .
c
      end
c
c***********************************************************************
c
      subroutine ssyr2 ( uplo, n, alpha, x, incx, y, incy, a, lda )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      real*8             alpha
      integer            incx, incy, lda, n
      character*1        uplo
c     .. array arguments ..
      real*8             a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  ssyr2  performs the symmetric rank 2 operation
c
c     a := alpha*x*y' + alpha*y*x' + a,
c
c  where alpha is a scalar, x and y are n element vectors and a is an n
c  by n symmetric matrix.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the array a is to be referenced as
c           follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of a
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of a
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y.
c           unchanged on exit.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular part of the symmetric matrix and the strictly
c           lower triangular part of a is not referenced. on exit, the
c           upper triangular part of the array a is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular part of the symmetric matrix and the strictly
c           upper triangular part of a is not referenced. on exit, the
c           lower triangular part of the array a is overwritten by the
c           lower triangular part of the updated matrix.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      else if( lda.lt.max( 1, n ) )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'ssyr2 ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     set up the start points in x and y if the increments are not both
c     unity.
c
      if( ( incx.ne.1 ).or.( incy.ne.1 ) )then
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            ky = 1
         else
            ky = 1 - ( n - 1 )*incy
         end if
         jx = kx
         jy = ky
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the triangular part
c     of a.
c
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when a is stored in the upper triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 20, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1 = alpha*y( j )
                  temp2 = alpha*x( j )
                  do 10, i = 1, j
                     a( i, j ) = a( i, j ) + x( i )*temp1 + y( i )*temp2
   10             continue
               end if
   20       continue
         else
            do 40, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1 = alpha*y( jy )
                  temp2 = alpha*x( jx )
                  ix    = kx
                  iy    = ky
                  do 30, i = 1, j
                     a( i, j ) = a( i, j ) + x( ix )*temp1
     $                                     + y( iy )*temp2
                     ix        = ix        + incx
                     iy        = iy        + incy
   30             continue
               end if
               jx = jx + incx
               jy = jy + incy
   40       continue
         end if
      else
c
c        form  a  when a is stored in the lower triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1 = alpha*y( j )
                  temp2 = alpha*x( j )
                  do 50, i = j, n
                     a( i, j ) = a( i, j ) + x( i )*temp1 + y( i )*temp2
   50             continue
               end if
   60       continue
         else
            do 80, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1 = alpha*y( jy )
                  temp2 = alpha*x( jx )
                  ix    = jx
                  iy    = jy
                  do 70, i = j, n
                     a( i, j ) = a( i, j ) + x( ix )*temp1
     $                                     + y( iy )*temp2
                     ix        = ix        + incx
                     iy        = iy        + incy
   70             continue
               end if
               jx = jx + incx
               jy = jy + incy
   80       continue
         end if
      end if
c
      return
c
c     end of ssyr2 .
c
      end
c
c***********************************************************************
c
      subroutine stbmv ( uplo, trans, diag, n, k, a, lda, x, incx )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      integer            incx, k, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      real*8             a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  stbmv  performs one of the matrix-vector operations
c
c     x := a*x,   or   x := a'*x,
c
c  where x is an n element vector and  a is an n by n unit, or non-unit,
c  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   x := a*x.
c
c              trans = 't' or 't'   x := a'*x.
c
c              trans = 'c' or 'c'   x := a'*x.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with uplo = 'u' or 'u', k specifies the number of
c           super-diagonals of the matrix a.
c           on entry with uplo = 'l' or 'l', k specifies the number of
c           sub-diagonals of the matrix a.
c           k must satisfy  0 .le. k.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
c           by n part of the array a must contain the upper triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row
c           ( k + 1 ) of the array, the first super-diagonal starting at
c           position 2 in row k, and so on. the top left k by k triangle
c           of the array a is not referenced.
c           the following program segment will transfer an upper
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 do 20, j = 1, n
c                    m = k + 1 - j
c                    do 10, i = max( 1, j - k ), j
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           before entry with uplo = 'l' or 'l', the leading ( k + 1 )
c           by n part of the array a must contain the lower triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row 1 of
c           the array, the first sub-diagonal starting at position 1 in
c           row 2, and so on. the bottom right k by k triangle of the
c           array a is not referenced.
c           the following program segment will transfer a lower
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 do 20, j = 1, n
c                    m = 1 - j
c                    do 10, i = j, min( n, j + k )
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           note that when diag = 'u' or 'u' the elements of the array a
c           corresponding to the diagonal elements of the matrix are not
c           referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( k + 1 ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x. on exit, x is overwritten with the
c           tranformed vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, j, jx, kplus1, kx, l
      logical            nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max, min
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( k.lt.0 )then
         info = 5
      else if( lda.lt.( k + 1 ) )then
         info = 7
      else if( incx.eq.0 )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'stbmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      nounit = lsame( diag, 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx   too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c         form  x := a*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     l    = kplus1 - j
                     do 10, i = max( 1, j - k ), j - 1
                        x( i ) = x( i ) + temp*a( l + i, j )
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( kplus1, j )
                  end if
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     l    = kplus1  - j
                     do 30, i = max( 1, j - k ), j - 1
                        x( ix ) = x( ix ) + temp*a( l + i, j )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( kplus1, j )
                  end if
                  jx = jx + incx
                  if( j.gt.k )
     $               kx = kx + incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     l    = 1      - j
                     do 50, i = min( n, j + k ), j + 1, -1
                        x( i ) = x( i ) + temp*a( l + i, j )
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( 1, j )
                  end if
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     l    = 1       - j
                     do 70, i = min( n, j + k ), j + 1, -1
                        x( ix ) = x( ix ) + temp*a( l + i, j )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( 1, j )
                  end if
                  jx = jx - incx
                  if( ( n - j ).ge.k )
     $               kx = kx - incx
   80          continue
            end if
         end if
      else
c
c        form  x := a'*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 100, j = n, 1, -1
                  temp = x( j )
                  l    = kplus1 - j
                  if( nounit )
     $               temp = temp*a( kplus1, j )
                  do 90, i = j - 1, max( 1, j - k ), -1
                     temp = temp + a( l + i, j )*x( i )
   90             continue
                  x( j ) = temp
  100          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 120, j = n, 1, -1
                  temp = x( jx )
                  kx   = kx      - incx
                  ix   = kx
                  l    = kplus1  - j
                  if( nounit )
     $               temp = temp*a( kplus1, j )
                  do 110, i = j - 1, max( 1, j - k ), -1
                     temp = temp + a( l + i, j )*x( ix )
                     ix   = ix   - incx
  110             continue
                  x( jx ) = temp
                  jx      = jx   - incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = 1, n
                  temp = x( j )
                  l    = 1      - j
                  if( nounit )
     $               temp = temp*a( 1, j )
                  do 130, i = j + 1, min( n, j + k )
                     temp = temp + a( l + i, j )*x( i )
  130             continue
                  x( j ) = temp
  140          continue
            else
               jx = kx
               do 160, j = 1, n
                  temp = x( jx )
                  kx   = kx      + incx
                  ix   = kx
                  l    = 1       - j
                  if( nounit )
     $               temp = temp*a( 1, j )
                  do 150, i = j + 1, min( n, j + k )
                     temp = temp + a( l + i, j )*x( ix )
                     ix   = ix   + incx
  150             continue
                  x( jx ) = temp
                  jx      = jx   + incx
  160          continue
            end if
         end if
      end if
c
      return
c
c     end of stbmv .
c
      end
c
c***********************************************************************
c
      subroutine stbsv ( uplo, trans, diag, n, k, a, lda, x, incx )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      integer            incx, k, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      real*8             a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  stbsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular band matrix, with ( k + 1 )
c  diagonals.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   a'*x = b.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with uplo = 'u' or 'u', k specifies the number of
c           super-diagonals of the matrix a.
c           on entry with uplo = 'l' or 'l', k specifies the number of
c           sub-diagonals of the matrix a.
c           k must satisfy  0 .le. k.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
c           by n part of the array a must contain the upper triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row
c           ( k + 1 ) of the array, the first super-diagonal starting at
c           position 2 in row k, and so on. the top left k by k triangle
c           of the array a is not referenced.
c           the following program segment will transfer an upper
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 do 20, j = 1, n
c                    m = k + 1 - j
c                    do 10, i = max( 1, j - k ), j
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           before entry with uplo = 'l' or 'l', the leading ( k + 1 )
c           by n part of the array a must contain the lower triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row 1 of
c           the array, the first sub-diagonal starting at position 1 in
c           row 2, and so on. the bottom right k by k triangle of the
c           array a is not referenced.
c           the following program segment will transfer a lower
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 do 20, j = 1, n
c                    m = 1 - j
c                    do 10, i = j, min( n, j + k )
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           note that when diag = 'u' or 'u' the elements of the array a
c           corresponding to the diagonal elements of the matrix are not
c           referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( k + 1 ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, j, jx, kplus1, kx, l
      logical            nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max, min
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( k.lt.0 )then
         info = 5
      else if( lda.lt.( k + 1 ) )then
         info = 7
      else if( incx.eq.0 )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'stbsv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      nounit = lsame( diag, 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed by sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := inv( a )*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     l = kplus1 - j
                     if( nounit )
     $                  x( j ) = x( j )/a( kplus1, j )
                     temp = x( j )
                     do 10, i = j - 1, max( 1, j - k ), -1
                        x( i ) = x( i ) - temp*a( l + i, j )
   10                continue
                  end if
   20          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 40, j = n, 1, -1
                  kx = kx - incx
                  if( x( jx ).ne.zero )then
                     ix = kx
                     l  = kplus1 - j
                     if( nounit )
     $                  x( jx ) = x( jx )/a( kplus1, j )
                     temp = x( jx )
                     do 30, i = j - 1, max( 1, j - k ), -1
                        x( ix ) = x( ix ) - temp*a( l + i, j )
                        ix      = ix      - incx
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     l = 1 - j
                     if( nounit )
     $                  x( j ) = x( j )/a( 1, j )
                     temp = x( j )
                     do 50, i = j + 1, min( n, j + k )
                        x( i ) = x( i ) - temp*a( l + i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  kx = kx + incx
                  if( x( jx ).ne.zero )then
                     ix = kx
                     l  = 1  - j
                     if( nounit )
     $                  x( jx ) = x( jx )/a( 1, j )
                     temp = x( jx )
                     do 70, i = j + 1, min( n, j + k )
                        x( ix ) = x( ix ) - temp*a( l + i, j )
                        ix      = ix      + incx
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else
c
c        form  x := inv( a')*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 100, j = 1, n
                  temp = x( j )
                  l    = kplus1 - j
                  do 90, i = max( 1, j - k ), j - 1
                     temp = temp - a( l + i, j )*x( i )
   90             continue
                  if( nounit )
     $               temp = temp/a( kplus1, j )
                  x( j ) = temp
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  l    = kplus1  - j
                  do 110, i = max( 1, j - k ), j - 1
                     temp = temp - a( l + i, j )*x( ix )
                     ix   = ix   + incx
  110             continue
                  if( nounit )
     $               temp = temp/a( kplus1, j )
                  x( jx ) = temp
                  jx      = jx   + incx
                  if( j.gt.k )
     $               kx = kx + incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = n, 1, -1
                  temp = x( j )
                  l    = 1      - j
                  do 130, i = min( n, j + k ), j + 1, -1
                     temp = temp - a( l + i, j )*x( i )
  130             continue
                  if( nounit )
     $               temp = temp/a( 1, j )
                  x( j ) = temp
  140          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 160, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  l    = 1       - j
                  do 150, i = min( n, j + k ), j + 1, -1
                     temp = temp - a( l + i, j )*x( ix )
                     ix   = ix   - incx
  150             continue
                  if( nounit )
     $               temp = temp/a( 1, j )
                  x( jx ) = temp
                  jx      = jx   - incx
                  if( ( n - j ).ge.k )
     $               kx = kx - incx
  160          continue
            end if
         end if
      end if
c
      return
c
c     end of stbsv .
c
      end
c
c***********************************************************************
c
      subroutine stpmv ( uplo, trans, diag, n, ap, x, incx )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      integer            incx, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      real*8             ap( * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  stpmv  performs one of the matrix-vector operations
c
c     x := a*x,   or   x := a'*x,
c
c  where x is an n element vector and  a is an n by n unit, or non-unit,
c  upper or lower triangular matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   x := a*x.
c
c              trans = 't' or 't'   x := a'*x.
c
c              trans = 'c' or 'c'   x := a'*x.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  ap     - real             array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )
c           respectively, and so on.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )
c           respectively, and so on.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x. on exit, x is overwritten with the
c           tranformed vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, j, jx, k, kk, kx
      logical            nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( incx.eq.0 )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'stpmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      nounit = lsame( diag, 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of ap are
c     accessed sequentially with one pass through ap.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x:= a*x.
c
         if( lsame( uplo, 'u' ) )then
            kk =1
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     k    = kk
                     do 10, i = 1, j - 1
                        x( i ) = x( i ) + temp*ap( k )
                        k      = k      + 1
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*ap( kk + j - 1 )
                  end if
                  kk = kk + j
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 30, k = kk, kk + j - 2
                        x( ix ) = x( ix ) + temp*ap( k )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*ap( kk + j - 1 )
                  end if
                  jx = jx + incx
                  kk = kk + j
   40          continue
            end if
         else
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     k    = kk
                     do 50, i = n, j + 1, -1
                        x( i ) = x( i ) + temp*ap( k )
                        k      = k      - 1
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*ap( kk - n + j )
                  end if
                  kk = kk - ( n - j + 1 )
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 70, k = kk, kk - ( n - ( j + 1 ) ), -1
                        x( ix ) = x( ix ) + temp*ap( k )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*ap( kk - n + j )
                  end if
                  jx = jx - incx
                  kk = kk - ( n - j + 1 )
   80          continue
            end if
         end if
      else
c
c        form  x := a'*x.
c
         if( lsame( uplo, 'u' ) )then
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 100, j = n, 1, -1
                  temp = x( j )
                  if( nounit )
     $               temp = temp*ap( kk )
                  k = kk - 1
                  do 90, i = j - 1, 1, -1
                     temp = temp + ap( k )*x( i )
                     k    = k    - 1
   90             continue
                  x( j ) = temp
                  kk     = kk   - j
  100          continue
            else
               jx = kx + ( n - 1 )*incx
               do 120, j = n, 1, -1
                  temp = x( jx )
                  ix   = jx
                  if( nounit )
     $               temp = temp*ap( kk )
                  do 110, k = kk - 1, kk - j + 1, -1
                     ix   = ix   - incx
                     temp = temp + ap( k )*x( ix )
  110             continue
                  x( jx ) = temp
                  jx      = jx   - incx
                  kk      = kk   - j
  120          continue
            end if
         else
            kk = 1
            if( incx.eq.1 )then
               do 140, j = 1, n
                  temp = x( j )
                  if( nounit )
     $               temp = temp*ap( kk )
                  k = kk + 1
                  do 130, i = j + 1, n
                     temp = temp + ap( k )*x( i )
                     k    = k    + 1
  130             continue
                  x( j ) = temp
                  kk     = kk   + ( n - j + 1 )
  140          continue
            else
               jx = kx
               do 160, j = 1, n
                  temp = x( jx )
                  ix   = jx
                  if( nounit )
     $               temp = temp*ap( kk )
                  do 150, k = kk + 1, kk + n - j
                     ix   = ix   + incx
                     temp = temp + ap( k )*x( ix )
  150             continue
                  x( jx ) = temp
                  jx      = jx   + incx
                  kk      = kk   + ( n - j + 1 )
  160          continue
            end if
         end if
      end if
c
      return
c
c     end of stpmv .
c
      end
c
c***********************************************************************
c
      subroutine stpsv ( uplo, trans, diag, n, ap, x, incx )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      integer            incx, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      real*8             ap( * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  stpsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix, supplied in packed form.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   a'*x = b.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  ap     - real             array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )
c           respectively, and so on.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )
c           respectively, and so on.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, j, jx, k, kk, kx
      logical            nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( incx.eq.0 )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'stpsv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      nounit = lsame( diag, 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of ap are
c     accessed sequentially with one pass through ap.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := inv( a )*x.
c
         if( lsame( uplo, 'u' ) )then
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/ap( kk )
                     temp = x( j )
                     k    = kk     - 1
                     do 10, i = j - 1, 1, -1
                        x( i ) = x( i ) - temp*ap( k )
                        k      = k      - 1
   10                continue
                  end if
                  kk = kk - j
   20          continue
            else
               jx = kx + ( n - 1 )*incx
               do 40, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/ap( kk )
                     temp = x( jx )
                     ix   = jx
                     do 30, k = kk - 1, kk - j + 1, -1
                        ix      = ix      - incx
                        x( ix ) = x( ix ) - temp*ap( k )
   30                continue
                  end if
                  jx = jx - incx
                  kk = kk - j
   40          continue
            end if
         else
            kk = 1
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/ap( kk )
                     temp = x( j )
                     k    = kk     + 1
                     do 50, i = j + 1, n
                        x( i ) = x( i ) - temp*ap( k )
                        k      = k      + 1
   50                continue
                  end if
                  kk = kk + ( n - j + 1 )
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/ap( kk )
                     temp = x( jx )
                     ix   = jx
                     do 70, k = kk + 1, kk + n - j
                        ix      = ix      + incx
                        x( ix ) = x( ix ) - temp*ap( k )
   70                continue
                  end if
                  jx = jx + incx
                  kk = kk + ( n - j + 1 )
   80          continue
            end if
         end if
      else
c
c        form  x := inv( a' )*x.
c
         if( lsame( uplo, 'u' ) )then
            kk = 1
            if( incx.eq.1 )then
               do 100, j = 1, n
                  temp = x( j )
                  k    = kk
                  do 90, i = 1, j - 1
                     temp = temp - ap( k )*x( i )
                     k    = k    + 1
   90             continue
                  if( nounit )
     $               temp = temp/ap( kk + j - 1 )
                  x( j ) = temp
                  kk     = kk   + j
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  do 110, k = kk, kk + j - 2
                     temp = temp - ap( k )*x( ix )
                     ix   = ix   + incx
  110             continue
                  if( nounit )
     $               temp = temp/ap( kk + j - 1 )
                  x( jx ) = temp
                  jx      = jx   + incx
                  kk      = kk   + j
  120          continue
            end if
         else
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 140, j = n, 1, -1
                  temp = x( j )
                  k = kk
                  do 130, i = n, j + 1, -1
                     temp = temp - ap( k )*x( i )
                     k    = k    - 1
  130             continue
                  if( nounit )
     $               temp = temp/ap( kk - n + j )
                  x( j ) = temp
                  kk     = kk   - ( n - j + 1 )
  140          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 160, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  do 150, k = kk, kk - ( n - ( j + 1 ) ), -1
                     temp = temp - ap( k )*x( ix )
                     ix   = ix   - incx
  150             continue
                  if( nounit )
     $               temp = temp/ap( kk - n + j )
                  x( jx ) = temp
                  jx      = jx   - incx
                  kk      = kk   - (n - j + 1 )
  160          continue
            end if
         end if
      end if
c
      return
c
c     end of stpsv .
c
      end
c
c***********************************************************************
c
      subroutine strmv ( uplo, trans, diag, n, a, lda, x, incx )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      real*8             a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  strmv  performs one of the matrix-vector operations
c
c     x := a*x,   or   x := a'*x,
c
c  where x is an n element vector and  a is an n by n unit, or non-unit,
c  upper or lower triangular matrix.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   x := a*x.
c
c              trans = 't' or 't'   x := a'*x.
c
c              trans = 'c' or 'c'   x := a'*x.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x. on exit, x is overwritten with the
c           tranformed vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, j, jx, kx
      logical            nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call xerbla( 'strmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      nounit = lsame( diag, 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := a*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 10, i = 1, j - 1
                        x( i ) = x( i ) + temp*a( i, j )
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 30, i = 1, j - 1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx + incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 50, i = n, j + 1, -1
                        x( i ) = x( i ) + temp*a( i, j )
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 70, i = n, j + 1, -1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx - incx
   80          continue
            end if
         end if
      else
c
c        form  x := a'*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 100, j = n, 1, -1
                  temp = x( j )
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 90, i = j - 1, 1, -1
                     temp = temp + a( i, j )*x( i )
   90             continue
                  x( j ) = temp
  100          continue
            else
               jx = kx + ( n - 1 )*incx
               do 120, j = n, 1, -1
                  temp = x( jx )
                  ix   = jx
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 110, i = j - 1, 1, -1
                     ix   = ix   - incx
                     temp = temp + a( i, j )*x( ix )
  110             continue
                  x( jx ) = temp
                  jx      = jx   - incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = 1, n
                  temp = x( j )
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 130, i = j + 1, n
                     temp = temp + a( i, j )*x( i )
  130             continue
                  x( j ) = temp
  140          continue
            else
               jx = kx
               do 160, j = 1, n
                  temp = x( jx )
                  ix   = jx
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 150, i = j + 1, n
                     ix   = ix   + incx
                     temp = temp + a( i, j )*x( ix )
  150             continue
                  x( jx ) = temp
                  jx      = jx   + incx
  160          continue
            end if
         end if
      end if
c
      return
c
c     end of strmv .
c
      end
c
c***********************************************************************
c
      subroutine strsv ( uplo, trans, diag, n, a, lda, x, incx )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      real*8             a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  strsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   a'*x = b.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - real             array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      real*8             zero
      parameter        ( zero = 0.0d+0 )
c     .. local scalars ..
      real*8             temp
      integer            i, info, ix, j, jx, kx
      logical            nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call xerbla( 'strsv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      nounit = lsame( diag, 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := inv( a )*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 10, i = j - 1, 1, -1
                        x( i ) = x( i ) - temp*a( i, j )
   10                continue
                  end if
   20          continue
            else
               jx = kx + ( n - 1 )*incx
               do 40, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 30, i = j - 1, 1, -1
                        ix      = ix      - incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 50, i = j + 1, n
                        x( i ) = x( i ) - temp*a( i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 70, i = j + 1, n
                        ix      = ix      + incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else
c
c        form  x := inv( a' )*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 100, j = 1, n
                  temp = x( j )
                  do 90, i = 1, j - 1
                     temp = temp - a( i, j )*x( i )
   90             continue
                  if( nounit )
     $               temp = temp/a( j, j )
                  x( j ) = temp
  100          continue
            else
               jx = kx
               do 120, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  do 110, i = 1, j - 1
                     temp = temp - a( i, j )*x( ix )
                     ix   = ix   + incx
  110             continue
                  if( nounit )
     $               temp = temp/a( j, j )
                  x( jx ) = temp
                  jx      = jx   + incx
  120          continue
            end if
         else
            if( incx.eq.1 )then
               do 140, j = n, 1, -1
                  temp = x( j )
                  do 130, i = n, j + 1, -1
                     temp = temp - a( i, j )*x( i )
  130             continue
                  if( nounit )
     $               temp = temp/a( j, j )
                  x( j ) = temp
  140          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 160, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  do 150, i = n, j + 1, -1
                     temp = temp - a( i, j )*x( ix )
                     ix   = ix   - incx
  150             continue
                  if( nounit )
     $               temp = temp/a( j, j )
                  x( jx ) = temp
                  jx      = jx   - incx
  160          continue
            end if
         end if
      end if
c
      return
c
c     end of strsv .
c
      end
      subroutine xerbla ( srname, info )
c     ..    scalar arguments ..
      integer            info
      character*6        srname
c     ..
c
c  purpose
c  =======
c
c  xerbla  is an error handler for the level 2 blas routines.
c
c  it is called by the level 2 blas routines if an input parameter is
c  invalid.
c
c  installers should consider modifying the stop statement in order to
c  call system-specific exception-handling facilities.
c
c  parameters
c  ==========
c
c  srname - character*6.
c           on entry, srname specifies the name of the routine which
c           called xerbla.
c
c  info   - integer.
c           on entry, info specifies the position of the invalid
c           parameter in the parameter-list of the calling routine.
c
c
c  auxiliary routine for level 2 blas.
c
c  written on 20-july-1986.
c
c     .. executable statements ..
c
      write (*,99999) srname, info
c
      stop
c
99999 format ( ' ** on entry to ', a6, ' parameter number ', i2,
     $         ' had an illegal value' )
c
c     end of xerbla.
c
      end
      logical function lsame ( ca, cb )
c     .. scalar arguments ..
      character*1            ca, cb
c     ..
c
c  purpose
c  =======
c
c  lsame  tests if ca is the same letter as cb regardless of case.
c  cb is assumed to be an upper case letter. lsame returns .true. if
c  ca is either the same as cb or the equivalent lower case letter.
c
c  n.b. this version of the routine is only correct for ascii code.
c       installers must modify the routine for other character-codes.
c
c       for ebcdic systems the constant ioff must be changed to -64.
c       for cdc systems using 6-12 bit representations, the system-
c       specific code in comments must be activated.
c
c  parameters
c  ==========
c
c  ca     - character*1
c  cb     - character*1
c           on entry, ca and cb specify characters to be compared.
c           unchanged on exit.
c
c
c  auxiliary routine for level 2 blas.
c
c  -- written on 20-july-1986
c     richard hanson, sandia national labs.
c     jeremy du croz, nag central office.
c
c     .. parameters ..
      integer                ioff
      parameter            ( ioff=32 )
c     .. intrinsic functions ..
      intrinsic              ichar
c     .. executable statements ..
c
c     test if the characters are equal
c
      lsame = ca .eq. cb
c
c     now test for equivalence
c
      if ( .not.lsame ) then
         lsame = ichar(ca) - ioff .eq. ichar(cb)
      end if
c
      return
c
c     end of lsame.
c
      end
c***********************************************************************
c
      subroutine cgbmv ( trans, m, n, kl, ku, alpha, a, lda, x, incx,
     $                   beta, y, incy )
c     .. scalar arguments ..
      complex            alpha, beta
      integer            incx, incy, kl, ku, lda, m, n
      character*1        trans
c     .. array arguments ..
      complex            a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  cgbmv  performs one of the matrix-vector operations
c
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,   or
c
c     y := alpha*conjg( a' )*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
c
c  parameters
c  ==========
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
c
c              trans = 't' or 't'   y := alpha*a'*x + beta*y.
c
c              trans = 'c' or 'c'   y := alpha*conjg( a' )*x + beta*y.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  kl     - integer.
c           on entry, kl specifies the number of sub-diagonals of the
c           matrix a. kl must satisfy  0 .le. kl.
c           unchanged on exit.
c
c  ku     - integer.
c           on entry, ku specifies the number of super-diagonals of the
c           matrix a. ku must satisfy  0 .le. ku.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry, the leading ( kl + ku + 1 ) by n part of the
c           array a must contain the matrix of coefficients, supplied
c           column by column, with the leading diagonal of the matrix in
c           row ( ku + 1 ) of the array, the first super-diagonal
c           starting at position 2 in row ku, the first sub-diagonal
c           starting at position 1 in row ( ku + 2 ), and so on.
c           elements in the array a that do not correspond to elements
c           in the band matrix (such as the top left ku by ku triangle)
c           are not referenced.
c           the following program segment will transfer a band matrix
c           from conventional full matrix storage to band storage:
c
c                 do 20, j = 1, n
c                    k = ku + 1 - j
c                    do 10, i = max( 1, j - ku ), min( m, j + kl )
c                       a( k + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( kl + ku + 1 ).
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
c           before entry, the incremented array y must contain the
c           vector y. on exit, y is overwritten by the updated vector y.
c
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, iy, j, jx, jy, k, kup1, kx, ky,
     $                   lenx, leny
      logical            noconj
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, min
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( kl.lt.0 )then
         info = 4
      else if( ku.lt.0 )then
         info = 5
      else if( lda.lt.( kl + ku + 1 ) )then
         info = 8
      else if( incx.eq.0 )then
         info = 10
      else if( incy.eq.0 )then
         info = 13
      end if
      if( info.ne.0 )then
         call xerbla( 'cgbmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
      noconj = lsame( trans, 't' )
c
c     set  lenx  and  leny, the lengths of the vectors x and y, and set
c     up the start points in  x  and  y.
c
      if( lsame( trans, 'n' ) )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the band part of a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      kup1 = ku + 1
      if( lsame( trans, 'n' ) )then
c
c        form  y := alpha*a*x + y.
c
         jx = kx
         if( incy.eq.1 )then
            do 60, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  k    = kup1 - j
                  do 50, i = max( 1, j - ku ), min( m, j + kl )
                     y( i ) = y( i ) + temp*a( k + i, j )
   50             continue
               end if
               jx = jx + incx
   60       continue
         else
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy   = ky
                  k    = kup1 - j
                  do 70, i = max( 1, j - ku ), min( m, j + kl )
                     y( iy ) = y( iy ) + temp*a( k + i, j )
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
               if( j.gt.ku )
     $            ky = ky + incy
   80       continue
         end if
      else
c
c        form  y := alpha*a'*x + y  or  y := alpha*conjg( a' )*x + y.
c
         jy = ky
         if( incx.eq.1 )then
            do 110, j = 1, n
               temp = zero
               k    = kup1 - j
               if( noconj )then
                  do 90, i = max( 1, j - ku ), min( m, j + kl )
                     temp = temp + a( k + i, j )*x( i )
   90             continue
               else
                  do 100, i = max( 1, j - ku ), min( m, j + kl )
                     temp = temp + conjg( a( k + i, j ) )*x( i )
  100             continue
               end if
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  110       continue
         else
            do 140, j = 1, n
               temp = zero
               ix   = kx
               k    = kup1 - j
               if( noconj )then
                  do 120, i = max( 1, j - ku ), min( m, j + kl )
                     temp = temp + a( k + i, j )*x( ix )
                     ix   = ix   + incx
  120             continue
               else
                  do 130, i = max( 1, j - ku ), min( m, j + kl )
                     temp = temp + conjg( a( k + i, j ) )*x( ix )
                     ix   = ix   + incx
  130             continue
               end if
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
               if( j.gt.ku )
     $            kx = kx + incx
  140       continue
         end if
      end if
c
      return
c
c     end of cgbmv .
c
      end
c
c***********************************************************************
c
      subroutine cgemv ( trans, m, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )
c     .. scalar arguments ..
      complex            alpha, beta
      integer            incx, incy, lda, m, n
      character*1        trans
c     .. array arguments ..
      complex            a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  cgemv  performs one of the matrix-vector operations
c
c     y := alpha*a*x + beta*y,   or   y := alpha*a'*x + beta*y,   or
c
c     y := alpha*conjg( a' )*x + beta*y,
c
c  where alpha and beta are scalars, x and y are vectors and a is an
c  m by n matrix.
c
c  parameters
c  ==========
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   y := alpha*a*x + beta*y.
c
c              trans = 't' or 't'   y := alpha*a'*x + beta*y.
c
c              trans = 'c' or 'c'   y := alpha*conjg( a' )*x + beta*y.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( m - 1 )*abs( incx ) ) otherwise.
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( m - 1 )*abs( incy ) ) when trans = 'n' or 'n'
c           and at least
c           ( 1 + ( n - 1 )*abs( incy ) ) otherwise.
c           before entry with beta non-zero, the incremented array y
c           must contain the vector y. on exit, y is overwritten by the
c           updated vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, iy, j, jx, jy, kx, ky, lenx, leny
      logical            noconj
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 1
      else if( m.lt.0 )then
         info = 2
      else if( n.lt.0 )then
         info = 3
      else if( lda.lt.max( 1, m ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'cgemv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
      noconj = lsame( trans, 't' )
c
c     set  lenx  and  leny, the lengths of the vectors x and y, and set
c     up the start points in  x  and  y.
c
      if( lsame( trans, 'n' ) )then
         lenx = n
         leny = m
      else
         lenx = m
         leny = n
      end if
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( lenx - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( leny - 1 )*incy
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, leny
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, leny
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, leny
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, leny
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( lsame( trans, 'n' ) )then
c
c        form  y := alpha*a*x + y.
c
         jx = kx
         if( incy.eq.1 )then
            do 60, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  do 50, i = 1, m
                     y( i ) = y( i ) + temp*a( i, j )
   50             continue
               end if
               jx = jx + incx
   60       continue
         else
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*x( jx )
                  iy   = ky
                  do 70, i = 1, m
                     y( iy ) = y( iy ) + temp*a( i, j )
                     iy      = iy      + incy
   70             continue
               end if
               jx = jx + incx
   80       continue
         end if
      else
c
c        form  y := alpha*a'*x + y  or  y := alpha*conjg( a' )*x + y.
c
         jy = ky
         if( incx.eq.1 )then
            do 110, j = 1, n
               temp = zero
               if( noconj )then
                  do 90, i = 1, m
                     temp = temp + a( i, j )*x( i )
   90             continue
               else
                  do 100, i = 1, m
                     temp = temp + conjg( a( i, j ) )*x( i )
  100             continue
               end if
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  110       continue
         else
            do 140, j = 1, n
               temp = zero
               ix   = kx
               if( noconj )then
                  do 120, i = 1, m
                     temp = temp + a( i, j )*x( ix )
                     ix   = ix   + incx
  120             continue
               else
                  do 130, i = 1, m
                     temp = temp + conjg( a( i, j ) )*x( ix )
                     ix   = ix   + incx
  130             continue
               end if
               y( jy ) = y( jy ) + alpha*temp
               jy      = jy      + incy
  140       continue
         end if
      end if
c
      return
c
c     end of cgemv .
c
      end
c
c***********************************************************************
c
      subroutine cgerc ( m, n, alpha, x, incx, y, incy, a, lda )
c     .. scalar arguments ..
      complex            alpha
      integer            incx, incy, lda, m, n
c     .. array arguments ..
      complex            a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  cgerc  performs the rank 1 operation
c
c     a := alpha*x*conjg( y' ) + a,
c
c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and a is an m by n matrix.
c
c  parameters
c  ==========
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( m - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the m
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y.
c           unchanged on exit.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients. on exit, a is
c           overwritten by the updated matrix.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jy, kx
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( m.lt.0 )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      else if( lda.lt.max( 1, m ) )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'cgerc ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( incy.gt.0 )then
         jy = 1
      else
         jy = 1 - ( n - 1 )*incy
      end if
      if( incx.eq.1 )then
         do 20, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*conjg( y( jy ) )
               do 10, i = 1, m
                  a( i, j ) = a( i, j ) + x( i )*temp
   10          continue
            end if
            jy = jy + incy
   20    continue
      else
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( m - 1 )*incx
         end if
         do 40, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*conjg( y( jy ) )
               ix   = kx
               do 30, i = 1, m
                  a( i, j ) = a( i, j ) + x( ix )*temp
                  ix        = ix        + incx
   30          continue
            end if
            jy = jy + incy
   40    continue
      end if
c
      return
c
c     end of cgerc .
c
      end
c
c***********************************************************************
c
      subroutine cgeru ( m, n, alpha, x, incx, y, incy, a, lda )
c     .. scalar arguments ..
      complex            alpha
      integer            incx, incy, lda, m, n
c     .. array arguments ..
      complex            a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  cgeru  performs the rank 1 operation
c
c     a := alpha*x*y' + a,
c
c  where alpha is a scalar, x is an m element vector, y is an n element
c  vector and a is an m by n matrix.
c
c  parameters
c  ==========
c
c  m      - integer.
c           on entry, m specifies the number of rows of the matrix a.
c           m must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( m - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the m
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y.
c           unchanged on exit.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry, the leading m by n part of the array a must
c           contain the matrix of coefficients. on exit, a is
c           overwritten by the updated matrix.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jy, kx
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( m.lt.0 )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      else if( lda.lt.max( 1, m ) )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'cgeru ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( incy.gt.0 )then
         jy = 1
      else
         jy = 1 - ( n - 1 )*incy
      end if
      if( incx.eq.1 )then
         do 20, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               do 10, i = 1, m
                  a( i, j ) = a( i, j ) + x( i )*temp
   10          continue
            end if
            jy = jy + incy
   20    continue
      else
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( m - 1 )*incx
         end if
         do 40, j = 1, n
            if( y( jy ).ne.zero )then
               temp = alpha*y( jy )
               ix   = kx
               do 30, i = 1, m
                  a( i, j ) = a( i, j ) + x( ix )*temp
                  ix        = ix        + incx
   30          continue
            end if
            jy = jy + incy
   40    continue
      end if
c
      return
c
c     end of cgeru .
c
      end
c
c***********************************************************************
c
      subroutine chbmv ( uplo, n, k, alpha, a, lda, x, incx,
     $                   beta, y, incy )
c     .. scalar arguments ..
      complex            alpha, beta
      integer            incx, incy, k, lda, n
      character*1        uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  chbmv  performs the matrix-vector  operation
c
c     y := alpha*a*x + beta*y,
c
c  where alpha and beta are scalars, x and y are n element vectors and
c  a is an n by n hermitian band matrix, with k super-diagonals.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the band matrix a is being supplied as
c           follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  being supplied.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  being supplied.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry, k specifies the number of super-diagonals of the
c           matrix a. k must satisfy  0 .le. k.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
c           by n part of the array a must contain the upper triangular
c           band part of the hermitian matrix, supplied column by
c           column, with the leading diagonal of the matrix in row
c           ( k + 1 ) of the array, the first super-diagonal starting at
c           position 2 in row k, and so on. the top left k by k triangle
c           of the array a is not referenced.
c           the following program segment will transfer the upper
c           triangular part of a hermitian band matrix from conventional
c           full matrix storage to band storage:
c
c                 do 20, j = 1, n
c                    m = k + 1 - j
c                    do 10, i = max( 1, j - k ), j
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           before entry with uplo = 'l' or 'l', the leading ( k + 1 )
c           by n part of the array a must contain the lower triangular
c           band part of the hermitian matrix, supplied column by
c           column, with the leading diagonal of the matrix in row 1 of
c           the array, the first sub-diagonal starting at position 1 in
c           row 2, and so on. the bottom right k by k triangle of the
c           array a is not referenced.
c           the following program segment will transfer the lower
c           triangular part of a hermitian band matrix from conventional
c           full matrix storage to band storage:
c
c                 do 20, j = 1, n
c                    m = 1 - j
c                    do 10, i = j, min( n, j + k )
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           note that the imaginary parts of the diagonal elements need
c           not be set and are assumed to be zero.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( k + 1 ).
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the
c           vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the
c           vector y. on exit, y is overwritten by the updated vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kplus1, kx, ky, l
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, min, real
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( k.lt.0 )then
         info = 3
      else if( lda.lt.( k + 1 ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      else if( incy.eq.0 )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'chbmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set up the start points in  x  and  y.
c
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( n - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( n - 1 )*incy
      end if
c
c     start the operations. in this version the elements of the array a
c     are accessed sequentially with one pass through a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, n
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, n
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, n
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, n
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( lsame( uplo, 'u' ) )then
c
c        form  y  when upper triangle of a is stored.
c
         kplus1 = k + 1
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               l     = kplus1 - j
               do 50, i = max( 1, j - k ), j - 1
                  y( i ) = y( i ) + temp1*a( l + i, j )
                  temp2  = temp2  + conjg( a( l + i, j ) )*x( i )
   50          continue
               y( j ) = y( j ) + temp1*real( a( kplus1, j ) )
     $                         + alpha*temp2
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               l     = kplus1 - j
               do 70, i = max( 1, j - k ), j - 1
                  y( iy ) = y( iy ) + temp1*a( l + i, j )
                  temp2   = temp2   + conjg( a( l + i, j ) )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*real( a( kplus1, j ) )
     $                           + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
               if( j.gt.k )then
                  kx = kx + incx
                  ky = ky + incy
               end if
   80       continue
         end if
      else
c
c        form  y  when lower triangle of a is stored.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j ) + temp1*real( a( 1, j ) )
               l      = 1      - j
               do 90, i = j + 1, min( n, j + k )
                  y( i ) = y( i ) + temp1*a( l + i, j )
                  temp2  = temp2  + conjg( a( l + i, j ) )*x( i )
   90          continue
               y( j ) = y( j ) + alpha*temp2
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy ) + temp1*real( a( 1, j ) )
               l       = 1       - j
               ix      = jx
               iy      = jy
               do 110, i = j + 1, min( n, j + k )
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*a( l + i, j )
                  temp2   = temp2   + conjg( a( l + i, j ) )*x( ix )
  110          continue
               y( jy ) = y( jy ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
  120       continue
         end if
      end if
c
      return
c
c     end of chbmv .
c
      end
c
c***********************************************************************
c
      subroutine chemv ( uplo, n, alpha, a, lda, x, incx,
     $                   beta, y, incy )
c     .. scalar arguments ..
      complex            alpha, beta
      integer            incx, incy, lda, n
      character*1        uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  chemv  performs the matrix-vector  operation
c
c     y := alpha*a*x + beta*y,
c
c  where alpha and beta are scalars, x and y are n element vectors and
c  a is an n by n hermitian matrix.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the array a is to be referenced as
c           follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of a
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of a
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular part of the hermitian matrix and the strictly
c           lower triangular part of a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular part of the hermitian matrix and the strictly
c           upper triangular part of a is not referenced.
c           note that the imaginary parts of the diagonal elements need
c           not be set and are assumed to be zero.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y. on exit, y is overwritten by the updated
c           vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, real
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( lda.lt.max( 1, n ) )then
         info = 5
      else if( incx.eq.0 )then
         info = 7
      else if( incy.eq.0 )then
         info = 10
      end if
      if( info.ne.0 )then
         call xerbla( 'chemv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set up the start points in  x  and  y.
c
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( n - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( n - 1 )*incy
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the triangular part
c     of a.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, n
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, n
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, n
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, n
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      if( lsame( uplo, 'u' ) )then
c
c        form  y  when a is stored in upper triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               do 50, i = 1, j - 1
                  y( i ) = y( i ) + temp1*a( i, j )
                  temp2  = temp2  + conjg( a( i, j ) )*x( i )
   50          continue
               y( j ) = y( j ) + temp1*real( a( j, j ) ) + alpha*temp2
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               do 70, i = 1, j - 1
                  y( iy ) = y( iy ) + temp1*a( i, j )
                  temp2   = temp2   + conjg( a( i, j ) )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*real( a( j, j ) ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
   80       continue
         end if
      else
c
c        form  y  when a is stored in lower triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j ) + temp1*real( a( j, j ) )
               do 90, i = j + 1, n
                  y( i ) = y( i ) + temp1*a( i, j )
                  temp2  = temp2  + conjg( a( i, j ) )*x( i )
   90          continue
               y( j ) = y( j ) + alpha*temp2
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy ) + temp1*real( a( j, j ) )
               ix      = jx
               iy      = jy
               do 110, i = j + 1, n
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*a( i, j )
                  temp2   = temp2   + conjg( a( i, j ) )*x( ix )
  110          continue
               y( jy ) = y( jy ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
  120       continue
         end if
      end if
c
      return
c
c     end of chemv .
c
      end
c
c***********************************************************************
c
      subroutine cher  ( uplo, n, alpha, x, incx, a, lda )
c     .. scalar arguments ..
      real               alpha
      integer            incx, lda, n
      character*1        uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  cher   performs the hermitian rank 1 operation
c
c     a := alpha*x*conjg( x' ) + a,
c
c  where alpha is a real scalar, x is an n element vector and a is an
c  n by n hermitian matrix.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the array a is to be referenced as
c           follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of a
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of a
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular part of the hermitian matrix and the strictly
c           lower triangular part of a is not referenced. on exit, the
c           upper triangular part of the array a is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular part of the hermitian matrix and the strictly
c           upper triangular part of a is not referenced. on exit, the
c           lower triangular part of the array a is overwritten by the
c           lower triangular part of the updated matrix.
c           note that the imaginary parts of the diagonal elements need
c           not be set, they are assumed to be zero, and on exit they
c           are set to zero.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, kx
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, real
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( lda.lt.max( 1, n ) )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'cher  ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.real( zero ) ) )
     $   return
c
c     set the start point in x if the increment is not unity.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the triangular part
c     of a.
c
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when a is stored in upper triangle.
c
         if( incx.eq.1 )then
            do 20, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*conjg( x( j ) )
                  do 10, i = 1, j - 1
                     a( i, j ) = a( i, j ) + x( i )*temp
   10             continue
                  a( j, j ) = real( a( j, j ) ) + real( x( j )*temp )
               else
                  a( j, j ) = real( a( j, j ) )
               end if
   20       continue
         else
            jx = kx
            do 40, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*conjg( x( jx ) )
                  ix   = kx
                  do 30, i = 1, j - 1
                     a( i, j ) = a( i, j ) + x( ix )*temp
                     ix        = ix        + incx
   30             continue
                  a( j, j ) = real( a( j, j ) ) + real( x( jx )*temp )
               else
                  a( j, j ) = real( a( j, j ) )
               end if
               jx = jx + incx
   40       continue
         end if
      else
c
c        form  a  when a is stored in lower triangle.
c
         if( incx.eq.1 )then
            do 60, j = 1, n
               if( x( j ).ne.zero )then
                  temp      = alpha*conjg( x( j ) )
                  a( j, j ) = real( a( j, j ) ) + real( temp*x( j ) )
                  do 50, i = j + 1, n
                     a( i, j ) = a( i, j ) + x( i )*temp
   50             continue
               else
                  a( j, j ) = real( a( j, j ) )
               end if
   60       continue
         else
            jx = kx
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp      = alpha*conjg( x( jx ) )
                  a( j, j ) = real( a( j, j ) ) + real( temp*x( jx ) )
                  ix        = jx
                  do 70, i = j + 1, n
                     ix        = ix        + incx
                     a( i, j ) = a( i, j ) + x( ix )*temp
   70             continue
               else
                  a( j, j ) = real( a( j, j ) )
               end if
               jx = jx + incx
   80       continue
         end if
      end if
c
      return
c
c     end of cher  .
c
      end
c
c***********************************************************************
c
      subroutine cher2 ( uplo, n, alpha, x, incx, y, incy, a, lda )
c     .. scalar arguments ..
      complex            alpha
      integer            incx, incy, lda, n
      character*1        uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  cher2  performs the hermitian rank 2 operation
c
c     a := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + a,
c
c  where alpha is a scalar, x and y are n element vectors and a is an n
c  by n hermitian matrix.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the array a is to be referenced as
c           follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of a
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of a
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y.
c           unchanged on exit.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular part of the hermitian matrix and the strictly
c           lower triangular part of a is not referenced. on exit, the
c           upper triangular part of the array a is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular part of the hermitian matrix and the strictly
c           upper triangular part of a is not referenced. on exit, the
c           lower triangular part of the array a is overwritten by the
c           lower triangular part of the updated matrix.
c           note that the imaginary parts of the diagonal elements need
c           not be set, they are assumed to be zero, and on exit they
c           are set to zero.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, real
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      else if( lda.lt.max( 1, n ) )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'cher2 ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     set up the start points in x and y if the increments are not both
c     unity.
c
      if( ( incx.ne.1 ).or.( incy.ne.1 ) )then
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            ky = 1
         else
            ky = 1 - ( n - 1 )*incy
         end if
         jx = kx
         jy = ky
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through the triangular part
c     of a.
c
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when a is stored in the upper triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 20, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1 = alpha*conjg( y( j ) )
                  temp2 = conjg( alpha*x( j ) )
                  do 10, i = 1, j - 1
                     a( i, j ) = a( i, j ) + x( i )*temp1 + y( i )*temp2
   10             continue
                  a( j, j ) = real( a( j, j ) ) +
     $                        real( x( j )*temp1 + y( j )*temp2 )
               else
                  a( j, j ) = real( a( j, j ) )
               end if
   20       continue
         else
            do 40, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1 = alpha*conjg( y( jy ) )
                  temp2 = conjg( alpha*x( jx ) )
                  ix    = kx
                  iy    = ky
                  do 30, i = 1, j - 1
                     a( i, j ) = a( i, j ) + x( ix )*temp1
     $                                     + y( iy )*temp2
                     ix        = ix        + incx
                     iy        = iy        + incy
   30             continue
                  a( j, j ) = real( a( j, j ) ) +
     $                        real( x( jx )*temp1 + y( jy )*temp2 )
               else
                  a( j, j ) = real( a( j, j ) )
               end if
               jx = jx + incx
               jy = jy + incy
   40       continue
         end if
      else
c
c        form  a  when a is stored in the lower triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1     = alpha*conjg( y( j ) )
                  temp2     = conjg( alpha*x( j ) )
                  a( j, j ) = real( a( j, j ) ) +
     $                        real( x( j )*temp1 + y( j )*temp2 )
                  do 50, i = j + 1, n
                     a( i, j ) = a( i, j ) + x( i )*temp1 + y( i )*temp2
   50             continue
               else
                  a( j, j ) = real( a( j, j ) )
               end if
   60       continue
         else
            do 80, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1     = alpha*conjg( y( jy ) )
                  temp2     = conjg( alpha*x( jx ) )
                  a( j, j ) = real( a( j, j ) ) +
     $                        real( x( jx )*temp1 + y( jy )*temp2 )
                  ix        = jx
                  iy        = jy
                  do 70, i = j + 1, n
                     ix        = ix        + incx
                     iy        = iy        + incy
                     a( i, j ) = a( i, j ) + x( ix )*temp1
     $                                     + y( iy )*temp2
   70             continue
               else
                  a( j, j ) = real( a( j, j ) )
               end if
               jx = jx + incx
               jy = jy + incy
   80       continue
         end if
      end if
c
      return
c
c     end of cher2 .
c
      end
c
c***********************************************************************
c
      subroutine chpmv ( uplo, n, alpha, ap, x, incx, beta, y, incy )
c     .. scalar arguments ..
      complex            alpha, beta
      integer            incx, incy, n
      character*1        uplo
c     .. array arguments ..
      complex            ap( * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  chpmv  performs the matrix-vector operation
c
c     y := alpha*a*x + beta*y,
c
c  where alpha and beta are scalars, x and y are n element vectors and
c  a is an n by n hermitian matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the matrix a is supplied in the packed
c           array ap as follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  supplied in ap.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  supplied in ap.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  ap     - complex          array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with uplo = 'u' or 'u', the array ap must
c           contain the upper triangular part of the hermitian matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
c           and a( 2, 2 ) respectively, and so on.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular part of the hermitian matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
c           and a( 3, 1 ) respectively, and so on.
c           note that the imaginary parts of the diagonal elements need
c           not be set and are assumed to be zero.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry, beta specifies the scalar beta. when beta is
c           supplied as zero then y need not be set on input.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y. on exit, y is overwritten by the updated
c           vector y.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, real
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 6
      else if( incy.eq.0 )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'chpmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     set up the start points in  x  and  y.
c
      if( incx.gt.0 )then
         kx = 1
      else
         kx = 1 - ( n - 1 )*incx
      end if
      if( incy.gt.0 )then
         ky = 1
      else
         ky = 1 - ( n - 1 )*incy
      end if
c
c     start the operations. in this version the elements of the array ap
c     are accessed sequentially with one pass through ap.
c
c     first form  y := beta*y.
c
      if( beta.ne.one )then
         if( incy.eq.1 )then
            if( beta.eq.zero )then
               do 10, i = 1, n
                  y( i ) = zero
   10          continue
            else
               do 20, i = 1, n
                  y( i ) = beta*y( i )
   20          continue
            end if
         else
            iy = ky
            if( beta.eq.zero )then
               do 30, i = 1, n
                  y( iy ) = zero
                  iy      = iy   + incy
   30          continue
            else
               do 40, i = 1, n
                  y( iy ) = beta*y( iy )
                  iy      = iy           + incy
   40          continue
            end if
         end if
      end if
      if( alpha.eq.zero )
     $   return
      kk = 1
      if( lsame( uplo, 'u' ) )then
c
c        form  y  when ap contains the upper triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               temp1 = alpha*x( j )
               temp2 = zero
               k     = kk
               do 50, i = 1, j - 1
                  y( i ) = y( i ) + temp1*ap( k )
                  temp2  = temp2  + conjg( ap( k ) )*x( i )
                  k      = k      + 1
   50          continue
               y( j ) = y( j ) + temp1*real( ap( kk + j - 1 ) )
     $                         + alpha*temp2
               kk     = kk     + j
   60       continue
         else
            jx = kx
            jy = ky
            do 80, j = 1, n
               temp1 = alpha*x( jx )
               temp2 = zero
               ix    = kx
               iy    = ky
               do 70, k = kk, kk + j - 2
                  y( iy ) = y( iy ) + temp1*ap( k )
                  temp2   = temp2   + conjg( ap( k ) )*x( ix )
                  ix      = ix      + incx
                  iy      = iy      + incy
   70          continue
               y( jy ) = y( jy ) + temp1*real( ap( kk + j - 1 ) )
     $                           + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
               kk      = kk      + j
   80       continue
         end if
      else
c
c        form  y  when ap contains the lower triangle.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 100, j = 1, n
               temp1  = alpha*x( j )
               temp2  = zero
               y( j ) = y( j ) + temp1*real( ap( kk ) )
               k      = kk     + 1
               do 90, i = j + 1, n
                  y( i ) = y( i ) + temp1*ap( k )
                  temp2  = temp2  + conjg( ap( k ) )*x( i )
                  k      = k      + 1
   90          continue
               y( j ) = y( j ) + alpha*temp2
               kk     = kk     + ( n - j + 1 )
  100       continue
         else
            jx = kx
            jy = ky
            do 120, j = 1, n
               temp1   = alpha*x( jx )
               temp2   = zero
               y( jy ) = y( jy ) + temp1*real( ap( kk ) )
               ix      = jx
               iy      = jy
               do 110, k = kk + 1, kk + n - j
                  ix      = ix      + incx
                  iy      = iy      + incy
                  y( iy ) = y( iy ) + temp1*ap( k )
                  temp2   = temp2   + conjg( ap( k ) )*x( ix )
  110          continue
               y( jy ) = y( jy ) + alpha*temp2
               jx      = jx      + incx
               jy      = jy      + incy
               kk      = kk      + ( n - j + 1 )
  120       continue
         end if
      end if
c
      return
c
c     end of chpmv .
c
      end
c
c***********************************************************************
c
      subroutine chpr  ( uplo, n, alpha, x, incx, ap )
c     .. scalar arguments ..
      real               alpha
      integer            incx, n
      character*1        uplo
c     .. array arguments ..
      complex            ap( * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  chpr    performs the hermitian rank 1 operation
c
c     a := alpha*x*conjg( x' ) + a,
c
c  where alpha is a real scalar, x is an n element vector and a is an
c  n by n hermitian matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the matrix a is supplied in the packed
c           array ap as follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  supplied in ap.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  supplied in ap.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  ap     - complex          array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular part of the hermitian matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
c           and a( 2, 2 ) respectively, and so on. on exit, the array
c           ap is overwritten by the upper triangular part of the
c           updated matrix.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular part of the hermitian matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
c           and a( 3, 1 ) respectively, and so on. on exit, the array
c           ap is overwritten by the lower triangular part of the
c           updated matrix.
c           note that the imaginary parts of the diagonal elements need
c           not be set, they are assumed to be zero, and on exit they
c           are set to zero.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, k, kk, kx
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, real
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      end if
      if( info.ne.0 )then
         call xerbla( 'chpr  ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.real( zero ) ) )
     $   return
c
c     set the start point in x if the increment is not unity.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of the array ap
c     are accessed sequentially with one pass through ap.
c
      kk = 1
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when upper triangle is stored in ap.
c
         if( incx.eq.1 )then
            do 20, j = 1, n
               if( x( j ).ne.zero )then
                  temp = alpha*conjg( x( j ) )
                  k    = kk
                  do 10, i = 1, j - 1
                     ap( k ) = ap( k ) + x( i )*temp
                     k       = k       + 1
   10             continue
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
     $                               + real( x( j )*temp )
               else
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
               end if
               kk = kk + j
   20       continue
         else
            jx = kx
            do 40, j = 1, n
               if( x( jx ).ne.zero )then
                  temp = alpha*conjg( x( jx ) )
                  ix   = kx
                  do 30, k = kk, kk + j - 2
                     ap( k ) = ap( k ) + x( ix )*temp
                     ix      = ix      + incx
   30             continue
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
     $                               + real( x( jx )*temp )
               else
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
               end if
               jx = jx + incx
               kk = kk + j
   40       continue
         end if
      else
c
c        form  a  when lower triangle is stored in ap.
c
         if( incx.eq.1 )then
            do 60, j = 1, n
               if( x( j ).ne.zero )then
                  temp     = alpha*conjg( x( j ) )
                  ap( kk ) = real( ap( kk ) ) + real( temp*x( j ) )
                  k        = kk               + 1
                  do 50, i = j + 1, n
                     ap( k ) = ap( k ) + x( i )*temp
                     k       = k       + 1
   50             continue
               else
                  ap( kk ) = real( ap( kk ) )
               end if
               kk = kk + n - j + 1
   60       continue
         else
            jx = kx
            do 80, j = 1, n
               if( x( jx ).ne.zero )then
                  temp    = alpha*conjg( x( jx ) )
                  ap( kk ) = real( ap( kk ) ) + real( temp*x( jx ) )
                  ix      = jx
                  do 70, k = kk + 1, kk + n - j
                     ix      = ix      + incx
                     ap( k ) = ap( k ) + x( ix )*temp
   70             continue
               else
                  ap( kk ) = real( ap( kk ) )
               end if
               jx = jx + incx
               kk = kk + n - j + 1
   80       continue
         end if
      end if
c
      return
c
c     end of chpr  .
c
      end
c
c***********************************************************************
c
      subroutine chpr2 ( uplo, n, alpha, x, incx, y, incy, ap )
c     .. scalar arguments ..
      complex            alpha
      integer            incx, incy, n
      character*1        uplo
c     .. array arguments ..
      complex            ap( * ), x( * ), y( * )
c     ..
c
c  purpose
c  =======
c
c  chpr2  performs the hermitian rank 2 operation
c
c     a := alpha*x*conjg( y' ) + conjg( alpha )*y*conjg( x' ) + a,
c
c  where alpha is a scalar, x and y are n element vectors and a is an
c  n by n hermitian matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the upper or lower
c           triangular part of the matrix a is supplied in the packed
c           array ap as follows:
c
c              uplo = 'u' or 'u'   the upper triangular part of a is
c                                  supplied in ap.
c
c              uplo = 'l' or 'l'   the lower triangular part of a is
c                                  supplied in ap.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x.
c           unchanged on exit.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c  y      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incy ) ).
c           before entry, the incremented array y must contain the n
c           element vector y.
c           unchanged on exit.
c
c  incy   - integer.
c           on entry, incy specifies the increment for the elements of
c           y. incy must not be zero.
c           unchanged on exit.
c
c  ap     - complex          array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular part of the hermitian matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 1, 2 )
c           and a( 2, 2 ) respectively, and so on. on exit, the array
c           ap is overwritten by the upper triangular part of the
c           updated matrix.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular part of the hermitian matrix
c           packed sequentially, column by column, so that ap( 1 )
c           contains a( 1, 1 ), ap( 2 ) and ap( 3 ) contain a( 2, 1 )
c           and a( 3, 1 ) respectively, and so on. on exit, the array
c           ap is overwritten by the lower triangular part of the
c           updated matrix.
c           note that the imaginary parts of the diagonal elements need
c           not be set, they are assumed to be zero, and on exit they
c           are set to zero.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp1, temp2
      integer            i, info, ix, iy, j, jx, jy, k, kk, kx, ky
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, real
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo, 'u' ).and.
     $         .not.lsame( uplo, 'l' )      )then
         info = 1
      else if( n.lt.0 )then
         info = 2
      else if( incx.eq.0 )then
         info = 5
      else if( incy.eq.0 )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'chpr2 ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.( alpha.eq.zero ) )
     $   return
c
c     set up the start points in x and y if the increments are not both
c     unity.
c
      if( ( incx.ne.1 ).or.( incy.ne.1 ) )then
         if( incx.gt.0 )then
            kx = 1
         else
            kx = 1 - ( n - 1 )*incx
         end if
         if( incy.gt.0 )then
            ky = 1
         else
            ky = 1 - ( n - 1 )*incy
         end if
         jx = kx
         jy = ky
      end if
c
c     start the operations. in this version the elements of the array ap
c     are accessed sequentially with one pass through ap.
c
      kk = 1
      if( lsame( uplo, 'u' ) )then
c
c        form  a  when upper triangle is stored in ap.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 20, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1 = alpha*conjg( y( j ) )
                  temp2 = conjg( alpha*x( j ) )
                  k     = kk
                  do 10, i = 1, j - 1
                     ap( k ) = ap( k ) + x( i )*temp1 + y( i )*temp2
                     k       = k       + 1
   10             continue
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) ) +
     $                               real( x( j )*temp1 + y( j )*temp2 )
               else
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
               end if
               kk = kk + j
   20       continue
         else
            do 40, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1 = alpha*conjg( y( jy ) )
                  temp2 = conjg( alpha*x( jx ) )
                  ix    = kx
                  iy    = ky
                  do 30, k = kk, kk + j - 2
                     ap( k ) = ap( k ) + x( ix )*temp1 + y( iy )*temp2
                     ix      = ix      + incx
                     iy      = iy      + incy
   30             continue
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) ) +
     $                               real( x( jx )*temp1 +
     $                                     y( jy )*temp2 )
               else
                  ap( kk + j - 1 ) = real( ap( kk + j - 1 ) )
               end if
               jx = jx + incx
               jy = jy + incy
               kk = kk + j
   40       continue
         end if
      else
c
c        form  a  when lower triangle is stored in ap.
c
         if( ( incx.eq.1 ).and.( incy.eq.1 ) )then
            do 60, j = 1, n
               if( ( x( j ).ne.zero ).or.( y( j ).ne.zero ) )then
                  temp1   = alpha*conjg( y( j ) )
                  temp2   = conjg( alpha*x( j ) )
                  ap( kk ) = real( ap( kk ) ) +
     $                       real( x( j )*temp1 + y( j )*temp2 )
                  k        = kk               + 1
                  do 50, i = j + 1, n
                     ap( k ) = ap( k ) + x( i )*temp1 + y( i )*temp2
                     k       = k       + 1
   50             continue
               else
                  ap( kk ) = real( ap( kk ) )
               end if
               kk = kk + n - j + 1
   60       continue
         else
            do 80, j = 1, n
               if( ( x( jx ).ne.zero ).or.( y( jy ).ne.zero ) )then
                  temp1    = alpha*conjg( y( jy ) )
                  temp2    = conjg( alpha*x( jx ) )
                  ap( kk ) = real( ap( kk ) ) +
     $                       real( x( jx )*temp1 + y( jy )*temp2 )
                  ix       = jx
                  iy       = jy
                  do 70, k = kk + 1, kk + n - j
                     ix      = ix      + incx
                     iy      = iy      + incy
                     ap( k ) = ap( k ) + x( ix )*temp1 + y( iy )*temp2
   70             continue
               else
                  ap( kk ) = real( ap( kk ) )
               end if
               jx = jx + incx
               jy = jy + incy
               kk = kk + n - j + 1
   80       continue
         end if
      end if
c
      return
c
c     end of chpr2 .
c
      end
c
c***********************************************************************
c
      subroutine ctbmv ( uplo, trans, diag, n, k, a, lda, x, incx )
c     .. scalar arguments ..
      integer            incx, k, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  ctbmv  performs one of the matrix-vector operations
c
c     x := a*x,   or   x := a'*x,   or   x := conjg( a' )*x,
c
c  where x is an n element vector and  a is an n by n unit, or non-unit,
c  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   x := a*x.
c
c              trans = 't' or 't'   x := a'*x.
c
c              trans = 'c' or 'c'   x := conjg( a' )*x.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with uplo = 'u' or 'u', k specifies the number of
c           super-diagonals of the matrix a.
c           on entry with uplo = 'l' or 'l', k specifies the number of
c           sub-diagonals of the matrix a.
c           k must satisfy  0 .le. k.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
c           by n part of the array a must contain the upper triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row
c           ( k + 1 ) of the array, the first super-diagonal starting at
c           position 2 in row k, and so on. the top left k by k triangle
c           of the array a is not referenced.
c           the following program segment will transfer an upper
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 do 20, j = 1, n
c                    m = k + 1 - j
c                    do 10, i = max( 1, j - k ), j
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           before entry with uplo = 'l' or 'l', the leading ( k + 1 )
c           by n part of the array a must contain the lower triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row 1 of
c           the array, the first sub-diagonal starting at position 1 in
c           row 2, and so on. the bottom right k by k triangle of the
c           array a is not referenced.
c           the following program segment will transfer a lower
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 do 20, j = 1, n
c                    m = 1 - j
c                    do 10, i = j, min( n, j + k )
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           note that when diag = 'u' or 'u' the elements of the array a
c           corresponding to the diagonal elements of the matrix are not
c           referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( k + 1 ).
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x. on exit, x is overwritten with the
c           tranformed vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, kplus1, kx, l
      logical            noconj, nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, min
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( k.lt.0 )then
         info = 5
      else if( lda.lt.( k + 1 ) )then
         info = 7
      else if( incx.eq.0 )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'ctbmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      noconj = lsame( trans, 't' )
      nounit = lsame( diag , 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx   too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c         form  x := a*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     l    = kplus1 - j
                     do 10, i = max( 1, j - k ), j - 1
                        x( i ) = x( i ) + temp*a( l + i, j )
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( kplus1, j )
                  end if
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     l    = kplus1  - j
                     do 30, i = max( 1, j - k ), j - 1
                        x( ix ) = x( ix ) + temp*a( l + i, j )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( kplus1, j )
                  end if
                  jx = jx + incx
                  if( j.gt.k )
     $               kx = kx + incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     l    = 1      - j
                     do 50, i = min( n, j + k ), j + 1, -1
                        x( i ) = x( i ) + temp*a( l + i, j )
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( 1, j )
                  end if
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     l    = 1       - j
                     do 70, i = min( n, j + k ), j + 1, -1
                        x( ix ) = x( ix ) + temp*a( l + i, j )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( 1, j )
                  end if
                  jx = jx - incx
                  if( ( n - j ).ge.k )
     $               kx = kx - incx
   80          continue
            end if
         end if
      else
c
c        form  x := a'*x  or  x := conjg( a' )*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 110, j = n, 1, -1
                  temp = x( j )
                  l    = kplus1 - j
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( kplus1, j )
                     do 90, i = j - 1, max( 1, j - k ), -1
                        temp = temp + a( l + i, j )*x( i )
   90                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( kplus1, j ) )
                     do 100, i = j - 1, max( 1, j - k ), -1
                        temp = temp + conjg( a( l + i, j ) )*x( i )
  100                continue
                  end if
                  x( j ) = temp
  110          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 140, j = n, 1, -1
                  temp = x( jx )
                  kx   = kx      - incx
                  ix   = kx
                  l    = kplus1  - j
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( kplus1, j )
                     do 120, i = j - 1, max( 1, j - k ), -1
                        temp = temp + a( l + i, j )*x( ix )
                        ix   = ix   - incx
  120                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( kplus1, j ) )
                     do 130, i = j - 1, max( 1, j - k ), -1
                        temp = temp + conjg( a( l + i, j ) )*x( ix )
                        ix   = ix   - incx
  130                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   - incx
  140          continue
            end if
         else
            if( incx.eq.1 )then
               do 170, j = 1, n
                  temp = x( j )
                  l    = 1      - j
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( 1, j )
                     do 150, i = j + 1, min( n, j + k )
                        temp = temp + a( l + i, j )*x( i )
  150                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( 1, j ) )
                     do 160, i = j + 1, min( n, j + k )
                        temp = temp + conjg( a( l + i, j ) )*x( i )
  160                continue
                  end if
                  x( j ) = temp
  170          continue
            else
               jx = kx
               do 200, j = 1, n
                  temp = x( jx )
                  kx   = kx      + incx
                  ix   = kx
                  l    = 1       - j
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( 1, j )
                     do 180, i = j + 1, min( n, j + k )
                        temp = temp + a( l + i, j )*x( ix )
                        ix   = ix   + incx
  180                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( 1, j ) )
                     do 190, i = j + 1, min( n, j + k )
                        temp = temp + conjg( a( l + i, j ) )*x( ix )
                        ix   = ix   + incx
  190                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   + incx
  200          continue
            end if
         end if
      end if
c
      return
c
c     end of ctbmv .
c
      end
c
c***********************************************************************
c
      subroutine ctbsv ( uplo, trans, diag, n, k, a, lda, x, incx )
c     .. scalar arguments ..
      integer            incx, k, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  ctbsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,   or   conjg( a' )*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular band matrix, with ( k + 1 )
c  diagonals.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   conjg( a' )*x = b.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with uplo = 'u' or 'u', k specifies the number of
c           super-diagonals of the matrix a.
c           on entry with uplo = 'l' or 'l', k specifies the number of
c           sub-diagonals of the matrix a.
c           k must satisfy  0 .le. k.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry with uplo = 'u' or 'u', the leading ( k + 1 )
c           by n part of the array a must contain the upper triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row
c           ( k + 1 ) of the array, the first super-diagonal starting at
c           position 2 in row k, and so on. the top left k by k triangle
c           of the array a is not referenced.
c           the following program segment will transfer an upper
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 do 20, j = 1, n
c                    m = k + 1 - j
c                    do 10, i = max( 1, j - k ), j
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           before entry with uplo = 'l' or 'l', the leading ( k + 1 )
c           by n part of the array a must contain the lower triangular
c           band part of the matrix of coefficients, supplied column by
c           column, with the leading diagonal of the matrix in row 1 of
c           the array, the first sub-diagonal starting at position 1 in
c           row 2, and so on. the bottom right k by k triangle of the
c           array a is not referenced.
c           the following program segment will transfer a lower
c           triangular band matrix from conventional full matrix storage
c           to band storage:
c
c                 do 20, j = 1, n
c                    m = 1 - j
c                    do 10, i = j, min( n, j + k )
c                       a( m + i, j ) = matrix( i, j )
c              10    continue
c              20 continue
c
c           note that when diag = 'u' or 'u' the elements of the array a
c           corresponding to the diagonal elements of the matrix are not
c           referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           ( k + 1 ).
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, kplus1, kx, l
      logical            noconj, nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, min
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( k.lt.0 )then
         info = 5
      else if( lda.lt.( k + 1 ) )then
         info = 7
      else if( incx.eq.0 )then
         info = 9
      end if
      if( info.ne.0 )then
         call xerbla( 'ctbsv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      noconj = lsame( trans, 't' )
      nounit = lsame( diag , 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed by sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := inv( a )*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     l = kplus1 - j
                     if( nounit )
     $                  x( j ) = x( j )/a( kplus1, j )
                     temp = x( j )
                     do 10, i = j - 1, max( 1, j - k ), -1
                        x( i ) = x( i ) - temp*a( l + i, j )
   10                continue
                  end if
   20          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 40, j = n, 1, -1
                  kx = kx - incx
                  if( x( jx ).ne.zero )then
                     ix = kx
                     l  = kplus1 - j
                     if( nounit )
     $                  x( jx ) = x( jx )/a( kplus1, j )
                     temp = x( jx )
                     do 30, i = j - 1, max( 1, j - k ), -1
                        x( ix ) = x( ix ) - temp*a( l + i, j )
                        ix      = ix      - incx
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     l = 1 - j
                     if( nounit )
     $                  x( j ) = x( j )/a( 1, j )
                     temp = x( j )
                     do 50, i = j + 1, min( n, j + k )
                        x( i ) = x( i ) - temp*a( l + i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  kx = kx + incx
                  if( x( jx ).ne.zero )then
                     ix = kx
                     l  = 1  - j
                     if( nounit )
     $                  x( jx ) = x( jx )/a( 1, j )
                     temp = x( jx )
                     do 70, i = j + 1, min( n, j + k )
                        x( ix ) = x( ix ) - temp*a( l + i, j )
                        ix      = ix      + incx
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else
c
c        form  x := inv( a' )*x  or  x := inv( conjg( a') )*x.
c
         if( lsame( uplo, 'u' ) )then
            kplus1 = k + 1
            if( incx.eq.1 )then
               do 110, j = 1, n
                  temp = x( j )
                  l    = kplus1 - j
                  if( noconj )then
                     do 90, i = max( 1, j - k ), j - 1
                        temp = temp - a( l + i, j )*x( i )
   90                continue
                     if( nounit )
     $                  temp = temp/a( kplus1, j )
                  else
                     do 100, i = max( 1, j - k ), j - 1
                        temp = temp - conjg( a( l + i, j ) )*x( i )
  100                continue
                     if( nounit )
     $                  temp = temp/conjg( a( kplus1, j ) )
                  end if
                  x( j ) = temp
  110          continue
            else
               jx = kx
               do 140, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  l    = kplus1  - j
                  if( noconj )then
                     do 120, i = max( 1, j - k ), j - 1
                        temp = temp - a( l + i, j )*x( ix )
                        ix   = ix   + incx
  120                continue
                     if( nounit )
     $                  temp = temp/a( kplus1, j )
                  else
                     do 130, i = max( 1, j - k ), j - 1
                        temp = temp - conjg( a( l + i, j ) )*x( ix )
                        ix   = ix   + incx
  130                continue
                     if( nounit )
     $                  temp = temp/conjg( a( kplus1, j ) )
                  end if
                  x( jx ) = temp
                  jx      = jx   + incx
                  if( j.gt.k )
     $               kx = kx + incx
  140          continue
            end if
         else
            if( incx.eq.1 )then
               do 170, j = n, 1, -1
                  temp = x( j )
                  l    = 1      - j
                  if( noconj )then
                     do 150, i = min( n, j + k ), j + 1, -1
                        temp = temp - a( l + i, j )*x( i )
  150                continue
                     if( nounit )
     $                  temp = temp/a( 1, j )
                  else
                     do 160, i = min( n, j + k ), j + 1, -1
                        temp = temp - conjg( a( l + i, j ) )*x( i )
  160                continue
                     if( nounit )
     $                  temp = temp/conjg( a( 1, j ) )
                  end if
                  x( j ) = temp
  170          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 200, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  l    = 1       - j
                  if( noconj )then
                     do 180, i = min( n, j + k ), j + 1, -1
                        temp = temp - a( l + i, j )*x( ix )
                        ix   = ix   - incx
  180                continue
                     if( nounit )
     $                  temp = temp/a( 1, j )
                  else
                     do 190, i = min( n, j + k ), j + 1, -1
                        temp = temp - conjg( a( l + i, j ) )*x( ix )
                        ix   = ix   - incx
  190                continue
                     if( nounit )
     $                  temp = temp/conjg( a( 1, j ) )
                  end if
                  x( jx ) = temp
                  jx      = jx   - incx
                  if( ( n - j ).ge.k )
     $               kx = kx - incx
  200          continue
            end if
         end if
      end if
c
      return
c
c     end of ctbsv .
c
      end
c
c***********************************************************************
c
      subroutine ctpmv ( uplo, trans, diag, n, ap, x, incx )
c     .. scalar arguments ..
      integer            incx, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      complex            ap( * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  ctpmv  performs one of the matrix-vector operations
c
c     x := a*x,   or   x := a'*x,   or   x := conjg( a' )*x,
c
c  where x is an n element vector and  a is an n by n unit, or non-unit,
c  upper or lower triangular matrix, supplied in packed form.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   x := a*x.
c
c              trans = 't' or 't'   x := a'*x.
c
c              trans = 'c' or 'c'   x := conjg( a' )*x.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  ap     - complex          array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )
c           respectively, and so on.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )
c           respectively, and so on.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x. on exit, x is overwritten with the
c           tranformed vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, k, kk, kx
      logical            noconj, nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( incx.eq.0 )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'ctpmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      noconj = lsame( trans, 't' )
      nounit = lsame( diag , 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of ap are
c     accessed sequentially with one pass through ap.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x:= a*x.
c
         if( lsame( uplo, 'u' ) )then
            kk = 1
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     k    = kk
                     do 10, i = 1, j - 1
                        x( i ) = x( i ) + temp*ap( k )
                        k      = k      + 1
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*ap( kk + j - 1 )
                  end if
                  kk = kk + j
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 30, k = kk, kk + j - 2
                        x( ix ) = x( ix ) + temp*ap( k )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*ap( kk + j - 1 )
                  end if
                  jx = jx + incx
                  kk = kk + j
   40          continue
            end if
         else
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     k    = kk
                     do 50, i = n, j + 1, -1
                        x( i ) = x( i ) + temp*ap( k )
                        k      = k      - 1
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*ap( kk - n + j )
                  end if
                  kk = kk - ( n - j + 1 )
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 70, k = kk, kk - ( n - ( j + 1 ) ), -1
                        x( ix ) = x( ix ) + temp*ap( k )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*ap( kk - n + j )
                  end if
                  jx = jx - incx
                  kk = kk - ( n - j + 1 )
   80          continue
            end if
         end if
      else
c
c        form  x := a'*x  or  x := conjg( a' )*x.
c
         if( lsame( uplo, 'u' ) )then
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 110, j = n, 1, -1
                  temp = x( j )
                  k    = kk     - 1
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*ap( kk )
                     do 90, i = j - 1, 1, -1
                        temp = temp + ap( k )*x( i )
                        k    = k    - 1
   90                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( ap( kk ) )
                     do 100, i = j - 1, 1, -1
                        temp = temp + conjg( ap( k ) )*x( i )
                        k    = k    - 1
  100                continue
                  end if
                  x( j ) = temp
                  kk     = kk   - j
  110          continue
            else
               jx = kx + ( n - 1 )*incx
               do 140, j = n, 1, -1
                  temp = x( jx )
                  ix   = jx
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*ap( kk )
                     do 120, k = kk - 1, kk - j + 1, -1
                        ix   = ix   - incx
                        temp = temp + ap( k )*x( ix )
  120                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( ap( kk ) )
                     do 130, k = kk - 1, kk - j + 1, -1
                        ix   = ix   - incx
                        temp = temp + conjg( ap( k ) )*x( ix )
  130                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   - incx
                  kk      = kk   - j
  140          continue
            end if
         else
            kk = 1
            if( incx.eq.1 )then
               do 170, j = 1, n
                  temp = x( j )
                  k    = kk     + 1
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*ap( kk )
                     do 150, i = j + 1, n
                        temp = temp + ap( k )*x( i )
                        k    = k    + 1
  150                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( ap( kk ) )
                     do 160, i = j + 1, n
                        temp = temp + conjg( ap( k ) )*x( i )
                        k    = k    + 1
  160                continue
                  end if
                  x( j ) = temp
                  kk     = kk   + ( n - j + 1 )
  170          continue
            else
               jx = kx
               do 200, j = 1, n
                  temp = x( jx )
                  ix   = jx
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*ap( kk )
                     do 180, k = kk + 1, kk + n - j
                        ix   = ix   + incx
                        temp = temp + ap( k )*x( ix )
  180                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( ap( kk ) )
                     do 190, k = kk + 1, kk + n - j
                        ix   = ix   + incx
                        temp = temp + conjg( ap( k ) )*x( ix )
  190                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   + incx
                  kk      = kk   + ( n - j + 1 )
  200          continue
            end if
         end if
      end if
c
      return
c
c     end of ctpmv .
c
      end
c
c***********************************************************************
c
      subroutine ctpsv ( uplo, trans, diag, n, ap, x, incx )
c     .. scalar arguments ..
      integer            incx, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      complex            ap( * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  ctpsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,   or   conjg( a' )*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix, supplied in packed form.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   conjg( a' )*x = b.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  ap     - complex          array of dimension at least
c           ( ( n*( n + 1 ) )/2 ).
c           before entry with  uplo = 'u' or 'u', the array ap must
c           contain the upper triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 1, 2 ) and a( 2, 2 )
c           respectively, and so on.
c           before entry with uplo = 'l' or 'l', the array ap must
c           contain the lower triangular matrix packed sequentially,
c           column by column, so that ap( 1 ) contains a( 1, 1 ),
c           ap( 2 ) and ap( 3 ) contain a( 2, 1 ) and a( 3, 1 )
c           respectively, and so on.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced, but are assumed to be unity.
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, k, kk, kx
      logical            noconj, nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( incx.eq.0 )then
         info = 7
      end if
      if( info.ne.0 )then
         call xerbla( 'ctpsv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      noconj = lsame( trans, 't' )
      nounit = lsame( diag , 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of ap are
c     accessed sequentially with one pass through ap.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := inv( a )*x.
c
         if( lsame( uplo, 'u' ) )then
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/ap( kk )
                     temp = x( j )
                     k    = kk     - 1
                     do 10, i = j - 1, 1, -1
                        x( i ) = x( i ) - temp*ap( k )
                        k      = k      - 1
   10                continue
                  end if
                  kk = kk - j
   20          continue
            else
               jx = kx + ( n - 1 )*incx
               do 40, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/ap( kk )
                     temp = x( jx )
                     ix   = jx
                     do 30, k = kk - 1, kk - j + 1, -1
                        ix      = ix      - incx
                        x( ix ) = x( ix ) - temp*ap( k )
   30                continue
                  end if
                  jx = jx - incx
                  kk = kk - j
   40          continue
            end if
         else
            kk = 1
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/ap( kk )
                     temp = x( j )
                     k    = kk     + 1
                     do 50, i = j + 1, n
                        x( i ) = x( i ) - temp*ap( k )
                        k      = k      + 1
   50                continue
                  end if
                  kk = kk + ( n - j + 1 )
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/ap( kk )
                     temp = x( jx )
                     ix   = jx
                     do 70, k = kk + 1, kk + n - j
                        ix      = ix      + incx
                        x( ix ) = x( ix ) - temp*ap( k )
   70                continue
                  end if
                  jx = jx + incx
                  kk = kk + ( n - j + 1 )
   80          continue
            end if
         end if
      else
c
c        form  x := inv( a' )*x  or  x := inv( conjg( a' ) )*x.
c
         if( lsame( uplo, 'u' ) )then
            kk = 1
            if( incx.eq.1 )then
               do 110, j = 1, n
                  temp = x( j )
                  k    = kk
                  if( noconj )then
                     do 90, i = 1, j - 1
                        temp = temp - ap( k )*x( i )
                        k    = k    + 1
   90                continue
                     if( nounit )
     $                  temp = temp/ap( kk + j - 1 )
                  else
                     do 100, i = 1, j - 1
                        temp = temp - conjg( ap( k ) )*x( i )
                        k    = k    + 1
  100                continue
                     if( nounit )
     $                  temp = temp/conjg( ap( kk + j - 1 ) )
                  end if
                  x( j ) = temp
                  kk     = kk   + j
  110          continue
            else
               jx = kx
               do 140, j = 1, n
                  temp = x( jx )
                  ix   = kx
                  if( noconj )then
                     do 120, k = kk, kk + j - 2
                        temp = temp - ap( k )*x( ix )
                        ix   = ix   + incx
  120                continue
                     if( nounit )
     $                  temp = temp/ap( kk + j - 1 )
                  else
                     do 130, k = kk, kk + j - 2
                        temp = temp - conjg( ap( k ) )*x( ix )
                        ix   = ix   + incx
  130                continue
                     if( nounit )
     $                  temp = temp/conjg( ap( kk + j - 1 ) )
                  end if
                  x( jx ) = temp
                  jx      = jx   + incx
                  kk      = kk   + j
  140          continue
            end if
         else
            kk = ( n*( n + 1 ) )/2
            if( incx.eq.1 )then
               do 170, j = n, 1, -1
                  temp = x( j )
                  k    = kk
                  if( noconj )then
                     do 150, i = n, j + 1, -1
                        temp = temp - ap( k )*x( i )
                        k    = k    - 1
  150                continue
                     if( nounit )
     $                  temp = temp/ap( kk - n + j )
                  else
                     do 160, i = n, j + 1, -1
                        temp = temp - conjg( ap( k ) )*x( i )
                        k    = k    - 1
  160                continue
                     if( nounit )
     $                  temp = temp/conjg( ap( kk - n + j ) )
                  end if
                  x( j ) = temp
                  kk     = kk   - ( n - j + 1 )
  170          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 200, j = n, 1, -1
                  temp = x( jx )
                  ix   = kx
                  if( noconj )then
                     do 180, k = kk, kk - ( n - ( j + 1 ) ), -1
                        temp = temp - ap( k )*x( ix )
                        ix   = ix   - incx
  180                continue
                     if( nounit )
     $                  temp = temp/ap( kk - n + j )
                  else
                     do 190, k = kk, kk - ( n - ( j + 1 ) ), -1
                        temp = temp - conjg( ap( k ) )*x( ix )
                        ix   = ix   - incx
  190                continue
                     if( nounit )
     $                  temp = temp/conjg( ap( kk - n + j ) )
                  end if
                  x( jx ) = temp
                  jx      = jx   - incx
                  kk      = kk   - ( n - j + 1 )
  200          continue
            end if
         end if
      end if
c
      return
c
c     end of ctpsv .
c
      end
c
c***********************************************************************
c
      subroutine ctrmv ( uplo, trans, diag, n, a, lda, x, incx )
c     .. scalar arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  ctrmv  performs one of the matrix-vector operations
c
c     x := a*x,   or   x := a'*x,   or   x := conjg( a' )*x,
c
c  where x is an n element vector and  a is an n by n unit, or non-unit,
c  upper or lower triangular matrix.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   x := a*x.
c
c              trans = 't' or 't'   x := a'*x.
c
c              trans = 'c' or 'c'   x := conjg( a' )*x.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element vector x. on exit, x is overwritten with the
c           tranformed vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, kx
      logical            noconj, nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call xerbla( 'ctrmv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      noconj = lsame( trans, 't' )
      nounit = lsame( diag , 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := a*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 20, j = 1, n
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 10, i = 1, j - 1
                        x( i ) = x( i ) + temp*a( i, j )
   10                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   20          continue
            else
               jx = kx
               do 40, j = 1, n
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 30, i = 1, j - 1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      + incx
   30                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx + incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     temp = x( j )
                     do 50, i = n, j + 1, -1
                        x( i ) = x( i ) + temp*a( i, j )
   50                continue
                     if( nounit )
     $                  x( j ) = x( j )*a( j, j )
                  end if
   60          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 80, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     temp = x( jx )
                     ix   = kx
                     do 70, i = n, j + 1, -1
                        x( ix ) = x( ix ) + temp*a( i, j )
                        ix      = ix      - incx
   70                continue
                     if( nounit )
     $                  x( jx ) = x( jx )*a( j, j )
                  end if
                  jx = jx - incx
   80          continue
            end if
         end if
      else
c
c        form  x := a'*x  or  x := conjg( a' )*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 110, j = n, 1, -1
                  temp = x( j )
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( j, j )
                     do 90, i = j - 1, 1, -1
                        temp = temp + a( i, j )*x( i )
   90                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( j, j ) )
                     do 100, i = j - 1, 1, -1
                        temp = temp + conjg( a( i, j ) )*x( i )
  100                continue
                  end if
                  x( j ) = temp
  110          continue
            else
               jx = kx + ( n - 1 )*incx
               do 140, j = n, 1, -1
                  temp = x( jx )
                  ix   = jx
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( j, j )
                     do 120, i = j - 1, 1, -1
                        ix   = ix   - incx
                        temp = temp + a( i, j )*x( ix )
  120                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( j, j ) )
                     do 130, i = j - 1, 1, -1
                        ix   = ix   - incx
                        temp = temp + conjg( a( i, j ) )*x( ix )
  130                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   - incx
  140          continue
            end if
         else
            if( incx.eq.1 )then
               do 170, j = 1, n
                  temp = x( j )
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( j, j )
                     do 150, i = j + 1, n
                        temp = temp + a( i, j )*x( i )
  150                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( j, j ) )
                     do 160, i = j + 1, n
                        temp = temp + conjg( a( i, j ) )*x( i )
  160                continue
                  end if
                  x( j ) = temp
  170          continue
            else
               jx = kx
               do 200, j = 1, n
                  temp = x( jx )
                  ix   = jx
                  if( noconj )then
                     if( nounit )
     $                  temp = temp*a( j, j )
                     do 180, i = j + 1, n
                        ix   = ix   + incx
                        temp = temp + a( i, j )*x( ix )
  180                continue
                  else
                     if( nounit )
     $                  temp = temp*conjg( a( j, j ) )
                     do 190, i = j + 1, n
                        ix   = ix   + incx
                        temp = temp + conjg( a( i, j ) )*x( ix )
  190                continue
                  end if
                  x( jx ) = temp
                  jx      = jx   + incx
  200          continue
            end if
         end if
      end if
c
      return
c
c     end of ctrmv .
c
      end
c
c***********************************************************************
c
      subroutine ctrsv ( uplo, trans, diag, n, a, lda, x, incx )
c     .. scalar arguments ..
      integer            incx, lda, n
      character*1        diag, trans, uplo
c     .. array arguments ..
      complex            a( lda, * ), x( * )
c     ..
c
c  purpose
c  =======
c
c  ctrsv  solves one of the systems of equations
c
c     a*x = b,   or   a'*x = b,   or   conjg( a' )*x = b,
c
c  where b and x are n element vectors and a is an n by n unit, or
c  non-unit, upper or lower triangular matrix.
c
c  no test for singularity or near-singularity is included in this
c  routine. such tests must be performed before calling this routine.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry, trans specifies the equations to be solved as
c           follows:
c
c              trans = 'n' or 'n'   a*x = b.
c
c              trans = 't' or 't'   a'*x = b.
c
c              trans = 'c' or 'c'   conjg( a' )*x = b.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit
c           triangular as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the order of the matrix a.
c           n must be at least zero.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, n ).
c           before entry with  uplo = 'u' or 'u', the leading n by n
c           upper triangular part of the array a must contain the upper
c           triangular matrix and the strictly lower triangular part of
c           a is not referenced.
c           before entry with uplo = 'l' or 'l', the leading n by n
c           lower triangular part of the array a must contain the lower
c           triangular matrix and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u', the diagonal elements of
c           a are not referenced either, but are assumed to be unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. lda must be at least
c           max( 1, n ).
c           unchanged on exit.
c
c  x      - complex          array of dimension at least
c           ( 1 + ( n - 1 )*abs( incx ) ).
c           before entry, the incremented array x must contain the n
c           element right-hand side vector b. on exit, x is overwritten
c           with the solution vector x.
c
c  incx   - integer.
c           on entry, incx specifies the increment for the elements of
c           x. incx must not be zero.
c           unchanged on exit.
c
c
c  level 2 blas routine.
c
c  -- written on 22-october-1986.
c     jack dongarra, argonne national lab.
c     jeremy du croz, nag central office.
c     sven hammarling, nag central office.
c     richard hanson, sandia national labs.
c
c
c     .. parameters ..
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     .. local scalars ..
      complex            temp
      integer            i, info, ix, j, jx, kx
      logical            noconj, nounit
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      info = 0
      if     ( .not.lsame( uplo , 'u' ).and.
     $         .not.lsame( uplo , 'l' )      )then
         info = 1
      else if( .not.lsame( trans, 'n' ).and.
     $         .not.lsame( trans, 't' ).and.
     $         .not.lsame( trans, 'c' )      )then
         info = 2
      else if( .not.lsame( diag , 'u' ).and.
     $         .not.lsame( diag , 'n' )      )then
         info = 3
      else if( n.lt.0 )then
         info = 4
      else if( lda.lt.max( 1, n ) )then
         info = 6
      else if( incx.eq.0 )then
         info = 8
      end if
      if( info.ne.0 )then
         call xerbla( 'ctrsv ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
      noconj = lsame( trans, 't' )
      nounit = lsame( diag , 'n' )
c
c     set up the start point in x if the increment is not unity. this
c     will be  ( n - 1 )*incx  too small for descending loops.
c
      if( incx.le.0 )then
         kx = 1 - ( n - 1 )*incx
      else if( incx.ne.1 )then
         kx = 1
      end if
c
c     start the operations. in this version the elements of a are
c     accessed sequentially with one pass through a.
c
      if( lsame( trans, 'n' ) )then
c
c        form  x := inv( a )*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 20, j = n, 1, -1
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 10, i = j - 1, 1, -1
                        x( i ) = x( i ) - temp*a( i, j )
   10                continue
                  end if
   20          continue
            else
               jx = kx + ( n - 1 )*incx
               do 40, j = n, 1, -1
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 30, i = j - 1, 1, -1
                        ix      = ix      - incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   30                continue
                  end if
                  jx = jx - incx
   40          continue
            end if
         else
            if( incx.eq.1 )then
               do 60, j = 1, n
                  if( x( j ).ne.zero )then
                     if( nounit )
     $                  x( j ) = x( j )/a( j, j )
                     temp = x( j )
                     do 50, i = j + 1, n
                        x( i ) = x( i ) - temp*a( i, j )
   50                continue
                  end if
   60          continue
            else
               jx = kx
               do 80, j = 1, n
                  if( x( jx ).ne.zero )then
                     if( nounit )
     $                  x( jx ) = x( jx )/a( j, j )
                     temp = x( jx )
                     ix   = jx
                     do 70, i = j + 1, n
                        ix      = ix      + incx
                        x( ix ) = x( ix ) - temp*a( i, j )
   70                continue
                  end if
                  jx = jx + incx
   80          continue
            end if
         end if
      else
c
c        form  x := inv( a' )*x  or  x := inv( conjg( a' ) )*x.
c
         if( lsame( uplo, 'u' ) )then
            if( incx.eq.1 )then
               do 110, j = 1, n
                  temp = x( j )
                  if( noconj )then
                     do 90, i = 1, j - 1
                        temp = temp - a( i, j )*x( i )
   90                continue
                     if( nounit )
     $                  temp = temp/a( j, j )
                  else
                     do 100, i = 1, j - 1
                        temp = temp - conjg( a( i, j ) )*x( i )
  100                continue
                     if( nounit )
     $                  temp = temp/conjg( a( j, j ) )
                  end if
                  x( j ) = temp
  110          continue
            else
               jx = kx
               do 140, j = 1, n
                  ix   = kx
                  temp = x( jx )
                  if( noconj )then
                     do 120, i = 1, j - 1
                        temp = temp - a( i, j )*x( ix )
                        ix   = ix   + incx
  120                continue
                     if( nounit )
     $                  temp = temp/a( j, j )
                  else
                     do 130, i = 1, j - 1
                        temp = temp - conjg( a( i, j ) )*x( ix )
                        ix   = ix   + incx
  130                continue
                     if( nounit )
     $                  temp = temp/conjg( a( j, j ) )
                  end if
                  x( jx ) = temp
                  jx      = jx   + incx
  140          continue
            end if
         else
            if( incx.eq.1 )then
               do 170, j = n, 1, -1
                  temp = x( j )
                  if( noconj )then
                     do 150, i = n, j + 1, -1
                        temp = temp - a( i, j )*x( i )
  150                continue
                     if( nounit )
     $                  temp = temp/a( j, j )
                  else
                     do 160, i = n, j + 1, -1
                        temp = temp - conjg( a( i, j ) )*x( i )
  160                continue
                     if( nounit )
     $                  temp = temp/conjg( a( j, j ) )
                  end if
                  x( j ) = temp
  170          continue
            else
               kx = kx + ( n - 1 )*incx
               jx = kx
               do 200, j = n, 1, -1
                  ix   = kx
                  temp = x( jx )
                  if( noconj )then
                     do 180, i = n, j + 1, -1
                        temp = temp - a( i, j )*x( ix )
                        ix   = ix   - incx
  180                continue
                     if( nounit )
     $                  temp = temp/a( j, j )
                  else
                     do 190, i = n, j + 1, -1
                        temp = temp - conjg( a( i, j ) )*x( ix )
                        ix   = ix   - incx
  190                continue
                     if( nounit )
     $                  temp = temp/conjg( a( j, j ) )
                  end if
                  x( jx ) = temp
                  jx      = jx   - incx
  200          continue
            end if
         end if
      end if
c
      return
c
c     end of ctrsv .
c
      end
c***********************************************************************
c                          blas3
c***********************************************************************
c
      subroutine sgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      character*1        transa, transb
      integer            m, n, k, lda, ldb, ldc
      real*8             alpha, beta
c     .. array arguments ..
      real*8             a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  sgemm  performs one of the matrix-matrix operations
c
c     c := alpha*op( a )*op( b ) + beta*c,
c
c  where  op( x ) is one of
c
c     op( x ) = x   or   op( x ) = x',
c
c  alpha and beta are scalars, and a, b and c are matrices, with op( a )
c  an m by k matrix,  op( b )  a  k by n matrix and  c an m by n matrix.
c
c  parameters
c  ==========
c
c  transa - character*1.
c           on entry, transa specifies the form of op( a ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'n' or 'n',  op( a ) = a.
c
c              transa = 't' or 't',  op( a ) = a'.
c
c              transa = 'c' or 'c',  op( a ) = a'.
c
c           unchanged on exit.
c
c  transb - character*1.
c           on entry, transb specifies the form of op( b ) to be used in
c           the matrix multiplication as follows:
c
c              transb = 'n' or 'n',  op( b ) = b.
c
c              transb = 't' or 't',  op( b ) = b'.
c
c              transb = 'c' or 'c',  op( b ) = b'.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry,  m  specifies  the number  of rows  of the  matrix
c           op( a )  and of the  matrix  c.  m  must  be at least  zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n  specifies the number  of columns of the matrix
c           op( b ) and the number of columns of the matrix c. n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry,  k  specifies  the number of columns of the matrix
c           op( a ) and the number of rows of the matrix op( b ). k must
c           be at least  zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, ka ), where ka is
c           k  when  transa = 'n' or 'n',  and is  m  otherwise.
c           before entry with  transa = 'n' or 'n',  the leading  m by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by m  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. when  transa = 'n' or 'n' then
c           lda must be at least  max( 1, m ), otherwise  lda must be at
c           least  max( 1, k ).
c           unchanged on exit.
c
c  b      - real             array of dimension ( ldb, kb ), where kb is
c           n  when  transb = 'n' or 'n',  and is  k  otherwise.
c           before entry with  transb = 'n' or 'n',  the leading  k by n
c           part of the array  b  must contain the matrix  b,  otherwise
c           the leading  n by k  part of the array  b  must contain  the
c           matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in the calling (sub) program. when  transb = 'n' or 'n' then
c           ldb must be at least  max( 1, k ), otherwise  ldb must be at
c           least  max( 1, n ).
c           unchanged on exit.
c
c  beta   - real            .
c           on entry,  beta  specifies the scalar  beta.  when  beta  is
c           supplied as zero then c need not be set on input.
c           unchanged on exit.
c
c  c      - real             array of dimension ( ldc, n ).
c           before entry, the leading  m by n  part of the array  c must
c           contain the matrix  c,  except when  beta  is zero, in which
c           case c need not be set on entry.
c           on exit, the array  c  is overwritten by the  m by n  matrix
c           ( alpha*op( a )*op( b ) + beta*c ).
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            nota, notb
      integer            i, info, j, l, ncola, nrowa, nrowb
      real*8             temp
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     ..
c     .. executable statements ..
c
c     set  nota  and  notb  as  true if  a  and  b  respectively are not
c     transposed and set  nrowa, ncola and  nrowb  as the number of rows
c     and  columns of  a  and the  number of  rows  of  b  respectively.
c
      nota  = lsame( transa, 'n' )
      notb  = lsame( transb, 'n' )
      if( nota )then
         nrowa = m
         ncola = k
      else
         nrowa = k
         ncola = m
      end if
      if( notb )then
         nrowb = k
      else
         nrowb = n
      end if
c
c     test the input parameters.
c
      info = 0
      if(      ( .not.nota                 ).and.
     $         ( .not.lsame( transa, 'c' ) ).and.
     $         ( .not.lsame( transa, 't' ) )      )then
         info = 1
      else if( ( .not.notb                 ).and.
     $         ( .not.lsame( transb, 'c' ) ).and.
     $         ( .not.lsame( transb, 't' ) )      )then
         info = 2
      else if( m  .lt.0               )then
         info = 3
      else if( n  .lt.0               )then
         info = 4
      else if( k  .lt.0               )then
         info = 5
      else if( lda.lt.max( 1, nrowa ) )then
         info = 8
      else if( ldb.lt.max( 1, nrowb ) )then
         info = 10
      else if( ldc.lt.max( 1, m     ) )then
         info = 13
      end if
      if( info.ne.0 )then
         call xerbla( 'sgemm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
c
c     and if  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( beta.eq.zero )then
            do 20, j = 1, n
               do 10, i = 1, m
                  c( i, j ) = zero
   10          continue
   20       continue
         else
            do 40, j = 1, n
               do 30, i = 1, m
                  c( i, j ) = beta*c( i, j )
   30          continue
   40       continue
         end if
         return
      end if
c
c     start the operations.
c
      if( notb )then
         if( nota )then
c
c           form  c := alpha*a*b + beta*c.
c
            do 90, j = 1, n
               if( beta.eq.zero )then
                  do 50, i = 1, m
                     c( i, j ) = zero
   50             continue
               else if( beta.ne.one )then
                  do 60, i = 1, m
                     c( i, j ) = beta*c( i, j )
   60             continue
               end if
               do 80, l = 1, k
                  if( b( l, j ).ne.zero )then
                     temp = alpha*b( l, j )
                     do 70, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
   70                continue
                  end if
   80          continue
   90       continue
         else
c
c           form  c := alpha*a'*b + beta*c
c
            do 120, j = 1, n
               do 110, i = 1, m
                  temp = zero
                  do 100, l = 1, k
                     temp = temp + a( l, i )*b( l, j )
  100             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  110          continue
  120       continue
         end if
      else
         if( nota )then
c
c           form  c := alpha*a*b' + beta*c
c
            do 170, j = 1, n
               if( beta.eq.zero )then
                  do 130, i = 1, m
                     c( i, j ) = zero
  130             continue
               else if( beta.ne.one )then
                  do 140, i = 1, m
                     c( i, j ) = beta*c( i, j )
  140             continue
               end if
               do 160, l = 1, k
                  if( b( j, l ).ne.zero )then
                     temp = alpha*b( j, l )
                     do 150, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  150                continue
                  end if
  160          continue
  170       continue
         else
c
c           form  c := alpha*a'*b' + beta*c
c
            do 200, j = 1, n
               do 190, i = 1, m
                  temp = zero
                  do 180, l = 1, k
                     temp = temp + a( l, i )*b( j, l )
  180             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  190          continue
  200       continue
         end if
      end if
c
      return
c
c     end of sgemm .
c
      end
c
c***********************************************************************
c
      subroutine ssymm ( side, uplo, m, n, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      character*1        side, uplo
      integer            m, n, lda, ldb, ldc
      real*8             alpha, beta
c     .. array arguments ..
      real*8             a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  ssymm  performs one of the matrix-matrix operations
c
c     c := alpha*a*b + beta*c,
c
c  or
c
c     c := alpha*b*a + beta*c,
c
c  where alpha and beta are scalars,  a is a symmetric matrix and  b and
c  c are  m by n matrices.
c
c  parameters
c  ==========
c
c  side   - character*1.
c           on entry,  side  specifies whether  the  symmetric matrix  a
c           appears on the  left or right  in the  operation as follows:
c
c              side = 'l' or 'l'   c := alpha*a*b + beta*c,
c
c              side = 'r' or 'r'   c := alpha*b*a + beta*c,
c
c           unchanged on exit.
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  symmetric  matrix   a  is  to  be
c           referenced as follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of the
c                                  symmetric matrix is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of the
c                                  symmetric matrix is to be referenced.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry,  m  specifies the number of rows of the matrix  c.
c           m  must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix c.
c           n  must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, ka ), where ka is
c           m  when  side = 'l' or 'l'  and is  n otherwise.
c           before entry  with  side = 'l' or 'l',  the  m by m  part of
c           the array  a  must contain the  symmetric matrix,  such that
c           when  uplo = 'u' or 'u', the leading m by m upper triangular
c           part of the array  a  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
c           the leading  m by m  lower triangular part  of the  array  a
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  a  is not
c           referenced.
c           before entry  with  side = 'r' or 'r',  the  n by n  part of
c           the array  a  must contain the  symmetric matrix,  such that
c           when  uplo = 'u' or 'u', the leading n by n upper triangular
c           part of the array  a  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
c           the leading  n by n  lower triangular part  of the  array  a
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  a  is not
c           referenced.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program.  when  side = 'l' or 'l'  then
c           lda must be at least  max( 1, m ), otherwise  lda must be at
c           least  max( 1, n ).
c           unchanged on exit.
c
c  b      - real             array of dimension ( ldb, n ).
c           before entry, the leading  m by n part of the array  b  must
c           contain the matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   ldb  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c  beta   - real            .
c           on entry,  beta  specifies the scalar  beta.  when  beta  is
c           supplied as zero then c need not be set on input.
c           unchanged on exit.
c
c  c      - real             array of dimension ( ldc, n ).
c           before entry, the leading  m by n  part of the array  c must
c           contain the matrix  c,  except when  beta  is zero, in which
c           case c need not be set on entry.
c           on exit, the array  c  is overwritten by the  m by n updated
c           matrix.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            upper
      integer            i, info, j, k, nrowa
      real*8             temp1, temp2
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     ..
c     .. executable statements ..
c
c     set nrowa as the number of rows of a.
c
      if( lsame( side, 'l' ) )then
         nrowa = m
      else
         nrowa = n
      end if
      upper = lsame( uplo, 'u' )
c
c     test the input parameters.
c
      info = 0
      if(      ( .not.lsame( side, 'l' ) ).and.
     $         ( .not.lsame( side, 'r' ) )      )then
         info = 1
      else if( ( .not.upper              ).and.
     $         ( .not.lsame( uplo, 'l' ) )      )then
         info = 2
      else if( m  .lt.0               )then
         info = 3
      else if( n  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldb.lt.max( 1, m     ) )then
         info = 9
      else if( ldc.lt.max( 1, m     ) )then
         info = 12
      end if
      if( info.ne.0 )then
         call xerbla( 'ssymm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( beta.eq.zero )then
            do 20, j = 1, n
               do 10, i = 1, m
                  c( i, j ) = zero
   10          continue
   20       continue
         else
            do 40, j = 1, n
               do 30, i = 1, m
                  c( i, j ) = beta*c( i, j )
   30          continue
   40       continue
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( side, 'l' ) )then
c
c        form  c := alpha*a*b + beta*c.
c
         if( upper )then
            do 70, j = 1, n
               do 60, i = 1, m
                  temp1 = alpha*b( i, j )
                  temp2 = zero
                  do 50, k = 1, i - 1
                     c( k, j ) = c( k, j ) + temp1    *a( k, i )
                     temp2     = temp2     + b( k, j )*a( k, i )
   50             continue
                  if( beta.eq.zero )then
                     c( i, j ) = temp1*a( i, i ) + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           temp1*a( i, i ) + alpha*temp2
                  end if
   60          continue
   70       continue
         else
            do 100, j = 1, n
               do 90, i = m, 1, -1
                  temp1 = alpha*b( i, j )
                  temp2 = zero
                  do 80, k = i + 1, m
                     c( k, j ) = c( k, j ) + temp1    *a( k, i )
                     temp2     = temp2     + b( k, j )*a( k, i )
   80             continue
                  if( beta.eq.zero )then
                     c( i, j ) = temp1*a( i, i ) + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           temp1*a( i, i ) + alpha*temp2
                  end if
   90          continue
  100       continue
         end if
      else
c
c        form  c := alpha*b*a + beta*c.
c
         do 170, j = 1, n
            temp1 = alpha*a( j, j )
            if( beta.eq.zero )then
               do 110, i = 1, m
                  c( i, j ) = temp1*b( i, j )
  110          continue
            else
               do 120, i = 1, m
                  c( i, j ) = beta*c( i, j ) + temp1*b( i, j )
  120          continue
            end if
            do 140, k = 1, j - 1
               if( upper )then
                  temp1 = alpha*a( k, j )
               else
                  temp1 = alpha*a( j, k )
               end if
               do 130, i = 1, m
                  c( i, j ) = c( i, j ) + temp1*b( i, k )
  130          continue
  140       continue
            do 160, k = j + 1, n
               if( upper )then
                  temp1 = alpha*a( j, k )
               else
                  temp1 = alpha*a( k, j )
               end if
               do 150, i = 1, m
                  c( i, j ) = c( i, j ) + temp1*b( i, k )
  150          continue
  160       continue
  170    continue
      end if
c
      return
c
c     end of ssymm .
c
      end
c
c***********************************************************************
c
      subroutine ssyrk ( uplo, trans, n, k, alpha, a, lda,
     $                   beta, c, ldc )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      character*1        uplo, trans
      integer            n, k, lda, ldc
      real*8             alpha, beta
c     .. array arguments ..
      real*8             a( lda, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  ssyrk  performs one of the symmetric rank k operations
c
c     c := alpha*a*a' + beta*c,
c
c  or
c
c     c := alpha*a'*a + beta*c,
c
c  where  alpha and beta  are scalars, c is an  n by n  symmetric matrix
c  and  a  is an  n by k  matrix in the first case and a  k by n  matrix
c  in the second case.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  c  is to be  referenced  as
c           follows:
c
c              uplo = 'u' or 'u'   only the  upper triangular part of  c
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the  lower triangular part of  c
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry,  trans  specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   c := alpha*a*a' + beta*c.
c
c              trans = 't' or 't'   c := alpha*a'*a + beta*c.
c
c              trans = 'c' or 'c'   c := alpha*a'*a + beta*c.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n specifies the order of the matrix c.  n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with  trans = 'n' or 'n',  k  specifies  the number
c           of  columns   of  the   matrix   a,   and  on   entry   with
c           trans = 't' or 't' or 'c' or 'c',  k  specifies  the  number
c           of rows of the matrix  a.  k must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, ka ), where ka is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by n  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  lda must be at least  max( 1, n ), otherwise  lda must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  c      - real             array of dimension ( ldc, n ).
c           before entry  with  uplo = 'u' or 'u',  the leading  n by n
c           upper triangular part of the array c must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of c is not referenced.  on exit, the
c           upper triangular part of the array  c is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry  with  uplo = 'l' or 'l',  the leading  n by n
c           lower triangular part of the array c must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of c is not referenced.  on exit, the
c           lower triangular part of the array  c is overwritten by the
c           lower triangular part of the updated matrix.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            upper
      integer            i, info, j, l, nrowa
      real*8             temp
c     .. parameters ..
      real*8             one ,         zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      if( lsame( trans, 'n' ) )then
         nrowa = n
      else
         nrowa = k
      end if
      upper = lsame( uplo, 'u' )
c
      info = 0
      if(      ( .not.upper               ).and.
     $         ( .not.lsame( uplo , 'l' ) )      )then
         info = 1
      else if( ( .not.lsame( trans, 'n' ) ).and.
     $         ( .not.lsame( trans, 't' ) ).and.
     $         ( .not.lsame( trans, 'c' ) )      )then
         info = 2
      else if( n  .lt.0               )then
         info = 3
      else if( k  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldc.lt.max( 1, n     ) )then
         info = 10
      end if
      if( info.ne.0 )then
         call xerbla( 'ssyrk ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( upper )then
            if( beta.eq.zero )then
               do 20, j = 1, n
                  do 10, i = 1, j
                     c( i, j ) = zero
   10             continue
   20          continue
            else
               do 40, j = 1, n
                  do 30, i = 1, j
                     c( i, j ) = beta*c( i, j )
   30             continue
   40          continue
            end if
         else
            if( beta.eq.zero )then
               do 60, j = 1, n
                  do 50, i = j, n
                     c( i, j ) = zero
   50             continue
   60          continue
            else
               do 80, j = 1, n
                  do 70, i = j, n
                     c( i, j ) = beta*c( i, j )
   70             continue
   80          continue
            end if
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( trans, 'n' ) )then
c
c        form  c := alpha*a*a' + beta*c.
c
         if( upper )then
            do 130, j = 1, n
               if( beta.eq.zero )then
                  do 90, i = 1, j
                     c( i, j ) = zero
   90             continue
               else if( beta.ne.one )then
                  do 100, i = 1, j
                     c( i, j ) = beta*c( i, j )
  100             continue
               end if
               do 120, l = 1, k
                  if( a( j, l ).ne.zero )then
                     temp = alpha*a( j, l )
                     do 110, i = 1, j
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  110                continue
                  end if
  120          continue
  130       continue
         else
            do 180, j = 1, n
               if( beta.eq.zero )then
                  do 140, i = j, n
                     c( i, j ) = zero
  140             continue
               else if( beta.ne.one )then
                  do 150, i = j, n
                     c( i, j ) = beta*c( i, j )
  150             continue
               end if
               do 170, l = 1, k
                  if( a( j, l ).ne.zero )then
                     temp      = alpha*a( j, l )
                     do 160, i = j, n
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  160                continue
                  end if
  170          continue
  180       continue
         end if
      else
c
c        form  c := alpha*a'*a + beta*c.
c
         if( upper )then
            do 210, j = 1, n
               do 200, i = 1, j
                  temp = zero
                  do 190, l = 1, k
                     temp = temp + a( l, i )*a( l, j )
  190             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  200          continue
  210       continue
         else
            do 240, j = 1, n
               do 230, i = j, n
                  temp = zero
                  do 220, l = 1, k
                     temp = temp + a( l, i )*a( l, j )
  220             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  230          continue
  240       continue
         end if
      end if
c
      return
c
c     end of ssyrk .
c
      end
c
c***********************************************************************
c
      subroutine ssyr2k( uplo, trans, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      character*1        uplo, trans
      integer            n, k, lda, ldb, ldc
      real*8             alpha, beta
c     .. array arguments ..
      real*8             a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  ssyr2k  performs one of the symmetric rank 2k operations
c
c     c := alpha*a*b' + alpha*b*a' + beta*c,
c
c  or
c
c     c := alpha*a'*b + alpha*b'*a + beta*c,
c
c  where  alpha and beta  are scalars, c is an  n by n  symmetric matrix
c  and  a and b  are  n by k  matrices  in the  first  case  and  k by n
c  matrices in the second case.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  c  is to be  referenced  as
c           follows:
c
c              uplo = 'u' or 'u'   only the  upper triangular part of  c
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the  lower triangular part of  c
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry,  trans  specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   c := alpha*a*b' + alpha*b*a' +
c                                        beta*c.
c
c              trans = 't' or 't'   c := alpha*a'*b + alpha*b'*a +
c                                        beta*c.
c
c              trans = 'c' or 'c'   c := alpha*a'*b + alpha*b'*a +
c                                        beta*c.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n specifies the order of the matrix c.  n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with  trans = 'n' or 'n',  k  specifies  the number
c           of  columns  of the  matrices  a and b,  and on  entry  with
c           trans = 't' or 't' or 'c' or 'c',  k  specifies  the  number
c           of rows of the matrices  a and b.  k must be at least  zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, ka ), where ka is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by n  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  lda must be at least  max( 1, n ), otherwise  lda must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  b      - real             array of dimension ( ldb, kb ), where kb is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  b  must contain the matrix  b,  otherwise
c           the leading  k by n  part of the array  b  must contain  the
c           matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  ldb must be at least  max( 1, n ), otherwise  ldb must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  c      - real             array of dimension ( ldc, n ).
c           before entry  with  uplo = 'u' or 'u',  the leading  n by n
c           upper triangular part of the array c must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of c is not referenced.  on exit, the
c           upper triangular part of the array  c is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry  with  uplo = 'l' or 'l',  the leading  n by n
c           lower triangular part of the array c must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of c is not referenced.  on exit, the
c           lower triangular part of the array  c is overwritten by the
c           lower triangular part of the updated matrix.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            upper
      integer            i, info, j, l, nrowa
      real*8             temp1, temp2
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      if( lsame( trans, 'n' ) )then
         nrowa = n
      else
         nrowa = k
      end if
      upper = lsame( uplo, 'u' )
c
      info = 0
      if(      ( .not.upper               ).and.
     $         ( .not.lsame( uplo , 'l' ) )      )then
         info = 1
      else if( ( .not.lsame( trans, 'n' ) ).and.
     $         ( .not.lsame( trans, 't' ) ).and.
     $         ( .not.lsame( trans, 'c' ) )      )then
         info = 2
      else if( n  .lt.0               )then
         info = 3
      else if( k  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldb.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldc.lt.max( 1, n     ) )then
         info = 12
      end if
      if( info.ne.0 )then
         call xerbla( 'ssyr2k', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( upper )then
            if( beta.eq.zero )then
               do 20, j = 1, n
                  do 10, i = 1, j
                     c( i, j ) = zero
   10             continue
   20          continue
            else
               do 40, j = 1, n
                  do 30, i = 1, j
                     c( i, j ) = beta*c( i, j )
   30             continue
   40          continue
            end if
         else
            if( beta.eq.zero )then
               do 60, j = 1, n
                  do 50, i = j, n
                     c( i, j ) = zero
   50             continue
   60          continue
            else
               do 80, j = 1, n
                  do 70, i = j, n
                     c( i, j ) = beta*c( i, j )
   70             continue
   80          continue
            end if
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( trans, 'n' ) )then
c
c        form  c := alpha*a*b' + alpha*b*a' + c.
c
         if( upper )then
            do 130, j = 1, n
               if( beta.eq.zero )then
                  do 90, i = 1, j
                     c( i, j ) = zero
   90             continue
               else if( beta.ne.one )then
                  do 100, i = 1, j
                     c( i, j ) = beta*c( i, j )
  100             continue
               end if
               do 120, l = 1, k
                  if( ( a( j, l ).ne.zero ).or.
     $                ( b( j, l ).ne.zero )     )then
                     temp1 = alpha*b( j, l )
                     temp2 = alpha*a( j, l )
                     do 110, i = 1, j
                        c( i, j ) = c( i, j ) +
     $                              a( i, l )*temp1 + b( i, l )*temp2
  110                continue
                  end if
  120          continue
  130       continue
         else
            do 180, j = 1, n
               if( beta.eq.zero )then
                  do 140, i = j, n
                     c( i, j ) = zero
  140             continue
               else if( beta.ne.one )then
                  do 150, i = j, n
                     c( i, j ) = beta*c( i, j )
  150             continue
               end if
               do 170, l = 1, k
                  if( ( a( j, l ).ne.zero ).or.
     $                ( b( j, l ).ne.zero )     )then
                     temp1 = alpha*b( j, l )
                     temp2 = alpha*a( j, l )
                     do 160, i = j, n
                        c( i, j ) = c( i, j ) +
     $                              a( i, l )*temp1 + b( i, l )*temp2
  160                continue
                  end if
  170          continue
  180       continue
         end if
      else
c
c        form  c := alpha*a'*b + alpha*b'*a + c.
c
         if( upper )then
            do 210, j = 1, n
               do 200, i = 1, j
                  temp1 = zero
                  temp2 = zero
                  do 190, l = 1, k
                     temp1 = temp1 + a( l, i )*b( l, j )
                     temp2 = temp2 + b( l, i )*a( l, j )
  190             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp1 + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           alpha*temp1 + alpha*temp2
                  end if
  200          continue
  210       continue
         else
            do 240, j = 1, n
               do 230, i = j, n
                  temp1 = zero
                  temp2 = zero
                  do 220, l = 1, k
                     temp1 = temp1 + a( l, i )*b( l, j )
                     temp2 = temp2 + b( l, i )*a( l, j )
  220             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp1 + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           alpha*temp1 + alpha*temp2
                  end if
  230          continue
  240       continue
         end if
      end if
c
      return
c
c     end of ssyr2k.
c
      end
c
c***********************************************************************
c
      subroutine strmm ( side, uplo, transa, diag, m, n, alpha, a, lda,
     $                   b, ldb )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      character*1        side, uplo, transa, diag
      integer            m, n, lda, ldb
      real*8             alpha
c     .. array arguments ..
      real*8             a( lda, * ), b( ldb, * )
c     ..
c
c  purpose
c  =======
c
c  strmm  performs one of the matrix-matrix operations
c
c     b := alpha*op( a )*b,   or   b := alpha*b*op( a ),
c
c  where  alpha  is a scalar,  b  is an m by n matrix,  a  is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( a )  is one  of
c
c     op( a ) = a   or   op( a ) = a'.
c
c  parameters
c  ==========
c
c  side   - character*1.
c           on entry,  side specifies whether  op( a ) multiplies b from
c           the left or right as follows:
c
c              side = 'l' or 'l'   b := alpha*op( a )*b.
c
c              side = 'r' or 'r'   b := alpha*b*op( a ).
c
c           unchanged on exit.
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix a is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  transa - character*1.
c           on entry, transa specifies the form of op( a ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'n' or 'n'   op( a ) = a.
c
c              transa = 't' or 't'   op( a ) = a'.
c
c              transa = 'c' or 'c'   op( a ) = a'.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit triangular
c           as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of b. m must be at
c           least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of b.  n must be
c           at least zero.
c           unchanged on exit.
c
c  alpha  - real*8          .
c           on entry,  alpha specifies the scalar  alpha. when  alpha is
c           zero then  a is not referenced and  b need not be set before
c           entry.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, k ), where k is m
c           when  side = 'l' or 'l'  and is  n  when  side = 'r' or 'r'.
c           before entry  with  uplo = 'u' or 'u',  the  leading  k by k
c           upper triangular part of the array  a must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           a is not referenced.
c           before entry  with  uplo = 'l' or 'l',  the  leading  k by k
c           lower triangular part of the array  a must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u',  the diagonal elements of
c           a  are not referenced either,  but are assumed to be  unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program.  when  side = 'l' or 'l'  then
c           lda  must be at least  max( 1, m ),  when  side = 'r' or 'r'
c           then lda must be at least max( 1, n ).
c           unchanged on exit.
c
c  b      - real             array of dimension ( ldb, n ).
c           before entry,  the leading  m by n part of the array  b must
c           contain the matrix  b,  and  on exit  is overwritten  by the
c           transformed matrix.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   ldb  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            lside, nounit, upper
      integer            i, info, j, k, nrowa
      real*8             temp
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      lside  = lsame( side  , 'l' )
      if( lside )then
         nrowa = m
      else
         nrowa = n
      end if
      nounit = lsame( diag  , 'n' )
      upper  = lsame( uplo  , 'u' )
c
      info   = 0
      if(      ( .not.lside                ).and.
     $         ( .not.lsame( side  , 'r' ) )      )then
         info = 1
      else if( ( .not.upper                ).and.
     $         ( .not.lsame( uplo  , 'l' ) )      )then
         info = 2
      else if( ( .not.lsame( transa, 'n' ) ).and.
     $         ( .not.lsame( transa, 't' ) ).and.
     $         ( .not.lsame( transa, 'c' ) )      )then
         info = 3
      else if( ( .not.lsame( diag  , 'u' ) ).and.
     $         ( .not.lsame( diag  , 'n' ) )      )then
         info = 4
      else if( m  .lt.0               )then
         info = 5
      else if( n  .lt.0               )then
         info = 6
      else if( lda.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldb.lt.max( 1, m     ) )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'strmm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         do 20, j = 1, n
            do 10, i = 1, m
               b( i, j ) = zero
   10       continue
   20    continue
         return
      end if
c
c     start the operations.
c
      if( lside )then
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*a*b.
c
            if( upper )then
               do 50, j = 1, n
                  do 40, k = 1, m
                     if( b( k, j ).ne.zero )then
                        temp = alpha*b( k, j )
                        do 30, i = 1, k - 1
                           b( i, j ) = b( i, j ) + temp*a( i, k )
   30                   continue
                        if( nounit )
     $                     temp = temp*a( k, k )
                        b( k, j ) = temp
                     end if
   40             continue
   50          continue
            else
               do 80, j = 1, n
                  do 70 k = m, 1, -1
                     if( b( k, j ).ne.zero )then
                        temp      = alpha*b( k, j )
                        b( k, j ) = temp
                        if( nounit )
     $                     b( k, j ) = b( k, j )*a( k, k )
                        do 60, i = k + 1, m
                           b( i, j ) = b( i, j ) + temp*a( i, k )
   60                   continue
                     end if
   70             continue
   80          continue
            end if
         else
c
c           form  b := alpha*b*a'.
c
            if( upper )then
               do 110, j = 1, n
                  do 100, i = m, 1, -1
                     temp = b( i, j )
                     if( nounit )
     $                  temp = temp*a( i, i )
                     do 90, k = 1, i - 1
                        temp = temp + a( k, i )*b( k, j )
   90                continue
                     b( i, j ) = alpha*temp
  100             continue
  110          continue
            else
               do 140, j = 1, n
                  do 130, i = 1, m
                     temp = b( i, j )
                     if( nounit )
     $                  temp = temp*a( i, i )
                     do 120, k = i + 1, m
                        temp = temp + a( k, i )*b( k, j )
  120                continue
                     b( i, j ) = alpha*temp
  130             continue
  140          continue
            end if
         end if
      else
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*b*a.
c
            if( upper )then
               do 180, j = n, 1, -1
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 150, i = 1, m
                     b( i, j ) = temp*b( i, j )
  150             continue
                  do 170, k = 1, j - 1
                     if( a( k, j ).ne.zero )then
                        temp = alpha*a( k, j )
                        do 160, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  160                   continue
                     end if
  170             continue
  180          continue
            else
               do 220, j = 1, n
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 190, i = 1, m
                     b( i, j ) = temp*b( i, j )
  190             continue
                  do 210, k = j + 1, n
                     if( a( k, j ).ne.zero )then
                        temp = alpha*a( k, j )
                        do 200, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  200                   continue
                     end if
  210             continue
  220          continue
            end if
         else
c
c           form  b := alpha*b*a'.
c
            if( upper )then
               do 260, k = 1, n
                  do 240, j = 1, k - 1
                     if( a( j, k ).ne.zero )then
                        temp = alpha*a( j, k )
                        do 230, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  230                   continue
                     end if
  240             continue
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( k, k )
                  if( temp.ne.one )then
                     do 250, i = 1, m
                        b( i, k ) = temp*b( i, k )
  250                continue
                  end if
  260          continue
            else
               do 300, k = n, 1, -1
                  do 280, j = k + 1, n
                     if( a( j, k ).ne.zero )then
                        temp = alpha*a( j, k )
                        do 270, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  270                   continue
                     end if
  280             continue
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( k, k )
                  if( temp.ne.one )then
                     do 290, i = 1, m
                        b( i, k ) = temp*b( i, k )
  290                continue
                  end if
  300          continue
            end if
         end if
      end if
c
      return
c
c     end of strmm .
c
      end
c
c***********************************************************************
c
      subroutine strsm ( side, uplo, transa, diag, m, n, alpha, a, lda,
     $                   b, ldb )
      implicit real*8 (a-h, o-z)
c     .. scalar arguments ..
      character*1        side, uplo, transa, diag
      integer            m, n, lda, ldb
      real*8             alpha
c     .. array arguments ..
      real*8             a( lda, * ), b( ldb, * )
c     ..
c
c  purpose
c  =======
c
c  strsm  solves one of the matrix equations
c
c     op( a )*x = alpha*b,   or   x*op( a ) = alpha*b,
c
c  where alpha is a scalar, x and b are m by n matrices, a is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( a )  is one  of
c
c     op( a ) = a   or   op( a ) = a'.
c
c  the matrix x is overwritten on b.
c
c  parameters
c  ==========
c
c  side   - character*1.
c           on entry, side specifies whether op( a ) appears on the left
c           or right of x as follows:
c
c              side = 'l' or 'l'   op( a )*x = alpha*b.
c
c              side = 'r' or 'r'   x*op( a ) = alpha*b.
c
c           unchanged on exit.
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix a is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  transa - character*1.
c           on entry, transa specifies the form of op( a ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'n' or 'n'   op( a ) = a.
c
c              transa = 't' or 't'   op( a ) = a'.
c
c              transa = 'c' or 'c'   op( a ) = a'.
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit triangular
c           as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of b. m must be at
c           least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of b.  n must be
c           at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry,  alpha specifies the scalar  alpha. when  alpha is
c           zero then  a is not referenced and  b need not be set before
c           entry.
c           unchanged on exit.
c
c  a      - real             array of dimension ( lda, k ), where k is m
c           when  side = 'l' or 'l'  and is  n  when  side = 'r' or 'r'.
c           before entry  with  uplo = 'u' or 'u',  the  leading  k by k
c           upper triangular part of the array  a must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           a is not referenced.
c           before entry  with  uplo = 'l' or 'l',  the  leading  k by k
c           lower triangular part of the array  a must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u',  the diagonal elements of
c           a  are not referenced either,  but are assumed to be  unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program.  when  side = 'l' or 'l'  then
c           lda  must be at least  max( 1, m ),  when  side = 'r' or 'r'
c           then lda must be at least max( 1, n ).
c           unchanged on exit.
c
c  b      - real             array of dimension ( ldb, n ).
c           before entry,  the leading  m by n part of the array  b must
c           contain  the  right-hand  side  matrix  b,  and  on exit  is
c           overwritten by the solution matrix  x.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   ldb  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            lside, nounit, upper
      integer            i, info, j, k, nrowa
      real*8             temp
c     .. parameters ..
      real*8             one         , zero
      parameter        ( one = 1.0d+0, zero = 0.0d+0 )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      lside  = lsame( side  , 'l' )
      if( lside )then
         nrowa = m
      else
         nrowa = n
      end if
      nounit = lsame( diag  , 'n' )
      upper  = lsame( uplo  , 'u' )
c
      info   = 0
      if(      ( .not.lside                ).and.
     $         ( .not.lsame( side  , 'r' ) )      )then
         info = 1
      else if( ( .not.upper                ).and.
     $         ( .not.lsame( uplo  , 'l' ) )      )then
         info = 2
      else if( ( .not.lsame( transa, 'n' ) ).and.
     $         ( .not.lsame( transa, 't' ) ).and.
     $         ( .not.lsame( transa, 'c' ) )      )then
         info = 3
      else if( ( .not.lsame( diag  , 'u' ) ).and.
     $         ( .not.lsame( diag  , 'n' ) )      )then
         info = 4
      else if( m  .lt.0               )then
         info = 5
      else if( n  .lt.0               )then
         info = 6
      else if( lda.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldb.lt.max( 1, m     ) )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'strsm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         do 20, j = 1, n
            do 10, i = 1, m
               b( i, j ) = zero
   10       continue
   20    continue
         return
      end if
c
c     start the operations.
c
      if( lside )then
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*inv( a )*b.
c
            if( upper )then
               do 60, j = 1, n
                  if( alpha.ne.one )then
                     do 30, i = 1, m
                        b( i, j ) = alpha*b( i, j )
   30                continue
                  end if
                  do 50, k = m, 1, -1
                     if( b( k, j ).ne.zero )then
                        if( nounit )
     $                     b( k, j ) = b( k, j )/a( k, k )
                        do 40, i = 1, k - 1
                           b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
   40                   continue
                     end if
   50             continue
   60          continue
            else
               do 100, j = 1, n
                  if( alpha.ne.one )then
                     do 70, i = 1, m
                        b( i, j ) = alpha*b( i, j )
   70                continue
                  end if
                  do 90 k = 1, m
                     if( b( k, j ).ne.zero )then
                        if( nounit )
     $                     b( k, j ) = b( k, j )/a( k, k )
                        do 80, i = k + 1, m
                           b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
   80                   continue
                     end if
   90             continue
  100          continue
            end if
         else
c
c           form  b := alpha*inv( a' )*b.
c
            if( upper )then
               do 130, j = 1, n
                  do 120, i = 1, m
                     temp = alpha*b( i, j )
                     do 110, k = 1, i - 1
                        temp = temp - a( k, i )*b( k, j )
  110                continue
                     if( nounit )
     $                  temp = temp/a( i, i )
                     b( i, j ) = temp
  120             continue
  130          continue
            else
               do 160, j = 1, n
                  do 150, i = m, 1, -1
                     temp = alpha*b( i, j )
                     do 140, k = i + 1, m
                        temp = temp - a( k, i )*b( k, j )
  140                continue
                     if( nounit )
     $                  temp = temp/a( i, i )
                     b( i, j ) = temp
  150             continue
  160          continue
            end if
         end if
      else
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*b*inv( a ).
c
            if( upper )then
               do 210, j = 1, n
                  if( alpha.ne.one )then
                     do 170, i = 1, m
                        b( i, j ) = alpha*b( i, j )
  170                continue
                  end if
                  do 190, k = 1, j - 1
                     if( a( k, j ).ne.zero )then
                        do 180, i = 1, m
                           b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  180                   continue
                     end if
  190             continue
                  if( nounit )then
                     temp = one/a( j, j )
                     do 200, i = 1, m
                        b( i, j ) = temp*b( i, j )
  200                continue
                  end if
  210          continue
            else
               do 260, j = n, 1, -1
                  if( alpha.ne.one )then
                     do 220, i = 1, m
                        b( i, j ) = alpha*b( i, j )
  220                continue
                  end if
                  do 240, k = j + 1, n
                     if( a( k, j ).ne.zero )then
                        do 230, i = 1, m
                           b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  230                   continue
                     end if
  240             continue
                  if( nounit )then
                     temp = one/a( j, j )
                     do 250, i = 1, m
                       b( i, j ) = temp*b( i, j )
  250                continue
                  end if
  260          continue
            end if
         else
c
c           form  b := alpha*b*inv( a' ).
c
            if( upper )then
               do 310, k = n, 1, -1
                  if( nounit )then
                     temp = one/a( k, k )
                     do 270, i = 1, m
                        b( i, k ) = temp*b( i, k )
  270                continue
                  end if
                  do 290, j = 1, k - 1
                     if( a( j, k ).ne.zero )then
                        temp = a( j, k )
                        do 280, i = 1, m
                           b( i, j ) = b( i, j ) - temp*b( i, k )
  280                   continue
                     end if
  290             continue
                  if( alpha.ne.one )then
                     do 300, i = 1, m
                        b( i, k ) = alpha*b( i, k )
  300                continue
                  end if
  310          continue
            else
               do 360, k = 1, n
                  if( nounit )then
                     temp = one/a( k, k )
                     do 320, i = 1, m
                        b( i, k ) = temp*b( i, k )
  320                continue
                  end if
                  do 340, j = k + 1, n
                     if( a( j, k ).ne.zero )then
                        temp = a( j, k )
                        do 330, i = 1, m
                           b( i, j ) = b( i, j ) - temp*b( i, k )
  330                   continue
                     end if
  340             continue
                  if( alpha.ne.one )then
                     do 350, i = 1, m
                        b( i, k ) = alpha*b( i, k )
  350                continue
                  end if
  360          continue
            end if
         end if
      end if
c
      return
c
c     end of strsm .
c
      end
c***********************************************************************
c
      subroutine cgemm ( transa, transb, m, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
c     .. scalar arguments ..
      character*1        transa, transb
      integer            m, n, k, lda, ldb, ldc
      complex            alpha, beta
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  cgemm  performs one of the matrix-matrix operations
c
c     c := alpha*op( a )*op( b ) + beta*c,
c
c  where  op( x ) is one of
c
c     op( x ) = x   or   op( x ) = x'   or   op( x ) = conjg( x' ),
c
c  alpha and beta are scalars, and a, b and c are matrices, with op( a )
c  an m by k matrix,  op( b )  a  k by n matrix and  c an m by n matrix.
c
c  parameters
c  ==========
c
c  transa - character*1.
c           on entry, transa specifies the form of op( a ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'n' or 'n',  op( a ) = a.
c
c              transa = 't' or 't',  op( a ) = a'.
c
c              transa = 'c' or 'c',  op( a ) = conjg( a' ).
c
c           unchanged on exit.
c
c  transb - character*1.
c           on entry, transb specifies the form of op( b ) to be used in
c           the matrix multiplication as follows:
c
c              transb = 'n' or 'n',  op( b ) = b.
c
c              transb = 't' or 't',  op( b ) = b'.
c
c              transb = 'c' or 'c',  op( b ) = conjg( b' ).
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry,  m  specifies  the number  of rows  of the  matrix
c           op( a )  and of the  matrix  c.  m  must  be at least  zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n  specifies the number  of columns of the matrix
c           op( b ) and the number of columns of the matrix c. n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry,  k  specifies  the number of columns of the matrix
c           op( a ) and the number of rows of the matrix op( b ). k must
c           be at least  zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
c           k  when  transa = 'n' or 'n',  and is  m  otherwise.
c           before entry with  transa = 'n' or 'n',  the leading  m by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by m  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program. when  transa = 'n' or 'n' then
c           lda must be at least  max( 1, m ), otherwise  lda must be at
c           least  max( 1, k ).
c           unchanged on exit.
c
c  b      - complex          array of dimension ( ldb, kb ), where kb is
c           n  when  transb = 'n' or 'n',  and is  k  otherwise.
c           before entry with  transb = 'n' or 'n',  the leading  k by n
c           part of the array  b  must contain the matrix  b,  otherwise
c           the leading  n by k  part of the array  b  must contain  the
c           matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in the calling (sub) program. when  transb = 'n' or 'n' then
c           ldb must be at least  max( 1, k ), otherwise  ldb must be at
c           least  max( 1, n ).
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry,  beta  specifies the scalar  beta.  when  beta  is
c           supplied as zero then c need not be set on input.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
c           before entry, the leading  m by n  part of the array  c must
c           contain the matrix  c,  except when  beta  is zero, in which
c           case c need not be set on entry.
c           on exit, the array  c  is overwritten by the  m by n  matrix
c           ( alpha*op( a )*op( b ) + beta*c ).
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c     .. local scalars ..
      logical            conja, conjb, nota, notb
      integer            i, info, j, l, ncola, nrowa, nrowb
      complex            temp
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     ..
c     .. executable statements ..
c
c     set  nota  and  notb  as  true if  a  and  b  respectively are not
c     conjugated or transposed, set  conja and conjb  as true if  a  and
c     b  respectively are to be  transposed but  not conjugated  and set
c     nrowa, ncola and  nrowb  as the number of rows and  columns  of  a
c     and the number of rows of  b  respectively.
c
      nota  = lsame( transa, 'n' )
      notb  = lsame( transb, 'n' )
      conja = lsame( transa, 'c' )
      conjb = lsame( transb, 'c' )
      if( nota )then
         nrowa = m
         ncola = k
      else
         nrowa = k
         ncola = m
      end if
      if( notb )then
         nrowb = k
      else
         nrowb = n
      end if
c
c     test the input parameters.
c
      info = 0
      if(      ( .not.nota                 ).and.
     $         ( .not.conja                ).and.
     $         ( .not.lsame( transa, 't' ) )      )then
         info = 1
      else if( ( .not.notb                 ).and.
     $         ( .not.conjb                ).and.
     $         ( .not.lsame( transb, 't' ) )      )then
         info = 2
      else if( m  .lt.0               )then
         info = 3
      else if( n  .lt.0               )then
         info = 4
      else if( k  .lt.0               )then
         info = 5
      else if( lda.lt.max( 1, nrowa ) )then
         info = 8
      else if( ldb.lt.max( 1, nrowb ) )then
         info = 10
      else if( ldc.lt.max( 1, m     ) )then
         info = 13
      end if
      if( info.ne.0 )then
         call xerbla( 'cgemm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( beta.eq.zero )then
            do 20, j = 1, n
               do 10, i = 1, m
                  c( i, j ) = zero
   10          continue
   20       continue
         else
            do 40, j = 1, n
               do 30, i = 1, m
                  c( i, j ) = beta*c( i, j )
   30          continue
   40       continue
         end if
         return
      end if
c
c     start the operations.
c
      if( notb )then
         if( nota )then
c
c           form  c := alpha*a*b + beta*c.
c
            do 90, j = 1, n
               if( beta.eq.zero )then
                  do 50, i = 1, m
                     c( i, j ) = zero
   50             continue
               else if( beta.ne.one )then
                  do 60, i = 1, m
                     c( i, j ) = beta*c( i, j )
   60             continue
               end if
               do 80, l = 1, k
                  if( b( l, j ).ne.zero )then
                     temp = alpha*b( l, j )
                     do 70, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
   70                continue
                  end if
   80          continue
   90       continue
         else if( conja )then
c
c           form  c := alpha*conjg( a' )*b + beta*c.
c
            do 120, j = 1, n
               do 110, i = 1, m
                  temp = zero
                  do 100, l = 1, k
                     temp = temp + conjg( a( l, i ) )*b( l, j )
  100             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  110          continue
  120       continue
         else
c
c           form  c := alpha*a'*b + beta*c
c
            do 150, j = 1, n
               do 140, i = 1, m
                  temp = zero
                  do 130, l = 1, k
                     temp = temp + a( l, i )*b( l, j )
  130             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  140          continue
  150       continue
         end if
      else if( nota )then
         if( conjb )then
c
c           form  c := alpha*a*conjg( b' ) + beta*c.
c
            do 200, j = 1, n
               if( beta.eq.zero )then
                  do 160, i = 1, m
                     c( i, j ) = zero
  160             continue
               else if( beta.ne.one )then
                  do 170, i = 1, m
                     c( i, j ) = beta*c( i, j )
  170             continue
               end if
               do 190, l = 1, k
                  if( b( j, l ).ne.zero )then
                     temp = alpha*conjg( b( j, l ) )
                     do 180, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  180                continue
                  end if
  190          continue
  200       continue
         else
c
c           form  c := alpha*a*b'          + beta*c
c
            do 250, j = 1, n
               if( beta.eq.zero )then
                  do 210, i = 1, m
                     c( i, j ) = zero
  210             continue
               else if( beta.ne.one )then
                  do 220, i = 1, m
                     c( i, j ) = beta*c( i, j )
  220             continue
               end if
               do 240, l = 1, k
                  if( b( j, l ).ne.zero )then
                     temp = alpha*b( j, l )
                     do 230, i = 1, m
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  230                continue
                  end if
  240          continue
  250       continue
         end if
      else if( conja )then
         if( conjb )then
c
c           form  c := alpha*conjg( a' )*conjg( b' ) + beta*c.
c
            do 280, j = 1, n
               do 270, i = 1, m
                  temp = zero
                  do 260, l = 1, k
                     temp = temp + conjg( a( l, i ) )*conjg( b( j, l ) )
  260             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  270          continue
  280       continue
         else
c
c           form  c := alpha*conjg( a' )*b' + beta*c
c
            do 310, j = 1, n
               do 300, i = 1, m
                  temp = zero
                  do 290, l = 1, k
                     temp = temp + conjg( a( l, i ) )*b( j, l )
  290             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  300          continue
  310       continue
         end if
      else
         if( conjb )then
c
c           form  c := alpha*a'*conjg( b' ) + beta*c
c
            do 340, j = 1, n
               do 330, i = 1, m
                  temp = zero
                  do 320, l = 1, k
                     temp = temp + a( l, i )*conjg( b( j, l ) )
  320             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  330          continue
  340       continue
         else
c
c           form  c := alpha*a'*b' + beta*c
c
            do 370, j = 1, n
               do 360, i = 1, m
                  temp = zero
                  do 350, l = 1, k
                     temp = temp + a( l, i )*b( j, l )
  350             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  360          continue
  370       continue
         end if
      end if
c
      return
c
c     end of cgemm .
c
      end
c
c***********************************************************************
c
      subroutine csymm ( side, uplo, m, n, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
c     .. scalar arguments ..
      character*1        side, uplo
      integer            m, n, lda, ldb, ldc
      complex            alpha, beta
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  csymm  performs one of the matrix-matrix operations
c
c     c := alpha*a*b + beta*c,
c
c  or
c
c     c := alpha*b*a + beta*c,
c
c  where  alpha and beta are scalars, a is a symmetric matrix and  b and
c  c are m by n matrices.
c
c  parameters
c  ==========
c
c  side   - character*1.
c           on entry,  side  specifies whether  the  symmetric matrix  a
c           appears on the  left or right  in the  operation as follows:
c
c              side = 'l' or 'l'   c := alpha*a*b + beta*c,
c
c              side = 'r' or 'r'   c := alpha*b*a + beta*c,
c
c           unchanged on exit.
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  symmetric  matrix   a  is  to  be
c           referenced as follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of the
c                                  symmetric matrix is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of the
c                                  symmetric matrix is to be referenced.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry,  m  specifies the number of rows of the matrix  c.
c           m  must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix c.
c           n  must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
c           m  when  side = 'l' or 'l'  and is n  otherwise.
c           before entry  with  side = 'l' or 'l',  the  m by m  part of
c           the array  a  must contain the  symmetric matrix,  such that
c           when  uplo = 'u' or 'u', the leading m by m upper triangular
c           part of the array  a  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
c           the leading  m by m  lower triangular part  of the  array  a
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  a  is not
c           referenced.
c           before entry  with  side = 'r' or 'r',  the  n by n  part of
c           the array  a  must contain the  symmetric matrix,  such that
c           when  uplo = 'u' or 'u', the leading n by n upper triangular
c           part of the array  a  must contain the upper triangular part
c           of the  symmetric matrix and the  strictly  lower triangular
c           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
c           the leading  n by n  lower triangular part  of the  array  a
c           must  contain  the  lower triangular part  of the  symmetric
c           matrix and the  strictly upper triangular part of  a  is not
c           referenced.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the  calling (sub) program. when  side = 'l' or 'l'  then
c           lda must be at least  max( 1, m ), otherwise  lda must be at
c           least max( 1, n ).
c           unchanged on exit.
c
c  b      - complex          array of dimension ( ldb, n ).
c           before entry, the leading  m by n part of the array  b  must
c           contain the matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   ldb  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry,  beta  specifies the scalar  beta.  when  beta  is
c           supplied as zero then c need not be set on input.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
c           before entry, the leading  m by n  part of the array  c must
c           contain the matrix  c,  except when  beta  is zero, in which
c           case c need not be set on entry.
c           on exit, the array  c  is overwritten by the  m by n updated
c           matrix.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            upper
      integer            i, info, j, k, nrowa
      complex            temp1, temp2
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     ..
c     .. executable statements ..
c
c     set nrowa as the number of rows of a.
c
      if( lsame( side, 'l' ) )then
         nrowa = m
      else
         nrowa = n
      end if
      upper = lsame( uplo, 'u' )
c
c     test the input parameters.
c
      info = 0
      if(      ( .not.lsame( side, 'l' ) ).and.
     $         ( .not.lsame( side, 'r' ) )      )then
         info = 1
      else if( ( .not.upper              ).and.
     $         ( .not.lsame( uplo, 'l' ) )      )then
         info = 2
      else if( m  .lt.0               )then
         info = 3
      else if( n  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldb.lt.max( 1, m     ) )then
         info = 9
      else if( ldc.lt.max( 1, m     ) )then
         info = 12
      end if
      if( info.ne.0 )then
         call xerbla( 'csymm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( beta.eq.zero )then
            do 20, j = 1, n
               do 10, i = 1, m
                  c( i, j ) = zero
   10          continue
   20       continue
         else
            do 40, j = 1, n
               do 30, i = 1, m
                  c( i, j ) = beta*c( i, j )
   30          continue
   40       continue
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( side, 'l' ) )then
c
c        form  c := alpha*a*b + beta*c.
c
         if( upper )then
            do 70, j = 1, n
               do 60, i = 1, m
                  temp1 = alpha*b( i, j )
                  temp2 = zero
                  do 50, k = 1, i - 1
                     c( k, j ) = c( k, j ) + temp1    *a( k, i )
                     temp2     = temp2     + b( k, j )*a( k, i )
   50             continue
                  if( beta.eq.zero )then
                     c( i, j ) = temp1*a( i, i ) + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           temp1*a( i, i ) + alpha*temp2
                  end if
   60          continue
   70       continue
         else
            do 100, j = 1, n
               do 90, i = m, 1, -1
                  temp1 = alpha*b( i, j )
                  temp2 = zero
                  do 80, k = i + 1, m
                     c( k, j ) = c( k, j ) + temp1    *a( k, i )
                     temp2     = temp2     + b( k, j )*a( k, i )
   80             continue
                  if( beta.eq.zero )then
                     c( i, j ) = temp1*a( i, i ) + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           temp1*a( i, i ) + alpha*temp2
                  end if
   90          continue
  100       continue
         end if
      else
c
c        form  c := alpha*b*a + beta*c.
c
         do 170, j = 1, n
            temp1 = alpha*a( j, j )
            if( beta.eq.zero )then
               do 110, i = 1, m
                  c( i, j ) = temp1*b( i, j )
  110          continue
            else
               do 120, i = 1, m
                  c( i, j ) = beta*c( i, j ) + temp1*b( i, j )
  120          continue
            end if
            do 140, k = 1, j - 1
               if( upper )then
                  temp1 = alpha*a( k, j )
               else
                  temp1 = alpha*a( j, k )
               end if
               do 130, i = 1, m
                  c( i, j ) = c( i, j ) + temp1*b( i, k )
  130          continue
  140       continue
            do 160, k = j + 1, n
               if( upper )then
                  temp1 = alpha*a( j, k )
               else
                  temp1 = alpha*a( k, j )
               end if
               do 150, i = 1, m
                  c( i, j ) = c( i, j ) + temp1*b( i, k )
  150          continue
  160       continue
  170    continue
      end if
c
      return
c
c     end of csymm .
c
      end
c
c***********************************************************************
c
      subroutine chemm ( side, uplo, m, n, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
c     .. scalar arguments ..
      character*1        side, uplo
      integer            m, n, lda, ldb, ldc
      complex            alpha, beta
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  chemm  performs one of the matrix-matrix operations
c
c     c := alpha*a*b + beta*c,
c
c  or
c
c     c := alpha*b*a + beta*c,
c
c  where alpha and beta are scalars, a is an hermitian matrix and  b and
c  c are m by n matrices.
c
c  parameters
c  ==========
c
c  side   - character*1.
c           on entry,  side  specifies whether  the  hermitian matrix  a
c           appears on the  left or right  in the  operation as follows:
c
c              side = 'l' or 'l'   c := alpha*a*b + beta*c,
c
c              side = 'r' or 'r'   c := alpha*b*a + beta*c,
c
c           unchanged on exit.
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of  the  hermitian  matrix   a  is  to  be
c           referenced as follows:
c
c              uplo = 'u' or 'u'   only the upper triangular part of the
c                                  hermitian matrix is to be referenced.
c
c              uplo = 'l' or 'l'   only the lower triangular part of the
c                                  hermitian matrix is to be referenced.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry,  m  specifies the number of rows of the matrix  c.
c           m  must be at least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of the matrix c.
c           n  must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
c           m  when  side = 'l' or 'l'  and is n  otherwise.
c           before entry  with  side = 'l' or 'l',  the  m by m  part of
c           the array  a  must contain the  hermitian matrix,  such that
c           when  uplo = 'u' or 'u', the leading m by m upper triangular
c           part of the array  a  must contain the upper triangular part
c           of the  hermitian matrix and the  strictly  lower triangular
c           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
c           the leading  m by m  lower triangular part  of the  array  a
c           must  contain  the  lower triangular part  of the  hermitian
c           matrix and the  strictly upper triangular part of  a  is not
c           referenced.
c           before entry  with  side = 'r' or 'r',  the  n by n  part of
c           the array  a  must contain the  hermitian matrix,  such that
c           when  uplo = 'u' or 'u', the leading n by n upper triangular
c           part of the array  a  must contain the upper triangular part
c           of the  hermitian matrix and the  strictly  lower triangular
c           part of  a  is not referenced,  and when  uplo = 'l' or 'l',
c           the leading  n by n  lower triangular part  of the  array  a
c           must  contain  the  lower triangular part  of the  hermitian
c           matrix and the  strictly upper triangular part of  a  is not
c           referenced.
c           note that the imaginary parts  of the diagonal elements need
c           not be set, they are assumed to be zero.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the  calling (sub) program. when  side = 'l' or 'l'  then
c           lda must be at least  max( 1, m ), otherwise  lda must be at
c           least max( 1, n ).
c           unchanged on exit.
c
c  b      - complex          array of dimension ( ldb, n ).
c           before entry, the leading  m by n part of the array  b  must
c           contain the matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   ldb  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry,  beta  specifies the scalar  beta.  when  beta  is
c           supplied as zero then c need not be set on input.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
c           before entry, the leading  m by n  part of the array  c must
c           contain the matrix  c,  except when  beta  is zero, in which
c           case c need not be set on entry.
c           on exit, the array  c  is overwritten by the  m by n updated
c           matrix.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, real
c     .. local scalars ..
      logical            upper
      integer            i, info, j, k, nrowa
      complex            temp1, temp2
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     ..
c     .. executable statements ..
c
c     set nrowa as the number of rows of a.
c
      if( lsame( side, 'l' ) )then
         nrowa = m
      else
         nrowa = n
      end if
      upper = lsame( uplo, 'u' )
c
c     test the input parameters.
c
      info = 0
      if(      ( .not.lsame( side, 'l' ) ).and.
     $         ( .not.lsame( side, 'r' ) )      )then
         info = 1
      else if( ( .not.upper              ).and.
     $         ( .not.lsame( uplo, 'l' ) )      )then
         info = 2
      else if( m  .lt.0               )then
         info = 3
      else if( n  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldb.lt.max( 1, m     ) )then
         info = 9
      else if( ldc.lt.max( 1, m     ) )then
         info = 12
      end if
      if( info.ne.0 )then
         call xerbla( 'chemm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( m.eq.0 ).or.( n.eq.0 ).or.
     $    ( ( alpha.eq.zero ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( beta.eq.zero )then
            do 20, j = 1, n
               do 10, i = 1, m
                  c( i, j ) = zero
   10          continue
   20       continue
         else
            do 40, j = 1, n
               do 30, i = 1, m
                  c( i, j ) = beta*c( i, j )
   30          continue
   40       continue
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( side, 'l' ) )then
c
c        form  c := alpha*a*b + beta*c.
c
         if( upper )then
            do 70, j = 1, n
               do 60, i = 1, m
                  temp1 = alpha*b( i, j )
                  temp2 = zero
                  do 50, k = 1, i - 1
                     c( k, j ) = c( k, j ) + temp1*a( k, i )
                     temp2     = temp2     +
     $                           b( k, j )*conjg(  a( k, i ) )
   50             continue
                  if( beta.eq.zero )then
                     c( i, j ) = temp1*real( a( i, i ) ) +
     $                           alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j )         +
     $                           temp1*real( a( i, i ) ) +
     $                           alpha*temp2
                  end if
   60          continue
   70       continue
         else
            do 100, j = 1, n
               do 90, i = m, 1, -1
                  temp1 = alpha*b( i, j )
                  temp2 = zero
                  do 80, k = i + 1, m
                     c( k, j ) = c( k, j ) + temp1*a( k, i )
                     temp2     = temp2     +
     $                           b( k, j )*conjg(  a( k, i ) )
   80             continue
                  if( beta.eq.zero )then
                     c( i, j ) = temp1*real( a( i, i ) ) +
     $                           alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j )         +
     $                           temp1*real( a( i, i ) ) +
     $                           alpha*temp2
                  end if
   90          continue
  100       continue
         end if
      else
c
c        form  c := alpha*b*a + beta*c.
c
         do 170, j = 1, n
            temp1 = alpha*real( a( j, j ) )
            if( beta.eq.zero )then
               do 110, i = 1, m
                  c( i, j ) = temp1*b( i, j )
  110          continue
            else
               do 120, i = 1, m
                  c( i, j ) = beta*c( i, j ) + temp1*b( i, j )
  120          continue
            end if
            do 140, k = 1, j - 1
               if( upper )then
                  temp1 = alpha*a( k, j )
               else
                  temp1 = alpha*conjg( a( j, k ) )
               end if
               do 130, i = 1, m
                  c( i, j ) = c( i, j ) + temp1*b( i, k )
  130          continue
  140       continue
            do 160, k = j + 1, n
               if( upper )then
                  temp1 = alpha*conjg( a( j, k ) )
               else
                  temp1 = alpha*a( k, j )
               end if
               do 150, i = 1, m
                  c( i, j ) = c( i, j ) + temp1*b( i, k )
  150          continue
  160       continue
  170    continue
      end if
c
      return
c
c     end of chemm .
c
      end
c
c***********************************************************************
c
      subroutine csyrk ( uplo, trans, n, k, alpha, a, lda,
     $                   beta, c, ldc )
c     .. scalar arguments ..
      character*1        uplo, trans
      integer            n, k, lda, ldc
      complex            alpha, beta
c     .. array arguments ..
      complex            a( lda, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  csyrk  performs one of the symmetric rank k operations
c
c     c := alpha*a*a' + beta*c,
c
c  or
c
c     c := alpha*a'*a + beta*c,
c
c  where  alpha and beta  are scalars,  c is an  n by n symmetric matrix
c  and  a  is an  n by k  matrix in the first case and a  k by n  matrix
c  in the second case.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  c  is to be  referenced  as
c           follows:
c
c              uplo = 'u' or 'u'   only the  upper triangular part of  c
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the  lower triangular part of  c
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry,  trans  specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   c := alpha*a*a' + beta*c.
c
c              trans = 't' or 't'   c := alpha*a'*a + beta*c.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n specifies the order of the matrix c.  n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with  trans = 'n' or 'n',  k  specifies  the number
c           of  columns   of  the   matrix   a,   and  on   entry   with
c           trans = 't' or 't',  k  specifies  the number of rows of the
c           matrix a.  k must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by n  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  lda must be at least  max( 1, n ), otherwise  lda must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
c           before entry  with  uplo = 'u' or 'u',  the leading  n by n
c           upper triangular part of the array c must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of c is not referenced.  on exit, the
c           upper triangular part of the array  c is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry  with  uplo = 'l' or 'l',  the leading  n by n
c           lower triangular part of the array c must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of c is not referenced.  on exit, the
c           lower triangular part of the array  c is overwritten by the
c           lower triangular part of the updated matrix.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            upper
      integer            i, info, j, l, nrowa
      complex            temp
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      if( lsame( trans, 'n' ) )then
         nrowa = n
      else
         nrowa = k
      end if
      upper = lsame( uplo, 'u' )
c
      info = 0
      if(      ( .not.upper               ).and.
     $         ( .not.lsame( uplo , 'l' ) )      )then
         info = 1
      else if( ( .not.lsame( trans, 'n' ) ).and.
     $         ( .not.lsame( trans, 't' ) )      )then
         info = 2
      else if( n  .lt.0               )then
         info = 3
      else if( k  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldc.lt.max( 1, n     ) )then
         info = 10
      end if
      if( info.ne.0 )then
         call xerbla( 'csyrk ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( upper )then
            if( beta.eq.zero )then
               do 20, j = 1, n
                  do 10, i = 1, j
                     c( i, j ) = zero
   10             continue
   20          continue
            else
               do 40, j = 1, n
                  do 30, i = 1, j
                     c( i, j ) = beta*c( i, j )
   30             continue
   40          continue
            end if
         else
            if( beta.eq.zero )then
               do 60, j = 1, n
                  do 50, i = j, n
                     c( i, j ) = zero
   50             continue
   60          continue
            else
               do 80, j = 1, n
                  do 70, i = j, n
                     c( i, j ) = beta*c( i, j )
   70             continue
   80          continue
            end if
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( trans, 'n' ) )then
c
c        form  c := alpha*a*a' + beta*c.
c
         if( upper )then
            do 130, j = 1, n
               if( beta.eq.zero )then
                  do 90, i = 1, j
                     c( i, j ) = zero
   90             continue
               else if( beta.ne.one )then
                  do 100, i = 1, j
                     c( i, j ) = beta*c( i, j )
  100             continue
               end if
               do 120, l = 1, k
                  if( a( j, l ).ne.zero )then
                     temp = alpha*a( j, l )
                     do 110, i = 1, j
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  110                continue
                  end if
  120          continue
  130       continue
         else
            do 180, j = 1, n
               if( beta.eq.zero )then
                  do 140, i = j, n
                     c( i, j ) = zero
  140             continue
               else if( beta.ne.one )then
                  do 150, i = j, n
                     c( i, j ) = beta*c( i, j )
  150             continue
               end if
               do 170, l = 1, k
                  if( a( j, l ).ne.zero )then
                     temp      = alpha*a( j, l )
                     do 160, i = j, n
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  160                continue
                  end if
  170          continue
  180       continue
         end if
      else
c
c        form  c := alpha*a'*a + beta*c.
c
         if( upper )then
            do 210, j = 1, n
               do 200, i = 1, j
                  temp = zero
                  do 190, l = 1, k
                     temp = temp + a( l, i )*a( l, j )
  190             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  200          continue
  210       continue
         else
            do 240, j = 1, n
               do 230, i = j, n
                  temp = zero
                  do 220, l = 1, k
                     temp = temp + a( l, i )*a( l, j )
  220             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  230          continue
  240       continue
         end if
      end if
c
      return
c
c     end of csyrk .
c
      end
c
c***********************************************************************
c
      subroutine cherk ( uplo, trans, n, k, alpha, a, lda,
     $                   beta, c, ldc )
c     .. scalar arguments ..
      character*1        uplo, trans
      integer            n, k, lda, ldc
      real               alpha, beta
c     .. array arguments ..
      complex            a( lda, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  cherk  performs one of the hermitian rank k operations
c
c     c := alpha*a*conjg( a' ) + beta*c,
c
c  or
c
c     c := alpha*conjg( a' )*a + beta*c,
c
c  where  alpha and beta  are  real scalars,  c is an  n by n  hermitian
c  matrix and  a  is an  n by k  matrix in the  first case and a  k by n
c  matrix in the second case.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  c  is to be  referenced  as
c           follows:
c
c              uplo = 'u' or 'u'   only the  upper triangular part of  c
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the  lower triangular part of  c
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry,  trans  specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'   c := alpha*a*conjg( a' ) + beta*c.
c
c              trans = 'c' or 'c'   c := alpha*conjg( a' )*a + beta*c.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n specifies the order of the matrix c.  n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with  trans = 'n' or 'n',  k  specifies  the number
c           of  columns   of  the   matrix   a,   and  on   entry   with
c           trans = 'c' or 'c',  k  specifies  the number of rows of the
c           matrix a.  k must be at least zero.
c           unchanged on exit.
c
c  alpha  - real            .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by n  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  lda must be at least  max( 1, n ), otherwise  lda must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
c           before entry  with  uplo = 'u' or 'u',  the leading  n by n
c           upper triangular part of the array c must contain the upper
c           triangular part  of the  hermitian matrix  and the strictly
c           lower triangular part of c is not referenced.  on exit, the
c           upper triangular part of the array  c is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry  with  uplo = 'l' or 'l',  the leading  n by n
c           lower triangular part of the array c must contain the lower
c           triangular part  of the  hermitian matrix  and the strictly
c           upper triangular part of c is not referenced.  on exit, the
c           lower triangular part of the array  c is overwritten by the
c           lower triangular part of the updated matrix.
c           note that the imaginary parts of the diagonal elements need
c           not be set,  they are assumed to be zero,  and on exit they
c           are set to zero.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          cmplx, conjg, max, real
c     .. local scalars ..
      logical            upper
      integer            i, info, j, l, nrowa
      real               rtemp
      complex            temp
c     .. parameters ..
      real               one ,         zero
      parameter        ( one = 1.0e+0, zero = 0.0e+0 )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      if( lsame( trans, 'n' ) )then
         nrowa = n
      else
         nrowa = k
      end if
      upper = lsame( uplo, 'u' )
c
      info = 0
      if(      ( .not.upper               ).and.
     $         ( .not.lsame( uplo , 'l' ) )      )then
         info = 1
      else if( ( .not.lsame( trans, 'n' ) ).and.
     $         ( .not.lsame( trans, 'c' ) )      )then
         info = 2
      else if( n  .lt.0               )then
         info = 3
      else if( k  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldc.lt.max( 1, n     ) )then
         info = 10
      end if
      if( info.ne.0 )then
         call xerbla( 'cherk ', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( upper )then
            if( beta.eq.zero )then
               do 20, j = 1, n
                  do 10, i = 1, j
                     c( i, j ) = zero
   10             continue
   20          continue
            else
               do 40, j = 1, n
                  do 30, i = 1, j - 1
                     c( i, j ) = beta*c( i, j )
   30             continue
                  c( j, j ) = beta*real( c( j, j ) )
   40          continue
            end if
         else
            if( beta.eq.zero )then
               do 60, j = 1, n
                  do 50, i = j, n
                     c( i, j ) = zero
   50             continue
   60          continue
            else
               do 80, j = 1, n
                  c( j, j ) = beta*real( c( j, j ) )
                  do 70, i = j + 1, n
                     c( i, j ) = beta*c( i, j )
   70             continue
   80          continue
            end if
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( trans, 'n' ) )then
c
c        form  c := alpha*a*conjg( a' ) + beta*c.
c
         if( upper )then
            do 130, j = 1, n
               if( beta.eq.zero )then
                  do 90, i = 1, j
                     c( i, j ) = zero
   90             continue
               else if( beta.ne.one )then
                  do 100, i = 1, j - 1
                     c( i, j ) = beta*c( i, j )
  100             continue
                  c( j, j ) = beta*real( c( j, j ) )
               end if
               do 120, l = 1, k
                  if( a( j, l ).ne.cmplx( zero ) )then
                     temp = alpha*conjg( a( j, l ) )
                     do 110, i = 1, j - 1
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  110                continue
                     c( j, j ) = real( c( j, j )      ) +
     $                           real( temp*a( i, l ) )
                  end if
  120          continue
  130       continue
         else
            do 180, j = 1, n
               if( beta.eq.zero )then
                  do 140, i = j, n
                     c( i, j ) = zero
  140             continue
               else if( beta.ne.one )then
                  c( j, j ) = beta*real( c( j, j ) )
                  do 150, i = j + 1, n
                     c( i, j ) = beta*c( i, j )
  150             continue
               end if
               do 170, l = 1, k
                  if( a( j, l ).ne.cmplx( zero ) )then
                     temp      = alpha*conjg( a( j, l ) )
                     c( j, j ) = real( c( j, j )      )   +
     $                           real( temp*a( j, l ) )
                     do 160, i = j + 1, n
                        c( i, j ) = c( i, j ) + temp*a( i, l )
  160                continue
                  end if
  170          continue
  180       continue
         end if
      else
c
c        form  c := alpha*conjg( a' )*a + beta*c.
c
         if( upper )then
            do 220, j = 1, n
               do 200, i = 1, j - 1
                  temp = zero
                  do 190, l = 1, k
                     temp = temp + conjg( a( l, i ) )*a( l, j )
  190             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  200          continue
               rtemp = zero
               do 210, l = 1, k
                  rtemp = rtemp + conjg( a( l, j ) )*a( l, j )
  210          continue
               if( beta.eq.zero )then
                  c( j, j ) = alpha*rtemp
               else
                  c( j, j ) = alpha*rtemp + beta*real( c( j, j ) )
               end if
  220       continue
         else
            do 260, j = 1, n
               rtemp = zero
               do 230, l = 1, k
                  rtemp = rtemp + conjg( a( l, j ) )*a( l, j )
  230          continue
               if( beta.eq.zero )then
                  c( j, j ) = alpha*rtemp
               else
                  c( j, j ) = alpha*rtemp + beta*real( c( j, j ) )
               end if
               do 250, i = j + 1, n
                  temp = zero
                  do 240, l = 1, k
                     temp = temp + conjg( a( l, i ) )*a( l, j )
  240             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp
                  else
                     c( i, j ) = alpha*temp + beta*c( i, j )
                  end if
  250          continue
  260       continue
         end if
      end if
c
      return
c
c     end of cherk .
c
      end
c
c***********************************************************************
c
      subroutine csyr2k( uplo, trans, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
c     .. scalar arguments ..
      character*1        uplo, trans
      integer            n, k, lda, ldb, ldc
      complex            alpha, beta
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  csyr2k  performs one of the symmetric rank 2k operations
c
c     c := alpha*a*b' + alpha*b*a' + beta*c,
c
c  or
c
c     c := alpha*a'*b + alpha*b'*a + beta*c,
c
c  where  alpha and beta  are scalars,  c is an  n by n symmetric matrix
c  and  a and b  are  n by k  matrices  in the  first  case  and  k by n
c  matrices in the second case.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  c  is to be  referenced  as
c           follows:
c
c              uplo = 'u' or 'u'   only the  upper triangular part of  c
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the  lower triangular part of  c
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry,  trans  specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'    c := alpha*a*b' + alpha*b*a' +
c                                         beta*c.
c
c              trans = 't' or 't'    c := alpha*a'*b + alpha*b'*a +
c                                         beta*c.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n specifies the order of the matrix c.  n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with  trans = 'n' or 'n',  k  specifies  the number
c           of  columns  of the  matrices  a and b,  and on  entry  with
c           trans = 't' or 't',  k  specifies  the number of rows of the
c           matrices  a and b.  k must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by n  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  lda must be at least  max( 1, n ), otherwise  lda must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  b      - complex          array of dimension ( ldb, kb ), where kb is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  b  must contain the matrix  b,  otherwise
c           the leading  k by n  part of the array  b  must contain  the
c           matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  ldb must be at least  max( 1, n ), otherwise  ldb must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  beta   - complex         .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
c           before entry  with  uplo = 'u' or 'u',  the leading  n by n
c           upper triangular part of the array c must contain the upper
c           triangular part  of the  symmetric matrix  and the strictly
c           lower triangular part of c is not referenced.  on exit, the
c           upper triangular part of the array  c is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry  with  uplo = 'l' or 'l',  the leading  n by n
c           lower triangular part of the array c must contain the lower
c           triangular part  of the  symmetric matrix  and the strictly
c           upper triangular part of c is not referenced.  on exit, the
c           lower triangular part of the array  c is overwritten by the
c           lower triangular part of the updated matrix.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          max
c     .. local scalars ..
      logical            upper
      integer            i, info, j, l, nrowa
      complex            temp1, temp2
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      if( lsame( trans, 'n' ) )then
         nrowa = n
      else
         nrowa = k
      end if
      upper = lsame( uplo, 'u' )
c
      info = 0
      if(      ( .not.upper               ).and.
     $         ( .not.lsame( uplo , 'l' ) )      )then
         info = 1
      else if( ( .not.lsame( trans, 'n' ) ).and.
     $         ( .not.lsame( trans, 't' ) )      )then
         info = 2
      else if( n  .lt.0               )then
         info = 3
      else if( k  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldb.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldc.lt.max( 1, n     ) )then
         info = 12
      end if
      if( info.ne.0 )then
         call xerbla( 'csyr2k', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( upper )then
            if( beta.eq.zero )then
               do 20, j = 1, n
                  do 10, i = 1, j
                     c( i, j ) = zero
   10             continue
   20          continue
            else
               do 40, j = 1, n
                  do 30, i = 1, j
                     c( i, j ) = beta*c( i, j )
   30             continue
   40          continue
            end if
         else
            if( beta.eq.zero )then
               do 60, j = 1, n
                  do 50, i = j, n
                     c( i, j ) = zero
   50             continue
   60          continue
            else
               do 80, j = 1, n
                  do 70, i = j, n
                     c( i, j ) = beta*c( i, j )
   70             continue
   80          continue
            end if
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( trans, 'n' ) )then
c
c        form  c := alpha*a*b' + alpha*b*a' + c.
c
         if( upper )then
            do 130, j = 1, n
               if( beta.eq.zero )then
                  do 90, i = 1, j
                     c( i, j ) = zero
   90             continue
               else if( beta.ne.one )then
                  do 100, i = 1, j
                     c( i, j ) = beta*c( i, j )
  100             continue
               end if
               do 120, l = 1, k
                  if( ( a( j, l ).ne.zero ).or.
     $                ( b( j, l ).ne.zero )     )then
                     temp1 = alpha*b( j, l )
                     temp2 = alpha*a( j, l )
                     do 110, i = 1, j
                        c( i, j ) = c( i, j ) + a( i, l )*temp1 +
     $                                          b( i, l )*temp2
  110                continue
                  end if
  120          continue
  130       continue
         else
            do 180, j = 1, n
               if( beta.eq.zero )then
                  do 140, i = j, n
                     c( i, j ) = zero
  140             continue
               else if( beta.ne.one )then
                  do 150, i = j, n
                     c( i, j ) = beta*c( i, j )
  150             continue
               end if
               do 170, l = 1, k
                  if( ( a( j, l ).ne.zero ).or.
     $                ( b( j, l ).ne.zero )     )then
                     temp1 = alpha*b( j, l )
                     temp2 = alpha*a( j, l )
                     do 160, i = j, n
                        c( i, j ) = c( i, j ) + a( i, l )*temp1 +
     $                                          b( i, l )*temp2
  160                continue
                  end if
  170          continue
  180       continue
         end if
      else
c
c        form  c := alpha*a'*b + alpha*b'*a + c.
c
         if( upper )then
            do 210, j = 1, n
               do 200, i = 1, j
                  temp1 = zero
                  temp2 = zero
                  do 190, l = 1, k
                     temp1 = temp1 + a( l, i )*b( l, j )
                     temp2 = temp2 + b( l, i )*a( l, j )
  190             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp1 + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           alpha*temp1 + alpha*temp2
                  end if
  200          continue
  210       continue
         else
            do 240, j = 1, n
               do 230, i = j, n
                  temp1 = zero
                  temp2 = zero
                  do 220, l = 1, k
                     temp1 = temp1 + a( l, i )*b( l, j )
                     temp2 = temp2 + b( l, i )*a( l, j )
  220             continue
                  if( beta.eq.zero )then
                     c( i, j ) = alpha*temp1 + alpha*temp2
                  else
                     c( i, j ) = beta *c( i, j ) +
     $                           alpha*temp1 + alpha*temp2
                  end if
  230          continue
  240       continue
         end if
      end if
c
      return
c
c     end of csyr2k.
c
      end
c
c***********************************************************************
c
      subroutine cher2k( uplo, trans, n, k, alpha, a, lda, b, ldb,
     $                   beta, c, ldc )
c     .. scalar arguments ..
      character*1        uplo, trans
      integer            n, k, lda, ldb, ldc
      real               beta
      complex            alpha
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * ), c( ldc, * )
c     ..
c
c  purpose
c  =======
c
c  cher2k  performs one of the hermitian rank 2k operations
c
c     c := alpha*a*conjg( b' ) + conjg( alpha )*b*conjg( a' ) + beta*c,
c
c  or
c
c     c := alpha*conjg( a' )*b + conjg( alpha )*conjg( b' )*a + beta*c,
c
c  where  alpha and beta  are scalars with  beta  real,  c is an  n by n
c  hermitian matrix and  a and b  are  n by k matrices in the first case
c  and  k by n  matrices in the second case.
c
c  parameters
c  ==========
c
c  uplo   - character*1.
c           on  entry,   uplo  specifies  whether  the  upper  or  lower
c           triangular  part  of the  array  c  is to be  referenced  as
c           follows:
c
c              uplo = 'u' or 'u'   only the  upper triangular part of  c
c                                  is to be referenced.
c
c              uplo = 'l' or 'l'   only the  lower triangular part of  c
c                                  is to be referenced.
c
c           unchanged on exit.
c
c  trans  - character*1.
c           on entry,  trans  specifies the operation to be performed as
c           follows:
c
c              trans = 'n' or 'n'    c := alpha*a*conjg( b' )          +
c                                         conjg( alpha )*b*conjg( a' ) +
c                                         beta*c.
c
c              trans = 'c' or 'c'    c := alpha*conjg( a' )*b          +
c                                         conjg( alpha )*conjg( b' )*a +
c                                         beta*c.
c
c           unchanged on exit.
c
c  n      - integer.
c           on entry,  n specifies the order of the matrix c.  n must be
c           at least zero.
c           unchanged on exit.
c
c  k      - integer.
c           on entry with  trans = 'n' or 'n',  k  specifies  the number
c           of  columns  of the  matrices  a and b,  and on  entry  with
c           trans = 'c' or 'c',  k  specifies  the number of rows of the
c           matrices  a and b.  k must be at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry, alpha specifies the scalar alpha.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, ka ), where ka is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  a  must contain the matrix  a,  otherwise
c           the leading  k by n  part of the array  a  must contain  the
c           matrix a.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  lda must be at least  max( 1, n ), otherwise  lda must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  b      - complex          array of dimension ( ldb, kb ), where kb is
c           k  when  trans = 'n' or 'n',  and is  n  otherwise.
c           before entry with  trans = 'n' or 'n',  the  leading  n by k
c           part of the array  b  must contain the matrix  b,  otherwise
c           the leading  k by n  part of the array  b  must contain  the
c           matrix b.
c           unchanged on exit.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   when  trans = 'n' or 'n'
c           then  ldb must be at least  max( 1, n ), otherwise  ldb must
c           be at least  max( 1, k ).
c           unchanged on exit.
c
c  beta   - real            .
c           on entry, beta specifies the scalar beta.
c           unchanged on exit.
c
c  c      - complex          array of dimension ( ldc, n ).
c           before entry  with  uplo = 'u' or 'u',  the leading  n by n
c           upper triangular part of the array c must contain the upper
c           triangular part  of the  hermitian matrix  and the strictly
c           lower triangular part of c is not referenced.  on exit, the
c           upper triangular part of the array  c is overwritten by the
c           upper triangular part of the updated matrix.
c           before entry  with  uplo = 'l' or 'l',  the leading  n by n
c           lower triangular part of the array c must contain the lower
c           triangular part  of the  hermitian matrix  and the strictly
c           upper triangular part of c is not referenced.  on exit, the
c           lower triangular part of the array  c is overwritten by the
c           lower triangular part of the updated matrix.
c           note that the imaginary parts of the diagonal elements need
c           not be set,  they are assumed to be zero,  and on exit they
c           are set to zero.
c
c  ldc    - integer.
c           on entry, ldc specifies the first dimension of c as declared
c           in  the  calling  (sub)  program.   ldc  must  be  at  least
c           max( 1, n ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max, real
c     .. local scalars ..
      logical            upper
      integer            i, info, j, l, nrowa
      complex            temp1, temp2
c     .. parameters ..
      real               one
      parameter        ( one  = 1.0e+0 )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      if( lsame( trans, 'n' ) )then
         nrowa = n
      else
         nrowa = k
      end if
      upper = lsame( uplo, 'u' )
c
      info = 0
      if(      ( .not.upper               ).and.
     $         ( .not.lsame( uplo , 'l' ) )      )then
         info = 1
      else if( ( .not.lsame( trans, 'n' ) ).and.
     $         ( .not.lsame( trans, 'c' ) )      )then
         info = 2
      else if( n  .lt.0               )then
         info = 3
      else if( k  .lt.0               )then
         info = 4
      else if( lda.lt.max( 1, nrowa ) )then
         info = 7
      else if( ldb.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldc.lt.max( 1, n     ) )then
         info = 12
      end if
      if( info.ne.0 )then
         call xerbla( 'cher2k', info )
         return
      end if
c
c     quick return if possible.
c
      if( ( n.eq.0 ).or.
     $    ( ( ( alpha.eq.zero ).or.( k.eq.0 ) ).and.( beta.eq.one ) ) )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         if( upper )then
            if( beta.eq.real( zero ) )then
               do 20, j = 1, n
                  do 10, i = 1, j
                     c( i, j ) = zero
   10             continue
   20          continue
            else
               do 40, j = 1, n
                  do 30, i = 1, j - 1
                     c( i, j ) = beta*c( i, j )
   30             continue
                  c( j, j ) = beta*real( c( j, j ) )
   40          continue
            end if
         else
            if( beta.eq.real( zero ) )then
               do 60, j = 1, n
                  do 50, i = j, n
                     c( i, j ) = zero
   50             continue
   60          continue
            else
               do 80, j = 1, n
                  c( j, j ) = beta*real( c( j, j ) )
                  do 70, i = j + 1, n
                     c( i, j ) = beta*c( i, j )
   70             continue
   80          continue
            end if
         end if
         return
      end if
c
c     start the operations.
c
      if( lsame( trans, 'n' ) )then
c
c        form  c := alpha*a*conjg( b' ) + conjg( alpha )*b*conjg( a' ) +
c                   c.
c
         if( upper )then
            do 130, j = 1, n
               if( beta.eq.real( zero ) )then
                  do 90, i = 1, j
                     c( i, j ) = zero
   90             continue
               else if( beta.ne.one )then
                  do 100, i = 1, j - 1
                     c( i, j ) = beta*c( i, j )
  100             continue
                  c( j, j ) = beta*real( c( j, j ) )
               end if
               do 120, l = 1, k
                  if( ( a( j, l ).ne.zero ).or.
     $                ( b( j, l ).ne.zero )     )then
                     temp1 = alpha*conjg( b( j, l ) )
                     temp2 = conjg( alpha*a( j, l ) )
                     do 110, i = 1, j - 1
                        c( i, j ) = c( i, j ) + a( i, l )*temp1 +
     $                                          b( i, l )*temp2
  110                continue
                     c( j, j ) = real( c( j, j ) )         +
     $                           real( a( j, l )*temp1 +
     $                                 b( j, l )*temp2   )
                  end if
  120          continue
  130       continue
         else
            do 180, j = 1, n
               if( beta.eq.real( zero ) )then
                  do 140, i = j, n
                     c( i, j ) = zero
  140             continue
               else if( beta.ne.one )then
                  do 150, i = j + 1, n
                     c( i, j ) = beta*c( i, j )
  150             continue
                  c( j, j ) = beta*real( c( j, j ) )
               end if
               do 170, l = 1, k
                  if( ( a( j, l ).ne.zero ).or.
     $                ( b( j, l ).ne.zero )     )then
                     temp1 = alpha*conjg( b( j, l ) )
                     temp2 = conjg( alpha*a( j, l ) )
                     do 160, i = j + 1, n
                        c( i, j ) = c( i, j ) + a( i, l )*temp1 +
     $                                          b( i, l )*temp2
  160                continue
                     c( j, j ) = real( c( j, j ) )         +
     $                           real( a( j, l )*temp1 +
     $                                 b( j, l )*temp2   )
                  end if
  170          continue
  180       continue
         end if
      else
c
c        form  c := alpha*conjg( a' )*b + conjg( alpha )*conjg( b' )*a +
c                   c.
c
         if( upper )then
            do 210, j = 1, n
               do 200, i = 1, j
                  temp1 = zero
                  temp2 = zero
                  do 190, l = 1, k
                     temp1 = temp1 + conjg( a( l, i ) )*b( l, j )
                     temp2 = temp2 + conjg( b( l, i ) )*a( l, j )
  190             continue
                  if( i.eq.j )then
                     if( beta.eq.real( zero ) )then
                        c( j, j ) = real(        alpha  *temp1 +
     $                                    conjg( alpha )*temp2   )
                     else
                        c( j, j ) = beta*real( c( j, j ) )         +
     $                              real(        alpha  *temp1 +
     $                                    conjg( alpha )*temp2   )
                     end if
                  else
                     if( beta.eq.real( zero ) )then
                        c( i, j ) = alpha*temp1 + conjg( alpha )*temp2
                     else
                        c( i, j ) = beta *c( i, j ) +
     $                              alpha*temp1 + conjg( alpha )*temp2
                     end if
                  end if
  200          continue
  210       continue
         else
            do 240, j = 1, n
               do 230, i = j, n
                  temp1 = zero
                  temp2 = zero
                  do 220, l = 1, k
                     temp1 = temp1 + conjg( a( l, i ) )*b( l, j )
                     temp2 = temp2 + conjg( b( l, i ) )*a( l, j )
  220             continue
                  if( i.eq.j )then
                     if( beta.eq.real( zero ) )then
                        c( j, j ) = real(        alpha  *temp1 +
     $                                    conjg( alpha )*temp2   )
                     else
                        c( j, j ) = beta*real( c( j, j ) )         +
     $                              real(        alpha  *temp1 +
     $                                    conjg( alpha )*temp2   )
                     end if
                  else
                     if( beta.eq.real( zero ) )then
                        c( i, j ) = alpha*temp1 + conjg( alpha )*temp2
                     else
                        c( i, j ) = beta *c( i, j ) +
     $                              alpha*temp1 + conjg( alpha )*temp2
                     end if
                  end if
  230          continue
  240       continue
         end if
      end if
c
      return
c
c     end of cher2k.
c
      end
c
c***********************************************************************
c
      subroutine ctrmm ( side, uplo, transa, diag, m, n, alpha, a, lda,
     $                   b, ldb )
c     .. scalar arguments ..
      character*1        side, uplo, transa, diag
      integer            m, n, lda, ldb
      complex            alpha
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * )
c     ..
c
c  purpose
c  =======
c
c  ctrmm  performs one of the matrix-matrix operations
c
c     b := alpha*op( a )*b,   or   b := alpha*b*op( a )
c
c  where  alpha  is a scalar,  b  is an m by n matrix,  a  is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( a )  is one  of
c
c     op( a ) = a   or   op( a ) = a'   or   op( a ) = conjg( a' ).
c
c  parameters
c  ==========
c
c  side   - character*1.
c           on entry,  side specifies whether  op( a ) multiplies b from
c           the left or right as follows:
c
c              side = 'l' or 'l'   b := alpha*op( a )*b.
c
c              side = 'r' or 'r'   b := alpha*b*op( a ).
c
c           unchanged on exit.
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix a is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  transa - character*1.
c           on entry, transa specifies the form of op( a ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'n' or 'n'   op( a ) = a.
c
c              transa = 't' or 't'   op( a ) = a'.
c
c              transa = 'c' or 'c'   op( a ) = conjg( a' ).
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit triangular
c           as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of b. m must be at
c           least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of b.  n must be
c           at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry,  alpha specifies the scalar  alpha. when  alpha is
c           zero then  a is not referenced and  b need not be set before
c           entry.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, k ), where k is m
c           when  side = 'l' or 'l'  and is  n  when  side = 'r' or 'r'.
c           before entry  with  uplo = 'u' or 'u',  the  leading  k by k
c           upper triangular part of the array  a must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           a is not referenced.
c           before entry  with  uplo = 'l' or 'l',  the  leading  k by k
c           lower triangular part of the array  a must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u',  the diagonal elements of
c           a  are not referenced either,  but are assumed to be  unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program.  when  side = 'l' or 'l'  then
c           lda  must be at least  max( 1, m ),  when  side = 'r' or 'r'
c           then lda must be at least max( 1, n ).
c           unchanged on exit.
c
c  b      - complex          array of dimension ( ldb, n ).
c           before entry,  the leading  m by n part of the array  b must
c           contain the matrix  b,  and  on exit  is overwritten  by the
c           transformed matrix.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   ldb  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c     .. local scalars ..
      logical            lside, noconj, nounit, upper
      integer            i, info, j, k, nrowa
      complex            temp
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      lside  = lsame( side  , 'l' )
      if( lside )then
         nrowa = m
      else
         nrowa = n
      end if
      noconj = lsame( transa, 't' )
      nounit = lsame( diag  , 'n' )
      upper  = lsame( uplo  , 'u' )
c
      info   = 0
      if(      ( .not.lside                ).and.
     $         ( .not.lsame( side  , 'r' ) )      )then
         info = 1
      else if( ( .not.upper                ).and.
     $         ( .not.lsame( uplo  , 'l' ) )      )then
         info = 2
      else if( ( .not.lsame( transa, 'n' ) ).and.
     $         ( .not.lsame( transa, 't' ) ).and.
     $         ( .not.lsame( transa, 'c' ) )      )then
         info = 3
      else if( ( .not.lsame( diag  , 'u' ) ).and.
     $         ( .not.lsame( diag  , 'n' ) )      )then
         info = 4
      else if( m  .lt.0               )then
         info = 5
      else if( n  .lt.0               )then
         info = 6
      else if( lda.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldb.lt.max( 1, m     ) )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'ctrmm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         do 20, j = 1, n
            do 10, i = 1, m
               b( i, j ) = zero
   10       continue
   20    continue
         return
      end if
c
c     start the operations.
c
      if( lside )then
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*a*b.
c
            if( upper )then
               do 50, j = 1, n
                  do 40, k = 1, m
                     if( b( k, j ).ne.zero )then
                        temp = alpha*b( k, j )
                        do 30, i = 1, k - 1
                           b( i, j ) = b( i, j ) + temp*a( i, k )
   30                   continue
                        if( nounit )
     $                     temp = temp*a( k, k )
                        b( k, j ) = temp
                     end if
   40             continue
   50          continue
            else
               do 80, j = 1, n
                  do 70 k = m, 1, -1
                     if( b( k, j ).ne.zero )then
                        temp      = alpha*b( k, j )
                        b( k, j ) = temp
                        if( nounit )
     $                     b( k, j ) = b( k, j )*a( k, k )
                        do 60, i = k + 1, m
                           b( i, j ) = b( i, j ) + temp*a( i, k )
   60                   continue
                     end if
   70             continue
   80          continue
            end if
         else
c
c           form  b := alpha*b*a'   or   b := alpha*b*conjg( a' ).
c
            if( upper )then
               do 120, j = 1, n
                  do 110, i = m, 1, -1
                     temp = b( i, j )
                     if( noconj )then
                        if( nounit )
     $                     temp = temp*a( i, i )
                        do 90, k = 1, i - 1
                           temp = temp + a( k, i )*b( k, j )
   90                   continue
                     else
                        if( nounit )
     $                     temp = temp*conjg( a( i, i ) )
                        do 100, k = 1, i - 1
                           temp = temp + conjg( a( k, i ) )*b( k, j )
  100                   continue
                     end if
                     b( i, j ) = alpha*temp
  110             continue
  120          continue
            else
               do 160, j = 1, n
                  do 150, i = 1, m
                     temp = b( i, j )
                     if( noconj )then
                        if( nounit )
     $                     temp = temp*a( i, i )
                        do 130, k = i + 1, m
                           temp = temp + a( k, i )*b( k, j )
  130                   continue
                     else
                        if( nounit )
     $                     temp = temp*conjg( a( i, i ) )
                        do 140, k = i + 1, m
                           temp = temp + conjg( a( k, i ) )*b( k, j )
  140                   continue
                     end if
                     b( i, j ) = alpha*temp
  150             continue
  160          continue
            end if
         end if
      else
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*b*a.
c
            if( upper )then
               do 200, j = n, 1, -1
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 170, i = 1, m
                     b( i, j ) = temp*b( i, j )
  170             continue
                  do 190, k = 1, j - 1
                     if( a( k, j ).ne.zero )then
                        temp = alpha*a( k, j )
                        do 180, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  180                   continue
                     end if
  190             continue
  200          continue
            else
               do 240, j = 1, n
                  temp = alpha
                  if( nounit )
     $               temp = temp*a( j, j )
                  do 210, i = 1, m
                     b( i, j ) = temp*b( i, j )
  210             continue
                  do 230, k = j + 1, n
                     if( a( k, j ).ne.zero )then
                        temp = alpha*a( k, j )
                        do 220, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  220                   continue
                     end if
  230             continue
  240          continue
            end if
         else
c
c           form  b := alpha*b*a'   or   b := alpha*b*conjg( a' ).
c
            if( upper )then
               do 280, k = 1, n
                  do 260, j = 1, k - 1
                     if( a( j, k ).ne.zero )then
                        if( noconj )then
                           temp = alpha*a( j, k )
                        else
                           temp = alpha*conjg( a( j, k ) )
                        end if
                        do 250, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  250                   continue
                     end if
  260             continue
                  temp = alpha
                  if( nounit )then
                     if( noconj )then
                        temp = temp*a( k, k )
                     else
                        temp = temp*conjg( a( k, k ) )
                     end if
                  end if
                  if( temp.ne.one )then
                     do 270, i = 1, m
                        b( i, k ) = temp*b( i, k )
  270                continue
                  end if
  280          continue
            else
               do 320, k = n, 1, -1
                  do 300, j = k + 1, n
                     if( a( j, k ).ne.zero )then
                        if( noconj )then
                           temp = alpha*a( j, k )
                        else
                           temp = alpha*conjg( a( j, k ) )
                        end if
                        do 290, i = 1, m
                           b( i, j ) = b( i, j ) + temp*b( i, k )
  290                   continue
                     end if
  300             continue
                  temp = alpha
                  if( nounit )then
                     if( noconj )then
                        temp = temp*a( k, k )
                     else
                        temp = temp*conjg( a( k, k ) )
                     end if
                  end if
                  if( temp.ne.one )then
                     do 310, i = 1, m
                        b( i, k ) = temp*b( i, k )
  310                continue
                  end if
  320          continue
            end if
         end if
      end if
c
      return
c
c     end of ctrmm .
c
      end
c
c***********************************************************************
c
      subroutine ctrsm ( side, uplo, transa, diag, m, n, alpha, a, lda,
     $                   b, ldb )
c     .. scalar arguments ..
      character*1        side, uplo, transa, diag
      integer            m, n, lda, ldb
      complex            alpha
c     .. array arguments ..
      complex            a( lda, * ), b( ldb, * )
c     ..
c
c  purpose
c  =======
c
c  ctrsm  solves one of the matrix equations
c
c     op( a )*x = alpha*b,   or   x*op( a ) = alpha*b,
c
c  where alpha is a scalar, x and b are m by n matrices, a is a unit, or
c  non-unit,  upper or lower triangular matrix  and  op( a )  is one  of
c
c     op( a ) = a   or   op( a ) = a'   or   op( a ) = conjg( a' ).
c
c  the matrix x is overwritten on b.
c
c  parameters
c  ==========
c
c  side   - character*1.
c           on entry, side specifies whether op( a ) appears on the left
c           or right of x as follows:
c
c              side = 'l' or 'l'   op( a )*x = alpha*b.
c
c              side = 'r' or 'r'   x*op( a ) = alpha*b.
c
c           unchanged on exit.
c
c  uplo   - character*1.
c           on entry, uplo specifies whether the matrix a is an upper or
c           lower triangular matrix as follows:
c
c              uplo = 'u' or 'u'   a is an upper triangular matrix.
c
c              uplo = 'l' or 'l'   a is a lower triangular matrix.
c
c           unchanged on exit.
c
c  transa - character*1.
c           on entry, transa specifies the form of op( a ) to be used in
c           the matrix multiplication as follows:
c
c              transa = 'n' or 'n'   op( a ) = a.
c
c              transa = 't' or 't'   op( a ) = a'.
c
c              transa = 'c' or 'c'   op( a ) = conjg( a' ).
c
c           unchanged on exit.
c
c  diag   - character*1.
c           on entry, diag specifies whether or not a is unit triangular
c           as follows:
c
c              diag = 'u' or 'u'   a is assumed to be unit triangular.
c
c              diag = 'n' or 'n'   a is not assumed to be unit
c                                  triangular.
c
c           unchanged on exit.
c
c  m      - integer.
c           on entry, m specifies the number of rows of b. m must be at
c           least zero.
c           unchanged on exit.
c
c  n      - integer.
c           on entry, n specifies the number of columns of b.  n must be
c           at least zero.
c           unchanged on exit.
c
c  alpha  - complex         .
c           on entry,  alpha specifies the scalar  alpha. when  alpha is
c           zero then  a is not referenced and  b need not be set before
c           entry.
c           unchanged on exit.
c
c  a      - complex          array of dimension ( lda, k ), where k is m
c           when  side = 'l' or 'l'  and is  n  when  side = 'r' or 'r'.
c           before entry  with  uplo = 'u' or 'u',  the  leading  k by k
c           upper triangular part of the array  a must contain the upper
c           triangular matrix  and the strictly lower triangular part of
c           a is not referenced.
c           before entry  with  uplo = 'l' or 'l',  the  leading  k by k
c           lower triangular part of the array  a must contain the lower
c           triangular matrix  and the strictly upper triangular part of
c           a is not referenced.
c           note that when  diag = 'u' or 'u',  the diagonal elements of
c           a  are not referenced either,  but are assumed to be  unity.
c           unchanged on exit.
c
c  lda    - integer.
c           on entry, lda specifies the first dimension of a as declared
c           in the calling (sub) program.  when  side = 'l' or 'l'  then
c           lda  must be at least  max( 1, m ),  when  side = 'r' or 'r'
c           then lda must be at least max( 1, n ).
c           unchanged on exit.
c
c  b      - complex          array of dimension ( ldb, n ).
c           before entry,  the leading  m by n part of the array  b must
c           contain  the  right-hand  side  matrix  b,  and  on exit  is
c           overwritten by the solution matrix  x.
c
c  ldb    - integer.
c           on entry, ldb specifies the first dimension of b as declared
c           in  the  calling  (sub)  program.   ldb  must  be  at  least
c           max( 1, m ).
c           unchanged on exit.
c
c
c  level 3 blas routine.
c
c  -- written on 8-february-1989.
c     jack dongarra, argonne national laboratory.
c     iain duff, aere harwell.
c     jeremy du croz, numerical algorithms group ltd.
c     sven hammarling, numerical algorithms group ltd.
c
c
c     .. external functions ..
      logical            lsame
      external           lsame
c     .. external subroutines ..
      external           xerbla
c     .. intrinsic functions ..
      intrinsic          conjg, max
c     .. local scalars ..
      logical            lside, noconj, nounit, upper
      integer            i, info, j, k, nrowa
      complex            temp
c     .. parameters ..
      complex            one
      parameter        ( one  = ( 1.0e+0, 0.0e+0 ) )
      complex            zero
      parameter        ( zero = ( 0.0e+0, 0.0e+0 ) )
c     ..
c     .. executable statements ..
c
c     test the input parameters.
c
      lside  = lsame( side  , 'l' )
      if( lside )then
         nrowa = m
      else
         nrowa = n
      end if
      noconj = lsame( transa, 't' )
      nounit = lsame( diag  , 'n' )
      upper  = lsame( uplo  , 'u' )
c
      info   = 0
      if(      ( .not.lside                ).and.
     $         ( .not.lsame( side  , 'r' ) )      )then
         info = 1
      else if( ( .not.upper                ).and.
     $         ( .not.lsame( uplo  , 'l' ) )      )then
         info = 2
      else if( ( .not.lsame( transa, 'n' ) ).and.
     $         ( .not.lsame( transa, 't' ) ).and.
     $         ( .not.lsame( transa, 'c' ) )      )then
         info = 3
      else if( ( .not.lsame( diag  , 'u' ) ).and.
     $         ( .not.lsame( diag  , 'n' ) )      )then
         info = 4
      else if( m  .lt.0               )then
         info = 5
      else if( n  .lt.0               )then
         info = 6
      else if( lda.lt.max( 1, nrowa ) )then
         info = 9
      else if( ldb.lt.max( 1, m     ) )then
         info = 11
      end if
      if( info.ne.0 )then
         call xerbla( 'ctrsm ', info )
         return
      end if
c
c     quick return if possible.
c
      if( n.eq.0 )
     $   return
c
c     and when  alpha.eq.zero.
c
      if( alpha.eq.zero )then
         do 20, j = 1, n
            do 10, i = 1, m
               b( i, j ) = zero
   10       continue
   20    continue
         return
      end if
c
c     start the operations.
c
      if( lside )then
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*inv( a )*b.
c
            if( upper )then
               do 60, j = 1, n
                  if( alpha.ne.one )then
                     do 30, i = 1, m
                        b( i, j ) = alpha*b( i, j )
   30                continue
                  end if
                  do 50, k = m, 1, -1
                     if( b( k, j ).ne.zero )then
                        if( nounit )
     $                     b( k, j ) = b( k, j )/a( k, k )
                        do 40, i = 1, k - 1
                           b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
   40                   continue
                     end if
   50             continue
   60          continue
            else
               do 100, j = 1, n
                  if( alpha.ne.one )then
                     do 70, i = 1, m
                        b( i, j ) = alpha*b( i, j )
   70                continue
                  end if
                  do 90 k = 1, m
                     if( b( k, j ).ne.zero )then
                        if( nounit )
     $                     b( k, j ) = b( k, j )/a( k, k )
                        do 80, i = k + 1, m
                           b( i, j ) = b( i, j ) - b( k, j )*a( i, k )
   80                   continue
                     end if
   90             continue
  100          continue
            end if
         else
c
c           form  b := alpha*inv( a' )*b
c           or    b := alpha*inv( conjg( a' ) )*b.
c
            if( upper )then
               do 140, j = 1, n
                  do 130, i = 1, m
                     temp = alpha*b( i, j )
                     if( noconj )then
                        do 110, k = 1, i - 1
                           temp = temp - a( k, i )*b( k, j )
  110                   continue
                        if( nounit )
     $                     temp = temp/a( i, i )
                     else
                        do 120, k = 1, i - 1
                           temp = temp - conjg( a( k, i ) )*b( k, j )
  120                   continue
                        if( nounit )
     $                     temp = temp/conjg( a( i, i ) )
                     end if
                     b( i, j ) = temp
  130             continue
  140          continue
            else
               do 180, j = 1, n
                  do 170, i = m, 1, -1
                     temp = alpha*b( i, j )
                     if( noconj )then
                        do 150, k = i + 1, m
                           temp = temp - a( k, i )*b( k, j )
  150                   continue
                        if( nounit )
     $                     temp = temp/a( i, i )
                     else
                        do 160, k = i + 1, m
                           temp = temp - conjg( a( k, i ) )*b( k, j )
  160                   continue
                        if( nounit )
     $                     temp = temp/conjg( a( i, i ) )
                     end if
                     b( i, j ) = temp
  170             continue
  180          continue
            end if
         end if
      else
         if( lsame( transa, 'n' ) )then
c
c           form  b := alpha*b*inv( a ).
c
            if( upper )then
               do 230, j = 1, n
                  if( alpha.ne.one )then
                     do 190, i = 1, m
                        b( i, j ) = alpha*b( i, j )
  190                continue
                  end if
                  do 210, k = 1, j - 1
                     if( a( k, j ).ne.zero )then
                        do 200, i = 1, m
                           b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  200                   continue
                     end if
  210             continue
                  if( nounit )then
                     temp = one/a( j, j )
                     do 220, i = 1, m
                        b( i, j ) = temp*b( i, j )
  220                continue
                  end if
  230          continue
            else
               do 280, j = n, 1, -1
                  if( alpha.ne.one )then
                     do 240, i = 1, m
                        b( i, j ) = alpha*b( i, j )
  240                continue
                  end if
                  do 260, k = j + 1, n
                     if( a( k, j ).ne.zero )then
                        do 250, i = 1, m
                           b( i, j ) = b( i, j ) - a( k, j )*b( i, k )
  250                   continue
                     end if
  260             continue
                  if( nounit )then
                     temp = one/a( j, j )
                     do 270, i = 1, m
                       b( i, j ) = temp*b( i, j )
  270                continue
                  end if
  280          continue
            end if
         else
c
c           form  b := alpha*b*inv( a' )
c           or    b := alpha*b*inv( conjg( a' ) ).
c
            if( upper )then
               do 330, k = n, 1, -1
                  if( nounit )then
                     if( noconj )then
                        temp = one/a( k, k )
                     else
                        temp = one/conjg( a( k, k ) )
                     end if
                     do 290, i = 1, m
                        b( i, k ) = temp*b( i, k )
  290                continue
                  end if
                  do 310, j = 1, k - 1
                     if( a( j, k ).ne.zero )then
                        if( noconj )then
                           temp = a( j, k )
                        else
                           temp = conjg( a( j, k ) )
                        end if
                        do 300, i = 1, m
                           b( i, j ) = b( i, j ) - temp*b( i, k )
  300                   continue
                     end if
  310             continue
                  if( alpha.ne.one )then
                     do 320, i = 1, m
                        b( i, k ) = alpha*b( i, k )
  320                continue
                  end if
  330          continue
            else
               do 380, k = 1, n
                  if( nounit )then
                     if( noconj )then
                        temp = one/a( k, k )
                     else
                        temp = one/conjg( a( k, k ) )
                     end if
                     do 340, i = 1, m
                        b( i, k ) = temp*b( i, k )
  340                continue
                  end if
                  do 360, j = k + 1, n
                     if( a( j, k ).ne.zero )then
                        if( noconj )then
                           temp = a( j, k )
                        else
                           temp = conjg( a( j, k ) )
                        end if
                        do 350, i = 1, m
                           b( i, j ) = b( i, j ) - temp*b( i, k )
  350                   continue
                     end if
  360             continue
                  if( alpha.ne.one )then
                     do 370, i = 1, m
                        b( i, k ) = alpha*b( i, k )
  370                continue
                  end if
  380          continue
            end if
         end if
      end if
c
      return
c
c     end of ctrsm .
c
      end
