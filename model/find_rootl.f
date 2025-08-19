      subroutine find_rootl( f, y, xlow, xhigh, delta, xnew, status )
c
c   find solution to 
c
c     log( f(x) ) - log( y ) = 0 
c
c   for x between x1 and x2 to within log( delta ) using secant method
c
c   if found, return:
c     status = 0
c     xnew   = root 
c   else if maximum number of iterations reached before convergence:
c     status = 1
c     xnew   = last estimate of root
c   else (no convergence):
c     status = -1
c     xnew   = last estimate of root
c
c   A. Ackerman Feb-2000
c
      implicit none


c   declare externals
    
      double precision f
      external f


c   declare local storage

      double precision x1, x2, f1, f2, slope, fnew
      double precision y, xlow, xhigh, delta, xnew
      integer status, iter


c   define maximum number of iterations
      integer MAX_ITER
      parameter( MAX_ITER = 50 )


c   define convergence criteria for independent variable
      double precision EPSILON
      parameter( EPSILON = 1d-7 )


c   copy input range into local variables
      x1 = xlow
      x2 = xhigh


c   abort if root not bracketed by initial guesses 

      f1 = log(f(x1)) - log(y)
      f2 = log(f(x2)) - log(y)

c   default return values
      xnew = x1 
      fnew = f1 

      if( f1*f2 .gt. 0 )then 
        status = -1
        return
      endif

c   take logarithm of precision
      delta = log( delta )

c   iterate until root is found or maximum number of iterations are taken

      iter = 0
      do while( ( abs(fnew) .gt. delta ) .and. ( iter .lt. MAX_ITER )
     $  .and. ( x2/x1-1. .gt. EPSILON ) )

c   estimate root from secant between endpoints

        slope = ( f2 - f1 ) / ( x2 - x1 )
        xnew = x1 - f1/slope
        fnew = log(f(xnew)) - log(y)

c   pick new x1 and x2 such that f(x) crosses zero between them
      
        if( fnew*f1 .le. 0 .and. x2/xnew-1. .gt. EPSILON )then
          x2 = xnew
          f2 = fnew
        else 
          x1 = xnew
          f1 = fnew
        endif

        iter = iter + 1

      enddo


c   set status flag accordingly

      if( iter .eq. MAX_ITER .and. x2/x1-1. .gt. EPSILON )then 
        status = 1
      else
        status = 0
      endif

      return
      end
