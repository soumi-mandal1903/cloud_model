      subroutine dblank(s, ns)
c
c
c  @(#) dblank.f  McKie  Jun-1988
c  This routine finds index of last nonblank char in a string.
c
c  Argument list input:
c    s = string to be examined
c
c  Argument list output:
c    ns = # chars in s, up to & including last non-blank
c
c
c  Declare subprogram arg(s)
c
      character*(*) s
      integer ns
c
c
c  Find last non-blank char in string (beyond 1st char), return to caller
c
      ns = len(s)
 2100 if( ( s(ns:ns) .eq. ' ' ) .and. ( ns .gt. 1 ) )then
       ns = ns - 1
       goto 2100
      endif
      return
      end
