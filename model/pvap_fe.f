      double precision function pvap_fe( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over MgSiO3
c     (Lunine et al., 1988, ApJ, 338, 314, citing Barshay and Lewis 1976)
c
c     input temperature in K
c
c     A. Ackerman Feb-2000 modifications for iron by M. Kress March-2000
c
      implicit none

      double precision t

      goto 9999   ! use new Fe expression

      if( t .gt. 1800. )then
            pvap_fe = dexp(-37120./t + 9.86)
      else 
            pvap_fe = dexp(-47664./t + 15.71)
      endif
   
      pvap_fe = pvap_fe*1e6    ! convert from bars to dyne/cm^2

9999  continue

c    OLD EXPRESSION BELOW
c     fit to Fe vapor pressure from email from Frank Ferguson
c     dated 3/25/03
c       pvap_fe = 10.d0**(-19140./t + 8.641)
c       pvap_fe = pvap_fe*1.332*1e3  ! convert from Torr to dyne/cm^2

c    NEW EXPRESSION from Channon Visscher, correspondance on 6/3/11, added 7/27/11 (cvm)

       pvap_fe = 10.d0**(7.09-20833./t)
       pvap_fe = pvap_fe * 1e6   ! convert from bars to dyne/cm^2
 
      return
      end
