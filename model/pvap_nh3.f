c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function pvap_nh3( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over NH3
c
c     input temperature in K
c
c     Marley/Beckmann Jul-2000 
c

      implicit none

      double precision t

      pvap_nh3 = dexp(-86596./t**2 - 2161./t + 10.53)
   
      pvap_nh3 = pvap_nh3*1e6    ! convert from bars to dyne/cm^2

      return
      end
