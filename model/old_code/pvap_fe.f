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

      if( t .gt. 1800. )then
            pvap_fe = exp(-37120./t + 9.86)
      else 
            pvap_fe = exp(-47664./t + 15.71)
      endif
   
      pvap_fe = pvap_fe*1e6    ! convert from bars to dyne/cm^2
 
      return
      end
