      double precision function pvap_mgsio3( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over MgSiO3
c     (Lunine et al., 1988, ApJ, 338, 314, citing Barshay and Lewis 1976)
c
c     input temperature in K
c
c     A. Ackerman Feb-2000 modifications for silicates by M. Kress March-2000
c
      implicit none

      double precision t

      pvap_mgsio3 = exp(-58663./t + 25.37)
    
      pvap_mgsio3 = pvap_mgsio3*1e6    ! convert from bars to dyne/cm^2

      return
      end
