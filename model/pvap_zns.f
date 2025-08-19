c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      double precision function pvap_zns( t )

c     calculate saturation vapor pressure (dyne/cm^2) of ZnS
c     (Channon's Email, 2011/06/21)
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_zns_bars, metallicityMH

c     Zn vapor pressure above cloud 
      metallicityMH=0.0
      pvap_zns_bars = 10.d0**(12.8117-15873./t - metallicityMH)
c     Then convert from bars to dynes/cm^2    
      pvap_zns = pvap_zns_bars*1e6   

      return
      end
