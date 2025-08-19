c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      double precision function pvap_mns( t )

c     calculate saturation vapor pressure (dyne/cm^2) of MnS
c     (Channon's Email, 2011/06/21)
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_mns_bars, metallicityMH

c     Mn vapor pressure above cloud 
      metallicityMH=0.0
      pvap_mns_bars = 10.d0**(11.5315-23810./t - metallicityMH)
c     Then convert from bars to dynes/cm^2    
      pvap_mns = pvap_mns_bars*1e6   

      return
      end
