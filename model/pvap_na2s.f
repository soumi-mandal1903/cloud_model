c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      double precision function pvap_na2s( t )

c     calculate saturation vapor pressure (dyne/cm^2) of Na2S
c     (Channon's Email, 2011/06/03
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_na2s_bars, metallicityMH

c     Na vapor pressure above cloud 
      metallicityMH=0.0
      pvap_na2s_bars = 10.d0**(8.5497-13889./t-0.5*metallicityMH)
c     Then convert from bars to dynes/cm^2    
      pvap_na2s = pvap_na2s_bars*1e6   

      return
      end
