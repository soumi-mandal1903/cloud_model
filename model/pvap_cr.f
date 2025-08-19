c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc


      double precision function pvap_cr( t )

c     calculate saturation vapor pressure (dyne/cm^2) of Cr
c     (Fit from data in email - Fortney 4/22/11) (cvm)
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_cr_bars

c     Cr vapor pressure above cloud 
      pvap_cr_bars = 10.d0**(7.2688-20353./t)
c     Then convert from bars to dynes/cm^2    
      pvap_cr = pvap_cr_bars*1e6   

      return
      end
