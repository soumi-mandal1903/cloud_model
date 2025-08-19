cc     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cc     ALTERNATE EXPRESSION TO CALCULATE PVAP_KCL
      double precision function pvap_kcl( t )

c     calculate saturation vapor pressure (dyne/cm^2) of KCl
c     (Fit from Lodders data)
c     input temperature in K
c     (cvm)
      implicit none
      double precision t, pvap_kcl_bars

      pvap_kcl_bars = 10.d0**(7.6106 - 11382./t)
c     Then convert from bars to dynes/cm^2    
      pvap_kcl = pvap_kcl_bars*1e6   

      return
      end
