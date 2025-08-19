c     cccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      double precision function pvap_mg2sio4( t , p )
c
c     calculate saturation vapor pressure (dyne/cm^2) over Mg2SiO4
c     Kozasa et al. Ap J. 344 325
c
c     input temperature in K
c
c     A. Ackerman Feb-2000 modifications for silicates by M. Kress March-2000
c
      implicit none

      double precision t, p, metallicityMH

c     OLD EXPRESSION BELOW 
c       pvap_mg2sio4 = exp(-62279./(t) + 20.94)
c       pvap_mg2sio4 = pvap_mg2sio4*1e6    ! convert from bars to dyne/cm^2

c    NEW EXPRESSION from Channon Visscher, correspondance on 7/27/11, added 8/1/11 (cvm)
c       metallicityMH = 0.0
c       pvap_mg2sio4 = 10.d0**(20.04 - 40540./t - 2*metallicityMH)
c       pvap_mg2sio4 = pvap_mg2sio4 * 1e6  ! convert from bars to dyne/cm^2

c Another new expression from Channon Visscher, correspondance on 10/6/11, includes total pressure dependence and met dep. 
        metallicityMH = 0.0
        pvap_mg2sio4 = 10.d0**(-32488./t + 14.88 - 0.2*log10(p*1e6) 
     $                  - 1.4*metallicityMH) * 1e6 !convered from bars to dynes/cm2

      return
      end
