      double precision function pvap_mgsio3( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over MgSiO3
c     (Lunine et al., 1986, ApJ, 338, 314, citing Barshay and Lewis 1976)
c
c     input temperature in K
c
c     A. Ackerman Feb-2000 modifications for silicates by M. Kress March-2000
c
      implicit none

      double precision t, metallicityMH

c      pvap_mgsio3 = exp(-58663./t + 25.37)   
c      pvap_mgsio3 = pvap_mgsio3*1e6    ! convert from bars to dyne/cm^2
c      NEW EXPRESSIONS FROM CV
       metallicityMH=0.0
       pvap_mgsio3 = 10.d0**(11.83 - 27250./t - metallicityMH)
       pvap_mgsio3 = 1e6 * pvap_mgsio3 !convert bars -> dynes/cm^2

c Another new expression from Channon Visscher, correspondance on 10/6/11, includes total pressure dependence and met dep. 
c This is the saturation vapor pressure of SiO (which is now the limiting factor, instead of Mg), assuming that the forsterite has 
c already formed underneath. You would need a different expression if only the enstatite cloud existed. 
c        pvap_mgsio3 = 10.d0**(-28665./t + 13.43) * 1e6
       
      return
      end
