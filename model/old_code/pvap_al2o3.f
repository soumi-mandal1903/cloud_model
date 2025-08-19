      double precision function pvap_al2o3( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over Al2O3
c     Kozasa et al. Ap J. 344 325
c
c     input temperature in K
c
c     A. Ackerman Feb-2000 modifications for al2o3 M Marley Nov. 2001
c
      implicit none

      double precision t

c     pvap_al2o3 = 0.0259*exp(-73503./t + 22.005)*1e6
      pvap_al2o3 = exp(-120000./t + 48.78)*1e6

      return
      end
