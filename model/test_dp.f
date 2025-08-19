      program main

c   --------------------------------------------------------------------
c
c   find dewpoint temperature
c  
c   A. Ackerman 2001
c
c   --------------------------------------------------------------------

      implicit none

      integer status_t

      double precision mw_atmos, mw_cloud, gas_mmr, dp
      double precision tlo, thi, delta_p, rhs, press, eps

      double precision pvap_ch4, pvap_nh3, pvap_h2o
      double precision pvap_fe, pvap_kcl, pvap_mgsio3

      external pvap_ch4, pvap_nh3, pvap_h2o
      external pvap_fe, pvap_kcl, pvap_mgsio3

c   --------------------------------------------------------------------
      mw_atmos = 29.66

c   water

      mw_cloud = 18.
      gas_mmr = 10.e-3
      press = 1.e6

c   ammonia

c     mw_cloud = 17.
c     gas_mmr = 3e-5*mw_cloud/mw_atmos
c     press = 42250.

c     eps = mw_cloud/mw_atmos

c   range of temperatures to search [K]
      tlo = 100.
      thi = 300.

c   rhs in pvap(t_dew) = rhs [dyne/cm^2]
      rhs = press*gas_mmr*mw_atmos/mw_cloud

c   precision of search [dyne/cm^2]
      delta_p = rhs/100.

      call find_root( pvap_h2o, rhs, tlo, thi, delta_p, 
     $                dp, status_t )

      print*,'delta_p,dewp=',delta_p,dp
 
      if( status_t .ne. 0 )then
        print*,'find_root(pvap) status = ',status_t
      endif

      end
