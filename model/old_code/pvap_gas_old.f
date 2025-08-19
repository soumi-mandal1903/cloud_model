      double precision function pvap_gas( gas_name, t_layer )
c
c   calculate vapor pressure for a gas
c
c   A. Ackerman Nov-2001
c
      implicit none

c  Declare externals

      double precision pvap_ch4, pvap_nh3, pvap_h2o
      double precision pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3

      external pvap_ch4, pvap_nh3, pvap_h2o
      external pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3

c  Declare passed arguments

      double precision t_layer
      character*(*) gas_name

c  Evaluate the vapor pressure using the appropriate function

      if( gas_name .eq. 'CH4' )then
        pvap_gas = pvap_ch4( t_layer ) 
      elseif( gas_name .eq. 'NH3' )then
        pvap_gas = pvap_nh3( t_layer )
      elseif( gas_name .eq. 'H2O' )then
        pvap_gas = pvap_h2o( t_layer )
      elseif( gas_name .eq. 'Fe' )then
        pvap_gas = pvap_fe( t_layer )
      elseif( gas_name .eq. 'KCl' )then
        pvap_gas = pvap_kcl( t_layer )
      elseif( gas_name .eq. 'MgSiO3' )then
        pvap_gas = pvap_mgsio3( t_layer )
      elseif( gas_name .eq. 'Al2O3' )then
        pvap_gas = pvap_al2o3( t_layer )
      else
        print*,'stop in pvap_gas(), bad gas_name = ',gas_name
      endif

      return
      end
