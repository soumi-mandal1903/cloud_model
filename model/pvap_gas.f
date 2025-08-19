       double precision function pvap_gas( gas_name, t_layer , p_layer)
c
c   calculate vapor pressure for a gas
c
c   A. Ackerman Nov-2001
c
      implicit none

c  Declare externals

      double precision pvap_ch4, pvap_nh3, pvap_h2o
      double precision pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3, 
     $      pvap_mg2sio4
      double precision pvap_mns, pvap_zns, pvap_na2s, pvap_cr

      external pvap_ch4, pvap_nh3, pvap_h2o
      external pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3, pvap_mg2sio4
      external pvap_mns, pvap_zns, pvap_na2s, pvap_cr

c  Declare passed arguments

      double precision t_layer, p_layer
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
      elseif( gas_name .eq. 'Mg2SiO4' )then
        pvap_gas = pvap_mg2sio4( t_layer , p_layer )
      elseif( gas_name .eq. 'Al2O3' )then
        pvap_gas = pvap_al2o3( t_layer )
      elseif( gas_name .eq. 'MnS' )then
	        pvap_gas = pvap_mns( t_layer )
	  elseif( gas_name .eq. 'Na2S' )then
	        pvap_gas = pvap_na2s( t_layer )
	  elseif( gas_name .eq. 'ZnS' )then
            pvap_gas = pvap_zns( t_layer )
      elseif( gas_name .eq. 'Cr' )then
		    pvap_gas = pvap_cr( t_layer )
      else
        print*,'stop in pvap_gas(), bad gas_name = ',gas_name
      endif
      
      return
      end











