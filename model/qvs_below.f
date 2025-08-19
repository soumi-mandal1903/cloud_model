      double precision function qvs_below( p_test )
c
c   calculate saturation mixing ratio for a gas, 
c   extrapolated below model domain
c
c   A. Ackerman Nov-2001
c
      implicit none

c   Declare externals

      external pvap_gas
      double precision pvap_gas

c   Declare common storage 

      double precision qv_dtdlnp, qv_p, qv_t, qv_factor
      character*10 qv_gas_name

      common / qvs_below_block /
     $  qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name

c  Declare passed arguments and local storage

      double precision p_test, t_test, pvap_test


c  Extrapolate temperature lapse rate to test pressure

      t_test = qv_t + log( qv_p / p_test )*qv_dtdlnp
 
c  Compute saturation mixing ratio

      pvap_test = pvap_gas( qv_gas_name, t_test, p_test )
      qvs_below = qv_factor * pvap_test / p_test

      return
      end
