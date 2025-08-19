      subroutine calc_qc( gas_name, rainf_layer, rho_p, mw_cloud, 
     $           q_below, supsat, w_convect, mixl,
     $           dz_layer, grav, mw_atmos, mfp, visc, t_layer, p_layer,
     $           sig_layer, qc_layer, qt_layer, rg_layer, reff_layer, 
     $           ndz_layer, qt_top, status_r, status_q )

c   --------------------------------------------------------------------
c
c   Calculate condensate optical depth and effective radius for a layer,
c   assuming geometric scatterers.
c
c   input:
c
c     gas_name    name of condensing vapor
c     rainf_layer rain factor for layer
c     rho_p       density of condensed vapor (g/cm^3)
c     mw_cloud    molecular weight of condensing vapor (g/mol)
c     q_below     total mixing ratio (vapor+condensate) below layer (g/g)
c     supsat      fractional supersaturation persisting after condensation
c     w_convect   convective velocity scale (cm/s)
c     mixl        convective mixing length scale (cm)
c     dz_layer    thickness of layer (cm)
c     grav        gravitational acceleration (cm/s^2)
c     mw_atmos    molecular weight of atmosphere (g/mol)
c     mfp         atmospheric mean free path (cm)
c     visc        atmospheric dynamic viscosity (dyne s/cm^2)
c     t_layer     temperature of layer mid-pt (K)
c     p_layer     air pressure (dyne/cm^2)
c     sig_layer   geometric std. dev. of lognormal size distribution
c
c   output:
c
c     qc_layer    condensate mixing ratio (g/g)
c     qt_layer    gas + condensate mixing ratio (g/g)
c     rg_layer    geometric mean radius of condensate (cm)
c     reff_layer  effective (area-weighted) radius of condensate (cm)
c     ndz_layer   column of particle concentration in layer (#/cm^2)
c     qt_top      top of layer
c     status_r    error status for finding rw
c     status_q    error status for finding qt
c   
c
c   A. Ackerman Feb-2001
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'


c   Declare externals
    
      double precision pvap_gas, vfall, advdiff
      external pvap_gas, vfall, advdiff
     

c   Declare common storage for find_root( vfall )

      double precision vf_grav, vf_mw_atmos, vf_mfp, vf_visc
      double precision vf_t, vf_p, vf_rhop

      common / vfall_block /
     $  vf_grav, vf_mw_atmos, vf_mfp, vf_visc,
     $  vf_t, vf_p, vf_rhop


c   Declare common storage for find_root( advdiff )

      double precision ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf

      common / advdiff_block /
     $  ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf


c   Declare local storage

      integer status_r, status_q
      character*(*) gas_name
      double precision rainf_layer, rho_p, supsat, fs, pvap, dz_layer
      double precision w_convect, qc_layer, qt_layer, rg_layer, mfp
      double precision reff_layer, ndz_layer, qvs, mw_cloud, q_below
      double precision qlo, qhi, delta_q, t_layer, p_layer, rho_atmos
      double precision rlo, rhi, delta_v, qt_top, mw_atmos, mixl, grav
      double precision rw_layer, lnsig2, sig_layer, sig_alpha, alpha
      double precision visc

c   vapor pressure of condensate (dyne/cm^2)
      pvap = pvap_gas( gas_name, t_layer, p_layer )

c   saturation factor
      fs = supsat + 1

c   atmospheric density (g/cm^3)
      rho_atmos = p_layer / ( R_GAS/mw_atmos * t_layer )

c   mass mixing ratio of saturated vapor (g/g)
      qvs = fs*pvap / ( (R_GAS/mw_cloud) * t_layer ) / rho_atmos

c   --------------------------------------------------------------------
c   Layer is cloud free

      if( q_below .lt. qvs )then
 
        qt_layer = q_below
        qt_top   = q_below
        qc_layer = 0.
        rg_layer = 0.
        reff_layer = 0.
        ndz_layer = 0.

      else

c   --------------------------------------------------------------------
c   Cloudy layer: first calculate qt and qc at top of layer,
c   then calculate layer averages

c   range of mixing ratios to search (g/g)
        qhi = q_below
        qlo = qhi / 1e3

c   precision of advective-diffusive solution (g/g)
        delta_q = q_below / 1000.

c   load parameters into advdiff common block

        ad_qbelow = q_below
        ad_qvs = qvs
        ad_mixl = mixl
        ad_dz = dz_layer
        ad_rainf = rainf_layer

c   Find total condensate mixing ratio at top of layer

        call find_root( advdiff, ZERO, qlo, qhi, delta_q, 
     $       qt_top, status_q )

c   Use trapezoid rule (for now) to calculate layer averages
c   -- should integrate exponential
        qt_layer = 0.5*( q_below + qt_top )

c   Diagnose condensate mixing ratio
        qc_layer = max( 0., qt_layer - qvs )

c   --------------------------------------------------------------------
c   Find <rw> corresponding to <w_convect> using function vfall()

c   load parameters into vfall common block
        vf_grav = grav
        vf_mw_atmos = mw_atmos
        vf_mfp = mfp
        vf_visc = visc
        vf_p = p_layer
        vf_t = t_layer
        vf_rhop = rho_p

c   range of particle radii to search (cm)
        rlo = 1.e-10
        rhi = 10.

c   precision of vfall solution (cm/s)
        delta_v = w_convect / 1000.

        call find_root( vfall, w_convect, rlo, rhi, delta_v, 
     $       rw_layer, status_r )

c   geometric std dev of lognormal size distribution
        lnsig2 = 0.5*log( sig_layer )**2


c   Compute exponent in vfall = w_convect r^alpha

c   sigma floor for the purpose of alpha calculation
        sig_alpha = max( 1.1, sig_layer )

        if( rainf_layer .gt. 1 )then

c   Bulk of precip at r > rw: exponent between rw and rw*sig
          alpha = log(
     $      vfall( rw_layer*sig_alpha ) / w_convect )
     $      / log( sig_alpha )

        else

c   Bulk of precip at r < rw: exponent between rw/sig and rw
          alpha = log(
     $      w_convect / vfall( rw_layer/sig_alpha) )
     $      / log( sig_alpha )

        endif


c   geometric mean radius of lognormal size distribution
        rg_layer = rainf_layer**(1./alpha) *
     $    rw_layer * exp( -(alpha+6)*lnsig2 )

c   droplet effective radius (cm)
        reff_layer = rg_layer*exp( 5*lnsig2 )

c   column droplet number concentration (cm^-2)
        ndz_layer = 3*rho_atmos*qc_layer*dz_layer /
     $    ( 4*PI*rho_p*rg_layer**3 ) * exp( -9*lnsig2 )

      endif

      return
      end
