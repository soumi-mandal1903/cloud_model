      subroutine layer( nsub_max, 
     $           gas_name, grav, mw_atmos, kz_min, cloudf_min,
     $           mw_cloud, rainf, 
     $           rho_p, supsat, sig_layer, 
     $           cloudf, q_below, 
     $           t_layer, p_layer, kz, chf, 
     $           t_top, t_bot, p_top, p_bot, 
     $           qc_layer, qt_layer, rg_layer, reff_layer, 
     $           ndz_layer, 
     $           report_status_r, report_status_q )

c   --------------------------------------------------------------------
c
c   Calculate layer condensate properties by iterating on optical depth
c   in one model layer (convering on optical depth over sublayers).
c
c   input:
c
c     nsub_max    maximum number of sublayers for mesh refinement
c     gas_name    name of condensing vapor
c     kz_min      minimum eddy diffusion coefficient (cm^2/s)
c     cloudf_min  minimum cloud fractional coverage
c     mw_cloud    molecular weight of condensate (g/mol)
c     rainf       rain factor
c     rho_p       density of condensed vapor (g/cm^3)
c     supsat      fractional supersaturation persisting after condensation
c     sig_layer   geometric std deviation of lognormal size distribution
c     t_layer     temperature at layer mid-pt (K)
c     p_layer     pressure at layer mid-pt (dyne/cm^2)
c     t_top,bot   temperature at top and bottom of layer (K)
c     p_top,bot   pressure at top and bottom of layer (dyne/cm^2)
c
c   output:
c
c     qc_layer    condensate mixing ratio (g/g)
c     qt_layer    gas + condensate mixing ratio (g/g)
c     rg_layer    geometric mean radius of condensate (cm)
c     reff_layer  effective (area-weighted) radius of condensate (cm)
c     ndz_layer   column of particle concentration in layer (#/cm^2)
c     opd_layer   optical depth for conservative geometric scatterers
c
c     report_status_r    report error status for finding rw
c     report_status_q    report error status for finding qt
c   
c
c   A. Ackerman Dec-2001
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'


c   Declare local storage

      logical report_status_r, report_status_q, converge
      integer status_r, status_q
      integer nsub_max, isub, nsub
      character*(*) gas_name
      double precision rainf, rho_p, supsat, sig_layer
      double precision dlnp, r_atmos, scale_h, mixl, scalef_kz, kz
      double precision p_top, t_top, chf, kz_min, cloudf, cloudf_min
      double precision p_bot, t_bot, n_atmos, q_below
      double precision qc_layer, qt_layer, rg_layer, dz_layer
      double precision reff_layer, ndz_layer, qt_below
      double precision lnsig2, dtdlnp, dp_sub, p_bot_sub, p_top_sub
      double precision qc_sub, qt_sub, t_sub, p_sub, dz_sub
      double precision z_layer, d_molecule, eps_k, c_p, lapse_ratio
      double precision rg_sub, reff_sub, ndz_sub, w_convect
      double precision qt_bot_sub, qt_top
      double precision opd_layer, opd_test, dp_layer, rho_atmos
      double precision grav, mw_atmos, mfp, visc
      double precision t_layer, p_layer, mw_cloud

c   Set error return codes to zero
      status_r = 0
      status_q = 0

c   Number of levels of grid refinement used 
      nsub = 1

c   diameter of atmospheric molecule (cm) (Rosner, 2000)
c   (3.711e-8 for air, 3.798e-8 for N2, 2.827e-8 for H2)
      d_molecule = 2.827e-8

c   parameter in Lennard-Jones potential for viscosity (K) (Rosner, 2000)
c   (78.6 for air, 71.4 for N2, 59.7 for H2)
      eps_k = 59.7

c   specific gas constant for atmosphere (erg/K/g)
      r_atmos = R_GAS / mw_atmos

c   specific heat of atmosphere (erg/K/g)
      c_p = 7./2. * r_atmos

c   pressure thickness of layer
      dp_layer = p_bot - p_top
      dlnp = log( p_bot/p_top )

c   temperature gradient 
      dtdlnp = ( t_top - t_bot ) / dlnp
      lapse_ratio = ( t_bot - t_top ) / dlnp / ( 2./7.*t_layer )

c   atmospheric density (g/cm^3)
      rho_atmos = p_layer / ( r_atmos * t_layer )

c   atmospheric scale height (cm)
      scale_h = r_atmos * t_layer / grav

c   convective mixing length scale (cm): no less than 1/10 scale height
      mixl = max( 0.1d0, lapse_ratio ) * scale_h

c   mixing length = scale height matches Lunine (1989) model
c     mixl = scale_h

c   scale factor for eddy diffusion: 1/3 is baseline
      scalef_kz = 1./3.

c   vertical eddy diffusion coefficient (cm^2/s)
c   from Gierasch and Conrath (1985)
      kz = scalef_kz * scale_h * (mixl/scale_h)**(4./3.) *
     $  ( ( r_atmos*chf ) / ( rho_atmos*c_p ) )**(1./3.)

c   no less than minimum value (for radiative regions)
      kz = max( kz, kz_min )

c   convective velocity scale (cm/s)
      w_convect = kz / mixl

c   cloud fractional coverage
      cloudf = cloudf_min +
     $  max( 0., min( 1., 1.-lapse_ratio )) *
     $  ( 1. - cloudf_min )

c   atmospheric number density (molecules/cm^3)
      n_atmos = p_layer / ( K_BOLTZ*t_layer )

c   atmospheric mean free path (cm)
      mfp = 1. / ( sqrt(2.)*n_atmos*PI*d_molecule**2 )

c   atmospheric viscosity (dyne s/cm^2)
      visc = 5./16.*sqrt( PI*K_BOLTZ*t_layer*(mw_atmos/AVOGADRO)) /
     $  ( PI*d_molecule**2 ) /
     $  ( 1.22 * ( t_layer / eps_k )**(-0.16) )


c   --------------------------------------------------------------------
c   Top of convergence loop

      converge = .false.
      do while ( .not. converge )

c   Zero cumulative values

        qc_layer = 0.
        qt_layer = 0.
        ndz_layer = 0.
        opd_layer = 0.

c   total mixing ratio and pressure at bottom of sub-layer

        qt_bot_sub = q_below
        p_bot_sub = p_bot

c   Loop over sub-layers

        dp_sub = dp_layer / nsub
        do isub = 1, nsub

          qt_below = qt_bot_sub
          p_top_sub = p_bot_sub - dp_sub
          dz_sub = scale_h * log( p_bot_sub/p_top_sub )
          p_sub = 0.5*( p_bot_sub + p_top_sub )
          t_sub = t_bot + log( p_bot/p_sub )*dtdlnp

c   Calculate condensate mixing ratio etc for sub-layer

          call calc_qc( gas_name, rainf, rho_p, mw_cloud,
     $         qt_below, supsat, w_convect, mixl,
     $         dz_sub, grav, mw_atmos, mfp, visc, t_sub, p_sub,
     $         sig_layer, qc_sub, qt_sub, rg_sub, reff_sub, 
     $         ndz_sub, qt_top, status_r, status_q )

c   vertical sums

          qc_layer = qc_layer + qc_sub*dp_sub/grav
          qt_layer = qt_layer + qt_sub*dp_sub/grav
          ndz_layer = ndz_layer + ndz_sub

          if( reff_sub .gt. 0. )then
            opd_layer = opd_layer + 
     $        1.5*qc_sub*dp_sub/grav/(rho_p*reff_sub)
          endif

c   Increment values at bottom of sub-layer

          qt_bot_sub = qt_top
          p_bot_sub = p_top_sub

        enddo

c    Check convergence on optical depth

        if( nsub_max .eq. 1 )then
          converge = .true.
        elseif( nsub .eq. 1 )then
          opd_test = opd_layer
        elseif( opd_layer .eq. 0. .or. nsub .ge. nsub_max )then
          converge = .true.
        elseif( abs( 1. - opd_test/opd_layer ) .le. 1e-2 )then
          converge = .true.
        else
          opd_test = opd_layer
        endif

        nsub = nsub * 2

      enddo
c   --------------------------------------------------------------------
c   Bottom of convergence loop


c     Report problems finding root the first time it happens

      if( status_r .ne. 0 .and. report_status_r )then
        print*, 'layer():'
        write(*,'(a,i3,3a,1pe10.2)')
     $    ' find_root(vfall) status = ',status_r,
     $    ' for ',gas_name,' at p = ',p_layer/1e6
        print*,' there may be more instances not reported'
        print*,''
        print*,'status_r = ',status_r
        report_status_r = .false.
      endif

      if( status_r .ne. 0 .and. report_status_q )then
        print*, 'layer():'
        write(*,'(a,i3,3a,1pe10.2)')
     $    ' find_root(advdiff) status = ',status_q,
     $    ' for ',gas_name,' at iz,p = ',p_layer/1e6
        print*,' there may be more instances not reported'
        print*,''
        print*,'status_q = ',status_q
        report_status_q = .false.
      endif


c   Update properties at bottom of next layer

      q_below = qt_top

c   Get layer averages

      if( opd_layer .gt. 0. )then
        reff_layer = 1.5*qc_layer / (rho_p*opd_layer)
        lnsig2 = 0.5*log( sig_layer )**2
        rg_layer = reff_layer*exp( -5*lnsig2 )
      else
        reff_layer = 0.
        rg_layer = 0.
      endif

      qc_layer = qc_layer*grav / dp_layer
      qt_layer = qt_layer*grav / dp_layer

      return
      end    
