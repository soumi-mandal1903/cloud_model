      subroutine eddysed( grav_in, teff, sig, rainf, 
     $  nz, z, p, t, dlnp, chf, lapse_ratio,
     $  ngas, gas_name, gas_mmr, mw_atmos_in, 
     $  kz, lhf, qt, qc, ndz, rg, reff )

c   --------------------------------------------------------------------
c
c   Given an atmosphere and condensates, calculate size and concentration
c   of condensates in balance between eddy diffusion and sedimentation.
c
c   input scalars:
c
c     grav_in      gravitational acceleration (cm/s^2)
c     teff         effective temperature (K)
c     nz           number of layers
c     ngas         number of condensing gases
c
c   input vectors:
c
c     sig          geometric std deviation of lognormal size distribution
c     rainf        ratio of microphysical sed flux to eddy sed flux
c     z            layer altitude (cm)
c     p            layer pressure (dyne/cm^2)
c     t            layer temperature (cm)
c     dlnp         ln(bot/top) pressures at bottom and top of layer
c     chf          layer convective heat flux (erg/s/cm^2)
c     lapse_ratio  ratio of lapse rate w/in layer to adiabat
c     gas_name     names of condensing gases
c     gas_mmr      mass mixing ratio of gas below cloud base (g/g)
c     mw_atmos_in  molecular weight of atmosphere (g/mol)
c
c
c   output vectors:
c
c     kz           eddy diffusion coefficient (cm^2/s)
c     lhf          total latent heat flux (erg/s/cm^2)
c     qt           total (gas+condensed) mixing ratio of condensate (g/g)
c     qc           mixing ratio of condensed condensate (g/g)
c     ndz          number column density of condensate (cm^-3)
c     rg           geometric mean radius of lognormal size distribution
c     reff         droplet effective radius (second moment of size distrib, cm)
c   
c
c   A. Ackerman Feb-2000
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'


c   Declare externals
    
      double precision pvap_ch4, pvap_nh3, pvap_h2o
      double precision pvap_fe, pvap_kcl, pvap_mgsio3
      double precision vfall, advdiff

      external pvap_ch4, pvap_nh3, pvap_h2o
      external pvap_fe, pvap_kcl, pvap_mgsio3
      external vfall, advdiff
     

c   Declare common storage for find_root( vfall )

      double precision grav, mw_atmos, mfp, visc
      double precision t_layer, p_layer, rho_p, mw_cloud

      common / vfall_block /
     $  grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p, mw_cloud
 

c   Declare common storage for find_root( advdiff )

      double precision qbelow, qvs, mixl_layer, zdiff, rainf_layer

      common / advdiff_block /
     $  qbelow, qvs, mixl_layer, zdiff, rainf_layer


c   Declare local storage

      integer cloud_base(MAXNGAS)
      integer nz, ngas, iz, igas, status
      integer itop, ibot, incr
      character*(*) gas_name(ngas)
      double precision z(MAXNZ)
      double precision p(MAXNZ)
      double precision t(MAXNZ)
      double precision dlnp(MAXNZ)
      double precision chf(MAXNZ)
      double precision lapse_ratio(MAXNZ)
      double precision kz(MAXNZ)
      double precision lhf(MAXNZ)
      double precision mixl(MAXNZ)
      double precision gas_mmr(MAXNGAS)
      double precision qt(MAXNZ,MAXNGAS)
      double precision qc(MAXNZ,MAXNGAS)
      double precision ndz(MAXNZ,MAXNGAS)
      double precision sig(MAXNZ,MAXNGAS)
      double precision rg(MAXNZ,MAXNGAS)
      double precision rainf(MAXNZ,MAXNGAS)
      double precision reff(MAXNZ,MAXNGAS)
      double precision rw(MAXNZ,MAXNGAS)
      double precision dz(MAXNZ)
      double precision rho_atmos, scale_h, n_atmos, mw_atmos_in
      double precision grav_in, teff, r_atmos, supsat, w_convect
      double precision alpha, lnsig2, opd_expect, sig_alpha
      double precision c_p, d_molecule, rlo, rhi, delta_v, pvap
      double precision qlo, qhi, delta_q, cond_path, fs, eps_k
      double precision lh, lh_sum, coeff_a, coeff_b, cubic_arg

c   Get values from argument list (also load into common vfall_block)
      grav = grav_in
      mw_atmos = mw_atmos_in

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

c   cloud base indices 
      do igas = 1, ngas
        cloud_base(igas) = 0
      enddo

c   saturation factor
      supsat = 0
      fs = supsat + 1
    
c   --------------------------------------------------------------------
c   Loop over atmospheric layers from the bottom up

      if( z(2) .gt. z(1) )then
        ibot = 1
        itop = nz
        incr = 1
      else
        ibot = nz
        itop = 1 
        incr = -1 
      endif

      do iz = ibot, itop, incr

c   atmospheric density (g/cm^3)
        rho_atmos = p(iz) / ( r_atmos * t(iz) )

c   atmospheric scale height (cm)
        scale_h = r_atmos * t(iz) / grav

c   physical thickness of model layer (cm)
        dz(iz) = scale_h*dlnp(iz)

c   convective mixing length scale (cm): no less than 1/10 scale height
        mixl(iz) = max( 0.1, lapse_ratio(iz) ) * scale_h
c hack
c       mixl(iz) = scale_h

c   atmospheric number density (molecules/cm^3)
        n_atmos = p(iz) / ( K_BOLTZ*t(iz) )

c   atmospheric mean free path (cm)
        mfp = 1. / ( sqrt(2.)*n_atmos*PI*d_molecule**2 )

c   atmospheric viscosity (dyne s/cm^2)
        visc = 5./16.*sqrt( PI*K_BOLTZ*t(iz)*(mw_atmos/AVOGADRO)) /
     $    ( PI*d_molecule**2 ) /
     $    ( 1.22 * ( t(iz) / eps_k )**(-0.16) )

c   Load layer temperature and pressure into common vfall_block 
        t_layer = t(iz)
        p_layer = p(iz)

c   range of particle radii to search (cm)
        rlo = 1.e-10
        rhi = 10.

c   sum of latent heat fluxes / w_convect
        lh_sum = 0.

c   --------------------------------------------------------------------
c   Loop over condensates
c     <mw_cloud> is molecular weight of condensate (g/mol)
c     <lh> is latent heat of sublimation (erg/mol)
c     <pvap> is vapor pressure of condensate (dyne/cm^2)
c     <rho_p> is density when condensed (g/cm^3)

        do igas = 1, ngas 

          if( gas_name(igas) .eq. 'CH4' )then
            pvap = pvap_ch4( t(iz) ) 
            mw_cloud = 16.
            lh = 0.                    ! need a value
c   V.G. Manzhelii and A.M. Tolkachev, Sov. Phys. Solid State 5, 2506 (1964)
c           rho_p = 0.42               ! liquid
            rho_p = 0.49               ! solid

          elseif( gas_name(igas) .eq. 'NH3' )then
            pvap = pvap_nh3( t(iz) ) 
            mw_cloud = 17.
c   Atkins, Physical Chemistry, 6th Ed. (1998)
            lh = 28.6e10
c   V.G. Manzhelii and A.M. Tolkachev, Sov. Phys. Solid State 5, 2506 (1964)
c           rho_p = 0.68               ! liquid
            rho_p = 0.84               ! solid

          elseif( gas_name(igas) .eq. 'H2O' )then
            pvap = pvap_h2o( t(iz) ) 
            mw_cloud = 18.
            lh = 46.7e10
c   Pruppacher and Klett (p. 86 of 1st edition, 1978)
c           rho_p = 1.0                ! liquid
            rho_p = 0.93               ! solid (T = 213 K)

          elseif( gas_name(igas) .eq. 'Fe' )then
            pvap = pvap_fe( t(iz) ) 
            mw_cloud = 55.845
            lh = 0.                    ! need a value
c   Lodders and Fegley (1998)
            rho_p = 7.875

          elseif( gas_name(igas) .eq. 'KCl' )then
            pvap = pvap_kcl( t(iz) ) 
            mw_cloud = 74.5
            lh = 0.                    ! need a value
c   from where?
            rho_p = 1.99

          elseif( gas_name(igas) .eq. 'MgSiO3' )then
            pvap = pvap_mgsio3( t(iz) ) 
            mw_cloud = 100.4
            lh = 0.                    ! need a value
c   Lodders and Fegley (1998)
            rho_p = 3.192

          else
            print*, 'eddysed(): bad igas = ', igas
            stop 1
          endif

c   mass mixing ratio of (super)saturated vapor (g/g)
          qvs = fs*pvap / ( (R_GAS/mw_cloud) * t(iz) ) / rho_atmos


c   Evaluate vapor mass mixing ratio in underlying layer and difference in
c   altitude between mid-pt of layer and mid-pt of underlying layer
c   (ambiguous for bottom layer)

          if( iz .eq. ibot )then
            qbelow = gas_mmr(igas)
            zdiff = dz(iz)/2.
          else
            qbelow = qt(iz-incr,igas)
            zdiff = z(iz) - z(iz-incr)
          endif

c   --------------------------------------------------------------------
c   Layer is cloud free

          if( qbelow .lt. qvs )then
 
            qt(iz,igas) = qbelow
            qc(iz,igas) = 0.
            rg(iz,igas) = 0.

          else

c   --------------------------------------------------------------------
c   Find microphysics/dynamics solution for cloudy layers

c   Load values into common advdiff_block 

            rainf_layer = rainf(iz,igas)
            mixl_layer = mixl(iz)

c   range of mixing ratios to search (g/g)
            qhi = qbelow
            qlo = qhi / 1e3

c   precision of advective-diffusive solution (g/g)
            delta_q = qbelow / 1000.

c   Find total condensate mixing ratio for layer

            call find_root( advdiff, ZERO, qlo, qhi, delta_q, 
     $        qt(iz,igas), status )

            if( status .ne. 0 )then
              write(*,'(a,i3,3a,i4,1pe10.3)')
     $        'find_root(advdiff) status = ',status,
     $        ' for ',gas_name(igas),' at iz,z = ',iz,z(iz)/1e5
              write(*,*)
            endif

            qc(iz,igas) = max( 0., qt(iz,igas) - qvs )

c   Sum latent heat fluxes due to precipitation 
c   (divided by w_convect, which is not yet known)
            lh_sum = lh_sum 
     $        + rainf(iz,igas) * qc(iz,igas)
     $        * lh / ( mw_cloud * c_p / r_atmos )

c   Add latent heat flux due to turbulent mixing
c           if( iz .ne. ibot )then
c             lh_sum = lh_sum 
c    $          - mixl(iz-incr)*(qt(iz,igas)-qt(iz-incr,igas)) / zdiff
cc   $          * lh / ( mw_cloud * c_p / r_atmos )
c    $          * lh / mw_cloud
c           endif

          endif
        enddo

c   Calculate convective velocity scale (cm/s) such that the convective
c   heat flux is the sum of dissipation plus the latent heat flux.
c   For cloudy layers, w_convect is solved from a cubic equation.

        coeff_a = lh_sum
        coeff_b = - chf(iz) / ( rho_atmos * c_p / r_atmos )

        if( lh_sum .eq. 0. )then
          w_convect = (-coeff_b)**(1./3.)
        else
          cubic_arg = coeff_b**2 / 4. + coeff_a**3 / 27.
      print*,'cubic_arg=',cubic_arg
          w_convect = ( -coeff_b/2. + sqrt( cubic_arg ) )**(1./3.)
     $              + ( -coeff_b/2. - sqrt( cubic_arg ) )**(1./3.)
        endif
c dbg
      if( lh_sum .gt. 0. )then
      igas = 1
      print*,'root=',w_convect**3 + coeff_a*w_convect + coeff_b
      print*,'iz, precip,turb LH flux=',iz,
     $        1e-3*lh_sum*w_convect*rho_atmos*c_p/r_atmos
c    $          rainf(iz,igas) * qc(iz,igas)
c    $        * lh / ( mw_cloud )*1e-3*w_convect*rho_atmos
c    $      , - mixl(iz)*(qt(iz,igas)-qt(iz-incr,igas)) / zdiff
c    $        * lh / ( mw_cloud )*1e-3*w_convect*rho_atmos
      endif


c   latent heat flux
        lhf(iz) = w_convect * lh_sum * rho_atmos * c_p / r_atmos

c   vertical eddy diffusion coefficient (cm^2/s) for diagnostic purposes only
        kz(iz) = w_convect * mixl(iz)

c   --------------------------------------------------------------------
c   Calculate optical properties for cloudy layers for each condensate

        do igas = 1, ngas

          if( qc(iz,igas) .gt. 0. )then


c   Calculate vertical index of cloud base

            if( cloud_base(igas) .eq. 0 )then
              cloud_base(igas) = iz
            endif


c   Find <rw> corresponding to <w_convect> using function vfall()

c   precision of vfall solution (cm/s)
            delta_v = w_convect / 1000.

            call find_root( vfall, w_convect, rlo, rhi, delta_v, 
     $        rw(iz,igas), status )

            if( status .ne. 0 )then
              write(*,'(a,i3,3a,i4,1p2e10.3)')
     $          'find_root(vfall) status = ',status,
     $          ' for ',gas_name(igas),' at iz,z,r = ',
     $          iz,z(iz)/1e5, rw(iz,igas)*1e4
              print*,''
            endif


c   geometric std dev of lognormal size distribution
            lnsig2 = 0.5*log( sig(iz,igas) )**2


c   Compute exponent in vfall = w_convect r^alpha

c   sigma floor for the purpose of alpha calculation
            sig_alpha = max( 1.1, sig(iz,igas) )

            if( rainf(iz,igas) .gt. 1 )then

c   Bulk of precip at r > rw: exponent between rw and rw*sig
              alpha = log(
     $          vfall( rw(iz,igas)*sig_alpha ) / w_convect )
     $          / log( sig_alpha )

            else

c   Bulk of precip at r < rw: exponent between rw/sig and rw
              alpha = log(
     $          w_convect / vfall( rw(iz,igas)/sig_alpha) )
     $          / log( sig_alpha )

            endif

c   geometric mean radius of lognormal size distribution
            rg(iz,igas) = rainf(iz,igas)**(1./alpha) *
     $        rw(iz,igas) * exp( -(alpha+6)*lnsig2 )

c   droplet effective radius (cm)
            reff(iz,igas) = rg(iz,igas)*exp( 5*lnsig2 )

c   column droplet number concentration (cm^-2)
            ndz(iz,igas) = 3*rho_atmos*qc(iz,igas)*dz(iz) /
     $        ( 4*PI*rho_p*rg(iz,igas)**3 ) * exp( -9*lnsig2 )

c dbg
c           if( iz .eq. cloud_base(igas) )then
c             print*,'alpha,rg,reff=',
c     $         alpha,rg(iz,igas)*1e4,reff(iz,igas)*1e4
c           endif

          else

            rg(iz,igas) = 0.
            reff(iz,igas) = 0.
            ndz(iz,igas) = 0.

          endif
 
c   Calculate and print analytic optical depth -- assumes cloud
c   extends upwards indefinitely

          if( iz .eq. cloud_base(igas) )then

            opd_expect = 1.5*gas_mmr(igas)*p(iz) /
     $        ( grav*( rainf(iz,igas) + 1 )*reff(iz,igas) )

            print*, ' '
            print*, 'eddysed(): condensing gas = ', 
     $       gas_name(igas)
            print*, 'T(cloud_base) = ', t( cloud_base(igas) )
            print*, 'opd_expect = ', opd_expect
            print*, ' '

          endif

        enddo                          ! ngas
      enddo                            ! nz


c   condensate path

      do igas = 1, ngas
        cond_path = 0.
        do iz = 1, nz
          cond_path = cond_path + qc(iz,igas)*dlnp(iz)*p(iz)/grav
        enddo
        if( cond_path .gt. 0. )then
          print*, ' '
          print*, 'eddysed(): condensing gas = ', 
     $     gas_name(igas)
          print*,' condensate path (g/m^2) = ', cond_path*1e4
          print*, ' '
        endif
      enddo

      return
      end
