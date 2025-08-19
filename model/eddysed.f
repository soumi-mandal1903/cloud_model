      subroutine eddysed( 
     $  grav, teff, kz_min, cloudf_min, nsub_max, 
     $  supsat, mw_atmos, do_virtual,
     $  nz, z, z_top, p, p_top, t, t_top, chf, 
     $  ngas, gas_name, gas_mmr, gas_mw, rho_p, sig, rainf, 
     $  kz, qt, qc, ndz, rg, reff, cloudf )

c   --------------------------------------------------------------------
c
c   Given an atmosphere and condensates, calculate size and concentration
c   of condensates in balance between eddy diffusion and sedimentation.
c
c   input scalars:
c
c     grav         gravitational acceleration (cm/s^2)
c     teff         effective temperature (K)
c     kz_min       minimum eddy diffusion coefficient (cm^2/s)
c     cloudf_min   minimum cloud fractional coverage
c     nsub_max     maximum number of sublayers for adaptive mesh refinement
c     supsat       supersaturation after condensation (fraction)
c     mw_atmos     molecular weight of atmosphere (g/mol)
c     do_virtual   reduce mixing ratio due to decrease below cloud base
c     nz           number of layers
c     ngas         number of condensing gases
c
c   input vectors:
c
c     z            altitude at layer mid-pt (by pressure) (cm)
c     z_top        altitude at top of layer (cm)
c     p            pressure at layer mid-pt (by pressure) (dyne/cm^2)
c     p_top        pressure at top of layer (dyne/cm^2)
c     t            temperature at layer mid-pt (by pressure) (K)
c     t_top        temperature at top of layer (K)
c     chf          layer convective heat flux (erg/s/cm^2)
c
c     gas_name     names of condensing gases
c     gas_mmr      mass mixing ratio of gas below cloud base (g/g)
c     gas_mw       molecular weight of gas (g/mol)
c     rho_p        density of condensed vapor (g/cm^3)
c     sig          geometric std deviation of lognormal size distribution
c     rainf        ratio of microphysical sed flux to eddy sed flux
c
c   output vectors:
c
c     kz           eddy diffusion coefficient (cm^2/s)
c     qt           total (gas+condensed) mixing ratio of condensate (g/g)
c     qc           mixing ratio of condensed condensate (g/g)
c     ndz          number column density of condensate (cm^-3)
c     rg           geometric mean radius of lognormal size distribution
c     reff         droplet effective radius (second moment of size distrib, cm)
c     cloudf       cloud fractional coverage (not applied to clouds)
c   
c
c   A. Ackerman Feb-2000
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'


c   Declare externals

      double precision qvs_below, pvap_gas
      external qvs_below, pvap_gas


c   Declare common storage for qvs_below

      double precision qv_dtdlnp, qv_p, qv_t, qv_factor
      character*10 qv_gas_name

      common / qvs_below_block /
     $  qv_dtdlnp, qv_p, qv_t, qv_factor, qv_gas_name
 

c   Declare local storage

      logical report_status_r(MAXNGAS)
      logical report_status_q(MAXNGAS)
      integer cloud_base(MAXNGAS)
      logical do_virtual, converge
      integer status_p
      integer nsub_max
      integer nz, ngas, iz, igas
      integer itop, ibot, incr
      character*(*) gas_name(ngas)
      double precision z(MAXNZ)
      double precision p(MAXNZ)
      double precision t(MAXNZ)
      double precision z_top(MAXNZ+1)
      double precision p_top(MAXNZ+1)
      double precision t_top(MAXNZ+1)
      double precision chf(MAXNZ)
      double precision kz(MAXNZ)
      double precision gas_mmr(MAXNGAS)
      double precision gas_mw(MAXNGAS)
      double precision rho_p(MAXNGAS)
      double precision qt(MAXNZ,MAXNGAS)
      double precision qc(MAXNZ,MAXNGAS)
      double precision ndz(MAXNZ,MAXNGAS)
      double precision sig(MAXNZ,MAXNGAS)
      double precision rg(MAXNZ,MAXNGAS)
      double precision rainf(MAXNZ,MAXNGAS)
      double precision reff(MAXNZ,MAXNGAS)
      double precision cloudf(MAXNZ)
      double precision rw(MAXNZ,MAXNGAS)
      double precision qc_path(MAXNGAS)
      double precision q_below(MAXNGAS)
      double precision rho_atmos, n_atmos, mw_atmos
      double precision teff, kz_min
      double precision lnsig2, cloudf_min
      double precision c_p, d_molecule
      double precision qlo, qhi, delta_q, eps_k
      double precision scalef_kz, supsat, qvs_factor
      double precision qc_layer, qt_layer, rg_layer, reff_layer
      double precision ndz_layer, kz_layer, dtdlnp
      double precision pvap, p_lo, p_hi, p_base, t_base
      double precision t_bot, p_bot, qvs, t_layer, p_layer, grav
      double precision opd_below


c   Initialize cloud base indices and logical flags that shut off 
c   reporting root problems after the first time

      do igas = 1, ngas
        cloud_base(igas) = 0
        report_status_r(igas) = .true.
        report_status_q(igas) = .true.
      enddo

c   Calculate indices of bottom and top of domain

      if( z(2) .gt. z(1) )then
        ibot = 1
        itop = nz
        incr = 1
      else
        ibot = nz
        itop = 1 
        incr = -1 
      endif

c   --------------------------------------------------------------------
c   Loop over condensates

      do igas = 1, ngas 

c   Start at bottom of domain: get p, T, and qt 

        t_bot = t_top(ibot-incr)
        p_bot = p_top(ibot-incr)

        q_below(igas) = gas_mmr(igas)

c   Initialize vertical path of condensate <qc_path> to zero
        qc_path(igas) = 0.


c   --------------------------------------------------------------------
c   Adjust mixing ratio at bottom of model domain if bottom layer 
c   is saturated

        if( do_virtual )then

          qvs_factor = (supsat+1)*gas_mw(igas)/mw_atmos
          pvap = pvap_gas( gas_name(igas), t_bot, p_bot )
          qvs = qvs_factor*pvap/p_bot

          if( qvs .le. q_below(igas) )then

c   Find pressure at cloud base

c   parameters for finding root 
            p_lo = p_bot
            p_hi = p_bot * 1e2
            delta_q = q_below(igas) / 1e2

c   temperature gradient
            dtdlnp = ( t_top(ibot) - t_bot ) / log( p_bot/p_top(ibot) )

c   load parameters into qvs_below common block
              
            qv_dtdlnp = dtdlnp
            qv_p = p_bot
            qv_t = t_bot
            qv_gas_name = gas_name(igas) 
            qv_factor = qvs_factor

            call find_rootl( qvs_below, q_below(igas), p_lo, p_hi,
     $           delta_q, p_base, status_p )
              
            if( status_p .ne. 0 )then
              print*, ''
              print*, 'unable to find cloud base pressure in eddysed():'
              print*, ' find_rootl(qvs_below) status = ',status_p,
     $          ' for ', gas_name(igas)
              print*,''
            endif

            t_base = t_bot + log( p_bot/p_base )*dtdlnp

c   Calculate temperature and pressure below bottom layer
c   by adding a virtual layer 

            p_layer = 0.5*( p_bot + p_base )
            t_layer = t_bot + log( p_bot/p_layer )*dtdlnp

c   Calculate qc, qt, rg, reff, and ndz for virtual layer

            call layer( nsub_max, 
     $           gas_name(igas), grav, mw_atmos, kz_min, cloudf_min,
     $           gas_mw(igas), rainf(ibot,igas),
     $           rho_p(igas), supsat, sig(ibot,igas), 
     $           cloudf(ibot), q_below(igas), 
     $           t_layer, p_layer, kz_layer, chf(ibot), 
     $           t_bot, t_base, p_bot, p_base,
     $           qc_layer, qt_layer, rg_layer, reff_layer, 
     $           ndz_layer, 
     $           report_status_r(igas), report_status_q(igas) )

c   Report optical depth below domain

            print*, ''
            print*, 'eddysed(): condensing gas = ', gas_name(igas)
            print*, ' cloud base at p, T = ', p_base/1e6, t_base
            print*, ' optical depth below domain = ',
     $          1.5*qc_layer*( p_base - p_bot )/grav /
     $          ( rho_p(igas)*reff_layer )
            print*, ''

          endif
        endif

c   --------------------------------------------------------------------
c   Loop over atmospheric layers from the bottom up
        do iz = ibot, itop, incr

c   Calculate layer qc, qt, rg, reff, and ndz for each layer

          call layer( nsub_max, 
     $         gas_name(igas), grav, mw_atmos, kz_min, cloudf_min,
     $         gas_mw(igas), rainf(iz,igas),
     $         rho_p(igas), supsat, sig(iz,igas), 
     $         cloudf(iz), q_below(igas), 
     $         t(iz), p(iz), kz(iz), chf(iz), 
     $         t_top(iz), t_top(iz-incr), p_top(iz), p_top(iz-incr),
     $         qc(iz,igas), qt(iz,igas), rg(iz,igas), reff(iz,igas), 
     $         ndz(iz,igas), 
     $         report_status_r(igas), report_status_q(igas) )

c   Accumulate vertical path of condensate 

          qc_path(igas) = qc_path(igas) + qc(iz,igas)*
     $       ( p_top(iz-incr) - p_top(iz) )/grav

        enddo                          ! nz

c   Print some diagnostics

      print*, ''
      print*, 'eddysed(): condensing gas = ', gas_name(igas)
      print*, ' condensate path = ', qc_path(igas)*1e4, ' g/m^2'
      print*, ''

      enddo                            ! ngas

      return
      end


