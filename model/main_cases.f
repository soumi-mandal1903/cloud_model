      program main

c   --------------------------------------------------------------------
c
c   Defines atmosphere and condensibles;
c   calculates Mie single-scattering efficiencies over a fixed 
c     particle radius grid;
c   calculates eddy-sedimentation vertical profiles of lognormal size
c     distributions of condensates (eddysed);
c   calculates optical properties integrated over particle size distributions.
c  
c   A. Ackerman Feb-2000
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'


c   Declare local storage

      integer NCASE, ncases, icase
      parameter( NCASE = 20 )

      logical do_ttyo, do_fileo, do_subcloud, do_optics, read_mie
      logical do_cases

      integer nz, ngas, iz, igas, itop, ibot, incr
      integer nwave, nrad, moment_rw, iz_max, nsub_max

      character*10 gas_name(MAXNGAS)
      character*80 input_fname

      double precision z(MAXNZ)
      double precision p(MAXNZ)
      double precision t(MAXNZ)
      double precision z_top(MAXNZ+1)
      double precision p_top(MAXNZ+1)
      double precision t_top(MAXNZ+1)
      double precision dlnp(MAXNZ)
      double precision chf(MAXNZ)
      double precision lapse_ratio(MAXNZ)
      double precision kz(MAXNZ)
      double precision gas_mmr(MAXNGAS)
      double precision qt(MAXNZ,MAXNGAS)
      double precision qc(MAXNZ,MAXNGAS)
      double precision ndz(MAXNZ,MAXNGAS)
      double precision sig(MAXNZ,MAXNGAS)
      double precision rainf(MAXNZ,MAXNGAS)
      double precision cloudf(MAXNZ)
      double precision cloudf_gas(MAXNGAS)
      double precision rg(MAXNZ,MAXNGAS)
      double precision reff(MAXNZ,MAXNGAS)
      double precision wave(MAXNWAVE)
      double precision radius(MAXNRAD,MAXNGAS)
      double precision dr(MAXNRAD,MAXNGAS)
      double precision qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision cos_qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision qext(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision opd(MAXNZ,MAXNWAVE)
      double precision w0(MAXNZ,MAXNWAVE)
      double precision g0(MAXNZ,MAXNWAVE)
      double precision opd_gas(MAXNZ,MAXNGAS)
      double precision rainf_case(NCASE)
      double precision reff_qc(NCASE,MAXNGAS)
      double precision reff_lwp(NCASE,MAXNGAS)
      double precision opd_case(NCASE,MAXNGAS)
      double precision grav, teff, opd_tot, mw_atmos, kz_min, cloudf_min
      double precision sig_all, rainf_all, rhoc, rhoc_max, lwp
      double precision lwp_layer, opd_eff
      double precision totpath, totpath_gas(NCASE,MAXNGAS)

c   --------------------------------------------------------------------
c   Define molecular weight of atmosphere (g/mol)
      mw_atmos = 2.2

c   Define number of condensing gases and parameters specific to each:
c     <gas_name> is name of gas
c     <gas_mmr> is total mass mixing ratio below cloud (g/g)

      ngas = 1
      gas_name(1) = 'NH3'

c     ngas = 3
c     gas_name(1) = 'H2O'
c     gas_name(2) = 'Fe'
c     gas_name(3) = 'MgSiO3'

c     ngas = 6
c     gas_name(1) = 'H2O'
c     gas_name(2) = 'Fe'
c     gas_name(3) = 'MgSiO3'
c     gas_name(4) = 'KCl'
c     gas_name(5) = 'NH3'
c     gas_name(6) = 'CH4'

      if( ngas .gt. MAXNGAS )then 
        print*, 'main(): ngas > MAXNGAS'
        stop 1
      endif

      do igas = 1, ngas

        if( gas_name(igas) .eq. 'CH4' )then
          gas_mmr(igas) = 5.9e-4 * (16./mw_atmos)

        elseif( gas_name(igas) .eq. 'NH3' )then
c WST86 (West, Strobel, and Tomasko, 1986) solar
c         gas_mmr(igas) = 1.5e-4 * (17./mw_atmos)
c Voyager IRIS data from Kunde et al (1982) @ 0.6 bar
c also corroborated by Brooke et al. (1998), and very close
c to Carlson et al. (1993) value of 2.8e-5 @ 0.6 bar
          gas_mmr(igas) = 3.0e-5 * (17./mw_atmos)

        elseif( gas_name(igas) .eq. 'H2O' )then
          gas_mmr(igas) = 1.06e-3 * (18./mw_atmos)

        elseif( gas_name(igas) .eq. 'Fe' )then
          gas_mmr(igas) = 1.30e-3

        elseif( gas_name(igas) .eq. 'KCl' )then
          gas_mmr(igas) = 2.2e-7 * (74.5/mw_atmos)

        elseif( gas_name(igas) .eq. 'MgSiO3' )then
          gas_mmr(igas) = 2.75e-3

        else
          print*, 'main(): bad gas_name = ', gas_name(igas)
          stop 1
        endif

      enddo

c   --------------------------------------------------------------------
c   User-defined model parameters for dynamics and microphysics

c   minimum eddy diffusion coefficient (cm^2/s)
      kz_min = 1e5      ! baseline

c   geometric standard deviation for lognormal size distributions
c     sig_all = 2.0     ! baseline
      sig_all = 2.0

c   rain factor 
c     rainf_all = 3.    ! baseline
      rainf_all = 3.

c   minimum cloud coverage (a diagnostic not applied to the clouds)
c     cloudf_min = 1.   ! old baseline
      cloudf_min = 0.75

c   minimum cloud coverage
c     nsub_max = 64     ! baseline
      nsub_max = 64

c   extra cloudy layers below cloud base in calc_optics()
      do_subcloud = .false.

c   optionally (when do_cases = .true.) calculate cloud profiles for
c   a number of different values of rainf; otherwise, just calculate
c   profiles for one value of rainf

      do_cases = .true.
      if( do_cases )then
        rainf_case(1)   = 0.1
        rainf_case(2)   = 0.2
        rainf_case(3)   = 0.3
        rainf_case(4)   = 0.4
        rainf_case(5)   = 0.5
        rainf_case(6)   = 0.6
        rainf_case(7)   = 0.7
        rainf_case(8)   = 0.8
        rainf_case(9)   = 0.9
        rainf_case(10)  = 1.0
        rainf_case(11)  = 1.5
        rainf_case(12)  = 2.
        rainf_case(13)  = 3.
        rainf_case(14)  = 4.
        rainf_case(15)  = 5.
        rainf_case(16)  = 6.
        rainf_case(17)  = 7.
        rainf_case(18)  = 8.
        rainf_case(19)  = 9.
        rainf_case(20)  = 10.
        ncases = NCASE
      else
        ncases = 1
      endif

c   --------------------------------------------------------------------
c   Define and initialize the optics arrays 
    
c     do_optics = .false.
c     read_mie = .true.
      do_optics = .true.
      read_mie = .false.

      call init_optics( do_optics, read_mie, ngas, gas_name,
     $  nwave, wave, nrad, radius, dr, qscat, qext, cos_qscat )

c   --------------------------------------------------------------------
c   Read the atmospheric profiles
c
c   Voyager profiles used in Ackerman and Marley, 2001
c
      input_fname='input/profiles/voyager.input'
      call read_voyager( input_fname, 
     $     grav, teff, nz, z, z_top, p, p_top, t, t_top, 
     $     dlnp, chf, lapse_ratio )
      
c     input_fname='input/profiles/t2000_g1e5.input'
c     input_fname='input/profiles/t1500_g1e5.input'
c     input_fname='input/profiles/t1000_g1e5.input'
c     input_fname='input/profiles/t900_g1e5.input'
c     input_fname='input/profiles/t500_g1e4.input'
c     input_fname='input/profiles/t128_g1e4.input'
c     call read_marley( input_fname, 
c    $     grav, teff, nz, z, z_top, p, p_top, t, t_top, 
c    $     dlnp, chf, lapse_ratio )

c   --------------------------------------------------------------------
c   loop over case

      do icase = 1, ncases

c   --------------------------------------------------------------------
c   Define geometric std deviations and rain factors for
c   each condensates at each altitude
      do igas = 1, ngas
        do iz = 1, nz
          sig(iz,igas) = sig_all
          if( do_cases )then
            rainf(iz,igas) = rainf_case(icase)
          else
            rainf(iz,igas) = rainf_all
          endif
        enddo
      enddo

c   --------------------------------------------------------------------
c   Calculate advective-diffusive balance for all condensing species
c   over all altitudes

      call eddysed( 
     $  grav, teff, kz_min, cloudf_min, nsub_max, sig, rainf,
     $  nz, z, z_top, p, p_top, t, t_top, dlnp, chf, lapse_ratio,
     $  ngas, gas_name, gas_mmr, mw_atmos, 
     $  kz, qt, qc, ndz, rg, reff, cloudf )

c   --------------------------------------------------------------------
c   Calculate optical properties 

      call calc_optics( do_subcloud, nz, ngas, nwave, nrad,
     $  z, gas_name, radius, dr, qscat, qext, cos_qscat,
     $  ndz, sig, rg, opd, opd_gas, w0, g0, opd_tot )

c   --------------------------------------------------------------------
c   Compute diagnostics:
c     reff_qc is eff. radius at condensate conc. maximum (um)
c     reff_av is eff. radius from total optical depth and condensate path (um)
c   Intermediates:
c     rhoc_max is condensate conc. maximum (g/cm^3)
c     lwp is vertical path of condensate (liquid water path) (g/cm^2)
c     opd_eff is optical depth / (1.5 * particle density) for
c     geometric scatterers

      if( z(2) .gt. z(1) )then
        itop = nz
        ibot = 1
        incr = -1
      else
        itop = 1
        ibot = nz
        incr = 1
      endif

      do igas = 1, ngas

        rhoc_max = 0.
        lwp = 0.
        totpath = 0.
        opd_eff = 0.

        do iz = itop, ibot, incr

          lwp_layer = qc(iz,igas) * (p_top(iz+incr)-p_top(iz)) / grav
          lwp = lwp + lwp_layer
          totpath = totpath + (qc(iz,igas)+qt(iz,igas)) * 
     $      (p_top(iz+incr)-p_top(iz)) / grav

          if( reff(iz,igas) .gt. 0. )then
            opd_eff = opd_eff + lwp_layer/reff(iz,igas)
          endif

          rhoc = qc(iz,igas) * p(iz)/(R_GAS/mw_atmos*t(iz))

          if( rhoc .gt. rhoc_max )then
            rhoc_max = rhoc
            iz_max = iz
          endif

        enddo

        reff_qc(icase,igas) = reff(iz_max,igas)*1e4
        opd_case(icase,igas) = opd_gas(ibot,igas)

        if( opd_eff .gt. 0. )then
          reff_lwp(icase,igas) = lwp/opd_eff*1e4
          cloudf_gas(igas) = cloudf(iz_max)
        else
          reff_lwp(icase,igas) = 0.
          cloudf_gas(igas) = 0.
        endif

        totpath_gas(icase,igas) = totpath

      enddo

c   --------------------------------------------------------------------
c   Select output formats
 
      if( do_cases )then
        do_ttyo = .false.
      else
        do_ttyo = .true.
      endif

      do_fileo = .true.

c   --------------------------------------------------------------------
c   terminal output

      if( do_ttyo )then

        print*, 'main():'
        print*, ''
        print*, 'Teff = ', int(teff+0.5), ' K'
        print*, 'grav = ', int(grav/1e2+0.5), ' m/s^2'

        do igas = 1, ngas

          print*, ''
          print*, 'condensing gas = ', gas_name(igas)
          print*, ''
          print*, ' P(bar)   z(m)     T(K) kz(cm^2/s) ' 
     $      // ' reff(um) qc(g/g) col_opd cloudf'     

          do iz = itop, ibot, incr
            write( *,
     $        '( 1pe9.1, 1pe9.2, 0pf7.0, 1p9e9.1)' )
     $        p(iz)/1e6, z(iz)/1e2, t(iz), 
     $        kz(iz), 
     $        reff(iz,igas)*1e4, qc(iz,igas), opd_gas(iz,igas),
     $        cloudf(iz), lapse_ratio(iz)
          enddo
  
          print*,''
          print*,' geometric scattering optical depth = ',
     $      opd_gas(ibot,igas)
          print*,' effective radius (um) at condensate conc. max. = ',
     $      reff_qc(icase,igas)
          print*,' effective radius (um) from lwp, opd = ',
     $      reff_lwp(icase,igas)
          print*,' fractional coverage at cloud base = ',
     $      cloudf_gas(igas)

        enddo

        print*, ''
        print*, ' total optical depth for all ',
     $    'geometric scatterers = ', opd_tot

      endif

c   --------------------------------------------------------------------
c   file output

      if( do_fileo )then

        open( LUNIO, file='eddysed.out', status='replace' )

        write( LUNIO, * ) teff, grav, nz, ngas, sig_all, rainf_all
      
        do igas = 1, ngas
          write( LUNIO, * ) gas_name(igas)
          write( LUNIO, * ) gas_mmr(igas), opd_gas(ibot,igas)
        enddo

        do iz = 1, nz
          write( LUNIO, '(1p20e11.3)' ) z(iz), p(iz), t(iz), 
     $      ( qt(iz,igas), igas=1,ngas ),
     $      ( qc(iz,igas), igas=1,ngas ),
     $      ( reff(iz,igas), igas=1,ngas ),
     $      ( rg(iz,igas), igas=1,ngas ),
     $      ( opd_gas(iz,igas), igas=1,ngas )
        enddo

        close( LUNIO )

      endif

c   --------------------------------------------------------------------
c   bottom of icase loop
c   --------------------------------------------------------------------

      enddo

      if( do_cases )then
        igas = 1
        write(*,'("rainf    = ",1p,20(e8.2,", "),e8.2)') 
     $     (rainf_case(icase),icase=1,NCASE)
        write(*,'("reff lwp = ",1p,20(e8.2,", "),e8.2)') 
     $     (reff_lwp(icase,igas),icase=1,NCASE)
        write(*,'("opd      = ",1p,20(e8.2,", "),e8.2)') 
     $     (opd_case(icase,igas),icase=1,NCASE)
        write(*,'("     totpath = ",1p,20(e8.2,", "),e8.2)') 
     $     (totpath_gas(icase,igas),icase=1,NCASE)
      endif

      end
