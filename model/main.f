      subroutine main( grav, 
     $   teff, nz, t_top, p_top, qt, qc, wave, opd, opd_gas, w0, g0, 
     $   gas_mmr, gas_mw, mw_atmos, kz)

      implicit none

      include 'globals.h'

c ============================================================
c External vapor pressure functions
c ============================================================

      double precision pvap_ch4, pvap_nh3, pvap_h2o
      double precision pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3
      double precision pvap_mg2sio4, pvap_mns, pvap_zns, pvap_na2s, pvap_cr

      external pvap_ch4, pvap_nh3, pvap_h2o
      external pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3, pvap_mg2sio4
      external pvap_mns, pvap_zns, pvap_na2s, pvap_cr

c ============================================================
c Local variables
c ============================================================

      logical do_ttyo, do_fileo, do_subcloud, do_optics, read_mie
      logical do_cases, do_virtual

      integer NCASE, ncases, icase
      parameter( NCASE = 20 )

      integer nz, ngas, iz, igas, itop, ibot, incr
      integer nwave, nrad, iz_max, nsub_max, status_t

      character*10 gas_name(MAXNGAS)
      character*80 input_fname

      double precision z(MAXNZ), p(MAXNZ), t(MAXNZ)
      double precision z_top(MAXNZ+1), p_top(MAXNZ+1), t_top(MAXNZ+1)
      double precision chf(MAXNZ), kz(MAXNZ)
      double precision gas_mmr(MAXNGAS), gas_mw(MAXNGAS)
      double precision qt(MAXNZ,MAXNGAS), qc(MAXNZ,MAXNGAS)
      double precision ndz(MAXNZ,MAXNGAS), sig(MAXNZ,MAXNGAS)
      double precision rainf(MAXNZ,MAXNGAS), cloudf(MAXNZ)
      double precision rho_p(MAXNGAS), rg(MAXNZ,MAXNGAS)
      double precision reff(MAXNZ,MAXNGAS)
      double precision wave(MAXNWAVE), radius(MAXNRAD,MAXNGAS)
      double precision dr(MAXNRAD,MAXNGAS)
      double precision qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision cos_qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision qext(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision opd(MAXNZ,MAXNWAVE), w0(MAXNZ,MAXNWAVE)
      double precision g0(MAXNZ,MAXNWAVE), opd_gas(MAXNZ,MAXNGAS)
      double precision dewp(MAXNZ,MAXNGAS)

      double precision grav, teff, opd_tot, mw_atmos
      double precision kz_min, cloudf_min, sig_all, rainf_all
      double precision rhoc, rhoc_max, lwp, lwp_layer, opd_eff
      double precision supsat, totpath, tlo, thi, delta_p, rhs, eps

c ============================================================
c Jovian atmosphere parameters (AM01)
c ============================================================

      mw_atmos = 2.2d0

c Only NH3 cloud
      ngas = 1
      gas_name(1) = 'NH3'

c NH3 molecular weight and deep mixing ratio
      gas_mw(1)  = 17.0d0
      gas_mmr(1) = 1.34d-4 * ( gas_mw(1) / mw_atmos )

c Condensate density (solid NH3)
      rho_p(1) = 0.84d0

c ============================================================
c Microphysics parameters (AM01 baseline)
c ============================================================

      kz_min    = 1.0d5
      sig_all   = 2.0d0
      rainf_all = 3.0d0
      cloudf_min = 0.75d0
      nsub_max  = 64
      do_virtual = .true.
      do_subcloud = .true.
      supsat = 0.0d0

      do_cases = .false.
      ncases = 1

c ============================================================
c Optics initialization
c ============================================================

      do_optics = .true.
      read_mie  = .true.

      call init_optics( do_optics, read_mie, ngas, gas_name,
     $   nwave, wave, nrad, radius, dr, qscat, qext, cos_qscat )

c ============================================================
c Read Voyager T-P profile (AM01)
c ============================================================

      input_fname = 'input/profiles/voyager.input'
      nz = 0

      call read_voyager( input_fname,
     $     grav, teff, nz, z, z_top, p, p_top, t, t_top, chf )

c ============================================================
c Cloud microphysics loop
c ============================================================

      do icase = 1, ncases

      do igas = 1, ngas
        do iz = 1, nz
          sig(iz,igas)   = sig_all
          rainf(iz,igas) = rainf_all
        enddo
      enddo

      call eddysed(
     $  grav, teff, kz_min, cloudf_min, nsub_max,
     $  supsat, mw_atmos, do_virtual,
     $  nz, z, z_top, p, p_top, t, t_top, chf,
     $  ngas, gas_name, gas_mmr, gas_mw, rho_p, sig, rainf,
     $  kz, qt, qc, ndz, rg, reff, cloudf )

      call calc_optics( do_subcloud, nz, ngas, nwave, nrad,
     $  z, gas_name, radius, dr, qscat, qext, cos_qscat,
     $  ndz, sig, rg, opd, opd_gas, w0, g0, opd_tot )

c ============================================================
c Diagnostics
c ============================================================

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

        rhoc_max = 0.d0
        lwp = 0.d0
        totpath = 0.d0
        opd_eff = 0.d0

        do iz = itop, ibot, incr

          lwp_layer = qc(iz,igas) * (p_top(iz+incr)-p_top(iz)) / grav
          lwp = lwp + lwp_layer

          if( reff(iz,igas) .gt. 0.d0 )then
            opd_eff = opd_eff + lwp_layer / reff(iz,igas)
          endif

          rhoc = qc(iz,igas) * p(iz)/(R_GAS/mw_atmos*t(iz))

          if( rhoc .gt. rhoc_max )then
            rhoc_max = rhoc
            iz_max = iz
          endif

        enddo

      enddo

c ============================================================
c Output
c ============================================================

      open( LUNIO, file='eddysed.out', status='replace' )

      write( LUNIO, * ) teff, grav, nz, ngas, sig_all, rainf_all
      write( LUNIO, * ) gas_name(1)
      write( LUNIO, * ) gas_mmr(1), opd_gas(ibot,1)

      do iz = 1, nz
        write( LUNIO, '(1p10e12.4)' ) z(iz), p(iz), t(iz),
     $      qt(iz,1), qc(iz,1), reff(iz,1), rg(iz,1), opd_gas(iz,1)
      enddo

      close( LUNIO )

      enddo

      return
      end
