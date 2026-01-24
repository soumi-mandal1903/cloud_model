
      subroutine main( grav, 
     $   teff, nz, t_top, p_top, qt, qc, wave, opd, opd_gas, w0, g0, 
     $   gas_mmr, gas_mw, mw_atmos, kz)
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


c   Declare externals

c      external pvap_ch4, pvap_nh3, pvap_h2o
c      external pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3

c      double precision pvap_ch4, pvap_nh3, pvap_h2o
c      double precision pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3

      double precision pvap_ch4, pvap_nh3, pvap_h2o
      double precision pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3, 
     $      pvap_mg2sio4
      double precision pvap_mns, pvap_zns, pvap_na2s, pvap_cr

      external pvap_ch4, pvap_nh3, pvap_h2o
      external pvap_fe, pvap_kcl, pvap_mgsio3, pvap_al2o3, pvap_mg2sio4
      external pvap_mns, pvap_zns, pvap_na2s, pvap_cr


c   Declare local storage

      logical do_ttyo, do_fileo, do_subcloud, do_optics, read_mie
      logical do_cases, do_virtual

      integer NCASE, ncases, icase
      parameter( NCASE = 20 )

      integer nz, ngas, iz, igas, itop, ibot, incr
      integer nwave, nrad, iz_max, nsub_max, status_t

      character*10 gas_name(MAXNGAS)
      character*80 input_fname

      double precision z(MAXNZ)
      double precision p(MAXNZ)
      double precision t(MAXNZ)
      double precision z_top(MAXNZ+1)
      double precision p_top(MAXNZ+1)
      double precision t_top(MAXNZ+1)
      double precision chf(MAXNZ)
      double precision lapse_ratio(MAXNZ)
      double precision kz(MAXNZ)
C     double precision gas_mmr(MAXNGAS)
      double precision gas_mmr(MAXNGAS)
      double precision gas_mw(MAXNGAS)
      double precision qt(MAXNZ,MAXNGAS)
      double precision qc(MAXNZ,MAXNGAS)
      double precision ndz(MAXNZ,MAXNGAS)
      double precision sig(MAXNZ,MAXNGAS)
      double precision rainf(MAXNZ,MAXNGAS)
      double precision cloudf(MAXNZ)
      double precision rho_p(MAXNGAS)
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
      double precision dewp(MAXNZ,MAXNGAS)
      double precision rainf_case(NCASE)
      double precision reff_qc(NCASE,MAXNGAS)
      double precision reff_lwp(NCASE,MAXNGAS)
      double precision totpath_gas(NCASE,MAXNGAS)
      double precision opd_case(NCASE,MAXNGAS)
      double precision grav, teff, opd_tot, mw_atmos, kz_min, cloudf_min
      double precision sig_all, rainf_all, rhoc, rhoc_max, lwp
      double precision lwp_layer, opd_eff, supsat
      double precision totpath, tlo, thi, delta_p, rhs, eps
      
      !! Some Print statements to debug ATMO Coupling if necessary
C     print *, gas_mmr(1)
C     print *, grav
C     print *, teff
C     print *, nz
C     print *, t_top
      print *, p_top
C     print *, mw_atmos
C     stop
c   --------------------------------------------------------------------
c   Define molecular weight of atmosphere (g/mol)
      !mw_atmos = 2.2
c   Define molecular weight of atmosphere (g/mol)
      mw_atmos = 2.2d0
c   Define number of condensing gases and parameters specific to each:
c     <gas_name> is name of gas
c     <gas_mw> is molecular weight of gas (g/mol)
c     <gas_mmr> is total mass mixing ratio below cloud (g/g)
c     <rho_p> is density of condensed phase (g/cm^3)

c    ngas = 1
c    gas_name(1) = 'NH3'
c   ---------------------------------------------------------
c   Define condensing gases
c   ---------------------------------------------------------

      ngas = 1
      gas_name(1) = 'NH3'

c      ngas = 5
c      gas_name(1) = 'H2O'
c      gas_name(2) = 'Fe'
c      gas_name(3) = 'MgSiO3'
c      gas_name(4) = 'Al2O3'
c      gas_name(5) = 'NH3'

c     ngas = 7
c     gas_name(1) = 'H2O'
c     gas_name(2) = 'Fe'
c     gas_name(3) = 'MgSiO3'
c     gas_name(4) = 'Al2O3'
c     gas_name(5) = 'KCl'
c     gas_name(6) = 'NH3'
c     gas_name(7) = 'CH4'

c       ngas=7
c       gas_name(1) = 'Al2O3'         
c       gas_name(2) = 'Fe'   
c       gas_name(3) = 'Na2S'  
c       gas_name(4) = 'KCl'  
c       gas_name(5) = 'Cr'
c C      gas_name(6) = 'MgSiO3' 
c       gas_name(6) = 'Mg2SiO4'
c C     gas_name(4) = 'NH3' 
c C     gas_name(6) = 'MnS'   
c       gas_name(7) = 'ZnS'  
C       
C     gas_name(9) = 'Mg2SiO4'      

      if( ngas .gt. MAXNGAS )then 
        print*, 'main(): ngas > MAXNGAS'
        stop 1
      endif
      
            do igas = 1, ngas

        if( gas_name(igas) .eq. 'CH4' )then
c         rho_p(igas) =  0.42               ! liquid
          rho_p(igas) =  0.49               ! solid

        elseif( gas_name(igas) .eq. 'NH3' )then
c         rho_p(igas) =  0.68               ! liquid
          rho_p(igas) =  0.84               ! solid

        elseif( gas_name(igas) .eq. 'H2O' )then
c         rho_p(igas) =  1.0                ! liquid
          rho_p(igas) =  0.93               ! solid (T = 213 K)

        elseif( gas_name(igas) .eq. 'Fe' )then
          rho_p(igas) =  7.875
  
        elseif( gas_name(igas) .eq. 'KCl' )then
          rho_p(igas) =  1.99

c ENSTATITE
        elseif( gas_name(igas) .eq. 'MgSiO3' )then
          rho_p(igas) =  3.192

c FORSTERITE
        elseif( gas_name(igas) .eq. 'Mg2SiO4' )then
          rho_p(igas) =  3.214
          
c CORUNDUM
        elseif( gas_name(igas) .eq. 'Al2O3' )then
          rho_p(igas) =  3.987
          
        elseif( gas_name(igas) .eq. 'Na2S' )then
          rho_p(igas) =  1.856
          
        elseif( gas_name(igas) .eq. 'MnS' )then
          rho_p(igas) =  4.0
          
        elseif( gas_name(igas) .eq. 'ZnS' )then
          rho_p(igas) =  4.04
          
c     (ZnS: values for rho-- "spahlerite")
        elseif( gas_name(igas) .eq. 'Cr' )then
          rho_p(igas) =  7.15

        else
          print*, 'main(): bad gas_name = ', gas_name(igas)
          stop 1
        endif

      enddo
c   ---------------------------------------------------------
c   Define molecular weights and deep mixing ratios
c   ---------------------------------------------------------
      do igas = 1, ngas

        if( gas_name(igas) .eq. 'NH3' )then

          gas_mw(igas) = 17.0
          gas_mmr(igas) = 1.34e-4 * ( gas_mw(igas) / mw_atmos )

        else

          gas_mw(igas)  = 0.0
          gas_mmr(igas) = 0.0

        endif

      enddo

C     do igas = 1, ngas
C
C       if( gas_name(igas) .eq. 'CH4' )then
C
C         gas_mw(igas) = 16.
Cc          gas_mmr(igas) = 5.9e-4 * (gas_mw(igas)/mw_atmos)
C         gas_mmr(igas) = 4.9e-4 * (gas_mw(igas)/mw_atmos)
Cc	Changed from 5.9e-4 to 4.9e-4 by JJF, to reflect Lodders(2003) abunds
Cc   V.G. Manzhelii and A.M. Tolkachev, Sov. Phys. Solid State 5, 2506 (1964)
Cc         rho_p(igas) =  0.42               ! liquid
C         rho_p(igas) =  0.49               ! solid
C
C       elseif( gas_name(igas) .eq. 'NH3' )then
C
C         gas_mw(igas) = 17.
C
Cc WST86 (West, Strobel, and Tomasko, 1986) solar
Cc         gas_mmr(igas) = 1.5e-4 * (gas_mw(igas)/mw_atmos)
Cc Voyager IRIS data from Kunde et al (1982) @ 0.6 bar
Cc also corroborated by Brooke et al. (1998), and very close
Cc to Carlson et al. (1993) value of 2.8e-5 @ 0.6 bar
Cc          gas_mmr(igas) = 3.0e-5 * (gas_mw(igas)/mw_atmos)
C         gas_mmr(igas) = 1.34e-4 * (gas_mw(igas)/mw_atmos) *1.d0
Cc	Changed from 3.0e-5 to 1.34e-4 by JJF, to reflect Lodders(2003) abunds
C
C
Cc   Atkins, Physical Chemistry, 6th Ed. (1998)
Cc   V.G. Manzhelii and A.M. Tolkachev, Sov. Phys. Solid State 5, 2506 (1964)
Cc         rho_p(igas) =  0.68               ! liquid
C         rho_p(igas) =  0.84               ! solid
C
C       elseif( gas_name(igas) .eq. 'H2O' )then
C
C         gas_mw(igas) = 18.
Cc          gas_mmr(igas) = 1.06e-3 * (gas_mw(igas)/mw_atmos)
C         gas_mmr(igas) = 7.54e-4 * (gas_mw(igas)/mw_atmos)  *1.d0
Cc	Changed from 1.03e-4 to 7.54e-4 by JJF, to reflect Lodders(2003) abunds	 
Cc   Pruppacher and Klett (p. 86 of 1st edition, 1978)
Cc         rho_p(igas) =  1.0                ! liquid
C         rho_p(igas) =  0.93               ! solid (T = 213 K)
C
C       elseif( gas_name(igas) .eq. 'Fe' )then
C
C         gas_mw(igas) = 55.845
C         gas_mmr(igas) = 1.30e-3
Cc   Lodders and Fegley (1998)
C         rho_p(igas) =  7.875
C 
C       elseif( gas_name(igas) .eq. 'KCl' )then
C
Cc          gas_mw(igas) = 74.5
Cc          gas_mmr(igas) = 2.2e-7 * (gas_mw(igas)/mw_atmos)
C         gas_mw(igas) = 74.5
Cc          gas_mmr(igas) = 2.2e-7 * (gas_mw(igas)/mw_atmos)
Cc   SOLAR METALLICITY (abunds tables, 900K, 1 bar)
C           gas_mmr(igas) = 2.2627E-07 * (gas_mw(igas)/mw_atmos)
Cc   10x SOLAR METALLICITY (abunds tables, 900K, 1 bar)
Cc            gas_mmr(igas) = 2.1829E-06 * (gas_mw(igas)/mw_atmos)
Cc   50x SOLAR METALLICITY (abunds tables, 900K, 1 bar)
Cc           gas_mmr(igas) = 8.1164E-06 * (gas_mw(igas)/mw_atmos)
Cc   source unknown
C         rho_p(igas) =  1.99
C
Cc ENSTATITE
C       elseif( gas_name(igas) .eq. 'MgSiO3' )then
C
C         gas_mw(igas) = 100.4
C         gas_mmr(igas) = 2.75e-3
Cc ADDING ENSTATITE as 2nd silicate cloud. Forsterite cloud has already formed below. 
Cc Evaluate saturation vap pressure of sio over forst at T=1900K, P=1 bar (same as forsterite mmr) to get mixing ratio
Cc log10(Psio over forst) = -9746/T +1.52 + 0.64*P_T +0.28*Fe/H
Cc   evaluated for SOLAR METALLICITY:
Cc         gas_mmr(igas) = 1.71e-5 * (gas_mw(igas)/mw_atmos)
Cc   Lodders and Fegley (1998)
C         rho_p(igas) =  3.192
C
Cc FORSTERITE
C       elseif( gas_name(igas) .eq. 'Mg2SiO4' )then
C
C         gas_mw(igas) = 140.7
Cc      OLD FORSTERITE
Cc           gas_mmr(igas) = 1.38e-3
Cc      NEW FORSTERITE (from Lodders et al. table, 1000mbar, 1900K)
C         gas_mmr(igas) = 7.1625e-05/2 * (gas_mw(igas)/mw_atmos)  
Cc   Lodders and Fegley (1998)
C         rho_p(igas) =  3.214
C
C
Cc CORUNDUM
C       elseif( gas_name(igas) .eq. 'Al2O3' )then
C
C         gas_mw(igas) = 101.961
C         gas_mmr(igas) = 2.51e-6 * (gas_mw(igas)/mw_atmos)
Cc   Lodders and Fegley (1998)
C         rho_p(igas) =  3.987
C
Cc Na2S
C       elseif( gas_name(igas) .eq. 'Na2S' )then
C
C         gas_mw(igas) = 78.05
Cc   Solar met.
C            gas_mmr(igas) = 1.97e-6 * (gas_mw(igas)/mw_atmos) 
Cc   10x solar met. 
Cc           gas_mmr(igas) =  0.5*3.6112E-05  * (gas_mw(igas)/mw_atmos) 
Cc   50x solar met. 
Cc             gas_mmr(igas) = 0.5*1.1041E-04 * (gas_mw(igas)/mw_atmos) 
Cc   Lodders and Fegley (2003) (cvm)
C         rho_p(igas) =  1.856
C
Cc MnS
C       elseif( gas_name(igas) .eq. 'MnS' )then
C
C         gas_mw(igas) = 87.00
Cc   Solar met.
C           gas_mmr(igas) = 6.37e-7 * (gas_mw(igas)/mw_atmos) 
Cc   10x solar met. 
Cc            gas_mmr(igas) =  10*6.37e-7  * (gas_mw(igas)/mw_atmos) 
Cc   50x solar met. 
Cc            gas_mmr(igas) = 50*6.37e-7 * (gas_mw(igas)/mw_atmos) 
Cc   Lodders and Fegley (2003) (cvm)
C         rho_p(igas) =  4.0
C
Cc ZnS
C       elseif( gas_name(igas) .eq. 'ZnS' )then
C
C         gas_mw(igas) = 97.46
Cc Solar met. 
C             gas_mmr(igas) = 8.40e-8 * (gas_mw(igas)/mw_atmos) 
Cc 10x solar met. 
Cc           gas_mmr(igas) = 10*8.40e-8 * (gas_mw(igas)/mw_atmos) 
Cc 50x solar met. 
Cc           gas_mmr(igas) = 50*8.40e-8 * (gas_mw(igas)/mw_atmos) 
C
Cc   Lodders and Fegley (2003) (cvm)
C         rho_p(igas) =  4.04
Cc     (ZnS: values for rho-- "spahlerite")
C
C
Cc CHROMIUM
C       elseif( gas_name(igas) .eq. 'Cr' )then
C
C         gas_mw(igas) = 51.996
Cc Solar met. 
C          gas_mmr(igas) = 8.80e-7 * (gas_mw(igas)/mw_atmos) 
Cc 10x solar met. 
Cc          gas_mmr(igas) = 8.6803E-06 * (gas_mw(igas)/mw_atmos) 
Cc 50x solar met. 
Cc          gas_mmr(igas) = 4.1308E-05 * (gas_mw(igas)/mw_atmos) 
Cc   Lodders and Fegley (2003) (cvm)
C         rho_p(igas) =  7.15
C
C       else
C         print*, 'main(): bad gas_name = ', gas_name(igas)
C         stop 1
C       endif
C
C     enddo

c   --------------------------------------------------------------------
c   User-defined model parameters for dynamics and microphysics

c   minimum eddy diffusion coefficient (cm^2/s)
      kz_min = 1e5      ! baseline

c   geometric standard deviation for lognormal size distributions
c     sig_all = 2.0     ! baseline
      sig_all = 2.0

c   rain factor 
c     rainf_all = 3.    ! baseline
c      rainf_all = 0.1
      rainf_all = 3     ! f_sed = 1 (later try 3, 10)
      cloudf_min = 1

c   minimum cloud coverage (a diagnostic not applied to the clouds)
c     cloudf_min = 1.   ! old baseline
c      cloudf_min = 0.75

c   minimum cloud coverage
c     nsub_max = 64     ! baseline
      nsub_max = 64

c   include decrease in condensate mixing ratio below model domain
      do_virtual = .true. !! was true but switched off while chemistry coupling testing
 
c   ramp up optical depth below cloud base in calc_optics()
      do_subcloud = .true.

c   saturation factor (after condensation)
      supsat = 0

c   optionally (when do_cases = .true.) calculate cloud profiles for
c   a number of different values of rainf; otherwise, just calculate
c   profiles for one value of rainf

      do_cases = .false.
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
    
      do_optics = .true.
      read_mie = .true.

       call init_optics( do_optics, read_mie, ngas, gas_name,
     $   nwave, wave, nrad, radius, dr, qscat, qext, cos_qscat )

c   --------------------------------------------------------------------
c   Read the atmospheric profiles
c
c   Voyager profile used in Ackerman and Marley, 2001
c
c     input_fname='input/profiles/voyager.input'

c   irradiated type IV roaster (Sudarsky et al 2000 according to Tim Brown)
c     input_fname='input/profiles/roaster1400.input'

c     call read_voyager( input_fname, 
c    $     grav, teff, nz, z, z_top, p, p_top, t, t_top, chf )

      input_fname = 'input/profiles/voyager.input'

      nz = 0
      call read_voyager( input_fname,
     $     grav, teff, nz, z, z_top, p, p_top, t, t_top, chf )

c     input_fname='t500g300f5.txt'
c     input_fname='input/profiles/t2000_g1e5.input'
c      input_fname='input/profiles/t1500_g1e5.input'
c     input_fname='input/profiles/t1000_g1e5.input'
c     input_fname='input/profiles/t900_g1e5.input'
c     input_fname='input/profiles/t500_g1e4.input'
c     input_fname='input/profiles/t128_g1e4.input'

      !input_fname='input/profiles/hd189_day_line.input'

C     nz=0
C     call read_lewis( input_fname, 
C    $     grav, teff, nz, z, z_top, p, p_top, t, t_top, chf )
c      call read_goyal( grav, 
c     $     teff, nz, z, z_top, p, p_top, t, t_top, chf )

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
     $  grav, teff, kz_min, cloudf_min, nsub_max, 
     $  supsat, mw_atmos, do_virtual,
     $  nz, z, z_top, p, p_top, t, t_top, chf, 
     $  ngas, gas_name, gas_mmr, gas_mw, rho_p, sig, rainf,
     $  kz, qt, qc, ndz, rg, reff, cloudf )

C     print *, "jm1 = ", qt(:,1)
C     stop
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
 
Cc  Zero diagnostics 
C
        rhoc_max = 0.
        lwp = 0.
        totpath = 0.
        opd_eff = 0.
 
        do iz = itop, ibot, incr
C
Cc  Increment vertical diagnostics
 
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
 
Cc   --------------------------------------------------------------------------
Cc   Dewpoint calculation
 
Cc   ratio of molecular weights
          eps = gas_mw(igas) / mw_atmos
 
Cc   right hand side in root expression:  pvap(t_dew) = p*q/(eps+q) [dyne/cm^2]
          rhs = p(iz)*gas_mmr(igas) / ( eps + gas_mmr(igas) )
 
Cc   delta_p is precision of vapor pressure root [dyne/cm^2]
          delta_p = rhs / 100.
 
Cc   tlo and thi are lower and upper guesses of dewpoint [K]
 
          if( gas_name(igas) .eq. 'CH4' )then
            tlo = 10.
            thi = 200.
            call find_rootl( pvap_ch4, rhs, tlo, thi, delta_p,
     $           dewp(iz,igas), status_t )
          elseif( gas_name(igas) .eq. 'NH3' )then
            tlo = 100.
            thi = 300.
            call find_rootl( pvap_nh3, rhs, tlo, thi, delta_p,
     $           dewp(iz,igas), status_t )
          elseif( gas_name(igas) .eq. 'H2O' )then
            tlo = 100.
            thi = 400.
            call find_rootl( pvap_h2o, rhs, tlo, thi, delta_p,
     $           dewp(iz,igas), status_t )
          elseif( gas_name(igas) .eq. 'Fe' )then
            tlo = 1000.
            thi = 3000.
            call find_rootl( pvap_fe, rhs, tlo, thi, delta_p,
     $           dewp(iz,igas), status_t )
          elseif( gas_name(igas) .eq. 'KCl' )then
            tlo = 1000.
            thi = 2000.
            call find_rootl( pvap_kcl, rhs, tlo, thi, delta_p,
     $           dewp(iz,igas), status_t )
          elseif( gas_name(igas) .eq. 'MgSiO3' )then
            tlo = 1000.
            thi = 3000.
            call find_rootl( pvap_mgsio3, rhs, tlo, thi, delta_p,
     $           dewp(iz,igas), status_t )
          elseif( gas_name(igas) .eq. 'Al2O3' )then
            tlo = 1000.
            thi = 4000.
            call find_rootl( pvap_al2o3, rhs, tlo, thi, delta_p,
     $           dewp(iz,igas), status_t )
          else
Cc            print*,'main(): for dew point need to add'
Cc            print*,' find_rootl(', gas_name(igas), ')'
            dewp(iz,igas) = 0.
          endif
 
          if( status_t .ne. 0 )then
Cc            print*,'main(): dew point not found'
Cc            write(*,'(a,i3,3a,i4,1pe10.2)')
Cc     $        ' find_rootl(pvap) status = ',status_t,
Cc     $        ' for ',gas_name(igas),' at iz,p = ', iz,p(iz)/1e6
          endif
 
        enddo

Cc  Cloud diagnostics
C
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
C
Cc   --------------------------------------------------------------------
Cc   Select output formats
C
      if( do_cases )then
        do_ttyo = .false.
      else
        do_ttyo = .true.
      endif
 
      do_fileo = .true.
 
Cc   --------------------------------------------------------------------
Cc   terminal output
C
      if( do_ttyo )then
 
        print*, 'main():'
        print*, ''
        print*, 'Teff = ', int(teff+0.5), ' K'
        print*, 'grav = ', int(grav/1e2+0.5), ' m/s^2'
        print*, 'nz = ', int(nz)
 
        do igas = 1, ngas
 
          print*, ''
          print*, 'condensing gas = ', gas_name(igas)
          print*, ''
          print*, 'iz  P(bar)   z(m)     T(K) kz(cm^2/s) ' 
     $      // ' reff(um) qc(g/g) col_opd cloudf dewp(K)'     
 
          do iz = itop, ibot, incr
            write( *,
     $        '( i3, 1pe9.1, 1pe9.2, 0pf7.0, 1p5e9.1, 0pf7.0)' )
     $        iz, p(iz)/1e6, z(iz)/1e2, t(iz),
     $        kz(iz),
     $        reff(iz,igas)*1e4, qc(iz,igas), opd_gas(iz,igas),
     $        cloudf(iz), dewp(iz,igas)
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
 
Cc   --------------------------------------------------------------------
Cc   file output
C
      if( do_fileo )then
 
        open( LUNIO, file='eddysed.out', status='replace' )
 
        write( LUNIO, * ) teff, grav, nz, ngas, sig_all, rainf_all
      
        do igas = 1, ngas
          write( LUNIO, * ) gas_name(igas)
          write( LUNIO, * ) gas_mmr(igas), opd_gas(ibot,igas)
        enddo
 
        do iz = 1, nz
          write( LUNIO, '(1p27e11.3)' ) z(iz), p(iz), t(iz), 
     $      ( qt(iz,igas), igas=1,ngas ),
     $      ( qc(iz,igas), igas=1,ngas ),
     $      ( reff(iz,igas), igas=1,ngas ),
     $      ( rg(iz,igas), igas=1,ngas ),
     $      ( opd_gas(iz,igas), igas=1,ngas ),
     $      ( dewp(iz,igas), igas=1,ngas )
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
