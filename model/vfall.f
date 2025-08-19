      double precision function vfall( r )
c
c   calculate fallspeed for a spherical particle at one layer in an
c   atmosphere, depending on Reynolds number for Stokes flow.
c
c   For Re_Stokes < 1, use Stokes velocity with slip correction
c   For Re_Stokes > 1, use fit to Re = exp( b1*x + b2*x^2 )
c     where x = log( Cd Re^2 / 24 )
c     where b2 = -0.1 (curvature term) and 
c     b1 from fit between Stokes at Re=1, Cd=24 and Re=1e3, Cd=0.45
c
c   and Precipitation, Reidel, Holland, 1978) and Carlson, Rossow, and
c   Orton (J. Atmos. Sci. 45, p. 2066, 1988)
c
c   all units are cgs
c
c   input parameters passed through argument list:
c
c     r         particle radius (cm)
c
c   input parameters passed through common vfall_block
c
c     vf_grav        acceleration of gravity (cm/s^2)
c     vf_mw_atmos    atmospheric molecular weight (g/mol)
c     vf_mfp         atmospheric molecular mean free path (cm)
c     vf_visc        atmospheric dynamic viscosity (dyne s/cm^2)
c     vf_t           atmospheric temperature (K)
c     vf_p           atmospheric pressure (dyne/cm^2)
c     vf_rhop        density of particle (g/cm^3)
c
c   A. Ackerman Feb-2000
c
      implicit none


c   declare common storage 

      double precision vf_grav, vf_mw_atmos, vf_mfp, vf_visc
      double precision vf_t, vf_p, vf_rhop

      common / vfall_block /
     $  vf_grav, vf_mw_atmos, vf_mfp, vf_visc, 
     $  vf_t, vf_p, vf_rhop

c   declare and define local storage

      double precision r, knudsen, rho_atmos, drho
      double precision reynolds, slip, x, y, b1, b2, cdrag

      parameter( b1 = 0.8 )            ! Ackerman
c     parameter( b1 = 0.86 )           ! Rossow
c     parameter( b1 = 0.72 )           ! Carlson
      parameter( b2 = -0.01 ) 
      parameter( cdrag = 0.45 )        ! Ackerman
c     parameter( cdrag = 0.2 )         ! Rossow
c     parameter( cdrag = 2.0 )         ! Carlson


c   universal gas constant (erg/mol/K)
      double precision R_GAS
      parameter( R_GAS = 8.3143e7 )


c   calculate vfall based on Knudsen and Reynolds numbers

      knudsen = vf_mfp / r
      rho_atmos = vf_p / ( (R_GAS/vf_mw_atmos) * vf_t )
      drho = vf_rhop - rho_atmos

c   Cunningham correction (slip factor for gas kinetic effects)
      slip = 1. + 1.26*knudsen

c   Stokes terminal velocity (low Reynolds number)
      vfall = slip*(2./9.)*drho*vf_grav*r**2 / vf_visc
      reynolds = 2.*r*rho_atmos*vfall / vf_visc

      if( reynolds .gt. 1. )then

c   correct drag coefficient for turbulence (Re = Cd Re^2 / 24)

        x = log( reynolds )
        y = b1*x + b2*x**2
        reynolds = exp(y)
        vfall = vf_visc*reynolds / (2.*r*rho_atmos)

        if( reynolds .gt. 1e3 )then

c   drag coefficient independent of Reynolds number

          vfall = slip*sqrt( 8.*drho*r*vf_grav / (3.*cdrag*rho_atmos) )

        endif

      endif

      return
      end
