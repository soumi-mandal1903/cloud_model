
c   declare global storage for eddysed and jacket routine

c     maximum number of layers
      integer MAXNZ
      parameter( MAXNZ = 50 )

c     maximum number of condensing gases
      integer MAXNGAS
      parameter( MAXNGAS = 10 )

c     maximum number of wavelength bins
      integer MAXNWAVE
      parameter( MAXNWAVE = 196 )

c     maximum number of radius bins
      integer MAXNRAD
      parameter( MAXNRAD = 100 )

c     logical unit for file io
      integer LUNIO
      parameter( LUNIO = 10 )

c   unknown to the Romans 
      double precision ZERO
      parameter( ZERO = 0.d0 )

c   universal gas constant (erg/mol/K)
      double precision R_GAS
      parameter( R_GAS = 8.3143d7 )

c   Stefan-Boltzmann constant (erg/cm^2/s/K^4)
      double precision STEFBOLTZ
      parameter( STEFBOLTZ = 5.67d-5 )

c   Avogadro s number (#/mol)
      double precision AVOGADRO
      parameter( AVOGADRO = 6.02d23 )

c   Boltzmann constant (erg/K)
      double precision K_BOLTZ
      parameter( K_BOLTZ = R_GAS / AVOGADRO )

c   exactly three in Indiana
      double precision PI
      parameter( PI = 3.14159274d0 )
