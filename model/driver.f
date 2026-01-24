      program driver
      implicit none

      include 'globals.h'

c     declare variables
      double precision grav, teff, mw_atmos
      integer nz

      double precision t_top(MAXNZ+1), p_top(MAXNZ+1)
      double precision qt(MAXNZ,MAXNGAS), qc(MAXNZ,MAXNGAS)
      double precision wave(MAXNWAVE), opd(MAXNZ,MAXNWAVE)
      double precision opd_gas(MAXNZ,MAXNGAS)
      double precision w0(MAXNZ,MAXNWAVE), g0(MAXNZ,MAXNWAVE)
      double precision gas_mmr(MAXNGAS), gas_mw(MAXNGAS)
      double precision kz(MAXNZ)

      integer i, j

c     initialize parameters
      grav = 1.0d3
      teff = 1000.d0
      mw_atmos = 2.3d0
      nz = 0

c     initialize arrays
      do 10 i = 1, MAXNZ+1
         t_top(i) = 0.d0
         p_top(i) = 0.d0
   10 continue

      do 20 i = 1, MAXNZ
         kz(i) = 0.d0
         do 21 j = 1, MAXNGAS
            qt(i,j) = 0.d0
            qc(i,j) = 0.d0
            opd_gas(i,j) = 0.d0
   21    continue
         do 22 j = 1, MAXNWAVE
            opd(i,j) = 0.d0
            w0(i,j) = 0.d0
            g0(i,j) = 0.d0
   22    continue
   20 continue

      do 30 i = 1, MAXNWAVE
         wave(i) = 0.d0
   30 continue

      do 40 i = 1, MAXNGAS
         gas_mmr(i) = 0.d0
         gas_mw(i)  = 0.d0
   40 continue

c     call cloud model subroutine
      call main( grav, teff, nz, t_top, p_top, qt, qc, wave, opd,
     &           opd_gas, w0, g0, gas_mmr, gas_mw, mw_atmos, kz )

      print *, 'EddySed cloud model finished successfully.'

      end
