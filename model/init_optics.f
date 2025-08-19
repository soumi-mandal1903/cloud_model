      subroutine init_optics( do_optics, read_mie, ngas, gas_name, 
     $  nwave, wave, nrad, radius, dr, qscat, qext, cos_qscat )

c   --------------------------------------------------------------------
c
c   Setup up a particle size grid and calculate single-particle scattering
c   and absorption efficiencies and other parameters to be used by
c   calc_optics()
c
c   input scalars:
c
c     do_optics    .false. means do nothing
c     read_mie     .true. means read Mie coefficients from file
c                  'gas_name'.mieff
c     ngas         number of condensing gases
c
c   input vectors:
c
c     gas_name     names of condensing gases
c
c
c ( intermediate vectors defined by this routine:
c
c     rup          upper bounds of radius bins (cm)
c     rmin         minimum radius of grid (cm)
c )
c
c   output vectors:
c
c     wave         wavelength bin centers (cm)
c     radius       radius bin centers (cm)
c     dr           widths of radius bins (cm)
c     qscat        scattering efficiency
c     qext         scattering efficiency
c     cos_qscat    qscat * average <cos (scattering angle)> 
c
c
c   A. Ackerman Feb-2000
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'


c   Declare local storage

      logical okay, read_mie, bad_mie, do_optics
      integer ngas, igas, nwave, iwave, nrad, irad
      integer mwave, mrad, iskip, isub, istatus
      integer n_thetd, nskip, ns, idummy
      character*(*) gas_name(ngas)
      character*(80) this_gas, filename
      character*1 dumstr
      double precision wave(MAXNWAVE)
      double precision radius(MAXNRAD,MAXNGAS)
      double precision dr(MAXNRAD,MAXNGAS)
      double precision rup(MAXNRAD,MAXNGAS)
      double precision qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision qext(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision cos_qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision rmin(MAXNGAS)
      double precision wvno, qs, qe, c_qs, corerad, corereal, coreimag
      double precision vrat, pw, f1, f2, r_in, wave_in
      double precision dr5, rr, nn, kk, thetd
      double precision qs_pass, qe_pass, c_qs_pass

c   --------------------------------------------------------------------

c   Report input controls

      print*,'init_optics(): do_optics, read_mie = ',do_optics,read_mie

c   If do_optics is false, zero nrad and nwave and return
c   (geometric scatterers only)

      if( .not. do_optics )then
        nrad = 0
        nwave = 0
        return
      endif

c   Define radius grids (could define separately for each species)

      nrad = 40
      if( nrad .gt. MAXNRAD )then
        print*, 'init_optics(): nrad > MAXNRAD'
        stop 1
      endif

      do igas = 1, ngas

        vrat = 2.2
        rmin(igas) = 1e-5

c       if( gas_name(igas) .eq. 'CH4' )then
c       elseif( gas_name(igas) .eq. 'NH3' )then
c       elseif( gas_name(igas) .eq. 'H2O' )then
c       elseif( gas_name(igas) .eq. 'Fe' )then
c       elseif( gas_name(igas) .eq. 'KCl' )then
c       elseif( gas_name(igas) .eq. 'MgSiO3' )then
c       elseif( gas_name(igas) .eq. 'Al2O3' )then
c       endif

        pw = 1. / 3.
        f1 = ( 2*vrat / ( 1 + vrat) )**pw
        f2 = ( 2 / ( 1 + vrat ) )**pw * (vrat**pw-1) !! Bug corrected after suggestion from Mark Marley

        do irad = 1, nrad
          radius(irad,igas) = rmin(igas) * vrat**(float(irad-1)/3.)
          rup(irad,igas) = f1*radius(irad,igas)
          dr(irad,igas) = f2*radius(irad,igas)
        enddo

c hack: Banfield et al
c       radius(1,igas) = 0.2e-4
c       dr(1,igas) = f2*radius(1,igas)
c       rup(1,igas) = f1*radius(1,igas)

c hack: Brooke et al
c       radius(1,igas) = 1e-4
c       dr(1,igas) = f2*radius(1,igas)
c       rup(1,igas) = f1*radius(1,igas)
c       radius(2,igas) = 10e-4
c       dr(2,igas) = f2*radius(2,igas)
c       rup(2,igas) = f1*radius(2,igas)

      enddo

c   --------------------------------------------------------------------
c   Define optical properties

      nwave = 196
c hack
c     nwave = 1

      if( nwave .gt. MAXNWAVE )then
        print*, 'init_optics(): nwave > MAXNWAVE'
        stop 1
      endif

      print*, 'init_optics(): nwave, nrad = ', nwave, nrad

c   --------------------------------------------------------------------
c   Read extinction and scattering coefficients for each condensing vapor

      if( read_mie )then

        do igas = 1, ngas
 
          this_gas = gas_name(igas)
          call dblank( this_gas, ns )
          filename = 'input/optics/' // this_gas(1:ns) // '.mieff'
          open( LUNIO, file=filename, form='formatted', status='old' )
 

c   Check that input file is consistent with radius and wavelength grids
 
          read( LUNIO, * ) mwave, mrad
          okay = (mwave .eq. nwave) .and. (mrad .eq. nrad) 
   
          if ( .not. okay )then
            print*,'init_optics(): input grid bad ',
     $        'for gas(', igas, ') = ', this_gas(1:ns), ':'
            print*,' mwave, mrad = ', mwave, mrad
            stop 1
          endif
 

          do irad = 1, nrad


c   Read and check input radii

            read( LUNIO, * ) r_in
            okay = abs( 1 - radius(irad,igas)/r_in ) .lt. 1d-6
            if ( .not. okay )then
              print*,'init_optics(): input radius grid bad ',
     $          'for gas(', igas, ') = ', this_gas(1:ns), ':'
              print*,' irad, radius, r_in = ', 
     $          irad, radius(irad,igas), r_in
              stop 1
            endif


c   Read and check wavelength and scattering efficiencies etc.

            do iwave = 1, nwave

              read( LUNIO ,* ) wave_in, qscat(iwave,irad,igas),
     $          qext(iwave,irad,igas), cos_qscat(iwave,irad,igas)

              if( igas .eq. 1 )then

                wave(iwave) = wave_in

              else

                okay = abs( 1 - wave(iwave)/wave_in ) .lt. 1d-7
                if ( .not. okay )then
                  print*,'init_optics(): input wavelength grid bad ',
     $              'for gas(', igas, ') = ', this_gas(1:ns), ':'
                  print*,' iwave, wave, wave_in = ', 
     $              iwave, wave(iwave), wave_in
                  stop 1
                endif

              endif

            enddo
          enddo

          close( LUNIO )

        enddo

      else 

c   --------------------------------------------------------------------
c   Calculate single-scattering efficiencies etc from refractive indices
c   for each condensing vapor

c   Mie parameters:
c   thetd is angle between incident and scattered radiation
c   n_thetd is number of thetd values to consider

        thetd = 0.0
        n_thetd = 1

        do igas = 1, ngas

          this_gas = gas_name(igas)
          call dblank( this_gas, ns )
          filename = 'input/optics/' // this_gas(1:ns) // '.refrind'
          open( LUNIO, file=filename, form='formatted', status='old' )


c   Skip filler and header lines

          nskip = 0
          do iskip = 1, nskip
            read( LUNIO, * ) dumstr
          enddo


c   Read wavelengths (convert um to cm) and refractive indices

          do iwave = 1, nwave

            read( LUNIO, * ) idummy, wave_in, nn, kk

            if( igas .eq. 1 )then

              wave(iwave) = wave_in*1e-4
              wvno = 2*PI / wave(iwave)

            else

c   Consistency check

              okay = abs( 1 - wave(iwave)/(wave_in*1e-4) ) .lt. 1d-7
              if( .not. okay )then
                print*,'init_optics(): input wavelength grid bad ',
     $            'for gas(', igas, ') = ', gas_name(igas), ':'
                print*,' iwave, wave, wave_in = ', 
     $            iwave, wave(iwave), wave_in*1e-4
                stop 1
              endif

            endif
   
            do irad = 1, nrad

c   --------------------------------------------------------------------
c   Subdivide radius grid into 6 bins (to avg out oscillations) and
c   call Mie code

              if( irad .eq. 1 )then
                dr5 = ( rup(1,igas) - radius(1,igas) ) / 5.
                rr  = radius(1,igas)
              else
                dr5 = ( rup(irad,igas) - rup(irad-1,igas) ) / 5.
                rr  = rup(irad-1,igas)
              endif

              qext(iwave,irad,igas) = 0.
              qscat(iwave,irad,igas) = 0.
              cos_qscat(iwave,irad,igas) = 0.

              corerad = 0.
              corereal = 1.
              coreimag = 0.

c    Only want one warning (not 6) message per radius bin
              bad_mie = .false.

              do isub = 1, 6

                call mie_calc( rr, nn, kk, thetd, n_thetd, 
     $            qe_pass, qs_pass, c_qs_pass, 
     $            corerad, corereal, coreimag, wvno,
     $            istatus )

                if( istatus .eq. 0 )then
                  qe = qe_pass
                  qs = qs_pass
                  c_qs = c_qs_pass
                else
                  if( .not. bad_mie )then
                    bad_mie = .true.
                    print*,'init_optics(): no Mie solution for '
     $                // 'irad, r(um), iwave, wave(um), n, k, gas = '
                  endif
                  print*,
     $              irad,rr*1e4,iwave,wave(iwave)*1e4,nn,kk,
     $              ' ',gas_name(igas)
                endif

                qext(iwave,irad,igas) = qext(iwave,irad,igas) + qe
                qscat(iwave,irad,igas) = qscat(iwave,irad,igas) + qs
                cos_qscat(iwave,irad,igas) = 
     $            cos_qscat(iwave,irad,igas) + c_qs

                rr = rr + dr5

              enddo
              
              qext(iwave,irad,igas) = qext(iwave,irad,igas) / 6.
              qscat(iwave,irad,igas) = qscat(iwave,irad,igas) / 6.
              cos_qscat(iwave,irad,igas) = 
     $          cos_qscat(iwave,irad,igas) / 6.

            enddo
          enddo
        enddo

        close(LUNIO)


c   --------------------------------------------------------------------
c   Write extinction and scattering coefficients
 
        do igas = 1, ngas

          this_gas = gas_name(igas)
          call dblank( this_gas, ns )
          filename = 'input/optics/' // this_gas(1:ns) // '.mieff'
          open( LUNIO, file=filename, form='formatted', 
     $      status='replace' )
 
          write( LUNIO, * ) nwave, nrad
 
          do irad = 1, nrad
            write( LUNIO, * ) radius(irad,igas)
            do iwave = 1, nwave
              write( LUNIO ,* ) wave(iwave), qscat(iwave,irad,igas),
     $          qext(iwave,irad,igas), cos_qscat(iwave,irad,igas)
            enddo
          enddo
 
          close( LUNIO )

        enddo
      endif

      return
      end
