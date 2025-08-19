      subroutine calc_optics( do_subcloud, nz, ngas, nwave, nrad,
     $  z, gas_name, radius, dr, qscat, qext, cos_qscat, 
     $  ndz, sig, rg, opd, opd_gas, w0, g0, opd_tot )

c   --------------------------------------------------------------------
c
c   Calculate spectrally-resolved profiles of optical depth, single-scattering
c   albedo, and asymmetry parameter.
c
c   input scalars:
c
c     do_subcloud  logical flag for Marley subcloud kludge 
c     nz           number of layers
c     ngas         number of condensing gases
c     nwave        number of wavelength bins
c     nrad         number of radius bins
c
c   input vectors:
c
c     z            layer altitude (cm)
c     gas_name     names of condensing gases
c     radius       radius bin centers (cm)
c     dr           width of radius bins (cm)
c     qscat        scattering efficiency
c     qext         scattering efficiency
c     cos_qscat    qscat-weighted <cos (scattering angle)>
c     ndz          number column density of condensate (cm^-3)
c     sig          geometric standard deviation of lognormal size distribution
c     rg           geometric mean radius of lognormal size distribution
c
c   output scalars:
c
c     opd_tot      total optical depth for geometric conservative scatterers
c
c   output vectors:
c
c     opd          extinction optical depth due to all condensates in layer
c     opd_gas      cumulative (from top) opd by condensing vapor as
c                  geometric conservative scatterers
c     w0           single scattering albedo
c     g0           asymmetry parameter = Q_scat wtd avg of <cos theta>
c
c
c   A. Ackerman Apr-2000
c
c   --------------------------------------------------------------------

      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'


c   Declare local storage

      logical do_subcloud
      integer nz, ngas, iz, igas, nwave, nrad, iwave, irad
      integer itop, ibot, incr, ibot_subcloud, ibot_gas
      integer ibot_cloud(MAXNGAS)

      character*(*) gas_name(MAXNGAS)
      double precision z(MAXNZ)
      double precision radius(MAXNRAD,MAXNGAS)
      double precision dr(MAXNRAD,MAXNGAS)
      double precision qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision qext(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision cos_qscat(MAXNWAVE,MAXNRAD,MAXNGAS)
      double precision scat_gas(MAXNZ,MAXNWAVE,MAXNGAS)
      double precision ext_gas(MAXNZ,MAXNWAVE,MAXNGAS)
      double precision cqs_gas(MAXNZ,MAXNWAVE,MAXNGAS)
      double precision ndz(MAXNZ,MAXNGAS)
      double precision sig(MAXNZ,MAXNGAS)
      double precision rg(MAXNZ,MAXNGAS)
      double precision opd(MAXNZ,MAXNWAVE)
      double precision w0(MAXNZ,MAXNWAVE)
      double precision g0(MAXNZ,MAXNWAVE)
      double precision opd_gas(MAXNZ,MAXNGAS)
      double precision opd_layer(MAXNZ,MAXNGAS)
      double precision opd_tot, rr, r2, norm
      double precision opd_scat, opd_ext, cos_qs
      double precision rsig, pir2ndz, arg1, arg2


c   Determine indices of top and bottom layers

      if( z(2) .gt. z(1) )then
        itop = nz
        ibot = 1
        incr = -1
      else
        itop = 1
        ibot = nz
        incr = 1
      endif


c   Initialize indices of bottoms of cloud layers for subcloud kludge

      if( do_subcloud )then
        do igas = 1, ngas
          ibot_cloud(igas) = ibot
        enddo
      endif

c   --------------------------------------------------------------------
c   Loop over layers and condensing vapors

      do iz = itop, ibot, incr
        do igas = 1, ngas

          opd_layer(iz,igas) = 0.


c   Zero spectral sums

          do iwave = 1, nwave
            scat_gas(iz,iwave,igas) = 0.
            ext_gas(iz,iwave,igas) = 0.
            cqs_gas(iz,iwave,igas) = 0.
          enddo

  
c   Optical depth for conservative geometric scatterers 

          if( ndz(iz,igas) .gt. 0. )then

            r2 = rg(iz,igas)**2 * exp( 2*log( sig(iz,igas) )**2 )
            opd_layer(iz,igas) = 2.*PI*r2*ndz(iz,igas)

c  Calculate normalization factor (forces lognormal sum = 1.0)

            rsig = sig(iz,igas)
            norm = 0.

            do irad = 1,nrad
              rr = radius(irad,igas)
              arg1 = dr(irad,igas) / ( sqrt(2.*PI)*rr*log(rsig) )
              arg2 = -log( rr/rg(iz,igas) )**2 / ( 2*log(rsig)**2 )
              norm = norm + arg1*exp( arg2 )
            enddo

            norm = ndz(iz,igas) / norm


c   Loop over wavelength and radius

            do irad = 1, nrad

              rr = radius(irad,igas)
              arg1 = dr(irad,igas) / ( sqrt(2.*PI)*log(rsig) )
              arg2 = -log( rr/rg(iz,igas) )**2 / ( 2*log(rsig)**2 )
              pir2ndz = norm*PI*rr*arg1*exp( arg2 )

              do iwave = 1, nwave
                scat_gas(iz,iwave,igas) = scat_gas(iz,iwave,igas) + 
     $            qscat(iwave,irad,igas)*pir2ndz
                ext_gas(iz,iwave,igas) = ext_gas(iz,iwave,igas) + 
     $            qext(iwave,irad,igas)*pir2ndz
                cqs_gas(iz,iwave,igas) = cqs_gas(iz,iwave,igas) + 
     $            cos_qscat(iwave,irad,igas)*pir2ndz
              enddo

            enddo   ! irad:1,nrad


c   index of bottom of cloud layer for subcloud kludge

            if( do_subcloud )then
              ibot_cloud(igas) = iz
            endif

          endif     ! ndz > 0
        enddo       ! igas:1,ngas
      enddo         ! iz:itop,ibot,incr

c   --------------------------------------------------------------------
c   subcloud kludge to soften discontinuity at cloud base
c   (10% in first layer below, 5% in second layer)

      if( do_subcloud )then
        do igas = 1, ngas
          ibot_gas = ibot_cloud(igas)
          if( ibot_gas .ne. ibot )then


c    Choose lower index to be within grid

            if( ibot_gas+incr .eq. ibot )then
              ibot_subcloud = ibot_gas + incr
            else
              ibot_subcloud = ibot_gas + 2*incr
            endif


            do iz = ibot_gas+incr, ibot_subcloud, incr

              if( iz .eq. ibot_gas+incr )then
                norm = 0.10
              else
                norm = 0.05
              endif

              opd_layer(iz,igas) = opd_layer(iz,igas) + 
     $          opd_layer(ibot_gas,igas)*norm

              do iwave = 1, nwave
                scat_gas(iz,iwave,igas) = scat_gas(iz,iwave,igas) + 
     $            scat_gas(ibot_gas,iwave,igas)*norm
                ext_gas(iz,iwave,igas) = ext_gas(iz,iwave,igas) + 
     $            ext_gas(ibot_gas,iwave,igas)*norm
                cqs_gas(iz,iwave,igas) = cqs_gas(iz,iwave,igas) + 
     $            cqs_gas(ibot_gas,iwave,igas)*norm
              enddo

            enddo
          endif
        enddo 
      endif 

c   --------------------------------------------------------------------
c   Sum over gases and compute spectral optical depth profile etc

      do iz = itop, ibot, incr
        do iwave = 1, nwave

          opd_scat = 0.
          opd_ext = 0.
          cos_qs = 0.

          do igas = 1, ngas
            opd_scat = opd_scat + scat_gas(iz,iwave,igas)
            opd_ext = opd_ext + ext_gas(iz,iwave,igas)
            cos_qs = cos_qs + cqs_gas(iz,iwave,igas)
          enddo

          if( opd_scat .gt. 0. )then
            opd(iz,iwave) = opd_ext
            w0(iz,iwave) = opd_scat / opd_ext
            g0(iz,iwave) = cos_qs / opd_scat
          else 
            opd(iz,iwave) = 0.
            w0(iz,iwave) = 0.
            g0(iz,iwave) = 0.
          endif

        enddo
      enddo


c   cumulative optical depths for conservative geometric scatterers

      opd_tot = 0.
      do igas = 1, ngas
        opd_gas(itop,igas) = opd_layer(itop,igas)
        do iz = itop+incr, ibot, incr
          opd_gas(iz,igas) = opd_gas(iz-incr,igas) + opd_layer(iz,igas)
        enddo
        opd_tot = opd_tot + opd_gas(ibot,igas)
      enddo

      return
      end
