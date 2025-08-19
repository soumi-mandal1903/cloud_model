      subroutine read_goyal( grav, 
     $   teff, nz, z, z_top, p, p_top, t, t_top, chf )
c
c   Read (and define) atmospheric properties for Marley input format
c   (variables defined in calling program)
c
c   A. Ackerman Mar-2000
c
      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'


c   Declare local storage

      integer nz, iz, nskip, idummy, nzp1, itop, ibot
      double precision z(MAXNZ)
      double precision p(MAXNZ)
      double precision t(MAXNZ)
      double precision chf(MAXNZ)
      double precision z_top(MAXNZ+1)
      double precision p_top(MAXNZ+1)
      double precision t_top(MAXNZ+1)
      double precision grav, teff, r_atmos, dum1, dum2, dum3, dlnp
      double precision scale_h, dz_layer, dz_pmid, dtdz, dtdlnp
C     character*(*) input_file
      character*1 dumstr


c   Open input file for reading atmospheric profiles
C     open( LUNIO, file=input_file, status='old' )
C
Cc   Read effective temperature (teff, in K), 
Cc   gravitational acceleration (grav, cm/s^2),
Cc   and number latitude (nx), longitude (ny), and layer edges (nzp1)
C
C     read(LUNIO,*) teff
C     read(LUNIO,*) grav
C     read(LUNIO,*) nzp1
C
Cc   Skip blank line
C
C     nskip = 1
C     do iz = 1, nskip
C       read(LUNIO,*) dumstr
C     enddo
C
Cc   Read the data -- column 2 with zeros means column 3 has pressures
Cc   (convert pressure in bars to dyne/cm^2)
C
C     if( nzp1-1 .gt. MAXNZ )then 
C       print*, 'read_lewis(): nz > MAXNZ'
C       stop 1
C     endif
C
C     do iz = 1, nzp1
C       read(LUNIO,*) idummy, dum1, dum2, dum3
C       t_top(iz) = dum2
C       if( dum1 .eq. 0. )then
C         p_top(iz) = dum3*1e6
C       else
C         p_top(iz) = dum1*1e6
C       endif
C     enddo
C
C     close(LUNIO)

c   Get pressure thickness and P, T at pressure mid-pt of layer
c   and hydrostatically calculate altitudes

c   should pass atmospheric molecular weight through argument list
      r_atmos = 8.3143e7 / 2.2
      
      z_top(nz+1) = 0.
      
      nzp1 = nz + 1
      do iz = nz, 1, -1
        itop = iz
        ibot = iz + 1
        dlnp = log( p_top(ibot)/p_top(itop) )
        p(iz) = 0.5*( p_top(itop) + p_top(ibot) )
c new
        dtdlnp = ( t_top(itop) - t_top(ibot) ) / dlnp
        t(iz) = t_top(ibot) + log( p_top(ibot)/p(iz) )*dtdlnp
        scale_h = r_atmos * t(iz) / grav
        dz_pmid = scale_h * log( p_top(ibot)/p(iz) )
        dz_layer = scale_h * dlnp
c old
c       scale_h = r_atmos * t_top(ibot) / grav
c       dz_layer = scale_h * dlnp
c       dtdz = ( t_top(itop) - t_top(ibot) ) / dz_layer
c       dz_pmid = scale_h * log( p_top(ibot)/p(iz) )
c       t(iz) = t_top(ibot) + dtdz*dz_pmid

        z(iz) = z_top(ibot) + dz_pmid
        z_top(iz) = z_top(ibot) + dz_layer
c       print*,'read_marley iz,dlnp,p,t,z=',iz,dlnp,p(iz),t(iz),z(iz)
      enddo


c   convective heat flux 

      do iz = 1, nz
        chf(iz) = STEFBOLTZ * teff**4
      enddo

      return
      end
