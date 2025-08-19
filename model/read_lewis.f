      subroutine read_lewis( input_file, 
     $   grav, teff, nz, z, z_top, p, p_top, t, t_top, chf )
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
      character*(*) input_file
      character*1 dumstr


c   Open input file for reading atmospheric profiles
      open( LUNIO, file=input_file, status='old' )

c   Read effective temperature (teff, in K), 
c   gravitational acceleration (grav, cm/s^2),
c   and number latitude (nx), longitude (ny), and layer edges (nzp1)

      read(LUNIO,*) teff
      read(LUNIO,*) grav
      read(LUNIO,*) nzp1

c   Skip blank line

      nskip = 1
      do iz = 1, nskip
        read(LUNIO,*) dumstr
      enddo

c   Read the data -- column 2 with zeros means column 3 has pressures
c   (convert pressure in bars to dyne/cm^2)

      if( nzp1-1 .gt. MAXNZ )then 
        print*, 'read_lewis(): nz > MAXNZ'
        stop 1
      endif
 
      do iz = 1, nzp1
        read(LUNIO,*) idummy, dum1, dum2, dum3
        t_top(iz) = dum2
        if( dum1 .eq. 0. )then
          p_top(iz) = dum3*1e6
        else
          p_top(iz) = dum1*1e6
        endif
      enddo

      close(LUNIO)

c   Get pressure thickness and P, T at pressure mid-pt of layer
c   and hydrostatically calculate altitudes

c   should pass atmospheric molecular weight through argument list
      r_atmos = 8.3143e7 / 2.2
      
      z_top(nz+1) = 0.
      
      nz = nzp1 - 1
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
