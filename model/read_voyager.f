      subroutine read_voyager( input_file, 
     $   grav, teff, nz, z, z_top, p, p_top, t, t_top, chf )
c
c   Read (and define) atmospheric properties for input format of
c   Galileo probe sounding (variables defined in calling program)
c
c   A. Ackerman Aug-2000
c
      implicit none


c   Include common data shared with eddysed()

      include 'globals.h'


c   Declare local storage

      integer nz, iz, nskip, idummy, nzp1, nmod, itop, ibot, k
      double precision z(MAXNZ)
      double precision p(MAXNZ)
      double precision t(MAXNZ)
      double precision z_top(MAXNZ+1)
      double precision p_top(MAXNZ+1)
      double precision t_top(MAXNZ+1)
      double precision chf(MAXNZ)
      double precision grav, teff, r_atmos, t_in, p_in, dlnp
      double precision scale_h, dz_layer, dz_pmid, dtdz
      character*(*) input_file
      character*1 dumstr


c   Open input file for reading atmospheric profiles

      open( LUNIO, file=input_file, status='old' )


c   Read effective temperature (teff, in K),
c   gravitational acceleration (grav, cm/s^2),
c   and number of layer edges (nzp1)

      read(LUNIO,*) teff
      read(LUNIO,*) grav
      read(LUNIO,*) nzp1

c   Skip text column headings

      nskip = 1
      do iz = 1, nskip
        read(LUNIO,*) dumstr
      enddo


c   Read the data

      if( nzp1-1 .gt. MAXNZ )then 
        print*, 'read_voyager(): nz > MAXNZ'
        stop 1
      endif

c   number of layers to average over
      nmod = 1

      read(LUNIO,*) t_in, p_in
      t_top(1) = t_in
      p_top(1) = p_in*1e6           ! convert bars to dyne/cm^2

      do iz = 1, nzp1-1
        read(LUNIO,*) t_in, p_in
        if( mod(iz,nmod) .eq. 0 )then 
          k = iz / nmod + 1
          t_top(k) = t_in
          p_top(k) = p_in*1e6
        endif
      enddo
      nz = (nzp1 - 1) / nmod

      close(LUNIO)


c   Get pressure thickness and P, T at pressure mid-pt of layer
c   and hydrostatically calculate altitudes

c   should pass atmospheric molecular weight through argument list
      r_atmos = 8.3143e7 / 2.2

      z_top(nz+1) = 0.

      do iz = nz, 1, -1
        itop = iz 
        ibot = iz + 1
        dlnp = log( p_top(ibot)/p_top(itop) )
        p(iz) = 0.5*( p_top(itop) + p_top(ibot) )
        scale_h = r_atmos * t_top(ibot) / grav
        dz_layer = scale_h * dlnp
        z_top(iz) = z_top(ibot) + dz_layer
        dtdz = ( t_top(itop) - t_top(ibot) ) / dz_layer
        dz_pmid = scale_h * log( p_top(ibot)/p(iz) )
        z(iz) = z_top(ibot) + dz_pmid
        t(iz) = t_top(ibot) + dtdz*dz_pmid
      enddo


c   convective heat flux 

      do iz = 1, nz
        chf(iz) = STEFBOLTZ * teff**4
      enddo


      return
      end
