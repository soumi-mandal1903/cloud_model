      program main
      
      implicit none

      integer irad
      double precision grav, mw_atmos, mfp, visc, rad(4), r
      double precision t_layer, p_layer, rho_p, mw_cloud

      data rad / 10., 100., 300., 1000. /

      common / vfall_block /
     $  grav, mw_atmos, mfp, visc, t_layer, p_layer, rho_p, mw_cloud

      double precision vfall_save, vfall_rossow, vfall
      external vfall_save, vfall_rossow, vfall

      grav = 980.
      mw_atmos = 29.
      mfp = 6.2e-6
      visc = 1.6e-4
      t_layer = 288.
      p_layer = 1e6
      rho_p = 1.

      do irad = 1, 4
        r = rad(irad)*1e-4
        print*,'r, vfall_old, rossow, new=',
     $    r*1e4,
     $    vfall_save(r),
     $    vfall_rossow(r),
     $    vfall(r)
       enddo

       end
