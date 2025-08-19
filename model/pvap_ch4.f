      double precision function pvap_ch4( t )
c
c     calculate saturation vapor pressure (dyne/cm^2) over CH4
c
c     input temperature in K
c
c     Rages Jul-2000 
c
      implicit none

      integer ic
      double precision a(2),b(2),c(2)

c     quantities to calculate methane vapor pressure     
c     (WHAT ARE THE UNITS?)
c
c     AMR   -- molecular weight / ideal gas constant
c     TCRIT -- triple point temperature
c     PCRIT --    "     "   pressure
c     AS    -- specific heat at constant pressure ( gas - solid )
c     AL    --    "      "   "     "        "     ( gas - liquid )
c     ALS   -- latent heat of sublimation
c     ALV   --   "     "   "  vaporization
c
      double precision t,amr,tcrit,pcrit,as,al,als,alv
      parameter( AMR = 16.043 / 8.3143, TCRIT = 90.68, PCRIT = .11719 )
      parameter( AS = 2.213 - 2.650, AL = 2.213 - 3.370 )
      parameter( ALS = 611.10, ALV = 552.36 )
c
c     ic=1: temperature below triple point
c     ic=2: temperature above triple point
c
      ic = 1
      if (t.gt.tcrit) ic = 2

      C(1) = - AMR * AS
      C(2) = - AMR * AL
      B(1) = - AMR * ( ALS + AS * TCRIT )
      B(2) = - AMR * ( ALV + AL * TCRIT )
      A(1) = PCRIT * TCRIT ** ( -C(1) ) * EXP( -B(1) / TCRIT )
      A(2) = PCRIT * TCRIT ** ( -C(2) ) * EXP( -B(2) / TCRIT )

      pvap_ch4 = A(IC) * t**C(IC) * EXP( B(IC) / t )
   
      pvap_ch4= pvap_ch4*1e6    ! convert from bars to dyne/cm^2
c     pvap_ch4= pvap_ch4*1e1    ! convert from Pa to dyne/cm^2: T ~ 240 K
c     pvap_ch4= pvap_ch4*1e3    ! convert from mb to dyne/cm^2
c     pvap_ch4= pvap_ch4*1e2    ! nonsense, but gives cloud base T ~ 140 K

      return
      end
