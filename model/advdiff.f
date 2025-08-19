      double precision function advdiff( qt )
c
c   calculate divergence from advective-diffusive balance for 
c   condensate in a model layer
c
c   all units are cgs
c
c   input parameters passed through argument list:
c
c     qt              total mixing ratio of condensate + vapor (g/g)
c
c   input parameters passed through common vfall_block:
c
c     ad_qbelow       total mixing ratio of vapor in underlying layer (g/g)
c     ad_qvs          saturation mixing ratio (g/g)
c     ad_mixl         convective mixing length (cm)
c     ad_dz           layer thickness (cm) 
c     ad_rainf        rain efficiency factor 
c
c   output parameters passed through common vfall_block:
c
c     ad_qc           mixing ratio of condensed condensate (g/g)
c   
c
c   A. Ackerman Feb-2000
c
      implicit none


c   Declare common storage 

      double precision ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf

      common / advdiff_block /
     $  ad_qbelow, ad_qvs, ad_mixl, ad_dz, ad_rainf


c   Declare local storage

      double precision qt, ad_qc


c   All vapor in excess of saturation condenses (supsat=0)
      ad_qc = max( 0., qt - ad_qvs )

c   Difference from advective-diffusive balance 
      advdiff = 
     $   ad_qbelow*exp( - ad_rainf*ad_qc*ad_dz / ( qt*ad_mixl ) ) - qt

      return
      end
