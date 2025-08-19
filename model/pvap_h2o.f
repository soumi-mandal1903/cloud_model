      double precision function pvap_h2o( t )

c   --------------------------------------------------------------------
c
c   calculate saturation vapor pressure (dyne/cm^2) over water:
c   liquid when T > 273 K
c   ice when colder
c
c   input temperature in K
c
c   A. Ackerman Feb-2000
c
c   --------------------------------------------------------------------

      implicit none

c   --------------------------------------------------------------------
c   define constants used in Buck's expressions
c   Buck, 1981 (J. Atmos. Sci., 20, p. 1527)

      double precision BAL, BBL, BCL, BDL
      double precision BAI, BBI, BCI, BDI

      parameter( BAL = 6.1121e3 )
      parameter( BBL = 18.729 )
      parameter( BCL = 257.87 )
      parameter( BDL = 227.3 )

      parameter( BAI = 6.1115e3 )
      parameter( BBI = 23.036 )
      parameter( BCI = 279.82 )
      parameter( BDI = 333.7 )

c   --------------------------------------------------------------------
c   define constants used in Wexler formulas
c   (see Flatau et al., 1992, J. Appl. Meteor. p. 1507)

      double precision GG0, GG1, GG2, GG3, GG4, GG5, GG6, GG7
      double precision HH0, HH1, HH2, HH3, HH4, HH5

      parameter( GG0 =-0.29912729e+4 )
      parameter( GG1 =-0.60170128e+4 )
      parameter( GG2 = 0.1887643854e+2 )
      parameter( GG3 =-0.28354721e-1 )
      parameter( GG4 = 0.17838301e-4 )
      parameter( GG5 =-0.84150417e-9 )
      parameter( GG6 = 0.44412543e-12 )
      parameter( GG7 = 0.28584870e+1 )

      parameter( HH0 = -0.58653696e+4 )
      parameter( HH1 =  0.2224103300e+2 )
      parameter( HH2 =  0.13749042e-1 )
      parameter( HH3 = -0.34031775e-4 )
      parameter( HH4 =  0.26967687e-7 )
      parameter( HH5 =  0.6918651 )

c   --------------------------------------------------------------------
c   define DO_BUCK: .true.  means use Buck's expression,
c                   .false. means use Wexler

      logical DO_BUCK
      parameter( DO_BUCK = .true. )


c   declare local storage

      double precision t, tc

c   --------------------------------------------------------------------
c   branch on temperature for liquid or ice
c   --------------------------------------------------------------------

      if( t .lt. 273.16 )then

c   --------------------------------------------------------------------
c   saturation vapor pressure over ice

        if( DO_BUCK )then 
          tc = t - 273.16
          pvap_h2o = BAI * exp( (BBI - tc/BDI)*tc / (tc + BCI) )
        else 
          pvap_h2o = 10*exp( 1.0/t* 
     $        ( HH0+(HH1+HH5*log(t)+
     $        ( HH2+(HH3+HH4*t)*t)*t)*t )  )
        endif

      else
    
c   --------------------------------------------------------------------
c   saturation vapor pressure over water
c   for T > 1050 K, fix at 600 bars 

        if( DO_BUCK )then 
          if( t .lt. 1048. )then 
            tc = t - 273.16
            pvap_h2o = BAL * exp( (BBL - tc/BDL)*tc / (tc + BCL) )
          else
            pvap_h2o = 600.e6
          endif
        else 
          pvap_h2o = 10*exp( (1.0/(t*t))* 
     $        ( GG0+(GG1+(GG2+GG7*log(t)+
     $        ( GG3+(GG4+(GG5+GG6*t)*t)*t)*t)*t)*t ) )

        endif

      endif

      return
      end
