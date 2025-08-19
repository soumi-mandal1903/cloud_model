       double precision function pvap_kcl( t )

       implicit none

       double precision t, pvaplog
       double precision a1,b1,c1
       double precision a2,b2,c2

* data from NIST web site (http://webbook.nist.gov)
*  A B and C for T = 1170-1466 K
       data a1,b1,c1 / 4.61668, 6910.833, -176.083 /
*  A B and C for T = 1094-1680 K
       data a2,b2,c2 / 4.78236, 7440.691, -122.709 /

* This is the first way to calculate the KCl vapor pressure from the NIST
* database. Valid from T = 1094 - 1680 K.
c         if (t .lt. 1094.) pvap_kcl = 0.0
c         if (t .ge. 1094. .and. t .le. 1680.) then
c         if (t .ge. 1094)  then
            pvaplog = a2 - b2/(t+c2)
            pvap_kcl = 10**(pvaplog) 
c         endif

       return 
       end
