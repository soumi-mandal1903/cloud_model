#
# Makefile for (1) creating executable "../eddy" (default)
#              (2) creating tar file of source, input, and output
#              (3) cleaning up object files and executable
#
# To make a tar file, type "make tar" at command line
#

# rule for creating .o from .f files

.SUFFIXES:
.SUFFIXES:   .h .f .o

#.f.o:	;	f77 -N32 -f -q -g -c $(<)
# %.o: %.f	;	
# #		g77 -ff90 -g -c $*.f
# 	gfortran -g -c $*.f
%.o: %.f ;
	gfortran -g -c $*.f

# define a string identifying version of model

VERSION = '3.2'

# object files that are linked together to create executable "eddy"
# in parent directory

OFILES=		driver.o \
		main.o \
		advdiff.o \
		calc_optics.o \
		calc_qc.o \
		dblank.o \
		eddysed.o \
		find_root.o \
		find_rootl.o \
		init_optics.o \
		layer.o \
		mie_calc.o \
		pvap_al2o3.o \
		pvap_ch4.o \
                pvap_cr.o \
		pvap_fe.o \
		pvap_h2o.o \
		pvap_gas.o \
		pvap_kcl.o \
		pvap_mgsio3.o \
		pvap_mg2sio4.o \
		pvap_mns.o \
		pvap_na2s.o \
		pvap_zns.o \
		pvap_nh3.o \
		qvs_below.o \
		read_marley.o \
		read_voyager.o \
                read_lewis.o \
		vfall.o

# create executable (default target) 

eddy:		$(OFILES)
#		g77 -g -o ../eddy $(OFILES)
#		g77 -o ../eddy $(OFILES) -L/usr/local/linux.shared/lib -lfpe
#		f77 -O -o ../eddy $(OFILES)
	gfortran -g -o ../eddy $(OFILES)

# create tar file

tar:		;
		( cd .. ; tar cvfz eddysed_$(VERSION).tgz \
		  *out* \
		  model/01.REVISIONS \
		  model/Makefile \
		  model/*.f \
		  model/*.h \
		  input/ \
		  idl/*.pro )


# clean things up

# clean:		;
# 		rm *.o ../eddy
clean: ;
		rm -f *.o ../eddy


# dependency list for object files
driver.o:	driver.f	globals.h
main.o:		main.f		globals.h
advdiff.o:	advdiff.f	globals.h
calc_optics.o:	calc_optics.f	globals.h
calc_qc.o:	calc_qc.f	globals.h
dblank.o:	dblank.f
eddysed.o:	eddysed.f	globals.h
find_root.o:	find_root.f
find_rootl.o:	find_rootl.f
init_optics.o:	init_optics.f	globals.h
layer.o:	layer.f		globals.h
mie_calc.o:	mie_calc.f
pvap_al2o3.o:	pvap_al2o3.f
pvap_ch4.o:	pvap_ch4.f
pvap_fe.o:	pvap_fe.f
pvap_h2o.o:	pvap_h2o.f
pvap_gas.o:	pvap_gas.f
pvap_kcl.o:	pvap_kcl.f
pvap_mgsio3.o:	pvap_mgsio3.f
pvap_nh3.o:	pvap_nh3.f
pvap_cr.o:	pvap_cr.f
pvap_mg2sio4.o:	pvap_mg2sio4.f
pvap_mns.o:	pvap_mns.f
pvap_na2s.o:	pvap_na2s.f
pvap_zns.o:	pvap_zns.f
qvs_below.o:	qvs_below.f
read_ack.o:	read_ack.f	globals.h
read_galileo.o:	read_galileo.f	globals.h
read_lunine.o:	read_lunine.f	globals.h
read_marley.o:	read_marley.f	globals.h
read_voyager.o:	read_voyager.f	globals.h
read_lewis.o:   read_lewis.f    globals.h
vfall.o:	vfall.f
