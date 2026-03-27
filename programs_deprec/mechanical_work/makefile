# Macros:
FC = ifx
# FFLAGS = -qopenmp -g -traceback -check all -debug all -I/usr/local/hdf5/intel-16.0/openmpi-1.10.2/1.8.16/include -I/usr/local/netcdf/intel-16.0/hdf5-1.8.16/openmpi-1.10.2/4.4.0/include
# FFLAGS = -Ofast -qopenmp -I/usr/local/hdf5/intel-2021.1/openmpi-4.1.0/1.10.6/include -I/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/openmpi-4.1.0/4.7.4/include
FFLAGS = -Ofast -qopenmp $$(nf-config --fflags)
# FFLAGS = -Ofast -I/usr/local/hdf5/intel-16.0/openmpi-1.10.2/1.8.16/include -I/usr/local/netcdf/intel-16.0/hdf5-1.8.16/openmpi-1.10.2/4.4.0/include
# LDLIBS = -L/usr/local/hdf5/intel-2021.1/openmpi-4.1.0/1.10.6/include/lib64:/usr/local/netcdf/intel-2021.1/hdf5-1.10.6/openmpi-4.1.0/4.7.4/lib64 -lnetcdff -lnetcdf
LDLIBS = $$(nf-config --flibs)
LDFLAGS = -qopenmp
sources = work.f90 libinterp.f90
execut = work

#--------------------------------------

objects = $(sources:.f90=.o)

${execut}: ${objects}
	$(FC) $(LDFLAGS) $^ $(LDLIBS) -o $@
# ${execut}: ${objects}
# 	$(FC) $^ $(LDLIBS) -o $@	

%.o: %.f90
	$(FC) $(FFLAGS) -o $@ -c $< 

.PHONY: clean

clean:
	rm -rf ${execut} ${objects}

work.o: libinterp.o
