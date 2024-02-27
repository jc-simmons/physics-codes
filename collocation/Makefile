FC=ifort
LINKER=ifort
FFLAGS = -c -O3 -xHost 
LFLAGS=	 -mkl -L./../ARPACK -larpack_frontenac	 
SRC1= dneupd2.f hp_sort2.f
SRC2 = sortmod.f90 indmod.f90 funcmod.f90 mapmod.f90 buildmatmod.f90 mvpmod.f90 coll_lu.f90
OBJ= ${SRC1:.f=.o} ${SRC2:.f90=.o}


%.o: %.f 
	$(FC) $(FFLAGS) $*.f -c

%.o: %.f90
	$(FC) $(FFLAGS) $*.f90 -c

coll: $(OBJ)
	$(FC) $(OBJ) -o $@ $(LFLAGS)  

clean:
	rm *.o
	rm *.mod
