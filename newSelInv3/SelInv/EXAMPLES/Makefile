include ../Makefile.in
include makefile_allo

SUPLDLT = ../LIB/libsupldlt.a 

LIBS = $(METIS) $(SUPLDLT) $(LAPACK) $(BLAS) -lm 

OBJS = C2Finterface.o 

all: lap2d5pt.x selinv.x

.c.o:
	$(CC) $(CFLAGS) -fPIC -shared -c $<

.f.o:
	${FC} ${FFLAGS} -c $<

lap2d5pt.x: lap2d5pt.o $(OBJS) $(SUPLDLT) $(METIS)
	$(LD) $(LDFLAGS) lap2d5pt.o $(OBJS) $(LIBS) -o $@

selinv.x: selinv.o $(OBJS) $(SUPLDLT) $(METIS) 
	$(LD) $(LDFLAGS) selinv.o $(OBJS) $(LIBS) -o $@


$(SUPLDLT):
	(cd ../LIB; make lib)

julia:
	gcc -g -fPIC -shared -c C2Finterface1.c
	gcc -g -fPIC -shared -c selinvJulia.c
	gfortran  -g -DTIMING -fPIC -shared -llapack -lblas -lm C2Finterface1.o  selinvJulia.o $(allObjs) -o julia.so 

clean:
	rm -f *.o selinv.x lap2d5pt.x *.so


forR:
	clang++ -I/Library/Frameworks/R.framework/Resources/include -DNDEBUG  -I/usr/local/include -I/usr/local/include/freetype2 -I/opt/X11/include  -I"/Library/Frameworks/R.framework/Versions/3.3/Resources/library/Rcpp/include" -I"/Users/shiyuanhe/Dropbox/project/SIinR/newSelInv3/SelInv/EXAMPLES"   -fPIC  -Wall -mtune=core2 -g -O2  -c selinvR.cpp -o selinvR.o
	clang++ -dynamiclib -Wl,-headerpad_max_install_names -undefined dynamic_lookup -single_module -multiply_defined suppress -L/Library/Frameworks/R.framework/Resources/lib -L/usr/local/lib  -llapack -lblas -lm -o sourceCpp_4.so  C2Finterface1.o selinvR.o $(allObjs) -F/Library/Frameworks/R.framework/.. -framework R -Wl,-framework -Wl,CoreFoundation
