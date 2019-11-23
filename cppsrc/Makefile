INSTALL_DIR=$(HOME)/bin
CC = gcc
LIB = -lstlport_gcc -lstdc++ -larrin -lpthread
WARN = \
 -Wall \
 -Wtraditional -Wpointer-arith \
 -Wcast-qual \
 -Wcast-align \
 -Wwrite-strings  \
 -Wstrict-prototypes \
 -Wmissing-prototypes \
 -Wmissing-declarations \
 -Wnested-externs \
 -Wconversion \
-Wno-ctor-dtor-privacy
########################################################
# Uncomment following line to turn off debugging
#
OPT1 = -DNDEBUG    
#
# Uncomment following line to turn on optimization
#
OPT2 = -O3
#
############# Don't edit below this line ###############
OPT = $(OPT1) $(OPT2)
CPPFLAGS = -g $(WARN) $(OPT)
TARGETS = hka xtndata
VERSION = C++_0_4

HKA=hka.o
hka: $(HKA)
	$(CC) -o $@ $(HKA) $(LIB)

hka.o xtndata.o: hkaversion.hpp

hkaversion.hpp : Makefile
	echo "#ifndef HKA_VERSION" > hkaversion.hpp
	echo "#define HKA_VERSION \"$(VERSION)\"" >> hkaversion.hpp
	echo "#endif // HKA_VERSION" >> hkaversion.hpp

XTNDATA_OBJ=xtndata.o  
xtndata: $(XTNDATA_OBJ)
	$(CC) $(PROFILE)  -pg -o $@ $(XTNDATA_OBJ) $(LIB)

XWIBTW_OBJ=xwibtw.o  
xwibtw: $(XWIBTW_OBJ)
	$(CC) $(PROFILE)  -pg -o $@ $(XWIBTW_OBJ) $(LIB)

ISIMTREE_OBJ=xisimtree.o  
xisimtree: $(ISIMTREE_OBJ)
	$(CC) $(PROFILE)  -pg -o $@ $(ISIMTREE_OBJ) $(LIB)

install : hka
	cp hka $(INSTALL_DIR)

clean :
	rm -rf *.o *~ try core $(TARGETS) *.tmp

dep :
	@ echo 'Generating ' `pwd`'/depend'
	-$(CC) -MM *.cpp > depend

.SUFFIXES:
.SUFFIXES: .cpp .o
.PHONY: dep depend clean

include depend
