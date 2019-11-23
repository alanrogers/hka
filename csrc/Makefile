INSTALL_DIR=$(HOME)/bin
CC = gcc
LIB =
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
 -Wconversion
########################################################
# Uncomment following line to turn off debugging
#
# OPT1 = -DNDEBUG    
#
# Uncomment following line to turn on optimization
#
# OPT2 = -O3
#
############# Don't edit below this line ###############
OPT = $(OPT1) $(OPT2)
CFLAGS = -g $(WARN) $(OPT)
TARGETS = hka
VERSION = 0_0

HKA = main.o hka.o getcic.o eprintf.o
hka: $(HKA)
	$(CC) -o $@ $(HKA) $(LIB)

main.o: hkaversion.h

hkaversion.h : Makefile
	echo "#ifndef HKA_VERSION" > hkaversion.h
	echo "#define HKA_VERSION \"$(VERSION)\"" >> hkaversion.h
	echo "#endif // HKA_VERSION" >> hkaversion.h

install : $(TARGETS)
	cp $(TARGETS) $(INSTALL_DIR)

clean :
	rm -rf *.o *~ try core $(TARGETS)

dep :
	@ echo 'Generating ' `pwd`'/depend'
	-$(CC) -MM *.c > depend

.SUFFIXES:
.SUFFIXES: .c .o
.PHONY: dep depend clean install

include depend
