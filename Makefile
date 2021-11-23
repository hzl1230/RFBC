BASEPATH=/mnt/c/Users/dell/Desktop/Plasma/eduPIC-main/RFBC
CXX=g++
CFLAGS=-std=c++11 -Wall -g -O2
PROG=main

OBJS=main.o espic_info.o mesh.o poisson.o scalar_field.o

EIGEN_PATH=${BASEPATH}/ThirdParty
EIGEN=${EIGEN_PATH}/Eigen3.3.7

INCLUDES=-I$(EIGEN)

LIBOBJDIR=Object 
LIBOBJ=libobject.a
LIBS=-L$(LIBOBJDIR)
LINKOPTS=-lobject -lm

all : libobject $(OBJS)
	$(CXX) $(CFLAGS) $(LIBS) $(OBJS) -o $(PROG) $(LINKOPTS)

libobject :
	@cd $(LIBOBJDIR); \
	make -f Makefile.object $(LIBOBJ) \
	CXX='$(CXX)' CFLAGS='$(CFLAGS)' \
	INCLUDES='$(INCLUDES)' LIBINJ='$(LIBOBJ)'

.cpp.o :
	$(CXX) $(CFLAGS) $(INCLUDES) -c $<

run :
	./$(PROG)

object_clean :
	cd $(LIBOBJDIR); make -f Makefile.object clean

clean : object_clean
	/bin/rm -f *.o

distclean: clean
