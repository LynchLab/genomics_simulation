#Testing with different standards. 
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y
#Right now I need c++11 to support constructor delegation, which I like a lot.

MODE=-O3
STD=-std=c++11
DESTDIR=/usr/local/bin
DOCPATH=\"$(shell pwd)/../../docs/man/\"
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
VERSION=-D VERSION=\"0.0\"
INCLUDE=-I ./ 

LDLIBS=

CXXFLAGS=$(VERSION) $(STD) $(MODE) $(WARNINGS) -fno-exceptions -pedantic $(INCLUDE) -fopenmp

DEBUG=-g -D DEBUG
OBJDIR=obj
TARGET=./pedigree_sim

ifeq ($(HAVE_MAN), yes)
	CXXFLAGS += -D USE_MAN -D PATH=$(DOCPATH)
endif

SUBDIR=
#>! Just a list of flags to try at some later date.
# -oFast -ffloat-store -fexcess-precision=style -ffast-math -fno-rounding-math -fno-signaling-nans
# -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math
# -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math
# -frounding-math -fsingle-precision-constant -fcx-limited-range -fcx-fortran-rules

#>! TODO ORGANIZE THESE!

MAIN=pedigree_sim.o

#These all depend on interface.h
#Depend on db_api

OBJ=$(MAIN) 

all: $(OBJ)
	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDLIBS)

clean:
	rm -f *.o
	for DIR in $(SUBDIR); do \
		cd $$DIR; \
		rm -f *.o; \
		cd ..; \
	done
.PHONY: clean
