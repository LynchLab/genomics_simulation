#Testing with different standards. 
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y
#Right now I need c++11 to support constructor delegation, which I like a lot.

MODE=-O2
STD=-std=c++11
DESTDIR=/usr/local/bin
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
VERSION=-D VERSION=\"0.0\"
INCLUDE=-I ./ -I ../common/include/ -I ~/src/lz4/lib  -I ~/src/eigen-eigen-5a0156e40feb
LDLIBS=-lgsl -lblas -L ../common/ -lmapgd -L ~/src/lz4/lib -llz4 -fopenmp 

DEBUG=-g -D DEBUG

#CXXFLAGS=$(VERSION) $(STD) $(MODE) $(WARNINGS) -fno-exceptions -pedantic $(INCLUDE) -fopenmp -DNDEBUG #$(DEBUG) #-mpopcnt
CXXFLAGS=$(VERSION) $(STD) $(MODE) $(WARNINGS) -fno-exceptions -pedantic $(INCLUDE) -fopenmp  #$(DEBUG) #-mpopcnt

#CXXFLAGS=$(VERSION) $(STD) $(DEBUG) $(WARNINGS) -fno-exceptions -pedantic $(INCLUDE) -fopenmp

OBJDIR=obj
TARGET=call_relatedness

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

MAIN=gsl_correl.o
#MAIN=gsl_correl.o

#These all depend on interface.h
#Depend on db_api

#OBJ=$(MAIN) 

#all: $(OBJ)
#	$(CXX) $(OBJ) $(CXXFLAGS) -o $(TARGET) $(LDLIBS)

known: gsl_correl_add.o
	$(CXX) gsl_correl_add.o $(CXXFLAGS) -o $(TARGET)_add $(LDLIBS)

old: gsl_correl_old.o
	$(CXX) gsl_correl_old.o $(CXXFLAGS) -o $(TARGET)_old $(LDLIBS)

new: gsl_correl.o
	$(CXX) gsl_correl.o $(CXXFLAGS) -o $(TARGET)_new $(LDLIBS)

inter: gsl_correl_interact.o
	$(CXX) gsl_correl_interact.o $(CXXFLAGS) -o $(TARGET)_interact $(LDLIBS)

lambda: gsl_correl_lambda.o
	$(CXX) gsl_correl_lambda.o $(CXXFLAGS) -o $(TARGET)_lambda $(LDLIBS)

count: gsl_correl_count_switch.o
	$(CXX) gsl_correl_count_switch.o $(CXXFLAGS) -o $(TARGET)_count $(LDLIBS)
	

clean:
	rm -f *.o
	for DIR in $(SUBDIR); do \
		cd $$DIR; \
		rm -f *.o; \
		cd ..; \
	done
.PHONY: clean
