#Testing with different standards. 
#c90 c89 c99 c9x c11 c1x gnu90 gnu89 gnu99 gnu9x gnu11 gnu1x c++98 c++03 gnu++98 gnu++03 c++11 c++0x gnu++11 gnu++0x c++1y gnu++1y
#Right now I need c++11 to support constructor delegation, which I like a lot.

MODE=-O3
DEBUG=-g -D DEBUG
STD=-std=c++11
WARNINGS=-Wall -Werror=return-type -Wno-write-strings
VERSION=-D VERSION=\"0.0\"
CXX=g++
INCLUDE=-I ./include -I ~/src/lz4/lib -L ~/src/lz4/lib
CXXFLAGS=$(VERSION) $(MODE) $(STD) $(WARNINGS)  -pedantic $(INCLUDE) -fopenmp -fpic -llz4

OBJDIR=obj

ifeq ($(HAVE_MAN), yes)
	CXXFLAGS += -D USE_MAN -D PATH=$(DOCPATH)
endif

SUBDIR=

VPATH = ./src
#>! Just a list of flags to try at some later date.
# -oFast -ffloat-store -fexcess-precision=style -ffast-math -fno-rounding-math -fno-signaling-nans
# -fcx-limited-range -fno-math-errno -funsafe-math-optimizations -fassociative-math -freciprocal-math
# -fassociative-math -freciprocal-math -ffinite-math-only -fno-signed-zeros -fno-trapping-math
# -frounding-math -fsingle-precision-constant -fcx-limited-range -fcx-fortran-rules

#>! TODO ORGANIZE THESE!

MAIN=mutate.o

#These all depend on interface.h
#Depend on db_api
TARGET=libmapgd.so

OBJ=obj/interface.o obj/stream_tools.o obj/fasta.o obj/sample_name.o obj/map_file.o obj/data.o obj/tmp_buffer.o obj/file_index.o obj/relatedness_data.o obj/state.o obj/correl_data.o

#all: $(OBJ)
#	$(CXX) -c $(OBJ) $(CXXFLAGS) 

all: $(OBJ)
	$(CXX) -shared -o $(TARGET) $(OBJ) $(CXXFALGS) $(INCLUDE)

obj/%.o: %.cc
	$(CXX) $(CXXFLAGS) -c -o $@ $< $(INCLUDE)

clean:
	rm -f *.o
	for DIR in $(SUBDIR); do \
		cd $$DIR; \
		rm -f *.o; \
		cd ..; \
	done
.PHONY: clean
