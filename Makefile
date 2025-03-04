all: dp-test

EIGEN=/usr/include/eigen3
LIBGEOM=../libgeom
MARCHING=../marching
IMC=../implicit_marching_cubes

INCLUDES=-I$(LIBGEOM) -I$(MARCHING) -I$(IMC) -I$(EIGEN)
LIBS=   -L$(LIBGEOM)/release -lgeom \
	-L$(MARCHING)/build -lmarching \
	-L$(IMC)/build -lmarching_cubes \
	-lasan

# Debug
CXXFLAGS=-std=c++20 -Wall -pedantic -g -O0 -fsanitize=address $(INCLUDES)
# Release
#CXXFLAGS=-std=c++20 -Wall -pedantic -O3 $(INCLUDES)

dp-test: dp-test.o dual-primal.o
	$(CXX) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	$(RM) *.o db-test
