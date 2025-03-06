all: dp-test

DC=../dual-contouring
EIGEN=/usr/include/eigen3
IMC=../implicit_marching_cubes
LIBGEOM=../libgeom
MARCHING=../marching

INCLUDES=-I$(LIBGEOM) -I$(MARCHING) -I$(IMC) -I$(EIGEN) -I$(DC)
LIBS=   -L$(LIBGEOM)/release -lgeom \
	-L$(MARCHING)/build -lmarching \
	-L$(IMC)/build -lmarching_cubes \
	-L$(DC)/build -ldualcontour \
	-lasan -lomp

# Debug
#CXXFLAGS=-std=c++20 -Wall -pedantic -g -O0 -fsanitize=address $(INCLUDES)
# Release
CXXFLAGS=-std=c++20 -Wall -pedantic -O3 $(INCLUDES)

dp-test: dp-test.o dual-primal.o solver.o
	$(CXX) -o $@ $^ $(LIBS)

.PHONY: clean
clean:
	$(RM) *.o db-test
