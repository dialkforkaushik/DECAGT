CC = clang++ -std=c++11 -O2 -fopenmp
OBJ1 = simplicial_complex.o utilities.o geometry.o discrete_exterior_calculus.o finite_element_exterior_calculus.o
HEADER = src/core/simplicial_complex.h src/core/definitions.h src/core/core_utils.h src/core/geometry.h src/core/discrete_exterior_calculus.h src/core/finite_element_exterior_calculus.h lib/Eigen
CFLAGS = -c -Wall -Isrc/core -Ilib
DOXYGEN = doxygen

doxyfile = decagt.inc

all: doc test1

doc: $(doxyfile)
	$(DOXYGEN) $(doxyfile)

test1: $(OBJ1)
	$(CC) $(OBJ1) -o $@

test1.o: tests/test1.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -o $@

simplicial_complex.o: src/core/simplicial_complex.cc src/core/definitions.h src/core/core_utils.h
	$(CC) $(CFLAGS) $< -o $@

utilities.o: src/core/core_utils.cc src/core/definitions.h src/core/core_utils.h
	$(CC) $(CFLAGS) $< -o $@

geometry.o: src/core/geometry.cc src/core/definitions.h src/core/core_utils.h src/core/simplicial_complex.h src/core/geometry.h
	$(CC) $(CFLAGS) $< -o $@

discrete_exterior_calculus.o: src/core/discrete_exterior_calculus.cc src/core/definitions.h src/core/core_utils.h src/core/simplicial_complex.h src/core/geometry.h src/core/discrete_exterior_calculus.h
	$(CC) $(CFLAGS) $< -o $@

finite_element_exterior_calculus.o: src/core/finite_element_exterior_calculus.cc src/core/definitions.h src/core/core_utils.h src/core/simplicial_complex.h src/core/geometry.h src/core/discrete_exterior_calculus.h src/core/finite_element_exterior_calculus.h
	$(CC) $(CFLAGS) $< -o $@

clean: 
	rm -rf *.o

cleanall:
	rm -rf test[0-9] *.o doc/latex/ doc/html/
