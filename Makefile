MULTICORE = true

CC = g++ -std=c++14 -O3 -fopenmp
OBJ1 = test1.o simplicial_complex.o utilities.o geometry.o discrete_exterior_calculus.o finite_element_exterior_calculus.o
OBJ2 = test2.o simplicial_complex.o utilities.o geometry.o discrete_exterior_calculus.o finite_element_exterior_calculus.o
OBJ3 = test3.o simplicial_complex.o utilities.o geometry.o discrete_exterior_calculus.o finite_element_exterior_calculus.o
OBJ4 = test4.o simplicial_complex.o utilities.o geometry.o discrete_exterior_calculus.o finite_element_exterior_calculus.o
OBJ5 = test5.o simplicial_complex.o utilities.o geometry.o discrete_exterior_calculus.o finite_element_exterior_calculus.o
OBJ6 = test6.o simplicial_complex.o utilities.o geometry.o discrete_exterior_calculus.o finite_element_exterior_calculus.o
HEADER = src/core/simplicial_complex.h src/core/definitions.h src/core/core_utils.h src/core/geometry.h src/core/discrete_exterior_calculus.h src/core/finite_element_exterior_calculus.h lib/Eigen

ifeq ($(MULTICORE), true)
	CFLAGS = -c -Wall -Isrc/core -Ilib -DMULTICORE
else
	CFLAGS = -c -Wall -Isrc/core -Ilib
endif

DOXYGEN = doxygen

doxyfile = decagt.inc

all: doc test1 test2 test3 test4 test5 test6

doc: $(doxyfile)
	$(DOXYGEN) $(doxyfile)

test1: $(OBJ1)
	$(CC) $(OBJ1) -o $@

test2: $(OBJ2)
	$(CC) $(OBJ2) -o $@

test3: $(OBJ3)
	$(CC) $(OBJ3) -o $@

test4: $(OBJ4)
	$(CC) $(OBJ4) -o $@

test5: $(OBJ5)
	$(CC) $(OBJ5) -o $@

test6: $(OBJ6)
	$(CC) $(OBJ6) -o $@

test1.o: tests/analytical.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -o $@

test2.o: tests/index_sets.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -o $@

test3.o: tests/bb_error.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -o $@

test4.o: tests/global_mass_matrix.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -o $@

test5.o: tests/random.cpp $(HEADER)
	$(CC) $(CFLAGS) $< -o $@

test6.o: tests/stiffness.cpp $(HEADER)
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
