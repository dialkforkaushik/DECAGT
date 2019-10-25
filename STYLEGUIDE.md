# DECAGT style guide #

# APIs #

Essential dictums (adapted from The Little Manual of API Design, Jasmin
Blanchette, 2008):

	1. Shall be easy to functionally recognize
	2. Shall lead to clean and readable code
	3. Shall provide mathematical flow of notions
	4. Shall be reusable and extendible

We regard extension here to mean method, performance, and parallelization.

# Coding #

    1. Readability, hence maintainability, above all
    2. Crisp and as efficient as possible 

## Files ##

    1. Header files use extension `.h`, C++ source has extension `.cc`, and C
       source has extension `.c`.

    2. Header files are protected by an include guard, which should be the
       name of the file converted to preprocessor symbol style. Do not
       include leading or trailing underscores on the include guard symbol.
    
    3. Files should not exceed 1000 lines, and lines should not exceed
       80 characters.
