# Discretizations of Exterior Calculus for Analysis, Geometry and Topology (DECAGT)

# Description #

Exterior calculus or calculus on manifolds is a significant generalization of vector calculus. It provides an elegant, coordinate-independent description of physical phenomena like electromagnetic fields as well as an intrinsic representation for operators of multivariable calculus like gradient, divergence and Laplacians. Built into it is also a philosophical separation of operators into purely topological and geometrical ones.

DECAGT provides a general, extendable software framework for discretizations of the objects and operators of exterior calculus. In this public beta, we provide two discretizations of exterior calculus, namely, discrete exterior calculus (DEC) and finite element exterior calculus (FEEC). In addition, this package will provide support for ancillary differential geometric and topological data analysis computations which can reuse the underlying simplicial discretization structure for spaces on which objects and operators are constructed.


## Version: 0.5.0 (Public beta) ##

# Credits #

Authors (alphabetically sorted by last name):

Pranav Jain

Kaushik Kalyanaraman

The concept and design of the software architecture is provided by Kaushik Kalyanaraman  with additional inputs from Pranav Jain. The development of code is due to Pranav Jain with inputs from Kaushik Kalyanaraman.


# Acknowledgements #

The origin of this project is from trying to rewrite the Python codebase of PyDEC: A Python Library for Discretizations of Exterior Calculus (https://github.com/hirani/pydec) in C++ with high performance compute abilities. Consequently, this project would not have been easily realized without the algorithmic heavy lifting already performed by PyDEC.
