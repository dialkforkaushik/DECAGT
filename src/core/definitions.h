#ifndef DEFINITIONS_H
#define DEFINITIONS_H

#include <vector>
#include <Eigen/Eigen>
#include <tuple>
#include <map>


#define SUCCESS 0
#define FAILURE 1
#define True true
#define False false

typedef std::vector < std::string > VectorS;
typedef std::vector < int > VectorI;
typedef std::vector < float > VectorF;
typedef std::vector < double > VectorD;
typedef std::vector < VectorI > Vector2I;
typedef std::vector < VectorD > Vector2D;
typedef std::vector < Vector2I > Vector3I;
typedef std::vector < Vector2D > Vector3D;

typedef std::map < int, int > MapI;
typedef std::map < Vector2I, int > MapVector2I;
typedef std::map < VectorI, int > MapVectorI;
typedef std::vector < std::map < int, int > > VectorMapI;

typedef Eigen::SparseMatrix < double > SpMatD;
typedef Eigen::MatrixXd DenMatD;
typedef Eigen::SparseMatrix < int, Eigen::RowMajor > SpMatI;
typedef Eigen::SparseMatrix < int, Eigen::ColMajor > SpMatIC;
typedef Eigen::MatrixXi DenMatI;
typedef Eigen::Map < Eigen::VectorXd > MapEigVectorD;
typedef Eigen::Map < Eigen::VectorXi > MapEigVectorI;
typedef Eigen::VectorXd EigVectorD;
typedef Eigen::VectorXi EigVectorI;

typedef std::tuple < VectorD, double > TupleVD;

typedef std::vector < SpMatI > VectorSpmatI;
typedef std::vector < SpMatIC > VectorSpmatIC;
typedef std::vector < SpMatD > VectorSpmatD;
typedef std::vector < DenMatD > VectorDenMatD;
typedef std::vector < VectorMapI > VectorMap2I;
typedef std::vector < std::vector < VectorMapI > > VectorMap3I;
typedef std::vector < EigVectorD > VectorEigVectorD;

typedef Eigen::Triplet < double > TripletD;
typedef Eigen::Triplet < int > TripletI;
typedef std::vector < TripletD > VectorTripletD;
typedef std::vector < TripletI > VectorTripletI;

typedef double(*function_d)(VectorD);
typedef VectorD(*function_v)(double*);

#endif




