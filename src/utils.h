#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cassert>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Int1e;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Int2e;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

int    read_nao_from_file(std::string file_name);
double get_int2e_element(Int2e int2e, int mu, int nu, int lm, int sg);

Int1e read_int1e_from_file(std::string file_name, int nao);
Int2e read_int2e_from_file(std::string file_name, int nao);