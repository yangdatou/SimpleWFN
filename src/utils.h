#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cassert>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Int1e;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Int2e;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

void print_matrix(const Matrix& mat, std::string title);

int read_nao_from_file(std::string file_name);

Int1e read_int1e_from_file(std::string file_name, int nao);
Int2e read_int2e_from_file(std::string file_name, int nao);

double get_eri_element(const Int2e& eri, int mu, int nu, int lm, int sg);