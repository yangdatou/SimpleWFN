#include <iostream>
#include <fstream>
#include <string>
#include <cstdio>
#include <cassert>

#include "Eigen/Core"
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Int1e;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Int1eAO;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Int1eMO;

typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Int2e;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Int2eAO;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Int2eMO;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MOCoeff;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> MOEnergy;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> DensityMatrixAO;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> Matrix;
typedef Eigen::Matrix<double, Eigen::Dynamic, 1> Vector;

void print_matrix(const Matrix& mat, std::string title);

int read_nao_from_file(std::string file_name);

Int1eAO read_int1e_ao_from_file(std::string file_name, int nao);
Int2eAO read_int2e_ao_from_file(std::string file_name, int nao);

double get_eri_ao_element(const Int2eAO& eri_ao, int mu, int nu, int lm, int sg);
double get_eri_mo_element(const Int2eMO& eri_mo, int p,  int q,  int r,  int s);

Int2eMO make_eri_mo(const Int2eAO& eri, MOCoeff mo_coeff);