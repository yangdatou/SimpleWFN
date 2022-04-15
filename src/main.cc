#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdio>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

#include "utils.h"

int main(int argc, char *argv[])
{
    std::string path{argv[1]}; 

    int nao      = read_nao_from_file(path + "s.dat");

    Int1e s1e   = read_int1e_from_file(path + "s.dat", nao);
    print_matrix(s1e, "Overlap matrix");
    
    Matrix s1e_sqrt = sqrt_matrix(s1e);
    print_matrix(s1e_sqrt, "Sqrt of overlap matrix");

    return 0;
}