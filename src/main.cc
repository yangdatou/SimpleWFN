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
    Int1e t1e   = read_int1e_from_file(path + "t.dat", nao);
    Int1e v1e   = read_int1e_from_file(path + "v.dat", nao);
    Int1e h1e   = t1e + v1e;
    Int2e eri   = read_int2e_from_file(path + "eri.dat", nao);

    std::cout << "Number of basis functions: " << nao << std::endl;
    std::cout << "Integrals:" << std::endl;
    std::cout << s1e << std::endl;
    std::cout << h1e << std::endl;

    double int2e_element = get_int2e_element(eri, 6, 5, 1, 0);
    std::cout << "Integral element: " << int2e_element << std::endl;
    int2e_element = get_int2e_element(eri, 5, 6, 1, 0);
    std::cout << "Integral element: " << int2e_element << std::endl;
    int2e_element = get_int2e_element(eri, 5, 6, 0, 1);
    std::cout << "Integral element: " << int2e_element << std::endl;
    std::cout << "Reference value: "  << 0.024847613930679 << std::endl;
    return 0;
}