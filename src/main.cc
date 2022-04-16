#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdio>
#include <cassert>

#include "Eigen/Dense"
#include "Eigen/Eigenvalues"
#include "Eigen/Core"

#include "utils.h"

const double   TOL = 1e-8;
const int MAX_ITER = 100;

int main(int argc, char *argv[])
{
    std::string path{argv[1]}; 
    int nelec_alpha = atoi(argv[2]);
    int nelec_beta  = atoi(argv[3]);
    assert(nelec_alpha == nelec_beta);

    printf("nelec_alpha = %d\n", nelec_alpha);
    printf("nelec_beta  = %d\n", nelec_beta);

    int nao      = read_nao_from_file(path + "s.dat");
    int nocc     = nelec_alpha; // restricted hartree fock

    Int1eAO s1e    = read_int1e_ao_from_file(path + "s.dat", nao);
    Int1eAO t1e    = read_int1e_ao_from_file(path + "t.dat", nao);
    Int1eAO v1e    = read_int1e_ao_from_file(path + "v.dat", nao);
    Int1eAO h1e    = t1e + v1e;
    Int2eAO eri_ao = read_int2e_ao_from_file(path + "eri.dat", nao);
    
    Eigen::GeneralizedSelfAdjointEigenSolver<Matrix> solver;
    MOCoeff  mo_coeff ;
    MOEnergy mo_energy;

    MOCoeff orbo;
    MOCoeff orbv;

    DensityMatrixAO dm;
    Int1eAO       fock;

    solver.compute(h1e, s1e);
    mo_coeff   = solver.eigenvectors();
    mo_energy  = solver.eigenvalues();
    int nmo = mo_energy.rows();

    orbo = mo_coeff(Eigen::all, Eigen::seq(0, nocc - 1));  
    orbv = mo_coeff(Eigen::all, Eigen::seq(nocc, nmo - 1));

    dm  = orbo * orbo.transpose();

    int    iter = 0;
    double cur_energy = 0.0;
    double pre_energy = 0.0;
    double err_energy = 1.0;
    double eri_mu_nu_lm_sg, eri_mu_lm_nu_sg;

    bool is_converged = false;
    bool is_max_iter  = false;

    while (not is_converged and not is_max_iter){
        fock = h1e;

        // TODO: Check the locality of the memory access
        for (int mu = 0; mu < nao; ++mu) {
            for (int nu = 0; nu < nao; ++nu) {
                for (int lm = 0; lm < nao; ++lm) {
                    for (int sg = 0; sg < nao; ++sg) {
                        eri_mu_nu_lm_sg  = get_eri_ao_element(eri_ao, mu, nu, lm, sg);
                        eri_mu_lm_nu_sg  = get_eri_ao_element(eri_ao, mu, lm, nu, sg);
                        fock(mu, nu)    += dm(lm, sg) * (2 * eri_mu_nu_lm_sg - eri_mu_lm_nu_sg);
                    }
                }
            }
        }

        solver.compute(fock, s1e);
        mo_coeff   = solver.eigenvectors();
        mo_energy  = solver.eigenvalues();

        orbo = mo_coeff(Eigen::all, Eigen::seq(0, nocc - 1));  
        orbv = mo_coeff(Eigen::all, Eigen::seq(nocc, nmo - 1));
        dm   = orbo * orbo.transpose();

        cur_energy = (dm * (h1e + fock)).trace();
        err_energy = std::abs(cur_energy - pre_energy);

        is_converged = (err_energy < TOL);
        is_max_iter  = (iter > MAX_ITER);

        printf("iter = % 3d, elec_energy = %.12f, err_energy = % 6.4e\n", iter, cur_energy, err_energy);

        pre_energy = cur_energy;
        iter += 1;
    }

    Int2eMO eri_mo = make_eri_mo(eri_ao, mo_coeff);

    double e_mp2 = 0.0;

    double eri_iajb = 0.0;
    double eri_ibja = 0.0;
    
    for(int i = 0; i < nocc; ++i){
        for(int a = nocc; a < nmo; ++a){
            for(int j = 0; j < nocc; ++j){
                for(int b = nocc; b < nmo; ++b){
                    eri_iajb = get_eri_mo_element(eri_mo, i, a, j, b);
                    eri_ibja = get_eri_mo_element(eri_mo, i, b, j, a);
                    e_mp2 += eri_iajb * (2 * eri_iajb - eri_ibja) / (mo_energy(i) + mo_energy(j) - mo_energy(a) - mo_energy(b));
                }
            }
        }
    }

    printf("MP2 energy = %.12f\n", e_mp2);

    return 0;
}