// #include "unit_test_tools.h"

#include "../ccsd_imds.h"

#define MAX_ITER 100
#define TOL      1e-8

// TEST(test_h2o_sto3g) {
int main()
{
    std::string path{"./input/h2o/STO-3G/"};
    int nelec_alpha = 5;
    int nelec_beta  = 5;
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
    int nvir = nmo - nocc;

    orbo = mo_coeff(Eigen::all, Eigen::seq(0, nocc - 1));  
    orbv = mo_coeff(Eigen::all, Eigen::seq(nocc, nmo - 1));

    dm  = orbo * orbo.transpose();

    int    iter = 0;
    double cur_energy = 0.0;
    double pre_energy = 0.0;
    double e_scf      = 0.0;
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

        printf("iter = % 3d, elec_energy = % 12.8f, err_energy = % 6.4e\n", iter, cur_energy, err_energy);

        pre_energy = cur_energy;
        iter += 1;
    }

    if (is_converged) {
        printf("SCF converged!\n");
        e_scf = cur_energy;
        printf("SCF energy = % 12.8f\n", e_scf);
    } else {
        printf("SCF not converged!\n");
        exit(1);
    }

    Int1eMO fock_mo = mo_coeff.transpose() * fock * mo_coeff;
    Int2eMO eri_mo  = make_eri_mo(eri_ao, mo_coeff);

    OV   t1(nocc, nvir);
    t1 << OV::Zero(nocc, nvir);
    OOVV t2 = OOVV(nocc, nvir);

    double e_mp2 = 0.0;

    double eri_iajb = 0.0;
    double eri_ibja = 0.0;
    double t2_ijab  = 0.0;
    
    for(int i = 0; i < nocc; ++i){
        for(int a = nocc; a < nmo; ++a){
            for(int j = 0; j < nocc; ++j){
                for(int b = nocc; b < nmo; ++b){
                    eri_iajb = get_eri_mo_element(eri_mo, i, a, j, b);
                    eri_ibja = get_eri_mo_element(eri_mo, i, b, j, a);
                    t2_ijab  =  eri_iajb / (mo_energy(i) + mo_energy(j) - mo_energy(a) - mo_energy(b));
                    e_mp2   += t2_ijab * (2 * eri_iajb - eri_ibja);
                    t2.set_element(i, j, a, b, t2_ijab);
                }
            }
        }
    }

    printf("MP2 energy = % 12.8f\n", e_mp2);

    auto foo = make_imds_foo(t1, t2, fock_mo, eri_mo);
    print_matrix(foo, "foo");
    
    auto fvv = make_imds_fvv(t1, t2, fock_mo, eri_mo);
    print_matrix(fvv, "fvv");

    auto fov = make_imds_fov(t1, t2, fock_mo, eri_mo);
    print_matrix(fov, "fov");

    auto loo = make_imds_loo(t1, t2, fock_mo, eri_mo);
    print_matrix(loo, "loo");

    auto lvv = make_imds_lvv(t1, t2, fock_mo, eri_mo);
    print_matrix(lvv, "lvv");

    auto woooo = make_imds_woooo(t1, t2, fock_mo, eri_mo);

    // OccIndex i, j, k, l;

    // FOR_OCC(i, nocc, nvir) {
    //     FOR_OCC(j, nocc, nvir) {
    //         FOR_OCC(k, nocc, nvir) {
    //             FOR_OCC(l, nocc, nvir) {
    //                 printf("- %d, %d, %d, %d, % 8.6f\n", i, j, k, l, woooo.get_element(i, j, k, l));
    //             }
    //         }
    //     }
    // }

    auto wvvvv = make_imds_wvvvv(t1, t2, fock_mo, eri_mo);

    // VirIndex a, b, c, d;

    // FOR_VIR(a, nocc, nvir) {
    //     FOR_VIR(b, nocc, nvir) {
    //         FOR_VIR(c, nocc, nvir) {
    //             FOR_VIR(d, nocc, nvir) {
    //                 printf("- %d, %d, %d, %d, % 8.6f\n", a, b, c, d, wvvvv.get_element(a, b, c, d));
    //             }
    //         }
    //     }
    // }

    auto wvoov = make_imds_wvoov(t1, t2, fock_mo, eri_mo);

    VirIndex a, b;
    OccIndex i, j;

    FOR_VIR(a, nocc, nvir) {
        FOR_VIR(b, nocc, nvir) {
            FOR_OCC(i, nocc, nvir) {
                FOR_OCC(j, nocc, nvir) {
                    printf("- %d, %d, %d, %d, % 8.6f\n", a, i, j, b, wvoov.get_element(a, i, j, b));
                }
            }
        }
    }

    return 0;
}

// TEST_MAIN()