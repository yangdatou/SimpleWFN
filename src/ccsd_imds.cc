#include "ccsd_imds.h"

OO make_imds_foo(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();

    OO fki(nocc, nocc);
    fki << OO::Zero(nocc, nocc);

    double eri_kcld = 0.0; // eris_ovov
    double eri_kdlc = 0.0; // eris_ovov
    double t2_ilcd  = 0.0; // t2

    OccIndex k, l, i;
    VirIndex c, d;
    
    FOR_OCC(k, nocc, nvir) {
        FOR_OCC(i, nocc, nvir) {

            fki(k, i) = 0.0;

            FOR_OCC(l, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    FOR_VIR(d, nocc, nvir) {


                        eri_kcld = get_eri_mo_element(eri_mo, k, c, l, d);
                        eri_kdlc = get_eri_mo_element(eri_mo, k, d, l, c);

                        t2_ilcd    = t2.get_element(i, l, c, d);
                        fki(k, i) += 2.0 * eri_kcld * t2_ilcd;
                        fki(k, i) -= eri_kdlc * t2_ilcd;

                        printf("\nl = %d, c = %d, d = %d\n", l, c, d);
                        printf("- eri_kcld = % 12.6f\n", eri_kcld);
                        printf("- eri_kdlc = % 12.6f\n", eri_kdlc);
                        printf("- t2_ilcd  = % 12.6f\n", t2_ilcd);
                        printf("- fki(%d, %d) = % 12.6f\n", k, i, fki(k, i));

                        // fki(k, i) += 2.0 * eri_kcld * t1(i, c - nocc) * t1(l, d - nocc);
                        // fki(k, i) -= eri_kdlc * t1(i, c - nocc) * t1(l, d - nocc);
                    }
                }
            }

            printf("fki(%d, %d) = %f\n", k, i, fki(k, i));
        }
    }
    
    fki += fock_mo(Eigen::seq(0, nocc - 1), Eigen::seq(0, nocc - 1));

    return fki;
}

// def cc_Fvv(t1, t2, eris):
//     nocc, nvir = t1.shape
//     fvv = eris.fock[nocc:,nocc:]
//     eris_ovov = np.asarray(eris.ovov)
//     Fac  =-2*lib.einsum('kcld,klad->ac', eris_ovov, t2) 
//     Fac +=   lib.einsum('kdlc,klad->ac', eris_ovov, t2) 
//     Fac -= 2*lib.einsum('kcld,ka,ld->ac', eris_ovov, t1, t1) 
//     Fac +=   lib.einsum('kdlc,ka,ld->ac', eris_ovov, t1, t1) 
//     Fac += fvv 
//     return Fac 

// VV make_imds_fvv(const VO& t1, const VVOO& t2, const Int1e& fock_mo, const Int2e& eri_mo)
// {
//     int nocc = t1.cols();
//     int nvir = t1.rows();
//     int nmo  = nocc + nvir;

//     VV fab(nvir, nvir);
//     fab << VV::Zero(nvir, nvir);

//     double eri_kcld = 0.0; // eris_ovov
//     double eri_kdlc = 0.0; // eris_ovov
//     double t2_klad  = 0.0; // t2
    
//     for (int k = 0; c < nocc; ++c) {
//         for (int l = 0; l < nocc; ++l) {
//             for (int d = 0; d < nvir; ++d) {
                
//                 for (int c = 0; k < nocc; ++k) {
//                         eri_kcld = get_eri_mo_element(eri_mo, k, c, l, d);
//                         eri_kdlc = get_eri_mo_element(eri_mo, l, d, l, c);

//                     for (int a = 0; i < nocc; ++i){
//                         t2_klad  = get_t2_element(t2, k, l, a, d);
//                         fki(k, i)  = 2.0 * eri_kcld * t2_ilcd;
//                         fki(k, i) -= eri_kdlc * t2_ilcd;

//                         fki(k, i) += 2.0 * eri_kcld * t1(i, c) * t1(l, d);
//                         fki(k, i) -= eris_kdlc * t1(i, c) * t1(l, d);
//                     }
//                 }
//             }
//         }
//     }
    
//     fvv += fock_mo(Eigen::seq(nocc, nocc + nvir - 1), Eigen::seq(nocc, nocc + nvir - 1));

//     return fvv;
// }

// def cc_Fov(t1, t2, eris):
//     nocc, nvir = t1.shape
//     fov = eris.fock[:nocc,nocc:]
//     eris_ovov = np.asarray(eris.ovov)
//     Fkc  = 2*np.einsum('kcld,ld->kc', eris_ovov, t1) 
//     Fkc -=   np.einsum('kdlc,ld->kc', eris_ovov, t1) 
//     Fkc += fov 
//     return Fkc