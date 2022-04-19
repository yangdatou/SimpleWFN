#include "ccsd_imds.h"

OO make_imds_foo(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();

    OO fki(nocc, nocc);
    fki << OO::Zero(nocc, nocc);

    double eri_kcld = 0.0; // eris_ovov
    double eri_kdlc = 0.0; // eris_ovov
    double t1_ic = 0.0;    // t1
    double t1_ld = 0.0;    // t1
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

                        t1_ic   = t1(i, c - nocc);
                        t1_ld   = t1(l, d - nocc);
                        t2_ilcd = t2.get_element(i, l, c, d);

                        fki(k, i) += 2.0 * eri_kcld * t2_ilcd;
                        fki(k, i) -= eri_kdlc * t2_ilcd;

                        fki(k, i) += 2.0 * eri_kcld * t1_ic * t1_ld;
                        fki(k, i) -= eri_kdlc * t1_ic * t1_ld;
                    }
                }
            }

        }
    }
    
    fki += fock_mo(Eigen::seq(0, nocc - 1), Eigen::seq(0, nocc - 1));

    return fki;
}

VV make_imds_fvv(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();
    int nmo  = nocc + nvir;

    VV fac(nvir, nvir);
    fac << VV::Zero(nvir, nvir);

    double eri_kcld = 0.0; // eris_ovov
    double eri_kdlc = 0.0; // eris_ovov
    double t1_ka    = 0.0; // t1
    double t1_ld    = 0.0; // t1
    double t2_klad  = 0.0; // t2
    
    OccIndex k, l;
    VirIndex a, c, d;
    
    FOR_VIR(a, nocc, nvir) {
        FOR_VIR(c, nocc, nvir) {

            fac(a - nocc, c - nocc) = 0.0;

            FOR_OCC(k, nocc, nvir) {
                FOR_OCC(l, nocc, nvir) {
                    FOR_VIR(d, nocc, nvir) {

                        eri_kcld = get_eri_mo_element(eri_mo, k, c, l, d);
                        eri_kdlc = get_eri_mo_element(eri_mo, k, d, l, c);

                        t1_ka     = t1(k, a - nocc);
                        t1_ld     = t1(l, d - nocc);
                        t2_klad  = t2.get_element(k, l, a, d);

                        fac(a - nocc, c - nocc) -= 2.0 * eri_kcld * t2_klad;
                        fac(a - nocc, c - nocc) += eri_kdlc * t2_klad;
                        fac(a - nocc, c - nocc) -= 2.0 * eri_kcld * t1_ka * t1_ld;
                        fac(a - nocc, c - nocc) += eri_kdlc * t1_ka * t1_ld;
                    }
                }
            }

        }
    }
    
    fac += fock_mo(Eigen::seq(nocc, nmo - 1), Eigen::seq(nocc, nmo - 1));

    return fac;
}

OV make_imds_fov(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();
    int nmo = nocc + nvir;

    OV fkc(nocc, nvir);
    fkc << OV::Zero(nocc, nvir);

    double eri_kcld = 0.0; // eris_ovov
    double eri_kdlc = 0.0; // eris_ovov
    double t1_ld    = 0.0; // t1

    OccIndex k, l;
    VirIndex c, d;
    
    FOR_OCC(k, nocc, nvir) {
        FOR_VIR(c, nocc, nvir) {

            fkc(k, c - nocc) = 0.0;

            FOR_OCC(l, nocc, nvir) {
                FOR_VIR(d, nocc, nvir) {

                    eri_kcld = get_eri_mo_element(eri_mo, k, c, c, d);
                    eri_kdlc = get_eri_mo_element(eri_mo, k, d, l, c);

                    t1_ld  = t1(l, d - nocc);

                    fkc(k, c - nocc) += 2.0 * eri_kcld * t1_ld;
                    fkc(k, c - nocc) -= eri_kdlc * t1_ld;
                }
            }

        }
    }
    
    fkc += fock_mo(Eigen::seq(0, nocc - 1), Eigen::seq(nocc, nmo-1));

    return fkc;
}

OO make_imds_loo(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();

    OO lki(nocc, nocc);
    lki << OO::Zero(nocc, nocc);

    double t1_ic    = 0.0; // t1
    double t1_lc    = 0.0; // t1

    double fkc      = 0.0; // fock_ov

    double eri_lcki = 0.0; // eris_ovov
    double eri_kcli = 0.0; // eris_ovov

    OccIndex k, l, i;
    VirIndex c;

    lki = make_imds_foo(t1, t2, fock_mo, eri_mo);

    FOR_OCC(k, nocc, nvir) {
        FOR_OCC(i, nocc, nvir) {
            FOR_VIR(c, nocc, nvir) {
                t1_ic      = t1(i, c - nocc);
                fkc        = fock_mo(k, c);
                lki(k, i) += fkc * t1_ic;
            }
        }
    }

    FOR_OCC(k, nocc, nvir) {
        FOR_OCC(i, nocc, nvir) {
            FOR_OCC(l, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    eri_lcki = get_eri_mo_element(eri_mo, l, c, k, i);
                    eri_kcli = get_eri_mo_element(eri_mo, k, c, l, i);

                    t1_lc  = t1(l, c - nocc);

                    lki(k, i) += 2.0 * eri_lcki * t1_lc;
                    lki(k, i) -= eri_kcli * t1_lc;
                }
            }
        }
    }

    return lki;

}

VV make_imds_lvv(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();

    VV lac(nvir, nvir);
    lac << VV::Zero(nvir, nvir);

    double t1_ka    = 0.0; // t1
    double t1_kd    = 0.0; // t1

    double fkc      = 0.0; // fock_ov

    double eri_kdac = 0.0; // eris_ovov
    double eri_kcad = 0.0; // eris_ovov

    OccIndex k;
    VirIndex a, c, d;

    lac = make_imds_fvv(t1, t2, fock_mo, eri_mo);

    FOR_VIR(a, nocc, nvir) {
        FOR_VIR(c, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                t1_ka      = t1(k, a - nocc);
                fkc        = fock_mo(k, c);
                lac(a - nocc, c - nocc) -= fkc * t1_ka;
            }
        }
    }

    FOR_OCC(k, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_VIR(c, nocc, nvir) {
                FOR_VIR(d, nocc, nvir) {
                    eri_kdac = get_eri_mo_element(eri_mo, k, d, a, c);
                    eri_kcad = get_eri_mo_element(eri_mo, k, c, a, d);

                    t1_kd  = t1(k, d - nocc);

                    lac(a - nocc, c - nocc) += 2.0 * eri_kdac * t1_kd;
                    lac(a - nocc, c - nocc) -= eri_kcad * t1_kd;
                }
            }
        }
    }

    return lac;
}

def cc_Woooo(t1, t2, eris):
    eris_ovoo = np.asarray(eris.ovoo)
    Wklij  = lib.einsum('lcki,jc->klij', eris_ovoo, t1) 
    Wklij += lib.einsum('kclj,ic->klij', eris_ovoo, t1) 
    eris_ovov = np.asarray(eris.ovov)
    Wklij += lib.einsum('kcld,ijcd->klij', eris_ovov, t2) 
    Wklij += lib.einsum('kcld,ic,jd->klij', eris_ovov, t1, t1) 
    Wklij += np.asarray(eris.oooo).transpose(0,2,1,3)
    return Wklij

def cc_Wvvvv(t1, t2, eris):
    # Incore
    eris_ovvv = np.asarray(eris.get_ovvv())
    Wabcd  = lib.einsum('kdac,kb->abcd', eris_ovvv,-t1)
    Wabcd -= lib.einsum('kcbd,ka->abcd', eris_ovvv, t1) 
    Wabcd += np.asarray(_get_vvvv(eris)).transpose(0,2,1,3)
    return Wabcd

def cc_Wvoov(t1, t2, eris):
    eris_ovvv = np.asarray(eris.get_ovvv())
    eris_ovoo = np.asarray(eris.ovoo)
    Wakic  = lib.einsum('kcad,id->akic', eris_ovvv, t1) 
    Wakic -= lib.einsum('kcli,la->akic', eris_ovoo, t1) 
    Wakic += np.asarray(eris.ovvo).transpose(2,0,3,1)
    eris_ovov = np.asarray(eris.ovov)
    Wakic -= 0.5*lib.einsum('ldkc,ilda->akic', eris_ovov, t2) 
    Wakic -= 0.5*lib.einsum('lckd,ilad->akic', eris_ovov, t2) 
    Wakic -= lib.einsum('ldkc,id,la->akic', eris_ovov, t1, t1) 
    Wakic += lib.einsum('ldkc,ilad->akic', eris_ovov, t2) 
    return Wakic

def cc_Wvovo(t1, t2, eris):
    eris_ovvv = np.asarray(eris.get_ovvv())
    eris_ovoo = np.asarray(eris.ovoo)
    Wakci  = lib.einsum('kdac,id->akci', eris_ovvv, t1)
    Wakci -= lib.einsum('lcki,la->akci', eris_ovoo, t1)
    Wakci += np.asarray(eris.oovv).transpose(2,0,3,1)
    eris_ovov = np.asarray(eris.ovov)
    Wakci -= 0.5*lib.einsum('lckd,ilda->akci', eris_ovov, t2)
    Wakci -= lib.einsum('lckd,id,la->akci', eris_ovov, t1, t1)
    return Wakci