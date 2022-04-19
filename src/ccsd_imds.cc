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

OOOO make_imds_woooo(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();

    OOOO wijkl(nocc, nvir);
    
    double t1_jc    = 0.0; // t1
    double t1_ic    = 0.0; // t1
    double t1_jd    = 0.0; // t1
    double t2_ijcd  = 0.0; // t2

    double eri_lcki = 0.0; // ovoo
    double eri_kclj = 0.0; // ovoo
    double eri_kcld = 0.0; // ovov

    double w_klij   = 0.0;

    OccIndex i, j, k, l;
    OccIndex c, d;

    FOR_OCC(l, nocc, nvir) {
        FOR_OCC(k, nocc, nvir) {
            FOR_OCC(j, nocc, nvir) {
                FOR_OCC(i, nocc, nvir) {
                    w_klij = 0.0;

                    FOR_VIR(c, nocc, nvir) {
                        eri_lcki = get_eri_mo_element(eri_mo, l, c, k, i);
                        eri_kclj = get_eri_mo_element(eri_mo, k, c, l, j);

                        t1_jc    = t1(j, c - nocc);
                        t1_ic    = t1(i, c - nocc);

                        w_klij  += eri_lcki * t1_jc;
                        w_klij  += eri_kclj * t1_ic;
                    }


                    FOR_VIR(c, nocc, nvir) {
                        FOR_VIR(d, nocc, nvir) {
                            eri_kcld = get_eri_mo_element(eri_mo, k, c, l, d);

                            t1_ic    = t1(i, c - nocc);
                            t1_jd    = t1(j, d - nocc);
                            t2_ijcd  = t2.get_element(i, j, c, d);

                            w_klij += eri_kcld * t1_jd * t1_ic;
                            w_klij += eri_kcld * t2_ijcd;
                        }
                    }

                    w_klij += get_eri_mo_element(eri_mo, k, i, l, j);
                    wijkl.set_element(k, l, i, j, w_klij);
                }
            }
        }
    }

    return wijkl;
}

VVVV make_imds_wvvvv(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();

    VVVV wabcd(nocc, nvir);

    double t1_kb    = 0.0; // t1
    double t1_ka    = 0.0; // t1

    double eri_kdac = 0.0; // ovvv
    double eri_kcbd = 0.0; // ovvv

    double w_abcd   = 0.0;

    OccIndex a, b, c, d;
    OccIndex k;

    FOR_VIR(a, nocc, nvir) {
        FOR_VIR(b, nocc, nvir) {
            FOR_VIR(c, nocc, nvir) {
                FOR_VIR(d, nocc, nvir) {
                    w_abcd = 0.0;

                    FOR_OCC(k, nocc, nvir) {
                        eri_kdac = get_eri_mo_element(eri_mo, k, d, a, c);
                        eri_kcbd = get_eri_mo_element(eri_mo, k, c, b, d);

                        t1_kb    = t1(k, b - nocc);
                        t1_ka    = t1(k, a - nocc);

                        w_abcd  -= eri_kdac * t1_kb;
                        w_abcd  -= eri_kcbd * t1_ka;
                    }

                    w_abcd += get_eri_mo_element(eri_mo, a, c, b, d);
                    wabcd.set_element(a, b, c, d, w_abcd);
                }
            }
        }
    }

    return wabcd;
}

VOOV make_imds_wvoov(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();

    VOOV wakic(nocc, nvir);

    double t1_id   = 0.0; // t1
    double t1_la   = 0.0; // t1
    double t2_ilda = 0.0; // t2
    double t2_ilad = 0.0; // t2

    double eri_kcad = 0.0; // ovvv
    double eri_kcli = 0.0; // ovoo
    double eri_iack = 0.0; // ovov
    double eri_ldkc = 0.0; // ovov
    double eri_lckd = 0.0; // ovov

    double w_akic   = 0.0;

    VirIndex a, c, d;
    OccIndex k, i, l;

    FOR_VIR(a, nocc, nvir) {
        FOR_OCC(k, nocc, nvir) {
            FOR_OCC(i, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    w_akic = 0.0;

                        FOR_OCC(l, nocc, nvir) {
                            FOR_VIR(d, nocc, nvir) {
                                eri_kcad = get_eri_ao_element(eri_mo, k, c, a, d);
                                eri_kcli = get_eri_ao_element(eri_mo, k, c, i, l);
                                eri_ldkc = get_eri_ao_element(eri_mo, l, d, k, c);
                                eri_lckd = get_eri_ao_element(eri_mo, l, c, k, d);

                                t2_ilda  = t2.get_element(i, l, d, a);
                                t2_ilad  = t2.get_element(i, l, a, d);

                                t1_id    = t1(i, d - nocc);
                                t1_la    = t1(l, a - nocc);

                                w_akic  += eri_kcad * t1_id;
                                w_akic  -= eri_kcli * t1_la;

                                w_akic  += eri_iack;

                                w_akic  -= 0.5 * eri_ldkc * t2_ilda;
                                w_akic  -= 0.5 * eri_lckd * t2_ilad;

                                w_akic  -= eri_ldkc * t1_id * t1_la;
                                w_akic  += eri_ldkc * t2_ilad;
                            }
                        }

                    w_akic += get_eri_mo_element(eri_mo, i, a, c, k);
                    wakic.set_element(a, k, i, c, w_akic);
                }
            }
        }
    }

    return wakic;

}

// def cc_Wvovo(t1, t2, eris):
//     eris_ovvv = np.asarray(eris.get_ovvv())
//     eris_ovoo = np.asarray(eris.ovoo)
//     Wakci  = lib.einsum('kdac,id->akci', eris_ovvv, t1)
//     Wakci -= lib.einsum('lcki,la->akci', eris_ovoo, t1)
//     Wakci += np.asarray(eris.oovv).transpose(2,0,3,1)
//     eris_ovov = np.asarray(eris.ovov)
//     Wakci -= 0.5*lib.einsum('lckd,ilda->akci', eris_ovov, t2)
//     Wakci -= lib.einsum('lckd,id,la->akci', eris_ovov, t1, t1)
//     return Wakci

VOVO make_imds_wvovo(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    int nocc = t1.rows();
    int nvir = t1.cols();

    VOVO wakci(nocc, nvir);
    
    double t1_id   = 0.0; // t1
    double t1_la   = 0.0; // t1
    double t2_ilda = 0.0; // t2

    double eri_kdac = 0.0; // ovvv
    double eri_lcki = 0.0; // ovoo
    double eri_lckd = 0.0; // ovov
    double eri_caik = 0.0; // ovov

    double w_akci   = 0.0;
}