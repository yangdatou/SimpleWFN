#include "ccsd.h"

void update_amps(const OV& pre_t1, const OOVV& pre_t2, OV& cur_t1, OOVV& cur_t2, const Int1e& fock_mo, const Int2eMO& eri_mo)
{
    OccIndex i, j, k, l;
    VirIndex a, b, c, d;

    int nocc = pre_t1.rows();
    int nvir = pre_t1.cols();
    int nmo  = nocc + nvir;

    auto foo   = make_imds_foo(pre_t1, pre_t2, fock_mo, eri_mo);
    auto fvv   = make_imds_fvv(pre_t1, pre_t2, fock_mo, eri_mo);
    auto fov   = make_imds_fov(pre_t1, pre_t2, fock_mo, eri_mo);
    auto loo   = make_imds_loo(pre_t1, pre_t2, fock_mo, eri_mo);
    auto lvv   = make_imds_lvv(pre_t1, pre_t2, fock_mo, eri_mo);

    auto woooo = make_imds_woooo(pre_t1, pre_t2, fock_mo, eri_mo);
    auto wvvvv = make_imds_wvvvv(pre_t1, pre_t2, fock_mo, eri_mo);
    auto wvoov = make_imds_wvoov(pre_t1, pre_t2, fock_mo, eri_mo);
    auto wvovo = make_imds_wvovo(pre_t1, pre_t2, fock_mo, eri_mo);

    cur_t1 << OV::Zero(nocc, nvir);
    cur_t2.set_zero();

    OV fock_ov = fock_mo(Eigen::seq(0, nocc - 1), Eigen::seq(nocc, nmo-1));
    OO fock_oo = fock_mo(Eigen::seq(0, nocc - 1), Eigen::seq(0, nocc - 1));
    VV fock_vv = fock_mo(Eigen::seq(nocc, nmo-1), Eigen::seq(nocc, nmo-1));

    // # Move energy terms to the other side
    foo -= fock_oo;
    fvv -= fock_vv;
    loo -= fock_oo;
    lvv -= fock_vv;

    // # T1 equation
    // t1new  =-2*np.einsum('kc,ka,ic->ia', fock_ov, t1, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    cur_t1(i, a - nocc) -= 2.0 * fock_mo(k, c) * pre_t1(k, a - nocc) * pre_t1(i, c - nocc);
                }
            }
        }
    }

    // t1new +=   np.einsum('ac,ic->ia', Fvv, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_VIR(c, nocc, nvir) {
                cur_t1(i, a - nocc) += fvv(a - nocc, c - nocc) * pre_t1(i, c - nocc);
            }
        }
    }

    // t1new +=  -np.einsum('ki,ka->ia', Foo, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                cur_t1(i, a - nocc) -= foo(k, i) * pre_t1(k, a - nocc);
            }
        }
    }

    // t1new += 2*np.einsum('kc,kica->ia', Fov, t2)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    cur_t1(i, a - nocc) += 2.0 * fov(k, c - nocc) * pre_t2.get_element(k, i, c, a);
                }
            }
        }
    } 

    // t1new +=  -np.einsum('kc,ikca->ia', Fov, t2)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    cur_t1(i, a - nocc) -= fov(k, c - nocc) * pre_t2.get_element(i, k, c, a);
                }
            }
        }
    }

    // t1new +=   np.einsum('kc,ic,ka->ia', Fov, t1, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    cur_t1(i, a - nocc) += fov(k, c - nocc) * pre_t1(i, c - nocc) * pre_t1(k, a - nocc);
                }
            }
        }
    }

    // t1new += fov
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            cur_t1(i, a - nocc) += fov(i, a - nocc);
        }
    }

    // t1new += 2*np.einsum('kcai,kc->ia', eris.ovvo, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    cur_t1(i, a - nocc) += 2.0 * get_eri_mo_element(eri_mo, k, c, a, i) * pre_t1(k, c - nocc);
                }
            }
        }
    }

    // t1new +=  -np.einsum('kiac,kc->ia', eris.oovv, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    cur_t1(i, a - nocc) -= get_eri_mo_element(eri_mo, k, i, a, c) * pre_t1(k, c - nocc);
                }
            }
        }
    }

    // eris_ovvv = np.asarray(eris.get_ovvv())
    // t1new += 2*lib.einsum('kdac,ikcd->ia', eris_ovvv, t2)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(d, nocc, nvir) {
                    FOR_VIR(c, nocc, nvir) {
                        cur_t1(i, a - nocc) += 2.0 * get_eri_mo_element(eri_mo, k, d, a, c) * pre_t2.get_element(i, k, c, d);
                    }
                }
            }
        }
    }

    // t1new +=  -lib.einsum('kcad,ikcd->ia', eris_ovvv, t2)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(d, nocc, nvir) {
                    FOR_VIR(c, nocc, nvir) {
                        cur_t1(i, a - nocc) -= get_eri_mo_element(eri_mo, k, c, a, d) * pre_t2.get_element(i, k, c, d);
                    }
                }
            }
        }
    }

    // t1new += 2*lib.einsum('kdac,kd,ic->ia', eris_ovvv, t1, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(d, nocc, nvir) {
                    FOR_VIR(c, nocc, nvir) {
                        cur_t1(i, a - nocc) += 2.0 * get_eri_mo_element(eri_mo, k, d, a, c) * pre_t1(k, d - nocc) * pre_t1(i, c - nocc);
                    }
                }
            }
        }
    }

    // t1new +=  -lib.einsum('kcad,kd,ic->ia', eris_ovvv, t1, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_OCC(k, nocc, nvir) {
                FOR_VIR(d, nocc, nvir) {
                    FOR_VIR(c, nocc, nvir) {
                        cur_t1(i, a - nocc) -= get_eri_mo_element(eri_mo, k, c, a, d) * pre_t1(k, d - nocc) * pre_t1(i, c - nocc);
                    }
                }
            }
        }
    }

    // eris_ovoo = np.asarray(eris.ovoo, order='C')
    // t1new +=-2*lib.einsum('lcki,klac->ia', eris_ovoo, t2)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_VIR(c, nocc, nvir) {
                FOR_OCC(k, nocc, nvir) {
                    FOR_OCC(l, nocc, nvir) {
                        cur_t1(i, a - nocc) -= 2.0 * get_eri_mo_element(eri_mo, l, c, k, i) * pre_t2.get_element(k, l, a, c);
                    }
                }
            }
        }
    }

    // t1new +=   lib.einsum('kcli,klac->ia', eris_ovoo, t2)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_VIR(c, nocc, nvir) {
                FOR_OCC(k, nocc, nvir) {
                    FOR_OCC(l, nocc, nvir) {
                        cur_t1(i, a - nocc) += get_eri_mo_element(eri_mo, k, c, l, i) * pre_t2.get_element(k, l, a, c);
                    }
                }
            }
        }
    }

    // t1new +=-2*lib.einsum('lcki,lc,ka->ia', eris_ovoo, t1, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_VIR(c, nocc, nvir) {
                FOR_OCC(k, nocc, nvir) {
                    FOR_OCC(l, nocc, nvir) {
                        cur_t1(i, a - nocc) -= 2.0 * get_eri_mo_element(eri_mo, k, c, l, i) * pre_t1(l, c - nocc) * pre_t1(k, a - nocc);
                    }
                }
            }
        }
    }

    // t1new +=   lib.einsum('kcli,lc,ka->ia', eris_ovoo, t1, t1)
    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            FOR_VIR(c, nocc, nvir) {
                FOR_OCC(k, nocc, nvir) {
                    FOR_OCC(l, nocc, nvir) {
                        cur_t1(i, a - nocc) += get_eri_mo_element(eri_mo, k, c, l, i) * pre_t1(l, c - nocc) * pre_t1(k, a - nocc);
                    }
                }
            }
        }
    }

    FOR_OCC(i, nocc, nvir) {
        FOR_VIR(a, nocc, nvir) {
            cur_t1(i, a - nocc) /= fock_mo(i, i) - fock_mo(a, a);
        }
    }

    // T2 equation
    // tmp2  = lib.einsum('kibc,ka->abic', eris.oovv, -t1)
    // tmp2 += np.asarray(eris_ovvv).conj().transpose(1,3,0,2)
    double eri_kibc  = 0.0;
    double eri_iacb  = 0.0;
    VVOV  tmp_vvov(nocc, nvir);
    OOVV  tmp_oovv(nocc, nvir);

    double tmp_abic = 0.0;

    FOR_VIR(a, nocc, nvir) {
        FOR_VIR(b, nocc, nvir) {
            FOR_OCC(i, nocc, nvir) {
                FOR_VIR(c, nocc, nvir) {
                    tmp_abic = 0.0;

                    FOR_OCC(k, nocc, nvir) {
                        eri_kibc = get_eri_mo_element(eri_mo, k, i, b, c);
                        eri_iacb = get_eri_mo_element(eri_mo, i, a, c, b);
                        tmp_abic -= eri_kibc * pre_t1(k, a - nocc);
                        tmp_abic += eri_iacb;
                    }

                    tmp_vvov.set_element(a, b, i, c, tmp_abic);
                }

            }
        }
    }

    // tmp = lib.einsum('abic,jc->ijab', tmp2, t1)
    double tmp_ijab = 0.0;
    FOR_OCC(i, nocc, nvir) {
        FOR_OCC(j, nocc, nvir) {
            FOR_VIR(a, nocc, nvir) {
                FOR_VIR(b, nocc, nvir) {
                    tmp_ijab = 0.0;

                    FOR_VIR(c, nocc, nvir) {
                        tmp_abic  = tmp_vvov.get_element(a, b, i, c);
                        tmp_ijab += tmp_abic * pre_t1(j, c - nocc);
                    }

                    tmp_oovv.set_element(i, j, a, b, tmp_ijab);
                }
            }
        }
    }

    // t2new = tmp + tmp.transpose(1,0,3,2)
    double t2_ijab = 0.0;
    FOR_OCC(i, nocc, nvir) {
        FOR_OCC(j, nocc, nvir) {
            FOR_VIR(a, nocc, nvir) {
                FOR_VIR(b, nocc, nvir) {
                    t2_ijab = tmp_oovv.get_element(i, j, a, b) + tmp_oovv.get_element(j, i, b, a);
                    cur_t2.set_element(i, j, a, b, t2_ijab);
                }
            }
        }
    }

    // tmp2  = lib.einsum('kcai,jc->akij', eris.ovvo, t1)
    tmp_vooo.set_zero();

    // tmp2 += eris_ovoo.transpose(1,3,0,2).conj()
    // tmp = lib.einsum('akij,kb->ijab', tmp2, t1)
    // t2new -= tmp + tmp.transpose(1,0,3,2)
    // t2new += np.asarray(eris.ovov).conj().transpose(0,2,1,3)
}