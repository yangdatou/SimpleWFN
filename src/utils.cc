#include "utils.h"

#define INDEX(mu, nu) mu > nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu

void print_matrix(const Matrix& mat, std::string title)
{   
    printf("\n%s\n", title.c_str());
    for (int j = 0; j < mat.cols(); ++j) {
        if (j == 0) {
            printf("% 13d", j);
        }
        else{
            printf("% 9d", j);
        }
        
    }
    printf("\n");

    for (int i = 0; i < mat.rows(); ++i) {
        printf(" -%2d ", i);
        for (int j = 0; j < mat.cols(); ++j) {
            
            printf(" % 8.4f", mat(i, j));
        }
        printf("\n");
    }
}

int read_nao_from_file(std::string file_name)
{   
    // Open the file
    std::ifstream input{file_name};

    // Check if the file is open
    assert(input.good());

    double val;
    int mu, nu;

    while (input >> mu >> nu >> val) {
        // do nothing
    }

    input.close();

    return mu;
}

Int1eAO read_int1e_ao_from_file(std::string file_name, int nao)
{
    // open filename
    std::ifstream input{file_name};
    assert(input.good());

    double val;
    int mu, nu; 

    Int1eAO int1e_ao(nao, nao);
    int1e_ao << Int1eAO::Zero(nao, nao);

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        int1e_ao(mu, nu) = val;
        int1e_ao(nu, mu) = val;
    }

    input.close();

    return int1e_ao;
}

Int2eAO read_int2e_ao_from_file(std::string file_name, int nao)
{
    // open filename
    std::ifstream input{file_name};
    assert(input.good());

    int npair_ao = nao * (nao + 1) / 2;
    int neri_ao  = npair_ao * (npair_ao + 1) / 2;

    double val;
    int mu, nu, lm, sg; 
    int munu; 
    int lmsg;
    int munulmsg;

    Int2eAO int2e_ao(neri_ao);
    int2e_ao << Int2eAO::Zero(neri_ao);

    while (input >> mu >> nu >> lm >> sg >> val) {
        mu = mu - 1;
        nu = nu - 1;
        lm = lm - 1;
        sg = sg - 1;
        munu     = INDEX(mu, nu);
        lmsg     = INDEX(lm, sg);
        munulmsg = INDEX(munu, lmsg);
        int2e_ao(munulmsg) = val;
    }

    input.close();

    return int2e_ao;
}

double get_eri_ao_element(const Int2eAO& eri_ao, int mu, int nu, int lm, int sg)
{
    int munu     = INDEX(mu, nu);
    int lmsg     = INDEX(lm, sg);
    int munulmsg = INDEX(munu, lmsg);
    return eri_ao(munulmsg);
}

double get_eri_mo_element(const Int2eMO& eri_mo, int p, int q, int r, int s)
{
    int pq   = INDEX(p, q);
    int rs   = INDEX(r, s);
    int pqrs = INDEX(pq, rs);
    return eri_mo(pqrs);
}

Int2eMO make_eri_mo(const Int2eAO& eri_ao, MOCoeff mo_coeff)
{   
    int p, q, r, s;
    int pq, rs, pqrs;

    int mu, nu, lm, sg;
    int munu, lmsg;

    int nao   = mo_coeff.rows();
    int nmo   = mo_coeff.cols();

    int npair_ao = nao * (nao + 1) / 2;
    int npair_mo = nmo * (nmo + 1) / 2;
    int neri_mo  = npair_mo * (npair_mo + 1) / 2;

    double eri_mu_nu_lm_sg;

    Matrix tmp(npair_mo, npair_ao);
    Int1eAO xx_ao(nao, nao);
    Int1eMO xx_mo(nmo, nmo);

    tmp    << Matrix::Zero(npair_mo, npair_ao);
    xx_ao  << Int1eAO::Zero(nao, nao);
    xx_mo  << Int1eMO::Zero(nmo, nmo);

    for(mu = 0, munu = 0; mu < nao; mu++) {
        for(nu = 0; nu <= mu; nu++, munu++) {

            for(lm = 0, lmsg = 0; lm < nao; lm++) {
                for(sg = 0; sg <= lm; sg++, lmsg++) {
                    eri_mu_nu_lm_sg = get_eri_ao_element(eri_ao, mu, nu, lm, sg);
                    xx_ao(lm, sg) = xx_ao(sg, lm) = eri_mu_nu_lm_sg;
                }
            }

            xx_mo = mo_coeff.transpose() * xx_ao * mo_coeff;

            for(p = 0, pq = 0; p < nmo; p++) {
                for(q = 0; q <= p; q++, pq++) {
                    tmp(pq, munu) = xx_mo(p, q);
                }
            }
        }
    }

    Int2eMO int2e_mo(neri_mo);
    int2e_mo << Int2eMO::Zero(neri_mo);

    for(p = 0, pq = 0; p < nmo; p++) {
        for(q = 0; q <= p; q++, pq++) {

            for(lm = 0, lmsg = 0; lm < nao; lm++) {
                for(sg = 0; sg <= lm; sg++, lmsg++) {
                    xx_ao(lm, sg) = xx_ao(sg, lm) = tmp(pq, lmsg);
                }
            }

            xx_mo = mo_coeff.transpose() * xx_ao * mo_coeff;

            for (r = 0, rs = 0; r < nmo; r++) {
                for (s = 0; s <= r; s++, rs++) {
                    pqrs = INDEX(pq, rs);
                    int2e_mo(pqrs) = xx_mo(r, s);
                }
            }
        }
    }
                                
    return int2e_mo;
}