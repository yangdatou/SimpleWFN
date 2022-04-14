#include "utils.h"

int read_nao_from_file(std::string file_name)
{
    std::ifstream input{file_name};
    assert(input.good());

    int a, b;
    double val;

    while (input >> a >> b >> val) {
        // do nothing
    }

    return a;
}

Int1e read_int1e_from_file(std::string file_name, int nao)
{
    // open filename
    std::ifstream input{file_name};
    assert(input.good());

    double val;
    int mu, nu; 

    Int1e int1e(nao, nao);

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        int1e(mu, nu) = val;
        int1e(nu, mu) = val;
    }

    input.close();

    return int1e;
}

Int2e read_int2e_from_file(std::string file_name, int nao)
{
    // open filename
    std::ifstream input{file_name};
    assert(input.good());

    int npair = nao * (nao + 1) / 2;
    int neri  = npair * (npair + 1) / 2;

    double val;
    int mu, nu, lm, sg; 
    int munu; 
    int lmsg;
    int munulmsg;

    Int2e int2e(neri);

    while (input >> mu >> nu >> lm >> sg >> val) {
        mu = mu - 1;
        nu = nu - 1;
        lm = lm - 1;
        sg = sg - 1;
        munu = mu > nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu;
        lmsg = lm > sg ? lm*(lm+1)/2 + sg : sg*(sg+1)/2 + lm;
        munulmsg = munu > lmsg ? munu*(munu+1)/2 + lmsg : lmsg*(lmsg+1)/2 + munu;
        int2e(munulmsg) = val;
    }

    input.close();

    return int2e;
}

double get_int2e_element(const Int2e int2e, int mu, int nu, int lm, int sg)
{
    int munu = mu > nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu;
    int lmsg = lm > sg ? lm*(lm+1)/2 + sg : sg*(sg+1)/2 + lm;
    int munulmsg = munu > lmsg ? munu*(munu+1)/2 + lmsg : lmsg*(lmsg+1)/2 + munu;
    return int2e(munulmsg);
}