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
        printf("%d %d %f\n", a, b, val);
        int1e(a, b) = val;
        int1e(b, a) = val;
    }

    input.close();

    return int1e;
}

Int2e read_int2e_from_file(std::string file_name, int nao)
{
    // open filename
    std::ifstream input{file_name};
    assert(input.good());

    double val;
    int mu, nu, lm, sg; 
    int munu; // = mu > nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu;
    int lmsg; // = lm > sg ? lm*(lm+1)/2 + sg : sg*(sg+1)/2 + lm;

    Int2e int2e(nao * (nao + 1) / 2, nao * (nao + 1) / 2);

    while (input >> mu >> nu >> lm >> sg >> val) {
        mu = mu - 1;
        nu = nu - 1;
        lm = lm - 1;
        sg = sg - 1;
        munu = mu > nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu;
        lmsg = lm > sg ? lm*(lm+1)/2 + sg : sg*(sg+1)/2 + lm;
        printf("%d %d %d %d %d %d\n", mu, nu, lm, sg, munu, lmsg);
        int2e(munu, lmsg) = val;
        int2e(lmsg, munu) = val;
    }

    input.close();

    return int2e;
}

double get_int2e_element(Int2e int2e, int mu, int nu, int lm, int sg)
{
    int munu = mu > nu ? mu*(mu+1)/2 + nu : nu*(nu+1)/2 + mu;
    int lmsg = lm > sg ? lm*(lm+1)/2 + sg : sg*(sg+1)/2 + lm;
    return int2e(munu, lmsg);
}