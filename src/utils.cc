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

Int1e read_int1e_from_file(std::string file_name, int nao)
{
    // open filename
    std::ifstream input{file_name};
    assert(input.good());

    double val;
    int mu, nu; 

    Int1e int1e(nao, nao);
    int1e << Int1e::Zero(nao, nao);

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
    int2e << Int2e::Zero(neri);

    while (input >> mu >> nu >> lm >> sg >> val) {
        mu = mu - 1;
        nu = nu - 1;
        lm = lm - 1;
        sg = sg - 1;
        munu     = INDEX(mu, nu);
        lmsg     = INDEX(lm, sg);
        munulmsg = INDEX(munu, lmsg);
        int2e(munulmsg) = val;
    }

    input.close();

    return int2e;
}

double get_eri_element(const Int2e& eri, int mu, int nu, int lm, int sg)
{
    int munu     = INDEX(mu, nu);
    int lmsg     = INDEX(lm, sg);
    int munulmsg = INDEX(munu, lmsg);
    return eri(munulmsg);
}