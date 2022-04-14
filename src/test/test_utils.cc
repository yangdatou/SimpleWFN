#include "unit_test_tools.h"

#include "../utils.h"

TEST(test_read_nao_from_file) {

    int nao = read_nao_from_file("./input/h2o/STO-3G/s.dat");
    ASSERT_EQUAL(nao, 7);

}

TEST(test_read_int1e_from_file) {

    std::string s_path = "./input/h2o/STO-3G/s.dat";
    std::string t_path = "./input/h2o/STO-3G/t.dat";
    std::string v_path = "./input/h2o/STO-3G/v.dat";

    int nao  = read_nao_from_file(s_path);

    Int1e s1 = read_int1e_from_file(s_path, nao);
    Int1e t1 = read_int1e_from_file(t_path, nao);
    Int1e v1 = read_int1e_from_file(v_path, nao);
    
    std::ifstream s_input{s_path};
    assert(s_input.good());

    double val;
    int mu, nu;

    while (s_input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(s1(mu, nu), val);
        ASSERT_EQUAL(s1(nu, mu), val);
    }

    std::ifstream t_input{t_path};
    assert(t_input.good());

    while (t_input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(t1(mu, nu), val);
        ASSERT_EQUAL(t1(nu, mu), val);
    }

    std::ifstream input2{v_path};
    assert(input2.good());

    while (input2 >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(v1(mu, nu), val);
        ASSERT_EQUAL(v1(nu, mu), val);
    }
}

TEST(test_read_int2e_from_file) {

    std::string s_path   = "./input/h2o/STO-3G/s.dat";
    std::string eri_path = "./input/h2o/STO-3G/eri.dat";

    int nao  = read_nao_from_file(s_path);

    Int2e eri = read_int2e_from_file(eri_path, nao);
    
    std::ifstream eri_input{eri_path};
    assert(eri_input.good());

    double val;
    int mu, nu, lm, sg;

    while (eri_input >> mu >> nu >> lm >> sg >> val) {
        mu = mu - 1;
        nu = nu - 1;
        lm = lm - 1;
        sg = sg - 1;
        ASSERT_EQUAL(get_int2e_element(eri, mu, nu, lm, sg), val);
        ASSERT_EQUAL(get_int2e_element(eri, nu, mu, lm, sg), val);
        ASSERT_EQUAL(get_int2e_element(eri, mu, nu, sg, lm), val);
        ASSERT_EQUAL(get_int2e_element(eri, nu, mu, sg, lm), val);
        ASSERT_EQUAL(get_int2e_element(eri, lm, sg, mu, nu), val);
        ASSERT_EQUAL(get_int2e_element(eri, lm, sg, nu, mu), val);
        ASSERT_EQUAL(get_int2e_element(eri, sg, lm, mu, nu), val);
        ASSERT_EQUAL(get_int2e_element(eri, sg, lm, nu, mu), val);
    }

}

TEST_MAIN()