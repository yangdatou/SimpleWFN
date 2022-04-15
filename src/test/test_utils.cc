#include "unit_test_tools.h"

#include "../utils.h"

TEST(test_h2o_sto3g) {

    std::string path{"./input/h2o/STO-3G/"};
    std::string ovlp_path = path + "s.dat";
    std::string kine_path = path + "t.dat";
    std::string potn_path = path + "v.dat";
    std::string eri_path  = path + "eri.dat";

    int nao  = read_nao_from_file(ovlp_path);
    ASSERT_EQUAL(nao, 7);

    Int1e s1 = read_int1e_from_file(ovlp_path, nao);
    Int1e t1 = read_int1e_from_file(kine_path, nao);
    Int1e v1 = read_int1e_from_file(potn_path, nao);
    Int2e eri = read_int2e_from_file(eri_path, nao);
    
    // Temporary variables for reading the integrals
    double val;
    int mu, nu;
    int lm, sg;

    std::ifstream input;

    // Check the overlap matrix elements
    input.open(ovlp_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(s1(mu, nu), val);
        ASSERT_EQUAL(s1(nu, mu), val);
    }
    input.close();

    // Check the kinetic energy matrix elements
    input.open(kine_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(t1(mu, nu), val);
        ASSERT_EQUAL(t1(nu, mu), val);
    }
    input.close();
    
    // Check the potential energy matrix elements
    input.open(potn_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(v1(mu, nu), val);
        ASSERT_EQUAL(v1(nu, mu), val);
    }
    input.close();

    // Check the electron repulsion integral elements
    input.open(eri_path);
    assert(input.good());

    while (input >> mu >> nu >> lm >> sg >> val) {
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
    input.close();
}


TEST(test_h2o_dzp) {

    std::string path{"./input/h2o/DZP/"};
    std::string ovlp_path = path + "s.dat";
    std::string kine_path = path + "t.dat";
    std::string potn_path = path + "v.dat";
    std::string eri_path  = path + "eri.dat";

    int nao  = read_nao_from_file(ovlp_path);
    ASSERT_EQUAL(nao, 26);

    Int1e s1 = read_int1e_from_file(ovlp_path, nao);
    Int1e t1 = read_int1e_from_file(kine_path, nao);
    Int1e v1 = read_int1e_from_file(potn_path, nao);
    Int2e eri = read_int2e_from_file(eri_path, nao);
    
    // Temporary variables for reading the integrals
    double val;
    int mu, nu;
    int lm, sg;

    std::ifstream input;

    // Check the overlap matrix elements
    input.open(ovlp_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(s1(mu, nu), val);
        ASSERT_EQUAL(s1(nu, mu), val);
    }
    input.close();

    // Check the kinetic energy matrix elements
    input.open(kine_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(t1(mu, nu), val);
        ASSERT_EQUAL(t1(nu, mu), val);
    }
    input.close();
    
    // Check the potential energy matrix elements
    input.open(potn_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(v1(mu, nu), val);
        ASSERT_EQUAL(v1(nu, mu), val);
    }
    input.close();

    // Check the electron repulsion integral elements
    input.open(eri_path);
    assert(input.good());

    while (input >> mu >> nu >> lm >> sg >> val) {
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
    input.close();
}

TEST(test_h2o_dz) {

    std::string path{"./input/h2o/DZ/"};
    std::string ovlp_path = path + "s.dat";
    std::string kine_path = path + "t.dat";
    std::string potn_path = path + "v.dat";
    std::string eri_path  = path + "eri.dat";

    int nao  = read_nao_from_file(ovlp_path);
    ASSERT_EQUAL(nao, 14);

    Int1e s1 = read_int1e_from_file(ovlp_path, nao);
    Int1e t1 = read_int1e_from_file(kine_path, nao);
    Int1e v1 = read_int1e_from_file(potn_path, nao);
    Int2e eri = read_int2e_from_file(eri_path, nao);
    
    // Temporary variables for reading the integrals
    double val;
    int mu, nu;
    int lm, sg;

    std::ifstream input;

    // Check the overlap matrix elements
    input.open(ovlp_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(s1(mu, nu), val);
        ASSERT_EQUAL(s1(nu, mu), val);
    }
    input.close();

    // Check the kinetic energy matrix elements
    input.open(kine_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(t1(mu, nu), val);
        ASSERT_EQUAL(t1(nu, mu), val);
    }
    input.close();
    
    // Check the potential energy matrix elements
    input.open(potn_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(v1(mu, nu), val);
        ASSERT_EQUAL(v1(nu, mu), val);
    }
    input.close();

    // Check the electron repulsion integral elements
    input.open(eri_path);
    assert(input.good());

    while (input >> mu >> nu >> lm >> sg >> val) {
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
    input.close();
}

TEST(test_ch4_sto3g) {

    std::string path{"./input/ch4/STO-3G/"};
    std::string ovlp_path = path + "s.dat";
    std::string kine_path = path + "t.dat";
    std::string potn_path = path + "v.dat";
    std::string eri_path  = path + "eri.dat";

    int nao  = read_nao_from_file(ovlp_path);
    ASSERT_EQUAL(nao, 9);

    Int1e s1 = read_int1e_from_file(ovlp_path, nao);
    Int1e t1 = read_int1e_from_file(kine_path, nao);
    Int1e v1 = read_int1e_from_file(potn_path, nao);
    Int2e eri = read_int2e_from_file(eri_path, nao);
    
    // Temporary variables for reading the integrals
    double val;
    int mu, nu;
    int lm, sg;

    std::ifstream input;

    // Check the overlap matrix elements
    input.open(ovlp_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(s1(mu, nu), val);
        ASSERT_EQUAL(s1(nu, mu), val);
    }
    input.close();

    // Check the kinetic energy matrix elements
    input.open(kine_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(t1(mu, nu), val);
        ASSERT_EQUAL(t1(nu, mu), val);
    }
    input.close();
    
    // Check the potential energy matrix elements
    input.open(potn_path);
    assert(input.good());

    while (input >> mu >> nu >> val) {
        mu = mu - 1;
        nu = nu - 1;
        ASSERT_EQUAL(v1(mu, nu), val);
        ASSERT_EQUAL(v1(nu, mu), val);
    }
    input.close();

    // Check the electron repulsion integral elements
    input.open(eri_path);
    assert(input.good());

    while (input >> mu >> nu >> lm >> sg >> val) {
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
    input.close();
}

TEST_MAIN()