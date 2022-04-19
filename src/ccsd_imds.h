#include "utils.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> OV;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> OO;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> VV;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _OOVV;
// typedef Eigen::Matrix<double, Eigen::Dynamic, 1> OVOV;

typedef int OccIndex;
typedef int VirIndex;

#define FOR_OCC(i, nocc, nvir) for (i = 0;    i < nocc; ++i)
#define FOR_VIR(a, nocc, nvir) for (a = nocc; a < nocc + nvir;  ++a)

class OOVV {
    public:
        int nocc;
        int nvir;
        int nmo;

        OOVV(int nocc, int nvir) {
            this->nocc = nocc;
            this->nvir = nvir;
            this->nmo  = nocc + nvir;
            this->_data =  _OOVV(nocc*nocc, nvir*nvir);
            this->_data << _OOVV::Zero(nocc*nocc, nvir*nvir);
        }

        void set_element(const OccIndex i, const OccIndex j, const VirIndex aa, const VirIndex bb, double value) {
            VirIndex a = aa - this -> nocc;
            VirIndex b = bb - this -> nocc;
            this->_data(i * nocc + j, a * nvir + b) = value;
        }

        double const get_element(OccIndex i, OccIndex j, VirIndex aa, VirIndex bb) const {
            VirIndex a = aa - this -> nocc;
            VirIndex b = bb - this -> nocc;
            return this->_data(i * nocc + j, a * nvir + b);
        }

    private:
        _OOVV _data;
};

OO make_imds_foo(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2e& eri_mo);