#include "utils.h"

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> OV;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> OO;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> VV;

typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _OOOO;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _OOVV;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _VVVV;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _VOOV;
typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor> _VOVO;

typedef int OccIndex;
typedef int VirIndex;

#define FOR_OCC(i, nocc, nvir) for (i = 0;    i < nocc; ++i)
#define FOR_VIR(a, nocc, nvir) for (a = nocc; a < nocc + nvir;  ++a)

class OOOO {
    public:
        int nocc;
        int nvir;
        int nmo;

        OOOO(int nocc, int nvir) {
            this->nocc = nocc;
            this->nvir = nvir;
            this->nmo  = nocc + nvir;
            this->_data =  _OOOO(nocc*nocc, nocc*nocc);
            this->_data << _OOOO::Zero(nocc*nocc, nocc*nocc);
        }

        void set_element(const OccIndex i, const OccIndex j, const OccIndex k, const OccIndex l, double value) {
            this->_data(i * nocc + j, k * nocc + l) = value;
        }

        double const get_element(const OccIndex i, const OccIndex j, const OccIndex k, const OccIndex l) const {
            return this->_data(i * nocc + j, k * nocc + l);
        }

    private:
        _OOOO _data;
};

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
            VirIndex a = aa - nocc;
            VirIndex b = bb - nocc;
            this->_data(i * nocc + j, a * nvir + b) = value;
        }

        double const get_element(OccIndex i, OccIndex j, VirIndex aa, VirIndex bb) const {
            VirIndex a = aa - nocc;
            VirIndex b = bb - nocc;
            return this->_data(i * nocc + j, a * nvir + b);
        }

    private:
        _OOVV _data;
};

class VVVV {
    public:
        int nocc;
        int nvir;
        int nmo;

        VVVV(int nocc, int nvir) {
            this->nocc = nocc;
            this->nvir = nvir;
            this->nmo  = nocc + nvir;
            this->_data =  _VVVV(nvir*nvir, nvir*nvir);
            this->_data << _VVVV::Zero(nvir*nvir, nvir*nvir);
        }

        void set_element(const VirIndex aa, const VirIndex bb, const VirIndex cc, const VirIndex dd, double value) {
            VirIndex a = aa - nocc;
            VirIndex b = bb - nocc;
            VirIndex c = cc - nocc;
            VirIndex d = dd - nocc;
            this->_data(a * nvir + b, c * nvir + d) = value;
        }

        double const get_element(VirIndex aa, VirIndex bb, VirIndex cc, VirIndex dd) const {
            VirIndex a = aa - nocc;
            VirIndex b = bb - nocc;
            VirIndex c = cc - nocc;
            VirIndex d = dd - nocc;
            return this->_data(a * nvir + b, c * nvir + d);
        }

    private:
        _VVVV _data;
};

class VOOV {
    public:
        int nocc;
        int nvir;
        int nmo;

        VOOV(int nocc, int nvir) {
            this->nocc = nocc;
            this->nvir = nvir;
            this->nmo  = nocc + nvir;
            this->_data =  _VOOV(nvir*nocc, nocc*nvir);
            this->_data << _VOOV::Zero(nvir*nocc, nocc*nvir);
        }

        void set_element(const VirIndex aa, const OccIndex i, const OccIndex j, const OccIndex bb, double value) {
            VirIndex a = aa - nocc;
            VirIndex b = bb - nocc;
            this->_data(a * nocc + i, j * nvir + b) = value;
        }

        double const get_element(const VirIndex aa, const OccIndex i, const OccIndex j, const OccIndex bb) const {
            VirIndex a = aa - nocc;
            VirIndex b = bb - nocc;
            return this->_data(a * nocc + i, j * nvir + b);
        }

    private:
        _VOOV _data;
};

class VOVO {
    public:
        int nocc;
        int nvir;
        int nmo;

        VOVO(int nocc, int nvir) {
            this->nocc = nocc;
            this->nvir = nvir;
            this->nmo  = nocc + nvir;
            this->_data =  _VOVO(nvir*nocc, nvir*nocc);
            this->_data << _VOVO::Zero(nvir*nocc, nvir*nocc);
        }

        void set_element(const VirIndex aa, const OccIndex i, const VirIndex bb, const OccIndex j, double value) {
            VirIndex a = aa - nocc;
            VirIndex b = bb - nocc;
            this->_data(a * nocc + i, b * nocc + j) = value;
        }

        double const get_element(const VirIndex aa, const OccIndex i, const VirIndex bb, const OccIndex j) const {
            VirIndex a = aa - nocc;
            VirIndex b = bb - nocc;
            return this->_data(a * nocc + i, b * nocc + j);
        }

    private:
        _VOVO _data;
};

OO make_imds_foo(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2e& eri_mo);
VV make_imds_fvv(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo);
OV make_imds_fov(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo);
OO make_imds_loo(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo);
VV make_imds_lvv(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo);
OOOO make_imds_woooo(const OV& t1, const OOVV& t2, const Int1e& fock_mo, const Int2eMO& eri_mo);