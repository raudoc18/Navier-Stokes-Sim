//
// Created by Dominik on 11.09.26.
//

#ifndef CPP_PROTOTYPE_ADAMSBASHFORTH_H
#define CPP_PROTOTYPE_ADAMSBASHFORTH_H
#include "Matrix.h"

class AdamsBashforth {
    public:
    AdamsBashforth(datastruct::Matrix<double> &u, datastruct::Matrix<double>&un, datastruct::Matrix<double>&unm1,
        datastruct::Matrix<double>&v, datastruct::Matrix<double>&vn, datastruct::Matrix<double>&vnm1
        ) : u(u), un(un), unm1(unm1), v(v), vn(vn), vnm1(vnm1) {}

    void AB2u() const;
    void AB2v() const;
    void forwardEuleru() const;
    void forwardEulerv() const;

    private:
    datastruct::Matrix<double> &u;
    datastruct::Matrix<double> &un;
    datastruct::Matrix<double> &unm1;
    datastruct::Matrix<double> &v;
    datastruct::Matrix<double> &vn;
    datastruct::Matrix<double> &vnm1;
    void getIndices(size_t batch_idx);
};

#endif //CPP_PROTOTYPE_ADAMSBASHFORTH_H
