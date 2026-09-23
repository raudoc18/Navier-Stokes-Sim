//
// Created by Dominik on 11.09.26.
//

#ifndef CPP_PROTOTYPE_CRANKNICOLSON_H
#define CPP_PROTOTYPE_CRANKNICOLSON_H
#include "Matrix.h"

class CrankNicolson {
    public:
    CrankNicolson(datastruct::Matrix<double> &u, datastruct::Matrix<double>&un, datastruct::Matrix<double>&unm1,
       datastruct::Matrix<double>&v, datastruct::Matrix<double>&vn, datastruct::Matrix<double>&vnm1
       ) : u(u), un(un), unm1(unm1), v(v), vn(vn), vnm1(vnm1) {}

    void CN_Wrapper();

private:
    datastruct::Matrix<double> &u;
    datastruct::Matrix<double> &un;
    datastruct::Matrix<double> &unm1;
    datastruct::Matrix<double> &v;
    datastruct::Matrix<double> &vn;
    datastruct::Matrix<double> &vnm1;

    static datastruct::Matrix<double> AdotCN(const datastruct::Matrix<double> &m) ;

    datastruct::Matrix<double> arangeRK(datastruct::Matrix<double> &u) const;

    void CN(datastruct::Matrix<double> &m);
};

#endif //CPP_PROTOTYPE_CRANKNICOLSON_H
