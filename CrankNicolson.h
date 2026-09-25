//
// Created by Dominik on 11.09.26.
//

#ifndef CPP_PROTOTYPE_CRANKNICOLSON_H
#define CPP_PROTOTYPE_CRANKNICOLSON_H
#include "Matrix.h"

class CrankNicolson {
    public:
    CrankNicolson(datastruct::Matrix<double> &u, datastruct::Matrix<double>&v
       ) : u(u), v(v) {}

    void CN_Wrapper();

private:
    datastruct::Matrix<double> &u;
    datastruct::Matrix<double> &v;

    static datastruct::Matrix<double> AdotCN(const datastruct::Matrix<double> &m) ;

    datastruct::Matrix<double> arangeRK(datastruct::Matrix<double> &u) const;

    void CN(datastruct::Matrix<double> &m);
};

#endif //CPP_PROTOTYPE_CRANKNICOLSON_H
