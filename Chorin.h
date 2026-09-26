//
// Created by Dominik on 12.09.26.
//

#ifndef CPP_PROTOTYPE_CHORIN_H
#define CPP_PROTOTYPE_CHORIN_H
#include "Matrix.h"

class Chorin {
public:
    Chorin(datastruct::Matrix<double> &p) : p(p) {
        initBoundaries();
    }

    void projection(datastruct::Matrix<double> &u, datastruct::Matrix<double> &v);

private:
    datastruct::Matrix<double> &p;
    datastruct::Matrix<double> an = datastruct::Matrix<double>();
    datastruct::Matrix<double> as = datastruct::Matrix<double>();
    datastruct::Matrix<double> ae = datastruct::Matrix<double>();
    datastruct::Matrix<double> aw = datastruct::Matrix<double>();
    datastruct::Matrix<double> ap = datastruct::Matrix<double>();

    datastruct::Matrix<double> AdotChorin(const datastruct::Matrix<double> &x);

    datastruct::Matrix<double> arangeRK(datastruct::Matrix<double> &u, datastruct::Matrix<double> &v);

    void initBoundaries();
    void laplaceSolver(datastruct::Matrix<double> &u, datastruct::Matrix<double> &v);

    void gradient(datastruct::Matrix<double> &x, datastruct::Matrix<double> &u, datastruct::Matrix<double> &v);
};

#endif //CPP_PROTOTYPE_CHORIN_H
