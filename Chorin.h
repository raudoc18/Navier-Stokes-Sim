//
// Created by Dominik on 12.09.26.
//

#ifndef CPP_PROTOTYPE_CHORIN_H
#define CPP_PROTOTYPE_CHORIN_H
#include "Matrix.h"

class Chorin {
public:
    Chorin(datastruct::Matrix<double> &p, datastruct::Matrix<double>&ut, datastruct::Matrix<double>&vt
       ) : p(p), ut(ut), vt(vt) {
        initBoundaries();
    }

    void projection();

private:
    datastruct::Matrix<double> &p;
    datastruct::Matrix<double> &ut;
    datastruct::Matrix<double> &vt;
    datastruct::Matrix<double> an = datastruct::Matrix<double>();
    datastruct::Matrix<double> as = datastruct::Matrix<double>();
    datastruct::Matrix<double> ae = datastruct::Matrix<double>();
    datastruct::Matrix<double> aw = datastruct::Matrix<double>();
    datastruct::Matrix<double> ap = datastruct::Matrix<double>();

    datastruct::Matrix<double> AdotChorin(const datastruct::Matrix<double> &x);

    datastruct::Matrix<double> arangeRK();

    void initBoundaries();
    void laplaceSolver();

    void gradient(datastruct::Matrix<double> &x);
};

#endif //CPP_PROTOTYPE_CHORIN_H
