//
// Created by Dominik on 25.09.26.
//

#ifndef CPP_PROTOTYPE_STREAMFUNCTION_H
#define CPP_PROTOTYPE_STREAMFUNCTION_H
#include "Matrix.h"

class StreamFunction {
public:
    StreamFunction(datastruct::Matrix<double> &u): u(u) {}

    void computeStreamFunc();

    void extractLines();
private:
    datastruct::Matrix<double> &u;
    datastruct::Matrix<double> streamFunc;
};

#endif //CPP_PROTOTYPE_STREAMFUNCTION_H
