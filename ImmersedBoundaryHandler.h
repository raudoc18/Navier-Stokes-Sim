//
// Created by Dominik on 25.09.26.
//

#ifndef CPP_PROTOTYPE_IMMERSEDBOUNDARY_H
#define CPP_PROTOTYPE_IMMERSEDBOUNDARY_H
#include <vector>

#include "Obstacle.h"
//
// Created by Dominik on 25.09.26.
//

class ImmersedBoundaryHandler {
    public:
    ImmersedBoundaryHandler(std::vector<Obstacle*> &obstacles, datastruct::Matrix<double> &u, datastruct::Matrix<double> &v,  datastruct::Matrix<double> &ut, datastruct::Matrix<double> &vt,  datastruct::Matrix<double> &f) : obstacles(
        obstacles), u(u), v(v), ut(ut), vt(vt), f(f){
    }

    void computeForceTerms();
    private:
    static double delta(double posx, double posu, double posy, double posv);
    std::vector<Obstacle*> &obstacles;
    datastruct::Matrix<double> &u;
    datastruct::Matrix<double> &v;
    datastruct::Matrix<double> &ut;
    datastruct::Matrix<double> &vt;
    datastruct::Matrix<double> &f;
};
#endif //CPP_PROTOTYPE_IMMERSEDBOUNDARY_H
