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

class ImmersedBoundary {
    public:
    ImmersedBoundary(std::vector<Obstacle*> &obstacles, datastruct::Matrix<double> &u, datastruct::Matrix<double> &v,  datastruct::Matrix<double> &ut, datastruct::Matrix<double> &vt) : obstacles(obstacles), u(u), v(v), ut(ut), vt(vt) {
    }

    void computeForceTerms();
    private:
    static double delta(double posx, double posu, double posy, double posv);
    std::vector<Obstacle*> &obstacles;
    datastruct::Matrix<double> &u;
    datastruct::Matrix<double> &v;
    datastruct::Matrix<double> &ut;
    datastruct::Matrix<double> &vt;
};
#endif //CPP_PROTOTYPE_IMMERSEDBOUNDARY_H
