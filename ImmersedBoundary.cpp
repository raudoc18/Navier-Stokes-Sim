//
// Created by Dominik on 25.09.26.
//

#include "ImmersedBoundary.h"

#include <valarray>

#include "Constants.h"

double ImmersedBoundary::delta(double posx, double posy, double posux, double posuy) {
    double rx = std::abs(posx - posux) / dx;
    double ry = std::abs(posy - posuy) / dy;

    double phi_x = 0;
    double phi_y = 0;

    if (0.5 <= rx && rx <= 1.5) {
        phi_x = 1.0/6.0 * (5.0 - 3.0 * rx - std::sqrt(-3.0 * std::pow(1.0 - rx, 2) + 1));
    }
    else if (rx < 0.5) {
        phi_x = 1.0/3.0 * (1.0 + std::sqrt(-3.0 * std::pow(rx, 2) + 1));
    }
    if (0.5 <= ry && ry <= 1.5) {
        phi_y = 1.0/6.0 * (5.0 - 3.0 * ry - std::sqrt(-3.0 * std::pow(1.0 - ry, 2) + 1));
    }
    else if (ry < 0.5) {
        phi_y = 1.0/3.0 * (1.0 + std::sqrt(-3.0 * std::pow(ry, 2) + 1));
    }
    return 1.0/dx * phi_x * 1.0/dy * phi_y;
}

void ImmersedBoundary::computeForceTerms() {
    for (auto obstacle : obstacles) {
        for (int i = 0; i < obstacle->getNumVertices(); i++) {
            auto pos_x = obstacle->getVertices()[i];
            auto pos_y = obstacle->getVertices()[i + 1];

            int nearest_i = static_cast<int>(std::round(pos_x/dx));
            int nearest_j = static_cast<int>(std::round(pos_y/dy));

            double u_current = 0.0;
            double v_current = 0.0;

            for (int i = nearest_i - 1; i < nearest_i + 2; i++) {
                for (int j = nearest_j - 1; j < nearest_j + 2; j++) {
                    auto pos_u_x = i * dx - 0.5 * dx;
                    auto pos_u_y = j * dy;

                    auto pos_v_x = i * dx;
                    auto pos_v_y = j * dy - 0.5 * dy;

                    u_current += ut(j, i) * delta(pos_x, pos_y, pos_u_x, pos_u_y) * dx * dy;
                    v_current += vt(j, i) * delta(pos_x, pos_y, pos_v_x, pos_v_y) * dx * dy;
                }
            }

            double u_exp = obstacle->getSpeeds()[i];
            double v_exp = obstacle->getSpeeds()[i + 1];

            double f_x = u_exp - u_current;
            double f_y = v_exp - v_current;

            for (int i = nearest_i - 1; i < nearest_i + 2; i++) {
                for (int j = nearest_j - 1; j < nearest_j + 2; j++) {
                    auto pos_u_x = i * dx - 0.5 * dx;
                    auto pos_u_y = j * dy;

                    auto pos_v_x = i * dx;
                    auto pos_v_y = j * dy - 0.5 * dy;

                    u(j, i) += f_x * delta(pos_x, pos_y, pos_u_x, pos_u_y) * dx * dy;
                    v(j, i) += f_y * delta(pos_x, pos_y, pos_v_x, pos_v_y) * dx * dy;
                }
            }
        }
    }
}
