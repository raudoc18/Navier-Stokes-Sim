//
// Created by Dominik on 25.09.26.
//

#include "Obstacle.h"

#include <numbers>

void Obstacle::circle() {
    double r = 0.05;
    for (int i = 0; i < numVertices; i++) {
        const auto x = static_cast<float>(centerx + r * std::cos(static_cast<double>(i) * 2.0 * std::numbers::pi/numVertices));
        const auto y = static_cast<float>(centery + r * std::sin(static_cast<double>(i) * 2.0 * std::numbers::pi/numVertices));

        vertices.push_back(x);
        vertices.push_back(y);

        speeds.push_back(0.0);
        speeds.push_back(0.0);
    }
}


