//
// Created by Dominik on 25.09.26.
//

#ifndef CPP_PROTOTYPE_OBSTACLE_H
#define CPP_PROTOTYPE_OBSTACLE_H
#include "Matrix.h"

class Obstacle {
public:
    Obstacle(float centerx, float centery): centerx(centerx), centery(centery) {
        circle();
    }

    int getNumVertices() {
        return numVertices;
    }

    std::vector<float>& getVertices() {
        return vertices;
    }

    std::vector<float>& getSpeeds() {
        return speeds;
    }

    float getCenterx() {
        return centerx;
    }

    float getCentery() {
        return centery;
    }

private:
    std::vector<float> vertices;
    std::vector<float> speeds;
    std::vector<int> verticeIndices;
    int numVertices = 30;
    float centerx;
    float centery;

    void circle();
};

#endif //CPP_PROTOTYPE_OBSTACLE_H
