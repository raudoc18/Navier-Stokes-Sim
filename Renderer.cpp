//
// Created by Dominik on 23.09.26.
//

#include "Renderer.h"

#include "Matrix.h"

typedef  datastruct::Matrix<float> Matrix;

void Renderer::init() {
    Matrix posData = Matrix();

    glGenBuffers(1, &VBO);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(float)*elem_cnt, posData.data(), GL_STATIC_DRAW);
}

void Renderer::draw() {
    glClear(GL_COLOR_BUFFER_BIT);
    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 0, 0);
    glDrawArrays(GL_POINTS, 0, elem_cnt);
    glDisableVertexAttribArray(0);
}
