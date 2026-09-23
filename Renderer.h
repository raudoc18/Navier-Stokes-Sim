//
// Created by Dominik on 23.09.26.
//

#ifndef CPP_PROTOTYPE_RENDERER_H
#define CPP_PROTOTYPE_RENDERER_H

#include "external/glad/include/glad/glad.h"
#include <GLFW/glfw3.h>

class Renderer {
public:
    Renderer() {
        init();
    }

    void init();

    void draw();

private:
    GLuint VBO;
};

#endif //CPP_PROTOTYPE_RENDERER_H
