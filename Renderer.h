//
// Created by Dominik on 23.09.26.
//

#ifndef CPP_PROTOTYPE_RENDERER_H
#define CPP_PROTOTYPE_RENDERER_H

#include "Matrix.h"
#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include "Obstacle.h"

class Renderer {
public:
    Renderer(datastruct::Matrix<double> &u,  datastruct::Matrix<double> &v, std::vector<Obstacle*> &obstacles, datastruct::Matrix<double> &p): u(u), v(v), obstacles(obstacles), p(p) {
        obstacleVBOs = std::vector<GLuint>(obstacles.size());
        obstacleEBOs = std::vector<GLuint>(obstacles.size());
        obstacleVAOs = std::vector<GLuint>(obstacles.size());
        init();
    }

    ~Renderer();

    void render();

    GLFWwindow* getWindow() const {
        return window;
    }

private:
    GLuint velocityVBO;
    GLuint velocityVAO;
    GLuint velocityEBO;
    GLuint uBuffer;
    GLuint vBuffer;

    GLFWwindow* window;

    GLuint velocityShaderProgram;
    GLuint obstacleShaderProgram;

    datastruct::Matrix<double> &u;
    datastruct::Matrix<double> &v;
    datastruct::Matrix<double> &p;
    std::vector<Obstacle *> &obstacles;

    std::vector<GLuint> obstacleVBOs;
    std::vector<GLuint>  obstacleVAOs;
    std::vector<GLuint>  obstacleEBOs;

    void createWindow();

    void loadContext();

    void init();

    void createContainerBuffer();

    void createObstacleBuffer();

    void createBuffers();

    static void addShader(GLuint shaderProgram, const char *shaderCode, GLenum shaderType);

    void compileShaders();

    void createR32TextureForVelocity(GLuint &bufferID);

    void updateVelocityBuffer(GLuint &bufferID, datastruct::Matrix<double> &m);

};

#endif //CPP_PROTOTYPE_RENDERER_H
