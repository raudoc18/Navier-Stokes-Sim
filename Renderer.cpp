//
// Created by Dominik on 23.09.26.
//

#include "Renderer.h"

#include <iostream>

#include "Matrix.h"
#include "Utils.h"
#include "external/glad/include/glad/glad.h"

typedef  datastruct::Matrix<double> Matrix;

void Renderer::createWindow() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        exit(1);
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // Required on macOS

    window = glfwCreateWindow(windowWidth, windowHeight, "OpenGL 4.1 on macOS", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        exit(1);
    }

    glfwMakeContextCurrent(window);

    // Initialize GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD\n";
        exit(1);
    }

    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << "\n";

    velocityShaderProgram = glCreateProgram();
    obstacleShaderProgram = glCreateProgram();

    if (velocityShaderProgram == 0 || obstacleShaderProgram == 0) {
        std::cerr << "Error when creating shader program" << std::endl;
    }

    createBuffers();
    compileShaders();
}

void Renderer::init() {
    createWindow();
}

void Renderer::createContainerBuffer() {
    // create vertices positions and indices
    float quadVertices[] = {
        // Positions   // Texture Coordinates (UV)
        -1.0f, -1.0f,   0.0f, 0.0f, // Bottom-Left
        -1.0f,  1.0f,   0.0f, 1.0f, // Top-Left
         1.0f,  1.0f,   1.0f, 1.0f, // Top-Right
         1.0f, -1.0f,   1.0f, 0.0f  // Bottom-Right
    };

    unsigned int quadIndices[] = {
        0, 1, 2,
        0, 2, 3
    };

    // gen a vertex array
    glGenVertexArrays(1, &velocityVAO);
    glGenBuffers(1, &velocityVBO);
    glGenBuffers(1, &velocityEBO);

    glBindVertexArray(velocityVAO);

    // gen the vertex and index buffers for the container and fill values
    glBindBuffer(GL_ARRAY_BUFFER, velocityVBO);
    glBufferData(GL_ARRAY_BUFFER, sizeof(quadVertices), quadVertices, GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, velocityEBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(quadIndices), quadIndices, GL_STATIC_DRAW);

    GLint stride = 4 * sizeof(float);

    // Attribute 0: Coords of the container
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, nullptr);

    // Attribute 1: TexCoords (location = 1 in vertex shader)
    glEnableVertexAttribArray(1);
    glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, (void *) (2 * sizeof(float)));

    // unbind vertex array
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void Renderer::createObstacleBuffer() {
    int numObstacles = obstacles.size();

    glGenVertexArrays(numObstacles, obstacleVAOs.data());
    glGenBuffers(numObstacles, obstacleVBOs.data());
    glGenBuffers(numObstacles, obstacleEBOs.data());

    for (int i = 0; i < numObstacles; i++) {
        std::vector<float> vertices;
        std::vector<int> indices;

        auto obstacle = obstacles[i];
        vertices.push_back(obstacle->getCenterx()/lx * 2 - 1);
        vertices.push_back(obstacle->getCentery()/ly * 2 - 1);
        for (int j = 0; j < obstacle->getNumVertices()*2; j+=2) {
            vertices.push_back(obstacle->getVertices()[j]/lx * 2 - 1);
            vertices.push_back(obstacle->getVertices()[j + 1]/ly * 2 - 1);
        }
        for (int j = 0; j < obstacle->getNumVertices(); j++) {
            indices.push_back(0);
            indices.push_back(j + 1);
            indices.push_back((j + 1) % obstacle->getNumVertices() + 1);
        }

        glBindVertexArray(obstacleVAOs[i]);
        glBindBuffer(GL_ARRAY_BUFFER, obstacleVBOs[i]);
        glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(), GL_STATIC_DRAW);

        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, obstacleEBOs[i]);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(int), indices.data(), GL_STATIC_DRAW);

        GLint stride = 2 * sizeof(float);

        // Attribute 0: Coords of the container
        glEnableVertexAttribArray(0);
        glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE, stride, nullptr);
    }

    // unbind vertex array
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);
}

void Renderer::createBuffers() {
    createContainerBuffer();

    createObstacleBuffer();

    createR32TextureForVelocity(uBuffer);
    createR32TextureForVelocity(vBuffer);
}

void Renderer::addShader(GLuint shaderProgram, const char* shaderCode, GLenum shaderType) {
    GLuint shaderObj = glCreateShader(shaderType);

    const GLchar* p[1];
    p[0] = shaderCode;

    GLint len[1];
    len[0] = static_cast<GLint>(strlen(shaderCode));

    glShaderSource(shaderObj, 1, p, len);
    glCompileShader(shaderObj);

    GLint success;
    glGetShaderiv(shaderObj, GL_COMPILE_STATUS, &success);

    if (!success) {
        GLchar infoLog[1024];
        glGetShaderInfoLog(shaderObj, 1024, nullptr, infoLog);
        std::cerr << infoLog << std::endl;
        exit(1);
    }

    glAttachShader(shaderProgram, shaderObj);
}

void Renderer::compileShaders() {
    const char* VSFileName = "./shaders/velocity/vertex.vs";
    const char* FSFileName = "./shaders/velocity/fragment.fs";
    const char* obstacleVSFileName = "./shaders/obstacle/vertex.fs";
    const char* obstacleFSFileName = "./shaders/obstacle/fragment.fs";

    std::string vs, fs;

    if (!readFile(VSFileName, vs)) {
        std::cerr << "Error when loading vertex shader" << std::endl;
    }

    addShader(velocityShaderProgram, vs.c_str(), GL_VERTEX_SHADER);

    if (!readFile(FSFileName, fs)) {
        std::cerr << "Error when loading fragment shader" << std::endl;
    }

    addShader(velocityShaderProgram, fs.c_str(), GL_FRAGMENT_SHADER);

    GLint success;
    GLchar errorLog[1024] = { 0 };

    glLinkProgram(velocityShaderProgram);

    glGetProgramiv(velocityShaderProgram, GL_LINK_STATUS, &success);

    if (!success) {
        glGetProgramInfoLog(velocityShaderProgram, 1024, nullptr, errorLog);
        std::cerr << errorLog << std::endl;
        exit(1);
    }

    glBindVertexArray(velocityVAO);

    glValidateProgram(velocityShaderProgram);
    glGetProgramiv(velocityShaderProgram, GL_VALIDATE_STATUS, &success);

    if (!success) {
        glGetProgramInfoLog(velocityShaderProgram, 1024, nullptr, errorLog);
        std::cerr << errorLog << std::endl;
        exit(1);
    }

    glBindVertexArray(0);

    //obstacle shaders

    vs, fs = "";

    if (!readFile(obstacleVSFileName, vs)) {
        std::cerr << "Error when loading vertex shader" << std::endl;
    }

    addShader(obstacleShaderProgram, vs.c_str(), GL_VERTEX_SHADER);

    if (!readFile(obstacleFSFileName, fs)) {
        std::cerr << "Error when loading fragment shader" << std::endl;
    }

    addShader(obstacleShaderProgram, fs.c_str(), GL_FRAGMENT_SHADER);

    char errorLog2[1024] = { 0 };

    glLinkProgram(obstacleShaderProgram);

    glGetProgramiv(obstacleShaderProgram, GL_LINK_STATUS, &success);

    if (!success) {
        glGetProgramInfoLog(obstacleShaderProgram, 1024, nullptr, errorLog2);
        std::cerr << errorLog2 << std::endl;
        exit(1);
    }

    glBindVertexArray(obstacleVAOs[0]);

    glValidateProgram(obstacleShaderProgram);
    glGetProgramiv(obstacleShaderProgram, GL_VALIDATE_STATUS, &success);

    if (!success) {
        glGetProgramInfoLog(obstacleShaderProgram, 1024, nullptr, errorLog2);
        std::cerr << errorLog2 << std::endl;
        exit(1);
    }

    glBindVertexArray(0);
}

void Renderer::createR32TextureForVelocity(GLuint &bufferID) {

    // create new texture binder and bind it
    glGenTextures(1, &bufferID);
    glBindTexture(GL_TEXTURE_2D, bufferID);

    // interpolation settings
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    // specify the type of the texture
    glTexImage2D(GL_TEXTURE_2D, 0, GL_R32F, nx + 2, ny + 2, 0, GL_RED, GL_FLOAT, nullptr);
}

Renderer::~Renderer() {
    glDeleteProgram(velocityShaderProgram);
    glDeleteTextures(1, &uBuffer);
    glDeleteTextures(1, &vBuffer);
    glDeleteVertexArrays(1, &velocityVAO);
}

void Renderer::updateVelocityBuffer(GLuint &bufferID, Matrix &m) {
    // Ensure tight row alignment for single-channel floats
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);

    datastruct::Matrix<float> mFloat = datastruct::Matrix<float>(0.0);

    for (int i = 0; i < nx + 2; i++) {
        for (int j = 0; j < ny + 2; j++) {
            double val = m(j, i);
            mFloat(j, i) = static_cast<float>(val);
        }
    }

    mFloat.collectAsVector();

    glBindTexture(GL_TEXTURE_2D, bufferID);
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, nx + 2, ny + 2, GL_RED, GL_FLOAT, mFloat.data());
}

void Renderer::render() {
    // clear Buffers
    glClear(GL_COLOR_BUFFER_BIT);

    updateVelocityBuffer(uBuffer, p);
    updateVelocityBuffer(vBuffer, p);

    // bind buffers to corresponding attributes for shader
    // position 0 for u
    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, uBuffer);
    glUniform1i(glGetUniformLocation(velocityShaderProgram, "u"), 0);

    // position 1 for v
    glActiveTexture(GL_TEXTURE1);
    glBindTexture(GL_TEXTURE_2D, vBuffer);
    glUniform1i(glGetUniformLocation(velocityShaderProgram, "v"), 1);

    // init rendering process
    glUseProgram(velocityShaderProgram);
    int maxVelLocation = glGetUniformLocation(velocityShaderProgram, "maxVel");
    glUniform1f(maxVelLocation, p.data()[p.maxIdx()]);
    glBindVertexArray(velocityVAO);
    glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);


    glUseProgram(obstacleShaderProgram);
    for (int i = 0; i < obstacles.size(); i++) {
        glBindVertexArray(obstacleVAOs[i]);
        glDrawElements(GL_TRIANGLES, obstacles[i]->getNumVertices()*3, GL_UNSIGNED_INT, 0);
    }

    glBindVertexArray(0);
    glUseProgram(0);

    glfwSwapBuffers(window);
    glfwPollEvents();
}
