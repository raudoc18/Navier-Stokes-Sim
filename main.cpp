#include "Constants.h"
#include "Matrix.h"
#include <iostream>
#include "AdamsBashforth.h"
#include "Chorin.h"
#include "CrankNicolson.h"

typedef datastruct::Matrix<double> Matrix;

void intializeMatrices(Matrix &u, Matrix &un, Matrix &v, Matrix &vn, Matrix &p) {
    for (int i = 0; i < nx + 2; ++i) {
        for (int j = 0; j < ny + 2; ++j) {
            u(j, i) = 0.0;
            un(j, i) = 0.0;
            v(j, i) = 0.0;
            vn(j, i) = 0.0;
            p(j, i) = 0.0;
        }
    }
}

void boundaryConditions(Matrix &u, Matrix&v) {
    for (int i = 0; i < nx + 2; ++i) {
        for (int j = 0; j < ny + 2; ++j) {
            if (j == ny + 1) {
                u(j, i) = 2.0 * Ut - u(j - 1, i);
            }

            if (j == 0) {
                u(j, i) = -u(j + 1, i);
            }

            if (i == 1) {
                u(j, i) = 0.0;
            }
            if (i == nx + 1) {
                u(j, i) = 0.0;
            }

            if (i == 0) {
                v(j, i) = -v(j, i + 1);
            }
            if (i == nx + 1) {
                v(j, i) = -v(j, i - 1);
            }

            if (j == 0) {
                v(j, i) = 0;
            }

            if (j == ny) {
                v(j, i) = 0;
            }
        }
    }
}

int test() {
    if (!glfwInit()) {
        std::cerr << "Failed to initialize GLFW\n";
        return -1;
    }

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 1);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE); // Required on macOS

    GLFWwindow* window = glfwCreateWindow(800, 600, "OpenGL 4.1 on macOS", nullptr, nullptr);
    if (!window) {
        std::cerr << "Failed to create GLFW window\n";
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);

    // Initialize GLAD
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress)) {
        std::cerr << "Failed to initialize GLAD\n";
        return -1;
    }

    std::cout << "OpenGL Version: " << glGetString(GL_VERSION) << "\n";

    // Main render loop
    while (!glfwWindowShouldClose(window)) {
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT);

        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    glfwDestroyWindow(window);
    glfwTerminate();
    return 0;
}

int main() {
    int cnt = 0;
    double deltat = 0;

    auto u = Matrix();
    auto un = Matrix();
    auto unm1 = Matrix();
    auto uAB = Matrix();
    auto uCN = Matrix();

    auto v = Matrix();
    auto vn = Matrix();
    auto vnm1 = Matrix();
    auto vAB = Matrix();
    auto vCN = Matrix();

    auto p = Matrix();
    auto pn = Matrix();

    AdamsBashforth adm = AdamsBashforth(uAB, un, unm1, vAB, vn, vnm1);
    CrankNicolson cn = CrankNicolson(u, un, unm1, v, vn, vnm1);
    Chorin cho = Chorin(p, u, v);

    intializeMatrices(u, un, v, vn, p);

    test();

    while (cnt < ntimesteps && deltat < animduration) {

        boundaryConditions(u, v);

        if (cnt == 0) {
            adm.forwardEuleru();
            adm.forwardEulerv();
        }
        else {
            adm.AB2u();
            adm.AB2v();
        }

        cn.CN_Wrapper();

        u.add(uAB);
        v.add(vAB);

        cho.projection();

        unm1 = std::move(un);
        un = Matrix();
        un.clone(u);

        vnm1 = std::move(vn);
        vn = Matrix();
        vn.clone(v);

        // visualization

        deltat += dt;
        ++cnt;
    }

    return 0;
}