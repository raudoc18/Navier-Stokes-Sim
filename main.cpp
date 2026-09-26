#include "Constants.h"
#include "Matrix.h"
#include <iostream>
#include <thread>

#include "AdamsBashforth.h"
#include "Chorin.h"
#include "CrankNicolson.h"
#include "ImmersedBoundary.h"
#include "Renderer.h"
#include "StreamFunction.h"

typedef datastruct::Matrix<double> Matrix;

void boundaryConditions(Matrix &u, Matrix&v) {
    for (int i = 0; i < nx + 2; ++i) {
        for (int j = 0; j < ny + 2; ++j) {
            if (j == ny + 1) {
                u(j, i) = -u(j - 1, i);
            }

            if (j == 0) {
                u(j, i) = -u(j + 1, i);
            }

            if (i == 1) {
                u(j, i) = Ut;
            }
            if (i == nx + 1) {
                u(j, i) = u(j, i - 1);
            }

            if (i == 0) {
                v(j, i) = -v(j, i + 1);
            }
            if (i == nx + 1) {
                v(j, i) = v(j, i - 1);
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

void boundaryConditionsAfterProjection(Matrix &u, Matrix&v) {
    for (int i = 1; i < ny + 1; ++i) {
        // u(i, 1) = Ut;
        // u(i, nx + 1) = u(i, nx);
    }
}

int main() {
    int cnt = 0;
    double deltat = 0;

    auto u = Matrix(0.0);
    auto ut = Matrix(0.0);
    auto unm1 = Matrix(0.0);
    auto uAB = Matrix(0.0);
    auto uCN = Matrix(0.0);

    auto v = Matrix(0.0);
    auto vt = Matrix(0.0);
    auto vnm1 = Matrix(0.0);
    auto vAB = Matrix(0.0);
    auto vCN = Matrix(0.0);

    auto p = Matrix(0.0);
    auto pn = Matrix(0.0);

    Obstacle circle = Obstacle(0.5, 0.5);

    std::vector obstacles = { &circle };

    AdamsBashforth adm = AdamsBashforth(uAB, u, unm1, vAB, v, vnm1);
    CrankNicolson cn = CrankNicolson(u, v);
    Chorin cho = Chorin(p);
    StreamFunction stream = StreamFunction(u);
    ImmersedBoundary imb = ImmersedBoundary(obstacles, u, v, ut, vt);

    Renderer ren = Renderer(u, v, obstacles);

    std::cout << "dt: " << dt << std::endl;
    std::cout << "batch size: " << batchsize << std::endl;
    std::cout << "nx: " << nx << std::endl;
    std::cout << "ny: " << ny << std::endl;
    std::cout << "Re: " << Re << std::endl;

    while (!glfwWindowShouldClose(ren.getWindow())) {
        while (true) {
            boundaryConditions(u, v);
            boundaryConditions(unm1, vnm1);

            p = Matrix(0.0);
            unm1.clone(u);
            vnm1.clone(v);

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

            ut.clone(u);
            vt.clone(v);

            auto start = std::chrono::steady_clock::now();
            cho.projection(ut, vt);
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> duration_ms = end - start;
            std::cout << "Execution time of projection: " << duration_ms.count() << " ms\n";

            imb.computeForceTerms();

            start = std::chrono::steady_clock::now();
            p = Matrix();
            cho.projection(u, v);
            end = std::chrono::steady_clock::now();
            duration_ms = end - start;
            std::cout << "Execution time of second projection: " << duration_ms.count() << " ms\n";

            boundaryConditionsAfterProjection(u, v);

            // visualization
            if (deltat > 1.0/30.0) {
                deltat = 0;
                ren.render();
            }
            ren.render();
            deltat += dt;
        }
    }

    return 0;
}