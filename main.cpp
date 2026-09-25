#include "Constants.h"
#include "Matrix.h"
#include <iostream>
#include <thread>

#include "AdamsBashforth.h"
#include "Chorin.h"
#include "CrankNicolson.h"
#include "Renderer.h"

typedef datastruct::Matrix<double> Matrix;

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

int main() {
    int cnt = 0;
    double deltat = 0;

    auto u = Matrix(0.0);
    auto un = Matrix(0.0);
    auto unm1 = Matrix(0.0);
    auto uAB = Matrix(0.0);
    auto uCN = Matrix(0.0);

    auto v = Matrix(0.0);
    auto vn = Matrix(0.0);
    auto vnm1 = Matrix(0.0);
    auto vAB = Matrix(0.0);
    auto vCN = Matrix(0.0);

    auto p = Matrix(0.0);
    auto pn = Matrix(0.0);

    AdamsBashforth adm = AdamsBashforth(uAB, u, unm1, vAB, v, vnm1);
    CrankNicolson cn = CrankNicolson(u, v);
    Chorin cho = Chorin(p, u, v);

    Renderer ren = Renderer(u, v);

    std::cout << "dt: " << dt << std::endl;
    std::cout << "batch size: " << batchsize << std::endl;
    std::cout << "nx: " << nx << std::endl;
    std::cout << "ny: " << ny << std::endl;

    while (!glfwWindowShouldClose(ren.getWindow())) {
        while (cnt < ntimesteps && deltat < animduration) {
            boundaryConditions(u, v);
            boundaryConditions(unm1, vnm1);

            p = Matrix(0.0);

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

            auto start = std::chrono::steady_clock::now();
            cho.projection();
            auto end = std::chrono::steady_clock::now();
            std::chrono::duration<double, std::milli> duration_ms = end - start;
            std::cout << "Execution time of boundaries: " << duration_ms.count() << " ms\n";

            unm1.clone(u);
            vnm1.clone(v);

            // visualization
            ren.render();


            std::cout << cnt << std::endl;
        }
    }

    return 0;
}