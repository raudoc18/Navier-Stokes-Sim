#include "Constants.h"
#include "Matrix.h"
#include <iostream>
#include <thread>

#include "ConvectionHandler.h"
#include "PressureHandler.h"
#include "DiffusionHandler.h"
#include "ImmersedBoundaryHandler.h"
#include "RenderingHandler.h"
#include "StreamFunction.h"

typedef datastruct::Matrix<double> Matrix;

void initialConditions(Matrix &u, Matrix&v) {
    for (int i = 0; i < nx + 2; ++i) {
        for (int j = 0; j < ny + 2; ++j) {
            // if inside
            if (j != 0 && j != ny + 1) {
                u(j, i) = Ut;
            }
        }
    }
}

void boundaryConditions(Matrix &u, Matrix&v) {
    for (int i = 0; i < nx + 2; ++i) {
        for (int j = 0; j < ny + 2; ++j) {
            // inner y points
            if (j != 0 && j != ny + 1) {
                // u boundary conditions
                // left side inlet
                if (i == 1) {
                    u(j, i) = Ut;
                }
                // right side outlet
                else if (i == nx + 1) {
                    u(j, i) = u(j, i - 1);
                }
                // left side inlet but v = 0
                if (i == 0) {
                    v(j, i) = -v(j, i + 1);
                }
                // right side outlet
                else if (i == nx + 1) {
                    v(j, i) = v(j, i - 1);
                }
            }
            // inner x
            if (i != 0 && i != nx + 1) {
                // upper side
                if (j == 0) {
                    u(j, i) = -u(j + 1, i);
                }
                // lower side
                else if (j == ny + 1) {
                    u(j, i) = -u(j - 1, i);
                }
                // upper side
                if (j == 1) {
                    v(j, i) = 0;
                }
                // lower side
                else if (j == ny + 1) {
                    v(j, i) = 0;
                }
            }

        }
    }
}

int main() {
    int cnt = 0;
    double deltat = 0;

    auto u = Matrix(0.0);
    auto ut = Matrix(0.0);
    auto unm1 = Matrix(0.0);
    auto uCon = Matrix(0.0);
    auto uDiff = Matrix(0.0);

    auto v = Matrix(0.0);
    auto vt = Matrix(0.0);
    auto vnm1 = Matrix(0.0);
    auto vCon = Matrix(0.0);
    auto vDiff = Matrix(0.0);

    auto f = Matrix(0.0);

    auto p = Matrix(0.0);
    auto pn = Matrix(0.0);

    Obstacle circle = Obstacle(0.5, 0.5);

    std::vector obstacles = { &circle };

    ConvectionHandler convection = ConvectionHandler(u, u, unm1, v, v, vnm1);
    DiffusionHandler diffusion = DiffusionHandler(u, v);
    PressureHandler pressure = PressureHandler(p);
    ImmersedBoundaryHandler imb = ImmersedBoundaryHandler(obstacles, u, v, ut, vt, f);

    RenderingHandler ren = RenderingHandler(u, v, obstacles, p);

    std::cout << "dt: " << dt << std::endl;
    std::cout << "batch size: " << batchsize << std::endl;
    std::cout << "nx: " << nx << std::endl;
    std::cout << "ny: " << ny << std::endl;
    std::cout << "Re: " << Re << std::endl;
    u.collectAsVector();
    initialConditions(u, v);

    while (!glfwWindowShouldClose(ren.getWindow())) {
        while (true) {
            boundaryConditions(u, v);

            uCon = Matrix(0.0);
            vCon = Matrix(0.0);

            convection.forwardEuler(uCon, vCon);

            diffusion.explicitDiffusion(uDiff, vDiff);

            uCon.add(uDiff);
            vCon.add(vDiff);

            u = std::move(uCon);
            v = std::move(vCon);

            pressure.gradient(p, u, v);

            boundaryConditions(u, v);

            pressure.projection(u, v);

            Matrix divu = Matrix(0.0);
            Matrix divv = Matrix(0.0);

            for (int i = 1; i < nx; i++) {
                for (int j = 1; j < ny + 1; j++) {
                    divu(j, i) = ((u(j, i + 1) - u(j, i)) / dx);
                }
            }


            for (int i = 1; i < nx + 1; i++) {
                for (int j = 1; j < ny + 1; j++) {
                    divv(j, i) = ((v(j + 1, i) - v(j, i)) / dy);
                }
            }

            divu.add(divv);

            std::cout << "Maximum divergence: " << divu.data()[divu.maxIdx()] << std::endl;

            // visualization
            // if (deltat > 1.0/30.0) {
            //     deltat = 0;
            //     ren.render();
            // }

            ren.render();

            double sum_inlet = 0.0;
            double sum_outlet = 0.0;
            for (int j = 1; j < ny + 1; ++j) {
                sum_inlet += u(j, 2);
                sum_outlet += u(j, nx + 1);
            }

            std::cout << "Mass Difference: " << sum_inlet - sum_outlet << std::endl;
            deltat += dt;
            cnt++;
        }
    }

    return 0;
}