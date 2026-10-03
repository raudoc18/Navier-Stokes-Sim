//
// Created by Dominik on 12.09.26.
//

#include "Matrix.h"
#include "PressureHandler.h"

#include <iostream>

typedef datastruct::Matrix<double> Matrix;

Matrix PressureHandler::AdotChorin(const datastruct::Matrix<double> &x) {
    auto res = Matrix(0.0);
    auto* res_ptr = &res;

    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize + 1;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

        for (int j = start_idx; j < end_idx; j++) {
            for (int i = 1; i < nx + 1; i++) {
                (*res_ptr)(j, i) =
                    - x(j + 1, i)*an(j, i) - x(j - 1, i)*as(j, i)
                    - x(j, i + 1)*ae(j, i) - x(j, i - 1)*aw(j, i)
                    + ap(j, i) * x(j, i);
            }
        }
    });
    return res;
}

void PressureHandler::laplaceSolver(Matrix &p_new, Matrix &u, Matrix &v) {
    Matrix rk = arangeRK(u, v);

    Matrix ones = Matrix(1.0);

    rk.subtract(ones, rk.mean());

    Matrix pk = Matrix();
    pk.clone(rk);

    Matrix Apk = std::move(AdotChorin(pk));

    int cnt = 0;

    double rkdot = rk.dot(rk);

    while (std::sqrt(rkdot) > tol && cnt < maxit) {
        double ak = rkdot / pk.dot(Apk);

        p_new.add(pk, ak);
        rk.subtract(Apk, ak);

        rk.subtract(ones, rk.mean());

        double rkdot_next = rk.dot(rk);

        if (std::sqrt(rkdot_next) < tol) {
            p_new.subtract(ones, p_new.mean());
            return;
        }

        double bk = rkdot_next/rkdot;
        pk.add(rk, 1.0, bk);

        rkdot = rkdot_next;

        Apk = std::move(AdotChorin(pk));

        ++cnt;
    }
    p_new.subtract(ones, p_new.mean());
    std::cerr << "Warning: Did not converge in the specified timesteps! " << cnt << std::endl;
}

void PressureHandler::gradient(Matrix &x, Matrix &u, Matrix &v) {
    auto gradx = Matrix();
    auto* resx_ptr = &gradx;
    auto grady = Matrix();
    auto* resy_ptr = &grady;

    for (int j = 1; j < ny + 1; j++) {
        for (int i = 2; i < nx + 2; i++) {
            (*resx_ptr)(j, i) = (x(j, i) - x(j, i - 1))/dx;
        }
    }
    for (int j = 2; j < ny + 1; j++) {
        for (int i = 1; i < nx + 1; i++) {
            (*resy_ptr)(j, i) = (x(j, i) - x(j - 1, i))/dy;
        }
    }

    gradx.printMatrix();

    u.subtract(gradx, dt/rho);
    v.subtract(grady, dt/rho);

}

void PressureHandler::projection(Matrix &u, Matrix &v) {
    auto p_correction = Matrix(0.0);
    laplaceSolver(p_correction, u, v);
    gradient(p_correction, u, v);
    p.add(p_correction);
}

Matrix PressureHandler::arangeRK(Matrix &u, Matrix &v) {
    auto divu = Matrix(0.0);
    auto divv = Matrix(0.0);

    for (int i = 1; i < nx + 1; i++) {
        for (int j = 1; j < ny + 1; j++) {
            divu(j, i) = -rho/dt * ((u(j, i + 1) - u(j, i)) / dx);
        }
    }


    for (int i = 1; i < nx + 1; i++) {
        for (int j = 1; j < ny + 1; j++) {
            divv(j, i) = -rho/dt * ((v(j + 1, i) - v(j, i)) / dy);
        }
    }

    divu.add(divv);

    return divu;
}

void PressureHandler::naivSolver(Matrix &u, Matrix &v) {
    Matrix rhs = Matrix();
    Matrix divu = Matrix();
    Matrix divv = Matrix();
    Matrix ones = Matrix(1.0);

    // 1. Calculate divergence strictly for interior cells (1 to nx, 1 to ny)
    for (int i = 1; i < nx; i++) {
        for (int j = 1; j < ny + 1; j++) {
            divu(j, i) = 1.0/dt * ((u(j, i + 1) - u(j, i)) / dx);
        }
    }


    for (int i = 1; i < nx + 1; i++) {
        for (int j = 2; j < ny + 1; j++) {
            divv(j, i) = 1.0/dt * ((v(j, i) - v(j - 1, i)) / dy);
        }
    }

    divu.add(divv);
    rhs = std::move(divu);

    int cnt = 0;
    Matrix pprev = Matrix();
    Matrix pnext = Matrix();
    Matrix diff = Matrix();

    double norm = 1.0;

    while (norm > tol && cnt < maxit) {
        // 2. Jacobi relaxation on internal cells
        for (int i = 1; i <= nx; i++) {
            for (int j = 1; j <= ny; j++) {
                pnext(j, i) = 0.25 * (pprev(j, i + 1) + pprev(j, i - 1) +
                                      pprev(j + 1, i) + pprev(j - 1, i) -
                                      dx * dy * rhs(j, i));
            }
        }

        // 3. Boundary conditions
        // Top and bottom walls (Neumann: dp/dy = 0)
        for (int i = 1; i <= nx; i++) {
            pnext(0, i)      = pnext(1, i);
            pnext(ny + 1, i) = pnext(ny, i);
        }

        // Inlet (Neumann: dp/dx = 0) & Outlet (Dirichlet: p = 0 at face)
        for (int j = 1; j <= ny; j++) {
            pnext(j, 0)      = pnext(j, 1);       // Inlet
            pnext(j, nx + 1) = -pnext(j, nx);      // Outlet: (p_ghost + p_int)/2 = 0
        }

        diff.clone(pnext);
        diff.subtract(pprev);
        norm = diff.norm();

        pprev.subtract(ones, pprev.mean());

        pprev = std::move(pnext);
        pnext = Matrix();
        cnt++;
    }

    std::cout << "Iterations: " << cnt << " | Final Norm: " << norm << std::endl;
    p = std::move(pprev);
}

void PressureHandler::initBoundaries() {
    for (int i = 1; i < nx + 1; i++) {
        for (int j = 1; j < ny + 1; j++) {
            aw(j, i) = 1.0/dx/dx;
            ae(j, i) = 1.0/dx/dx;
            an(j, i) = 1.0/dy/dy;
            as(j, i) = 1.0/dy/dy;

            if (i == 1) {
                aw(j, i) = 0.0;
            }
            else if (i == nx) {
                ae(j, i) = 0.0;
                //ap(j, i) += 1.0/dx/dx;
            }
            if (j == 1) {
                as(j, i) = 0.0;
            }
            else if (j == ny) {
                an(j, i) = 0.0;
            }

            ap(j, i) += aw(j, i) + ae(j, i) + an(j, i) + as(j, i);
        }

    }
}