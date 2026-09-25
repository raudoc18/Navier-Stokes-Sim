//
// Created by Dominik on 12.09.26.
//

#include "Matrix.h"
#include "Chorin.h"

#include <iostream>

typedef datastruct::Matrix<double> Matrix;

Matrix Chorin::AdotChorin(const datastruct::Matrix<double> &x) {
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
                    - (x(j + 1, i)*an(j, i) + x(j - 1, i)*as(j, i))
                    - (x(j, i + 1)*ae(j, i) + x(j, i - 1)*aw(j, i))
                    + ap(j, i) * x(j, i);
            }
        }
    });
    return res;
}

void Chorin::laplaceSolver() {
    Matrix rk = arangeRK();

    Matrix ones = Matrix(1.0);

    rk.subtract(ones, rk.mean());

    Matrix pk = Matrix();
    pk.clone(rk);

    Matrix Apk = std::move(AdotChorin(pk));

    int cnt = 0;

    double rkdot = rk.dot(rk);

    while (std::sqrt(rkdot) > tol && cnt < maxit) {
        double ak = rkdot / pk.dot(Apk);

        p.add(pk, ak);
        rk.subtract(Apk, ak);

        rk.subtract(ones, rk.mean());

        double rkdot_next = rk.dot(rk);

        if (std::sqrt(rkdot_next) < tol) {
            p.subtract(ones, p.mean());
            return;
        }

        double bk = rkdot_next/rkdot;
        pk.add(rk, 1.0, bk);

        rkdot = rkdot_next;

        Apk = std::move(AdotChorin(pk));

        ++cnt;
    }
    p.subtract(ones, p.mean());
    std::cerr << "Warning: Did not converge in the specified timesteps! " << cnt << std::endl;
}

void Chorin::gradient(Matrix &x) {
    auto gradx = Matrix();
    auto* resx_ptr = &gradx;
    auto grady = Matrix();
    auto* resy_ptr = &grady;

    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize + 1;
        int end_idx = std::min(start_idx + batchsize, ny + 1);
        for (int j = start_idx; j < end_idx; j++) {
            for (int i = 2; i < nx + 1; i++) {
                (*resx_ptr)(j, i) = (x(j, i) - x(j, i - 1))/dx;
            }
        }
        for (int j = start_idx + 1; j < end_idx; j++) {
            for (int i = 1; i < nx + 1; i++) {
                (*resy_ptr)(j, i) = (x(j, i) - x(j - 1, i))/dy;
            }
        }
    });

    ut.subtract(gradx, dt/rho);
    vt.subtract(grady, dt/rho);

}

void Chorin::projection() {
    laplaceSolver();
    gradient(p);
}

Matrix Chorin::arangeRK() {
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
                (*res_ptr)(j, i) = - rho/dt * ((ut(j, i + 1) - ut(j, i))/dx + (vt(j + 1, i) - vt(j, i))/dy)
                    + (- (p(j + 1, i)*an(j, i) - p(j - 1, i)*as(j, i))
                    - (p(j, i + 1)*ae(j, i) - p(j, i - 1)*aw(j, i))
                    + ap(j, i) * p(j, i));
            }
        }
    });

    return res;
}

void Chorin::initBoundaries() {
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
            }
            if (j == 1) {
                as(j, i) = 0.0;
            }
            else if (j == ny) {
                an(j, i) = 0.0;
            }

            ap(j, i) = aw(j, i) + ae(j, i) + an(j, i) + as(j, i);
        }

    }
}