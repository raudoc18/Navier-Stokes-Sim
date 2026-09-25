//
// Created by Dominik on 11.09.26.
//
#include <Accelerate/Accelerate.h>

#include "CrankNicolson.h"

#include <iostream>

typedef datastruct::Matrix<double> Matrix;

Matrix CrankNicolson::AdotCN(const Matrix &m) {
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
                (*res_ptr)(j, i) = k / 2 * (m(j + 1, i) + m(j - 1, i) + m(j, i + 1) + m(j, i - 1))
                + (-2 * k - 1) * m(j,  i);
            }
        }
    });
    return res;
}

Matrix CrankNicolson::arangeRK(Matrix &x) const {
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
                auto adotcn = 0.5 * k * (x(j + 1, i) + x(j - 1, i) + x(j, i + 1) + x(j, i - 1)) + (-2.0 * k - 1) * x(j, i);
                auto b = - 0.5 * k * (x(j + 1, i) + x(j - 1, i) + x(j, i + 1) + x(j, i - 1)) + (2.0 * k - 1) * x(j, i);
                (*res_ptr)(j, i) = b - adotcn;
            }
        }
    });
    return res;
}

void CrankNicolson::CN_Wrapper() {
    CN(u);
    CN(v);
}

void CrankNicolson::CN(Matrix &m) {
    Matrix rk = arangeRK(m);

    Matrix pk = Matrix();
    pk.clone(rk);

    int cnt = 0;

    Matrix Apk = std::move(AdotCN(pk));

    double norm = rk.norm();

    while (norm > tol && cnt < maxit) {
        double rkdot = rk.dot(rk);
        double ak = rkdot / pk.dot(Apk);

        m.add(pk, ak);
        rk.subtract(Apk, ak);

        norm = rk.norm();

        if (norm < tol) {
            return;
        }

        double bk = rk.dot(rk)/rkdot;
        pk.add(rk, 1.0, bk);

        Apk = std::move(AdotCN(pk));

        ++cnt;
    }
    if (cnt > maxit) {
        std::cerr << "Warning: Did not converge in the specified timesteps! ";
    }
}
