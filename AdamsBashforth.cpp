//
// Created by Dominik on 11.09.26.
//

#include "AdamsBashforth.h"

#include <iostream>
#include <ostream>
#include <dispatch/dispatch.h>

void AdamsBashforth::AB2u() const {
    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + 2 + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

            for (int j = start_idx + 1; j < end_idx; j++) {
                for (int i = 1; i < nx + 1; i++) {
                    double uwn = 0.5 * (un(j, i) + un(j, i - 1));
                    double uen = 0.5 * (un(j, i) + un(j, i + 1));
                    double unn = 0.5 * (un(j, i) + un(j + 1, i));
                    double usn = 0.5 * (un(j, i) + un(j - 1, i - 1));
                    double vnn = 0.5 * (vn(j + 1 , i - 1) + vn(j + 1, i));
                    double vsn = 0.5 * (vn(j, i) + vn(j, i - 1));

                    double uwnm1 = 0.5 * (unm1(j, i) + unm1(j, i - 1));
                    double uenm1 = 0.5 * (unm1(j, i) + unm1(j, i + 1));
                    double unnm1 = 0.5 * (unm1(j, i) + unm1(j + 1, i));
                    double usnm1 = 0.5 * (unm1(j, i) + unm1(j - 1, i - 1));
                    double vnnm1 = 0.5 * (vnm1(j + 1 , i - 1) + vnm1(j + 1, i));
                    double vsnm1 = 0.5 * (vnm1(j, i) + vnm1(j, i - 1));

                    u(j, i) = 3.0 * dt/8.0 * (
                        (std::pow(uen, 2) - std::pow(uwn, 2))/dx +
                        (unn * vnn - usn * vsn)/dy
                        )
                        + dt/8.0 * (
                            (std::pow(uenm1, 2) - std::pow(uwnm1, 2))/dx +
                            (unnm1 * vnnm1 - usnm1 * vsnm1)/dy
                            );
            }
            }
        }
    );
}

void AdamsBashforth::forwardEuleru() const {
    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + 2 + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

            for (int j = start_idx + 1; j < end_idx; j++) {
                for (int i = 1; i < nx + 1; i++) {
                double uwn = 0.5 * (un(j, i) + un(j, i - 1));
                double uen = 0.5 * (un(j, i) + un(j, i + 1));
                double unn = 0.5 * (un(j, i) + un(j + 1, i));
                double usn = 0.5 * (un(j, i) + un(j - 1, i - 1));
                double vnn = 0.5 * (vn(j + 1 , i - 1) + vn(j + 1, i));
                double vsn = 0.5 * (vn(j, i) + vn(j, i - 1));

                u(j, i) = dt/4.0 * (
                    (std::pow(uen, 2) - std::pow(uwn, 2))/dx +
                    (unn * vnn - usn * vsn)/dy);
            }
            }
        }
    );
}

void AdamsBashforth::AB2v() const {

    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + 2 + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

            for (int j = start_idx + 1; j < end_idx; j++) {
                for (int i = 1; i < nx + 1; i++) {
                    double unn = 0.5 * (un(j, i) + un(j - 1, i));
                    double usn = 0.5 * (un(j, i + 1) + un(j - 1, i + 1));
                    double ven = 0.5 * (vn(j, i) + vn(j, i + 1));
                    double vwn = 0.5 * (vn(j, i) + vn(j, i - 1));
                    double vnn = 0.5 * (vn(j, i) + vn(j + 1, i));
                    double vsn = 0.5 * (vn(j, i) + vn(j - 1, i));

                    double unnm1 = 0.5 * (unm1(j, i) + unm1(j - 1, i));
                    double usnm1 = 0.5 * (unm1(j, i + 1) + unm1(j - 1, i + 1));
                    double venm1 = 0.5 * (vnm1(j, i) + vnm1(j, i + 1));
                    double vwnm1 = 0.5 * (vnm1(j, i) + vnm1(j, i - 1));
                    double vnnm1 = 0.5 * (vnm1(j, i) + vnm1(j + 1, i));
                    double vsnm1 = 0.5 * (vnm1(j, i) + vnm1(j - 1, i));

                    v(j, i) = 3.0 * dt/8.0 * (
                        (std::pow(ven, 2) - std::pow(vwn, 2))/dx +
                        (unn * vnn - usn * vsn)/dy
                        )
                        + dt/8.0 * (
                            (std::pow(venm1, 2) - std::pow(vwnm1, 2))/dx +
                            (unnm1 * vnnm1 - usnm1 * vsnm1)/dy
                            );
            }
            }
        }
    );
}

void AdamsBashforth::forwardEulerv() const {
    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + 2 + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

            for (int j = start_idx + 1; j < end_idx; j++) {
                for (int i = 1; i < nx + 1; i++) {
                    double unn = 0.5 * (un(j, i) + un(j - 1, i));
                        double usn = 0.5 * (un(j, i + 1) + un(j - 1, i + 1));
                        double ven = 0.5 * (vn(j, i) + vn(j, i + 1));
                        double vwn = 0.5 * (vn(j, i) + vn(j, i - 1));
                        double vnn = 0.5 * (vn(j, i) + vn(j + 1, i));
                        double vsn = 0.5 * (vn(j, i) + vn(j - 1, i));

                        v(j, i) = dt/4.0 * (
                            (std::pow(ven, 2) - std::pow(vwn, 2))/dx +
                            (unn * vnn - usn * vsn)/dy
                            );
            }
            }
        }
    );
}
