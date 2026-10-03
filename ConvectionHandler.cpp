//
// Created by Dominik on 11.09.26.
//

#include "ConvectionHandler.h"

#include <iostream>
#include <ostream>
#include <dispatch/dispatch.h>

typedef datastruct::Matrix<double> Matrix;

void ConvectionHandler::AB2u() const {
    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize + 1;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

            for (int j = start_idx; j < end_idx; j++) {
                for (int i = 1; i < nx + 1; i++) {
                    double uwn = 0.5 * (un(j, i) + un(j, i - 1));
                    double uen = 0.5 * (un(j, i) + un(j, i + 1));
                    double unn = 0.5 * (un(j, i) + un(j + 1, i));
                    double usn = 0.5 * (un(j, i) + un(j - 1, i));
                    double vnn = 0.5 * (vn(j + 1 , i - 1) + vn(j + 1, i));
                    double vsn = 0.5 * (vn(j, i) + vn(j, i - 1));

                    double uwnm1 = 0.5 * (unm1(j, i) + unm1(j, i - 1));
                    double uenm1 = 0.5 * (unm1(j, i) + unm1(j, i + 1));
                    double unnm1 = 0.5 * (unm1(j, i) + unm1(j + 1, i));
                    double usnm1 = 0.5 * (unm1(j, i) + unm1(j - 1, i));
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

double ConvectionHandler::PSI(double r) const {
    // implement SUPERBEE flux limiter
    return std::max(0.0, std::max(std::min(2.0*r, 1.0), std::min(r, 2.0)));
}

void ConvectionHandler::forwardEuler(Matrix &uCon, Matrix &vCon) {
    forwardEuleru(uCon);
    forwardEulerv(vCon);
}

void ConvectionHandler::forwardEuleru(Matrix &uCon) const {
    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize + 1;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

            for (int j = start_idx; j < end_idx; j++) {
                for (int i = 2; i < nx + 1; i++) {
                    double Fe = 0.5 * (u(j, i) + u(j, i + 1)) * dy;
                    double Fw = 0.5 * (u(j, i) + u(j, i - 1)) * dy;
                    double Fn = 0.5 * (v(j + 1, i) + v(j + 1, i - 1)) * dx;
                    double Fs = 0.5 * (v(j, i) + v(j, i - 1)) * dx;

                    double ue = 0.0;
                    double uw = 0.0;
                    double un = 0.0;
                    double us = 0.0;

                    if (Fe > 0.0) {
                        double den = u(j, i + 1) - u(j, i);
                        double r = std::abs(den) > 1e-12 ? (u(j, i) - u(j, i - 1))/den : 0.0;
                        ue = u(j, i) + 0.5 * PSI(r) * (u(j, i + 1) - u(j, i));
                    }
                    else if (Fe < 0.0){
                        if (i == nx) {
                           ue = u(j, i + 1);
                        }
                        else {
                            double den = u(j, i + 1) - u(j, i);
                            double r = std::abs(den) > 1e-12 ? (u(j, i + 2) - u(j, i + 1))/den : 0.0;
                            ue = u(j, i + 1) + 0.5 * PSI(r) * (u(j, i) - u(j, i + 1));
                        }
                    }

                    if (Fw > 0.0) {
                        if (i == 2) {
                            uw = u(j, i - 1);
                        }
                        else {
                            double den = u(j, i) - u(j, i - 1);
                            double r = std::abs(den) > 1e-12 ? (u(j, i - 1) - u(j, i - 2))/den : 0.0;
                            uw = u(j, i - 1) + 0.5 * PSI(r) * (u(j, i) - u(j, i - 1));
                        }
                    }
                    else if (Fw < 0.0){
                        double den = u(j, i) - u(j, i - 1);
                        double r = std::abs(den) > 1e-12 ? (u(j, i + 1) - u(j, i))/den : 0.0;
                        uw = u(j, i) + 0.5 * PSI(r) * (u(j, i - 1) - u(j, i));
                    }

                    if (Fn > 0.0) {
                        double den = u(j + 1, i) - u(j, i);
                        double r = std::abs(den) > 1e-12 ? (u(j, i) - u(j - 1, i))/den : 0.0;
                        un = u(j, i) + 0.5 * PSI(r) * (u(j + 1, i) - u(j, i));
                    }
                    else if (Fn < 0.0){
                        if (j == ny) {
                            un = u(j + 1, i);
                        }
                        else {
                            double den = u(j + 1, i) - u(j, i);
                            double r = std::abs(den) > 1e-12 ? (u(j + 2, i) - u(j + 1, i))/den : 0.0;
                            un = u(j + 1, i) + 0.5 * PSI(r) * (u(j, i) - u(j + 1, i));
                        }
                    }

                    if (Fs > 0.0) {
                        if (j== 1) {
                            us = u(j - 1, i);
                        }
                        else {
                            double den = u(j, i) - u(j - 1, i);
                            double r = std::abs(den) > 1e-12 ? (u(j - 1, i) - u(j - 2, i))/den : 0.0;
                            us = u(j - 1, i) + 0.5 * PSI(r) * (u(j, i) - u(j - 1, i));
                        }
                    }
                    else if (Fs < 0.0){
                        double den = u(j, i) - u(j - 1, i);
                        double r = std::abs(den) > 1e-12 ? (u(j + 1, i) - u(j, i))/den : 0.0;
                        us = u(j, i) + 0.5 * PSI(r) * (u(j - 1, i) - u(j, i));
                    }

                    uCon(j, i) = u(j, i) - dt/dx/dy * (ue * Fe - uw * Fw + un * Fn - us * Fs);
                }
            }
        }
    );
}

void ConvectionHandler::AB2v() const {

    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize + 2;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

            for (int j = start_idx; j < end_idx; j++) {
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

void ConvectionHandler::forwardEulerv(Matrix &vCon) const {
    dispatch_queue_t queue = dispatch_get_global_queue(QOS_CLASS_USER_INITIATED, 0);

    // get the number of workers
    int n_workers = (ny + batchsize - 1) / batchsize;

    dispatch_apply(n_workers, queue, ^(size_t batch_idx) {
        int start_idx = static_cast<int>(batch_idx) * batchsize + 1;
        int end_idx = std::min(start_idx + batchsize, ny + 1);

            for (int j = start_idx; j < end_idx; j++) {
                for (int i = 1; i < nx + 1; i++) {
                    double Fw = (u(j, i) + u(j - 1, i))/2 * dy;
                    double Fe = (u(j, i + 1) + u(j - 1, i + 1))/2 * dy;
                    double Fn = (v(j, i) + v(j + 1, i))/2 * dx;
                    double Fs = (v(j, i) + v(j - 1, i))/2 * dx;

                    double vw = 0.0;
                    double ve = 0.0;
                    double vn = 0.0;
                    double vs = 0.0;

                    if (Fw > 0.0) {
                        if (i == 1) {
                            vw = v(j, i - 1);
                        }
                        else {
                            double den = v(j, i) - v(j, i - 1);
                            double r = std::abs(den) > 1e-12 ? (v(j, i - 1) - v(j, i - 2))/den : 0.0;
                            vw = v(j, i - 1) + 0.5 * PSI(r) * (v(j, i) - v(j, i - 1));
                        }
                    }
                    else {
                        double den = v(j, i) - v(j, i - 1);
                        double r = std::abs(den) > 1e-12 ? (v(j, i + 1) - v(j, i))/den : 0.0;
                        vw = v(j, i) + 0.5 * PSI(r) * (v(j, i - 1) - v(j, i));
                    }

                    if (Fe > 0.0) {
                        double den = v(j, i + 1) - v(j, i);
                        double r = std::abs(den) > 1e-12 ? (v(j, i) - v(j, i - 1))/den : 0.0;
                        ve = v(j, i) + 0.5 * PSI(r) * (v(j, i + 1) - v(j, i));
                    }
                    else {
                        if (i == nx) {
                            ve = v(j, i + 1);
                        }
                        else {
                            double den = v(j, i + 1) - v(j, i);
                            double r = std::abs(den) > 1e-12 ? (v(j, i + 2) - v(j, i + 1))/den : 0.0;
                            ve = v(j, i + 1) + 0.5 * PSI(r) * (v(j, i) - v(j, i + 1));
                        }
                    }

                    if (Fn > 0.0) {
                        double den = v(j + 1, i) - v(j, i);
                        double r = std::abs(den) > 1e-12 ? (v(j, i) - v(j - 1, i))/den : 0.0;
                        vn = v(j, i) + 0.5 * PSI(r) * (v(j + 1, i) - v(j, i));
                    }
                    else {
                        if (j == ny) {
                            vn = v(j + 1, i);
                        }
                        else {
                            double den = v(j + 1, i) - v(j, i);
                            double r = std::abs(den) > 1e-12 ? (v(j + 2, i) - v(j + 1, i))/den : 0.0;
                            vn = v(j + 1, i) + 0.5 * PSI(r) * (v(j, i) - v(j + 1, i));
                        }
                    }

                    if (Fs > 0.0) {
                        if (j == 1) {
                            vs = v(j - 1, i);
                        }
                        else {
                            double den = v(j, i) - v(j - 1, i);
                            double r = std::abs(den) > 1e-12 ? (v(j - 1, i) - v(j - 2, i))/den : 0.0;
                            vs = v(j - 1, i) + 0.5 * PSI(r) * (v(j, i) - v(j - 1, i));
                        }
                    }
                    else {
                        double den = v(j, i) - v(j - 1, i);
                        double r = std::abs(den) > 1e-12 ? (v(j + 1, i) - v(j, i))/den : 0.0;
                        vs = v(j, i) + 0.5 * PSI(r) * (v(j - 1, i) - v(j, i));
                    }

                    vCon(j, i) = v(j, i) - dt/dx/dy * (ve * Fe - vw * Fw + vn * Fn - vs * Fs);
                }
            }
        }
    );
}