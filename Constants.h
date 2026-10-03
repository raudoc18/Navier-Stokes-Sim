//
// Created by Dominik on 04.09.26.
//

#ifndef CPP_CONSTANTS_H
#define CPP_CONSTANTS_H
#include <algorithm>

const int nx = 70;
const int ny = 14;

const double tol = 1e-5;
const int maxit = 10000;

const int elem_cnt = (nx + 2) * (ny + 2);

const double lx = 5; // m
const double ly = 1; // m

const double dx = lx/static_cast<double>(nx);
const double dy = ly/static_cast<double>(ny);

const double rho = 1.0;
const double mu = 0.01;
const double nu = mu/rho;

const int windowWidth = lx / ly * 400;
const int windowHeight = 400;

const double courant_number = 0.5;
const double Ut = 5.0;

const double Re = rho * Ut * lx/mu;

const double min_delta = std::min(dx, dy);

// Convective CFL
const double dt_conv = courant_number * min_delta / Ut;

// Viscous / Diffusion CFL (safety factor 0.25 for 2D explicit)
const double dt_visc = 0.25 * (dx * dx * dy * dy) / (nu * (dx * dx + dy * dy));

// True stable time step
// const double dt = std::min(dt_conv, dt_visc);
const double dt = 0.001;

const double k = dt * nu / dx /dx;

const int ntimesteps = 10;
const double animduration = 1;

const int batchsize = 32;

#endif //CPP_CONSTANTS_H
