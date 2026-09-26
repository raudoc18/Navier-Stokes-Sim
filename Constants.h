//
// Created by Dominik on 04.09.26.
//

#ifndef CPP_CONSTANTS_H
#define CPP_CONSTANTS_H

const int nx = 64;
const int ny = 64;

const double tol = 1e-5;
const int maxit = 10000;

const int elem_cnt = (nx + 2) * (ny + 2);

const double lx = 3; // m
const double ly = 1; // m

const double dx = lx/static_cast<double>(nx);
const double dy = ly/static_cast<double>(ny);

const double rho = 1.0;
const double mu = 0.1;
const double nu = mu/rho;

const int windowWidth = lx / ly * 400;
const int windowHeight = 400;

const double courant_number = 0.5;
const double Ut = 20.0;

const double Re = rho * Ut * lx/mu;

const double dt = Ut != 0.0 ? courant_number * dx / Ut : courant_number * dx;

const double k = dt * nu / dx /dx;

const int ntimesteps = 10;
const double animduration = 1;

const int batchsize = ny / 8;

#endif //CPP_CONSTANTS_H
