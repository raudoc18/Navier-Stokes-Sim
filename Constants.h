//
// Created by Dominik on 04.09.26.
//

#ifndef CPP_CONSTANTS_H
#define CPP_CONSTANTS_H

const int nx = 32;
const int ny = 32;

const int elem_cnt = (nx + 2) * (ny + 2);

const double lx = 10; // m
const double ly = 10; // m

const double dx = lx/static_cast<double>(nx);
const double dy = ly/static_cast<double>(ny);

const double rho = 1.0;
const double mu = 0.1;
const double nu = mu/rho;

const double courant_number = 0.5;
const double Ut = 1.0;

const double dt = Ut != 0.0 ? courant_number * dx / Ut : courant_number * dx;

const double k = dt * nu / dx /dx;

const int ntimesteps = 10;
const double animduration = 1;

const int batchsize = 64;

#endif //CPP_CONSTANTS_H
