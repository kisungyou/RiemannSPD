#ifndef CASKET_SPD_H_
#define CASKET_SPD_H_

#include "RcppArmadillo.h"

// AUXILIARY FUNCTIONS ---------------------------------------------------------
arma::cube aux_cube_logm(arma::cube &data3d);
arma::cube aux_cube_expm(arma::cube &data3d);

#endif