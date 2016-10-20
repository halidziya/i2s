#pragma once
#include "FastMat.h"

extern  double kappa0;
extern  double kappa1;
extern  double alpha;
extern  int MAXCOMP;
extern  Normal priormean;
extern  IWishart priorcov;

Vector stickBreaker(double ustar, double betastar = 1.0, double alpha = 1);
