#ifndef ADSYMATFINITET_H
#define ADSYMATFINITET_H

const int dim = 6;

const double stepsize = 1.0e-8;

const double sqrt_lambda_parameter = 596.6;		//MeV
const double lambda = sqrt_lambda_parameter*sqrt_lambda_parameter;
const double gamma_parameter = 0.9722;			//random value from earlier paper
const double k = 0.9963*sqrt_lambda_parameter;	//MeV
const double q = 0.999*k;						//MeV, k>=q

inline int solution_indexer(int iy, int iv)
{
	return ( iy * dim + iv );
}

#endif
