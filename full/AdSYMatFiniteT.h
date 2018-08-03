#ifndef ADSYMATFINITET_H
#define ADSYMATFINITET_H

const int dim = 6;	//number of variables
const int mode = 1;	//0 - use definite range in y and terminate
					//1 - just go until the first negative result
					//for f(y) and determine y_H


const double stepsize = (mode==0) ? 1.0 : 1.e-6;

//const double sqrt_lambda_parameter = 596.6;		//MeV
const double sqrt_lambda_parameter = 638.5;		//MeV
const double lambda = sqrt_lambda_parameter*sqrt_lambda_parameter;
const double gamma_parameter = 0.9722;			//random value from earlier paper
//const double k = 0.9963*sqrt_lambda_parameter;	//MeV
const double k = 1.2383*sqrt_lambda_parameter;	//MeV
const double q = 1.0*k;						//MeV, k>=q

inline int solution_indexer(int iy, int iv)
{
	return ( iy * dim + iv );
}


#endif
