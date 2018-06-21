#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "AdSYMatFiniteT.h"

using namespace std;

int mode = 1;	//0 - use definite range in y and terminate
				//1 - just go until the first negative result
				//for f(y) and determine y_H

void set_initial_conditions(double v[]);
int func (double y, const double v[], double f[],
	void *params);
int jac (double y, const double v[], double *dfdy, 
	double dfdt[], void *params);
void run_ode_solver(double y0, double h,
	vector<double> * ypts, vector<double> * solution);
double get_yH(vector<double> * ypts, vector<double> * solution);

int main (void)
{
	const int nypts = 101;
	vector<double> ypts, solution;
	
	//initialize
	if (mode == 0)
	{
		for (int iy = 0; iy < nypts; ++iy)
		{
			ypts.push_back(0.0 + 0.001*double(iy));
			solution.push_back(0.0);
		}
	}

	//solve
	run_ode_solver(0.0, 1.e-6, &ypts, &solution);

	//output results
	for (int iy = 0; iy < ypts.size(); ++iy)
		cout << ypts[iy] << "   " << solution[iy] << endl;

	double yH = get_yH(&ypts, &solution);
	
	cout << "y_H = " << yH << endl;

	return 0;
}

/////////////////////////////////////////////////
/////////////////////////////////////////////////
/////////////////////////////////////////////////

void set_initial_conditions(double v[])
{
	//variable order: 0-f
	v[0] = 1.0;	//f(y=0) == 1.0

	return;	
}
int func (double y, const double v[], double f[],
   void *params)
{
	//variable order: 0-f
	f[0] = 4.0*(q*q*v[0] - k*k)/q;

	return GSL_SUCCESS;
}

int jac (double y, const double v[], double *dfdy, 
  double dfdt[], void *params)
{
	double mu = *(double *)params;
	gsl_matrix_view dfdy_mat 
	 = gsl_matrix_view_array (dfdy, dim, dim);
	gsl_matrix * m = &dfdy_mat.matrix;

	//initialize matrix to zero
	for (int i = 0; i < dim; ++i)
	for (int j = 0; j < dim; ++j)
		gsl_matrix_set (m, i, j, 0.0);

	//variable order: 0-f
	gsl_matrix_set (m, 0, 0, 1.0);
	//

	for (int i = 0; i < dim; ++i)
		dfdt[i] = 0.0;

	return GSL_SUCCESS;
}


void run_ode_solver(double y0, double h,
	vector<double> * ypts, vector<double> * solution)
{
	const gsl_odeiv_step_type * T 
	 = gsl_odeiv_step_rk8pd;

	gsl_odeiv_step * s 
	 = gsl_odeiv_step_alloc (T, dim);
	gsl_odeiv_control * c 
	 = gsl_odeiv_control_y_new (1e-6, 0.0);
	gsl_odeiv_evolve * e 
	 = gsl_odeiv_evolve_alloc (dim);

	double mu = 10;
	gsl_odeiv_system sys = {func, jac, dim, &mu};

	double v[dim] = { 0.0 };

	if (mode == 0)
	{
		for (int iy = 0; iy < ypts->size(); ++iy)
		{
			double y = y0, y1 = ypts->at(iy);
			double h = 1e-10;

			//set initials conditions
			//variable order: 0-f
			set_initial_conditions(v);

			while (y < y1)
			{
				int status = gsl_odeiv_evolve_apply (e, c, s,
									                &sys, 
									                &y, y1,
									                &h, v);
				if (status != GSL_SUCCESS)
					break;
			}
			//printf ("%.5e %.5e\n", y, v[0]);
			solution->at(iy) = v[0];
		}
	}
	else if (mode == 1)
	{
		int count = 0;
		do
		{
			double y = y0, y1 = 0.0 + 0.001*double(count);
			double h = 1e-10;

			//set initials conditions
			//variable order: 0-f
			set_initial_conditions(v);

			while (y < y1)
			{
				int status = gsl_odeiv_evolve_apply (e, c, s,
									                &sys, 
									                &y, y1,
									                &h, v);
				if (status != GSL_SUCCESS)
					break;
			}
			//printf ("%.5e %.5e\n", y, v[0]);
			ypts->push_back(y1);
			solution->push_back(v[0]);
			count++;
		} while( v[0] >= 0.0 );
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	return;
}

double get_yH(vector<double> * ypts, vector<double> * solution)
{
	double yH_estimate = 0.0;
	int size = ypts->size();

	double A0 = solution->at(size-3), A1 = solution->at(size-2), A2 = solution->at(size-1);

	double x0 = ypts->at(size-3);
	double dx = ypts->at(size-2) - ypts->at(size-3);	//assumed constant

	double num1 = 2.0*(A0-2.0*A1+A2)*x0 + dx*(3.0*A0-4.0*A1+A2);
	double num2 = dx * sqrt( A0*A0
					+ (A2-4.0*A1)*(A2-4.0*A1)
					- 2.0*A0*(4.0*A1+A2) );
	double den = 2.0*(A0-2.0*A1+A2);

	yH_estimate = ( num1 + num2 ) / den;
	if ( yH_estimate < x0 || yH_estimate > x0+2.0*dx )
		yH_estimate = ( num1 - num2 ) / den;

	return ( yH_estimate );
}






