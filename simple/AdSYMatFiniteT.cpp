#include <stdio.h>
#include <iostream>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "AdSYMatFiniteT.h"

using namespace std;

void set_initial_conditions(double v[]);
int func (double y, const double v[], double f[],
	void *params);
int jac (double y, const double v[], double *dfdy, 
	double dfdt[], void *params);
void run_ode_solver(double y0, double h, int nypts,
	const double ypts[], double solution[]);

int main (void)
{
	const int nypts = 101;
	double ypts[nypts], solution[nypts];
	
	//initialize
	for (int iy = 0; iy < nypts; ++iy)
	{
		ypts[iy] = 0.0 + 0.001*double(iy);
		solution[iy] = 0.0;
	}

	//solve
	run_ode_solver(0.0, 1.e-6, nypts, ypts, solution);

	//output results
	for (int iy = 0; iy < nypts; ++iy)
		cout << ypts[iy] << "   " << solution[iy] << endl;

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


void run_ode_solver(double y0, double h, int nypts,
	const double ypts[], double solution[])
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

	for (int iy = 0; iy < nypts; ++iy)
	{
		double y = y0, y1 = ypts[iy];
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
		solution[iy] = v[0];
	}


	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	return;
}
