#include <stdio.h>
#include <cmath>
#include <iostream>
#include <vector>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv.h>

#include "AdSYMatFiniteT.h"
#include "Potentials.h"

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
double get_yH( vector<double> * ypts,
				vector<double> * solution,
				double * f_prime_at_yH_estimate );

int main (void)
{
	const int nypts = 101;
	vector<double> ypts, solution;
	
	//initialize
	if (mode == 0)
	{
		for (int iy = 0; iy < nypts; ++iy)
		{
			ypts.push_back(0.0 + stepsize*double(iy));
			for (int iv = 0; iv < dim; ++iv)
				solution.push_back(0.0);
		}
	}

	//solve
	double y0 = 0.0;
	double h = 1.e-10;
	run_ode_solver(y0, h, &ypts, &solution);

	//output results
	for (int iy = 0; iy < ypts.size(); ++iy)
	{
		cout << ypts[iy];
		for (int iv = 0; iv < dim; ++iv)
			cout << "   " << solution[solution_indexer(iy, iv)];
		cout << endl;
	}

	if (mode == 1)
	{
		double f_prime_at_yH = 0.0;
		double yH = get_yH(&ypts, &solution, &f_prime_at_yH);
	
		cout << "y_H = " << yH << endl;
		//cout << "z_H = " << sgn << endl;
		cout << "f'(y_H) = " << f_prime_at_yH << endl;
	}

	return 0;
}

///////////////////////////////////////////////////////////

int func (double y, const double v[], double f[],
   void *params)
{
	double Bpq = v[2]+q;
	double Vt = Vtilde(v[0], v[1], y);
	double dVtdphi = dVtildedphi(v[0], v[1], y);
	double dVtdG = dVtildedG(v[0], v[1], y);

	//variable order: 0-phi, 1-G, 2-B, 3-f, 4-psi, 5-gamma
	f[0] = v[4];
	f[1] = v[5];
	f[2] = ( v[4]*v[4] + v[5]*v[5] )/6.0;
	f[3] = ( Vt - 12.0*k*k + 12.0*v[3]*Bpq*Bpq - 3.0*v[3]*f[2] ) / ( 3.0*Bpq );
	f[4] = ( dVtdphi - f[3]*v[4] + 4.0*v[3]*Bpq*v[4] ) / v[3];
	f[5] = ( dVtdG - f[3]*v[5] + 4.0*v[3]*Bpq*v[5] ) / v[3];

	return GSL_SUCCESS;
}

int jac (double y, const double v[], double *dfdy, 
  double dfdt[], void *params)
{
	double Bpq = v[2]+q;
	double Vt = Vtilde(v[0], v[1], y);
	double dVtdphi = dVtildedphi(v[0], v[1], y);
	double dVtdG = dVtildedG(v[0], v[1], y);
	double mu = *(double *)params;
	gsl_matrix_view dfdy_mat 
	 = gsl_matrix_view_array (dfdy, dim, dim);
	gsl_matrix * m = &dfdy_mat.matrix;

	//initialize matrix to zero
	for (int i = 0; i < dim; ++i)
	for (int j = 0; j < dim; ++j)
		gsl_matrix_set (m, i, j, 0.0);

	//variable order: 0-phi, 1-G, 2-B, 3-f, 4-psi, 5-gamma
	gsl_matrix_set (m, 0, 4, 1.0);
	//
	gsl_matrix_set (m, 1, 5, 1.0);
	//
	gsl_matrix_set (m, 2, 4, v[4]/3.0);
	gsl_matrix_set (m, 2, 5, v[5]/3.0);
	//
	gsl_matrix_set (m, 3, 2, ( 24.0*k*k-2.0*Vt+v[3]*( 24.0*Bpq*Bpq+v[4]*v[4]+v[5]*v[5] ) )
								/ (6.0*Bpq*Bpq));
	gsl_matrix_set (m, 3, 3, 4.0*Bpq-(v[4]*v[4]+v[5]*v[5])/(6.0*Bpq));
	gsl_matrix_set (m, 3, 4, -3.0*v[3]*v[4]/(3.0*Bpq));
	gsl_matrix_set (m, 3, 5, -3.0*v[3]*v[5]/(3.0*Bpq));
	//
	gsl_matrix_set (m, 4, 2, -( 24.0*k*k-2.0*Vt+v[3]*( v[4]*v[4]+v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq*Bpq));
	gsl_matrix_set (m, 4, 3, ( v[4]*( Vt - 12.0*k*k ) - 3.0*dVtdphi*Bpq )
								/ (3.0*v[3]*v[3]*Bpq));
	gsl_matrix_set (m, 4, 4, ( 24.0*k*k-2.0*Vt+v[3]*( 3.0*v[4]*v[4]+v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq));
	gsl_matrix_set (m, 4, 5, v[4]*v[5] / (3.0*Bpq));
	//
	gsl_matrix_set (m, 5, 2, -( 24.0*k*k-2.0*Vt+v[3]*( v[4]*v[4]+v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq*Bpq));
	gsl_matrix_set (m, 5, 3, ( v[5]*( Vt - 12.0*k*k ) - 3.0*dVtdG*Bpq )
								/ (3.0*v[3]*v[3]*Bpq));
	gsl_matrix_set (m, 5, 4, ( 24.0*k*k-2.0*Vt+v[3]*( v[4]*v[4]+3.0*v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq));
	gsl_matrix_set (m, 5, 5, v[4]*v[5] / (3.0*Bpq));

	for (int i = 0; i < dim; ++i)
		dfdt[i] = 0.0;

	return GSL_SUCCESS;
}


void set_initial_conditions(double v[])
{
	//OLD INITIAL CONDITIONS
	/*
	//variable order: 0-phi, 1-G, 2-B, 3-f, 4-psi, 5-gamma
	v[0] = 0.0;							//phi(y=0) == 0.0
	v[1] = 0.0;							//G(y=0) == 0.0
	v[2] = W(v[0], v[1], 0.0);
	v[3] = 1.0;							//f(y=0) == 1.0, by definition
	v[4] = 6.0 * dWdphi(v[0], v[1], 0.0);
	v[5] = 6.0 * dWdG(v[0], v[1], 0.0);
	*/
	
	//NEW INITIAL CONDITIONS
	//variable order: 0-phi, 1-G, 2-B, 3-f, 4-psi, 5-gamma
	v[0] = sqrt(8.0/3.0)*lambda/(k*k);	//phi(y=0) == sqrt(8./3.)*lambda/k^2
	v[1] = sqrt(24.0*lambda)/k;			//G(y=0) == sqrt(24.*lambda)/k
	v[2] = W(v[0], v[1], 0.0);			//B = W
	v[3] = 1.0;							//f(y=0) == 1.0, by definition
	v[4] = 2.0*exp(a*v[0])*sqrt_lambda_parameter/k;
	v[5] = sqrt(24.0)*exp(a*v[0]);

	//for (int i = 0; i < dim; ++i)
	//	cout << "v[" << i << "] = " << v[i] << endl;

	return;	
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

	double v[dim] = { 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 };

	if (mode == 0)
	{
		for (int iy = 0; iy < ypts->size(); ++iy)
		{
			double y = y0, y1 = ypts->at(iy);

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

				//printf ("%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
				//		y, v[0], v[1], v[2], v[3], v[4], v[5]);
			}
			//printf ("%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
			//		y, v[0], v[1], v[2], v[3], v[4], v[5]);

			for (int iv = 0; iv < dim; ++iv)
				solution->at(solution_indexer(iy, iv)) = v[iv];
		}
	}
	else if (mode == 1)
	{
		int count = 0;
		do
		{
			double y = y0, y1 = 0.0 + stepsize*double(count);

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

				//printf ("%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
				//		y, v[0], v[1], v[2], v[3], v[4], v[5]);
			}
			//printf ("%.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
			//		y, v[0], v[1], v[2], v[3], v[4], v[5]);
			ypts->push_back(y1);
			for (int iv = 0; iv < dim; ++iv)
				solution->push_back(v[iv]);
			count++;
		} while( v[3] >= 0.0 );	//termination condition is f(y>yH) < 0.0
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	return;
}

/*double get_yH( vector<double> * ypts,
				vector<double> * solution,
				vector<double> * interpolated_results )
{
	double yH_estimate = 0.0;
	int size = ypts->size();

	double A0 = solution->at(solution_indexer(size-3, 3));
	double A1 = solution->at(solution_indexer(size-2, 3));
	double A2 = solution->at(solution_indexer(size-1, 3));	//this is probably nan

	double x0 = ypts->at(size-3);
	double x1 = ypts->at(size-2);
	double dx = x1 - x0;			//assumed constant

	double num1 = 2.0*(A0-2.0*A1+A2)*x0 + dx*(3.0*A0-4.0*A1+A2);
	double num2 = dx * sqrt( A0*A0
					+ (A2-4.0*A1)*(A2-4.0*A1)
					- 2.0*A0*(4.0*A1+A2) );
	double den = 2.0*(A0-2.0*A1+A2);

	yH_estimate = ( num1 + num2 ) / den;
	if ( yH_estimate < x0 || yH_estimate > x0+2.0*dx )
		yH_estimate = ( num1 - num2 ) / den;

	//just do linear interpolation for now
	for (int iv = 0; iv < dim; ++iv)
	{
		int sol1 = solution->at(solution_indexer(size-2, iv));
		int sol2 = solution->at(solution_indexer(size-1, iv));
		interpolated_results->at(iv) = sol1 + (yH_estimate - x1)*(sol2-sol1)/dx;
	}

	return ( yH_estimate );
}*/

double get_yH( vector<double> * ypts,
				vector<double> * solution,
				double * f_prime_at_yH_estimate )
{
	double yH_estimate = 0.0;
	int size = ypts->size();

	double A0 = solution->at(solution_indexer(size-4, 3));
	double A1 = solution->at(solution_indexer(size-3, 3));
	double A2 = solution->at(solution_indexer(size-2, 3));

	double x0 = ypts->at(size-4);
	double x1 = ypts->at(size-3);
	double dx = x1 - x0;			//assumed constant

	double num1 = 2.0*(A0-2.0*A1+A2)*x0 + dx*(3.0*A0-4.0*A1+A2);
	double num2 = dx * sqrt( A0*A0
					+ (A2-4.0*A1)*(A2-4.0*A1)
					- 2.0*A0*(4.0*A1+A2) );
	double den = 2.0*(A0-2.0*A1+A2);

	yH_estimate = ( num1 + num2 ) / den;		//choose one quadratic root
	if ( yH_estimate < x0 || yH_estimate > ypts->at(size-1) )
		yH_estimate = ( num1 - num2 ) / den;	//or the other quadratic root

	//finally, estimate f'(y_H) < 0
	(*f_prime_at_yH_estimate)
		= -0.5 * sqrt( A0*A0 + (A2-4.0*A1)*(A2-4.0*A1) - 2.0*A0*(4.0*A1+A2) ) / dx;

	return ( yH_estimate );
}
