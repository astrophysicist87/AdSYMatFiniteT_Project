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

//prototypes
void set_initial_conditions(double v[], double y0);
int func (double y, const double v[], double f[],
	void *params);
int jac (double y, const double v[], double *dfdy, 
	double dfdt[], void *params);
void run_ode_solver(double y0, double h,
	vector<double> * ypts, vector<double> * solution);
double get_yH( vector<double> * ypts,
				vector<double> * solution,
				double * f_prime_at_yH_estimate );
double get_zH(vector<double> * ypts,
				vector<double> * solution,
				double f_prime_at_yH_estimate,
				double * f_dot_at_zH_estimate);
double integrate(vector<double> * xpts, vector<double> * fpts);
void get_zpts(vector<double> * ypts, vector<double> * solution, vector<double> * zpts);

//main driver
int main (void)
{
	const int nypts = 101;
	vector<double> ypts, solution;
	
	//initialize
	if (mode == 0)
	{
		for (int iy = 0; iy < nypts; ++iy)
		{
			//ypts.push_back(0.0 + (1.471981825e-05)*stepsize*double(iy));
			ypts.push_back(0.0 + 0.01*stepsize*double(iy));
			for (int iv = 0; iv < dim; ++iv)
				solution.push_back(0.0);
		}
	}

	//solve
	//double y0 = 0.0;
	double y0 = 1.e-10;
	double h = 1.e-10;
	run_ode_solver(y0, h, &ypts, &solution);

	//convert ypts to zpts too
	vector<double> zpts;
	get_zpts(&ypts, &solution, &zpts);

	//output results
	for (int iy = 0; iy < ypts.size(); ++iy)
	{
		cout << ypts[iy] << "   " << zpts[iy];	//ypts and zpts should be same length
		for (int iv = 0; iv < dim; ++iv)
			cout << "   " << solution[solution_indexer(iy, iv)] + double(iv==2)*q;
		cout << "   "
				<< Vtilde( solution[solution_indexer(iy, 0)],
							solution[solution_indexer(iy, 1)],
							ypts[iy] ) - 12.0*k*k
				<< "   "
				<< dVtildedphi( solution[solution_indexer(iy, 0)],
								solution[solution_indexer(iy, 1)],
								ypts[iy] )
				<< "   "
				<< dVtildedG( solution[solution_indexer(iy, 0)],
								solution[solution_indexer(iy, 1)],
								ypts[iy] );
		cout << endl;
	}

	bool verbose = true;
	if (mode == 1 and verbose)
	{
		double f_prime_at_yH = 0.0;
		double yH = get_yH(&ypts, &solution, &f_prime_at_yH);

		cout << "q = " << q << ", k = " << k << endl;
		cout << "y_H = " << yH << endl;
		//cout << "z_H = " << sgn << endl;
		cout << "f'(y_H) = " << f_prime_at_yH << endl;

		double f_dot_at_zH = 0.0;
		double zH = get_zH(&ypts, &solution, f_prime_at_yH, &f_dot_at_zH);
	
		cout << "z_H = " << zH << endl;
		cout << "dot f(z_H) = " << f_dot_at_zH
				<< " ==> T = " << -f_dot_at_zH / (4.0*M_PI) << " MeV" << endl;
	}

	return 0;
}

///////////////////////////////////////////////////////////



void set_initial_conditions(double v[], double y0)
{
	//INITIAL CONDITIONS TO GET 1D CASE
	/*
	//variable order: 0-phi, 1-G, 2-B, 3-f, 4-psi, 5-gamma
	v[0] = 0.0;							//phi(y=0) == 0.0
	v[1] = 0.0;							//G(y=0) == 0.0
	v[2] = W(v[0], v[1], 0.0);
	v[3] = 1.0;							//f(y=0) == 1.0, by definition
	v[4] = 6.0 * dWdphi(v[0], v[1], 0.0);
	v[5] = 6.0 * dWdG(v[0], v[1], 0.0);
	*/
	/*
	//INITIAL CONDITIONS AT Y0==0
	//variable order: 0-phi, 1-G, 2-B, 3-f, 4-psi, 5-gamma
	v[0] = sqrt(8.0/3.0)*lambda/(k*k);
										//phi(y=0) == sqrt(8./3.)*lambda/k^2
	v[1] = sqrt(24.0*lambda)/k;
										//G(y=0) == sqrt(24.*lambda)/k
	v[2] = 0.0;							//B = W
	v[3] = 1.0;							//f(y=0) == 1.0, by definition
	v[4] = exp(a*v[0])*4.0*sqrt(2.0/3.0)*lambda/k;
	v[5] = exp(a*v[0])*sqrt(24.0)*sqrt_lambda_parameter;
	*/
	/*
	//INITIAL CONDITIONS AT Y0==EPS
	//variable order: 0-phi, 1-G, 2-B, 3-f, 4-psi, 5-gamma
	double phi0 = sqrt(8.0/3.0)*lambda/(k*k);
	double phiPrime0 = exp(a*phi0)*4.0*sqrt(2.0/3.0)*lambda/k;
	double zeps = y0 * exp(a*phi0);	//to leading order in y0
	v[0] = sqrt(8.0/3.0)*lambda*(zeps + 1.0/k)*(zeps + 1.0/k);
										//phi(y=0) == sqrt(8./3.)*lambda/k^2
	v[1] = sqrt(24.0*lambda)*(zeps + 1.0/k);
										//G(y=0) == sqrt(24.*lambda)/k
	//v[2] = W(v[0], v[1], y0);			//B = W
	v[3] = 1.0;							//f(y=0) == 1.0, by definition
	v[4] = exp(a*v[0])*(1.0-a*phiPrime0*zeps)*4.0*sqrt(2.0/3.0)*lambda*(zeps + 1.0/k);
	v[5] = exp(a*v[0])*(1.0-a*phiPrime0*zeps)*sqrt(24.0)*sqrt_lambda_parameter;
	v[2] = y0 * ( v[4]*v[4]+v[5]*v[5] ) / 6.0;
	*/
	///*
	//CORRECTED INITIAL CONDITIONS AT Y0==EPS
	//variable order: 0-phi, 1-G, 2-B, 3-f, 4-psi, 5-gamma
	double phi0 = sqrt(8.0/3.0)*lambda/(k*k);
	double phiDot0 = 4.0*sqrt(2.0/3.0)*lambda/k;
	double phiDotDot0 = 4.0*sqrt(2.0/3.0)*lambda;
	double G0 = sqrt(24.0*lambda)/k;
	double GDot0 = sqrt(24.0*lambda);
	double GDotDot0 = 0.0;
	double expaphi0 = exp(a*phi0);

	v[0] = phi0 + y0*expaphi0*phiDot0;			//phi
	v[1] = G0 + y0*expaphi0*GDot0;				//G
	v[2] = W(v[0], v[1], y0);					//B
	v[3] = 1.0;									//f(y=0) == 1.0, by definition
	v[4] = expaphi0*phiDot0*( 1.0 + a*y0*phiDot0 )
			+ expaphi0*expaphi0*y0*( k*phiDot0 + phiDotDot0 );
	v[5] = expaphi0*GDot0*( 1.0 + a*y0*phiDot0 )
			+ expaphi0*expaphi0*y0*( k*GDot0 + GDotDot0 );
	//v[2] = y0 * ( v[4]*v[4]+v[5]*v[5] ) / 6.0;
	//*/

	//for (int i = 0; i < dim; ++i)
	//	cout << "v[" << i << "] = " << v[i] << endl;
//if (1) exit(8);

	return;	
}





int func (double y, const double v[], double f[],
   void *params)
{
	double Bpq = v[2]+q*sgn(y);
	double Vt = Vtilde(v[0], v[1], y);
	double dVtdphi = dVtildedphi(v[0], v[1], y);
	double dVtdG = dVtildedG(v[0], v[1], y);
//cout << "Fields: " << Bpq << "   " << Vt << "   " << dVtdphi << "   " << dVtdG << endl;

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
	double Bpq = v[2]+q*sgn(y);
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
	gsl_matrix_set (m, 0, 4, 1.0);		//checked
	//
	gsl_matrix_set (m, 1, 5, 1.0);		//checked
	//
	gsl_matrix_set (m, 2, 4, v[4]/3.0);	//checked
	gsl_matrix_set (m, 2, 5, v[5]/3.0);	//checked
	//
	gsl_matrix_set (m, 3, 2, ( 24.0*k*k-2.0*Vt+v[3]*( 24.0*Bpq*Bpq+v[4]*v[4]+v[5]*v[5] ) )
								/ (6.0*Bpq*Bpq));//checked
	gsl_matrix_set (m, 3, 3, 4.0*Bpq-(v[4]*v[4]+v[5]*v[5])/(6.0*Bpq));//checked
	gsl_matrix_set (m, 3, 4, -3.0*v[3]*v[4]/(3.0*Bpq));//checked
	gsl_matrix_set (m, 3, 5, -3.0*v[3]*v[5]/(3.0*Bpq));//checked
	//
	//gsl_matrix_set (m, 4, 2, -( 24.0*k*k-2.0*Vt+v[3]*( v[4]*v[4]+v[5]*v[5] ) )
	//							/ (6.0*v[3]*Bpq*Bpq));
	gsl_matrix_set (m, 4, 2, -v[4]*( 24.0*k*k-2.0*Vt+v[3]*( v[4]*v[4]+v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq*Bpq));
	gsl_matrix_set (m, 4, 3, ( v[4]*( Vt - 12.0*k*k ) - 3.0*dVtdphi*Bpq )
								/ (3.0*v[3]*v[3]*Bpq));//checked
	gsl_matrix_set (m, 4, 4, ( 24.0*k*k-2.0*Vt+v[3]*( 3.0*v[4]*v[4]+v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq));//checked
	gsl_matrix_set (m, 4, 5, v[4]*v[5] / (3.0*Bpq));//checked
	//
	/*
	gsl_matrix_set (m, 5, 2, -( 24.0*k*k-2.0*Vt+v[3]*( v[4]*v[4]+v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq*Bpq));
	gsl_matrix_set (m, 5, 3, ( v[5]*( Vt - 12.0*k*k ) - 3.0*dVtdG*Bpq )
								/ (3.0*v[3]*v[3]*Bpq));
	gsl_matrix_set (m, 5, 4, ( 24.0*k*k-2.0*Vt+v[3]*( v[4]*v[4]+3.0*v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq));
	gsl_matrix_set (m, 5, 5, v[4]*v[5] / (3.0*Bpq));
	*/
	gsl_matrix_set (m, 5, 2, -v[5]*( 24.0*k*k-2.0*Vt+v[3]*( v[4]*v[4]+v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq*Bpq));
	gsl_matrix_set (m, 5, 3, ( v[5]*( Vt - 12.0*k*k ) - 3.0*dVtdG*Bpq )
								/ (3.0*v[3]*v[3]*Bpq));
	gsl_matrix_set (m, 5, 4, v[4]*v[5] / (3.0*Bpq));
	gsl_matrix_set (m, 5, 5, ( 24.0*k*k-2.0*Vt+v[3]*( v[4]*v[4]+3.0*v[5]*v[5] ) )
								/ (6.0*v[3]*Bpq));

	for (int i = 0; i < dim; ++i)
		dfdt[i] = 0.0;

	return GSL_SUCCESS;
}


void run_ode_solver(double y0, double h,
	vector<double> * ypts, vector<double> * solution)
{
	const gsl_odeiv_step_type * T 
	 = gsl_odeiv_step_rk8pd;
	//const gsl_odeiv_step_type * T 
	// = gsl_odeiv_step_bsimp;

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
			set_initial_conditions(v, y0);

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
		//long int countmax = 1000;
		do
		{
			double y = y0, y1 = y0 + stepsize*double(count);

			//set initials conditions
			set_initial_conditions(v, y0);
			//printf ("ICs: %.5e %.5e %.5e %.5e %.5e %.5e %.5e\n",
			//		y, v[0], v[1], v[2], v[3], v[4], v[5]);

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
		}	//termination condition is f(y>yH) < 0.0
		while( v[3] >= 0.0 /*and count <= countmax*/ );
		/*if (count >= countmax)
		{
			cerr << "Never terminated!  Exiting..." << endl;
			exit(8);
		}*/
	}

	gsl_odeiv_evolve_free (e);
	gsl_odeiv_control_free (c);
	gsl_odeiv_step_free (s);

	return;
}


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

	//cout << "CHECK: " << (A0 - 2.0*A1 + A2) / (2.0*dx*dx) << endl;

	return ( yH_estimate );
}
double get_zH(vector<double> * ypts,
				vector<double> * solution,
				double f_prime_at_yH_estimate,
				double * f_dot_at_zH_estimate)
{
	int length = ypts->size();
	vector<double> integrand(length);
	double phi_at_H_estimate = solution->at(solution_indexer(length-2, 0));
	for (int iy = 0; iy < length; ++iy)
	{
		double phi = solution->at(solution_indexer(iy, 0));
		integrand.at(iy) = exp(a*phi);
	}

	double arg = integrate(ypts, &integrand);
	double zH_estimate = (exp(k*arg) - 1.0)/k;

	//finally, estimate f'(z_H) < 0
	(*f_dot_at_zH_estimate)
		= exp(-a*phi_at_H_estimate)
			* f_prime_at_yH_estimate
			/ ( k*zH_estimate + 1.0 );

	return ( zH_estimate );
}
//naive trapezoid rule
double integrate(vector<double> * xpts, vector<double> * fpts)
{
	int length = fpts->size();
	double h = xpts->at(1) - xpts->at(0);	//assume uniform grid

	double result = 0.0;
	for (int ix = 0; ix < length - 1; ++ix)	//last value is probably NaN, so cut it out
		result += fpts->at(ix);
	
	result = 0.5*h*(2.0*result - fpts->at(0) - fpts->at(length-2));

	return ( result );
}


void get_zpts(vector<double> * ypts, vector<double> * solution, vector<double> * zpts)
{
	//int length = ypts->size()-2;	//last two steps generally contain NaNs
	int length = ypts->size();
	zpts->clear();
	
	vector<double> yrange, integrand;
	for (int iz = 1; iz < length; ++iz)
	{
		yrange.clear();
		integrand.clear();
		for (int iy = 0; iy <= iz; ++iy)
		{
			double phi = solution->at(solution_indexer(iy, 0));
			yrange.push_back(ypts->at(iy));
			integrand.push_back(exp(a*phi));
//if (iz==4) cout << yrange[iy] << "   " << integrand[iy] << endl;
		}

		double arg = integrate(&yrange, &integrand);
//if (iz==4) cout << "CHECK: " << q << "   " << arg << endl;
		zpts->push_back((exp(k*arg) - 1.0)/k);
	}

	return;
}





