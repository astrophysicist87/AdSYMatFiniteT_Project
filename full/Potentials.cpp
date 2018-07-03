#include <iostream>
#include <cmath>

#include "Potentials.h"

using namespace std;

//potentials for initializing at y==0

double U(double G)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	double result = 0.0;
	switch(U_MODE)
	{
		case 0:
			result = G*G/12.0;
			break;
		case 1:
			result = G*G/12.0 + log(1.0 + gamma_parameter*G*G) / (4.0*gamma_parameter);
			break;
		default:
			cout << "U_MODE = " << U_MODE << " not supported!" << endl;
			exit(1);
			break;
	}
	return ( result );
}

double dUdG(double G)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	double result = 0.0;
	switch(U_MODE)
	{
		case 0:
			result = G/6.0;
			break;
		case 1:
			result = G*(4.0 + gamma_parameter*G*G)/( 6.0*(1.0 + gamma_parameter*G*G) );
			break;
		default:
			cout << "U_MODE = " << U_MODE << " not supported!" << endl;
			exit(1);
			break;
	}
	return ( result );
}

double d2UdG2(double G)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	double result = 0.0;
	double gG2 = gamma_parameter*G*G;
	switch(U_MODE)
	{
		case 0:
			result = 1.0/6.0;
			break;
		case 1:
			result = ( 4.0 + gG2*(gG2-1.0) ) / ( 6.0 * (1.0+gG2) * (1.0+gG2) );
			break;
		default:
			cout << "U_MODE = " << U_MODE << " not supported!" << endl;
			exit(1);
			break;
	}
	return ( result );
}

double W(double phi, double G, double y)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	return ( k * sgn(y) * ( exp(a*phi) * ( U(G) + 1.0 - a*phi ) - 1.0 ) );
}

double dWdphi(double phi, double G, double y)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	return ( k * sgn(y) * a * exp(a*phi) * ( U(G) - a*phi ) );
}

double dWdG(double phi, double G, double y)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	return ( k * sgn(y) * exp(a*phi) * dUdG(G) );
}

double d2Wdphi2(double phi, double G, double y)
{
	return ( k * sgn(y) * a * a * exp(a*phi) * ( U(G) - a*phi - 1.0 ) );
}

double d2WdphidG(double phi, double G, double y)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	return ( k * sgn(y) * a * exp(a*phi) * dUdG(G) );
}

double d2WdG2(double phi, double G, double y)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	return ( k * sgn(y) * exp(a*phi) * d2UdG2(G) );
}

double Vtilde(double phi, double G, double y)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	double W_local = W(phi, G, y);
	double dWdphi_local = dWdphi(phi, G, y);
	double dWdG_local = dWdG(phi, G, y);

	//cout << "Vt(" << phi << ", " << G << ", " << y << ") = " <<
	//	18.0 * ( dWdphi_local*dWdphi_local + dWdG_local*dWdG_local )
	//			-12.0 * W_local*W_local
	//			-24.0*k*W_local*( Theta(y) - Theta(-y) ) << endl;

	return ( 18.0 * ( dWdphi_local*dWdphi_local + dWdG_local*dWdG_local )
				-12.0 * W_local*W_local
				-24.0*k*W_local*( Theta(y) - Theta(-y) ) );
}

double dVtildedphi(double phi, double G, double y)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	double W_local = W(phi, G, y);
	double dWdphi_local = dWdphi(phi, G, y);
	double dWdG_local = dWdG(phi, G, y);
	double d2Wdphi2_local = d2Wdphi2(phi, G, y);
	double d2WdphidG_local = d2WdphidG(phi, G, y);
	//double d2WdG2_local = d2WdG2(phi, G, y);

	//cout << "dVtildedphi(" << phi << ", " << G << ", " << y << ") = " <<
	//	36.0 * ( dWdphi_local*d2Wdphi2_local + dWdG_local*d2WdphidG_local )
	//			-24.0 * W_local*dWdphi_local
	//			-24.0*k*dWdphi_local*sgn(y) << endl;

	return ( 36.0 * ( dWdphi_local*d2Wdphi2_local + dWdG_local*d2WdphidG_local )
				-24.0 * W_local*dWdphi_local
				-24.0*k*dWdphi_local*sgn(y) );
}

double dVtildedG(double phi, double G, double y)
{
	if (not INCLUDE_SCALARS)
		return ( 0.0 );

	double W_local = W(phi, G, y);
	double dWdphi_local = dWdphi(phi, G, y);
	double dWdG_local = dWdG(phi, G, y);
	//double d2Wdphi2_local = d2Wdphi2(phi, G, y);
	double d2WdphidG_local = d2WdphidG(phi, G, y);
	double d2WdG2_local = d2WdG2(phi, G, y);

	//cout << "dVtildedG(" << phi << ", " << G << ", " << y << ") = " <<
	//	36.0 * ( dWdphi_local*d2WdphidG_local + dWdG_local*d2WdG2_local )
	//			-24.0 * W_local*dWdG_local
	//			-24.0*k*dWdG_local*sgn(y) << endl;

	return ( 36.0 * ( dWdphi_local*d2WdphidG_local + dWdG_local*d2WdG2_local )
				-24.0 * W_local*dWdG_local
				-24.0*k*dWdG_local*sgn(y) );
}


