#ifndef POTENTIALS_H
#define POTENTIALS_H

#include <cmath>
#include <cstdlib>

#include "AdSYMatFiniteT.h"

#define U_MODE				0		//which U(G) potential to use
#define INCLUDE_SCALARS		true	//whether or not to include phi, G, and B

using namespace std;

const double a = 1.0/sqrt(6.0);

inline double Theta(double y)
{
	return double( y>=0.0 );
}

inline double sgn(double y)
{
	return ( Theta(y) - Theta(-y) );
}

double U(double G);
double dUdG(double G);
double d2UdG2(double G);

double W(double phi, double G, double y);
double dWdphi(double phi, double G, double y);
double dWdG(double phi, double G, double y);
double d2Wdphi2(double phi, double G, double y);
double d2WdphidG(double phi, double G, double y);
double d2WdG2(double phi, double G, double y);

double Vtilde(double phi, double G, double y);
double dVtildedphi(double phi, double G, double y);
double dVtildedG(double phi, double G, double y);

#endif
