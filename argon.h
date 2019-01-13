#ifndef _argon_h_
#define _argon_h_
#define _USE_MATH_DEFINES
#include <cstring>
#include <iostream>
#include <cstdlib>
#include <math.h>
#include <fstream>
#include <cstdlib>
using namespace std;

class atom
{
	public:
	int index;
	double r[3];
	double E[3];
	double p[3];
	double F[3];
	double Fv[3];
	double Vv;
	double Vs;
	double V;
	double Fs[3];
	atom();
};

double randomNum(void);
int randomSign(void);
double dyst(double, double, double,double, double, double);
double dyst2(atom, atom);
double lenght(atom);
double SPotential(atom,double,double);
void SPotArray(double, double, atom[],int);
double vanDerWaalsPot(double, double, atom, atom);
void vanDerWaalsPotArray(double, double, atom[],int);
double vanDerWaalsForce(double, double, atom, atom, int);
void vanDerWaalsForceArray(double, double, atom[],int);
void SForce(atom ,double,double);
void SForceArray(atom[],double ,double ,int);
double getP(double,atom*, int);
void evolv(double, atom[], int,double);
void setR0(int, atom[], double[],double[],double[]);
void saveXYZ(string, int,atom[]);
void setE0(int,atom[], double,double);
void setP0(int,atom[],double);
void sumV(int,atom[]);
void sumF(int,atom[]);
double getEk(double, atom);
double getT(int, atom[], double,double);
#endif
