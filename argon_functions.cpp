#include "argon.h"

atom::atom(){}

double randomNum()
{
		double randomNumber =rand()% 100000;
	return randomNumber / 100000;
}

int randomSign()
{
	int tmp = rand()%2;
	if(tmp ==1 )
	{
		return 1;
	}
	else
	{
		return -1;
	}
}

double lenght(atom a1)
{
	return sqrt( pow(a1.r[0],2)+pow(a1.r[1],2)+pow(a1.r[2],2));
}

double dyst(double r1x, double r1y, double r1z,double r2x, double r2y, double r2z)
{
	return sqrt( pow(r1x-r2x,2)+pow(r1y-r2y,2)+pow(r1z-r2z,2) );

}

double dyst2(atom a1, atom a2)
{
	return sqrt(pow(a1.r[0]-a2.r[0],2)+pow(a1.r[1]-a2.r[1],2)+pow(a1.r[2]-a2.r[2],2));
}

double SPotential(atom a,double L,double f)
{
	double r = lenght(a);
	if(r>L)
	{
		return f*pow(r-L,2)/2;
	}
	else
	{
		return 0;
	}
}
void SPotArray(double L, double f, atom a[],int N)
{
	for(int ii=0;ii<N;ii++)
	{
		a[ii].Vs=SPotential(a[ii],L,f);

	}
}
void SForce(atom a,double L,double f)
{
	double r = lenght(a);
	if(r>L)
	{
		for(int ii=0;ii<3;ii++) a.Fs[ii]=f*(r-L)*a.r[ii]/r;
	}
	else {}
}
void SForceArray(atom a[],double L,double f,int N)
{
	for(int ii=0;ii<N;ii++)
	{
		SForce(a[ii],L,f);
	}
}

double vanDerWaalsPot(double e, double R, atom a1, atom a2)
{
	double r= dyst2(a1, a2);
	return e*(pow(R/r,12)-2*pow(R/r,6));
}

double vanDerWaalsForce(double e, double R, atom a1, atom a2, int dir)
{
	double r= dyst2(a1, a2);
	return 12*e*(pow(R/r,12)-2*pow(R/r,6))/r/r*(a1.r[dir]-a2.r[dir]);
}

double getP(double L,atom *a[], int N)
{
	double tmp =0;
	for (int ii=0; ii<N; ii++){tmp+=sqrt(pow(a[ii]->F[0],2)+pow(a[ii]->F[1],2)+pow(a[ii]->F[2],2));}
	return tmp/4/M_PI/L/L;

}
 void vanDerWaalsPotArray(double e, double R, atom a[],int N)
{
	for(int ii=0;ii<N;ii++)
	{
		a[ii].Vv=0;
		for(int jj=0;jj<N;jj++)
		{
			if(!(ii==jj))
			{
				a[ii].Vv+=vanDerWaalsPot(e,R,a[ii],a[jj]);
			}
		}
	}
}
 void vanDerWaalsForceArray(double e, double R, atom a[],int N)
{
	for(int ii=0;ii<N;ii++)
	{
		a[ii].Vv=0;
		for(int jj=0;jj<N;jj++)
		{
			if(!(ii==jj))
			{
				for(int kk=0;kk<3;kk++) a[ii].Fv[kk]+=vanDerWaalsForce(e,R,a[ii],a[jj],kk);
			}
		}
	}
}
void evolv(double t, atom a[], int N,double m)
{
	for(int ii=0;ii<N;ii++)
	{
		for(int jj=0;jj<3;jj++)
		{
			a[ii].p[jj]=a[ii].p[jj]+t*a[ii].F[jj]/2;
		//	cout<<"p= "<<a[ii].p[jj]<<"+"<<t<<"*"<<a[ii].F[jj]<<"/"<<2<<" = "<<a[ii].p[jj]<<endl;
			a[ii].r[jj]=a[ii].r[jj]+t*a[ii].p[jj]/m;
		//	cout<<"r= "<<a[ii].r[jj]<<"+"<<t<<"*"<<a[ii].p[jj]<<"/"<<m<<" = "<<a[ii].r[jj]<<endl;
			a[ii].p[jj]=a[ii].p[jj]+t*a[ii].F[jj]/2;
		//	cout<<"p2= "<<a[ii].p[jj]<<"+"<<t<<"*"<<a[ii].F[jj]<<"/"<<2<<" = "<<a[ii].p[jj]<<endl;
		}
	}
}

void setR0(int n, atom atoms[], double b0[],double b1[],double b2[])
{
	for(int i0=0;i0<n;i0++)
	{
		for(int i1=0;i1<n;i1++)
		{
			for(int i2=0;i2<n;i2++)
			{
				int i = i0*n*n+i1*n+i2;
				for(int ii=0;ii<3;ii++)
				{
					atoms[i].r[ii]= (i0-(n-1)/2)*b0[ii]+(i1-(n-1)/2)*b1[ii]+(i2-(n-1)/2)*b2[ii];
				}

			}
		}
	}

}

void saveXYZ(string name, int N,atom atoms[])
{
	fstream dane;
	dane.open(name.c_str(), ios::out | ios::trunc);
    if(dane.good() == true)
    {
		dane<<N<<endl;
		dane<<"Argon"<<endl;
		for(int ii=0;ii<N;ii++)
		dane<<"Ar "<<atoms[ii].r[0]<<" "<< atoms[ii].r[1]<<" "<<atoms[ii].r[2]<<endl;
	}
	dane.close();
}

void setE0(int N,atom atoms[], double k,double T0)
{
	double sumE[3]={0,0,0};
	for(int ii=0;ii<N;ii++)
	{
		for(int jj=0;jj<3;jj++)
		{
			atoms[ii].E[jj]=-T0/2.0 * k*log(randomNum());
			sumE[jj]+=atoms[ii].E[jj];
		}

	}

	for(int jj=0;jj<3;jj++)sumE[jj]=sumE[jj]/N;
	double mean = k* T0/2.0;

	double normE[3];
	if(T0!=0)
	{
		for(int jj=0;jj<3;jj++)normE[jj]=mean/sumE[jj];

		for(int ii=0;ii<N;ii++)
		{
			for(int jj=0;jj<3;jj++)
			{
				atoms[ii].E[jj]=atoms[ii].E[jj]*normE[jj];
			}
		}
	}
}
void setP0(int N,atom atoms[],double m)
{
	double sumP[3]={0,0,0};

	for(int ii=0;ii<N;ii++)
	{
		for(int jj=0;jj<3;jj++)
		{
			int tmp =randomSign();
			atoms[ii].p[jj]= tmp*sqrt(2*m*atoms[ii].E[jj]);
			//cout<<"p= "<<tmp<<"*"<<atoms[ii].E[jj]<<endl;
			sumP[jj]+=atoms[ii].p[jj];
		}
	}
	for(int jj=0;jj<3;jj++){sumP[jj]=sumP[jj]/N;}

	for(int ii=0;ii<N;ii++)
	{
		for(int jj=0;jj<3;jj++){atoms[ii].p[jj]=atoms[ii].p[jj]-sumP[jj];}
	}
}
void sumV(int N,atom atoms[])
{
	for(int ii=0;ii<N;ii++){atoms[ii].V=atoms[ii].Vv+atoms[ii].Vs;}
}
void sumF(int N,atom atoms[])
{
	for(int ii=0;ii<N;ii++){for(int jj=0;jj<3;jj++)atoms[ii].F[jj]=atoms[ii].Fv[jj]+atoms[ii].Fs[jj];}
}
double getEk(double m, atom a)
{
	double absP=sqrt(pow(a.p[0],2)+pow(a.p[1],2)+pow(a.p[2],2));
	return pow(absP,2)/2/m;
}
double getT(int N, atom atoms[], double k,double m)
{
	double sumE =0;
	for(int ii= 0;ii<N;ii++)sumE+=getEk(m,atoms[ii]);
	return sumE*2.0/3.0/N/k;
}

