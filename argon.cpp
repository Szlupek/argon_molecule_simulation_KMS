#include "argon.h"
int main()
{
	srand( time( NULL ) );
	double a = 0.38;
	double k= 0.00831;
	int n = 3;
	int N = n*n*n;
	double m =40;//39.948;
	double b0[3]={a,0.0,0.0};
	double b1[3]={a/2.0,a*sqrt(3)/2,0};
	double b2[3]={a/2,a*sqrt(3)/6,a*sqrt(2.0/3.0)};
	double T0 =0;
	double R = 0.38;
	double e =1;
	double t = 0.001;
	double L  = 2.3;//2.3;//1.2
	double f =10000.0;
	const double dt = 0.01;
	const unsigned int Steps = 100;
	const unsigned int Steps_0 = 1;
	const unsigned int Steps_xyz = 10;
	const unsigned int Steps_out = 10;
	atom atoms[N];
	//ifstream params;
	//params.open( "params.txt" );
	//params>>n;



	setR0(n,atoms,b0,b1,b2);
	//saveXYZ("test1.txt", N,atoms);
	setE0(N,atoms,k,T0);
	setP0(N,atoms,m);
	vanDerWaalsPotArray(e,R, atoms, N);
	SPotArray(L,f, atoms,N);
	sumV(N,atoms);
	double sumV=0;

	fstream dataOut;
	dataOut.open("data.txt", ios::out | ios::trunc);

	fstream dane;
	dane.open("dane.xyz", ios::out | ios::trunc);
    if(dane.good() == true)
    {
		//dane<<N<<endl;
		//dane<<"Argon"<<endl;
		for(int ii=0;ii<N;ii++)
		dane/*<<"Ar "*/<<atoms[ii].r[0]<<" "<< atoms[ii].r[1]<<" "<<atoms[ii].r[2]<<endl;
		dane<<endl;
	}

	for(int ii=0;ii<N;ii++)sumV+=atoms[ii].V;
	sumV=sumV/2;
	cout<<sumV<<endl;

	//for(int jj=0;jj<N;jj++){cout<<atoms[jj].p[0]<<" "<<atoms[jj].p[1]<<" "<<atoms[jj].p[2]<<endl;}

	dataOut<<"t= "<<0<<" T= "<<getT(N,atoms,k,m)<<endl;

	for(int ii =0;ii<Steps+Steps_0;ii++)
	{
		vanDerWaalsForceArray(e,R, atoms, N);
		SForceArray(atoms, L,f,N);
		sumF(N,atoms);
		evolv(dt, atoms,N, m);

		if(!(ii<Steps_0))
		{
			if((ii-Steps_0)%Steps_xyz==0)
			{
				for(int jj=0;jj<N;jj++)
				dane/*<<"Ar "*/<<atoms[jj].r[0]<<" "<< atoms[jj].r[1]<<" "<<atoms[jj].r[2]<<endl;
				dane<<endl<<endl;
			}
			if((ii-Steps_0)%Steps_out==0)
			{
				dataOut<<"t= "<<ii*dt<<" T= "<<getT(N,atoms,k,m)<<endl;
			}


		}


	}
	dane.close();
	dataOut.close();
	return 0;
}
//Vs = - 600 dla N =5
