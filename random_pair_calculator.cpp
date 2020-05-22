//Computes RR


#include <iostream>
#include <cmath>
#include <valarray>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include "/storage/data03/Users/rpandit/logistics/headers/randparams.h"
#include "/storage/data03/Users/rpandit/logistics/headers/nrutil.h"
#include "/storage/data03/Users/rpandit/logistics/headers/nrutil.cpp"

using namespace std;

int main (){



	cout << "Begin count RR ranging from "<< smin <<" to "<< smax <<" Mpc in "<< sbins <<"x"<< mubins <<" bins" << endl; 

	int a1, a2, a3, a4, a5, a6, aux; 	
	double half_length = L/2.;
	double  rx, ry, rz, rperp, s, mu; 	
	
	

	double *xr=dvector(0,m);	
	double *yr=dvector(0,m);
	double *zr=dvector(0,m);
       	

	double  RR[sbins+1]={0};	
	double  RR2[sbins+1][mubins+1]={0};
	double 	RR3[sbins+1][sbins+1]={0};


/*Open random*/	

	ifstream in ("/storage/data03/Users/rpandit/logistics/randoms/random_x_y_z_L1500_1M.dat");
	assert(in);
	for (int i=0; i<=m; i++){ xr[i]=0;yr[i]=0;zr[i]=0;}
	for (int i=1; i<=m; i++)in >> aux >> xr[i]>> yr[i] >> zr[i];
	in.close();

//RR
 	cout<<"Computing RR"<<endl;

	#pragma omp parallel for private(a1, a2, a3, a4, rx, ry, rz, s, mu, rperp)
	
	for(int i=1; i<=m; i++){	
		for(int j=i+1; j<=m; j++){
			double rx_p=abs(xr[i]-xr[j]), ry_p=abs(yr[i]-yr[j]), rz_p=abs(zr[i]-zr[j]); 
			       
			rx= ( rx_p > half_length ?  L-rx_p : rx_p);
			ry= ( ry_p > half_length ?  L-ry_p : ry_p);
			rz= ( rz_p > half_length ?  L-rz_p : rz_p);
				
			s=sqrt(pow(rx, 2)+pow(ry, 2)+pow(rz, 2));
                               
			mu=rx/s;
 
			rperp=sqrt(pow(ry,2)+pow(rz,2));			
				
	
			if ((rx<=smax && rx>=smin) && (rperp<=smax && rperp>=smin)){			
				a1=(int)floor((rperp-smin)/delta_s);
			 	a2=(int)floor((rx-smin)/delta_s);
				#pragma omp atomic
				RR3[a1+1][a2+1]++;
           		 }     

			//2D correlation function
			if (s<=smax && s>=smin){
              		 	a3=(int)floor((s-smin)/delta_s);    
              			a4=(int)floor((mu-mu_min)/delta_mu);
				#pragma omp atomic
 				RR[a3+1]++;
				#pragma omp atomic
		        	RR2[a3+1][a4+1]++;
		        }
	
			
		}
	
	}


	ofstream rr("/storage/data03/Users/rpandit/logistics/randoms/random_pairs_12x50_1M_200Mpc_1d.dat");
		for (int i=1; i<=sbins; i++)
			rr << RR[i] << endl;
		
		rr.close();

	ofstream rr2("/storage/data03/Users/rpandit/logistics/randoms/random_pairs_12x50_1M_200Mpc_s_mu.dat");
		for (int i=1; i<=sbins; i++)
		for (int j=1; j<=mubins; j++)
			rr2 << RR2[i][j] << endl;

		rr2.close();

	ofstream rr3("/storage/data03/Users/rpandit/logistics/randoms/random_pairs_12x50_1M_200Mpc_spar_sperp.dat");
		for (int i=1; i<=sbins; i++)
		for (int j=1; j<=sbins; j++)
			rr3 << RR3[i][j] << endl;

		rr3.close();


	return 0;

}
