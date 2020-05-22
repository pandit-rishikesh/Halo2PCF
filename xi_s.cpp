// Calculates DD RR and DR binned in s(redshift distacne and mu(cos \theta) in redshift space. 


//====================================== Calling Header files ======================================//
#include <iostream>
#include <cmath>
#include <valarray>
#include <fstream>
#include <string>
#include <sstream>
#include <cassert>
#include <iomanip>
#include <string.h>
#include <stdlib.h>
#include <omp.h>
#include <time.h>
#include <sys/time.h>
#include "/storage/DATA-03/astrorm3/Users/rpandit/logistics/headers/timer.h"
#include "/storage/DATA-03/astrorm3/Users/rpandit/logistics/headers/params.h"
#include "/storage/DATA-03/astrorm3/Users/rpandit/logistics/headers/nrutil.cpp"
#include "/storage/DATA-03/astrorm3/Users/rpandit/logistics/headers/nrutil.h"

using namespace std;





int main (int argc, char *argv[]){


	int startTime, endTime, totalTime;

	startTime = time(NULL);

//==================================== argument for given range run ===============================//

	if (argc != 3)
	{
		cout <<"Please provide a range e.g. : main <from> <to>" <<endl;
		return 0;
	}

	int from = atoi(argv[1]);
	int to = atoi(argv[2]);

	
	cout << "Calculating pairs for catfiles :"<< from << " to " << to << endl;
	cout << "==================================================" << endl;	
	



	double *xr=dvector(0,m);	
	double *yr=dvector(0,m);
	double *zr=dvector(0,m);

//===================================== Looping over realisations =================================//
	for (int num=from; num<=to; num++){

//	cout <<"Hubble constant H(z) in [h Mpc/s/km^2] at z=0 = " << Hofz << endl;
	
	cout << "Current run is "<< num << endl;	

	

	string number;
	string redshiftno;
	stringstream convert;
	stringstream convert_z;
	convert << num;
	convert_z << Z;
	number = convert.str();
	redshiftno = convert_z.str();

	cout << "\n"  <<"We are currently at redshift z = "+redshiftno << endl;
	// Here I define and initiate variables and arrays
	
	int w,a1,a2,a3,a4,a5,a6,a7,a8; 	
	double mu;
	double half_length=L/2.;
	double sx, sy, sz, s_perp, s; 	
	double aux;

/*	double *DD=dvector(0,sbins), *RR=dvector(0,sbins), *DR=dvector(0,sbins); 
	double **DD2=dmatrix(0, sbins, 0, sbins), **RR2=dmatrix(0, sbins, 0, sbins), **DR2=dmatrix(0, sbins, 0, sbins); 
	
	double **DD3=dmatrix(0, sbins, 0, sbins), **RR3=dmatrix(0, sbins, 0, sbins), **DR3=dmatrix(0, sbins, 0, sbins); 
*/
	double DD[sbins+1], RR[sbins+1], DR[sbins+1]; 
	double DD2[sbins+1][sbins+1], RR2[sbins+1][sbins+1], DR2[sbins+1][sbins+1]; 
	
	double DD3[sbins+1][sbins+1], RR3[sbins+1][sbins+1], DR3[sbins+1][sbins+1]; 

	for(int i=0; i<=sbins; i++){DD[i]=0; RR[i]=0; DR[i]=0;}

	for(int i=0; i<=sbins; i++)
		for(int j=0; j<=sbins; j++) { DD2[i][j]=0; RR2[i][j]=0; DR2[i][j]=0; DD3[i][j]=0; RR3[i][j]=0; DR3[i][j]=0; }


//======================================= Reading DATA ============================================//
	
	long numLines = 0;
	
	ifstream in1 (string("/storage/DATA-03/astrorm3/Users/rpandit/sim/zHorizon/gc/anukruti/z_0/41-200/run"+number+".dat").c_str());
	assert(in1);

	string unused;
	while ( std::getline(in1, unused) )
	++numLines;
	
 	long n = numLines;

	cout << string( "For run"+number+":").c_str() << "\t" << n << " halos found" << endl;

 	double *xd=dvector(0,n);
 	double *yd=dvector(0,n);
 	double *zd=dvector(0,n);
	double *vx=dvector(0,n);
	
	in1.close();

	//cout << "initialized and ready for computing" << "\t" << t << endl;

	/*Open data*/

	
	ifstream in (string("/storage/DATA-03/astrorm3/Users/rpandit/sim/zHorizon/gc/anukruti/z_0/41-200/run"+number+".dat").c_str());
	assert(in);
	
	for (long i=0; i<=n; i++) {xd[i]=0;yd[i]=0;zd[i]=0;vx[i]=0;} 
	for (long i=1; i<=n; i++) in >> xd[i]>> yd[i] >> zd[i] >> vx[i] >> aux >> aux >> aux;	
	in.close();


    cout << "Pec. Velcoity before gadget norm + distance units :" << "\t" << vx[100] <<" km/s" << endl;
	cout << string("H(z) at redshift in [h Mpc/s/km^2]  z = "+redshiftno+" is :").c_str() << "\t" << Hofz << endl;


	for (long i=1; i<=n; i++) vx[i]/=(gadget_normalize*Hofz);       // measuring peculiar velocity in distance units (h^-1 Mpc).

	cout << "Pec. Velcoity AFTER gadget norm + distance units :"<< "\t" << vx[100] << " Mpc/h" << endl;


	/*Open random*/	

	ifstream in2 ("/storage/DATA-03/astrorm3/Users/rpandit/logistics/randoms/random_x_y_z_L1500_1M.dat");
	assert(in2);
	for (long i=0; i<=m; i++){ xr[i]=0;yr[i]=0;zr[i]=0;}
	for (long i=1; i<=m; i++)in2 >> w >> xr[i]>> yr[i] >> zr[i];


	cout << "Random catfile with "<< m << " particles " << endl;
	in2.close();


	cout << "All arrays initialized and data is acquired" << endl;


//==================================    Computational Zone   ======================================//
	
	cout << string("Computing pairs for run"+number+", I need some time").c_str() << endl;
				/*            DD             */

	cout<<"Going for DD.. this one looks easier" << endl;

	#pragma omp parallel for private(a1,a2,a3,a4,sx,sy,sz,s,mu,s_perp) //parallel looping OpenMP

	for(long i=1; i<=n; i++){
		for(long j=i+1; j<=n; j++){
			
				double sx_p=abs((xd[i]+vx[i])-(xd[j]+vx[j])), 
				       sy_p=abs(yd[i]-yd[j]), 
				       sz_p=abs(zd[i]-zd[j]);	
				
				sx= ( sx_p >half_length ?  L-sx_p : sx_p);	//spar
				sy= ( sy_p >half_length ?  L-sy_p : sy_p);
				sz= ( sz_p >half_length ?  L-sz_p : sz_p);
			
				s=sqrt(pow(sx,2)+pow(sy,2)+pow(sz,2));
			
				mu=sx/s;

				s_perp=sqrt(pow(sy,2)+pow(sz,2));  // sperp

			if ((sx<=smax && sx>=smin) && (s_perp<=smax && s_perp>=smin)){			
			 	a1=(int)floor((sx-smin)/delta_s);
				a2=(int)floor((s_perp-smin)/delta_s);
				#pragma omp atomic
				DD3[a1+1][a2+1]++;		// for xi_(spar, sperp)

			}

			if (s<=smax && s>=smin){
				a3=(int)floor((s-smin)/delta_s);	
				a4=(int)floor((mu-mu_min)/delta_mu);
				#pragma omp atomic
				DD[a3+1]++;			// for xi_(s)
				#pragma omp atomic
				DD2[a3+1][a4+1]++;		// for xi_(s,mu)
			}
				
				 			
				
                   	
		}
			

	}
	



	cout << "DD is done " << endl;

	                       /*            RR            */
	

	cout << "reading RR from database! (please check once if the binning is same as DD and DR)" << endl;
		
	ifstream in3 ("/storage/DATA-03/astrorm3/Users/rpandit/logistics/randoms/random_pairs_50x50_1M_200Mpc_1d.dat");
	assert(in3);
	for (int i=1; i<=sbins; i++) in3 >> RR[i] ;
	in3.close();	
	
		
	ifstream in4 ("/storage/DATA-03/astrorm3/Users/rpandit/logistics/randoms/random_pairs_50x50_1M_200Mpc_s_mu.dat");
	assert(in4);
	for (int i=1; i<=sbins; i++)
		for (int j=1; j<=mubins; j++)
	    	in4 >> RR2[i][j] ;
		
	in4.close();
	
	ifstream in5 ("/storage/DATA-03/astrorm3/Users/rpandit/logistics/randoms/random_pairs_50x50_1M_200Mpc_spar_sperp.dat");
	assert(in5);
	for (int i=1; i<=sbins; i++)
		for (int j=1; j<=sbins; j++)
	    	in5 >> RR3[i][j] ;
		
	in5.close();

	cout << "RR is done"<< endl;

				/*          DR          */
	

	cout <<"Now going for DR... will take much longer"<<endl;
	

	#pragma omp parallel for private(a5,a6,a7,a8,sx,sy,sz,s,mu,s_perp)	

	for(long i=1; i<=n; i++){
		for(long j=1; j<=m; j++){
			double sx_p=abs((xd[i]+vx[i])-xr[j]), 
			       sy_p=abs(yd[i]-yr[j]), 
			       sz_p=abs(zd[i]-zr[j]);	
				
			       sx= ( sx_p >half_length ?  L-sx_p : sx_p);        // spar
			       sy= ( sy_p >half_length ?  L-sy_p : sy_p);
			       sz= ( sz_p >half_length ?  L-sz_p : sz_p);
			
			       s=sqrt(pow(sx,2)+pow(sy,2)+pow(sz,2));
			       
			       mu=sx/s;
		
			       s_perp=sqrt(pow(sy,2)+pow(sz,2));		// sperp	
				

			if ((sx<=smax && sx>=smin) && (s_perp<=smax && s_perp>=smin)){			
			 	a5=(int)floor((sx-smin)/delta_s);
				a6=(int)floor((s_perp-smin)/delta_s);
				#pragma omp atomic
				DR3[a5+1][a6+1]++;		// for xi_(spar, sperp)
        		}     

								
			if (s<=smax && s>=smin){
				a7=(int)floor((s-smin)/delta_s);	
				a8=(int)floor((mu-mu_min)/delta_mu);
				#pragma omp atomic
				DR[a7+1]++;			// for xi_(s)
		     		#pragma omp atomic
		      	 	DR2[a7+1][a8+1]++;		// for xi_(s,mu)
		    }
				
		
		}


	}

	cout << "YUP! DR is done."<< endl;

//======================================= Writing output =========================================//

	cout << "Finally.. Writing output files" << endl;

	
	ofstream ddr(string("/storage/DATA-03/astrorm3/Users/rpandit/output/2pcf/red/z_0/pairs/41-200/ddrrdr_s_12bins_50_200Mpc_run"+number+".dat").c_str());

	for(int i=1; i<=sbins; i++)
		  
	  ddr << setiosflags(ios::uppercase) << std::scientific  << smin+(i-0.5)*delta_s << "\t" << DD[i] << "\t" << RR[i] << "\t" << DR[i] << endl;
	

	ddr.close();
	
	ofstream ddr2(string("/storage/DATA-03/astrorm3/Users/rpandit/output/2pcf/red/z_0/pairs/41-200/ddrrdr_s_mu_50x50_50_200Mpc_run"+number+".dat").c_str());
	for(int i=1; i<=sbins; i++)

		for(int j=1.; j<=mubins; j++) 
		  
		ddr2 << setiosflags(ios::uppercase) << std::scientific << smin+(i-0.5)*delta_s << "\t" << mu_min+(j-0.5)*delta_mu << "\t" << DD2[i][j] << "\t" << RR2[i][j] << "\t" << DR2[i][j] << endl;
	

	ddr2.close();

	ofstream ddr3(string("/storage/DATA-03/astrorm3/Users/rpandit/output/2pcf/red/z_0/pairs/41-200/ddrrdr_spar_sperp_12x12_50_200Mpc_run"+number+".dat").c_str());

	for(int i=1; i<=sbins; i++)

		for(int j=1.; j<=sbins; j++) 
		  
		ddr3 << setiosflags(ios::uppercase) << std::scientific << smin+(i-0.5)*delta_s << "\t" << smin+(j-0.5)*delta_s << "\t" << DD3[i][j] << "\t" << RR3[i][j] << "\t" << DR3[i][j] << endl;
	

	ddr3.close();



//================================== Free-ing up memosy for arrays =============================//
 	free_dvector(xd,0,n);
	free_dvector(yd,0,n);
	free_dvector(zd,0,n);
    	free_dvector(vx,0,n);

	/*free_dvector(DD, 0, sbins);
	free_dvector(RR, 0, sbins);
	free_dvector(DR, 0, sbins);
	
	free_dmatrix(DD2, 0, sbins, 0, sbins);
	free_dmatrix(RR2, 0, sbins, 0, sbins);
	free_dmatrix(DR2, 0, sbins, 0, sbins);
	free_dmatrix(DD3, 0, sbins, 0, sbins);
	free_dmatrix(RR3, 0, sbins, 0, sbins);
	free_dmatrix(DR3, 0, sbins, 0, sbins);
*/
	
        cout <<  string("Computed pairs for run "+number+" successfully").c_str()  <<   " :) " << endl;



	cout << "----------------------------------------------------------------------------------------------------" << endl;

	/* relevant code to benchmark in here */

	endTime = time(NULL);

	totalTime = endTime - startTime;
	
	std::cout << "Runtime: " << totalTime << " seconds.";


	
	
	}	
	
	free_dvector(xr,0,m);
	free_dvector(yr,0,m);
	free_dvector(zr,0,m);

//============================================ Writing LogFiles ================================//

	string logdate;
	stringstream convert_date_time;
	convert_date_time << time(NULL);
	logdate = convert_date_time.str();
	
	ofstream logfile(string("/storage/DATA-03/astrorm3/Users/rpandit/logfiles/redcorr_z_0_41-200"+logdate+".log").c_str());

	logfile <<"Log File with parameters for redshift space pair calculator\n\n\n";
	logfile <<"We are at redshift z = :"<< Z << "\n";
	logfile <<"No. of sbins :"  << "\t" << Nbins << "\n";
	logfile <<"No. of mubins :" << "\t" << mubins << "\n";
	logfile <<"H(z) in [h Mpc/s/km^2] :" << "\t" << Hofz << "\n";
	logfile <<"Min separation in [h^{-1} Mpc] " << "\t" << smin << "\n";
	logfile <<"Min separation in [h^{-1} Mpc] " << "\t" << smax << "\n";
	 
	logfile <<"------------------------------------------------------------------ " << endl;
	
		  
	logfile.close();	





//============================================ DONE!! ==========================================//	
	return 0;

}
