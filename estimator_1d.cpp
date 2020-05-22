/* THis calculates 1D correlation function estimators from the pairs DD RR and DR*/
#include <cmath>
#include <iostream>
#include <valarray>
#include <fstream>
#include <cassert>
#include <iomanip>
#include <sstream>
#include <string>
#include <string.h>


using namespace std;

const int N=160;
const int Nbins=50; 
int m = 1e+6;

int from = 41;
int to = 200;

int main()
{

	// num is the number of box under analysis	

	for (int num=from; num<=to; num++){

	string number;

	stringstream convert;

	convert <<  num ;

    	number = convert.str();

	int numLines = 0;

	ifstream in1(string("/storage/DATA-03/astrorm3/Users/rpandit/sim/zHorizon/gc/anukruti/z_0.5/41-200/run"+number+".dat").c_str());

	assert(in1);
	std::string unused;
	while ( std::getline(in1, unused) )
	++numLines;
	in1.close();	
	
         int n = numLines; 

	// Here I do the normalisation of DD, RR and DR and measure the estimators.

	double ratio=(double) m /(double) n;
	double Npairs_data=0.5*((double) n)*((double)n-1.);
	double Npairs_random=0.5*((double) m)*((double)m-1.);
	double Npairs_data_random=((double) n)*((double) m);

	
	double r[Nbins+1]={0}, dd[Nbins+1]={0}, rr[Nbins+1]={0}, dr[Nbins+1]={0}, c1[Nbins+1]={0}, c2[Nbins+1]={0}, c3[Nbins+1]={0}, c4[Nbins+1]={0};


	ifstream in2 (string("/storage/DATA-03/astrorm3/Users/rpandit/output/2pcf/real/z_0.5/pairs/41-200/ddrrdr_r_50x50_0_200Mpc_run"+number+".dat").c_str());
		assert(in2);
		for (int i=1; i<=Nbins; i++)
		    in2 >> r[i] >> dd[i] >> rr[i] >> dr[i] ;
		
		in2.close();





	for(int i=1; i<=Nbins; i++){

				dd[i]/=Npairs_data;
				rr[i]/=Npairs_random;
				dr[i]/=Npairs_data_random;

			}
	



	for(int i=1; i<=Nbins; i++){

			    c1[i]=dd[i]/rr[i]-1.;				// Peebles & Hauser
				c2[i]=dd[i]/dr[i]-1.;				// Davis & Peebles
				c3[i]=(dd[i]*rr[i])/pow(dr[i],2.)-1.;		// Hamilton
				c4[i]=(dd[i]-2.*dr[i]+rr[i])/rr[i];		// Landy & Szalay
	}




	ofstream out (string("/storage/DATA-03/astrorm3/Users/rpandit/output/2pcf/real/z_0.5/estimators/41-200/estimator_1d_50bins_50_200Mpc_run"+number+".dat").c_str());

	for(int i=1; i<=Nbins; i++){
	
		out << std::scientific << setiosflags(ios::uppercase) <<  r[i] << "\t" << c1[i]<< "\t" << c2[i]<< "\t" << c3[i]<< "\t" << c4[i] << endl;
	
		}

	
	out.close();

}

return 0;

}
