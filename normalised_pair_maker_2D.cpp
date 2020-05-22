/* THis calculates correlation function estimators from the pairs DD RR and DR*/
#include <cmath>
#include <iostream>
#include <valarray>
#include <fstream>
#include <cassert>
#include "/storage/DATA-03/astrorm3/Users/rpandit/logistics/headers/params.h"
#include <iomanip>
#include <sstream>
#include <string>
#include <string.h>


using namespace std;

int from = 41;
int to = 200;


// Here I do the normalisation of DD, RR and DR and measure the estimators.
int main()
{

	for (int xxx=from; xxx<=to; xxx++){

	string number;

	stringstream convert;

	convert <<  xxx ;

    	number = convert.str();

	int numLines = 0;

	ifstream in1(string("/storage/DATA-03/astrorm3/Users/rpandit/sim/zHorizon/gc/anukruti/z_1/41-200/run"+number+".dat").c_str());

	assert(in1);
	std::string unused;
	while ( std::getline(in1, unused) )
	++numLines;
	in1.close();

         int n = numLines;

	double Npairs_data=0.5*((double) n)*((double)n-1.);
	double Npairs_random=0.5*((double) m)*((double)m-1.);
	double Npairs_data_random=((double) n)*((double) m);


	double sigma[Nbins+1][Nbins+1]={0}, pi[Nbins+1][Nbins+1]={0}, Np[Nbins+1][Nbins+1]={0}, dd[Nbins+1][Nbins+1]={0}, rr[Nbins+1][Nbins+1]={0}, dr[Nbins+1][Nbins+1]={0}, c1[Nbins+1][Nbins+1]={0}, c2[Nbins+1][Nbins+1]={0}, c3[Nbins+1][Nbins+1]={0}, c4[Nbins+1][Nbins+1]={0};



	ifstream in2 (string("/storage/DATA-03/astrorm3/Users/rpandit/output/2pcf/red/z_1/pairs/41-200/ddrrdr_spar_sperp_50x50_0_200Mpc_run"+number+".dat").c_str());
		assert(in2);
		for (int i=1; i<=Nbins; i++)
		  for (int j=1; j<=Nbins; j++)
		    in2 >> sigma[i][j] >> pi[i][j] >> dd[i][j] >> rr[i][j] >> dr[i][j] ;

		in2.close();


	//in.close();
	//for(int i=1; i<=nbins; i++)Np[i]=dd[i];
	for(int i=1; i<=Nbins; i++){
		for (int j=1; j<=Nbins; j++){

				dd[i][j]/=Npairs_data;
				rr[i][j]/=Npairs_random;
				dr[i][j]/=Npairs_data_random;

			}

	}




	ofstream pairs (string("/storage/DATA-03/astrorm3/Users/rpandit/output/2pcf/red/z_1/normalised_pairs/41-200/normalised_pairs_spar_sperp_50x50_0_200Mpc_run"+number+".dat").c_str());


	for(int i=1; i<=Nbins; i++){
		for (int j=1; j<=Nbins; j++){

		pairs << i << "\t" << j << "\t" << dd[i][j]<< "\t" <<dr[i][j]<< "\t" <<rr[i][j]<< endl;

		}

	}

	pairs.close();

}

return 0;

}
