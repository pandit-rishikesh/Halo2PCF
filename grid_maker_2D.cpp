#include<iostream>
#include<cmath>
#include<valarray>
#include<fstream>
#include<cassert>
#include "/storage/DATA-03/astrorm3/Users/rpandit/logistics/headers/params.h"
#include<sstream>
#include<string>
#include<iomanip>

using namespace std;


int from =41;
int to = 200;

int main (){

	double aux;

	double rx[Nbins+1][Nbins+1]={0}, ry[Nbins+1][Nbins+1]={0};

	double delta=(smax-smin)/Nbins;

	double Sum_DD[Nbins+1][Nbins+1]={0}, Sum_RR[Nbins+1][Nbins+1]={0}, Sum_DR[Nbins+1][Nbins+1]={0};


	double dd[Nbins+1][Nbins+1]={0}, rr[Nbins+1][Nbins+1]={0}, dr[Nbins+1][Nbins+1]={0};

	double c1[Nbins+1][Nbins+1]={0}, c2[Nbins+1][Nbins+1]={0}, c3[Nbins+1][Nbins+1]={0}, c4[Nbins+1][Nbins+1]={0};




for(int xxx=from; xxx<=to; xxx++){

 string number;

 stringstream convert;

 convert << xxx;

 number = convert.str();

 int numLines = 0;

 ifstream in(string("/storage/DATA-03/astrorm3/Users/rpandit/output/2pcf/red/z_0.5/normalised_pairs/41-200/normalised_pairs_spar_sperp_50x50_0_200Mpc_run"+number+".dat").c_str());
 assert(in);


	for (int i=1; i<=Nbins; i++){

		for(int j=1; j<=Nbins; j++){

			in >> i >> j >> dd[i][j] >> dr[i][j] >> rr[i][j];

		}
	}


 in.close();



	for(int i=1; i<=Nbins; i++){

 		for(int j=1; j<=Nbins; j++){

 			Sum_DD[i][j]+=dd[i][j];
 			Sum_RR[i][j]+=rr[i][j];
			Sum_DR[i][j]+=dr[i][j];

		}

	}

 }

//cout << Sum_DD[1][1] << endl;
//return 0;



cout << "Building estimators" << endl;


	for (int i=1; i<=Nbins; i++)
		  for (int j=1; j<=Nbins; j++){


		c1[i][j]=(double)Sum_DD[i][j]/(double)Sum_RR[i][j]-1.0;
		c2[i][j]=(double)Sum_DD[i][j]/(double)Sum_DR[i][j]-1.0;
		c3[i][j]=((double)Sum_DD[i][j])*((double)Sum_RR[i][j])/pow((double)Sum_DR[i][j], 2.)-1.0;
		c4[i][j]=((double)Sum_DD[i][j]-2.*(double)Sum_DR[i][j]+(double)Sum_RR[i][j])/((double)Sum_RR[i][j]);

		}



	cout<<"Writing output"<<endl;

	ofstream out1("/storage/DATA-03/astrorm3/Users/rpandit/output/2pcf/red/z_0.5/grid/grid_z1.dat");
	out1.precision(12);

	for(int i=1.; i<=Nbins; i++){
		for (int j=1.; j<=Nbins; j++)

	out1 << setiosflags(ios::uppercase) << std::scientific << i << "\t"<< j << "\t" << smin+(i-0.5)*delta << "\t" << smin+(j-0.5)*delta <<  "\t" << c1[i][j] << "\t" << c2[i][j] << "\t" << c3[i][j] << "\t" << c4[i][j] << endl;
	}

	out1.close();

	return 0 ;
}
