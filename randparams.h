//200Mpc

const int sbins=12; //number of lin-bins
const int mubins=50; // number of mubins 
double rmax=200.00, rmin=40.;
double smax=200.00, smin=40.;
double delta_s= (double) (smax-smin)/sbins;
double mu_max=1., mu_min=0.;
double delta_mu= (double) (mu_max-mu_min)/sbins;
const int L=1500;     //box side in Mpc/h
const int m=1e+06;   //number of random
const int H_0=100;   // hubble parameter h Mpc^{-1} km/s
const float omega_m=0.25; // matter density
