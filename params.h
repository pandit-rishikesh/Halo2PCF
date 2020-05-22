//200Mpc

const int Nbins=50; //number of lin-bins
const int rbins = 50; // real space bins
const int sbins = 50; // redshift space bins
const int mubins=50; // number of mubins
double rmax=200.00, rmin=0.;
double smax=200.00, smin=0.;
double delta_s= (double) (smax-smin)/sbins;
double delta_r = delta_s;
double mu_max=1., mu_min=0.;
double delta_mu= (double) (mu_max-mu_min)/mubins;
const int L=1500;     //box side in Mpc/h


const int m=1e+06;   //number of random at z = 0 , 0.5

const int H_0=100;   // h km s^-1 Mpc^-1
const int hubble = 0.7;  // hubble parameter
const float omega_m=0.25; // matter density
const float omega_lambda = 1 - omega_m;
const float Z = 0.5;
const float Hofz = H_0 * sqrt(omega_m*pow((1+Z), 3.)+omega_lambda);
const float gadget_normalize = sqrt(1/(1+Z));

