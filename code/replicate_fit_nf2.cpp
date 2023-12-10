//C++ standard headers
#include <iostream>
#include <climits>
#include <vector>

//ROOT headers for plotting
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TGraphErrors.h"

//ROOT headers for math
#include "TMath.h"
#include "TComplex.h"
#include "TMatrixDSym.h"
#include "TMatrixD.h"

//ROOT headers for minimizing
#include "TFitter.h"
#include "TVirtualFitter.h"
#include "TMinuit.h"

//in-house header files
#include "BetaFn.h"

using namespace std;

#define NFLAV 2

#define NOTPADE 0
#define PADE 1

#define DIRECT 0
#define RECIPROCAL 1

Double_t RATIO=1.00; //for testing the sensitivity to scale

#define NDATA 15
#define NBETA 4

//data from Table 1 of arxiv:hep-ph/0502212v2
Double_t beta[NDATA]         ={5.20    ,5.20    ,5.20    ,5.20    ,5.20    ,5.25    ,5.25    ,5.25    ,
                               5.29    ,5.29    ,5.29    ,5.29    ,5.40    ,5.40    ,5.40             };
Double_t r0_over_a[NDATA]    ={4.077   ,4.754   ,5.041   ,5.250   ,5.320   ,4.737   ,5.138   ,5.532   ,
                               4.813   ,5.227   ,5.566   ,5.880   ,6.092   ,6.381   ,6.714            };
Double_t r0_over_a_err[NDATA]={0.070   ,0.045   ,0.053   ,0.075   ,0.095   ,0.050   ,0.055   ,0.040   ,
                               0.082   ,0.075   ,0.064   ,0.100   ,0.067   ,0.053   ,0.064            };
Double_t plaquette[NDATA]    ={0.528994,0.533670,0.536250,0.537070,0.537670,0.538770,0.541150,0.543135,
                               0.542400,0.545520,0.547094,0.548286,0.559000,0.560246,0.561281         };
Double_t plaquette_err[NDATA]={0.000058,0.000040,0.000030,0.000100,0.000030,0.000041,0.000030,0.000015,
                               0.000050,0.000029,0.000023,0.000057,0.000019,0.000010,0.000008         };

//for nf=2 usage
Double_t kappa[NDATA]        ={0.1342  ,0.1350  ,0.1355  ,0.13565 ,0.1358  ,0.1346  ,0.1352  ,0.13575 ,
                               0.1340  ,0.1350  ,0.1355  ,0.1359  ,0.1350  ,0.1356  ,0.1361           };
Double_t kappa_c[NDATA]      ={0.136008,0.136008,0.136008,0.136008,0.136008,0.136250,0.136250,0.136250,
                               0.136410,0.136410,0.136410,0.136410,0.136690,0.136690,0.136690         };
Double_t kappa_c_err[NDATA]  ={0.000015,0.000015,0.000015,0.000015,0.000015,0.000007,0.000007,0.000007,
                               0.000009,0.000009,0.000009,0.000009,0.000022,0.000022,0.000022         };

// -- for chiral limit extrapolation
Double_t quarkmass[NDATA];  //acting as another name for am_q
Double_t quarkmass_err[NDATA];

// -- for continuum limit extrapolation
Double_t beta_extrap[NBETA]={5.20,5.25,5.29,5.40};
Double_t r0_over_a_extrap[NBETA];
Double_t r0_over_a_extrap_err[NBETA];
Double_t plaquette_extrap[NBETA];
Double_t plaquette_extrap_err[NBETA];
Double_t a_over_r0_sqr[NBETA];
Double_t a_over_r0_sqr_err[NBETA];

//results need to be global for minuit to work properly
Double_t methodI_r0Lambda[NBETA]; 
Double_t methodI_r0Lambda_err[NBETA];

//results need to be global for minuit to work properly
Double_t methodII_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodII_r0Lambda_err[NDATA];

//results need to be global for minuit to work properly
Double_t methodIIP_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodIIP_r0Lambda_err[NDATA];

//results need to be global for minuit to work properly
Double_t methodIII_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodIII_r0Lambda_err[NDATA];

//results need to be global for minuit to work properly
Double_t methodIIIP_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodIIIP_r0Lambda_err[NDATA];

//analytic formula for mu_over_lambda using beta function coefficients
Double_t mu_over_lambda(Double_t x, BetaFn bf){
  std::vector<TComplex> P(bf.nloops()-1,0.);
  std::vector<TComplex> roots=bf.nonzeroroots();
  for(Int_t i=0;i<bf.nloops()-1;++i){
    P[i]=(-1./(bf[bf.nloops()-1]*roots[i]*roots[i]));
    for(Int_t j=0;j<bf.nloops()-1;++j){
      if(i!=j){
        P[i]*=( 1./(roots[i]-roots[j]) );
      }
    }
  }
  
  TComplex term1=0.,term2=0.,ln_arg=1.;
  if(bf.variable()==2){ //using g as a variable
    term1=1./(2.*bf[0]*x*x);
  }
  else{ //using alpha or a as a variable
    term1=1./(2.*bf[0]*x);
  }
  if(bf.variable()==2){  //use g as a variable
    term2=( bf[1]/(2.*bf[0]*bf[0]) )*TMath::Log( bf[0]*x*x );
  }
  else{  //use alpha or a as a variable
    term2=( bf[1]/(2.*bf[0]*bf[0]) )*TMath::Log( bf[0]*x );
  }
  for(Int_t i=0;i<bf.nloops()-1;++i){
    if(bf.variable()>1){
      ln_arg*=TComplex::Power( ((TComplex)1.)-x*x/roots[i],P[i]/2. );
    }
    else{
      ln_arg*=TComplex::Power( ((TComplex)1.)-x/roots[i],P[i]/2. );
    }
  }
  if(bf.nloops()<2){
    return TMath::Exp(term1);
  }
  return TComplex::Exp( term1+term2+TComplex::Log(ln_arg)).Re();
};

// equation (13) from hep-ph/0502212
Double_t mu_over_lambda_Roger(Double_t g,BetaFn bf){
  Double_t fac1=TMath::Exp(1./(2.*bf[0]*g*g));
  Double_t fac2=TMath::Power( bf[0]*g*g/(1.+(bf[1]/bf[0]-bf[2]/bf[1])*g*g) , bf[1]/(2.*bf[0]*bf[0]) );
  return fac1*fac2;
};

// here follows many useful variables relating g_msbar to g0 and g_boost
Double_t u0(Int_t ind){
  if( (ind<0) || (ind>(NBETA-1)) ) {cout<<"\n u0: Bad index ["<<ind<<"]\n\n"; exit(-1);}
  return TMath::Power(plaquette_extrap[ind],0.25);
};
Double_t u0_err(Int_t ind){
  if( (ind<0) || (ind>(NBETA-1)) ){cout<<"\n u0_err: Bad index ["<<ind<<"]\n\n"; exit(-1);}
  Double_t partial_plaquette=0.25*TMath::Power(plaquette_extrap[ind],-0.75);
  Double_t term1=TMath::Power(partial_plaquette,2.)*TMath::Power(plaquette_extrap_err[ind],2.);
  return TMath::Sqrt(term1);
};

Double_t c_sw(Int_t ind){
  if( (ind<0) || (ind>(NBETA-1)) ){cout<<"\n c_sw: Bad index ["<<ind<<"]\n\n"; exit(-1);}
  Double_t gsqr=6./beta_extrap[ind];
  Double_t numer=1.-0.454*gsqr-0.175*gsqr*gsqr+0.012*gsqr*gsqr*gsqr+0.045*gsqr*gsqr*gsqr*gsqr;
  Double_t denom=1.-0.720*gsqr;
  return numer/denom;
};
Double_t c_sw_err(Int_t ind){
  return 0.;  //there is no error associated to c_sw, since there is no error associated with beta, and 
              //the rational form of c_sw is taken as the definition, see Advanced Lattice QCD by M. Luscher
};

Double_t c_sw_boosted(Int_t ind){
  return c_sw(ind)*TMath::Power( u0(ind) , 3. );
};
Double_t c_sw_boosted_err(Int_t ind){
  Double_t partial_c_sw=TMath::Power( u0(ind) , 3. );
  Double_t partial_u0=3.*c_sw(ind)*TMath::Power( u0(ind) , 2. );
  Double_t corr_c_sw__u0=0.;      // no correlations

  Double_t term1=TMath::Power(partial_c_sw,2.)*TMath::Power(c_sw_err(ind),2.);
  Double_t term2=TMath::Power(partial_u0,2.)*TMath::Power(u0_err(ind),2.);
  Double_t termcov=2.*partial_c_sw*partial_u0*c_sw_err(ind)*u0_err(ind)*corr_c_sw__u0;

  return TMath::Sqrt(term1+term2+termcov);
};

Double_t amq(Int_t ind){
  if( (ind<0) || (ind>(NDATA-1)) ){cout<<"\n c_sw: Bad index ["<<ind<<"]\n\n"; exit(-1);}
  return 0.5*( (1./kappa[ind])-(1./kappa_c[ind]) );
};
Double_t amq_err(Int_t ind){  // no correlations
  if( (ind<0) || (ind>(NDATA-1)) ){cout<<"\n u0: Bad index ["<<ind<<"]\n\n"; exit(-1);}
  Double_t partial_kappa_c=0.5*TMath::Power( kappa_c[ind] , -2. );
  Double_t term1=TMath::Power(partial_kappa_c,2.)*TMath::Power(kappa_c_err[ind],2.);
  return TMath::Sqrt(term1);
};

Double_t t1_lat_chirallimit(Int_t ind){
  Double_t csw=c_sw(ind);
  return ( 0.4682013-NFLAV*(0.0066960-0.0050467*csw+0.0298435*csw*csw) );
};
Double_t t1_lat_chirallimit_err(Int_t ind){
  //c_sw has no error. However, there should be truncation errors on the order of c_sw^3, but we don't know the coefficient,
  //so we'll leave that error to be soaked up in the systematic error
  Double_t csw=c_sw(ind);
  Double_t partial_c_sw=-NFLAV*(-0.0050467+2.*0.0298435*csw);
  Double_t term1=TMath::Power(partial_c_sw,2.)*TMath::Power(c_sw_err(ind),2.);
  return TMath::Sqrt(term1);
};

/*  //you shouldn't use this since you are taking the continuum limit *in* the chiral limit
Double_t t1_lat(Int_t ind){
  Double_t csw=c_sw(ind);
  Double_t _amq=amq(ind);
  return ( t1_lat_chirallimit(ind) - NFLAV*_amq*(-0.0272837+0.0223503*csw-0.0070667*csw*csw));
};
Double_t t1_lat_err(Int_t ind){
  //c_sw has no error. However, there should be truncation errors on the order of c_sw^3, but we don't know the coefficient,
  //so we'll leave that error to be soaked up in the systematic error
  Double_t csw=c_sw(ind);
  Double_t _amq=amq(ind);

  Double_t partial_c_sw=-NFLAV*(-0.0050467+2.*0.0298435*csw+_amq*(-0.0223503+2.*0.0070667*csw));
  Double_t partial_amq=-NFLAV*(-0.0272837+0.0223503*csw-0.0070667*csw*csw);
  Double_t corr_c_sw__amq=0.;      // no correlations

  Double_t term1=TMath::Power(partial_c_sw,2.)*TMath::Power(c_sw_err(ind),2.);
  Double_t term2=TMath::Power(partial_amq,2.)*TMath::Power(amq_err(ind),2.);
  Double_t termcov=2.*partial_c_sw*partial_amq*c_sw_err(ind)*amq_err(ind)*corr_c_sw__amq;

  return TMath::Sqrt(term1+term2+termcov);
};
*/

Double_t t2_lat(Int_t ind){
  Double_t csw=c_sw(ind);
  return ( 0.0556675 - NFLAV*(0.002600 + 0.000155*csw - 0.012834*csw*csw - 0.000474*csw*csw*csw - 0.000104*csw*csw*csw*csw) );
};
Double_t t2_lat_err(Int_t ind){
  //c_sw has no error. However, there should be trauncation errors on the order of c_sw^3, but we don't know the coefficients,
  //so we'll leave that error to be soaked up in the systematic error
  Double_t csw=c_sw(ind);
  Double_t partial_c_sw=-NFLAV*(0.000155- 2.*0.012834*csw - 3.*0.000474*csw*csw - 4.*0.000104*csw*csw*csw);
  Double_t term1=TMath::Power(partial_c_sw,2.)*TMath::Power(c_sw_err(ind),2.);
  return TMath::Sqrt(term1);   
};

Double_t t1_boost(Int_t ind,Double_t val=1.){
  Double_t cswb=val;
  if(ind>=0){cswb=c_sw_boosted(ind);}
  return ( 0.1348680 - NFLAV*(0.0066960-0.0050467*cswb+0.0298435*cswb*cswb ) );
};
Double_t t1_boost_err(Int_t ind,Double_t val=1.){
  Double_t cswb;
  if(ind<0){cswb=val; return 0.;}
  else{cswb=c_sw_boosted(ind);}
  Double_t partial_c_sw_boosted=-NFLAV*(-0.0050467+2.*0.0298435*cswb);
  Double_t term1=TMath::Power(partial_c_sw_boosted,2.)*TMath::Power(c_sw_boosted_err(ind),2.);
  return TMath::Sqrt(term1);
};

Double_t t2_boost(Int_t ind,Double_t val=1.){
  Double_t cswb;
  if(ind<0){cswb=val;}
  else{cswb=c_sw_boosted(ind);}
  return ( 0.0217565 - NFLAV*(0.000753-0.001053*cswb+0.000498*cswb*cswb-0.000474*cswb*cswb*cswb-0.000104*cswb*cswb*cswb*cswb) );
};
Double_t t2_boost_err(Int_t ind,Double_t val=1.){
  Double_t cswb;
  if(ind<0){cswb=val; return 0.;}
  else{cswb=c_sw_boosted(ind);}
  Double_t partial_c_sw_boosted=-NFLAV*(0.001053-2.*0.000498*cswb-3.*0.000474*cswb*cswb-4.*0.000104*cswb*cswb*cswb);
  Double_t term1=TMath::Power(partial_c_sw_boosted,2.)*TMath::Power(c_sw_boosted_err(ind),2.);
  return TMath::Sqrt(term1);
};

Double_t g_lat(Int_t ind){
  if( (ind<0) || (ind>(NBETA-1)) ){cout<<"\n g_boost: Bad index \n\n"; exit(-1);}
  return TMath::Sqrt( 6./beta_extrap[ind] );
};
Double_t g_lat_err(Int_t ind){
  return 0.;
};
Double_t g_boost(Int_t ind){
  if( (ind<0) || (ind>(NBETA-1)) ){cout<<"\n g_boost: Bad index \n\n"; exit(-1);}
  return TMath::Sqrt( 6./(beta_extrap[ind]*plaquette_extrap[ind]) );
};
Double_t g_boost_err(Int_t ind){
  if( (ind<0) || (ind>(NBETA-1)) ){cout<<"\n g_boost_err: Bad index \n\n"; exit(-1);}
  Double_t partial_plaquette=-6./(beta_extrap[ind]*plaquette_extrap[ind]*plaquette_extrap[ind]);
  Double_t term1=TMath::Power(partial_plaquette,2.)*TMath::Power(plaquette_extrap_err[ind],2.);
  return TMath::Sqrt( term1 );
};

// boosted beta coeffs for method II
Double_t b2_boost_II(Int_t ind){
  BetaFn bf_ms(4,NFLAV,2);
  return bf_ms[2]+bf_ms[1]*t1_boost(ind)-bf_ms[0]*t2_boost(ind);
};
Double_t b2_boost_II_err(Int_t ind){
  BetaFn bf_ms(4,NFLAV,2);
  Double_t partial_t1_boost=bf_ms[1];
  Double_t partial_t2_boost=-bf_ms[0];

  Double_t corr_t1_boost__t2_boost=1.;
  
  Double_t term1=TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.);
  Double_t term2=TMath::Power(partial_t2_boost,2.)*TMath::Power(t2_boost_err(ind),2.);

  Double_t termcorr1=2.*partial_t1_boost*partial_t1_boost*g_boost_err(ind)*t1_boost_err(ind)*corr_t1_boost__t2_boost;

  return TMath::Sqrt(term1+term2+termcorr1);
};
Double_t b3_boost_II(Int_t ind){
  BetaFn bf_ms(4,NFLAV,2);
  return b2_boost_II(ind)*b2_boost_II(ind)/bf_ms[1];;
};
Double_t b3_boost_II_err(Int_t ind){
  Double_t partial_b2_boost=2.*b3_boost_II(ind)/b2_boost_II(ind);
  Double_t term1=TMath::Power(partial_b2_boost,2.)*TMath::Power(b2_boost_II_err(ind),2.);
  return TMath::Sqrt(term1);
};

// boosted beta coeffs for method III
Double_t b2_boost_III(Int_t ind){
  BetaFn bf_ms(4,NFLAV,2);
  return bf_ms[2]+bf_ms[1]*t1_boost(-1,1.)-bf_ms[0]*t2_boost(-1,1.)-bf_ms[0]*(0.2659-0.25)*(-NFLAV*(-0.0050467+2.*0.0298435));
};
Double_t b2_boost_III_err(Int_t ind){
  BetaFn bf_ms(4,NFLAV,2);
  Double_t partial_t1_boost=bf_ms[1];
  Double_t partial_t2_boost=-bf_ms[0];

  Double_t corr_t1_boost__t2_boost=1.;
  
  Double_t term1=TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.);
  Double_t term2=TMath::Power(partial_t2_boost,2.)*TMath::Power(t2_boost_err(ind),2.);

  Double_t termcorr1=2.*partial_t1_boost*partial_t1_boost*g_boost_err(ind)*t1_boost_err(ind)*corr_t1_boost__t2_boost;

  return TMath::Sqrt(term1+term2+termcorr1);
};
Double_t b3_boost_III(Int_t ind){
  BetaFn bf_ms(4,NFLAV,2);
  return b2_boost_III(ind)*b2_boost_III(ind)/bf_ms[1];;
};
Double_t b3_boost_III_err(Int_t ind){
  Double_t partial_b2_boost=2.*b3_boost_III(ind)/b2_boost_III(ind);
  Double_t term1=TMath::Power(partial_b2_boost,2.)*TMath::Power(b2_boost_III_err(ind),2.);
  return TMath::Sqrt(term1);
};


//fitting r0Lambda vs (a/r0)^2
//with y=a+bx
Double_t fitline(Double_t x, Double_t a, Double_t b){
  return a+b*x;
}
Double_t TF1_fitline(Double_t* x, Double_t* par){
  return fitline(x[0],par[0],par[1]);
}

//fitting r0Lambda vs (a/r0)^2
//with y=a+bx+cx^2
Double_t fit_parabola(Double_t x, Double_t a, Double_t b, Double_t c){
  return a+b*x+c*x*x;
}
Double_t fit_parabola_deriv(Double_t x, Double_t a, Double_t b, Double_t c){
  return b+2.*c*x;
}
Double_t TF1_fit_parabola(Double_t* x, Double_t* par){
  return fit_parabola(x[0],par[0],par[1],par[2]);
}
Double_t fit_global(Double_t _amq, Double_t beta, Double_t* par){
  return TMath::Exp( (par[0]+par[1]*beta) 
                   + (par[2]+par[3]*beta+par[4]*beta*beta)*_amq
                   + (par[5]+par[6]*beta+par[7]*beta*beta)*_amq*_amq );
};
Double_t fit_global_deriv(Double_t _amq, Double_t beta, Double_t* par){
  return fit_global(_amq,beta,par)*( 
                  (par[2]+par[3]*beta+par[4]*beta*beta) 
             + 2.*(par[5]+par[6]*beta+par[7]*beta*beta)*_amq );
};
Double_t TF1_fit_global(Double_t* x, Double_t* par){
  return fit_global(x[0],par[8],par);
}

// ------------- chiral extrapolation for all methods ---------------

// you're using a chi-squared fit, which means you're ignoring the error
// on the x-variable. Try to find a better way

// *Okay, we're going to implement the "effective variance" approach. Here we use
// \chi^2 = sum_i^N [y_i - f(x_i)]^2/[\sigma_y^2+(\sigma_x*f'(x_i))^2]

// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_chiral_plaq520_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=0;i<=4;++i) { 
      num=TMath::Power( plaquette[i] - fit_parabola(quarkmass[i],p[0],p[1],p[2]) , 2. );
      den=TMath::Power(plaquette_err[i],2.)+TMath::Power(quarkmass_err[i]*fit_parabola_deriv(quarkmass[i],p[0],p[1],p[2]),2.);
      chi2 += num/den;
  }
  fval = chi2;
}
// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_chiral_plaq525_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=5;i<=7;++i) { 
      num=TMath::Power( plaquette[i] - fit_parabola(quarkmass[i],p[0],p[1],p[2]) , 2. );
      den=TMath::Power(plaquette_err[i],2.)+TMath::Power(quarkmass_err[i]*fit_parabola_deriv(quarkmass[i],p[0],p[1],p[2]),2.);
      chi2 += num/den;
  }
  fval = chi2;
}
// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_chiral_plaq529_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=8;i<=11;++i) { 
      num=TMath::Power( plaquette[i] - fit_parabola(quarkmass[i],p[0],p[1],p[2]) , 2. );
      den=TMath::Power(plaquette_err[i],2.)+TMath::Power(quarkmass_err[i]*fit_parabola_deriv(quarkmass[i],p[0],p[1],p[2]),2.);
      chi2 += num/den;
  }
  fval = chi2;
}
// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_chiral_plaq540_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=12;i<=14;++i) { 
      num=TMath::Power( plaquette[i] - fit_parabola(quarkmass[i],p[0],p[1],p[2]) , 2. );
      den=TMath::Power(plaquette_err[i],2.)+TMath::Power(quarkmass_err[i]*fit_parabola_deriv(quarkmass[i],p[0],p[1],p[2]),2.);
      chi2 += num/den;
  }
  fval = chi2;
}

// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_chiral_r0a_global_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=0;i<NDATA;++i) { 
      num=TMath::Power( r0_over_a[i] - fit_global(quarkmass[i],beta[i],p) , 2. );
      den=TMath::Power(r0_over_a_err[i],2.)+TMath::Power(quarkmass_err[i]*fit_global_deriv(quarkmass[i],beta[i],p),2.);
      chi2 += num/den;
  }
  fval = chi2;
}

// -------------------------- method I ---------------------
Double_t g_ms_I(Int_t ind,Int_t opt=0){
// opt==0 -- use direct expansion of g_MS ie g_MS^2=g_boost^2+d_1*g_boost^4+...
// opt==1 -- use reciprocal expansion of g_MS ie 1/g_MS^2=1/g_boost^2+...  -- this one agrees with the r0Lambda values in the paper
  BetaFn bf(4,NFLAV,2);
  if(opt==RECIPROCAL){
    return TMath::Power( TMath::Power( g_boost(ind) , -2. ) + (2.*bf[0]*TMath::Log(RATIO))
                        + (2.*bf[1]*TMath::Log(RATIO)-t2_boost(ind)+bf[1]*t1_boost(ind)/bf[0])*TMath::Power(g_boost(ind),2.) , -0.5 );
  }
  // opt==DIRECT
  return TMath::Sqrt( TMath::Power( g_boost(ind) , 2. ) - (2.*bf[0]*TMath::Log(RATIO))*TMath::Power( g_boost(ind) , 4. )
                        + (4.*bf[0]*bf[0]*TMath::Log(RATIO)*TMath::Log(RATIO)-2.*bf[1]*TMath::Log(RATIO)+t2_boost(ind)-bf[1]*t1_boost(ind)/bf[0])*TMath::Power(g_boost(ind),6.) );
};

Double_t g_ms_I_err(Int_t ind,Int_t opt=0){
  BetaFn bf(4,NFLAV,2);
  // opt==0 -- use direct expansion of g_MS ie g_MS^2=g_boost^2+d_1*g_boost^4+...
  // opt==1 -- use reciprocal expansion of g_MS ie 1/g_MS^2=1/g_boost^2+...  -- this one agrees with the r0Lambda values in the paper
  if(opt==RECIPROCAL){

    Double_t prefac=0.25*TMath::Power(g_ms_I(ind,opt),6.);

    Double_t partial_g_boost=-2.*TMath::Power(g_boost(ind),-3.)
                    +2.*(2.*bf[1]*TMath::Log(RATIO)-t2_boost(ind)+bf[1]*t1_boost(ind)/bf[0])*TMath::Power( g_boost(ind) , 1. );
    Double_t partial_t1_boost=bf[1]*TMath::Power( g_boost(ind) , 2. )/bf[0];
    Double_t partial_t2_boost=-TMath::Power( g_boost(ind) , 2. );
 
    Double_t corr_g_boost__t1_boost=1.;
    Double_t corr_g_boost__t2_boost=1.;
    Double_t corr_t1_boost__t2_boost=1.;
  
    Double_t term1=TMath::Power(partial_g_boost,2.)*TMath::Power(g_boost_err(ind),2.);
    Double_t term2=TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.);
    Double_t term3=TMath::Power(partial_t2_boost,2.)*TMath::Power(t2_boost_err(ind),2.);

    Double_t termcorr1=2.*partial_g_boost*partial_t1_boost*g_boost_err(ind)*t1_boost_err(ind)*corr_g_boost__t1_boost;
    Double_t termcorr2=2.*partial_g_boost*partial_t2_boost*g_boost_err(ind)*t2_boost_err(ind)*corr_g_boost__t2_boost;
    Double_t termcorr3=2.*partial_t1_boost*partial_t2_boost*t1_boost_err(ind)*t2_boost_err(ind)*corr_t1_boost__t2_boost;

    return TMath::Sqrt(prefac*(term1+term2+term3+termcorr1+termcorr2+termcorr3));
  }
  Double_t partial_g_boost=2.*g_boost(ind)+6.*(t2_boost(ind)-bf[1]*t1_boost(ind)/bf[0])*TMath::Power( g_boost(ind) , 5. );
  Double_t partial_t1_boost=bf[1]*TMath::Power( g_boost(ind) , 6. )/bf[0];
  Double_t partial_t2_boost=TMath::Power( g_boost(ind) , 6. );
 
  Double_t corr_g_boost__t1_boost=1.;
  Double_t corr_g_boost__t2_boost=1.;
  Double_t corr_t1_boost__t2_boost=1.;
  
  Double_t term1=TMath::Power(partial_g_boost,2.)*TMath::Power(g_boost_err(ind),2.)/TMath::Power(2.*g_ms_I(ind,opt),2.);
  Double_t term2=TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.)/TMath::Power(2.*g_ms_I(ind,opt),2.);
  Double_t term3=TMath::Power(partial_t2_boost,2.)*TMath::Power(t2_boost_err(ind),2.)/TMath::Power(2.*g_ms_I(ind,opt),2.);

  Double_t termcorr1=2.*partial_g_boost*partial_t1_boost*g_boost_err(ind)*t1_boost_err(ind)*corr_g_boost__t1_boost;
  Double_t termcorr2=2.*partial_g_boost*partial_t2_boost*g_boost_err(ind)*t2_boost_err(ind)*corr_g_boost__t2_boost;
  Double_t termcorr3=2.*partial_t1_boost*partial_t2_boost*t1_boost_err(ind)*t2_boost_err(ind)*corr_t1_boost__t2_boost;

  return TMath::Sqrt(term1+term2+term3+termcorr1+termcorr2+termcorr3);
};

Double_t methodI_getr0Lambda(Int_t ind,Int_t opt=0){
  BetaFn bf(4,NFLAV,2);
  Double_t fac1=r0_over_a_extrap[ind];
  Double_t fac2=TMath::Exp(t1_boost(ind)/(2.*bf[0]));
  Double_t fac3=1./mu_over_lambda(g_ms_I(ind,opt),bf);
  cout << "beta 0 , 1 , 2 , 3 :" << "  " << bf[0] << "  " << bf[1] << "  " << bf[2] << "  " << bf[3] << endl;
  return RATIO*fac1*fac2*fac3;
};
Double_t methodI_getr0Lambda_err(Int_t ind,Int_t opt=0){   //did you take the covariance between t1_boost and g_ms_I??? Yep!
  BetaFn bf(4,NFLAV,2);

  Double_t partial_r0_over_a=RATIO*methodI_getr0Lambda(ind)/r0_over_a_extrap[ind];
  Double_t partial_t1_boost=RATIO*methodI_getr0Lambda(ind)/(2.*bf[0]);
  Double_t partial_g_ms_I= -RATIO*methodI_getr0Lambda(ind)/bf.eval(g_ms_I(ind,opt));  

  Double_t corr_r0_over_a__t1_boost=0.;
  Double_t corr_r0_over_a__g_ms_I=0.;
  Double_t corr_t1_boost__g_ms_I=1.;

  Double_t term1=TMath::Power(partial_r0_over_a,2.)*TMath::Power(r0_over_a_extrap_err[ind],2.);
  Double_t term2=TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.);
  Double_t term3=TMath::Power(partial_g_ms_I,2.)*TMath::Power(g_ms_I_err(ind,opt),2.);

  Double_t termcorr1=2.*partial_r0_over_a*partial_t1_boost*r0_over_a_extrap_err[ind]*t1_boost_err(ind)*corr_r0_over_a__t1_boost;
  Double_t termcorr2=2.*partial_r0_over_a*partial_g_ms_I*r0_over_a_extrap_err[ind]*g_ms_I_err(ind,opt)*corr_r0_over_a__g_ms_I;
  Double_t termcorr3=2.*partial_t1_boost*partial_g_ms_I*t1_boost_err(ind)*g_ms_I_err(ind,opt)*corr_t1_boost__g_ms_I;

  return TMath::Sqrt(term1+term2+term3+termcorr1+termcorr2+termcorr3);
};

// *Okay, we're going to implement the "effective variance" approach. Here we use
// \chi^2 = sum_i^N [y_i - f(x_i)]^2/[\sigma_y^2+(\sigma_x*f'(x_i))^2]
// Minuit functions (chi2) p[0]=a  p[1]=b  p[2]=c
void minuit_methodI_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=0;i<NBETA;++i) { //skip the first two due to order(a^2) effects and skip the last
    num=TMath::Power( methodI_r0Lambda[i] - fitline(a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodI_r0Lambda_err[i],2.)+TMath::Power(a_over_r0_sqr_err[i]*p[1],2.);
    chi2 += num/den;
  }
  fval = chi2;
}

// ----------------------- method II --------------------------

Double_t g_ms_II(Int_t ind){
  BetaFn bf(4,NFLAV,2);
  //not using the simple "return g_boost" because we want to look at the influence of RATIO on the outcome of r0Lambda
  return TMath::Sqrt( TMath::Power( g_boost(ind) , 2. ) - (2.*bf[0]*TMath::Log(RATIO))*TMath::Power( g_boost(ind) , 4. )
                        + (4.*bf[0]*bf[0]*TMath::Log(RATIO)*TMath::Log(RATIO)-2.*bf[1]*TMath::Log(RATIO))*TMath::Power(g_boost(ind),6.) );
};
Double_t TF1_partial_b2_boost_II_factor_integrand(Double_t* x,Double_t* par){
// par[0]==ind, par[1]==pade?
  Int_t pade=(Int_t) par[1];
  BetaFn bf_boost(3+pade,NFLAV,2);
  bf_boost[2]=b2_boost_II(par[0]);
  if(pade==1){bf_boost[3]=b3_boost_II(par[0]);}
  return TMath::Power(x[0],7.) * TMath::Power( bf_boost.eval(x[0]), -2. );
}
Double_t partial_b2_boost_II_factor(Int_t ind,Int_t pade){
  TF1 scratch("scratch",TF1_partial_b2_boost_II_factor_integrand,0,100,2);
  scratch.SetParameter(0,ind);
  scratch.SetParameter(1,pade);
  return scratch.Integral(0.,g_boost(ind));
};
Double_t TF1_partial_b3_boost_II_factor_integrand(Double_t* x,Double_t* par){
// par[0]==ind
  BetaFn bf_boost(4,NFLAV,2);
  bf_boost[2]=b2_boost_II(par[0]);
  bf_boost[3]=bf_boost[2]*bf_boost[2]/bf_boost[1];
  return TMath::Power(x[0],9.) * TMath::Power( bf_boost.eval(x[0]), -2. );
}
Double_t partial_b3_boost_II_factor(Int_t ind,Int_t pade){
  if(pade!=1){return 0.;}
  TF1 scratch("scratch",TF1_partial_b3_boost_II_factor_integrand,0,100,1);
  scratch.SetParameter(0,ind);
  return scratch.Integral(0.,g_boost(ind));
};

Double_t methodII_getr0Lambda(Int_t ind,Int_t pade){
  BetaFn bf_ms(4,NFLAV,2); 
  BetaFn bf_boost(3+pade,NFLAV,2);
  bf_boost[2]=b2_boost_II(ind);
  if(pade){bf_boost[3]=bf_boost[2]*bf_boost[2]/bf_boost[1];}

  Double_t fac1=r0_over_a_extrap[ind];
  Double_t fac2=TMath::Exp(t1_boost(ind)/(2.*bf_ms[0]));
  Double_t fac3=mu_over_lambda(g_boost(ind),bf_ms)/(mu_over_lambda(g_ms_II(ind),bf_ms)*mu_over_lambda(g_boost(ind),bf_boost) );
  cout << "beta 0 , 1 , 2 , 3 :" << "  " << bf_boost[0] << "  " << bf_boost[1] << "  " << bf_boost[2] << "  " << bf_boost[3] << endl;
  return RATIO*fac1*fac2*fac3;
};
Double_t methodII_getr0Lambda_err(Int_t ind,Int_t pade){
  BetaFn bf_ms(4,NFLAV,2);
  BetaFn bf_boost(3,NFLAV,2);
  bf_boost[2]=b2_boost_II(ind);
  if(pade){bf_boost[3]=bf_boost[2]*bf_boost[2]/bf_boost[1];}

  Double_t partial_r0_over_a=methodII_getr0Lambda(ind,pade)/r0_over_a_extrap[ind];
  Double_t partial_t1_boost=methodII_getr0Lambda(ind,pade)/(2.*bf_ms[0]);
  Double_t partial_g_boost= -methodII_getr0Lambda(ind,pade)/bf_boost.eval(g_boost(ind)); 
  Double_t partial_b2_boost= -methodII_getr0Lambda(ind,pade)*partial_b2_boost_II_factor(ind,pade); 
  Double_t partial_b3_boost= -methodII_getr0Lambda(ind,pade)*partial_b3_boost_II_factor(ind,pade); 

  Double_t corr_r0_over_a__t1_boost=0.;
  Double_t corr_r0_over_a__g_boost=0.;
  Double_t corr_r0_over_a__b2_boost=0.;
  Double_t corr_r0_over_a__b3_boost=0.;
  Double_t corr_t1_boost__g_boost=1.;
  Double_t corr_t1_boost__b2_boost=1.;
  Double_t corr_t1_boost__b3_boost=1.;
  Double_t corr_g_boost__b2_boost=1.;
  Double_t corr_g_boost__b3_boost=1.;
  Double_t corr_b2_boost__b3_boost=1.;

  Double_t term1=TMath::Power(partial_r0_over_a,2.)*TMath::Power(r0_over_a_extrap_err[ind],2.);
  Double_t term2=TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.);
  Double_t term3=TMath::Power(partial_g_boost,2.)*TMath::Power(g_boost_err(ind),2.);
  Double_t term4=TMath::Power(partial_b2_boost,2.)*TMath::Power(b2_boost_II_err(ind),2.);
  Double_t term5=TMath::Power(partial_b3_boost,2.)*TMath::Power(b3_boost_II_err(ind),2.);

  Double_t termcorr1=2.*partial_r0_over_a*partial_t1_boost*r0_over_a_extrap_err[ind]*t1_boost_err(ind)*corr_r0_over_a__t1_boost;
  Double_t termcorr2=2.*partial_r0_over_a*partial_g_boost*r0_over_a_extrap_err[ind]*g_boost_err(ind)*corr_r0_over_a__g_boost;
  Double_t termcorr3=2.*partial_r0_over_a*partial_b2_boost*r0_over_a_extrap_err[ind]*b2_boost_II_err(ind)*corr_r0_over_a__b2_boost;
  Double_t termcorr4=2.*partial_r0_over_a*partial_b3_boost*r0_over_a_extrap_err[ind]*b3_boost_II_err(ind)*corr_r0_over_a__b3_boost;
  Double_t termcorr5=2.*partial_t1_boost*partial_g_boost*t1_boost_err(ind)*g_boost_err(ind)*corr_t1_boost__g_boost;
  Double_t termcorr6=2.*partial_t1_boost*partial_b2_boost*t1_boost_err(ind)*b2_boost_II_err(ind)*corr_t1_boost__b2_boost;
  Double_t termcorr7=2.*partial_t1_boost*partial_b3_boost*t1_boost_err(ind)*b3_boost_II_err(ind)*corr_t1_boost__b3_boost;
  Double_t termcorr8=2.*partial_g_boost*partial_b2_boost*g_boost_err(ind)*b2_boost_II_err(ind)*corr_g_boost__b2_boost;
  Double_t termcorr9=2.*partial_g_boost*partial_b3_boost*g_boost_err(ind)*b3_boost_II_err(ind)*corr_g_boost__b3_boost;
  Double_t termcorr10=2.*partial_b2_boost*partial_b3_boost*b2_boost_II_err(ind)*b3_boost_II_err(ind)*corr_b2_boost__b3_boost;

  return TMath::Sqrt(term1+term2+term3+term4+term5
                      +termcorr1+termcorr2+termcorr3+termcorr4+termcorr5+termcorr6+termcorr7+termcorr8+termcorr9+termcorr10);
};

// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_methodII_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=0;i<NBETA;++i) { //skip the first two due to order(a^2) effects and skip the last
    num=TMath::Power( methodII_r0Lambda[i] - fitline(a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodII_r0Lambda_err[i],2.)+TMath::Power(a_over_r0_sqr_err[i]*p[1],2.);
    chi2 += num/den;
  }
  fval = chi2;
}
// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_methodIIP_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=0;i<NBETA;++i) { //skip the first two due to order(a^2) effects and skip the last
    num=TMath::Power( methodIIP_r0Lambda[i] - fitline(a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodIIP_r0Lambda_err[i],2.)+TMath::Power(a_over_r0_sqr_err[i]*p[1],2.);
    chi2 += num/den;
  }
  fval = chi2;
}

// ----------------------- method III --------------------------

Double_t g_ms_III(Int_t ind){
  BetaFn bf(4,NFLAV,2);
  //not using the simple "return g_boost" because we want to look at the influence of RATIO on the outcome of r0Lambda
  return TMath::Sqrt( TMath::Power( g_boost(ind) , 2. ) - (2.*bf[0]*TMath::Log(RATIO))*TMath::Power( g_boost(ind) , 4. )
                        + (4.*bf[0]*bf[0]*TMath::Log(RATIO)*TMath::Log(RATIO)-2.*bf[1]*TMath::Log(RATIO))*TMath::Power(g_boost(ind),6.) );
};
Double_t TF1_partial_b2_boost_III_factor_integrand(Double_t* x,Double_t* par){
// par[0]==ind, par[1]==pade?
  Int_t pade=(Int_t) par[1];
  BetaFn bf_boost(3+pade,NFLAV,2);
  bf_boost[2]=b2_boost_III(par[0]);
  if(pade==PADE){bf_boost[3]=b3_boost_III(par[0]);}
  return TMath::Power(x[0],7.) * TMath::Power( bf_boost.eval(x[0]), -2. );
}
Double_t partial_b2_boost_III_factor(Int_t ind,Int_t pade){
  TF1 scratch("scratch",TF1_partial_b2_boost_III_factor_integrand,0,100,2);
  scratch.SetParameter(0,ind);
  scratch.SetParameter(1,pade);
  return scratch.Integral(0.,g_boost(ind));
};
Double_t TF1_partial_b3_boost_III_factor_integrand(Double_t* x,Double_t* par){
// par[0]==ind
  BetaFn bf_boost(4,NFLAV,2);
  bf_boost[2]=b2_boost_III(par[0]);
  bf_boost[3]=bf_boost[2]*bf_boost[2]/bf_boost[1];
  return TMath::Power(x[0],9.) * TMath::Power( bf_boost.eval(x[0]), -2. );
}
Double_t partial_b3_boost_III_factor(Int_t ind,Int_t pade){
  if(pade!=1){return 0.;}
  TF1 scratch("scratch",TF1_partial_b3_boost_III_factor_integrand,0,100,1);
  scratch.SetParameter(0,ind);
  return scratch.Integral(0.,g_boost(ind));
};

Double_t methodIII_getr0Lambda(Int_t ind,Int_t pade){
  BetaFn bf_ms(4,NFLAV,2); 
  BetaFn bf_boost(3+pade,NFLAV,2);
  bf_boost[2]=b2_boost_III(ind);
  if(pade==PADE){bf_boost[3]=bf_boost[2]*bf_boost[2]/bf_boost[1];}

  Double_t fac1=r0_over_a_extrap[ind];
//  Double_t fac2=TMath::Exp(t1_boost(ind)/(2.*bf_ms[0])); // HERIN LIES THE ERROR
  Double_t fac2=TMath::Exp(t1_boost(-1,1.)/(2.*bf_ms[0])); // AND HERIN LIES THE FIX
  Double_t fac3=mu_over_lambda(g_boost(ind),bf_ms)/(mu_over_lambda(g_ms_III(ind),bf_ms)*mu_over_lambda(g_boost(ind),bf_boost) );
  cout << "beta 0 , 1 , 2 , 3 :" << "  " << bf_boost[0] << "  " << bf_boost[1] << "  " << bf_boost[2] << "  " << bf_boost[3] << "  " << fac1 << "  " << fac2 << "  " << fac3 << "  " << fac1*fac2*fac3<< endl;
  return RATIO*fac1*fac2*fac3;
};
Double_t methodIII_getr0Lambda_err(Int_t ind,Int_t pade){
  BetaFn bf_ms(4,NFLAV,2);
  BetaFn bf_boost(3+pade,NFLAV,2);
  bf_boost[2]=b2_boost_III(ind);
  if(pade){bf_boost[3]=bf_boost[2]*bf_boost[2]/bf_boost[1];}

  Double_t partial_r0_over_a=methodIII_getr0Lambda(ind,pade)/r0_over_a_extrap[ind];
  Double_t partial_t1_boost=methodIII_getr0Lambda(ind,pade)/(2.*bf_ms[0]);
  Double_t partial_g_boost= -methodIII_getr0Lambda(ind,pade)/bf_boost.eval(g_boost(ind)); 
  Double_t partial_b2_boost= -methodIII_getr0Lambda(ind,pade)*partial_b2_boost_III_factor(ind,pade); 
  Double_t partial_b3_boost= -methodIII_getr0Lambda(ind,pade)*partial_b3_boost_III_factor(ind,pade); 

  Double_t corr_r0_over_a__t1_boost=0.;  // I think for nf=2, they're correlated since r0/a and P both come from am_q... hmmm...
  Double_t corr_r0_over_a__g_boost=0.;
  Double_t corr_r0_over_a__b2_boost=0.;
  Double_t corr_r0_over_a__b3_boost=0.;
  Double_t corr_t1_boost__g_boost=1.;
  Double_t corr_t1_boost__b2_boost=1.;
  Double_t corr_t1_boost__b3_boost=1.;
  Double_t corr_g_boost__b2_boost=1.;
  Double_t corr_g_boost__b3_boost=1.;
  Double_t corr_b2_boost__b3_boost=1.;

  Double_t term1=TMath::Power(partial_r0_over_a,2.)*TMath::Power(r0_over_a_extrap_err[ind],2.);
  Double_t term2=0.*TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.);  // we fixed c_sw_boost to be 1, so there is no error on t1_boost
  Double_t term3=TMath::Power(partial_g_boost,2.)*TMath::Power(g_boost_err(ind),2.);
  Double_t term4=TMath::Power(partial_b2_boost,2.)*TMath::Power(b2_boost_III_err(ind),2.);
  Double_t term5=TMath::Power(partial_b3_boost,2.)*TMath::Power(b3_boost_III_err(ind),2.);

  Double_t termcorr1=2.*partial_r0_over_a*partial_t1_boost*r0_over_a_extrap_err[ind]*t1_boost_err(ind)*corr_r0_over_a__t1_boost;
  Double_t termcorr2=2.*partial_r0_over_a*partial_g_boost*r0_over_a_extrap_err[ind]*g_boost_err(ind)*corr_r0_over_a__g_boost;
  Double_t termcorr3=2.*partial_r0_over_a*partial_b2_boost*r0_over_a_extrap_err[ind]*b2_boost_III_err(ind)*corr_r0_over_a__b2_boost;
  Double_t termcorr4=2.*partial_r0_over_a*partial_b3_boost*r0_over_a_extrap_err[ind]*b3_boost_III_err(ind)*corr_r0_over_a__b3_boost;
  Double_t termcorr5=2.*partial_t1_boost*partial_g_boost*t1_boost_err(ind)*g_boost_err(ind)*corr_t1_boost__g_boost;
  Double_t termcorr6=2.*partial_t1_boost*partial_b2_boost*t1_boost_err(ind)*b2_boost_III_err(ind)*corr_t1_boost__b2_boost;
  Double_t termcorr7=2.*partial_t1_boost*partial_b3_boost*t1_boost_err(ind)*b3_boost_III_err(ind)*corr_t1_boost__b3_boost;
  Double_t termcorr8=2.*partial_g_boost*partial_b2_boost*g_boost_err(ind)*b2_boost_III_err(ind)*corr_g_boost__b2_boost;
  Double_t termcorr9=2.*partial_g_boost*partial_b3_boost*g_boost_err(ind)*b3_boost_III_err(ind)*corr_g_boost__b3_boost;
  Double_t termcorr10=2.*partial_b2_boost*partial_b3_boost*b2_boost_III_err(ind)*b3_boost_III_err(ind)*corr_b2_boost__b3_boost;

  return TMath::Sqrt(term1+term2+term3+term4+term5
                      +termcorr1+termcorr2+termcorr3+termcorr4+termcorr5+termcorr6+termcorr7+termcorr8+termcorr9+termcorr10);
};

// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_methodIII_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=0;i<NBETA;++i) { //skip the first two due to order(a^2) effects and skip the last
    num=TMath::Power( methodIII_r0Lambda[i] - fitline(a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodIII_r0Lambda_err[i],2.)+TMath::Power(a_over_r0_sqr_err[i]*p[1],2.);
    chi2 += num/den;
  }
  fval = chi2;
}
// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_methodIIIP_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=0;i<NBETA;++i) { //skip the first two due to order(a^2) effects and skip the last
    num=TMath::Power( methodIIIP_r0Lambda[i] - fitline(a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodIIIP_r0Lambda_err[i],2.)+TMath::Power(a_over_r0_sqr_err[i]*p[1],2.);
    chi2 += num/den;
  }
  fval = chi2;
}

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//                                                  MAIN BEGIN                                                    //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

Int_t main(){

  TVirtualFitter::SetDefaultFitter("Minuit");
  Double_t arglist[100];

  cout << "\n\n ----- BEGIN CHIRAL LIMIT EXTRAPOLATION ------ \n\n";

  cout << "          amq              |      plaquette       |       r0/a          \n";
  for(Int_t i=0;i<NDATA;++i){
    quarkmass[i]=amq(i);
    quarkmass_err[i]=amq_err(i);
    cout << quarkmass[i]<<"+/-"<<quarkmass_err[i]<<"       "
         << plaquette[i]<<"+/-"<<plaquette_err[i]<<"       "
         << r0_over_a[i]<<"+/-"<<r0_over_a_err[i]<<"       "
         <<endl;
  }

  // ---------------------------------------------------------------- //
  //                                                                  //
  //                     Chiral Limit Extrapolation                   //
  //                                                                  //
  // ---------------------------------------------------------------- //


  // ++++++++++++++++++++++++++++ PLAQUETTE +++++++++++++++++++++++++++++++++

  //setting up the canvas and frame for plaquette data
  TCanvas* chiralcanv = new TCanvas("chiralcanv","chiralcanv",700,700);
  TH1F* chiralframe=new TH1F("chiralframe","Average Plaquette vs ?Quark Mass?  ",1000,200,600);
  chiralframe->SetStats(0);
  chiralframe->SetMinimum(0.52);
  chiralframe->SetMaximum(0.57);
  chiralframe->GetXaxis()->SetTitle("am_{q}");
  chiralframe->GetXaxis()->SetTickLength(0.02);
  chiralframe->GetXaxis()->SetLabelSize(0.020);
  chiralframe->GetXaxis()->SetLimits(-0.005,0.08); //minx/maxx limits
  chiralframe->GetYaxis()->SetTitle("P");
  chiralframe->GetYaxis()->SetTickLength(0.02);
  chiralframe->GetYaxis()->SetLabelSize(0.020);
  chiralframe->Draw(" ");

  //plotting the data points
  TGraphErrors* chiralgraph=new TGraphErrors(15,quarkmass,plaquette,quarkmass_err,plaquette_err);
  chiralgraph->SetMarkerColor(1);
  chiralgraph->SetMarkerStyle(1);
  chiralgraph->SetMarkerSize(1);
  chiralgraph->Draw("P");

  //line at x=0
  Double_t linex_chiralplaq[2]={0.,0.};
  Double_t liney_chiralplaq[2]={0.,1.};
  TGraph* chiralplaqline=new TGraph(2,linex_chiralplaq,liney_chiralplaq);
  chiralplaqline->SetLineColor(1);
  chiralplaqline->SetLineStyle(3);
  chiralplaqline->Draw("L");

  //  -------------- fit with minuit for beta=5.20 ---------------
  // Minimize the statistics with Minuit and thus find the best estimator for the y intercept
  TVirtualFitter* minuit_chiral_plaq520 = TVirtualFitter::Fitter(0,3);
  minuit_chiral_plaq520->SetParameter(0,"d0",0.001, 0.1, 0., 1.);
  minuit_chiral_plaq520->SetParameter(1,"d1",0.0, 0.1,-5., 5.);
  minuit_chiral_plaq520->SetParameter(2,"d2",0.0, 0.1,-5., 5.);
  minuit_chiral_plaq520->SetFCN(minuit_chiral_plaq520_chisqr);

  // minimize
  minuit_chiral_plaq520->ExecuteCommand("MIGRAD",arglist,0);
  minuit_chiral_plaq520->ExecuteCommand("MINOS",arglist,0);

  // local variables
  Double_t pval_chiral_plaq520[3],perr_chiral_plaq520[3],plo_chiral_plaq520[3],phi_chiral_plaq520[3],pco_chiral_plaq520[3];
  
  // Fit results
  pval_chiral_plaq520[0] = minuit_chiral_plaq520->GetParameter(0);
  pval_chiral_plaq520[1] = minuit_chiral_plaq520->GetParameter(1);
  pval_chiral_plaq520[2] = minuit_chiral_plaq520->GetParameter(2);
  minuit_chiral_plaq520->GetErrors(0,phi_chiral_plaq520[0],plo_chiral_plaq520[0],perr_chiral_plaq520[0],pco_chiral_plaq520[0]);
  minuit_chiral_plaq520->GetErrors(1,phi_chiral_plaq520[1],plo_chiral_plaq520[1],perr_chiral_plaq520[1],pco_chiral_plaq520[1]);
  minuit_chiral_plaq520->GetErrors(2,phi_chiral_plaq520[2],plo_chiral_plaq520[2],perr_chiral_plaq520[2],pco_chiral_plaq520[2]);
  cout << endl << pval_chiral_plaq520[0] << " +/- " << perr_chiral_plaq520[0] << " + " << phi_chiral_plaq520[0] << " " << plo_chiral_plaq520[0] << " " << pco_chiral_plaq520[0] << endl;
  cout << pval_chiral_plaq520[1] << " +/- " << perr_chiral_plaq520[1] << " + " << phi_chiral_plaq520[1] << " " << plo_chiral_plaq520[1] << " " << pco_chiral_plaq520[1] << endl;
  cout << pval_chiral_plaq520[2] << " +/- " << perr_chiral_plaq520[2] << " + " << phi_chiral_plaq520[2] << " " << plo_chiral_plaq520[2] << " " << pco_chiral_plaq520[2] << endl << endl;

  TF1* method_chiral_plaq520_fit=new TF1("method_chiral_plaq520_fit",TF1_fit_parabola,0.,0.07,3);
  method_chiral_plaq520_fit->SetParameter(0,pval_chiral_plaq520[0]);
  method_chiral_plaq520_fit->SetParameter(1,pval_chiral_plaq520[1]);
  method_chiral_plaq520_fit->SetParameter(2,pval_chiral_plaq520[2]);
  method_chiral_plaq520_fit->SetLineWidth(1.5);
  method_chiral_plaq520_fit->SetLineColor(1);
  method_chiral_plaq520_fit->SetLineStyle(2);
  method_chiral_plaq520_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_plaq520[1]={0.};
  Double_t extrap_y_chiral_plaq520[1]={pval_chiral_plaq520[0]};
  Double_t extrap_x_err_chiral_plaq520[1]={0.};
  Double_t extrap_y_err_chiral_plaq520[1]={perr_chiral_plaq520[0]};
  TGraphErrors* extrap_chiral_plaq520=new TGraphErrors(1,extrap_x_chiral_plaq520,extrap_y_chiral_plaq520,extrap_x_err_chiral_plaq520,extrap_y_err_chiral_plaq520);
  extrap_chiral_plaq520->SetMarkerColor(1);
  extrap_chiral_plaq520->SetMarkerStyle(1);
  extrap_chiral_plaq520->SetMarkerSize(1);
  extrap_chiral_plaq520->Draw("P");

  // Result:
  plaquette_extrap[0]=extrap_y_chiral_plaq520[0];
  plaquette_extrap_err[0]=extrap_y_err_chiral_plaq520[0];

  //  -------------- fit with minuit for beta=5.25 ---------------
  // Minimize the statistics with Minuit and thus find the best estimator for the y intercept
  TVirtualFitter* minuit_chiral_plaq525 = TVirtualFitter::Fitter(0,3);
  minuit_chiral_plaq525->SetParameter(0,"d0",0.001, 0.1, 0., 1.);
  minuit_chiral_plaq525->SetParameter(1,"d1",0.0, 0.1,-5., 5.);
  minuit_chiral_plaq525->SetParameter(2,"d2",0.0, 0.1,-5., 5.);
  minuit_chiral_plaq525->SetFCN(minuit_chiral_plaq525_chisqr);

  // minimize
  minuit_chiral_plaq525->ExecuteCommand("MIGRAD",arglist,0);
  minuit_chiral_plaq525->ExecuteCommand("MINOS",arglist,0);

  // local variables
  Double_t pval_chiral_plaq525[3],perr_chiral_plaq525[3],plo_chiral_plaq525[3],phi_chiral_plaq525[3],pco_chiral_plaq525[3];
  
  // Fit results
  pval_chiral_plaq525[0] = minuit_chiral_plaq525->GetParameter(0);
  pval_chiral_plaq525[1] = minuit_chiral_plaq525->GetParameter(1);
  pval_chiral_plaq525[2] = minuit_chiral_plaq525->GetParameter(2);
  minuit_chiral_plaq525->GetErrors(0,phi_chiral_plaq525[0],plo_chiral_plaq525[0],perr_chiral_plaq525[0],pco_chiral_plaq525[0]);
  minuit_chiral_plaq525->GetErrors(1,phi_chiral_plaq525[1],plo_chiral_plaq525[1],perr_chiral_plaq525[1],pco_chiral_plaq525[1]);
  minuit_chiral_plaq525->GetErrors(2,phi_chiral_plaq525[2],plo_chiral_plaq525[2],perr_chiral_plaq525[2],pco_chiral_plaq525[2]);
  cout << endl << pval_chiral_plaq525[0] << " +/- " << perr_chiral_plaq525[0] << " + " << phi_chiral_plaq525[0] << " " << plo_chiral_plaq525[0] << " " << pco_chiral_plaq525[0] << endl;
  cout << pval_chiral_plaq525[1] << " +/- " << perr_chiral_plaq525[1] << " + " << phi_chiral_plaq525[1] << " " << plo_chiral_plaq525[1] << " " << pco_chiral_plaq525[1] << endl;
  cout << pval_chiral_plaq525[2] << " +/- " << perr_chiral_plaq525[2] << " + " << phi_chiral_plaq525[2] << " " << plo_chiral_plaq525[2] << " " << pco_chiral_plaq525[2] << endl << endl;

  TF1* method_chiral_plaq525_fit=new TF1("method_chiral_plaq525_fit",TF1_fit_parabola,0.,0.07,3);
  method_chiral_plaq525_fit->SetParameter(0,pval_chiral_plaq525[0]);
  method_chiral_plaq525_fit->SetParameter(1,pval_chiral_plaq525[1]);
  method_chiral_plaq525_fit->SetParameter(2,pval_chiral_plaq525[2]);
  method_chiral_plaq525_fit->SetLineWidth(1.5);
  method_chiral_plaq525_fit->SetLineColor(1);
  method_chiral_plaq525_fit->SetLineStyle(2);
  method_chiral_plaq525_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_plaq525[1]={0.};
  Double_t extrap_y_chiral_plaq525[1]={pval_chiral_plaq525[0]};
  Double_t extrap_x_err_chiral_plaq525[1]={0.};
  Double_t extrap_y_err_chiral_plaq525[1]={perr_chiral_plaq525[0]};
  TGraphErrors* extrap_chiral_plaq525=new TGraphErrors(1,extrap_x_chiral_plaq525,extrap_y_chiral_plaq525,extrap_x_err_chiral_plaq525,extrap_y_err_chiral_plaq525);
  extrap_chiral_plaq525->SetMarkerColor(1);
  extrap_chiral_plaq525->SetMarkerStyle(1);
  extrap_chiral_plaq525->SetMarkerSize(1);
  extrap_chiral_plaq525->Draw("P");

  // Result:
  plaquette_extrap[1]=extrap_y_chiral_plaq525[0];
  plaquette_extrap_err[1]=extrap_y_err_chiral_plaq525[0];

  //  -------------- fit with minuit for beta=5.29 ---------------
  // Minimize the statistics with Minuit and thus find the best estimator for the y intercept
  TVirtualFitter* minuit_chiral_plaq529 = TVirtualFitter::Fitter(0,3);
  minuit_chiral_plaq529->SetParameter(0,"d0",0.001, 0.1, 0., 1.);
  minuit_chiral_plaq529->SetParameter(1,"d1",0.0, 0.1,-5., 5.);
  minuit_chiral_plaq529->SetParameter(2,"d2",0.0, 0.1,-5., 5.);
  minuit_chiral_plaq529->SetFCN(minuit_chiral_plaq529_chisqr);

  // minimize
  minuit_chiral_plaq529->ExecuteCommand("MIGRAD",arglist,0);
  minuit_chiral_plaq529->ExecuteCommand("MINOS",arglist,0);

  // local variables
  Double_t pval_chiral_plaq529[3],perr_chiral_plaq529[3],plo_chiral_plaq529[3],phi_chiral_plaq529[3],pco_chiral_plaq529[3];
  
  // Fit results
  pval_chiral_plaq529[0] = minuit_chiral_plaq529->GetParameter(0);
  pval_chiral_plaq529[1] = minuit_chiral_plaq529->GetParameter(1);
  pval_chiral_plaq529[2] = minuit_chiral_plaq529->GetParameter(2);
  minuit_chiral_plaq529->GetErrors(0,phi_chiral_plaq529[0],plo_chiral_plaq529[0],perr_chiral_plaq529[0],pco_chiral_plaq529[0]);
  minuit_chiral_plaq529->GetErrors(1,phi_chiral_plaq529[1],plo_chiral_plaq529[1],perr_chiral_plaq529[1],pco_chiral_plaq529[1]);
  minuit_chiral_plaq529->GetErrors(2,phi_chiral_plaq529[2],plo_chiral_plaq529[2],perr_chiral_plaq529[2],pco_chiral_plaq529[2]);
  cout << endl << pval_chiral_plaq529[0] << " +/- " << perr_chiral_plaq529[0] << " + " << phi_chiral_plaq529[0] << " " << plo_chiral_plaq529[0] << " " << pco_chiral_plaq529[0] << endl;
  cout << pval_chiral_plaq529[1] << " +/- " << perr_chiral_plaq529[1] << " + " << phi_chiral_plaq529[1] << " " << plo_chiral_plaq529[1] << " " << pco_chiral_plaq529[1] << endl;
  cout << pval_chiral_plaq529[2] << " +/- " << perr_chiral_plaq529[2] << " + " << phi_chiral_plaq529[2] << " " << plo_chiral_plaq529[2] << " " << pco_chiral_plaq529[2] << endl << endl;

  TF1* method_chiral_plaq529_fit=new TF1("method_chiral_plaq529_fit",TF1_fit_parabola,0.,0.07,3);
  method_chiral_plaq529_fit->SetParameter(0,pval_chiral_plaq529[0]);
  method_chiral_plaq529_fit->SetParameter(1,pval_chiral_plaq529[1]);
  method_chiral_plaq529_fit->SetParameter(2,pval_chiral_plaq529[2]);
  method_chiral_plaq529_fit->SetLineWidth(1.5);
  method_chiral_plaq529_fit->SetLineColor(1);
  method_chiral_plaq529_fit->SetLineStyle(2);
  method_chiral_plaq529_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_plaq529[1]={0.};
  Double_t extrap_y_chiral_plaq529[1]={pval_chiral_plaq529[0]};
  Double_t extrap_x_err_chiral_plaq529[1]={0.};
  Double_t extrap_y_err_chiral_plaq529[1]={perr_chiral_plaq529[0]};
  TGraphErrors* extrap_chiral_plaq529=new TGraphErrors(1,extrap_x_chiral_plaq529,extrap_y_chiral_plaq529,extrap_x_err_chiral_plaq529,extrap_y_err_chiral_plaq529);
  extrap_chiral_plaq529->SetMarkerColor(1);
  extrap_chiral_plaq529->SetMarkerStyle(1);
  extrap_chiral_plaq529->SetMarkerSize(1);
  extrap_chiral_plaq529->Draw("P");

  // Result:
  plaquette_extrap[2]=extrap_y_chiral_plaq529[0];
  plaquette_extrap_err[2]=extrap_y_err_chiral_plaq529[0];

  //  -------------- fit with minuit for beta=5.40 ---------------
  // Minimize the statistics with Minuit and thus find the best estimator for the y intercept
  TVirtualFitter* minuit_chiral_plaq540 = TVirtualFitter::Fitter(0,3);
  minuit_chiral_plaq540->SetParameter(0,"d0",0.001, 0.1, 0., 1.);
  minuit_chiral_plaq540->SetParameter(1,"d1",0.0, 0.1,-5., 5.);
  minuit_chiral_plaq540->SetParameter(2,"d2",0.0, 0.1,-5., 5.);
  minuit_chiral_plaq540->SetFCN(minuit_chiral_plaq540_chisqr);

  // minimize
  minuit_chiral_plaq540->ExecuteCommand("MIGRAD",arglist,0);
  minuit_chiral_plaq540->ExecuteCommand("MINOS",arglist,0);

  // local variables
  Double_t pval_chiral_plaq540[3],perr_chiral_plaq540[3],plo_chiral_plaq540[3],phi_chiral_plaq540[3],pco_chiral_plaq540[3];
  
  // Fit results
  pval_chiral_plaq540[0] = minuit_chiral_plaq540->GetParameter(0);
  pval_chiral_plaq540[1] = minuit_chiral_plaq540->GetParameter(1);
  pval_chiral_plaq540[2] = minuit_chiral_plaq540->GetParameter(2);
  minuit_chiral_plaq540->GetErrors(0,phi_chiral_plaq540[0],plo_chiral_plaq540[0],perr_chiral_plaq540[0],pco_chiral_plaq540[0]);
  minuit_chiral_plaq540->GetErrors(1,phi_chiral_plaq540[1],plo_chiral_plaq540[1],perr_chiral_plaq540[1],pco_chiral_plaq540[1]);
  minuit_chiral_plaq540->GetErrors(2,phi_chiral_plaq540[2],plo_chiral_plaq540[2],perr_chiral_plaq540[2],pco_chiral_plaq540[2]);
  cout << endl << pval_chiral_plaq540[0] << " +/- " << perr_chiral_plaq540[0] << " + " << phi_chiral_plaq540[0] << " " << plo_chiral_plaq540[0] << " " << pco_chiral_plaq540[0] << endl;
  cout << pval_chiral_plaq540[1] << " +/- " << perr_chiral_plaq540[1] << " + " << phi_chiral_plaq540[1] << " " << plo_chiral_plaq540[1] << " " << pco_chiral_plaq540[1] << endl;
  cout << pval_chiral_plaq540[2] << " +/- " << perr_chiral_plaq540[2] << " + " << phi_chiral_plaq540[2] << " " << plo_chiral_plaq540[2] << " " << pco_chiral_plaq540[2] << endl << endl;

  TF1* method_chiral_plaq540_fit=new TF1("method_chiral_plaq540_fit",TF1_fit_parabola,0.,0.07,3);
  method_chiral_plaq540_fit->SetParameter(0,pval_chiral_plaq540[0]);
  method_chiral_plaq540_fit->SetParameter(1,pval_chiral_plaq540[1]);
  method_chiral_plaq540_fit->SetParameter(2,pval_chiral_plaq540[2]);
  method_chiral_plaq540_fit->SetLineWidth(1.5);
  method_chiral_plaq540_fit->SetLineColor(1);
  method_chiral_plaq540_fit->SetLineStyle(2);
  method_chiral_plaq540_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_plaq540[1]={0.};
  Double_t extrap_y_chiral_plaq540[1]={pval_chiral_plaq540[0]};
  Double_t extrap_x_err_chiral_plaq540[1]={0.};
  Double_t extrap_y_err_chiral_plaq540[1]={perr_chiral_plaq540[0]};
  TGraphErrors* extrap_chiral_plaq540=new TGraphErrors(1,extrap_x_chiral_plaq540,extrap_y_chiral_plaq540,extrap_x_err_chiral_plaq540,extrap_y_err_chiral_plaq540);
  extrap_chiral_plaq540->SetMarkerColor(1);
  extrap_chiral_plaq540->SetMarkerStyle(1);
  extrap_chiral_plaq540->SetMarkerSize(1);
  extrap_chiral_plaq540->Draw("P");

  // Result:
  plaquette_extrap[3]=extrap_y_chiral_plaq540[0];
  plaquette_extrap_err[3]=extrap_y_err_chiral_plaq540[0];

  chiralcanv->SaveAs("Plaquette_Chiral_Extrap.pdf");

  // +++++++++++++++++++++++++++++++++ R0 OVER A  +++++++++++++++++++++++++++++++++++

  // ------------ setting up the canvas and frame for plaquette data  ----------
  TCanvas* chiralr0canv = new TCanvas("chiralr0canv","chiralr0canv",700,700);
  TH1F* chiralr0frame=new TH1F("chiralr0frame","r_{0}/a vs ?Quark Mass?  ",1000,200,600);
  chiralr0frame->SetStats(0);
  chiralr0frame->SetMinimum(3.);
  chiralr0frame->SetMaximum(8.);
  chiralr0frame->GetXaxis()->SetTitle("am_{q}");
  chiralr0frame->GetXaxis()->SetTickLength(0.02);
  chiralr0frame->GetXaxis()->SetLabelSize(0.020);
  chiralr0frame->GetXaxis()->SetLimits(-0.005,0.08); //minx/maxx limits
  chiralr0frame->GetYaxis()->SetTitle("r_{0}/a");
  chiralr0frame->GetYaxis()->SetTickLength(0.02);
  chiralr0frame->GetYaxis()->SetLabelSize(0.020);
  chiralr0frame->Draw(" ");

  //plotting the data points
  TGraphErrors* chiralr0graph=new TGraphErrors(15,quarkmass,r0_over_a,quarkmass_err,r0_over_a_err);
  chiralr0graph->SetMarkerColor(1);
  chiralr0graph->SetMarkerStyle(1);
  chiralr0graph->SetMarkerSize(1);
  chiralr0graph->Draw("P");

  //line at x=0
  Double_t linex_chiralr0[2]={0.,0.};
  Double_t liney_chiralr0[2]={0.,9.};
  TGraph* chiralr0line=new TGraph(2,linex_chiralr0,liney_chiralr0);
  chiralr0line->SetLineColor(1);
  chiralr0line->SetLineStyle(3);
  chiralr0line->Draw("L");

  // Minimize the statistics with Minuit and thus find the best estimator for the y intercept
  TVirtualFitter* minuit_chiral_r0_global = TVirtualFitter::Fitter(0,8);
  minuit_chiral_r0_global->SetParameter(0,"A00",0., 0.1, -10., 10.);
  minuit_chiral_r0_global->SetParameter(1,"A01",0.0, 100,-10000., 10000.);
  minuit_chiral_r0_global->SetParameter(2,"A10",0.0, 100,-10000., 10000.);
  minuit_chiral_r0_global->SetParameter(3,"A11",0.0, 100,-10000., 10000.);
  minuit_chiral_r0_global->SetParameter(4,"A12",0.0, 100,-10000., 10000.);
  minuit_chiral_r0_global->SetParameter(5,"A20",0.0, 100,-10000., 10000.);
  minuit_chiral_r0_global->SetParameter(6,"A21",0.0, 100,-10000., 10000.);
  minuit_chiral_r0_global->SetParameter(7,"A22",0.0, 100,-10000., 10000.);
  minuit_chiral_r0_global->SetFCN(minuit_chiral_r0a_global_chisqr);

  // minimize
  minuit_chiral_r0_global->ExecuteCommand("MIGRAD",arglist,0);
  minuit_chiral_r0_global->ExecuteCommand("MINOS",arglist,0);

  // local variables
  Double_t pval_chiral_r0_global[8],perr_chiral_r0_global[8],plo_chiral_r0_global[8],phi_chiral_r0_global[8],pco_chiral_r0_global[8];
  
  // Fit results
  pval_chiral_r0_global[0] = minuit_chiral_r0_global->GetParameter(0);
  pval_chiral_r0_global[1] = minuit_chiral_r0_global->GetParameter(1);
  pval_chiral_r0_global[2] = minuit_chiral_r0_global->GetParameter(2);
  pval_chiral_r0_global[3] = minuit_chiral_r0_global->GetParameter(3);
  pval_chiral_r0_global[4] = minuit_chiral_r0_global->GetParameter(4);
  pval_chiral_r0_global[5] = minuit_chiral_r0_global->GetParameter(5);
  pval_chiral_r0_global[6] = minuit_chiral_r0_global->GetParameter(6);
  pval_chiral_r0_global[7] = minuit_chiral_r0_global->GetParameter(7);

  minuit_chiral_r0_global->GetErrors(0,phi_chiral_r0_global[0],plo_chiral_r0_global[0],perr_chiral_r0_global[0],pco_chiral_r0_global[0]);
  minuit_chiral_r0_global->GetErrors(1,phi_chiral_r0_global[1],plo_chiral_r0_global[1],perr_chiral_r0_global[1],pco_chiral_r0_global[1]);
  minuit_chiral_r0_global->GetErrors(2,phi_chiral_r0_global[2],plo_chiral_r0_global[2],perr_chiral_r0_global[2],pco_chiral_r0_global[2]);
  minuit_chiral_r0_global->GetErrors(3,phi_chiral_r0_global[3],plo_chiral_r0_global[3],perr_chiral_r0_global[3],pco_chiral_r0_global[3]);
  minuit_chiral_r0_global->GetErrors(4,phi_chiral_r0_global[4],plo_chiral_r0_global[4],perr_chiral_r0_global[4],pco_chiral_r0_global[4]);
  minuit_chiral_r0_global->GetErrors(5,phi_chiral_r0_global[5],plo_chiral_r0_global[5],perr_chiral_r0_global[5],pco_chiral_r0_global[5]);
  minuit_chiral_r0_global->GetErrors(6,phi_chiral_r0_global[6],plo_chiral_r0_global[6],perr_chiral_r0_global[6],pco_chiral_r0_global[6]);
  minuit_chiral_r0_global->GetErrors(7,phi_chiral_r0_global[7],plo_chiral_r0_global[7],perr_chiral_r0_global[7],pco_chiral_r0_global[7]);

  cout << endl << pval_chiral_r0_global[0] << " +/- " << perr_chiral_r0_global[0] << " + " << phi_chiral_r0_global[0] << " " << plo_chiral_r0_global[0] << " " << pco_chiral_r0_global[0] << endl;
  cout << pval_chiral_r0_global[1] << " +/- " << perr_chiral_r0_global[1] << " + " << phi_chiral_r0_global[1] << " " << plo_chiral_r0_global[1] << " " << pco_chiral_r0_global[1] << endl;
  cout << pval_chiral_r0_global[2] << " +/- " << perr_chiral_r0_global[2] << " + " << phi_chiral_r0_global[2] << " " << plo_chiral_r0_global[2] << " " << pco_chiral_r0_global[2] << endl;
  cout << pval_chiral_r0_global[3] << " +/- " << perr_chiral_r0_global[3] << " + " << phi_chiral_r0_global[3] << " " << plo_chiral_r0_global[3] << " " << pco_chiral_r0_global[3] << endl;
  cout << pval_chiral_r0_global[4] << " +/- " << perr_chiral_r0_global[4] << " + " << phi_chiral_r0_global[4] << " " << plo_chiral_r0_global[4] << " " << pco_chiral_r0_global[4] << endl;
  cout << pval_chiral_r0_global[5] << " +/- " << perr_chiral_r0_global[5] << " + " << phi_chiral_r0_global[5] << " " << plo_chiral_r0_global[5] << " " << pco_chiral_r0_global[5] << endl;
  cout << pval_chiral_r0_global[6] << " +/- " << perr_chiral_r0_global[6] << " + " << phi_chiral_r0_global[6] << " " << plo_chiral_r0_global[6] << " " << pco_chiral_r0_global[6] << endl;
  cout << pval_chiral_r0_global[7] << " +/- " << perr_chiral_r0_global[7] << " + " << phi_chiral_r0_global[7] << " " << plo_chiral_r0_global[7] << " " << pco_chiral_r0_global[7] << endl << endl;

  TMatrixDSym cov(8); 
  cov.SetMatrixArray(minuit_chiral_r0_global->GetCovarianceMatrix());
  cov.Print();

  cout << cov(0,0) << endl;

  // ------------ for beta=5.20 -------------
  TF1* method_chiral_r0_520_fit=new TF1("method_chiral_r0_520_fit",TF1_fit_global,0.,0.07,9);
  method_chiral_r0_520_fit->SetParameter(0,pval_chiral_r0_global[0]);
  method_chiral_r0_520_fit->SetParameter(1,pval_chiral_r0_global[1]);
  method_chiral_r0_520_fit->SetParameter(2,pval_chiral_r0_global[2]);
  method_chiral_r0_520_fit->SetParameter(3,pval_chiral_r0_global[3]);
  method_chiral_r0_520_fit->SetParameter(4,pval_chiral_r0_global[4]);
  method_chiral_r0_520_fit->SetParameter(5,pval_chiral_r0_global[5]);
  method_chiral_r0_520_fit->SetParameter(6,pval_chiral_r0_global[6]);
  method_chiral_r0_520_fit->SetParameter(7,pval_chiral_r0_global[7]);
  method_chiral_r0_520_fit->SetParameter(8,beta_extrap[0]);

  method_chiral_r0_520_fit->SetLineWidth(1.5);
  method_chiral_r0_520_fit->SetLineColor(1);
  method_chiral_r0_520_fit->SetLineStyle(2);
  method_chiral_r0_520_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_r0_520[1]={0.};
  Double_t extrap_y_chiral_r0_520[1]={fit_global(0,beta_extrap[0],pval_chiral_r0_global)};
  Double_t extrap_x_err_chiral_r0_520[1]={0.};
  Double_t extrap_y_err_chiral_r0_520[1]={fit_global(0,beta_extrap[0],pval_chiral_r0_global)*TMath::Sqrt(
      TMath::Power(perr_chiral_r0_global[0],2.)+TMath::Power(beta_extrap[0]*perr_chiral_r0_global[1],2.)+2*beta_extrap[0]*cov(0,1) )};
  TGraphErrors* extrap_chiral_r0_520=new TGraphErrors(1,extrap_x_chiral_r0_520,extrap_y_chiral_r0_520,extrap_x_err_chiral_r0_520,extrap_y_err_chiral_r0_520);
  extrap_chiral_r0_520->SetMarkerColor(1);
  extrap_chiral_r0_520->SetMarkerStyle(1);
  extrap_chiral_r0_520->SetMarkerSize(1);
  extrap_chiral_r0_520->Draw("P");

  // Result:
  r0_over_a_extrap[0]=extrap_y_chiral_r0_520[0];
  r0_over_a_extrap_err[0]=extrap_y_err_chiral_r0_520[0];

  // ------------ for beta=5.25 -------------
  TF1* method_chiral_r0_525_fit=new TF1("method_chiral_r0_525_fit",TF1_fit_global,0.,0.07,9);
  method_chiral_r0_525_fit->SetParameter(0,pval_chiral_r0_global[0]);
  method_chiral_r0_525_fit->SetParameter(1,pval_chiral_r0_global[1]);
  method_chiral_r0_525_fit->SetParameter(2,pval_chiral_r0_global[2]);
  method_chiral_r0_525_fit->SetParameter(3,pval_chiral_r0_global[3]);
  method_chiral_r0_525_fit->SetParameter(4,pval_chiral_r0_global[4]);
  method_chiral_r0_525_fit->SetParameter(5,pval_chiral_r0_global[5]);
  method_chiral_r0_525_fit->SetParameter(6,pval_chiral_r0_global[6]);
  method_chiral_r0_525_fit->SetParameter(7,pval_chiral_r0_global[7]);
  method_chiral_r0_525_fit->SetParameter(8,beta_extrap[1]);

  method_chiral_r0_525_fit->SetLineWidth(1.5);
  method_chiral_r0_525_fit->SetLineColor(1);
  method_chiral_r0_525_fit->SetLineStyle(2);
  method_chiral_r0_525_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_r0_525[1]={0.};
  Double_t extrap_y_chiral_r0_525[1]={fit_global(0,beta_extrap[1],pval_chiral_r0_global)};
  Double_t extrap_x_err_chiral_r0_525[1]={0.};
  Double_t extrap_y_err_chiral_r0_525[1]={fit_global(0,beta_extrap[1],pval_chiral_r0_global)*TMath::Sqrt(
     TMath::Power(perr_chiral_r0_global[0],2.)+TMath::Power(beta_extrap[1]*perr_chiral_r0_global[1],2.)+2*beta_extrap[1]*cov(0,1) )};
  TGraphErrors* extrap_chiral_r0_525=new TGraphErrors(1,extrap_x_chiral_r0_525,extrap_y_chiral_r0_525,extrap_x_err_chiral_r0_525,extrap_y_err_chiral_r0_525);
  extrap_chiral_r0_525->SetMarkerColor(1);
  extrap_chiral_r0_525->SetMarkerStyle(1);
  extrap_chiral_r0_525->SetMarkerSize(1);
  extrap_chiral_r0_525->Draw("P");

  // Result:
  r0_over_a_extrap[1]=extrap_y_chiral_r0_525[0];
  r0_over_a_extrap_err[1]=extrap_y_err_chiral_r0_525[0];

  // ------------ for beta=5.29 -------------
  TF1* method_chiral_r0_529_fit=new TF1("method_chiral_r0_529_fit",TF1_fit_global,0.,0.07,9);
  method_chiral_r0_529_fit->SetParameter(0,pval_chiral_r0_global[0]);
  method_chiral_r0_529_fit->SetParameter(1,pval_chiral_r0_global[1]);
  method_chiral_r0_529_fit->SetParameter(2,pval_chiral_r0_global[2]);
  method_chiral_r0_529_fit->SetParameter(3,pval_chiral_r0_global[3]);
  method_chiral_r0_529_fit->SetParameter(4,pval_chiral_r0_global[4]);
  method_chiral_r0_529_fit->SetParameter(5,pval_chiral_r0_global[5]);
  method_chiral_r0_529_fit->SetParameter(6,pval_chiral_r0_global[6]);
  method_chiral_r0_529_fit->SetParameter(7,pval_chiral_r0_global[7]);
  method_chiral_r0_529_fit->SetParameter(8,beta_extrap[2]);

  method_chiral_r0_529_fit->SetLineWidth(1.5);
  method_chiral_r0_529_fit->SetLineColor(1);
  method_chiral_r0_529_fit->SetLineStyle(2);
  method_chiral_r0_529_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_r0_529[1]={0.};
  Double_t extrap_y_chiral_r0_529[1]={fit_global(0,beta_extrap[2],pval_chiral_r0_global)};
  Double_t extrap_x_err_chiral_r0_529[1]={0.};
  Double_t extrap_y_err_chiral_r0_529[1]={fit_global(0,beta_extrap[2],pval_chiral_r0_global)*TMath::Sqrt(
      TMath::Power(perr_chiral_r0_global[0],2.)+TMath::Power(beta_extrap[2]*perr_chiral_r0_global[1],2.)+2*beta_extrap[2]*cov(0,1) )};
  TGraphErrors* extrap_chiral_r0_529=new TGraphErrors(1,extrap_x_chiral_r0_529,extrap_y_chiral_r0_529,extrap_x_err_chiral_r0_529,extrap_y_err_chiral_r0_529);
  extrap_chiral_r0_529->SetMarkerColor(1);
  extrap_chiral_r0_529->SetMarkerStyle(1);
  extrap_chiral_r0_529->SetMarkerSize(1);
  extrap_chiral_r0_529->Draw("P");

  // Result:
  r0_over_a_extrap[2]=extrap_y_chiral_r0_529[0];
  r0_over_a_extrap_err[2]=extrap_y_err_chiral_r0_529[0];

  // ------------ for beta=5.40 -------------
  TF1* method_chiral_r0_540_fit=new TF1("method_chiral_r0_540_fit",TF1_fit_global,0.,0.07,9);
  method_chiral_r0_540_fit->SetParameter(0,pval_chiral_r0_global[0]);
  method_chiral_r0_540_fit->SetParameter(1,pval_chiral_r0_global[1]);
  method_chiral_r0_540_fit->SetParameter(2,pval_chiral_r0_global[2]);
  method_chiral_r0_540_fit->SetParameter(3,pval_chiral_r0_global[3]);
  method_chiral_r0_540_fit->SetParameter(4,pval_chiral_r0_global[4]);
  method_chiral_r0_540_fit->SetParameter(5,pval_chiral_r0_global[5]);
  method_chiral_r0_540_fit->SetParameter(6,pval_chiral_r0_global[6]);
  method_chiral_r0_540_fit->SetParameter(7,pval_chiral_r0_global[7]);
  method_chiral_r0_540_fit->SetParameter(8,beta_extrap[3]);

  method_chiral_r0_540_fit->SetLineWidth(1.5);
  method_chiral_r0_540_fit->SetLineColor(1);
  method_chiral_r0_540_fit->SetLineStyle(2);
  method_chiral_r0_540_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_r0_540[1]={0.};
  Double_t extrap_y_chiral_r0_540[1]={fit_global(0,beta_extrap[3],pval_chiral_r0_global)};
  Double_t extrap_x_err_chiral_r0_540[1]={0.};
  Double_t extrap_y_err_chiral_r0_540[1]={fit_global(0,beta_extrap[3],pval_chiral_r0_global)*TMath::Sqrt(
      TMath::Power(perr_chiral_r0_global[0],2.)+TMath::Power(beta_extrap[3]*perr_chiral_r0_global[1],2.)+2*beta_extrap[3]*cov(0,1) )};
  TGraphErrors* extrap_chiral_r0_540=new TGraphErrors(1,extrap_x_chiral_r0_540,extrap_y_chiral_r0_540,extrap_x_err_chiral_r0_540,extrap_y_err_chiral_r0_540);
  extrap_chiral_r0_540->SetMarkerColor(1);
  extrap_chiral_r0_540->SetMarkerStyle(1);
  extrap_chiral_r0_540->SetMarkerSize(1);
  extrap_chiral_r0_540->Draw("P");

  // Result:
  r0_over_a_extrap[3]=extrap_y_chiral_r0_540[0];
  r0_over_a_extrap_err[3]=extrap_y_err_chiral_r0_540[0];

  chiralr0canv->SaveAs("r0_Over_a_Chiral_Extrap.pdf");

  //  --------- looking at plaquette vs r0_over_a -----------------

  // ------------ setting up the canvas and frame for plaquette data  ----------
  TCanvas* testcanv = new TCanvas("testcanv","testcanv",700,700);
  TH1F* testframe=new TH1F("testframe","plaquette vs r_{0}/a",1000,200,600);
  testframe->SetStats(0);
  testframe->SetMinimum(0.52);
  testframe->SetMaximum(0.57);
  testframe->GetXaxis()->SetTitle("r_{0}/a");
  testframe->GetXaxis()->SetTickLength(0.02);
  testframe->GetXaxis()->SetLabelSize(0.020);
  testframe->GetXaxis()->SetLimits(5.,8.); //minx/maxx limits
  testframe->GetYaxis()->SetTitle("P");
  testframe->GetYaxis()->SetTickLength(0.02);
  testframe->GetYaxis()->SetLabelSize(0.020);
  testframe->Draw(" ");

  //plotting the data points
  TGraphErrors* testgraph=new TGraphErrors(4,r0_over_a_extrap,plaquette_extrap,r0_over_a_err,plaquette_err);
  testgraph->SetMarkerColor(1);
  testgraph->SetMarkerStyle(1);
  testgraph->SetMarkerSize(1);
  testgraph->Fit("pol1");
  testgraph->Draw("P");

  testcanv->SaveAs("Plaq_vs_r0a.pdf");

  for(Int_t i=0;i<NBETA;++i){
    a_over_r0_sqr[i]=TMath::Power(r0_over_a_extrap[i],-2.);
    a_over_r0_sqr_err[i]=2.*r0_over_a_extrap_err[i]/(r0_over_a_extrap[i]*r0_over_a_extrap[i]*r0_over_a_extrap[i]);
  }

  cout << "\n\n ----- BEGIN METHOD I ------ \n\n";

  for(Int_t i=0;i<NBETA;++i){
    methodI_r0Lambda[i]=methodI_getr0Lambda(i,DIRECT);
    methodI_r0Lambda_err[i]=methodI_getr0Lambda_err(i,DIRECT);
    cout<<a_over_r0_sqr[i]<<"+/-"<<a_over_r0_sqr_err[i]<<"    "<<methodI_r0Lambda[i]<<"+/-"<<methodI_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;
/*
  TCanvas* methodIcanv = new TCanvas("methodIcanv","methodIcanv",700,700);
  
  TH1F* methodIframe=new TH1F("methodIframe","Method I with N_{f} = 2",1000,200,600);
  methodIframe->SetStats(0);
  methodIframe->SetMinimum(0.45);
  methodIframe->SetMaximum(0.70);
  methodIframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIframe->GetXaxis()->SetTickLength(0.02);
  methodIframe->GetXaxis()->SetLabelSize(0.020);
  methodIframe->GetXaxis()->SetLimits(-0.002,0.05); //minx/maxx limits
  methodIframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIframe->GetYaxis()->SetTickLength(0.02);
  methodIframe->GetYaxis()->SetLabelSize(0.020);
  methodIframe->Draw(" ");

  TGraphErrors* methodIgraph=new TGraphErrors(4,a_over_r0_sqr,methodI_r0Lambda,a_over_r0_sqr_err,methodI_r0Lambda_err);
  methodIgraph->SetMarkerColor(1);
  methodIgraph->SetMarkerStyle(33);
  methodIgraph->SetMarkerSize(1);
  methodIgraph->Draw("P");

  Double_t linexI[2]={0.,0.};
  Double_t lineyI[2]={0.,1.};
  TGraph* methodIline=new TGraph(2,linexI,lineyI);
  methodIline->SetLineColor(1);
  methodIline->SetLineStyle(3);
  methodIline->Draw("L");

  delete minuit_chiral_r0_global ;

  // Minimize the statistics with Minuit and thus find the best estimator for x_0
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuitI = TVirtualFitter::Fitter(0,2);
  minuitI->SetParameter(0,"a",0.0, 0.1, 0., 5.);
  minuitI->SetParameter(1,"b",0.0, 0.1,-5., 5.);
  minuitI->SetFCN(minuit_methodI_chisqr);

  // minimize
  Double_t arglistI[100];
  minuitI->ExecuteCommand("MIGRAD",arglistI,0);
  minuitI->ExecuteCommand("MINOS",arglistI,0);

  // local variables
  Double_t pvalI[2],perrI[2],ploI[2],phiI[2],pcoI[2];
  
  // Fit results
  pvalI[0] = minuitI->GetParameter(0);
  pvalI[1] = minuitI->GetParameter(1);
  minuitI->GetErrors(0,phiI[0],ploI[0],perrI[0],pcoI[0]);
  minuitI->GetErrors(1,phiI[1],ploI[1],perrI[1],pcoI[1]);
  cout << endl << pvalI[0] << " +/- " << perrI[0] << " + " << phiI[0] << " " << ploI[0] << " " << pcoI[0] << endl;
  cout << pvalI[1] << " +/- " << perrI[1] << " + " << phiI[1] << " " << ploI[1] << " " << pcoI[1] << endl << endl;

  TF1* methodI_fit=new TF1("methodI_fit",TF1_fitline,0.,0.045,2);
  methodI_fit->SetParameter(0,pvalI[0]);
  methodI_fit->SetParameter(1,pvalI[1]);
  methodI_fit->SetLineWidth(1.5);
  methodI_fit->SetLineColor(1);
  methodI_fit->SetLineStyle(2);
  methodI_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_xI[1]={0.};
  Double_t extrap_yI[1]={pvalI[0]};
  Double_t extrap_x_errI[1]={0.};
  Double_t extrap_y_errI[1]={perrI[0]};
  TGraphErrors* extrapI=new TGraphErrors(1,extrap_xI,extrap_yI,extrap_x_errI,extrap_y_errI);
  extrapI->SetMarkerColor(1);
  extrapI->SetMarkerStyle(27);
  extrapI->SetMarkerSize(1);
  extrapI->Draw("P");

  methodIcanv->SaveAs("methodI_nf=2.pdf");
*/
  cout << "\n\n ----- BEGIN METHOD II ------ \n\n";

  for(Int_t i=0;i<NBETA;++i){
    methodII_r0Lambda[i]=methodII_getr0Lambda(i,NOTPADE);
    methodII_r0Lambda_err[i]=methodII_getr0Lambda_err(i,NOTPADE);
    cout<<a_over_r0_sqr[i]<<"+/-"<<a_over_r0_sqr_err[i]<<"    "<<methodII_r0Lambda[i]<<"+/-"<<methodII_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;
/*
  TCanvas* methodIIcanv = new TCanvas("methodIIcanv","methodIIcanv",700,700);
  
  TH1F* methodIIframe=new TH1F("methodIIframe","Method II with N_{f} = 2",1000,200,600);
  methodIIframe->SetStats(0);
  methodIIframe->SetMinimum(0.45);
  methodIIframe->SetMaximum(0.70);
  methodIIframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIIframe->GetXaxis()->SetTickLength(0.02);
  methodIIframe->GetXaxis()->SetLabelSize(0.020);
  methodIIframe->GetXaxis()->SetLimits(-0.002,0.05); //minx/maxx limits
  methodIIframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIIframe->GetYaxis()->SetTickLength(0.02);
  methodIIframe->GetYaxis()->SetLabelSize(0.020);
  methodIIframe->Draw(" ");

  TGraphErrors* methodIIgraph=new TGraphErrors(12,a_over_r0_sqr,methodII_r0Lambda,a_over_r0_sqr_err,methodII_r0Lambda_err);
  methodIIgraph->SetMarkerColor(1);
  methodIIgraph->SetMarkerStyle(33);
  methodIIgraph->SetMarkerSize(1);
  methodIIgraph->Draw("P");

  Double_t linexII[2]={0.,0.};
  Double_t lineyII[2]={0.,1.};
  TGraph* methodIIline=new TGraph(2,linexII,lineyII);
  methodIIline->SetLineColor(1);
  methodIIline->SetLineStyle(3);
  methodIIline->Draw("L");

  // Minimize the statistics with Minuit and thus find the best estimator for x_0
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuitII = TVirtualFitter::Fitter(0,2);
  minuitII->SetParameter(0,"a",0.0, 0.1, 0., 5.);
  minuitII->SetParameter(1,"b",0.0, 0.1,-5., 5.);
  minuitII->SetFCN(minuit_methodII_chisqr);

  // minimize
  Double_t arglistII[100];
  minuitII->ExecuteCommand("MIGRAD",arglistII,0);
  minuitII->ExecuteCommand("MINOS",arglistII,0);

  // local variables
  Double_t pvalII[2],perrII[2],ploII[2],phiII[2],pcoII[2];
  
  // Fit results
  pvalII[0] = minuitII->GetParameter(0);
  pvalII[1] = minuitII->GetParameter(1);
  minuitII->GetErrors(0,phiII[0],ploII[0],perrII[0],pcoII[0]);
  minuitII->GetErrors(1,phiII[1],ploII[1],perrII[1],pcoII[1]);
  cout << endl << pvalII[0] << " +/- " << perrII[0] << " + " << phiII[0] << " " << ploII[0] << " " << pcoII[0] << endl;
  cout << pvalII[1] << " +/- " << perrII[1] << " + " << phiII[1] << " " << ploII[1] << " " << pcoII[1] << endl << endl;

  TF1* methodII_fit=new TF1("methodI_fit",TF1_fitline,0.,0.045,2);
  methodII_fit->SetParameter(0,pvalII[0]);
  methodII_fit->SetParameter(1,pvalII[1]);
  methodII_fit->SetLineWidth(1.5);
  methodII_fit->SetLineColor(1);
  methodII_fit->SetLineStyle(2);
  methodII_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_xII[1]={0.};
  Double_t extrap_yII[1]={pvalII[0]};
  Double_t extrap_x_errII[1]={0.};
  Double_t extrap_y_errII[1]={perrII[0]};
  TGraphErrors* extrapII=new TGraphErrors(1,extrap_xII,extrap_yII,extrap_x_errII,extrap_y_errII);
  extrapII->SetMarkerColor(1);
  extrapII->SetMarkerStyle(27);
  extrapII->SetMarkerSize(1);
  extrapII->Draw("P");

  methodIIcanv->SaveAs("methodII_nf=2.pdf");
*/
  cout << "\n\n ----- BEGIN METHOD IIP ------ \n\n";

  for(Int_t i=0;i<NBETA;++i){
    methodIIP_r0Lambda[i]=methodII_getr0Lambda(i,PADE);
    methodIIP_r0Lambda_err[i]=methodII_getr0Lambda_err(i,PADE);
    cout<<a_over_r0_sqr[i]<<"+/-"<<a_over_r0_sqr_err[i]<<"    "<<methodIIP_r0Lambda[i]<<"+/-"<<methodIIP_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;
/*
  TCanvas* methodIIPcanv = new TCanvas("methodIIPcanv","methodIIPcanv",700,700);
  
  TH1F* methodIIPframe=new TH1F("methodIIPframe","Method IIP with N_{f} = 2",1000,200,600);
  methodIIPframe->SetStats(0);
  methodIIPframe->SetMinimum(0.45);
  methodIIPframe->SetMaximum(0.70);
  methodIIPframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIIPframe->GetXaxis()->SetTickLength(0.02);
  methodIIPframe->GetXaxis()->SetLabelSize(0.020);
  methodIIPframe->GetXaxis()->SetLimits(-0.002,0.05); //minx/maxx limits
  methodIIPframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIIPframe->GetYaxis()->SetTickLength(0.02);
  methodIIPframe->GetYaxis()->SetLabelSize(0.020);
  methodIIPframe->Draw(" ");

  TGraphErrors* methodIIPgraph=new TGraphErrors(12,a_over_r0_sqr,methodIIP_r0Lambda,a_over_r0_sqr_err,methodIIP_r0Lambda_err);
  methodIIPgraph->SetMarkerColor(1);
  methodIIPgraph->SetMarkerStyle(33);
  methodIIPgraph->SetMarkerSize(1);
  methodIIPgraph->Draw("P");

  Double_t linexIIP[2]={0.,0.};
  Double_t lineyIIP[2]={0.,1.};
  TGraph* methodIIPline=new TGraph(2,linexIIP,lineyIIP);
  methodIIPline->SetLineColor(1);
  methodIIPline->SetLineStyle(3);
  methodIIPline->Draw("L");

  // Minimize the statistics with Minuit and thus find the best estimator for x_0
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuitIIP = TVirtualFitter::Fitter(0,2);
  minuitIIP->SetParameter(0,"a",0.0, 0.1, 0., 5.);
  minuitIIP->SetParameter(1,"b",0.0, 0.1,-5., 5.);
  minuitIIP->SetFCN(minuit_methodIIP_chisqr);

  // minimize
  Double_t arglistIIP[100];
  minuitIIP->ExecuteCommand("MIGRAD",arglistIIP,0);
  minuitIIP->ExecuteCommand("MINOS",arglistIIP,0);

  // local variables
  Double_t pvalIIP[2],perrIIP[2],ploIIP[2],phIIPI[2],pcoIIP[2];
  
  // Fit results
  pvalIIP[0] = minuitIIP->GetParameter(0);
  pvalIIP[1] = minuitIIP->GetParameter(1);
  minuitIIP->GetErrors(0,phIIPI[0],ploIIP[0],perrIIP[0],pcoIIP[0]);
  minuitIIP->GetErrors(1,phIIPI[1],ploIIP[1],perrIIP[1],pcoIIP[1]);
  cout << endl << pvalIIP[0] << " +/- " << perrIIP[0] << " + " << phIIPI[0] << " " << ploIIP[0] << " " << pcoIIP[0] << endl;
  cout << pvalIIP[1] << " +/- " << perrIIP[1] << " + " << phIIPI[1] << " " << ploIIP[1] << " " << pcoIIP[1] << endl << endl;

  TF1* methodIIP_fit=new TF1("methodIIP_fit",TF1_fitline,0.,0.045,2);
  methodIIP_fit->SetParameter(0,pvalIIP[0]);
  methodIIP_fit->SetParameter(1,pvalIIP[1]);
  methodIIP_fit->SetLineWidth(1.5);
  methodIIP_fit->SetLineColor(1);
  methodIIP_fit->SetLineStyle(2);
  methodIIP_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_xIIP[1]={0.};
  Double_t extrap_yIIP[1]={pvalIIP[0]};
  Double_t extrap_x_errIIP[1]={0.};
  Double_t extrap_y_errIIP[1]={perrIIP[0]};
  TGraphErrors* extrapIIP=new TGraphErrors(1,extrap_xIIP,extrap_yIIP,extrap_x_errIIP,extrap_y_errIIP);
  extrapIIP->SetMarkerColor(1);
  extrapIIP->SetMarkerStyle(27);
  extrapIIP->SetMarkerSize(1);
  extrapIIP->Draw("P");

  methodIIPcanv->SaveAs("methodIIP_nf=2.pdf");
*/
  cout << "\n\n ----- BEGIN METHOD III ------ \n\n";

  for(Int_t i=0;i<NBETA;++i){
    methodIII_r0Lambda[i]=methodIII_getr0Lambda(i,NOTPADE);
    methodIII_r0Lambda_err[i]=methodIII_getr0Lambda_err(i,NOTPADE);
    cout<<a_over_r0_sqr[i]<<"+/-"<<a_over_r0_sqr_err[i]<<"    "<<methodIII_r0Lambda[i]<<"+/-"<<methodIII_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;


  TCanvas* methodIIIcanv = new TCanvas("methodIIIcanv","methodIIIcanv",700,700);
  
  TH1F* methodIIIframe=new TH1F("methodIIIframe","Method III with N_{f} = 2",1000,200,600);
  methodIIIframe->SetStats(0);
  methodIIIframe->SetMinimum(0.45);
  methodIIIframe->SetMaximum(0.70);
  methodIIIframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIIIframe->GetXaxis()->SetTickLength(0.02);
  methodIIIframe->GetXaxis()->SetLabelSize(0.020);
  methodIIIframe->GetXaxis()->SetLimits(-0.002,0.05); //minx/maxx limits
  methodIIIframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIIIframe->GetYaxis()->SetTickLength(0.02);
  methodIIIframe->GetYaxis()->SetLabelSize(0.020);
  methodIIIframe->Draw(" ");

  TGraphErrors* methodIIIgraph=new TGraphErrors(12,a_over_r0_sqr,methodIII_r0Lambda,a_over_r0_sqr_err,methodIII_r0Lambda_err);
  methodIIIgraph->SetMarkerColor(1);
  methodIIIgraph->SetMarkerStyle(33);
  methodIIIgraph->SetMarkerSize(1);
  methodIIIgraph->Draw("P");

  Double_t linexIII[2]={0.,0.};
  Double_t lineyIII[2]={0.,1.};
  TGraph* methodIIIline=new TGraph(2,linexIII,lineyIII);
  methodIIIline->SetLineColor(1);
  methodIIIline->SetLineStyle(3);
  methodIIIline->Draw("L");

  // Minimize the statistics with Minuit and thus find the best estimator for x_0
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuitIII = TVirtualFitter::Fitter(0,2);
  minuitIII->SetParameter(0,"a",0.0, 0.1, 0., 5.);
  minuitIII->SetParameter(1,"b",0.0, 0.1,-5., 5.);
  minuitIII->SetFCN(minuit_methodIII_chisqr);

  // minimize
  Double_t arglistIII[100];
  minuitIII->ExecuteCommand("MIGRAD",arglistIII,0);
  minuitIII->ExecuteCommand("MINOS",arglistIII,0);

  // local variables
  Double_t pvalIII[2],perrIII[2],ploIII[2],phiIII[2],pcoIII[2];
  
  // Fit results
  pvalIII[0] = minuitIII->GetParameter(0);
  pvalIII[1] = minuitIII->GetParameter(1);
  minuitIII->GetErrors(0,phiIII[0],ploIII[0],perrIII[0],pcoIII[0]);
  minuitIII->GetErrors(1,phiIII[1],ploIII[1],perrIII[1],pcoIII[1]);
  cout << endl << pvalIII[0] << " +/- " << perrIII[0] << " + " << phiIII[0] << " " << ploIII[0] << " " << pcoIII[0] << endl;
  cout << pvalIII[1] << " +/- " << perrIII[1] << " + " << phiIII[1] << " " << ploIII[1] << " " << pcoIII[1] << endl << endl;

  TF1* methodIII_fit=new TF1("methodIII_fit",TF1_fitline,0.,0.045,2);
  methodIII_fit->SetParameter(0,pvalIII[0]);
  methodIII_fit->SetParameter(1,pvalIII[1]);
  methodIII_fit->SetLineWidth(1.5);
  methodIII_fit->SetLineColor(1);
  methodIII_fit->SetLineStyle(2);
  methodIII_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_xIII[1]={0.};
  Double_t extrap_yIII[1]={pvalIII[0]};
  Double_t extrap_x_errIII[1]={0.};
  Double_t extrap_y_errIII[1]={perrIII[0]};
  TGraphErrors* extrapIII=new TGraphErrors(1,extrap_xIII,extrap_yIII,extrap_x_errIII,extrap_y_errIII);
  extrapIII->SetMarkerColor(1);
  extrapIII->SetMarkerStyle(27);
  extrapIII->SetMarkerSize(1);
  extrapIII->Draw("P");

  methodIIIcanv->SaveAs("methodIII_nf=2.pdf");

  cout << "\n\n ----- BEGIN METHOD IIIP ------ \n\n";

  for(Int_t i=0;i<NBETA;++i){
    methodIIIP_r0Lambda[i]=methodIII_getr0Lambda(i,PADE);
    methodIIIP_r0Lambda_err[i]=methodIII_getr0Lambda_err(i,PADE);
    cout<<a_over_r0_sqr[i]<<"+/-"<<a_over_r0_sqr_err[i]<<"    "<<methodIIIP_r0Lambda[i]<<"+/-"<<methodIIIP_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;

//  exit(1);

  TCanvas* methodIIIPcanv = new TCanvas("methodIIIPcanv","methodIIIPcanv",700,700);
  
  TH1F* methodIIIPframe=new TH1F("methodIIIPframe","Method IIIP with N_{f} = 2",1000,200,600);
  methodIIIPframe->SetStats(0);
  methodIIIPframe->SetMinimum(0.45);
  methodIIIPframe->SetMaximum(0.70);
  methodIIIPframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIIIPframe->GetXaxis()->SetTickLength(0.02);
  methodIIIPframe->GetXaxis()->SetLabelSize(0.020);
  methodIIIPframe->GetXaxis()->SetLimits(-0.002,0.05); //minx/maxx limits
  methodIIIPframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIIIPframe->GetYaxis()->SetTickLength(0.02);
  methodIIIPframe->GetYaxis()->SetLabelSize(0.020);
  methodIIIPframe->Draw(" ");

  TGraphErrors* methodIIIPgraph=new TGraphErrors(12,a_over_r0_sqr,methodIIIP_r0Lambda,a_over_r0_sqr_err,methodIIIP_r0Lambda_err);
  methodIIIPgraph->SetMarkerColor(1);
  methodIIIPgraph->SetMarkerStyle(33);
  methodIIIPgraph->SetMarkerSize(1);
  methodIIIPgraph->Draw("P");

  Double_t linexIIIP[2]={0.,0.};
  Double_t lineyIIIP[2]={0.,1.};
  TGraph* methodIIIPline=new TGraph(2,linexIIIP,lineyIIIP);
  methodIIIPline->SetLineColor(1);
  methodIIIPline->SetLineStyle(3);
  methodIIIPline->Draw("L");

  // Minimize the statistics with Minuit and thus find the best estimator for x_0
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuitIIIP = TVirtualFitter::Fitter(0,2);
  minuitIIIP->SetParameter(0,"a",0.0, 0.1, 0., 5.);
  minuitIIIP->SetParameter(1,"b",0.0, 0.1,-5., 5.);
  minuitIIIP->SetFCN(minuit_methodIIIP_chisqr);

  // minimize
  Double_t arglistIIIP[100];
  minuitIIIP->ExecuteCommand("MIGRAD",arglistIIIP,0);
  minuitIIIP->ExecuteCommand("MINOS",arglistIIIP,0);

  // local variables
  Double_t pvalIIIP[2],perrIIIP[2],ploIIIP[2],phIIIPI[2],pcoIIIP[2];
  
  // Fit results
  pvalIIIP[0] = minuitIIIP->GetParameter(0);
  pvalIIIP[1] = minuitIIIP->GetParameter(1);
  minuitIIIP->GetErrors(0,phIIIPI[0],ploIIIP[0],perrIIIP[0],pcoIIIP[0]);
  minuitIIIP->GetErrors(1,phIIIPI[1],ploIIIP[1],perrIIIP[1],pcoIIIP[1]);
  cout << endl << pvalIIIP[0] << " +/- " << perrIIIP[0] << " + " << phIIIPI[0] << " " << ploIIIP[0] << " " << pcoIIIP[0] << endl;
  cout << pvalIIIP[1] << " +/- " << perrIIIP[1] << " + " << phIIIPI[1] << " " << ploIIIP[1] << " " << pcoIIIP[1] << endl << endl;

  TF1* methodIIIP_fit=new TF1("methodIIIP_fit",TF1_fitline,0.,0.045,2);
  methodIIIP_fit->SetParameter(0,pvalIIIP[0]);
  methodIIIP_fit->SetParameter(1,pvalIIIP[1]);
  methodIIIP_fit->SetLineWidth(1.5);
  methodIIIP_fit->SetLineColor(1);
  methodIIIP_fit->SetLineStyle(2);
  methodIIIP_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_xIIIP[1]={0.};
  Double_t extrap_yIIIP[1]={pvalIIIP[0]};
  Double_t extrap_x_errIIIP[1]={0.};
  Double_t extrap_y_errIIIP[1]={perrIIIP[0]};
  TGraphErrors* extrapIIIP=new TGraphErrors(1,extrap_xIIIP,extrap_yIIIP,extrap_x_errIIIP,extrap_y_errIIIP);
  extrapIIIP->SetMarkerColor(1);
  extrapIIIP->SetMarkerStyle(27);
  extrapIIIP->SetMarkerSize(1);
  extrapIIIP->Draw("P");

  methodIIIPcanv->SaveAs("methodIIIP_nf=2.pdf");

  cout << " ---- " << b2_boost_III(-1) << " ---- " << endl;
  

/*
  BetaFn bf_deriv(3+PADE,NFLAV,2);
  bf_deriv[2]=b2_boost_III(0);
  bf_deriv[3]=b3_boost_III(0);

  for(Int_t i=0;i<NBETA;++i){
    Double_t deriv=-0.5*methodIIIP_r0Lambda[i]*r0_over_a_extrap[i]*r0_over_a_extrap[i];
    Double_t dpdr0a=0.0119537;
    Double_t prefac=-0.5*methodIIIP_r0Lambda[i]*r0_over_a_extrap[i]*r0_over_a_extrap[i]*r0_over_a_extrap[i]*dpdr0a;
    Double_t parenterm1=3.*c_sw_boosted(i)*NFLAV*(0.0050467-2.*0.0298435*c_sw_boosted(i))/(8.*bf_deriv[0]*g_boost(i)*g_boost(i)*plaquette_extrap[i]);
    Double_t parenterm2=t1_boost(i)/(2.*bf_deriv[0]*g_boost(i)*g_boost(i)*plaquette_extrap[i]);
    Double_t parenterm3=g_boost(i)/(2.*plaquette_extrap[i]*bf_deriv.eval(g_boost(i)));
    Double_t parenterm41=bf_deriv[1]*NFLAV*(0.00050467-2.*0.0298435)*3./(4.*plaquette_extrap[i]);
    Double_t parenterm42=-bf_deriv[0]*NFLAV*(0.001053-2.*0.000498+3.*0.000474+4.*0.000104)*3./(4.*plaquette_extrap[i]);
    Double_t parenterm43=-bf_deriv[0]*(0.2659-0.25)*(-NFLAV*(2.*0.0298435))*3./(4.*plaquette_extrap[i]);
    Double_t parenterm4=(parenterm41+parenterm42+parenterm43)*partial_b2_boost_III_factor(i,PADE);
    cout << deriv << "  " << prefac << "  " << parenterm1 << "  " << parenterm2 << "  " << parenterm3 << "  " << parenterm4 << endl;
    deriv+=(prefac*(parenterm1+parenterm2+parenterm3-parenterm4));
    cout << deriv << endl << endl;
  }
*/

  return 0;

};


