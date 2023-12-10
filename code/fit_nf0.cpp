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

#define NFLAV 0

#define NOTPADE 0
#define PADE 1

#define DIRECT 0
#define RECIPROCAL 1

Double_t RATIO=1.00; //for testing the sensitivity to scale

#define NDATA 11

//data from Table 1 of arxiv:hep-ph/0502212v2
Double_t beta[NDATA]         ={5.70    ,5.80    ,5.95    ,6.00    ,6.07    ,6.20    ,6.40    , 6.57    , 6.69    , 6.81    , 6.92    };
Double_t r0_over_a[NDATA]    ={2.922   ,3.673   ,4.898   ,5.368   ,6.033   ,7.380   ,9.740   ,12.18    ,14.20    ,16.54    ,19.13    };
Double_t r0_over_a_err[NDATA]={0.009   ,0.005   ,0.012   ,0.033   ,0.017   ,0.026   ,0.050   , 0.10    , 0.12    , 0.12    , 0.15    };
Double_t plaquette[NDATA]    ={0.549195,0.567651,0.588006,0.593679,0.601099,0.613633,0.630633, 0.643524, 0.651936, 0.659877, 0.666721};
Double_t plaquette_err[NDATA]={0.000025,0.000021,0.000020,0.000008,0.000018,0.000002,0.000004, 0.000015, 0.000015, 0.000013, 0.000012};

//for nf=2 usage
Double_t kappa[NDATA]        ={1.      ,1.      ,1.      ,1.      ,1.      ,1.      ,1.      , 1.      , 1.      , 1.      , 1.      };
Double_t kappa_c[NDATA]      ={1.      ,1.      ,1.      ,1.      ,1.      ,1.      ,1.      , 1.      , 1.      , 1.      , 1.      };  
Double_t kappa_c_err[NDATA]  ={0.      ,0.      ,0.      ,0.      ,0.      ,0.      ,0.      , 0.      , 0.      , 0.      , 0.      };  

//results need to be global for minuit to work properly
Double_t methodI_lat_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodI_lat_r0Lambda_err[NDATA];
Double_t methodI_lat_a_over_r0_sqr[NDATA];
Double_t methodI_lat_a_over_r0_sqr_err[NDATA];

//results need to be global for minuit to work properly
Double_t methodI_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodI_r0Lambda_err[NDATA];
Double_t methodI_a_over_r0_sqr[NDATA];
Double_t methodI_a_over_r0_sqr_err[NDATA];

//results need to be global for minuit to work properly
Double_t methodII_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodII_r0Lambda_err[NDATA];
Double_t methodII_a_over_r0_sqr[NDATA];
Double_t methodII_a_over_r0_sqr_err[NDATA];

//results need to be global for minuit to work properly
Double_t methodIIP_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodIIP_r0Lambda_err[NDATA];
Double_t methodIIP_a_over_r0_sqr[NDATA];
Double_t methodIIP_a_over_r0_sqr_err[NDATA];

//results need to be global for minuit to work properly
Double_t methodIII_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodIII_r0Lambda_err[NDATA];
Double_t methodIII_a_over_r0_sqr[NDATA];
Double_t methodIII_a_over_r0_sqr_err[NDATA];

//results need to be global for minuit to work properly
Double_t methodIIIP_r0Lambda[NDATA]; //don't use the first two datapoints since you start to see order(a^2) artifacts
Double_t methodIIIP_r0Lambda_err[NDATA];
Double_t methodIIIP_a_over_r0_sqr[NDATA];
Double_t methodIIIP_a_over_r0_sqr_err[NDATA];

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

// here follows many usefull variables relating g_msbar to g0 and g_boost
Double_t u0(Int_t ind){
  return TMath::Power(plaquette[ind],0.25);
};
Double_t u0_err(Int_t ind){
  Double_t partial_plaquette=0.25*TMath::Power(plaquette[ind],-0.75);
  Double_t term1=TMath::Power(partial_plaquette,2.)*TMath::Power(plaquette_err[ind],2.);
  return TMath::Sqrt(term1);
};

Double_t c_sw(Int_t ind){
  Double_t gsqr=6./beta[ind];
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
  return 0.5*( (1./kappa[ind])-(1./kappa_c[ind]) );
};
Double_t amq_err(Int_t ind){  // no correlations
  Double_t partial_kappa_c=0.5*TMath::Power( kappa_c[ind] , -2. );
  Double_t term1=TMath::Power(partial_kappa_c,2.)*TMath::Power(kappa_c_err[ind],2.);
  return TMath::Sqrt(term1);
};

Double_t t1_lat_chirallimit(Int_t ind){
  Double_t csw=c_sw(ind);
//  if(NFLAV>0){csw=c_sw(ind);}  // I don't see why we need to do this
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
  Double_t cswb;
  if(ind<0){cswb=val;}
  else{cswb=c_sw_boosted(ind);}
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
  return TMath::Sqrt( 6./beta[ind] );
};
Double_t g_lat_err(Int_t ind){
  return 0.;
};
Double_t g_boost(Int_t ind){
  return TMath::Sqrt( 6./(beta[ind]*plaquette[ind]) );
};
Double_t g_boost_err(Int_t ind){
  Double_t partial_plaquette=-6./(beta[ind]*plaquette[ind]*plaquette[ind]);
  Double_t term1=TMath::Power(partial_plaquette,2.)*TMath::Power(plaquette_err[ind],2.);
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

// -------------------------- method 0 ---------------------
// This method is like Roger's method I but without using the 
// boosted coupling, rather it just uses the bare lattice coupling
Double_t g_ms_I_lat(Int_t ind,Int_t opt=0){
  // opt==0 -- use direct expansion of g_MS ie g_MS^2=g_lat^2+d_1*g_lat^4+...
  // opt==1 -- use reciprocal expansion of g_MS ie 1/g_MS^2=1/g_lat^2+...  -- this one agrees with the r0Lambda values in the paper
  BetaFn bf(4,NFLAV,2);
  if(opt==RECIPROCAL){
    return TMath::Power( TMath::Power( g_lat(ind) , -2. ) + (2.*bf[0]*TMath::Log(RATIO))
                        + (2.*bf[1]*TMath::Log(RATIO)-t2_lat(ind)+bf[1]*t1_lat(ind)/bf[0])*TMath::Power(g_lat(ind),2.) , -0.5 );
  }
  // opt==DIRECT
  return TMath::Sqrt( TMath::Power( g_lat(ind) , 2. ) - (2.*bf[0]*TMath::Log(RATIO))*TMath::Power( g_lat(ind) , 4. )
                        + (4.*bf[0]*bf[0]*TMath::Log(RATIO)*TMath::Log(RATIO)-2.*bf[1]*TMath::Log(RATIO)+t2_lat(ind)-bf[1]*t1_lat(ind)/bf[0])*TMath::Power(g_lat(ind),6.) );
};

Double_t g_ms_I_lat_err(Int_t ind,Int_t opt=0){
  BetaFn bf(4,NFLAV,2);
  // opt==0 -- use direct expansion of g_MS ie g_MS^2=g_lat^2+d_1*g_lat^4+...
  // opt==1 -- use reciprocal expansion of g_MS ie 1/g_MS^2=1/g_lat^2+...  -- this one agrees with the r0Lambda values in the paper
  if(opt==RECIPROCAL){

    Double_t prefac=0.25*TMath::Power(g_ms_I_lat(ind,opt),6.);

    Double_t partial_g_lat=-2.*TMath::Power(g_lat(ind),-3.)
                    +2.*(2.*bf[1]*TMath::Log(RATIO)-t2_lat(ind)+bf[1]*t1_lat(ind)/bf[0])*TMath::Power( g_lat(ind) , 1. );
    Double_t partial_t1_lat=bf[1]*TMath::Power( g_lat(ind) , 2. )/bf[0];
    Double_t partial_t2_lat=-TMath::Power( g_lat(ind) , 2. );
 
    Double_t corr_g_lat__t1_lat=1.;
    Double_t corr_g_lat__t2_lat=1.;
    Double_t corr_t1_lat__t2_lat=1.;
  
    Double_t term1=TMath::Power(partial_g_lat,2.)*TMath::Power(g_lat_err(ind),2.);
    Double_t term2=TMath::Power(partial_t1_lat,2.)*TMath::Power(t1_lat_err(ind),2.);
    Double_t term3=TMath::Power(partial_t2_lat,2.)*TMath::Power(t2_lat_err(ind),2.);

    Double_t termcorr1=2.*partial_g_lat*partial_t1_lat*g_lat_err(ind)*t1_lat_err(ind)*corr_g_lat__t1_lat;
    Double_t termcorr2=2.*partial_g_lat*partial_t2_lat*g_lat_err(ind)*t2_lat_err(ind)*corr_g_lat__t2_lat;
    Double_t termcorr3=2.*partial_t1_lat*partial_t2_lat*t1_lat_err(ind)*t2_lat_err(ind)*corr_t1_lat__t2_lat;

    return TMath::Sqrt(prefac*(term1+term2+term3+termcorr1+termcorr2+termcorr3));
  }
  Double_t partial_g_lat=2.*g_lat(ind)+6.*(t2_lat(ind)-bf[1]*t1_lat(ind)/bf[0])*TMath::Power( g_lat(ind) , 5. );
  Double_t partial_t1_lat=bf[1]*TMath::Power( g_lat(ind) , 6. )/bf[0];
  Double_t partial_t2_lat=TMath::Power( g_lat(ind) , 6. );
 
  Double_t corr_g_lat__t1_lat=1.;
  Double_t corr_g_lat__t2_lat=1.;
  Double_t corr_t1_lat__t2_lat=1.;
  
  Double_t term1=TMath::Power(partial_g_lat,2.)*TMath::Power(g_lat_err(ind),2.)/TMath::Power(2.*g_ms_I_lat(ind,opt),2.);
  Double_t term2=TMath::Power(partial_t1_lat,2.)*TMath::Power(t1_lat_err(ind),2.)/TMath::Power(2.*g_ms_I_lat(ind,opt),2.);
  Double_t term3=TMath::Power(partial_t2_lat,2.)*TMath::Power(t2_lat_err(ind),2.)/TMath::Power(2.*g_ms_I_lat(ind,opt),2.);

  Double_t termcorr1=2.*partial_g_lat*partial_t1_lat*g_lat_err(ind)*t1_lat_err(ind)*corr_g_lat__t1_lat;
  Double_t termcorr2=2.*partial_g_lat*partial_t2_lat*g_lat_err(ind)*t2_lat_err(ind)*corr_g_lat__t2_lat;
  Double_t termcorr3=2.*partial_t1_lat*partial_t2_lat*t1_lat_err(ind)*t2_lat_err(ind)*corr_t1_lat__t2_lat;

  return TMath::Sqrt(term1+term2+term3+termcorr1+termcorr2+termcorr3);
};

Double_t methodI_lat_getr0Lambda(Int_t ind,Int_t opt=0){
  BetaFn bf(4,NFLAV,2);
  Double_t fac1=r0_over_a[ind];
  Double_t fac2=TMath::Exp(t1_lat(ind)/(2.*bf[0]));
  Double_t fac3=1./mu_over_lambda(g_ms_I_lat(ind,opt),bf);
  return RATIO*fac1*fac2*fac3;
};
Double_t methodI_lat_getr0Lambda_err(Int_t ind, Int_t opt=0){   //did you take the covariance between t1_lat and g_ms_I??? Yep!
  BetaFn bf(4,NFLAV,2);

  Double_t partial_r0_over_a=RATIO*methodI_lat_getr0Lambda(ind,opt)/r0_over_a[ind];
  Double_t partial_t1_lat=RATIO*methodI_lat_getr0Lambda(ind,opt)/(2.*bf[0]);
  Double_t partial_g_ms_I_lat= -RATIO*methodI_lat_getr0Lambda(ind,opt)/bf.eval(g_ms_I_lat(ind,opt));  

  Double_t corr_r0_over_a__t1_lat=0.;
  Double_t corr_r0_over_a__g_ms_I_lat=0.;
  Double_t corr_t1_lat__g_ms_I_lat=1.;

  Double_t term1=TMath::Power(partial_r0_over_a,2.)*TMath::Power(r0_over_a_err[ind],2.);
  Double_t term2=TMath::Power(partial_t1_lat,2.)*TMath::Power(t1_lat_err(ind),2.);
  Double_t term3=TMath::Power(partial_g_ms_I_lat,2.)*TMath::Power(g_ms_I_lat_err(ind,ind),2.);

  Double_t termcorr1=2.*partial_r0_over_a*partial_t1_lat*r0_over_a_err[ind]*t1_lat_err(ind)*corr_r0_over_a__t1_lat;
  Double_t termcorr2=2.*partial_r0_over_a*partial_g_ms_I_lat*r0_over_a_err[ind]*g_ms_I_lat_err(ind,opt)*corr_r0_over_a__g_ms_I_lat;
  Double_t termcorr3=2.*partial_t1_lat*partial_g_ms_I_lat*t1_lat_err(ind)*g_ms_I_lat_err(ind,opt)*corr_t1_lat__g_ms_I_lat;

  return TMath::Sqrt(term1+term2+term3+termcorr1+termcorr2+termcorr3);
};

// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_methodI_lat_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=2;i<NDATA;++i) { //skip the first two due to order(a^2) effects
    num=TMath::Power( methodI_lat_r0Lambda[i] - fit_parabola(methodI_lat_a_over_r0_sqr[i],p[0],p[1],p[2]) , 2. );
    den=TMath::Power(methodI_lat_r0Lambda_err[i],2.)
       +TMath::Power(methodI_lat_a_over_r0_sqr_err[i]*fit_parabola_deriv(methodI_lat_a_over_r0_sqr[i],p[0],p[1],p[2]),2.);
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
  Double_t fac1=RATIO*r0_over_a[ind];
  Double_t fac2=TMath::Exp(t1_boost(ind)/(2.*bf[0]));
  Double_t fac3=1./mu_over_lambda(g_ms_I(ind,opt),bf);
  return fac1*fac2*fac3;
};
Double_t methodI_getr0Lambda_err(Int_t ind,Int_t opt=0){   //did you take the covariance between t1_boost and g_ms_I??? Yep!
  BetaFn bf(4,NFLAV,2);

  Double_t partial_r0_over_a=RATIO*methodI_getr0Lambda(ind)/r0_over_a[ind];
  Double_t partial_t1_boost=RATIO*methodI_getr0Lambda(ind)/(2.*bf[0]);
  Double_t partial_g_ms_I= -RATIO*methodI_getr0Lambda(ind)/bf.eval(g_ms_I(ind,opt));  

  Double_t corr_r0_over_a__t1_boost=0.;
  Double_t corr_r0_over_a__g_ms_I=0.;
  Double_t corr_t1_boost__g_ms_I=1.;

  Double_t term1=TMath::Power(partial_r0_over_a,2.)*TMath::Power(r0_over_a_err[ind],2.);
  Double_t term2=TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.);
  Double_t term3=TMath::Power(partial_g_ms_I,2.)*TMath::Power(g_ms_I_err(ind,opt),2.);

  Double_t termcorr1=2.*partial_r0_over_a*partial_t1_boost*r0_over_a_err[ind]*t1_boost_err(ind)*corr_r0_over_a__t1_boost;
  Double_t termcorr2=2.*partial_r0_over_a*partial_g_ms_I*r0_over_a_err[ind]*g_ms_I_err(ind,opt)*corr_r0_over_a__g_ms_I;
  Double_t termcorr3=2.*partial_t1_boost*partial_g_ms_I*t1_boost_err(ind)*g_ms_I_err(ind,opt)*corr_t1_boost__g_ms_I;

  return TMath::Sqrt(term1+term2+term3+termcorr1+termcorr2+termcorr3);
};

// *Okay, we're going to implement the "effective variance" approach. Here we use
// \chi^2 = sum_i^N [y_i - f(x_i)]^2/[\sigma_y^2+(\sigma_x*f'(x_i))^2]
// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_methodI_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=2;i<NDATA;++i) { //skip the first two due to order(a^2) effects
//    if( (i==3) || (i==10) ){ continue;} // Roger's neglected points
    num=TMath::Power( methodI_r0Lambda[i] - fitline(methodI_a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodI_r0Lambda_err[i],2.)+TMath::Power(methodI_a_over_r0_sqr_err[i]*p[1],2.);
    chi2 += num/den;
  }
  fval = chi2;
}

// ----------------------- method I - dependence on scale setting --------------------------

Double_t TF1_r0Lambda_of_muratio_I(Double_t* x, Double_t* par){
  Double_t temprat=RATIO;
  RATIO=x[0];
  for(Int_t i=0;i<NDATA;++i){
    methodI_r0Lambda[i]=methodI_getr0Lambda(i,DIRECT);
    methodI_r0Lambda_err[i]=methodI_getr0Lambda_err(i,DIRECT);
    methodI_a_over_r0_sqr[i]=TMath::Power(r0_over_a[i],-2.);
    methodI_a_over_r0_sqr_err[i]=2.*r0_over_a_err[i]/(r0_over_a[i]*r0_over_a[i]*r0_over_a[i]);
  }
  // Minimize the statistics with Minuit and thus find the best estimator for x_0
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuit = TVirtualFitter::Fitter(0,2);
  minuit->SetParameter(0,"a",0.0, 0.001, 0., 10.);
  minuit->SetParameter(1,"b",0.0, 0.001, -10., 0.);
  minuit->SetFCN(minuit_methodI_chisqr);

  // minimize
  Double_t arglist[100];
  minuit->ExecuteCommand("MIGRAD",arglist,0);
  minuit->ExecuteCommand("MINOS",arglist,0);
  
  // Fit results
  RATIO=temprat;
  Double_t val=minuit->GetParameter(0);
  return val;
};
Double_t TF1_r0Lambda_of_muratio_I_deriv(Double_t* x, Double_t* par){
  TF1 scratch("scratch",TF1_r0Lambda_of_muratio_I,0.01,40.,0);
  return scratch.Derivative(x[0],par,0.0001);
};
// ----------------------- method II --------------------------

Double_t g_ms_II(Int_t ind){
  BetaFn bf(4,NFLAV,2);
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

  Double_t fac1=r0_over_a[ind];
  Double_t fac2=TMath::Exp(t1_boost(ind)/(2.*bf_ms[0]));
  Double_t fac3=mu_over_lambda(g_boost(ind),bf_ms)/(mu_over_lambda(g_ms_II(ind),bf_ms)*mu_over_lambda(g_boost(ind),bf_boost) );
  return RATIO*fac1*fac2*fac3;
};
Double_t methodII_getr0Lambda_err(Int_t ind,Int_t pade){
  BetaFn bf_ms(4,NFLAV,2);
  BetaFn bf_boost(3,NFLAV,2);
  bf_boost[2]=b2_boost_II(ind);
  if(pade){bf_boost[3]=bf_boost[2]*bf_boost[2]/bf_boost[1];}

  Double_t partial_r0_over_a=methodII_getr0Lambda(ind,pade)/r0_over_a[ind];
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

  Double_t term1=TMath::Power(partial_r0_over_a,2.)*TMath::Power(r0_over_a_err[ind],2.);
  Double_t term2=TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.);
  Double_t term3=TMath::Power(partial_g_boost,2.)*TMath::Power(g_boost_err(ind),2.);
  Double_t term4=TMath::Power(partial_b2_boost,2.)*TMath::Power(b2_boost_II_err(ind),2.);
  Double_t term5=TMath::Power(partial_b3_boost,2.)*TMath::Power(b3_boost_II_err(ind),2.);

  Double_t termcorr1=2.*partial_r0_over_a*partial_t1_boost*r0_over_a_err[ind]*t1_boost_err(ind)*corr_r0_over_a__t1_boost;
  Double_t termcorr2=2.*partial_r0_over_a*partial_g_boost*r0_over_a_err[ind]*g_boost_err(ind)*corr_r0_over_a__g_boost;
  Double_t termcorr3=2.*partial_r0_over_a*partial_b2_boost*r0_over_a_err[ind]*b2_boost_II_err(ind)*corr_r0_over_a__b2_boost;
  Double_t termcorr4=2.*partial_r0_over_a*partial_b3_boost*r0_over_a_err[ind]*b3_boost_II_err(ind)*corr_r0_over_a__b3_boost;
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
  for(Int_t i=2;i<NDATA;++i) { //skip the first two due to order(a^2) effects
    num=TMath::Power( methodII_r0Lambda[i] - fitline(methodII_a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodII_r0Lambda_err[i],2.)+TMath::Power(methodII_a_over_r0_sqr_err[i]*p[1],2.);
    chi2 += num/den;
  }
  fval = chi2;
}
// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_methodIIP_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=2;i<NDATA;++i) { //skip the first two due to order(a^2) effects
    num=TMath::Power( methodIIP_r0Lambda[i] - fitline(methodIIP_a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodIIP_r0Lambda_err[i],2.)+TMath::Power(methodIIP_a_over_r0_sqr_err[i]*p[1],2.);
    chi2 += num/den;
  }
  fval = chi2;
}

// ----------------------- method II - dependence on scale setting --------------------------

Double_t TF1_r0Lambda_of_muratio_II(Double_t* x, Double_t* par){
  Double_t temprat=RATIO;
  RATIO=x[0];
  for(Int_t i=0;i<NDATA;++i){
    methodII_r0Lambda[i]=methodII_getr0Lambda(i,NOTPADE);
    methodII_r0Lambda_err[i]=methodII_getr0Lambda_err(i,NOTPADE);
    methodII_a_over_r0_sqr[i]=TMath::Power(r0_over_a[i],-2.);
    methodII_a_over_r0_sqr_err[i]=2.*r0_over_a_err[i]/(r0_over_a[i]*r0_over_a[i]*r0_over_a[i]);
  }
  // Minimize the statistics with Minuit and thus find the best estimator for x_0
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuit = TVirtualFitter::Fitter(0,2);
  minuit->SetParameter(0,"a",0.0, 0.1, 0., 10.);
  minuit->SetParameter(1,"b",0.0, 0.1,-5., 5.);
  minuit->SetFCN(minuit_methodII_chisqr);

  // minimize
  Double_t arglist[100];
  minuit->ExecuteCommand("MIGRAD",arglist,0);
  minuit->ExecuteCommand("MINOS",arglist,0);
  
  // Fit results
  RATIO=temprat;
  return minuit->GetParameter(0);
};
Double_t TF1_r0Lambda_of_muratio_II_deriv(Double_t* x, Double_t* par){
  TF1 scratch("scratch",TF1_r0Lambda_of_muratio_II,0.01,40.,0);
  return scratch.Derivative(x[0],par,0.0001);
};

// ----------------------- method III --------------------------

Double_t g_ms_III(Int_t ind,Int_t opt=0){
  BetaFn bf(4,NFLAV,2);
  if(opt>0){
    //put in reciprocal stuff here
  }
  return TMath::Sqrt( TMath::Power( g_boost(ind) , 2. ) - (2.*bf[0]*TMath::Log(RATIO))*TMath::Power( g_boost(ind) , 4. )
                        + (4.*bf[0]*bf[0]*TMath::Log(RATIO)*TMath::Log(RATIO)-2.*bf[1]*TMath::Log(RATIO))*TMath::Power(g_boost(ind),6.) );
};
Double_t TF1_partial_b2_boost_III_factor_integrand(Double_t* x,Double_t* par){
// par[0]==ind, par[1]==pade?
  Int_t pade=(Int_t) par[1];
  BetaFn bf_boost(3+pade,NFLAV,2);
  bf_boost[2]=b2_boost_III(par[0]);
  if(pade==1){bf_boost[3]=b3_boost_III(par[0]);}
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
  if(pade){bf_boost[3]=bf_boost[2]*bf_boost[2]/bf_boost[1];}

  Double_t fac1=r0_over_a[ind];
  Double_t fac2=TMath::Exp(t1_boost(ind)/(2.*bf_ms[0]));
  Double_t fac3=mu_over_lambda(g_boost(ind),bf_ms)/(mu_over_lambda(g_ms_III(ind),bf_ms)*mu_over_lambda(g_boost(ind),bf_boost) );
  return RATIO*fac1*fac2*fac3;
};
Double_t methodIII_getr0Lambda_err(Int_t ind,Int_t pade){
  BetaFn bf_ms(4,NFLAV,2);
  BetaFn bf_boost(3,NFLAV,2);
  bf_boost[2]=b2_boost_III(ind);
  if(pade){bf_boost[3]=bf_boost[2]*bf_boost[2]/bf_boost[1];}

  Double_t partial_r0_over_a=methodIII_getr0Lambda(ind,pade)/r0_over_a[ind];
  Double_t partial_t1_boost=methodIII_getr0Lambda(ind,pade)/(2.*bf_ms[0]);
  Double_t partial_g_boost= -methodIII_getr0Lambda(ind,pade)/bf_boost.eval(g_boost(ind)); 
  Double_t partial_b2_boost= -methodIII_getr0Lambda(ind,pade)*partial_b2_boost_III_factor(ind,pade); 
  Double_t partial_b3_boost= -methodIII_getr0Lambda(ind,pade)*partial_b3_boost_III_factor(ind,pade); 

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

  Double_t term1=TMath::Power(partial_r0_over_a,2.)*TMath::Power(r0_over_a_err[ind],2.);
  Double_t term2=TMath::Power(partial_t1_boost,2.)*TMath::Power(t1_boost_err(ind),2.);
  Double_t term3=TMath::Power(partial_g_boost,2.)*TMath::Power(g_boost_err(ind),2.);
  Double_t term4=TMath::Power(partial_b2_boost,2.)*TMath::Power(b2_boost_III_err(ind),2.);
  Double_t term5=TMath::Power(partial_b3_boost,2.)*TMath::Power(b3_boost_III_err(ind),2.);

  Double_t termcorr1=2.*partial_r0_over_a*partial_t1_boost*r0_over_a_err[ind]*t1_boost_err(ind)*corr_r0_over_a__t1_boost;
  Double_t termcorr2=2.*partial_r0_over_a*partial_g_boost*r0_over_a_err[ind]*g_boost_err(ind)*corr_r0_over_a__g_boost;
  Double_t termcorr3=2.*partial_r0_over_a*partial_b2_boost*r0_over_a_err[ind]*b2_boost_III_err(ind)*corr_r0_over_a__b2_boost;
  Double_t termcorr4=2.*partial_r0_over_a*partial_b3_boost*r0_over_a_err[ind]*b3_boost_III_err(ind)*corr_r0_over_a__b3_boost;
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
  for(Int_t i=2;i<NDATA;++i) { //skip the first two due to order(a^2) effects
//    if( (i==3) || (i==10) ){continue;}  //Roger's skipped points
    num=TMath::Power( methodIII_r0Lambda[i] - fitline(methodIII_a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodIII_r0Lambda_err[i],2.)/*+TMath::Power(methodIII_a_over_r0_sqr_err[i]*p[1],2.)*/;
    chi2 += num/den;
  }
  fval = chi2;
}
// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_methodIIIP_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  for(Int_t i=2;i<NDATA;++i) { //skip the first two due to order(a^2) effects
//    if( (i==3) || (i==10) ){continue;}  //Roger's skipped points
    num=TMath::Power( methodIIIP_r0Lambda[i] - fitline(methodIIIP_a_over_r0_sqr[i],p[0],p[1]) , 2. );
    den=TMath::Power(methodIIIP_r0Lambda_err[i],2.)+TMath::Power(methodIIIP_a_over_r0_sqr_err[i]*p[1],2.);
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

// I want to show that using the boosted coupling actually improves convergence
// *NOTE* I used improved_convergence.cpp to check how badly the g_0 to g_ms series converges, and it's not bad at all
//  cout<<g_ms_I_lat(2,0)<<endl;
//  cout<<g_ms_I_lat(2,1)<<endl;
//  cout<<g_ms_I(2,0)<<endl;
//  cout<<g_ms_I(2,1)<<endl;

// NOTE -- If you want accurate results for the Method I, II, IIP, III, and IIIP, then you need to comment out
//            the method I lat method because it causes minuit to always have 3 minimization parameters

/*
  cout << "\n\n ----- BEGIN METHOD I Lattice ------ \n\n";

  for(Int_t i=0;i<NDATA;++i){
    methodI_lat_r0Lambda[i]=methodI_lat_getr0Lambda(i,DIRECT);
    methodI_lat_r0Lambda_err[i]=methodI_lat_getr0Lambda_err(i,DIRECT);
    methodI_lat_a_over_r0_sqr[i]=TMath::Power(r0_over_a[i],-2.);
    methodI_lat_a_over_r0_sqr_err[i]=2.*r0_over_a_err[i]/(r0_over_a[i]*r0_over_a[i]*r0_over_a[i]);
    cout<<methodI_lat_a_over_r0_sqr[i]<<"+/-"<<methodI_lat_a_over_r0_sqr_err[i]<<"    "<<methodI_lat_r0Lambda[i]<<"+/-"<<methodI_lat_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;

  TCanvas* methodI_latcanv = new TCanvas("methodI_latcanv","methodI_latcanv",700,700);
  
  TH1F* methodI_latframe=new TH1F("methodI_latframe"," ",1000,200,600);
  methodI_latframe->SetStats(0);
  methodI_latframe->SetMinimum(0.3);
  methodI_latframe->SetMaximum(0.85);
  methodI_latframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodI_latframe->GetXaxis()->SetTickLength(0.02);
  methodI_latframe->GetXaxis()->SetLabelSize(0.020);
  methodI_latframe->GetXaxis()->SetLimits(-0.005,0.12); //minx/maxx limits
  methodI_latframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodI_latframe->GetYaxis()->SetTickLength(0.02);
  methodI_latframe->GetYaxis()->SetLabelSize(0.020);
  methodI_latframe->Draw(" ");

  TGraphErrors* methodI_latgraph=new TGraphErrors(NDATA,methodI_lat_a_over_r0_sqr,methodI_lat_r0Lambda,methodI_lat_a_over_r0_sqr_err,methodI_lat_r0Lambda_err);
  methodI_latgraph->SetMarkerColor(1);
  methodI_latgraph->SetMarkerStyle(33);
  methodI_latgraph->SetMarkerSize(1);
  methodI_latgraph->Draw("P");

  Double_t linexI_lat[2]={0.,0.};
  Double_t lineyI_lat[2]={0.,1.};
  TGraph* methodI_latline=new TGraph(2,linexI_lat,lineyI_lat);
  methodI_latline->SetLineColor(1);
  methodI_latline->SetLineStyle(3);
  methodI_latline->Draw("L");

  // Minimize the statistics with Minuit and thus find the best estimator for x_0
  TVirtualFitter::SetDefaultFitter("Minuit");
  TVirtualFitter* minuitI_lat = TVirtualFitter::Fitter(0,3);
  minuitI_lat->SetParameter(0,"p0",0.0, 0.1, -1., 1.);
  minuitI_lat->SetParameter(1,"p1",0.0, 0.1,-50., 50.);
  minuitI_lat->SetParameter(2,"p2",0.0, 1.,-500., 500.);
  minuitI_lat->SetFCN(minuit_methodI_lat_chisqr);

  // minimize
  Double_t arglistI_lat[100];
  minuitI_lat->ExecuteCommand("MIGRAD",arglistI_lat,0);
  minuitI_lat->ExecuteCommand("MINOS",arglistI_lat,0);

  // local variables
  Double_t pvalI_lat[3],perrI_lat[3],ploI_lat[3],phiI_lat[3],pcoI_lat[3];
  
  // Fit results
  pvalI_lat[0] = minuitI_lat->GetParameter(0);
  pvalI_lat[1] = minuitI_lat->GetParameter(1);
  pvalI_lat[2] = minuitI_lat->GetParameter(2);
  minuitI_lat->GetErrors(0,phiI_lat[0],ploI_lat[0],perrI_lat[0],pcoI_lat[0]);
  minuitI_lat->GetErrors(1,phiI_lat[1],ploI_lat[1],perrI_lat[1],pcoI_lat[1]);
  minuitI_lat->GetErrors(2,phiI_lat[2],ploI_lat[2],perrI_lat[2],pcoI_lat[2]);
  cout << endl << pvalI_lat[0] << " +/- " << perrI_lat[0] << " + " << phiI_lat[0] << " " << ploI_lat[0] << " " << pcoI_lat[0] << endl;
  cout << pvalI_lat[1] << " +/- " << perrI_lat[1] << " + " << phiI_lat[1] << " " << ploI_lat[1] << " " << pcoI_lat[1] << endl;
  cout << pvalI_lat[2] << " +/- " << perrI_lat[2] << " + " << phiI_lat[2] << " " << ploI_lat[2] << " " << pcoI_lat[2] << endl << endl;

  TF1* methodI_lat_fit=new TF1("methodI_lat_fit",TF1_fit_parabola,0.,0.045,3);
  methodI_lat_fit->SetParameter(0,pvalI_lat[0]);
  methodI_lat_fit->SetParameter(1,pvalI_lat[1]);
  methodI_lat_fit->SetParameter(2,pvalI_lat[2]);
  methodI_lat_fit->SetLineWidth(1.5);
  methodI_lat_fit->SetLineColor(1);
  methodI_lat_fit->SetLineStyle(2);
  methodI_lat_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_xI_lat[1]={0.};
  Double_t extrap_yI_lat[1]={pvalI_lat[0]};
  Double_t extrap_x_errI_lat[1]={0.};
  Double_t extrap_y_errI_lat[1]={perrI_lat[0]};
  TGraphErrors* extrapI_lat=new TGraphErrors(1,extrap_xI_lat,extrap_yI_lat,extrap_x_errI_lat,extrap_y_errI_lat);
  extrapI_lat->SetMarkerColor(1);
  extrapI_lat->SetMarkerStyle(27);
  extrapI_lat->SetMarkerSize(1);
  extrapI_lat->Draw("P");

  methodI_latcanv->SaveAs("methodI_lat_quenched.pdf");

  exit(1);  // this has to be here, since the previous fit job screws up the following fit jobs
*/

  cout << "\n\n ----- BEGIN METHOD I Boosted ------ \n\n";

  for(Int_t i=0;i<NDATA;++i){
    methodI_r0Lambda[i]=methodI_getr0Lambda(i,DIRECT);
    methodI_r0Lambda_err[i]=methodI_getr0Lambda_err(i,DIRECT);
    methodI_a_over_r0_sqr[i]=TMath::Power(r0_over_a[i],-2.);
    methodI_a_over_r0_sqr_err[i]=2.*r0_over_a_err[i]/(r0_over_a[i]*r0_over_a[i]*r0_over_a[i]);
    cout<<methodI_a_over_r0_sqr[i]<<"+/-"<<methodI_a_over_r0_sqr_err[i]<<"    "<<methodI_r0Lambda[i]<<"+/-"<<methodI_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;

  TCanvas* methodIcanv = new TCanvas("methodIcanv","methodIcanv",700,700);
  
  TH1F* methodIframe=new TH1F("methodIframe"," ",1000,200,600);
  methodIframe->SetStats(0);
  methodIframe->SetMinimum(0.45);
  methodIframe->SetMaximum(0.70);
  methodIframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIframe->GetXaxis()->SetTickLength(0.02);
  methodIframe->GetXaxis()->SetLabelSize(0.020);
  methodIframe->GetXaxis()->SetLimits(-0.005,0.12); //minx/maxx limits
  methodIframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIframe->GetYaxis()->SetTickLength(0.02);
  methodIframe->GetYaxis()->SetLabelSize(0.020);
  methodIframe->Draw(" ");

  TGraphErrors* methodIgraph=new TGraphErrors(NDATA,methodI_a_over_r0_sqr,methodI_r0Lambda,methodI_a_over_r0_sqr_err,methodI_r0Lambda_err);
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

  methodIcanv->SaveAs("methodIquenched.pdf");
/*
  cout << "\n\n ----- BEGIN METHOD I - testing muratio dependence------ \n\n";

  TCanvas* methodImucanv = new TCanvas("methodImucanv","methodImucanv",700,700);
  methodImucanv->SetLogx();
  
  TH1F* methodImuframe=new TH1F("methodImuframe"," ",1000,200,600);
  methodImuframe->SetStats(0);
  methodImuframe->SetMinimum(0.);
  methodImuframe->SetMaximum(5.);
  methodImuframe->GetXaxis()->SetTitle("#mu/#mu_{I}");
  methodImuframe->GetXaxis()->SetTickLength(0.02);
  methodImuframe->GetXaxis()->SetLabelSize(0.020);
  methodImuframe->GetXaxis()->SetLimits(0.01,30.); //minx/maxx limits
  methodImuframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodImuframe->GetYaxis()->SetTickLength(0.02);
  methodImuframe->GetYaxis()->SetLabelSize(0.020);
  methodImuframe->Draw(" ");

  TF1* methodImu_fit=new TF1("methodImu_fit",TF1_r0Lambda_of_muratio_I,0.01,40.,0);
  methodImu_fit->SetLineWidth(1.5);
  methodImu_fit->SetLineColor(1);
  methodImu_fit->SetLineStyle(1);
  methodImu_fit->Draw("SAME");

  TF1* methodImu_deriv=new TF1("methodImu_fit",TF1_r0Lambda_of_muratio_I_deriv,0.01,40.,0);
  methodImu_deriv->SetLineWidth(1.5);
  methodImu_deriv->SetLineColor(1);
  methodImu_deriv->SetLineStyle(2);
  methodImu_deriv->Draw("SAME");

  methodImucanv->SaveAs("method_I_stability.pdf");
*/
  cout << "\n\n ----- BEGIN METHOD II ------ \n\n";

  for(Int_t i=0;i<NDATA;++i){
    methodII_r0Lambda[i]=methodII_getr0Lambda(i,NOTPADE);
    methodII_r0Lambda_err[i]=methodII_getr0Lambda_err(i,NOTPADE);
    methodII_a_over_r0_sqr[i]=TMath::Power(r0_over_a[i],-2.);
    methodII_a_over_r0_sqr_err[i]=2.*r0_over_a_err[i]/(r0_over_a[i]*r0_over_a[i]*r0_over_a[i]);
    cout<<methodII_a_over_r0_sqr[i]<<"+/-"<<methodII_a_over_r0_sqr_err[i]<<"    "<<methodII_r0Lambda[i]<<"+/-"<<methodII_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;

  TCanvas* methodIIcanv = new TCanvas("methodIIcanv","methodIIcanv",700,700);
  
  TH1F* methodIIframe=new TH1F("methodIIframe"," ",1000,200,600);
  methodIIframe->SetStats(0);
  methodIIframe->SetMinimum(0.45);
  methodIIframe->SetMaximum(0.70);
  methodIIframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIIframe->GetXaxis()->SetTickLength(0.02);
  methodIIframe->GetXaxis()->SetLabelSize(0.020);
  methodIIframe->GetXaxis()->SetLimits(-0.005,0.12); //minx/maxx limits
  methodIIframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIIframe->GetYaxis()->SetTickLength(0.02);
  methodIIframe->GetYaxis()->SetLabelSize(0.020);
  methodIIframe->Draw(" ");

  TGraphErrors* methodIIgraph=new TGraphErrors(12,methodII_a_over_r0_sqr,methodII_r0Lambda,methodII_a_over_r0_sqr_err,methodII_r0Lambda_err);
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

  methodIIcanv->SaveAs("methodIIquenched.pdf");
/*
  cout << "\n\n ----- BEGIN METHOD II - testing muratio dependence------ \n\n";

  TCanvas* methodIImucanv = new TCanvas("methodIImucanv","methodIImucanv",700,700);
  methodIImucanv->SetLogx();
  
  TH1F* methodIImuframe=new TH1F("methodIImuframe"," ",1000,200,600);
  methodIImuframe->SetStats(0);
  methodIImuframe->SetMinimum(0.);
  methodIImuframe->SetMaximum(5.);
  methodIImuframe->GetXaxis()->SetTitle("#mu/#mu_{II}");
  methodIImuframe->GetXaxis()->SetTickLength(0.02);
  methodIImuframe->GetXaxis()->SetLabelSize(0.020);
  methodIImuframe->GetXaxis()->SetLimits(0.01,30.); //minx/maxx limits
  methodIImuframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIImuframe->GetYaxis()->SetTickLength(0.02);
  methodIImuframe->GetYaxis()->SetLabelSize(0.020);
  methodIImuframe->Draw(" ");

  TF1* methodIImu_fit=new TF1("methodIImu_fit",TF1_r0Lambda_of_muratio_II,0.01,40.,0);
  methodIImu_fit->SetLineWidth(1.5);
  methodIImu_fit->SetLineColor(1);
  methodIImu_fit->SetLineStyle(1);
  methodIImu_fit->Draw("SAME");

  TF1* methodIImu_deriv=new TF1("methodIImu_fit",TF1_r0Lambda_of_muratio_II_deriv,0.01,40.,0);
  methodIImu_deriv->SetLineWidth(1.5);
  methodIImu_deriv->SetLineColor(1);
  methodIImu_deriv->SetLineStyle(2);
  methodIImu_deriv->Draw("SAME");

  methodIImucanv->SaveAs("method_II_stability.pdf");
*/
  cout << "\n\n ----- BEGIN METHOD IIP ------ \n\n";

  for(Int_t i=0;i<NDATA;++i){
    methodIIP_r0Lambda[i]=methodII_getr0Lambda(i,PADE);
    methodIIP_r0Lambda_err[i]=methodII_getr0Lambda_err(i,PADE);
    methodIIP_a_over_r0_sqr[i]=TMath::Power(r0_over_a[i],-2.);
    methodIIP_a_over_r0_sqr_err[i]=2.*r0_over_a_err[i]/(r0_over_a[i]*r0_over_a[i]*r0_over_a[i]);
    cout<<methodIIP_a_over_r0_sqr[i]<<"+/-"<<methodIIP_a_over_r0_sqr_err[i]<<"    "<<methodIIP_r0Lambda[i]<<"+/-"<<methodIIP_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;

  TCanvas* methodIIPcanv = new TCanvas("methodIIPcanv","methodIIPcanv",700,700);
  
  TH1F* methodIIPframe=new TH1F("methodIIPframe"," ",1000,200,600);
  methodIIPframe->SetStats(0);
  methodIIPframe->SetMinimum(0.45);
  methodIIPframe->SetMaximum(0.70);
  methodIIPframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIIPframe->GetXaxis()->SetTickLength(0.02);
  methodIIPframe->GetXaxis()->SetLabelSize(0.020);
  methodIIPframe->GetXaxis()->SetLimits(-0.005,0.12); //minx/maxx limits
  methodIIPframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIIPframe->GetYaxis()->SetTickLength(0.02);
  methodIIPframe->GetYaxis()->SetLabelSize(0.020);
  methodIIPframe->Draw(" ");

  TGraphErrors* methodIIPgraph=new TGraphErrors(12,methodIIP_a_over_r0_sqr,methodIIP_r0Lambda,methodIIP_a_over_r0_sqr_err,methodIIP_r0Lambda_err);
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

  methodIIPcanv->SaveAs("methodIIPquenched.pdf");

  cout << "\n\n ----- BEGIN METHOD III ------ \n\n";

  for(Int_t i=0;i<NDATA;++i){
    methodIII_r0Lambda[i]=methodIII_getr0Lambda(i,NOTPADE);
    methodIII_r0Lambda_err[i]=methodIII_getr0Lambda_err(i,NOTPADE);
    methodIII_a_over_r0_sqr[i]=TMath::Power(r0_over_a[i],-2.);
    methodIII_a_over_r0_sqr_err[i]=2.*r0_over_a_err[i]/(r0_over_a[i]*r0_over_a[i]*r0_over_a[i]);
    cout<<methodIII_a_over_r0_sqr[i]<<"+/-"<<methodIII_a_over_r0_sqr_err[i]<<"    "<<methodIII_r0Lambda[i]<<"+/-"<<methodIII_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;

  TCanvas* methodIIIcanv = new TCanvas("methodIIIcanv","methodIIIcanv",700,700);
  
  TH1F* methodIIIframe=new TH1F("methodIIIframe"," ",1000,200,600);
  methodIIIframe->SetStats(0);
  methodIIIframe->SetMinimum(0.45);
  methodIIIframe->SetMaximum(0.70);
  methodIIIframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIIIframe->GetXaxis()->SetTickLength(0.02);
  methodIIIframe->GetXaxis()->SetLabelSize(0.020);
  methodIIIframe->GetXaxis()->SetLimits(-0.005,0.12); //minx/maxx limits
  methodIIIframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIIIframe->GetYaxis()->SetTickLength(0.02);
  methodIIIframe->GetYaxis()->SetLabelSize(0.020);
  methodIIIframe->Draw(" ");

  TGraphErrors* methodIIIgraph=new TGraphErrors(12,methodIII_a_over_r0_sqr,methodIII_r0Lambda,methodIII_a_over_r0_sqr_err,methodIII_r0Lambda_err);
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

  methodIIIcanv->SaveAs("methodIIIquenched.pdf");

  cout << "\n\n ----- BEGIN METHOD IIIP ------ \n\n";

  for(Int_t i=0;i<NDATA;++i){
    methodIIIP_r0Lambda[i]=methodIII_getr0Lambda(i,PADE);
    methodIIIP_r0Lambda_err[i]=methodIII_getr0Lambda_err(i,PADE);
    methodIIIP_a_over_r0_sqr[i]=TMath::Power(r0_over_a[i],-2.);
    methodIIIP_a_over_r0_sqr_err[i]=2.*r0_over_a_err[i]/(r0_over_a[i]*r0_over_a[i]*r0_over_a[i]);
    cout<<methodIIIP_a_over_r0_sqr[i]<<"+/-"<<methodIIIP_a_over_r0_sqr_err[i]<<"    "<<methodIIIP_r0Lambda[i]<<"+/-"<<methodIIIP_r0Lambda_err[i]<<endl;
  }
  cout<<endl<<endl;

  TCanvas* methodIIIPcanv = new TCanvas("methodIIIPcanv","methodIIIPcanv",700,700);
  
  TH1F* methodIIIPframe=new TH1F("methodIIIPframe"," ",1000,200,600);
  methodIIIPframe->SetStats(0);
  methodIIIPframe->SetMinimum(0.45);
  methodIIIPframe->SetMaximum(0.70);
  methodIIIPframe->GetXaxis()->SetTitle("(a/r_{0})^{2}");
  methodIIIPframe->GetXaxis()->SetTickLength(0.02);
  methodIIIPframe->GetXaxis()->SetLabelSize(0.020);
  methodIIIPframe->GetXaxis()->SetLimits(-0.005,0.12); //minx/maxx limits
  methodIIIPframe->GetYaxis()->SetTitle("r_{0} #Lambda");
  methodIIIPframe->GetYaxis()->SetTickLength(0.02);
  methodIIIPframe->GetYaxis()->SetLabelSize(0.020);
  methodIIIPframe->Draw(" ");

  TGraphErrors* methodIIIPgraph=new TGraphErrors(12,methodIIIP_a_over_r0_sqr,methodIIIP_r0Lambda,methodIIIP_a_over_r0_sqr_err,methodIIIP_r0Lambda_err);
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

  methodIIIPcanv->SaveAs("methodIIIPquenched.pdf");

  return 0;

};


