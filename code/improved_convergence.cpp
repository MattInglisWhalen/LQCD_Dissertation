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

#define NDATA 11
#define NOTPADE 0
#define PADE 1

Double_t RATIO=1.00; //for testing the sensitivity to scale

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

//relating g_msbar to g0 and g_boost
Double_t u0(Int_t ind){
  return TMath::Power(plaquette[ind],0.25);
};
Double_t u0_err(Int_t ind){
  Double_t partial_plaquette=0.25*TMath::Power(plaquette[ind],-0.75);
  Double_t term1=TMath::Power(partial_plaquette,2.)*TMath::Power(plaquette_err[ind],2.);
  return TMath::Sqrt(term1);
};

Double_t c_sw(Int_t ind,Double_t _beta=1.){
  Double_t gsqr=6./_beta;
  if(ind>=0){gsqr=6./beta[ind];}
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
  if(ind<0){return 0.;}
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

Double_t t1_lat(Int_t ind,Double_t _beta=1.){
  Double_t csw=c_sw(ind,_beta);
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

Double_t t2_lat(Int_t ind,Double_t _beta=1.){
  Double_t csw=c_sw(ind,_beta);
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

Double_t t1_boost(Int_t ind,Double_t val=0.){
  Double_t cswb;
  if(ind<0){cswb=val;}
  else{cswb=c_sw_boosted(ind);}
  return ( 0.1348680 - NFLAV*(0.0066960-0.0050467*cswb+0.0298435*cswb*cswb ) );
};
Double_t t1_boost_err(Int_t ind,Double_t val=0.){
  Double_t cswb;
  if(ind<0){cswb=val; return 0.;}
  else{cswb=c_sw_boosted(ind);}
  Double_t partial_c_sw_boosted=-NFLAV*(-0.0050467+2.*0.0298435*cswb);
  Double_t term1=TMath::Power(partial_c_sw_boosted,2.)*TMath::Power(c_sw_boosted_err(ind),2.);
  return TMath::Sqrt(term1);
};

Double_t t2_boost(Int_t ind,Double_t val=0.){
  Double_t cswb;
  if(ind<0){cswb=val;}
  else{cswb=c_sw_boosted(ind);}
  return ( 0.0217565 - NFLAV*(0.000753-0.001053*cswb+0.000498*cswb*cswb-0.000474*cswb*cswb*cswb-0.000104*cswb*cswb*cswb*cswb) );
};
Double_t t2_boost_err(Int_t ind,Double_t val=0.){
  Double_t cswb;
  if(ind<0){cswb=val; return 0.;}
  else{cswb=c_sw_boosted(ind);}
  Double_t partial_c_sw_boosted=-NFLAV*(0.001053-2.*0.000498*cswb-3.*0.000474*cswb*cswb-4.*0.000104*cswb*cswb*cswb);
  Double_t term1=TMath::Power(partial_c_sw_boosted,2.)*TMath::Power(c_sw_boosted_err(ind),2.);
  return TMath::Sqrt(term1);
};

Double_t g_lat(Int_t ind,Double_t _beta=1.){
  if(ind<0){return TMath::Sqrt( 6./_beta );}
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

Double_t g_ms_I_lat_direct(Int_t ind,Double_t _beta=1.){
  BetaFn bf(4,NFLAV,2);
  return TMath::Sqrt( TMath::Abs( TMath::Power( g_lat(ind,_beta) , 2. ) - (2.*bf[0]*TMath::Log(RATIO))*TMath::Power( g_lat(ind,_beta) , 4. )
                        + (4.*bf[0]*bf[0]*TMath::Log(RATIO)*TMath::Log(RATIO)-2.*bf[1]*TMath::Log(RATIO)+t2_lat(ind,_beta)-bf[1]*t1_lat(ind,_beta)/bf[0])*TMath::Power(g_lat(ind,_beta),6.) ) );
};
Double_t g_ms_I_lat_recip(Int_t ind,Double_t _beta=1.){
  BetaFn bf(4,NFLAV,2);
  return TMath::Power( TMath::Abs( TMath::Power( g_lat(ind,_beta) , -2. ) + (2.*bf[0]*TMath::Log(RATIO))
                        + (2.*bf[1]*TMath::Log(RATIO)-t2_lat(ind,_beta)+bf[1]*t1_lat(ind,_beta)/bf[0])*TMath::Power(g_lat(ind,_beta),2.) ) , -0.5 );
};

// x[0] is g in the lattice scheme
Double_t TF1_g_ms_I_lat_direct(Double_t* x, Double_t* par){
  return g_ms_I_lat_direct(-1,6./(x[0]*x[0]));
}
Double_t TF1_g_ms_I_lat_recip(Double_t* x, Double_t* par){
  return g_ms_I_lat_recip(-1,6./(x[0]*x[0]));
}
Double_t TF1_rat(Double_t* x, Double_t* par){
  return TMath::Abs( 1-(g_ms_I_lat_recip(-1,6./(x[0]*x[0])))/(g_ms_I_lat_direct(-1,6./(x[0]*x[0]))) );
};

// -------------------------- method I ---------------------
Double_t g_ms_I(Int_t ind,Int_t opt=0){
// opt==0 -- use direct expansion of g_MS ie g_MS^2=g_boost^2+d_1*g_boost^4+...
// opt==1 -- use reciprocal expansion of g_MS ie 1/g_MS^2=1/g_boost^2+...  -- this one agrees with the r0Lambda values in the paper
  BetaFn bf(4,NFLAV,2);
  if(opt>0){
    return TMath::Power( TMath::Power( g_boost(ind) , -2. ) + (2.*bf[0]*TMath::Log(RATIO))
                        + (2.*bf[1]*TMath::Log(RATIO)-t2_boost(ind)+bf[1]*t1_boost(ind)/bf[0])*TMath::Power(g_boost(ind),2.) , -0.5 );
  }
  return TMath::Sqrt( TMath::Power( g_boost(ind) , 2. ) - (2.*bf[0]*TMath::Log(RATIO))*TMath::Power( g_boost(ind) , 4. )
                        + (4.*bf[0]*bf[0]*TMath::Log(RATIO)*TMath::Log(RATIO)-2.*bf[1]*TMath::Log(RATIO)+t2_boost(ind)-bf[1]*t1_boost(ind)/bf[0])*TMath::Power(g_boost(ind),6.) );
};

// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//                                                  MAIN BEGIN                                                    //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

Int_t main(){

  Double_t g_lat_data[NDATA];
  Double_t g_ms_boost_direct_data[NDATA];
  Double_t g_ms_boost_recip_data[NDATA];
  Double_t g_ms_boost_rat_data[NDATA];

  for(Int_t i=0;i<NDATA;++i){
    g_lat_data[i]=TMath::Sqrt(6./beta[i]);
    g_ms_boost_direct_data[i]=g_ms_I(i,0);
    g_ms_boost_recip_data[i]=g_ms_I(i,1);
    g_ms_boost_rat_data[i]=TMath::Abs(1.-g_ms_I(i,1)/g_ms_I(i,0));
    cout << g_lat_data[i] << " : " << g_ms_boost_direct_data[i] << "  " << g_ms_boost_recip_data[i] << "  " << g_ms_boost_rat_data[i] << endl;
  }

  TCanvas* canv = new TCanvas("canv","canv",700,700);
  canv->SetLogy();

  TH1F* frame=new TH1F("frame","g_{#bar{MS}} vs g_{0}",1000,200,600);
  frame->SetStats(0);
  frame->SetMinimum(0.0001);
  frame->SetMaximum(10.);
  frame->GetXaxis()->SetTitle("g_{0}");
  frame->GetXaxis()->SetTickLength(0.02);
  frame->GetXaxis()->SetLabelSize(0.020);
  frame->GetXaxis()->SetLimits(0.0,2.4); //minx/maxx limits
  frame->GetYaxis()->SetTitle("g_{#bar{MS}}");
  frame->GetYaxis()->SetTickLength(0.02);
  frame->GetYaxis()->SetLabelSize(0.020);
  frame->Draw(" ");

  TF1* func_direct=new TF1("func_recip",TF1_g_ms_I_lat_direct,0.,2.4);
  func_direct->SetLineWidth(1.5);
  func_direct->SetLineColor(1);
  func_direct->SetLineStyle(1);
  func_direct->Draw("SAME");

  TF1* func_recip=new TF1("func_direct",TF1_g_ms_I_lat_recip,0.,2.4);
  func_recip->SetLineWidth(1.5);
  func_recip->SetLineColor(2);
  func_recip->SetLineStyle(2);
  func_recip->Draw("SAME");

  TF1* func_rat=new TF1("func_direct",TF1_rat,0.,2.4);
  func_rat->SetLineWidth(1.5);
  func_rat->SetLineColor(4);
  func_rat->SetLineStyle(3);
  func_rat->Draw("SAME");

  TGraph* gr_direct=new TGraph(NDATA,g_lat_data,g_ms_boost_direct_data);
  gr_direct->SetMarkerColor(1);
  gr_direct->SetMarkerStyle(33);
  gr_direct->SetMarkerSize(1);
  gr_direct->Draw("P");

  TGraph* gr_recip=new TGraph(NDATA,g_lat_data,g_ms_boost_recip_data);
  gr_recip->SetMarkerColor(2);
  gr_recip->SetMarkerStyle(33);
  gr_recip->SetMarkerSize(0.5);
  gr_recip->Draw("P");

  TGraph* gr_rat=new TGraph(NDATA,g_lat_data,g_ms_boost_rat_data);
  gr_rat->SetMarkerColor(4);
  gr_rat->SetMarkerStyle(33);
  gr_rat->SetMarkerSize(1);
  gr_rat->Draw("P");

  canv->SaveAs("comparison.pdf");

// I want to show that using the boosted coupling actually improves convergence
//  cout<<g_ms_I_lat(2,0)<<endl;
//  cout<<g_ms_I_lat(2,1)<<endl;
//  cout<<g_ms_I(2,0)<<endl;
//  cout<<g_ms_I(2,1)<<endl;

  return 0;

};


