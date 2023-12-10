//C++ standard headers
#include <iostream>
#include <climits>
#include <vector>

//ROOT headers for math
#include "TMath.h"
#include "TComplex.h"

//ROOT headers for plotting
#include "TF1.h"

//in-house header files
#include "BetaFn.h"

Double_t Lambda0=-1.;  Double_t Lambda0_err=-1.;
Double_t Lambda1=-1.;  Double_t Lambda1_err=-1.;
Double_t Lambda2=-1.;  Double_t Lambda2_err=-1.;
Double_t Lambda3=-1.;  Double_t Lambda3_err=-1.;
Double_t Lambda4=-1.;  Double_t Lambda4_err=-1.;
Double_t Lambda5=-1.;  Double_t Lambda5_err=-1.;
Double_t Lambda6=-1.;  Double_t Lambda6_err=-1.;

//all masses in MeV
Double_t mass_up=2.;   Double_t mass_charm=1275.;  Double_t mass_top=175000.;
Double_t mass_down=2.; Double_t mass_strange=95.;  Double_t mass_bottom=4180.;
Double_t mass_z=91187.6;
Double_t mass_tau=1777.0;

//all scales in MeV
Double_t mu_up=mass_up;   Double_t mu_charm=mass_charm;  Double_t mu_top=mass_top;
Double_t mu_down=mass_down; Double_t mu_strange=mass_strange;  Double_t mu_bottom=mass_bottom;
Double_t mu_z=91187.6;

Double_t mass_up_err=-1.;  Double_t mass_charm_err=25.; Double_t mass_top_err=-1.;
Double_t mass_down_err=-1.;  Double_t mass_strange_err=-1.; Double_t mass_bottom_err=30.;

using namespace std;

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
Double_t TF1_mu_over_lambda(Double_t* x, Double_t* par){
  BetaFn bf( (Int_t)par[0] , (Int_t)par[1], (Int_t)par[2] );
  bf[2]=par[3];  bf[3]=par[4];
  return mu_over_lambda(x[0],bf);
};
Double_t alpha_of_mu(Double_t mu_over_lambda, BetaFn bf){
  TF1 scratch("scratch",TF1_mu_over_lambda,0.0001,1000,5);
  scratch.SetParameter(0,bf.nloops());
  scratch.SetParameter(1,bf.nflav());
  scratch.SetParameter(2,bf.variable());
  scratch.SetParameter(3,bf[2]);
  scratch.SetParameter(4,bf[3]);
  Double_t x=scratch.GetX(mu_over_lambda,0.001,10);
  if(bf.variable()==2){return ( x*x/(4.*TMath::Pi()) );}
  if(bf.variable()==0){return 4.*x*TMath::Pi();}
  return x;
};

//from http://arxiv.org/pdf/hep-ph/0004189v1.pdf 
// this function is for matching down from alpha(nf) to alpha(nf-1) 
//here we assume that the matching is done at threshold, ie \mu=m_h where m_h is the msbar threshold mass
Double_t zeta_g_sqr(Double_t alpha_nf,Int_t nf){
  Double_t coeff0=1.;
  Double_t coeff1=0.;
  Double_t coeff2=11./72.;
  Double_t coeff3=564731./124416.-82043.*RZ3/27648.-2633.*(nf-1.)/31104.;
  return ( coeff0*TMath::Power( alpha_nf/TMath::Pi() , 0. )
         + coeff1*TMath::Power( alpha_nf/TMath::Pi() , 1. )
         + coeff2*TMath::Power( alpha_nf/TMath::Pi() , 2. )
         + coeff3*TMath::Power( alpha_nf/TMath::Pi() , 3. )  );
};

//from http://arxiv.org/pdf/hep-ph/0004189v1.pdf 
// this function is for matching up from alpha(nf-1) to alpha(nf) 
//here we assume that the matching is done at threshold, ie \mu=m_h where m_h is the msbar threshold mass
Double_t recip_zeta_g_sqr(Double_t alpha_nl,Int_t nl){
  Double_t coeff0=1.;
  Double_t coeff1=0.;
  Double_t coeff2=-11./72.;
  Double_t coeff3=-564731./124416.+82043.*RZ3/27648.+2633.*(nl)/31104.;
  return ( coeff0*TMath::Power( alpha_nl/TMath::Pi() , 0. )
         + coeff1*TMath::Power( alpha_nl/TMath::Pi() , 1. )
         + coeff2*TMath::Power( alpha_nl/TMath::Pi() , 2. )
         + coeff3*TMath::Power( alpha_nl/TMath::Pi() , 3. )  );
};


// -------------------------------------------------------------------------------------------------------------- //
//                                                                                                                //
//                                                  MAIN BEGIN                                                    //
//                                                                                                                //
// -------------------------------------------------------------------------------------------------------------- //

Int_t main(){

  // 4-loop beta functions with alpha as the variable
  BetaFn bf_ms_0(4,0,1);
  BetaFn bf_ms_1(4,1,1);
  BetaFn bf_ms_2(4,2,1);
  BetaFn bf_ms_3(4,3,1);
  BetaFn bf_ms_4(4,4,1);
  BetaFn bf_ms_5(4,5,1);
  BetaFn bf_ms_6(4,6,1);

  Lambda3=307+17;     //testing what would happen if I use r0=0.472
  Lambda3_err=17.;
//  Lambda3=286.;     //in MeV 
//  Lambda3_err=16.; 
  Double_t alpha3_mc=alpha_of_mu(mu_charm/Lambda3,bf_ms_3);
  Double_t alpha4_mc=alpha3_mc*recip_zeta_g_sqr(alpha3_mc,3);

  Lambda4=mu_charm/mu_over_lambda(alpha4_mc,bf_ms_4);
  Double_t alpha4_mtau=alpha_of_mu(mass_tau/Lambda4,bf_ms_4);
  Double_t alpha4_mb=alpha_of_mu(mu_bottom/Lambda4,bf_ms_4);
  Double_t alpha5_mb=alpha4_mb*recip_zeta_g_sqr(alpha4_mb,4);

  Lambda5=mu_bottom/mu_over_lambda(alpha5_mb,bf_ms_5);
  Double_t alpha5_mz=alpha_of_mu(mu_z/Lambda5,bf_ms_5);

  cout <<endl << "alpha(mu=mc,nf=3) = " << alpha3_mc << endl;
  cout << "alpha(mu=mc,nf=4) = " << alpha4_mc << endl;
  cout << "alpha(mu=mb,nf=4) = " << alpha4_mb << endl;
  cout << "alpha(mu=mb,nf=5) = " << alpha5_mb << endl;
  cout << "alpha(mu=mz,nf=5) = " << alpha5_mz << endl << endl;
  cout << "Lambda(nf=3) = " << Lambda3 << endl;
  cout << "Lambda(nf=4) = " << Lambda4 << endl;
  cout << "Lambda(nf=5) = " << Lambda5 << endl;

  cout << "alpha(mu=mtau,nf=4) = " << alpha4_mtau << endl;

  exit(1);
  
  cout << endl << endl << endl;
  Double_t _Lambda3,_Lambda4,_Lambda5,_mass_charm,_mu_charm,_mass_bottom,_mu_bottom;
  // finding the systematic spread
  for(Double_t LL=(Lambda3-1.*Lambda3_err); LL<(1.+Lambda3+1.*Lambda3_err); LL=(LL+Lambda3_err)){
    _Lambda3=LL;
    for(Double_t MC=(mass_charm-1.*mass_charm_err); MC<(1.+mass_charm+1.*mass_charm_err); MC=(MC+mass_charm_err)){
      _mass_charm=MC;  _mu_charm=MC;
      for(Double_t MB=(mass_bottom-1.*mass_bottom_err); MB<(1.+mass_bottom+1.*mass_bottom_err); MB=(MB+mass_bottom_err)){
        _mass_bottom=MB;  _mu_bottom=MB;

        alpha3_mc=alpha_of_mu(_mu_charm/_Lambda3,bf_ms_3);
        alpha4_mc=alpha3_mc*recip_zeta_g_sqr(alpha3_mc,3);

        _Lambda4=_mu_charm/mu_over_lambda(alpha4_mc,bf_ms_4);
        alpha4_mb=alpha_of_mu(_mu_bottom/_Lambda4,bf_ms_4);
        alpha5_mb=alpha4_mb*recip_zeta_g_sqr(alpha4_mb,4);

        _Lambda5=_mu_bottom/mu_over_lambda(alpha5_mb,bf_ms_5);
        alpha5_mz=alpha_of_mu(mu_z/_Lambda5,bf_ms_5);

        cout << "********* " << LL << "  " << MC << "  " << MB << " ***************\n";
        cout <<  endl << "alpha(mu=mc,nf=3) = " << alpha3_mc << endl;
        cout << "alpha(mu=mc,nf=4) = " << alpha4_mc << endl;
        cout << "alpha(mu=mb,nf=4) = " << alpha4_mb << endl;
        cout << "alpha(mu=mb,nf=5) = " << alpha5_mb << endl;
        cout << "alpha(mu=mz,nf=5) = " << alpha5_mz << endl << endl;
        cout << "Lambda(nf=3) = " << _Lambda3 << endl;
        cout << "Lambda(nf=4) = " << _Lambda4 << endl;
        cout << "Lambda(nf=5) = " << _Lambda5 << endl;
      }
    }
  }

  return 0;

};


