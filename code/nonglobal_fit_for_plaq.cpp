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

#define NDATA 32
#define NBETA 4
#define NIGNORE 9

// * plaquette data from QCDSF database -- get reference *
// r0/a data extrapolated from Table 1 of arxiv:1206.7034
Double_t beta[NDATA]         ={5.20    ,5.20    ,5.20    ,5.20    ,5.20    ,
                               5.25    ,5.25    ,5.25    ,5.25    ,5.25    ,
                               5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,
                               5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,
                               5.40    ,5.40    ,5.40    ,5.40    ,5.40    ,5.40    ,5.40    ,5.40    };
Double_t r0_over_a[NDATA]    ={6.1139  ,6.11291 ,6.11269 ,6.11269 ,6.11271 ,
                               6.59568 ,6.59556 ,6.59498 ,6.59457 ,6.59418 ,
                               7.00588 ,7.00889 ,7.00885 ,7.00885 ,7.00885 ,7.00810 ,7.00810 ,
                               7.00810 ,7.00713 ,7.00664 ,7.00664 ,7.00664 ,7.00629 ,7.00629 ,
                               8.28189 ,8.28442 ,8.28304 ,8.28202 ,8.28073 ,8.28073 ,8.27857 ,8.27857 };
Double_t r0_over_a_err[NDATA]={0.024718,0.018622,0.014121,0.014447,0.017008,
                               0.022982,0.020002,0.015019,0.017242,0.023948,
                               0.052213,0.022759,0.019935,0.019935,0.019935,0.015948,0.015948,
                               0.015948,0.019393,0.023723,0.023723,0.023723,0.027428,0.027428,
                               0.031622,0.026577,0.020004,0.019148,0.020906,0.020906,0.028367,0.028367};

Double_t plaquette[NDATA]    ={0.528994,0.533670,0.536250,0.537070,0.537670,
                               0.538770,0.541150,0.543135,0.544038,0.544782,
                               0.542400,0.545520,0.547293,0.547104,0.547094,0.548569,0.548317,
                               0.548286,0.549187,0.549542,0.549541,0.549540,0.549546,0.549782,
                               0.559000,0.560246,0.561281,0.561550,0.561895,0.561887,0.562278,0.562274};
Double_t plaquette_err[NDATA]={0.000058,0.000040,0.000030,0.000100,0.000030,
                               0.000041,0.000030,0.000015,0.000009,0.000029,
                               0.000050,0.000029,0.000044,0.000044,0.000023,0.000059,0.000030,
                               0.000057,0.000016,0.000009,0.000008,0.000008,0.000012,0.000006,
                               0.000019,0.000010,0.000008,0.000007,0.000012,0.000010,0.000005,0.000005};

//for nf=2 usage
Double_t kappa[NDATA]        ={0.1342  ,0.1350  ,0.1355  ,0.13565 ,0.1358  ,
                               0.1346  ,0.1352  ,0.13575 ,0.1360  ,0.1362  ,
                               0.1340  ,0.1350  ,0.1355  ,0.1355  ,0.1355  ,0.1359  ,0.1359  ,
                               0.1359  ,0.1362  ,0.13632 ,0.13632 ,0.13632 ,0.1364  ,0.1364  ,
                               0.1350  ,0.1356  ,0.1361  ,0.13625 ,0.1364  ,0.1364  ,0.1366  ,0.1366  };
Double_t kappa_c[NDATA]      ={0.136008,0.136008,0.136008,0.136008,0.136008,
                               0.136250,0.136250,0.136250,0.136250,0.136250,
                               0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,
                               0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,
                               0.136690,0.136690,0.136690,0.136690,0.136690,0.136690,0.136690,0.136690};
Double_t kappa_c_err[NDATA]  ={0.000015,0.000015,0.000015,0.000015,0.000015,
                               0.000007,0.000007,0.000007,0.000007,0.000007,
                               0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,
                               0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,
                               0.000022,0.000022,0.000022,0.000022,0.000022,0.000022,0.000022,0.000022};

Double_t index_ignore_list[NIGNORE]={12,13,15,16,19,20,22,28,30};

// -- for chiral limit extrapolation
Double_t quarkmass[NDATA];  //acting as another name for am_q
Double_t quarkmass_err[NDATA];

// -- for continuum limit extrapolation
Double_t beta_extrap[NBETA]={5.20,5.25,5.29,5.40};
Double_t plaquette_extrap[NBETA];
Double_t plaquette_extrap_err[NBETA];

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

// ------------- chiral extrapolation for all methods ---------------

// you're using a chi-squared fit, which means you're ignoring the error
// on the x-variable. Try to find a better way

// *Okay, we're going to implement the "effective variance" approach. Here we use
// \chi^2 = sum_i^N [y_i - f(x_i)]^2/[\sigma_y^2+(\sigma_x*f'(x_i))^2]

// Minuit functions (chi2) p[0]=intercept  p[1]=slope 
void minuit_chiral_plaq520_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  Bool_t ignore_flag;
  for(Int_t i=0;i<=4;++i) { 
    ignore_flag=false;
    for(Int_t j=0;j<NIGNORE;++j){
      if(i==index_ignore_list[j]){
        ignore_flag=true;
      }
    }
    if(ignore_flag){continue;}
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
  Bool_t ignore_flag;
  for(Int_t i=5;i<=9;++i) { 
    ignore_flag=false;
    for(Int_t j=0;j<NIGNORE;++j){
      if(i==index_ignore_list[j]){
        ignore_flag=true;
      }
    }
    if(ignore_flag){continue;}
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
  Bool_t ignore_flag;
  for(Int_t i=10;i<=23;++i) { 
    ignore_flag=false;
    for(Int_t j=0;j<NIGNORE;++j){
      if(i==index_ignore_list[j]){
        ignore_flag=true;
      }
    }
    if(ignore_flag){continue;}
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
  Bool_t ignore_flag;
  for(Int_t i=24;i<=32;++i) { 
    ignore_flag=false;
    for(Int_t j=0;j<NIGNORE;++j){
      if(i==index_ignore_list[j]){
        ignore_flag=true;
      }
    }
    if(ignore_flag){continue;}
      num=TMath::Power( plaquette[i] - fit_parabola(quarkmass[i],p[0],p[1],p[2]) , 2. );
      den=TMath::Power(plaquette_err[i],2.)+TMath::Power(quarkmass_err[i]*fit_parabola_deriv(quarkmass[i],p[0],p[1],p[2]),2.);
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

  Double_t extrap_x_offset=0.00075;

  // ++++++++++++++++++++++++++++ PLAQUETTE +++++++++++++++++++++++++++++++++

  //setting up the canvas and frame for plaquette data
  TCanvas* chiralcanv = new TCanvas("chiralcanv","chiralcanv",700,700);
  TH1F* chiralframe=new TH1F("chiralframe"," ",1000,200,600);
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

  for(Int_t i=0;i<NDATA;++i) { 
    for(Int_t j=0;j<NIGNORE;++j){
      if(i==index_ignore_list[j]){
        plaquette[i]=100.;
      }
    }
  }

  //plotting the data points
  TGraphErrors* chiralgraph=new TGraphErrors(NDATA,quarkmass,plaquette,quarkmass_err,plaquette_err);
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
  Double_t extrap_x_chiral_plaq520[1]={-extrap_x_offset};
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
  Double_t extrap_x_chiral_plaq525[1]={-extrap_x_offset};
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
  Double_t extrap_x_chiral_plaq529[1]={-extrap_x_offset};
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
  Double_t extrap_x_chiral_plaq540[1]={-extrap_x_offset};
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

  chiralcanv->SaveAs("Plaquette_Chiral_Extrap_new.pdf");

  cout << endl;

  cout << "plaquette_extrap \n";
  for(Int_t i=0;i<NBETA;++i){
    cout << plaquette_extrap[i] << " ," << flush;
    if( (i==4) || (i==9) || (i==16) || (i==23) ) { cout << endl; }
  }

  cout << "\n\nplaquette_extrap_err \n";
  for(Int_t i=0;i<NBETA;++i){
    cout << plaquette_extrap_err[i] << " ," << flush;
    if( (i==4) || (i==9) || (i==16) || (i==23) ) { cout << endl; }
  }

  cout << endl << endl;

  return 0;

};


