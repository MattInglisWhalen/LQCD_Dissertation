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

using namespace std;

#define NFLAV 2

#define NOTPADE 0
#define PADE 1

#define DIRECT 0
#define RECIPROCAL 1

Double_t RATIO=1.00; //for testing the sensitivity to scale

#define NDATAR0A 27
#define NFITPLQ 32

//#define NIGNORE 1
#define NIGNORE 9
//#define NIGNORE 24

//data from Table 1 of arxiv:1206.7034
Double_t beta[NDATAR0A]          ={5.25    ,5.25    ,5.25    ,5.25    ,5.25    ,
                                   5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,
                                   5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,
                                   5.40    ,5.40    ,5.40    ,5.40    ,5.40    ,5.40    ,5.40    ,5.40    };

Double_t r0mpi[NDATAR0A]         ={3.256   ,2.523   ,1.687   ,1.215   ,0.658   ,
                                   4.039   ,2.946   ,2.525   ,2.329   ,2.290   ,2.360   ,1.763   ,
                                   1.677   ,1.087   ,0.779   ,0.750   ,0.735   ,0.463   ,0.399   ,
                                   3.339   ,2.588   ,1.829   ,1.576   ,1.274   ,1.246   ,0.700   ,0.660   };
Double_t r0mpi_err[NDATAR0A]     ={0.027   ,0.022   ,0.014   ,0.011   ,0.009   ,
                                   0.032   ,0.024   ,0.030   ,0.020   ,0.018   ,0.047   ,0.017   ,
                                   0.013   ,0.010   ,0.009   ,0.007   ,0.006   ,0.007   ,0.006   ,
                                   0.030   ,0.024   ,0.017   ,0.015   ,0.014   ,0.012   ,0.008   ,0.007   };
Double_t ampi[NDATAR0A]          ={0.4932  ,0.3821  ,0.2556  ,0.1840  ,0.0997  ,
                                   0.5767  ,0.4206  ,0.3605  ,0.3325  ,0.3270  ,0.3369  ,0.2518  ,
                                   0.2395  ,0.1552  ,0.1112  ,0.1070  ,0.1050  ,0.0660  ,0.0570  ,
                                   0.4030  ,0.3123  ,0.2208  ,0.1902  ,0.1538  ,0.1505  ,0.0845  ,0.0797  };
Double_t ampi_err[NDATAR0A]      ={0.0010  ,0.0013  ,0.0005  ,0.0007  ,0.0011  ,
                                   0.0011  ,0.0009  ,0.0032  ,0.0014  ,0.0006  ,0.0062  ,0.0015  ,
                                   0.0005  ,0.0006  ,0.0009  ,0.0005  ,0.0003  ,0.0008  ,0.0007  ,
                                   0.0004  ,0.0007  ,0.0007  ,0.0006  ,0.0010  ,0.0005  ,0.0006  ,0.0003  }; 

Double_t kappa[NDATAR0A]         ={0.1346  ,0.1352  ,0.13575 ,0.1360  ,0.1362  ,
                                   0.1340  ,0.1350  ,0.1355  ,0.1355  ,0.1355  ,0.1359  ,0.1359  ,
                                   0.1359  ,0.1362  ,0.13632 ,0.13632 ,0.13632 ,0.1364  ,0.1364  ,
                                   0.1350  ,0.1356  ,0.1361  ,0.13625 ,0.1364  ,0.1364  ,0.1366  ,0.1366  };
Double_t kappa_c[NDATAR0A]       ={0.136250,0.136250,0.136250,0.136250,0.136250,
                                   0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,
                                   0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,
                                   0.136690,0.136690,0.136690,0.136690,0.136690,0.136690,0.136690,0.136690};
Double_t kappa_c_err[NDATAR0A]   ={0.000007,0.000007,0.000007,0.000007,0.000007,
                                   0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,
                                   0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,
                                   0.000022,0.000022,0.000022,0.000022,0.000022,0.000022,0.000022,0.000022};

//Int_t index_ignore_list[NIGNORE]={100};
Int_t index_ignore_list[NIGNORE]={7,8,10,11,14,15,17,23,25};
//Int_t index_ignore_list[NIGNORE]={1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21,22,23,24,25,26};

// -- must be calculated first --
Double_t r0_over_a[NDATAR0A];
Double_t r0_over_a_err[NDATAR0A];
Double_t quarkmass[NDATAR0A];  //acting as another name for am_q
Double_t quarkmass_err[NDATAR0A];

// -- placeholders for the variables I want to calculate --
// * plaquette data from QCDSF database -- get reference *
Double_t beta_extrap[NFITPLQ]={5.20    ,5.20    ,5.20    ,5.20    ,5.20    ,
                               5.25    ,5.25    ,5.25    ,5.25    ,5.25    ,
                               5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,
                               5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,5.29    ,
                               5.40    ,5.40    ,5.40    ,5.40    ,5.40    ,5.40    ,5.40    ,5.40    };

Double_t kappa_extrap[NFITPLQ]          ={0.1342  ,0.1350  ,0.1355  ,0.13565 ,0.1358  ,
                                        0.1346  ,0.1352  ,0.13575 ,0.1360  ,0.1362  ,
                                        0.1340  ,0.1350  ,0.1355  ,0.1355  ,0.1355  ,0.1359  ,0.1359  ,
                                        0.1359  ,0.1362  ,0.13632 ,0.13632 ,0.13632 ,0.1364  ,0.1364  ,
                                        0.1350  ,0.1356  ,0.1361  ,0.13625 ,0.1364  ,0.1364  ,0.1366  ,0.1366  };
Double_t kappa_c_extrap[NFITPLQ]        ={0.136008,0.136008,0.136008,0.136008,0.136008,
                                        0.136250,0.136250,0.136250,0.136250,0.136250,
                                        0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,
                                        0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,0.136410,
                                        0.136690,0.136690,0.136690,0.136690,0.136690,0.136690,0.136690,0.136690};
Double_t kappa_c_err_extrap[NFITPLQ]    ={0.000015,0.000015,0.000015,0.000015,0.000015,
                                        0.000007,0.000007,0.000007,0.000007,0.000007,
                                        0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,
                                        0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,0.000009,
                                        0.000022,0.000022,0.000022,0.000022,0.000022,0.000022,0.000022,0.000022};

Double_t r0_over_a_extrap[NFITPLQ];
Double_t r0_over_a_err_extrap[NFITPLQ];

Double_t calc_r0_over_a(Int_t ind){
  return r0mpi[ind]/ampi[ind];
};
Double_t calc_r0_over_a_err(Int_t ind){
  Double_t partial_r0mpi=calc_r0_over_a(ind)/r0mpi[ind];
  Double_t partial_ampi=-calc_r0_over_a(ind)/ampi[ind];

  Double_t corr_r0mpi__ampi=0.;      // no correlations

  Double_t term1=TMath::Power(partial_r0mpi*r0mpi_err[ind],2.);
  Double_t term2=TMath::Power(partial_ampi*ampi_err[ind],2.);
  Double_t termcov=2.*partial_r0mpi*partial_ampi*r0mpi_err[ind]*ampi_err[ind]*corr_r0mpi__ampi;

  return TMath::Sqrt( term1+term2+termcov );
}

Double_t amq(Int_t ind,Int_t is_extrap=0){
  if(is_extrap==1){  
    if( (ind<0) || (ind>(NFITPLQ-1)) ){cout<<"\n c_sw: Bad index ["<<ind<<"]\n\n"; exit(-1);}
    return 0.5*( (1./kappa_extrap[ind])-(1./kappa_c_extrap[ind]) );
  }
  if( (ind<0) || (ind>(NDATAR0A-1)) ){cout<<"\n c_sw: Bad index ["<<ind<<"]\n\n"; exit(-1);}

  return 0.5*( (1./kappa[ind])-(1./kappa_c[ind]) );
};
Double_t amq_err(Int_t ind,Int_t is_extrap=0){  
  if(is_extrap==1){  
    if( (ind<0) || (ind>(NFITPLQ-1)) ){cout<<"\n u0: Bad index ["<<ind<<"]\n\n"; exit(-1);}
    Double_t partial_kappa_c=0.5*TMath::Power( kappa_c_extrap[ind] , -2. );
    Double_t term1=TMath::Power(partial_kappa_c,2.)*TMath::Power(kappa_c_err_extrap[ind],2.);
    return TMath::Sqrt(term1);
  }
  if( (ind<0) || (ind>(NDATAR0A-1)) ){cout<<"\n u0: Bad index ["<<ind<<"]\n\n"; exit(-1);}
  Double_t partial_kappa_c=0.5*TMath::Power( kappa_c[ind] , -2. );
  Double_t term1=TMath::Power(partial_kappa_c,2.)*TMath::Power(kappa_c_err[ind],2.);
  return TMath::Sqrt(term1);
};

// ------------- Functions for global fitting ------------
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
Double_t fit_global_err(Double_t _amq,Double_t _amq_err,Double_t beta, Double_t* par, TMatrixDSym cov){

  Double_t partial_amq=fit_global_deriv(_amq,beta,par);

  Double_t partial_par[8];  // for error on the parameters
  partial_par[0]=fit_global(_amq,beta,par);
  partial_par[1]=fit_global(_amq,beta,par)*beta;
  partial_par[2]=fit_global(_amq,beta,par)*_amq;
  partial_par[3]=fit_global(_amq,beta,par)*beta*_amq;
  partial_par[4]=fit_global(_amq,beta,par)*beta*beta*_amq;
  partial_par[5]=fit_global(_amq,beta,par)*_amq*_amq;
  partial_par[6]=fit_global(_amq,beta,par)*beta*_amq*_amq;
  partial_par[7]=fit_global(_amq,beta,par)*beta*beta*_amq*_amq;

  Double_t error_sqr=TMath::Power(partial_amq,2.)*TMath::Power(_amq_err,2.);

  for(Int_t i=0;i<8;++i){
    for(Int_t j=0;j<8;++j){
      error_sqr+=partial_par[i]*partial_par[j]*cov(i,j);
    }
  }

  return TMath::Sqrt( error_sqr );
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
void minuit_chiral_r0a_global_chisqr(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t chi2=0.0;
  Double_t num,den;
  Bool_t ignore_flag;
  for(Int_t i=0;i<NDATAR0A;++i) { 
    ignore_flag=false;
    for(Int_t j=0;j<NIGNORE;++j){
      if(i==index_ignore_list[j]){
        ignore_flag=true;
      }
    }
    if(ignore_flag){continue;}
    num=TMath::Power( r0_over_a[i] - fit_global(quarkmass[i],beta[i],p) , 2. );
    den=TMath::Power(r0_over_a_err[i],2.)+TMath::Power(quarkmass_err[i]*fit_global_deriv(quarkmass[i],beta[i],p),2.);
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

  cout << "              amq              |            r0/a         |          beta         \n";
  for(Int_t i=0;i<NDATAR0A;++i){
    r0_over_a[i]=calc_r0_over_a(i);
    r0_over_a_err[i]=calc_r0_over_a_err(i);
    quarkmass[i]=amq(i);
    quarkmass_err[i]=amq_err(i);
    cout << i << "   " 
         << quarkmass[i]<<"+/-"<<quarkmass_err[i]<<"       "
         << r0_over_a[i]<<"+/-"<<r0_over_a_err[i]<<"             "
         << beta[i]<<"+/-"<<0<<"       "
         <<endl;
  }

  // ---------------------------------------------------------------- //
  //                                                                  //
  //                     Chiral Limit Extrapolation                   //
  //                                                                  //
  // ---------------------------------------------------------------- //

  // +++++++++++++++++++++++++++++++++ R0 OVER A  +++++++++++++++++++++++++++++++++++

  // ------------ setting up the canvas and frame for plaquette data  ----------
  TCanvas* chiralr0canv = new TCanvas("chiralr0canv","chiralr0canv",700,700);
  TH1F* chiralr0frame=new TH1F("chiralr0frame"," ",1000,200,600);
  chiralr0frame->SetStats(0);
  chiralr0frame->SetMinimum(6.);
  chiralr0frame->SetMaximum(9.);
  chiralr0frame->GetXaxis()->SetTitle("am_{q}");
  chiralr0frame->GetXaxis()->SetTickLength(0.02);
  chiralr0frame->GetXaxis()->SetLabelSize(0.020);
  chiralr0frame->GetXaxis()->SetLimits(-0.005,0.08); //minx/maxx limits
  chiralr0frame->GetYaxis()->SetTitle("r_{0}/a");
  chiralr0frame->GetYaxis()->SetTickLength(0.02);
  chiralr0frame->GetYaxis()->SetLabelSize(0.020);
  chiralr0frame->Draw(" ");

  for(Int_t i=0;i<NDATAR0A;++i) { 
    for(Int_t j=0;j<NIGNORE;++j){
      if(i==index_ignore_list[j]){
        r0_over_a[i]=100.;
      }
    }
  }

  //plotting the data points
  TGraphErrors* chiralr0graph=new TGraphErrors(NDATAR0A,quarkmass,r0_over_a,quarkmass_err,r0_over_a_err);
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
  minuit_chiral_r0_global->SetParameter(0,"A00",-6.0, 0.01, -7., -5.);
  minuit_chiral_r0_global->SetParameter(1,"A01", 1.5, 0.01,  0., 2.);
  minuit_chiral_r0_global->SetParameter(2,"A10", 0.0, 0.01, -3., 3.);
  minuit_chiral_r0_global->SetParameter(3,"A11", 0.0, 0.01, -2., 2.);
  minuit_chiral_r0_global->SetParameter(4,"A12", 0.0, 0.01, -1., 1.);
  minuit_chiral_r0_global->SetParameter(5,"A20", 0.0, 0.1, -40., 40.);
  minuit_chiral_r0_global->SetParameter(6,"A21", 0.0, 0.01, -5., 5.);
  minuit_chiral_r0_global->SetParameter(7,"A22", 0.0, 0.01, -3., 3.);

//  minuit_chiral_r0_global->FixParameter(0);
//  minuit_chiral_r0_global->FixParameter(1);
  minuit_chiral_r0_global->FixParameter(2);
  minuit_chiral_r0_global->FixParameter(3);
  minuit_chiral_r0_global->FixParameter(4);
  minuit_chiral_r0_global->FixParameter(5);
  minuit_chiral_r0_global->FixParameter(6);
  minuit_chiral_r0_global->FixParameter(7);

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

  Double_t offset=0.00075;

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
  method_chiral_r0_525_fit->SetParameter(8,5.25);

  method_chiral_r0_525_fit->SetLineWidth(1.5);
  method_chiral_r0_525_fit->SetLineColor(1);
  method_chiral_r0_525_fit->SetLineStyle(2);
  method_chiral_r0_525_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_r0_525[1]={-offset};
  Double_t extrap_y_chiral_r0_525[1]={fit_global(0,5.25,pval_chiral_r0_global)};
  Double_t extrap_x_err_chiral_r0_525[1]={0.};
  Double_t extrap_y_err_chiral_r0_525[1]={fit_global(0,beta_extrap[1],pval_chiral_r0_global)*TMath::Sqrt(
     TMath::Power(perr_chiral_r0_global[0],2.)+TMath::Power(beta_extrap[1]*perr_chiral_r0_global[1],2.)+2*beta_extrap[1]*cov(0,1) )};
  TGraphErrors* extrap_chiral_r0_525=new TGraphErrors(1,extrap_x_chiral_r0_525,extrap_y_chiral_r0_525,extrap_x_err_chiral_r0_525,extrap_y_err_chiral_r0_525);
  extrap_chiral_r0_525->SetMarkerColor(1);
  extrap_chiral_r0_525->SetMarkerStyle(1);
  extrap_chiral_r0_525->SetMarkerSize(1);
  extrap_chiral_r0_525->Draw("P");

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
  method_chiral_r0_529_fit->SetParameter(8,5.29);

  method_chiral_r0_529_fit->SetLineWidth(1.5);
  method_chiral_r0_529_fit->SetLineColor(1);
  method_chiral_r0_529_fit->SetLineStyle(2);
  method_chiral_r0_529_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_r0_529[1]={-offset};
  Double_t extrap_y_chiral_r0_529[1]={fit_global(0,5.29,pval_chiral_r0_global)};
  Double_t extrap_x_err_chiral_r0_529[1]={0.};
  Double_t extrap_y_err_chiral_r0_529[1]={fit_global(0,beta_extrap[2],pval_chiral_r0_global)*TMath::Sqrt(
      TMath::Power(perr_chiral_r0_global[0],2.)+TMath::Power(beta_extrap[2]*perr_chiral_r0_global[1],2.)+2*beta_extrap[2]*cov(0,1) )};
  TGraphErrors* extrap_chiral_r0_529=new TGraphErrors(1,extrap_x_chiral_r0_529,extrap_y_chiral_r0_529,extrap_x_err_chiral_r0_529,extrap_y_err_chiral_r0_529);
  extrap_chiral_r0_529->SetMarkerColor(1);
  extrap_chiral_r0_529->SetMarkerStyle(1);
  extrap_chiral_r0_529->SetMarkerSize(1);
  extrap_chiral_r0_529->Draw("P");

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
  method_chiral_r0_540_fit->SetParameter(8,5.40);

  method_chiral_r0_540_fit->SetLineWidth(1.5);
  method_chiral_r0_540_fit->SetLineColor(1);
  method_chiral_r0_540_fit->SetLineStyle(2);
  method_chiral_r0_540_fit->Draw("SAME");

  //putting in the extrapolated value
  Double_t extrap_x_chiral_r0_540[1]={-offset};
  Double_t extrap_y_chiral_r0_540[1]={fit_global(0,5.40,pval_chiral_r0_global)};
  Double_t extrap_x_err_chiral_r0_540[1]={0.};
  Double_t extrap_y_err_chiral_r0_540[1]={fit_global(0,beta_extrap[3],pval_chiral_r0_global)*TMath::Sqrt(
      TMath::Power(perr_chiral_r0_global[0],2.)+TMath::Power(beta_extrap[3]*perr_chiral_r0_global[1],2.)+2*beta_extrap[3]*cov(0,1) )};
  TGraphErrors* extrap_chiral_r0_540=new TGraphErrors(1,extrap_x_chiral_r0_540,extrap_y_chiral_r0_540,extrap_x_err_chiral_r0_540,extrap_y_err_chiral_r0_540);
  extrap_chiral_r0_540->SetMarkerColor(1);
  extrap_chiral_r0_540->SetMarkerStyle(1);
  extrap_chiral_r0_540->SetMarkerSize(1);
  extrap_chiral_r0_540->Draw("P");

  chiralr0canv->SaveAs("r0_Over_a_Chiral_Extrap_new.pdf");

  cout << "\n\n Extrapolated values for r0/a \n";
  cout << "Beta=5.20 : " << fit_global(0.,5.20,pval_chiral_r0_global) << " +/- " << fit_global_err(0.,0.,5.20,pval_chiral_r0_global,cov) << endl;
  cout << "Beta=5.25 : " << fit_global(0.,5.25,pval_chiral_r0_global) << " +/- " << fit_global_err(0.,0.,5.25,pval_chiral_r0_global,cov) << endl;
  cout << "Beta=5.29 : " << fit_global(0.,5.29,pval_chiral_r0_global) << " +/- " << fit_global_err(0.,0.,5.29,pval_chiral_r0_global,cov) << endl;
  cout << "Beta=5.40 : " << fit_global(0.,5.40,pval_chiral_r0_global) << " +/- " << fit_global_err(0.,0.,5.40,pval_chiral_r0_global,cov) << endl << endl << endl;
  cout << "Beta=5.55 : " << fit_global(0.,5.55,pval_chiral_r0_global) << " +/- " << fit_global_err(0.,0.,5.55,pval_chiral_r0_global,cov) << endl << endl << endl;

  for(Int_t i=0;i<NFITPLQ;++i){
    r0_over_a_extrap[i]=fit_global( amq(i,1) , beta_extrap[i] , pval_chiral_r0_global);
    r0_over_a_err_extrap[i]=fit_global_err( amq(i,1) , amq_err(i,1) , beta_extrap[i] , pval_chiral_r0_global , cov );
  };

  cout << "r0_over_a_extrap \n";
  for(Int_t i=0;i<NFITPLQ;++i){
    cout << r0_over_a_extrap[i] << " ," << flush;
    if( (i==4) || (i==9) || (i==16) || (i==23) ) { cout << endl; }
  }

  cout << "\n\nr0_over_a_err_extrap \n";
  for(Int_t i=0;i<NFITPLQ;++i){
    cout << r0_over_a_err_extrap[i] << " ," << flush;
    if( (i==4) || (i==9) || (i==16) || (i==23) ) { cout << endl; }
  }

  cout << endl << endl;

  // Now we extrapolate to the beta and kappa values needed for the plaquette data

  return 0;

};


