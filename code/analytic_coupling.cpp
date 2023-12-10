#include <iostream>
#include <climits>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TLegend.h"
#include "TMath.h"
#include "TComplex.h"

#include "BetaFn.h"

Double_t mu_over_lambda(Double_t x, Int_t nloop, Int_t nflav, Int_t var, Int_t opt){
  BetaFn bf(nloop,nflav,var);
  std::vector<TComplex> P(nloop-1,0.);
  std::vector<TComplex> roots=bf.nonzeroroots();
  for(Int_t i=0;i<nloop-1;++i){
    P[i]=(-1./(bf[nloop-1]*roots[i]*roots[i]));
    for(Int_t j=0;j<nloop-1;++j){
      if(i!=j){
        P[i]*=( 1./(roots[i]-roots[j]) );
      }
    }
  }
  
  TComplex term1=0.,term2=0.,ln_arg=1.;
  if(var==2){ //using g as a variable
    term1=1./(2.*bf[0]*x*x);
  }
  else{ //using alpha or a as a variable
    term1=1./(2.*bf[0]*x);
  }
  if(var==2){  //use g as a variable
    if(opt==0){  //use kappa=b_0
      term2=( bf[1]/(2.*bf[0]*bf[0]) )*TMath::Log( bf[0]*x*x );
    }
    else{  //use kappa=b_1/b_0
      term2=( bf[1]/(2.*bf[0]*bf[0]) )*TMath::Log( bf[1]*x*x/bf[0] );
    } 
  }
  else{  //use alpha or a as a variable
    if(opt==0){  //use kappa=b_0
      term2=( bf[1]/(2.*bf[0]*bf[0]) )*TMath::Log( bf[0]*x );
    }
    else{  //use kappa=b_1/b_0
      term2=( bf[1]/(2.*bf[0]*bf[0]) )*TMath::Log( bf[1]*x/bf[0] );
    } 
  }
  for(Int_t i=0;i<nloop-1;++i){
    if(var>1){
      ln_arg*=TComplex::Power( ((TComplex)1.)-x*x/roots[i],P[i]/2. );
    }
    else{
      ln_arg*=TComplex::Power( ((TComplex)1.)-x/roots[i],P[i]/2. );
    }
  }
  if(nloop<2){
    return TMath::Exp(term1);
  }
  return TComplex::Exp( term1+term2+TComplex::Log(ln_arg)).Re();
};
Double_t TF1_mu_over_lambda(Double_t* x, Double_t* par){
  return mu_over_lambda(x[0],(Int_t)par[0],(Int_t)par[1],(Int_t)par[2],(Int_t)par[3]);
};

Double_t alpha_of_mu(Double_t mu, Int_t nloops, Int_t nflav, Int_t var, Int_t opt){
  TF1 scratch("scratch",TF1_mu_over_lambda,0.0001,1000,4);
  scratch.SetParameter(0,nloops);
  scratch.SetParameter(1,nflav);
  scratch.SetParameter(2,var);
  scratch.SetParameter(3,opt);
  Double_t x=scratch.GetX(mu,0.001,10);
  if(var==2){return ( x*x/(4.*TMath::Pi()) );}
  if(var==0){return 4.*x*TMath::Pi();}
  return x;
};
Double_t TF1_alpha_of_mu(Double_t* x, Double_t* par){
  return alpha_of_mu(x[0],(Int_t)par[0],(Int_t)par[1],(Int_t)par[2],(Int_t)par[3]);
};

Double_t Rogers_mu_over_lambda(Double_t g, Int_t nloop, Int_t nflav){
  BetaFn bf(nloop,nflav,2);
  TComplex A=bf[1]+TComplex::Sqrt(bf[1]*bf[1]-4.*bf[0]*bf[2]);
  TComplex B=bf[1]-TComplex::Sqrt(bf[1]*bf[1]-4.*bf[0]*bf[2]);
  TComplex pA=-((bf[1])/(4.*bf[0]*bf[0]))
              -((bf[1]*bf[1]-2.*bf[0]*bf[2])/(4.*bf[0]*bf[0]*TComplex::Sqrt(bf[1]*bf[1]-4.*bf[0]*bf[2])));
  TComplex pB=-((bf[1])/(4.*bf[0]*bf[0]))
              +((bf[1]*bf[1]-2.*bf[0]*bf[2])/(4.*bf[0]*bf[0]*TComplex::Sqrt(bf[1]*bf[1]-4.*bf[0]*bf[2])));
  TComplex fac1=TMath::Exp(1./(2.*bf[0]*g*g));
  TComplex fac2=TMath::Power(bf[0]*g*g,bf[1]/(2.*bf[0]*bf[0]));
  TComplex fac3=TComplex::Power( TComplex(1.)+A*g*g/(2*bf[0]),pA);
  TComplex fac4=TComplex::Power( TComplex(1.)+B*g*g/(2*bf[0]),pB);
  if(nloop<2){
    return fac1;
  }
  return (fac1*fac2*fac3*fac4).Re();
};

Int_t main(){


  std::cout << mu_over_lambda(0.2,1,2,1,0) << std::endl;
  std::cout << mu_over_lambda(0.2,2,2,1,0) << std::endl;
  std::cout << mu_over_lambda(0.2,3,2,1,0) << std::endl;

  std::cout << Rogers_mu_over_lambda(1.5853309,1,2) << std::endl;
  std::cout << Rogers_mu_over_lambda(1.5853309,2,2) << std::endl;
  std::cout << Rogers_mu_over_lambda(1.5853309,3,2) << std::endl;

  std::cout << alpha_of_mu(10,3,2,0,0) << std::endl;
  std::cout << alpha_of_mu(10,3,2,1,0) << std::endl;
  std::cout << alpha_of_mu(10,3,2,2,0) << std::endl;


  //exit(-1);

  // begin program

  TCanvas* myCanv1=new TCanvas("myCanv1","myCanv1 title",700,700);
  myCanv1->Divide(1,2);
  myCanv1->cd(1);
  gPad->SetLogy();
  gPad->SetTicky(1);

  TH1F* frame0=new TH1F("frame0","N_{f} = 0",1000,200,600);
  frame0->SetStats(0);
  frame0->SetMinimum(4);
  frame0->SetMaximum(1000);
  frame0->GetXaxis()->SetTitle("g");
  frame0->GetXaxis()->SetTickLength(0.02);
  frame0->GetXaxis()->SetLabelSize(0.020);
  frame0->GetXaxis()->SetLimits(1,2.5); //minx/maxx limits
  frame0->GetYaxis()->SetTitle("#mu / #Lambda");
  frame0->GetYaxis()->SetTickLength(0.02);
  frame0->GetYaxis()->SetLabelSize(0.020);
  frame0->Draw(" ");

  TF1* FNF0_NL1=new TF1("Ftest1",TF1_mu_over_lambda,1,10,4);
  FNF0_NL1->SetParameter(0,1);  //nloop
  FNF0_NL1->SetParameter(1,0);  //nflav
  FNF0_NL1->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL1->SetParameter(3,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL1->SetLineWidth(1.5);
  FNF0_NL1->SetLineColor(2);
  FNF0_NL1->SetLineStyle(2);
  FNF0_NL1->Draw("SAME");

  TF1* FNF0_NL2=new TF1("test2",TF1_mu_over_lambda,1,10,4);
  FNF0_NL2->SetParameter(0,2);  //nloop
  FNF0_NL2->SetParameter(1,0);  //nflav
  FNF0_NL2->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL2->SetParameter(3,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL2->SetLineWidth(1.5);
  FNF0_NL2->SetLineColor(3);
  FNF0_NL2->SetLineStyle(2);
  FNF0_NL2->Draw("SAME");

  TF1* FNF0_NL3=new TF1("test3",TF1_mu_over_lambda,1,10,4);
  FNF0_NL3->SetParameter(0,3);  //nloop
  FNF0_NL3->SetParameter(1,0);  //nflav
  FNF0_NL3->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL3->SetParameter(3,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL3->SetLineWidth(1.5);
  FNF0_NL3->SetLineColor(4);
  FNF0_NL3->SetLineStyle(2);
  FNF0_NL3->Draw("SAME");
/*
  TF1* FNF0_NL4=new TF1("test4",TF1_mu_over_lambda,0.001,1000,4);
  FNF0_NL4->SetParameter(0,4);  //nloop
  FNF0_NL4->SetParameter(1,0);  //nflav
  FNF0_NL4->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL4->SetParameter(3,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL4->SetLineWidth(1.5);
  FNF0_NL4->SetLineColor(1);
  FNF0_NL4->SetLineStyle(2);
  FNF0_NL4->Draw("SAME");
*/
  TF1* FNF0_NL1_opt1=new TF1("Ftest1",TF1_mu_over_lambda,1,10,4);
  FNF0_NL1_opt1->SetParameter(0,1);  //nloop
  FNF0_NL1_opt1->SetParameter(1,0);  //nflav
  FNF0_NL1_opt1->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL1_opt1->SetParameter(3,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL1_opt1->SetLineWidth(1.5);
  FNF0_NL1_opt1->SetLineColor(2);
  FNF0_NL1_opt1->SetLineStyle(1);
  FNF0_NL1_opt1->Draw("SAME");

  TF1* FNF0_NL2_opt1=new TF1("test2",TF1_mu_over_lambda,1,10,4);
  FNF0_NL2_opt1->SetParameter(0,2);  //nloop
  FNF0_NL2_opt1->SetParameter(1,0);  //nflav
  FNF0_NL2_opt1->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL2_opt1->SetParameter(3,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL2_opt1->SetLineWidth(1.5);
  FNF0_NL2_opt1->SetLineColor(3);
  FNF0_NL2_opt1->SetLineStyle(1);
  FNF0_NL2_opt1->Draw("SAME");

  TF1* FNF0_NL3_opt1=new TF1("test3",TF1_mu_over_lambda,1,10,4);
  FNF0_NL3_opt1->SetParameter(0,3);  //nloop
  FNF0_NL3_opt1->SetParameter(1,0);  //nflav
  FNF0_NL3_opt1->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL3_opt1->SetParameter(3,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL3_opt1->SetLineWidth(1.5);
  FNF0_NL3_opt1->SetLineColor(4);
  FNF0_NL3_opt1->SetLineStyle(1);
  FNF0_NL3_opt1->Draw("SAME");
/*
  TF1* FNF0_NL4_opt1=new TF1("test4",TF1_mu_over_lambda,0,10,4);
  FNF0_NL4_opt1->SetParameter(0,4);  //nloop
  FNF0_NL4_opt1->SetParameter(1,0);  //nflav
  FNF0_NL4_opt1->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL4_opt1->SetParameter(3,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL4_opt1->SetLineWidth(1.5);
  FNF0_NL4_opt1->SetLineColor(1);
  FNF0_NL4_opt1->SetLineStyle(1);
  FNF0_NL4_opt1->Draw("SAME");
*/

  TLegend* leg = new TLegend(0.75,0.6,0.95,0.8,"Loops");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0); 
  leg->AddEntry(FNF0_NL1_opt1,"1","l");
  leg->AddEntry(FNF0_NL2_opt1,"2","l");
  leg->AddEntry(FNF0_NL3_opt1,"3","l");
  leg->Draw();

  myCanv1->cd(2);
  gPad->SetLogy();
  gPad->SetTicky(1);

  // now with nflav=2

  TH1F* frame2=new TH1F("frame2","N_{f} = 2",1000,200,600);
  frame2->SetStats(0);
  frame2->SetMinimum(4);
  frame2->SetMaximum(1000);
  frame2->GetXaxis()->SetTitle("g");
  frame2->GetXaxis()->SetTickLength(0.02);
  frame2->GetXaxis()->SetLabelSize(0.020);
  frame2->GetXaxis()->SetLimits(1.,2.5); //minx/maxx limits
  frame2->GetYaxis()->SetTitle("#mu / #Lambda");
  frame2->GetYaxis()->SetTickLength(0.02);
  frame2->GetYaxis()->SetLabelSize(0.020);
  frame2->Draw(" ");

  TF1* FNF2_NL1=new TF1("Ftest1",TF1_mu_over_lambda,1,10,4);
  FNF2_NL1->SetParameter(0,1);  //nloop
  FNF2_NL1->SetParameter(1,2);  //nflav
  FNF2_NL1->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL1->SetParameter(3,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL1->SetLineWidth(1.5);
  FNF2_NL1->SetLineColor(2);
  FNF2_NL1->SetLineStyle(2);
  FNF2_NL1->Draw("SAME");

  TF1* FNF2_NL2=new TF1("test2",TF1_mu_over_lambda,1,10,4);
  FNF2_NL2->SetParameter(0,2);  //nloop
  FNF2_NL2->SetParameter(1,2);  //nflav
  FNF2_NL2->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL2->SetParameter(3,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL2->SetLineWidth(1.5);
  FNF2_NL2->SetLineColor(3);
  FNF2_NL2->SetLineStyle(2);
  FNF2_NL2->Draw("SAME");

  TF1* FNF2_NL3=new TF1("test3",TF1_mu_over_lambda,1,10,4);
  FNF2_NL3->SetParameter(0,3);  //nloop
  FNF2_NL3->SetParameter(1,2);  //nflav
  FNF2_NL3->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL3->SetParameter(3,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL3->SetLineWidth(1.5);
  FNF2_NL3->SetLineColor(4);
  FNF2_NL3->SetLineStyle(2);
  FNF2_NL3->Draw("SAME");
/*
  TF1* FNF2_NL4=new TF1("test4",TF1_mu_over_lambda,0.001,1000,4);
  FNF2_NL4->SetParameter(0,4);  //nloop
  FNF2_NL4->SetParameter(1,2);  //nflav
  FNF2_NL4->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL4->SetParameter(1,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL4->SetLineWidth(1.5);
  FNF2_NL4->SetLineColor(1);
  FNF2_NL4->SetLineStyle(2);
  FNF2_NL4->Draw("SAME");
*/
  TF1* FNF2_NL1_opt1=new TF1("Ftest1",TF1_mu_over_lambda,1,10,4);
  FNF2_NL1_opt1->SetParameter(0,1);  //nloop
  FNF2_NL1_opt1->SetParameter(1,2);  //nflav
  FNF2_NL1_opt1->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL1_opt1->SetParameter(3,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL1_opt1->SetLineWidth(1.5);
  FNF2_NL1_opt1->SetLineColor(2);
  FNF2_NL1_opt1->SetLineStyle(1);
  FNF2_NL1_opt1->Draw("SAME");

  TF1* FNF2_NL2_opt1=new TF1("test2",TF1_mu_over_lambda,1,10,4);
  FNF2_NL2_opt1->SetParameter(0,2);  //nloop
  FNF2_NL2_opt1->SetParameter(1,2);  //nflav
  FNF2_NL2_opt1->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL2_opt1->SetParameter(3,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL2_opt1->SetLineWidth(1.5);
  FNF2_NL2_opt1->SetLineColor(3);
  FNF2_NL2_opt1->SetLineStyle(1);
  FNF2_NL2_opt1->Draw("SAME");

  TF1* FNF2_NL3_opt1=new TF1("test3",TF1_mu_over_lambda,1,10,4);
  FNF2_NL3_opt1->SetParameter(0,3);  //nloop
  FNF2_NL3_opt1->SetParameter(1,2);  //nflav
  FNF2_NL3_opt1->SetParameter(2,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL3_opt1->SetParameter(3,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL3_opt1->SetLineWidth(1.5);
  FNF2_NL3_opt1->SetLineColor(4);
  FNF2_NL3_opt1->SetLineStyle(1);
  FNF2_NL3_opt1->Draw("SAME");
/*
  TF1* FNF2_NL4_opt1=new TF1("test4",TF1_mu_over_lambda,0.001,1000,4);
  FNF2_NL4_opt1->SetParameter(0,4);  //nloop
  FNF2_NL4_opt1->SetParameter(1,2);  //nflav
  FNF2_NL4_opt1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL4_opt1->SetParameter(1,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL4_opt1->SetLineWidth(1.5);
  FNF2_NL4_opt1->SetLineColor(1);
  FNF2_NL4_opt1->SetLineStyle(1);
  FNF2_NL4_opt1->Draw("SAME");
*/
  myCanv1->SaveAs("mu_over_lambda_analytic.pdf");

// --------------------------------------------------------------

  TCanvas* myCanv3=new TCanvas("myCanv3","myCanv3 title",1000,600);
  myCanv3->Divide(2,1);
  myCanv3->cd(1);
  gPad->SetLogx();
  gPad->SetTicky(1);
  
  TH1F* frame3=new TH1F("frame3"," ",1000,200,600);
  frame3->SetStats(0);
  frame3->SetMinimum(0.0);
  frame3->SetMaximum(0.5);
  frame3->GetXaxis()->SetTitle("#mu / #Lambda");
  frame3->GetXaxis()->SetTickLength(0.02);
  frame3->GetXaxis()->SetLabelSize(0.020);
  frame3->GetXaxis()->SetLimits(1,1000); //minx/maxx limits
  frame3->GetYaxis()->SetTitle("#alpha");
  frame3->GetYaxis()->SetTickLength(0.02);
  frame3->GetYaxis()->SetLabelSize(0.020);
  frame3->Draw(" ");

  TF1* NF0_NL1=new TF1("test1",TF1_alpha_of_mu,4,1000,4);
  NF0_NL1->SetParameter(0,1);
  NF0_NL1->SetParameter(1,0);
  NF0_NL1->SetParameter(2,1);
  NF0_NL1->SetParameter(3,0);
  NF0_NL1->SetLineWidth(1.5);
  NF0_NL1->SetLineColor(2);
  NF0_NL1->SetLineStyle(2);
  NF0_NL1->Draw("SAME");

  TF1* NF0_NL2=new TF1("test2",TF1_alpha_of_mu,4,1000,4);
  NF0_NL2->SetParameter(0,2);
  NF0_NL2->SetParameter(1,0);
  NF0_NL2->SetParameter(2,1);
  NF0_NL2->SetParameter(3,0);
  NF0_NL2->SetLineWidth(1.5);
  NF0_NL2->SetLineColor(3);
  NF0_NL2->SetLineStyle(3);
  NF0_NL2->Draw("SAME");

  TF1* NF0_NL3=new TF1("test3",TF1_alpha_of_mu,4,1000,4);
  NF0_NL3->SetParameter(0,3);
  NF0_NL3->SetParameter(1,0);
  NF0_NL3->SetLineWidth(1.5);
  NF0_NL3->SetParameter(2,1);
  NF0_NL1->SetParameter(3,0);
  NF0_NL3->SetLineColor(4);
  NF0_NL3->SetLineStyle(4);
  NF0_NL3->Draw("SAME");

  TF1* NF0_NL4=new TF1("test4",TF1_alpha_of_mu,4,1000,4);
  NF0_NL4->SetParameter(0,4);
  NF0_NL4->SetParameter(1,0);
  NF0_NL4->SetParameter(2,1);
  NF0_NL4->SetParameter(3,0);
  NF0_NL4->SetLineWidth(1.5);
  NF0_NL4->SetLineColor(1);
  NF0_NL4->SetLineStyle(1);
  NF0_NL4->Draw("SAME");

//  TCanvas* myCanv4=new TCanvas("myCanv4","myCanv4 title",700,700);
//  myCanv4->SetLogx();

  myCanv3->cd(2);
  gPad->SetLogx();
  gPad->SetTicky(1);

  TH1F* frame4=new TH1F("frame4"," ",1000,200,600);
  frame4->SetStats(0);
  frame4->SetMinimum(0.0);
  frame4->SetMaximum(0.5);
  frame4->GetXaxis()->SetTitle("#mu / #Lambda");
  frame4->GetXaxis()->SetTickLength(0.02);
  frame4->GetXaxis()->SetLabelSize(0.020);
  frame4->GetXaxis()->SetLimits(1,1000); //minx/maxx limits
  frame4->GetYaxis()->SetTitle("#alpha");
  frame4->GetYaxis()->SetTickLength(0.02);
  frame4->GetYaxis()->SetLabelSize(0.020);
  frame4->Draw(" ");

  TF1* NF2_NL1=new TF1("test11",TF1_alpha_of_mu,4,1000,2);
  NF2_NL1->SetParameter(0,1);
  NF2_NL1->SetParameter(1,2);
  NF2_NL1->SetParameter(2,1);
  NF2_NL1->SetParameter(3,0);
  NF2_NL1->SetLineWidth(1.5);
  NF2_NL1->SetLineColor(2);
  NF2_NL1->SetLineStyle(2);
  NF2_NL1->Draw("SAME");

  TF1* NF2_NL2=new TF1("test12",TF1_alpha_of_mu,4,1000,2);
  NF2_NL2->SetParameter(0,2);
  NF2_NL2->SetParameter(1,2);
  NF2_NL2->SetParameter(2,1);
  NF2_NL2->SetParameter(3,0);
  NF2_NL2->SetLineWidth(1.5);
  NF2_NL2->SetLineColor(3);
  NF2_NL2->SetLineStyle(3);
  NF2_NL2->Draw("SAME");

  TF1* NF2_NL3=new TF1("test13",TF1_alpha_of_mu,4,1000,2);
  NF2_NL3->SetParameter(0,3);
  NF2_NL3->SetParameter(1,2);
  NF2_NL3->SetParameter(2,1);
  NF2_NL3->SetParameter(3,0);
  NF2_NL3->SetLineWidth(1.5);
  NF2_NL3->SetLineColor(4);
  NF2_NL3->SetLineStyle(4);
  NF2_NL3->Draw("SAME");

  TF1* NF2_NL4=new TF1("test14",TF1_alpha_of_mu,4,1000,2);
  NF2_NL4->SetParameter(0,4);
  NF2_NL4->SetParameter(1,2);
  NF2_NL4->SetParameter(2,1);
  NF2_NL4->SetParameter(3,0);
  NF2_NL4->SetLineWidth(1.5);
  NF2_NL4->SetLineColor(1);
  NF2_NL4->SetLineStyle(1);
  NF2_NL4->Draw("SAME");

  TLegend* leg2 = new TLegend(0.7,0.65,1.1,0.85,"Loops");
  leg2->SetBorderSize(0);
  leg2->SetFillStyle(0); 
  leg2->AddEntry(NF0_NL1,"1","l");
  leg2->AddEntry(NF0_NL2,"2","l");
  leg2->AddEntry(NF0_NL3,"3","l");
  leg2->AddEntry(NF0_NL4,"4","l");
  leg2->Draw();
// alpha_of_mu(Double_t mu, Int_t nloops, Int_t nflav, Int_t var, Int_t opt){
  std::cout << alpha_of_mu(5.98316,4,4,2,0) << std::endl;

  myCanv3->SaveAs("alpha_of_mu_analytic.pdf");

  return 0;

}
