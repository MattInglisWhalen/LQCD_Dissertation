#include <iostream>
#include <climits>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TComplex.h"

#include "BetaFn.h"

Double_t TF1_betafn(Double_t* x, Double_t* par){
  BetaFn bf(par[0],par[1],par[2]); //par[0]==nloops, par[1]==nflav, par[2]==variable
  return bf.eval(x[0]);
};

Double_t integrand(Double_t x, Int_t nloops, Int_t nflav, Int_t var){
  if(x==0.){
    return 0.;
  }
  BetaFn bf(nloops,nflav,var); 
  if(var==2){
    return (   ( 1./( bf.eval(x)) ) + ( 1./(x*x*x*bf[0]) ) - ( bf[1]/(x*bf[0]*bf[0]) )   );
  }
  return 0.5*( 1./( bf.eval(x)) + 1./(x*x*bf[0]) - bf[1]/(x*bf[0]*bf[0]) );
};
Double_t TF1_integrand(Double_t* x, Double_t* par){
  return integrand(x[0],(Int_t)par[0],(Int_t)par[1],(Int_t)par[2]);
};

// numerically integrating and multiplying in the prefactors to get mu over lambda
Double_t mu_over_lambda(Double_t x, Int_t nloops, Int_t nflav, Int_t var, Int_t opt){
  std::cout<<nloops<<"  "<<nflav<<"  "<<var<<"  "<<opt<<std::endl;
  BetaFn bf(nloops,nflav,var);  
  TF1 scratch("scratch",TF1_integrand,0.001,1,3);
  scratch.SetParameter(0,nloops);
  scratch.SetParameter(1,nflav);
  scratch.SetParameter(2,var);
  Double_t term1,term2,term3;
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
  term3=scratch.Integral(0.,x);
  std::cout<<term1<<" + "<<term2<<" + "<<term3<<" = "<<term1+term2+term3<<std::endl;
  return TMath::Exp( term1+term2+term3 );
};
Double_t TF1_mu_over_lambda(Double_t* x, Double_t* par){
  return mu_over_lambda(x[0],(Int_t)par[0],(Int_t)par[1],(Int_t)par[2],(Int_t)par[3]);
};

Double_t alpha_of_mu(Double_t mu, Double_t nloops, Double_t nflav, Double_t var, Double_t opt){
  TF1 scratch("scratch",TF1_mu_over_lambda,0.0001,1000,4);
  scratch.SetParameter(0,nloops);
  scratch.SetParameter(1,nflav);
  scratch.SetParameter(2,var);
  scratch.SetParameter(3,opt);
  Double_t x=scratch.GetX(mu,0.0001,1000);
  if(var==2){return ( x*x/(4.*TMath::Pi()) );}
  if(var==0){return 4.*x*TMath::Pi();}
  return x;
};
Double_t TF1_alpha_of_mu(Double_t* x, Double_t* par){
  return alpha_of_mu(x[0],par[0],par[1],par[2],par[3]);
};

Int_t main(){

  // testing and debugging
  BetaFn bf_g(4,2,2);
  BetaFn bf_alpha(4,2,1);
  BetaFn bf_a(4,2,0);

  // ---

  for(Int_t i=0;i<4;++i){
    std::cout << bf_g[i] << "  ";
  }
  std::cout<<std::endl;
  for(Int_t i=0;i<4;++i){
    std::cout << bf_alpha[i] << "  ";
  }
  std::cout<<std::endl;
  for(Int_t i=0;i<4;++i){
    std::cout << bf_a[i] << "  ";
  }
  std::cout<<std::endl<<std::endl;

  // ---

  std::cout << bf_g.eval(0.5) << std::endl;
  std::cout << bf_alpha.eval(0.5) << std::endl;
  std::cout << bf_a.eval(0.5) << std::endl;

  std::cout << 1/bf_g.eval(0.5) << " + " << 1/(bf_g[0]*0.5*0.5*0.5) << " - " << bf_g[1]/(bf_g[0]*bf_g[0]*0.5)
            << " = " << 1/bf_g.eval(0.5)+1/(bf_g[0]*0.5*0.5*0.5)-bf_g[1]/(bf_g[0]*bf_g[0]*0.5) << std::endl << std::endl;

  // ---

  std::cout << mu_over_lambda(0.02,2,2,2,0) << std::endl;

  exit(-1);

  //chunk of code for plotting F(g) as a function of g

  TCanvas* myCanv1=new TCanvas("myCanv1","myCanv1 title",700,700);
  myCanv1->Divide(1,2);
  myCanv1->cd(1);
  gPad->SetLogy();
  gPad->SetTicky(1);

  TH1F* frame0=new TH1F("frame0","N_{f} = 0",1000,200,600);
  frame0->SetStats(0);
  frame0->SetMinimum(4);
  frame0->SetMaximum(1000);
  frame0->GetXaxis()->SetTitle("g^{#bar{MS}}(#mu)");
  frame0->GetXaxis()->SetTickLength(0.02);
  frame0->GetXaxis()->SetLabelSize(0.020);
  frame0->GetXaxis()->SetLimits(0.01,0.05); //minx/maxx limits
  frame0->GetYaxis()->SetTitle("#mu / #Lambda");
  frame0->GetYaxis()->SetTickLength(0.02);
  frame0->GetYaxis()->SetLabelSize(0.020);
  frame0->Draw(" ");

  TF1* FNF0_NL1=new TF1("Ftest1",TF1_mu_over_lambda,0.001,1,4);
  FNF0_NL1->SetParameter(0,1);  //nloop
  FNF0_NL1->SetParameter(1,0);  //nflav
  FNF0_NL1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL1->SetParameter(1,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL1->SetLineWidth(1.5);
  FNF0_NL1->SetLineColor(2);
  FNF0_NL1->SetLineStyle(2);
  FNF0_NL1->Draw("SAME");

  TF1* FNF0_NL2=new TF1("test2",TF1_mu_over_lambda,0.001,1,4);
  FNF0_NL2->SetParameter(0,2);  //nloop
  FNF0_NL2->SetParameter(1,0);  //nflav
  FNF0_NL2->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL2->SetParameter(1,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL2->SetLineWidth(1.5);
  FNF0_NL2->SetLineColor(3);
  FNF0_NL2->SetLineStyle(2);
  FNF0_NL2->Draw("SAME");

  TF1* FNF0_NL3=new TF1("test3",TF1_mu_over_lambda,0.001,1,4);
  FNF0_NL3->SetParameter(0,3);  //nloop
  FNF0_NL3->SetParameter(1,0);  //nflav
  FNF0_NL3->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL3->SetParameter(1,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL3->SetLineWidth(1.5);
  FNF0_NL3->SetLineColor(4);
  FNF0_NL3->SetLineStyle(2);
  FNF0_NL3->Draw("SAME");
/*
  TF1* FNF0_NL4=new TF1("test4",TF1_mu_over_lambda,0.001,1,4);
  FNF0_NL4->SetParameter(0,4);  //nloop
  FNF0_NL4->SetParameter(1,0);  //nflav
  FNF0_NL4->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL4->SetParameter(1,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL4->SetLineWidth(1.5);
  FNF0_NL4->SetLineColor(1);
  FNF0_NL4->SetLineStyle(2);
  FNF0_NL4->Draw("SAME");
*/
  TF1* FNF0_NL1_opt1=new TF1("Ftest1",TF1_mu_over_lambda,0.001,1,4);
  FNF0_NL1_opt1->SetParameter(0,1);  //nloop
  FNF0_NL1_opt1->SetParameter(1,0);  //nflav
  FNF0_NL1_opt1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL1_opt1->SetParameter(1,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL1_opt1->SetLineWidth(1.5);
  FNF0_NL1_opt1->SetLineColor(2);
  FNF0_NL1_opt1->SetLineStyle(1);
  FNF0_NL1_opt1->Draw("SAME");

  TF1* FNF0_NL2_opt1=new TF1("test2",TF1_mu_over_lambda,0.001,1,4);
  FNF0_NL2_opt1->SetParameter(0,2);  //nloop
  FNF0_NL2_opt1->SetParameter(1,0);  //nflav
  FNF0_NL2_opt1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL2_opt1->SetParameter(1,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL2_opt1->SetLineWidth(1.5);
  FNF0_NL2_opt1->SetLineColor(3);
  FNF0_NL2_opt1->SetLineStyle(1);
  FNF0_NL2_opt1->Draw("SAME");

  TF1* FNF0_NL3_opt1=new TF1("test3",TF1_mu_over_lambda,0.001,1,4);
  FNF0_NL3_opt1->SetParameter(0,3);  //nloop
  FNF0_NL3_opt1->SetParameter(1,0);  //nflav
  FNF0_NL3_opt1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL3_opt1->SetParameter(1,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL3_opt1->SetLineWidth(1.5);
  FNF0_NL3_opt1->SetLineColor(4);
  FNF0_NL3_opt1->SetLineStyle(1);
  FNF0_NL3_opt1->Draw("SAME");
/*
  TF1* FNF0_NL4_opt1=new TF1("test4",TF1_mu_over_lambda,0.001,1,4);
  FNF0_NL4_opt1->SetParameter(0,4);  //nloop
  FNF0_NL4_opt1->SetParameter(1,0);  //nflav
  FNF0_NL4_opt1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF0_NL4_opt1->SetParameter(1,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF0_NL4_opt1->SetLineWidth(1.5);
  FNF0_NL4_opt1->SetLineColor(1);
  FNF0_NL4_opt1->SetLineStyle(1);
  FNF0_NL4_opt1->Draw("SAME");
*/
  myCanv1->cd(2);
  gPad->SetLogy();
  gPad->SetTicky(1);

  // now with nflav=2

  TH1F* frame2=new TH1F("frame2","N_{f} = 2",1000,200,600);
  frame2->SetStats(0);
  frame2->SetMinimum(4);
  frame2->SetMaximum(1000);
  frame2->GetXaxis()->SetTitle("g^{#bar{MS}}(#mu)");
  frame2->GetXaxis()->SetTickLength(0.02);
  frame2->GetXaxis()->SetLabelSize(0.020);
  frame2->GetXaxis()->SetLimits(0.01,0.05); //minx/maxx limits
  frame2->GetYaxis()->SetTitle("#mu / #Lambda");
  frame2->GetYaxis()->SetTickLength(0.02);
  frame2->GetYaxis()->SetLabelSize(0.020);
  frame2->Draw(" ");

  TF1* FNF2_NL1=new TF1("Ftest1",TF1_mu_over_lambda,0.001,1,4);
  FNF2_NL1->SetParameter(0,1);  //nloop
  FNF2_NL1->SetParameter(1,2);  //nflav
  FNF2_NL1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL1->SetParameter(1,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL1->SetLineWidth(1.5);
  FNF2_NL1->SetLineColor(2);
  FNF2_NL1->SetLineStyle(2);
  FNF2_NL1->Draw("SAME");

  TF1* FNF2_NL2=new TF1("test2",TF1_mu_over_lambda,0.001,1,4);
  FNF2_NL2->SetParameter(0,2);  //nloop
  FNF2_NL2->SetParameter(1,2);  //nflav
  FNF2_NL2->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL2->SetParameter(1,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL2->SetLineWidth(1.5);
  FNF2_NL2->SetLineColor(3);
  FNF2_NL2->SetLineStyle(2);
  FNF2_NL2->Draw("SAME");

  TF1* FNF2_NL3=new TF1("test3",TF1_mu_over_lambda,0.001,1,4);
  FNF2_NL3->SetParameter(0,3);  //nloop
  FNF2_NL3->SetParameter(1,2);  //nflav
  FNF2_NL3->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL3->SetParameter(1,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL3->SetLineWidth(1.5);
  FNF2_NL3->SetLineColor(4);
  FNF2_NL3->SetLineStyle(2);
  FNF2_NL3->Draw("SAME");
/*
  TF1* FNF2_NL4=new TF1("test4",TF1_mu_over_lambda,0.001,1,4);
  FNF2_NL4->SetParameter(0,4);  //nloop
  FNF2_NL4->SetParameter(1,2);  //nflav
  FNF2_NL4->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL4->SetParameter(1,0);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL4->SetLineWidth(1.5);
  FNF2_NL4->SetLineColor(1);
  FNF2_NL4->SetLineStyle(2);
  FNF2_NL4->Draw("SAME");
*/
  TF1* FNF2_NL1_opt1=new TF1("Ftest1",TF1_mu_over_lambda,0.001,1,4);
  FNF2_NL1_opt1->SetParameter(0,1);  //nloop
  FNF2_NL1_opt1->SetParameter(1,2);  //nflav
  FNF2_NL1_opt1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL1_opt1->SetParameter(1,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL1_opt1->SetLineWidth(1.5);
  FNF2_NL1_opt1->SetLineColor(2);
  FNF2_NL1_opt1->SetLineStyle(1);
  FNF2_NL1_opt1->Draw("SAME");

  TF1* FNF2_NL2_opt1=new TF1("test2",TF1_mu_over_lambda,0.001,1,4);
  FNF2_NL2_opt1->SetParameter(0,2);  //nloop
  FNF2_NL2_opt1->SetParameter(1,2);  //nflav
  FNF2_NL2_opt1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL2_opt1->SetParameter(1,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL2_opt1->SetLineWidth(1.5);
  FNF2_NL2_opt1->SetLineColor(3);
  FNF2_NL2_opt1->SetLineStyle(1);
  FNF2_NL2_opt1->Draw("SAME");

  TF1* FNF2_NL3_opt1=new TF1("test3",TF1_mu_over_lambda,0.001,1,4);
  FNF2_NL3_opt1->SetParameter(0,3);  //nloop
  FNF2_NL3_opt1->SetParameter(1,2);  //nflav
  FNF2_NL3_opt1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL3_opt1->SetParameter(1,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL3_opt1->SetLineWidth(1.5);
  FNF2_NL3_opt1->SetLineColor(4);
  FNF2_NL3_opt1->SetLineStyle(1);
  FNF2_NL3_opt1->Draw("SAME");
/*
  TF1* FNF2_NL4_opt1=new TF1("test4",TF1_mu_over_lambda,0.001,1,4);
  FNF2_NL4_opt1->SetParameter(0,4);  //nloop
  FNF2_NL4_opt1->SetParameter(1,2);  //nflav
  FNF2_NL4_opt1->SetParameter(1,2);  //var (g=2, alpha=1, a=0)
  FNF2_NL4_opt1->SetParameter(1,1);  //opt (0 : kappa=b_0  --  1: kappa=b_1/b_0 )
  FNF2_NL4_opt1->SetLineWidth(1.5);
  FNF2_NL4_opt1->SetLineColor(1);
  FNF2_NL4_opt1->SetLineStyle(1);
  FNF2_NL4_opt1->Draw("SAME");
*/
  myCanv1->SaveAs("mu_over_lambda_numerical.pdf");


 //chunk of code for plotting alpha as a function of mu over lambda, numerically with infinity limits

  TCanvas* myCanv=new TCanvas("myCanv","myCanv title",700,700);
  myCanv->SetLogx();

  TH1F* frame=new TH1F("frame","N_{f} = 0",1000,200,600);
  frame->SetStats(0);
  frame->SetMinimum(0.0);
  frame->SetMaximum(0.5);
  frame->GetXaxis()->SetTitle("#mu / #Lambda^{#bar{MS}}");
  frame->GetXaxis()->SetTickLength(0.02);
  frame->GetXaxis()->SetLabelSize(0.020);
  frame->GetXaxis()->SetLimits(1,1000); //minx/maxx limits
  frame->GetYaxis()->SetTitle("#alpha_{s}^{#bar{MS}}(#mu)");
  frame->GetYaxis()->SetTickLength(0.02);
  frame->GetYaxis()->SetLabelSize(0.020);
  frame->Draw(" ");

  TF1* NF0_NL1=new TF1("test1",TF1_alpha_of_mu,4,1000,4);
  NF0_NL1->SetParameter(0,1);
  NF0_NL1->SetParameter(1,0);
  NF0_NL1->SetParameter(2,2);
  NF0_NL1->SetParameter(3,0);
  NF0_NL1->SetLineWidth(1.5);
  NF0_NL1->SetLineColor(2);
  NF0_NL1->SetLineStyle(2);
  NF0_NL1->Draw("SAME");

  TF1* NF0_NL2=new TF1("test2",TF1_alpha_of_mu,4,1000,4);
  NF0_NL2->SetParameter(0,2);
  NF0_NL2->SetParameter(1,0);
  NF0_NL2->SetParameter(2,2);
  NF0_NL2->SetParameter(3,0);
  NF0_NL2->SetLineWidth(1.5);
  NF0_NL2->SetLineColor(3);
  NF0_NL2->SetLineStyle(3);
  NF0_NL2->Draw("SAME");

  TF1* NF0_NL3=new TF1("test3",TF1_alpha_of_mu,4,1000,4);
  NF0_NL3->SetParameter(0,3);
  NF0_NL3->SetParameter(1,0);
  NF0_NL3->SetParameter(2,2);
  NF0_NL3->SetParameter(3,0);
  NF0_NL3->SetLineWidth(1.5);
  NF0_NL3->SetLineColor(4);
  NF0_NL3->SetLineStyle(4);
  NF0_NL3->Draw("SAME");

  TF1* NF0_NL4=new TF1("test4",TF1_alpha_of_mu,4,1000,4);
  NF0_NL4->SetParameter(0,4);
  NF0_NL4->SetParameter(1,0);
  NF0_NL4->SetParameter(2,2);
  NF0_NL4->SetParameter(3,0);
  NF0_NL4->SetLineWidth(1.5);
  NF0_NL4->SetLineColor(1);
  NF0_NL4->SetLineStyle(1);
  NF0_NL4->Draw("SAME");

  myCanv->SaveAs("alpha_of_mu_numerical.pdf");
/*
  TCanvas* myCanv2=new TCanvas("myCanv2","myCanv2 title",700,700);
  myCanv2->SetLogx();

  TH1F* frame2=new TH1F("frame2","N_{f} = 2",1000,200,600);
  frame2->SetStats(0);
  frame2->SetMinimum(0.0);
  frame2->SetMaximum(0.5);
  frame2->GetXaxis()->SetTitle("#mu / #Lambda^{#bar{MS}}");
  frame2->GetXaxis()->SetTickLength(0.02);
  frame2->GetXaxis()->SetLabelSize(0.020);
  frame2->GetXaxis()->SetLimits(1,1000); //minx/maxx limits
  frame2->GetYaxis()->SetTitle("#alpha_{s}^{#bar{MS}}(#mu)");
  frame2->GetYaxis()->SetTickLength(0.02);
  frame2->GetYaxis()->SetLabelSize(0.020);
  frame2->Draw(" ");

  TF1* NF2_NL1=new TF1("test11",TF1_alpha_of_mu,4,1000,2);
  NF2_NL1->SetParameter(0,2);
  NF2_NL1->SetParameter(1,1);
  NF2_NL1->SetLineWidth(1.5);
  NF2_NL1->SetLineColor(2);
  NF2_NL1->SetLineStyle(2);
  NF2_NL1->Draw("SAME");

  TF1* NF2_NL2=new TF1("test12",TF1_alpha_of_mu,4,1000,2);
  NF2_NL2->SetParameter(0,2);
  NF2_NL2->SetParameter(1,2);
  NF2_NL2->SetLineWidth(1.5);
  NF2_NL2->SetLineColor(3);
  NF2_NL2->SetLineStyle(3);
  NF2_NL2->Draw("SAME");

  TF1* NF2_NL3=new TF1("test13",TF1_alpha_of_mu,4,1000,2);
  NF2_NL3->SetParameter(0,2);
  NF2_NL3->SetParameter(1,3);
  NF2_NL3->SetLineWidth(1.5);
  NF2_NL3->SetLineColor(4);
  NF2_NL3->SetLineStyle(4);
  NF2_NL3->Draw("SAME");

  TF1* NF2_NL4=new TF1("test14",TF1_alpha_of_mu,4,1000,2);
  NF2_NL4->SetParameter(0,2);
  NF2_NL4->SetParameter(1,4);
  NF2_NL4->SetLineWidth(1.5);
  NF2_NL4->SetLineColor(1);
  NF2_NL4->SetLineStyle(1);
  NF2_NL4->Draw("SAME");

  myCanv2->SaveAs("NF2.pdf");


  //chunk of code for plotting alpha as a function of mu over lambda with roger's analytic formula

  TCanvas* myCanv3=new TCanvas("myCanv3","myCanv3 title",700,700);
  myCanv3->SetLogx();

  TH1F* frame3=new TH1F("frame3","Real Analytic N_{f} = 2",1000,200,600);
  frame3->SetStats(0);
  frame3->SetMinimum(0.0);
  frame3->SetMaximum(0.5);
  frame3->GetXaxis()->SetTitle("#mu / #Lambda^{#bar{MS}}");
  frame3->GetXaxis()->SetTickLength(0.02);
  frame3->GetXaxis()->SetLabelSize(0.020);
  frame3->GetXaxis()->SetLimits(1,1000); //minx/maxx limits
  frame3->GetYaxis()->SetTitle("#alpha_{s}^{#bar{MS}}(#mu)");
  frame3->GetYaxis()->SetTickLength(0.02);
  frame3->GetYaxis()->SetLabelSize(0.020);
  frame3->Draw(" ");

  TF1* RNF2_NL1_Re=new TF1("test11",TF1_Ralpha_of_mu,4,1000,4);
  RNF2_NL1_Re->SetParameter(0,2);
  RNF2_NL1_Re->SetParameter(1,1);
  RNF2_NL1_Re->SetParameter(2,1);
  RNF2_NL1_Re->SetParameter(3,0);
  RNF2_NL1_Re->SetLineWidth(1.5);
  RNF2_NL1_Re->SetLineColor(2);
  RNF2_NL1_Re->SetLineStyle(2);
  RNF2_NL1_Re->Draw("SAME");

  TF1* RNF2_NL2_Re=new TF1("test12",TF1_Ralpha_of_mu,4,1000,4);
  RNF2_NL2_Re->SetParameter(0,2);
  RNF2_NL2_Re->SetParameter(1,2);
  RNF2_NL2_Re->SetParameter(2,1);
  RNF2_NL2_Re->SetParameter(3,0);
  RNF2_NL2_Re->SetLineWidth(1.5);
  RNF2_NL2_Re->SetLineColor(3);
  RNF2_NL2_Re->SetLineStyle(3);
  RNF2_NL2_Re->Draw("SAME");

  TF1* RNF2_NL3_Re=new TF1("test13",TF1_Ralpha_of_mu,4,1000,4);
  RNF2_NL3_Re->SetParameter(0,2);
  RNF2_NL3_Re->SetParameter(1,3);
  RNF2_NL3_Re->SetParameter(2,1);
  RNF2_NL3_Re->SetParameter(3,0);
  RNF2_NL3_Re->SetLineWidth(1.5);
  RNF2_NL3_Re->SetLineColor(4);
  RNF2_NL3_Re->SetLineStyle(4);
  RNF2_NL3_Re->Draw("SAME");

  myCanv3->SaveAs("RNF2_Re_b0=1.pdf");

  TCanvas* myCanv4=new TCanvas("myCanv4","myCanv4 title",700,700);
  myCanv4->SetLogx();
  myCanv4->SetLogy();

  TH1F* frame4=new TH1F("frame4","Imaginary Analytic N_{f} = 2",1000,200,600);
  frame4->SetStats(0);
  frame4->SetMinimum(0.001);
  frame4->SetMaximum(100.);
  frame4->GetXaxis()->SetTitle("#mu / #Lambda^{#bar{MS}}");
  frame4->GetXaxis()->SetTickLength(0.02);
  frame4->GetXaxis()->SetLabelSize(0.020);
  frame4->GetXaxis()->SetLimits(1,1000); //minx/maxx limits
  frame4->GetYaxis()->SetTitle("#alpha_{s}^{#bar{MS}}(#mu)");
  frame4->GetYaxis()->SetTickLength(0.02);
  frame4->GetYaxis()->SetLabelSize(0.020);
  frame4->Draw(" ");

  TF1* RNF2_NL1_Im=new TF1("test11",TF1_Ralpha_of_mu,4,1000,4);
  RNF2_NL1_Im->SetParameter(0,2);
  RNF2_NL1_Im->SetParameter(1,1);
  RNF2_NL1_Im->SetParameter(2,0);
  RNF2_NL1_Im->SetParameter(3,0);
  RNF2_NL1_Im->SetLineWidth(1.5);
  RNF2_NL1_Im->SetLineColor(2);
  RNF2_NL1_Im->SetLineStyle(2);
  RNF2_NL1_Im->Draw("SAME");

  TF1* RNF2_NL2_Im=new TF1("test12",TF1_Ralpha_of_mu,4,1000,4);
  RNF2_NL2_Im->SetParameter(0,2);
  RNF2_NL2_Im->SetParameter(1,2);
  RNF2_NL2_Im->SetParameter(2,0);
  RNF2_NL2_Im->SetParameter(3,0);
  RNF2_NL2_Im->SetLineWidth(1.5);
  RNF2_NL2_Im->SetLineColor(3);
  RNF2_NL2_Im->SetLineStyle(3);
  RNF2_NL2_Im->Draw("SAME");

  TF1* RNF2_NL3_Im=new TF1("test13",TF1_Ralpha_of_mu,4,1000,4);
  RNF2_NL3_Im->SetParameter(0,2);
  RNF2_NL3_Im->SetParameter(1,3);
  RNF2_NL3_Im->SetParameter(2,0);
  RNF2_NL3_Im->SetParameter(3,0);
  RNF2_NL3_Im->SetLineWidth(1.5);
  RNF2_NL3_Im->SetLineColor(4);
  RNF2_NL3_Im->SetLineStyle(4);
  RNF2_NL3_Im->Draw("SAME");

  myCanv4->SaveAs("RNF2_Im_b0=1.pdf");
*/
/*
  std::cout << "TAG1\n";
  for(Int_t i=0;i<4;++i){
    std::cout << Ralpha_of_mu(TMath::Power(10,i) , 2 , 1 , 1) << std::endl;
  }
  for(Int_t i=0;i<4;++i){
    std::cout << Ralpha_of_mu(TMath::Power(10,i) , 2 , 2 , 1) << std::endl;
  }
  for(Int_t i=0;i<4;++i){
    std::cout << Ralpha_of_mu(TMath::Power(10,i) , 2 , 3 , 1) << std::endl;
  }
  std::cout << "TAG1\n";
  for(Int_t i=0;i<4;++i){
    std::cout << Ralpha_of_mu(TMath::Power(10,i) , 2 , 1 , 0) << std::endl;
  }
  for(Int_t i=0;i<4;++i){
    std::cout << Ralpha_of_mu(TMath::Power(10,i) , 2 , 2 , 0) << std::endl;
  }
  for(Int_t i=0;i<4;++i){
    std::cout << Ralpha_of_mu(TMath::Power(10,i) , 2 , 3 , 0) << std::endl;
  }
  std::cout << "TAG1\n";
  std::cout << "TAG2\n";
  for(Int_t i=0;i<4;++i){
    std::cout << RF_alpha( TMath::Power(10,i) , 2 , 1 , 1) << std::endl;
  }
  for(Int_t i=0;i<4;++i){
    std::cout << RF_alpha( TMath::Power(10,i) , 2 , 2 , 1) << std::endl;
  }
  for(Int_t i=0;i<4;++i){
    std::cout << RF_alpha( TMath::Power(10,i) , 2 , 3 , 1) << std::endl;
  }
  std::cout << "TAG2\n";
  for(Int_t i=0;i<4;++i){
    //std::cout << RF_alpha( TMath::Power(10,i) , 2 , 1 , 0) << std::endl;
  }
  for(Int_t i=0;i<4;++i){
    //std::cout << RF_alpha( TMath::Power(10,i) , 2 , 2 , 0) << std::endl;
  }
  for(Int_t i=-5;i<4;++i){
    std::cout << RF_alpha( TMath::Power(10,i) , 2 , 3 , 0) << std::endl;
  }
  std::cout << "TAG2\n";
*/

  return 0;

}
