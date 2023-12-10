#include <iostream>
#include <climits>

#include "TCanvas.h"
#include "TF1.h"
#include "TH1.h"
#include "TMath.h"
#include "TComplex.h"

// masses of the quarks in GeV
#define MASS_DOWN 0.003
#define MASS_UP 0.002
#define MASS_STRANGE 0.103
#define MASS_CHARM 1.5
#define MASS_BOTTOM 4.5
#define MASS_TOP 172.9

//beta function defined by
// \partial \alpha / \partial \ln \mu = \beta ( \alpha ) = - \beta_0 \alpha^2 - \beta_1 \alpha^3 - ...
//  with \alpha = g^2 / 4\pi 
Double_t betafn(Double_t alpha, Double_t Nf, Double_t order){

  Double_t RZ3=1.2020569031595942853997;  //Riemann zeta function at z=3.0

  Double_t beta=0.0;
  if(order>0){
    beta=beta-( 
                (11.) 
                      - Nf*(2./3.) 
                                   )*TMath::Power( ( alpha/(4.*TMath::Pi()) ) , 2 ) ;
  }
  if(order>1){
    beta=beta-(
                (102.) 
                      - Nf*(38./3.)
                                   )*TMath::Power( ( alpha/(4.*TMath::Pi()) ) , 3 ) ;
  }
  if(order>2){
    beta=beta-(
                (2857./2.)
                      -Nf*(5033./18.)
                                     + Nf*Nf*(325./54.)
                                   )*TMath::Power( ( alpha/(4.*TMath::Pi()) ) , 4 ) ;
  }
  if(order>3){
    beta=beta-(
               (149753./6. + 3564.*RZ3)
                      - Nf*(1078361./162.+6508.*RZ3/27.) 
                                     + Nf*Nf*(50065./162.+6472.*RZ3/81.)
                                                   + Nf*Nf*Nf*(1093./729.)
                                   )*TMath::Power( ( alpha/(4.*TMath::Pi()) ) , 5 ) ;
  }

  return 4*TMath::Pi()*beta;

};
Double_t TF1_betafn(Double_t* x, Double_t* par){
  return betafn(x[0],par[0],par[1]);
};

Double_t neg_recip_betafn(Double_t alpha, Double_t Nf, Double_t order){
  return -1/( 2*betafn(alpha,Nf,order) );
};
Double_t TF1_neg_recip_betafn(Double_t* x, Double_t* par){
  return neg_recip_betafn(x[0],par[0],par[1]);
};

// my formulae
Double_t F_alpha(Double_t alpha, Double_t Nf, Double_t order){
  TF1 scratch("scratch",TF1_neg_recip_betafn,0.0001,1000,2);
  scratch.SetParameter(0,Nf);
  scratch.SetParameter(1,order);
  return TMath::Exp( scratch.Integral(100,alpha) );
};
Double_t TF1_F_alpha(Double_t* x, Double_t* par){
  return F_alpha(x[0],par[0],par[1]);
};

Double_t alpha_of_mu(Double_t mu, Double_t Nf, Double_t order){
  TF1 scratch("scratch",TF1_F_alpha,0.0001,1000,2);
  scratch.SetParameter(0,Nf);
  scratch.SetParameter(1,order);
  return scratch.GetX(1/mu,0.0001,1000);
};
Double_t TF1_alpha_of_mu(Double_t* x, Double_t* par){
  return alpha_of_mu(x[0],par[0],par[1]);
};

//Roger's formulae
Double_t RF_alpha(Double_t g, Double_t Nf, Double_t order, Double_t isReal){

  Double_t b0=0.,b1=0.,b2=0.;

  if(order>0){
    b0= ( (11.) 
                        - Nf*(2./3.) 
                                     )*TMath::Power( 4.*TMath::Pi()  , -2 ) ;
  }
  if(order>1){
    b1=(  (102.) 
                        - Nf*(38./3.)
                                     )*TMath::Power( 4.*TMath::Pi() , -4 ) ;
  }
  if(order>2){
    b2=( (2857./2.)
                        - Nf*(5033./18.)
                                       + Nf*Nf*(325./54.)
                                     )*TMath::Power( 4.*TMath::Pi() , -6 ) ;
  }

  TComplex A=b1+TComplex::Sqrt(b1*b1-4*b0*b2);
  TComplex B=b1-TComplex::Sqrt(b1*b1-4*b0*b2);
  TComplex P_A=(-b1/(4.*b0*b0))-((b1*b1-2.*b0*b2)/(4.*b0*b0*TComplex::Sqrt(b1*b1-4.*b0*b2)));
  TComplex P_B=(-b1/(4.*b0*b0))+((b1*b1-2.*b0*b2)/(4.*b0*b0*TComplex::Sqrt(b1*b1-4.*b0*b2)));

  TComplex fac_1=TComplex::Exp(-1./(2.*b0*g*g));
  TComplex fac_2=TComplex::Power(TComplex(b0*g*g),TComplex(-b1/(2.*b0*b0)));
  TComplex fac_3=TComplex::Power(1.+A*g*g/(2.*b0),-P_A);
  TComplex fac_4=TComplex::Power(1.+B*g*g/(2.*b0),-P_B);

  std::cout << Nf << "  " << order << "  " << isReal << " : " << b0 << "  " << b1 << "  " << b2 << "  " << b1*b1-4*b0*b2 << std::endl;
  std::cout << A << "  " << B << "  " << P_A << "  " << P_B << " | " << fac_1 << "  " << fac_2 << "  " << fac_3 << "  " << fac_4 << std::endl;

  if(isReal<1.){
    if(order<2.){return fac_1.Im();}
    return (fac_1*fac_2*fac_3*fac_4).Im();
  }
  if(order<2.){return fac_1.Re();}
  return (fac_1*fac_2*fac_3*fac_4).Re();
};
Double_t TF1_RF_alpha(Double_t* x, Double_t* par){
  return RF_alpha(x[0],par[0],par[1],par[2]);
};

Double_t Ralpha_of_mu(Double_t mu, Double_t Nf, Double_t order, Double_t isReal){
  TF1 scratch("scratch",TF1_RF_alpha,0.0001,1000,3);
  scratch.SetParameter(0,Nf);
  scratch.SetParameter(1,order);
  scratch.SetParameter(2,isReal);
  Double_t g=scratch.GetX(1/mu,0.0001,1000);
  return (g*g/(4*TMath::Pi()));
};
Double_t TF1_Ralpha_of_mu(Double_t* x, Double_t* par){
  return Ralpha_of_mu(x[0],par[0],par[1],par[2]);
};

Int_t main(){
/*
  // testing how well the integration works
  TF1 scratch("scratch",TF1_neg_recip_betafn,0.0001,1000,2);
  scratch.SetParameter(0,2);
  scratch.SetParameter(1,4);
  Int_t np = 1000;
  Double_t x[1000];
  Double_t w[1000];
  scratch.CalcGaussLegendreSamplingPoints(np,x,w,1e-15);
  std::cout << scratch.IntegralFast(np,x,w,1,100) << std::endl;
  std::cout << scratch.IntegralFast(np,x,w,1,1000) << std::endl;
  std::cout << scratch.IntegralFast(np,x,w,1,10000) << std::endl;
  std::cout << scratch.Integral(1,100) << std::endl;
*/

  //chunk of code for plotting F(g) as a function of g

  TCanvas* c1=new TCanvas("c1","c1 title",700,700);
  c1->SetLogx();
  c1->SetLogy();

  TH1F* fr1=new TH1F("frame4","F(g)",1000,200,600);
  fr1->SetStats(0);
  fr1->SetMinimum(TMath::Power(10,-22) );
  fr1->SetMaximum(10000.);
  fr1->GetXaxis()->SetTitle("g");
  fr1->GetXaxis()->SetTickLength(0.02);
  fr1->GetXaxis()->SetLabelSize(0.020);
  fr1->GetXaxis()->SetLimits(0.001,10000); //minx/maxx limits
  fr1->GetYaxis()->SetTitle("F(g)");
  fr1->GetYaxis()->SetTickLength(0.02);
  fr1->GetYaxis()->SetLabelSize(0.020);
  fr1->Draw(" ");

  TF1* Fg_RNF2_NL1_Re=new TF1("test11",TF1_RF_alpha,0.001,10000,3);
  Fg_RNF2_NL1_Re->SetParameter(0,0);
  Fg_RNF2_NL1_Re->SetParameter(1,1);
  Fg_RNF2_NL1_Re->SetParameter(2,1);
  Fg_RNF2_NL1_Re->SetLineWidth(1.5);
  Fg_RNF2_NL1_Re->SetLineColor(2);
  Fg_RNF2_NL1_Re->SetLineStyle(2);
  Fg_RNF2_NL1_Re->Draw("SAME");

  TF1* Fg_RNF2_NL2_Re=new TF1("test12",TF1_RF_alpha,0.001,10000,3);
  Fg_RNF2_NL2_Re->SetParameter(0,0);
  Fg_RNF2_NL2_Re->SetParameter(1,2);
  Fg_RNF2_NL2_Re->SetParameter(2,1);
  Fg_RNF2_NL2_Re->SetLineWidth(1.5);
  Fg_RNF2_NL2_Re->SetLineColor(3);
  Fg_RNF2_NL2_Re->SetLineStyle(3);
  Fg_RNF2_NL2_Re->Draw("SAME");

  TF1* Fg_RNF2_NL3_Re=new TF1("test13",TF1_RF_alpha,0.001,10000,3);
  Fg_RNF2_NL3_Re->SetParameter(0,0);
  Fg_RNF2_NL3_Re->SetParameter(1,3);
  Fg_RNF2_NL3_Re->SetParameter(2,1);
  Fg_RNF2_NL3_Re->SetLineWidth(1.5);
  Fg_RNF2_NL3_Re->SetLineColor(4);
  Fg_RNF2_NL3_Re->SetLineStyle(4);
  Fg_RNF2_NL3_Re->Draw("SAME");

  TF1* Fg_RNF2_NL1_Im=new TF1("test11",TF1_RF_alpha,0.001,10000,3);
  Fg_RNF2_NL1_Im->SetParameter(0,0);
  Fg_RNF2_NL1_Im->SetParameter(1,1);
  Fg_RNF2_NL1_Im->SetParameter(2,0);
  Fg_RNF2_NL1_Im->SetLineWidth(1.5);
  Fg_RNF2_NL1_Im->SetLineColor(2);
  Fg_RNF2_NL1_Im->SetLineStyle(2);
  Fg_RNF2_NL1_Im->Draw("SAME");

  TF1* Fg_RNF2_NL2_Im=new TF1("test12",TF1_RF_alpha,0.001,10000,3);
  Fg_RNF2_NL2_Im->SetParameter(0,0);
  Fg_RNF2_NL2_Im->SetParameter(1,2);
  Fg_RNF2_NL2_Im->SetParameter(2,0);
  Fg_RNF2_NL2_Im->SetLineWidth(1.5);
  Fg_RNF2_NL2_Im->SetLineColor(3);
  Fg_RNF2_NL2_Im->SetLineStyle(3);
  Fg_RNF2_NL2_Im->Draw("SAME");

  TF1* Fg_RNF2_NL3_Im=new TF1("test13",TF1_RF_alpha,0.001,10000,3);
  Fg_RNF2_NL3_Im->SetParameter(0,0);
  Fg_RNF2_NL3_Im->SetParameter(1,3);
  Fg_RNF2_NL3_Im->SetParameter(2,0);
  Fg_RNF2_NL3_Im->SetLineWidth(1.5);
  Fg_RNF2_NL3_Im->SetLineColor(4);
  Fg_RNF2_NL3_Im->SetLineStyle(4);
  Fg_RNF2_NL3_Im->Draw("SAME");

  c1->SaveAs("Fg_re_and_im.pdf");

  //chunk of code for plotting alpha as a function of mu over lambda

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

  TF1* NF0_NL1=new TF1("test1",TF1_alpha_of_mu,4,1000,2);
  NF0_NL1->SetParameter(0,0);
  NF0_NL1->SetParameter(1,1);
  NF0_NL1->SetLineWidth(1.5);
  NF0_NL1->SetLineColor(2);
  NF0_NL1->SetLineStyle(2);
  NF0_NL1->Draw("SAME");

  TF1* NF0_NL2=new TF1("test2",TF1_alpha_of_mu,4,1000,2);
  NF0_NL2->SetParameter(0,0);
  NF0_NL2->SetParameter(1,2);
  NF0_NL2->SetLineWidth(1.5);
  NF0_NL2->SetLineColor(3);
  NF0_NL2->SetLineStyle(3);
  NF0_NL2->Draw("SAME");

  TF1* NF0_NL3=new TF1("test3",TF1_alpha_of_mu,4,1000,2);
  NF0_NL3->SetParameter(0,0);
  NF0_NL3->SetParameter(1,3);
  NF0_NL3->SetLineWidth(1.5);
  NF0_NL3->SetLineColor(4);
  NF0_NL3->SetLineStyle(4);
  NF0_NL3->Draw("SAME");

  TF1* NF0_NL4=new TF1("test4",TF1_alpha_of_mu,4,1000,2);
  NF0_NL4->SetParameter(0,0);
  NF0_NL4->SetParameter(1,4);
  NF0_NL4->SetLineWidth(1.5);
  NF0_NL4->SetLineColor(1);
  NF0_NL4->SetLineStyle(1);
  NF0_NL4->Draw("SAME");

  myCanv->SaveAs("NF0.pdf");

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

  TF1* RNF2_NL1_Re=new TF1("test11",TF1_Ralpha_of_mu,4,1000,3);
  RNF2_NL1_Re->SetParameter(0,2);
  RNF2_NL1_Re->SetParameter(1,1);
  RNF2_NL1_Re->SetParameter(2,1);
  RNF2_NL1_Re->SetLineWidth(1.5);
  RNF2_NL1_Re->SetLineColor(2);
  RNF2_NL1_Re->SetLineStyle(2);
  RNF2_NL1_Re->Draw("SAME");

  TF1* RNF2_NL2_Re=new TF1("test12",TF1_Ralpha_of_mu,4,1000,3);
  RNF2_NL2_Re->SetParameter(0,2);
  RNF2_NL2_Re->SetParameter(1,2);
  RNF2_NL2_Re->SetParameter(2,1);
  RNF2_NL2_Re->SetLineWidth(1.5);
  RNF2_NL2_Re->SetLineColor(3);
  RNF2_NL2_Re->SetLineStyle(3);
  RNF2_NL2_Re->Draw("SAME");

  TF1* RNF2_NL3_Re=new TF1("test13",TF1_Ralpha_of_mu,4,1000,3);
  RNF2_NL3_Re->SetParameter(0,2);
  RNF2_NL3_Re->SetParameter(1,3);
  RNF2_NL3_Re->SetParameter(2,1);
  RNF2_NL3_Re->SetLineWidth(1.5);
  RNF2_NL3_Re->SetLineColor(4);
  RNF2_NL3_Re->SetLineStyle(4);
  RNF2_NL3_Re->Draw("SAME");

  myCanv3->SaveAs("RNF2_Re.pdf");

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

  TF1* RNF2_NL1_Im=new TF1("test11",TF1_Ralpha_of_mu,4,1000,3);
  RNF2_NL1_Im->SetParameter(0,2);
  RNF2_NL1_Im->SetParameter(1,1);
  RNF2_NL1_Im->SetParameter(2,0);
  RNF2_NL1_Im->SetLineWidth(1.5);
  RNF2_NL1_Im->SetLineColor(2);
  RNF2_NL1_Im->SetLineStyle(2);
  RNF2_NL1_Im->Draw("SAME");

  TF1* RNF2_NL2_Im=new TF1("test12",TF1_Ralpha_of_mu,4,1000,3);
  RNF2_NL2_Im->SetParameter(0,2);
  RNF2_NL2_Im->SetParameter(1,2);
  RNF2_NL2_Im->SetParameter(2,0);
  RNF2_NL2_Im->SetLineWidth(1.5);
  RNF2_NL2_Im->SetLineColor(3);
  RNF2_NL2_Im->SetLineStyle(3);
  RNF2_NL2_Im->Draw("SAME");

  TF1* RNF2_NL3_Im=new TF1("test13",TF1_Ralpha_of_mu,4,1000,3);
  RNF2_NL3_Im->SetParameter(0,2);
  RNF2_NL3_Im->SetParameter(1,3);
  RNF2_NL3_Im->SetParameter(2,0);
  RNF2_NL3_Im->SetLineWidth(1.5);
  RNF2_NL3_Im->SetLineColor(4);
  RNF2_NL3_Im->SetLineStyle(4);
  RNF2_NL3_Im->Draw("SAME");

  myCanv4->SaveAs("RNF2_Im.pdf");

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
  return 0;

}
