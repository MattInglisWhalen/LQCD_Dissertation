//C++ standard headers
#include <iostream>
#include <climits>
#include <vector>

//ROOT headers for plotting
#include "TCanvas.h"
#include "TH1.h"
#include "TF1.h"
#include "TGraph.h"
#include "TLegend.h"

//ROOT headers for math
#include "TMath.h"
#include "TComplex.h"

//ROOT Minuit headers for minimization
#include "TFitter.h"
#include "TVirtualFitter.h"
#include "TMinuit.h"

//in-house header files
#include "BetaFn.h"

using namespace std;

#define ORDER_P1 5

#define NPOINT 20
#define STEPSIZE 0.025

#define NOTPADE 0
#define PADE 1

#define GAMMA_E 0.5772156649

Double_t xx1[NPOINT],yy1[NPOINT];
Double_t xx2[NPOINT],yy2[NPOINT];
Double_t xx3[NPOINT],yy3[NPOINT];
Double_t xx4[NPOINT],yy4[NPOINT];

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
    term2=( bf[1]/(2.*bf[0]*bf[0]) )*TComplex::Log( bf[0]*x );
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

Double_t a1_qq(Int_t nflav,Int_t var){
  Double_t prefac;
  if(var==2){prefac=TMath::Power(4.*TMath::Pi(),-2.);}
  else{prefac=TMath::Power(4.*TMath::Pi(),-1.);}
  Double_t coeff0=31./3.;
  Double_t coeff1=-10./9.;
  return prefac*(coeff0+coeff1*nflav);
}

Double_t a2_qq(Int_t nflav,Int_t var){
  Double_t prefac;
  if(var==2){prefac=TMath::Power(4.*TMath::Pi(),-4.);}
  else{prefac=TMath::Power(4.*TMath::Pi(),-2.);}
  Double_t coeff0=(4343./18. + 36.*TMath::Power(TMath::Pi(),2.) + 9.*TMath::Power(TMath::Pi(),4.)/4. + 66.*RZ3);
  Double_t coeff1=-(1229./27.+53.*RZ3);
  Double_t coeff2=100./81.;
  return prefac*(coeff0+coeff1*nflav+coeff2*nflav*nflav);
};

Double_t t1_qq(Int_t nflav,Int_t var){
  BetaFn bf(4,nflav,var);
  return (2.*bf[0]*(1.-GAMMA_E)-a1_qq(nflav,var) );
}

Double_t t2_qq(Int_t nflav,Int_t var){
  BetaFn bf(4,nflav,var);
  return (bf[0]*bf[0]*(4.-TMath::Power(TMath::Pi(),2.)/3.)+2.*bf[1]*(1.-GAMMA_E)+a1_qq(nflav,var)*a1_qq(nflav,var)-a2_qq(nflav,var) );
}

// transfer coefficients for the qqbar scheme
// be careful though, because you're expanding g_qq in powers of g_ms
// rather than the usual way
Double_t t1_qq_new(Int_t nflav,Int_t var){
  Double_t prefac=-1./TMath::Power(4.*TMath::Pi(),2.);
  Double_t coeff0=22.*(GAMMA_E-35./66.);
  Double_t coeff1=-4.*(GAMMA_E-1./6.)/3.;
  if(var==1){return TMath::Power(4.*TMath::Pi(),1.)*prefac*(coeff0+coeff1*nflav);}
  return prefac*(coeff0+coeff1*nflav);
}
Double_t t2_qq_new(Int_t nflav,Int_t var){
  Double_t prefac=1./TMath::Power(4.*TMath::Pi(),4.);
  Double_t coeff0=1107./2.-204.*GAMMA_E-229.*TMath::Power(TMath::Pi(),2.)/3.+9.*TMath::Power(TMath::Pi(),4.)/4.-66.*RZ3;
  Double_t coeff1=(-553./3.+76.*GAMMA_E+44.*TMath::Power(TMath::Pi(),2.)/3.+52.*RZ3)/3.;
  Double_t coeff2=4.*(12.-TMath::Power(TMath::Pi(),2.))/27.;
  if(var==1){return TMath::Power(4.*TMath::Pi(),2.)*prefac*(coeff0+coeff1*nflav+coeff2*nflav*nflav);}
  return prefac*(coeff0+coeff1*nflav+coeff2*nflav*nflav);
}

// beta function coefficients for the qqbar scheme
Double_t b2_qq(Int_t nflav,Int_t var){
  BetaFn bf_ms(4,nflav,var);  //nloop,nf,var
  return bf_ms[2]+bf_ms[1]*t1_qq(nflav,var)-bf_ms[0]*t2_qq(nflav,var);
};
Double_t b3_qq(Int_t nflav, Int_t var, Int_t pade){
  BetaFn bf_ms(4,nflav,var);
  if(pade==NOTPADE){return 0.;}
  return b2_qq(nflav,var)*b2_qq(nflav,var)/bf_ms[1];;
};

Double_t rLambda_qq(Double_t x,Int_t nflav,Int_t nloop,Int_t var){
  BetaFn bf_qq(nloop,nflav,var);
  if(nloop>=3){bf_qq[2]=b2_qq(nflav,var);}
  if(nloop>=4){bf_qq[3]=b3_qq(nflav,var,1);}
  return 1./mu_over_lambda(x,bf_qq);
};
Double_t rLambda_ms(Double_t x,Int_t nflav,Int_t nloop,Int_t var){
  BetaFn bf_ms(nloop,nflav,var);
  return TMath::Exp( t1_qq(nflav,var)/(2.*bf_ms[0]) )*rLambda_qq(x,nflav,nloop,var);
};

Double_t Lambda2_over_Lambda0(Double_t x,Int_t nloop,Int_t var){
  return rLambda_ms(x,2,nloop,var)/rLambda_ms(x,0,nloop,var);
};
Double_t TF1_Lambda2_over_Lambda0(Double_t* x, Double_t* par){
  return Lambda2_over_Lambda0(x[0],par[0],par[1]);
};

Double_t Lambda3_over_Lambda2(Double_t l2_over_l0,Int_t nloop,Int_t var){
  TF1 scratch("scratch",TF1_Lambda2_over_Lambda0,-1.,1.,2);
  scratch.SetParameter(0,nloop);
  scratch.SetParameter(1,var);
  Double_t x;
  if(var>1){x=scratch.GetX(l2_over_l0,0.01,2000.);}
  else{
    if(l2_over_l0<Lambda2_over_Lambda0(10000.,nloop,var)){
      x=scratch.GetX(l2_over_l0,0.001,100.);
    }
    else{
      x=scratch.GetX(l2_over_l0,-2000.,-0.001);
    }
  }
  return rLambda_ms(x,3,nloop,var)/rLambda_ms(x,2,nloop,var);
};

Double_t fit_pnomial(Double_t x,Double_t* p){
//  return p[0]+p[1]*TMath::Power(x,p[2])+p[3]*TMath::Power(x,p[4]);
  Double_t xx=x-1.1;
  Double_t y=0.;
  for(Int_t i=0;i<ORDER_P1;++i){
    y+=p[i]*TMath::Power(xx,i);
  }
  return y;
};

Double_t TF1_Lambda3_over_Lambda2(Double_t* x, Double_t* par){
  if( x[0]<Lambda2_over_Lambda0(.001,par[0],par[1]) ){
    return Lambda3_over_Lambda2(x[0],par[0],par[1]);
  }
  Double_t parnew[ORDER_P1];
  for(Int_t i=0;i<ORDER_P1;++i){
    parnew[i]=par[i+2];
  }
  return fit_pnomial(x[0],parnew);
};
// Minimize the r^2 value, should be exact because I have the exact function
void minuit_1loop(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t r2=0.0;
  for(Int_t i=0;i<NPOINT;++i) { 
    r2 += TMath::Power( yy1[i] - fit_pnomial(xx1[i],p) , 2. );
  }
  fval = r2;
}
// Minimize the r^2 value, should be exact because I have the exact function
void minuit_2loop(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t r2=0.0;
  for(Int_t i=0;i<NPOINT;++i) { 
    r2 += TMath::Power( yy2[i] - fit_pnomial(xx2[i],p) , 2. );
  }
  fval = r2;
}
// Minimize the r^2 value, should be exact because I have the exact function
void minuit_3loop(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t r2=0.0;
  for(Int_t i=0;i<NPOINT;++i) { 
    r2 += TMath::Power( yy3[i] - fit_pnomial(xx3[i],p) , 2. );
  }
  fval = r2;
}
// Minimize the r^2 value, should be exact because I have the exact function
void minuit_4loop(Int_t& /*nPar*/, Double_t* /*grad*/ , Double_t& fval, Double_t* p, Int_t /*iflag */  ){
  Double_t r2=0.0;
  for(Int_t i=0;i<NPOINT;++i) { 
    r2 += TMath::Power( yy4[i] - fit_pnomial(xx4[i],p) , 2. );
  }
  fval = r2;
}

Int_t main(){

  Int_t variable=1;

  // Minimize the statistics with Minuit and thus find the best estimator for x_0
  TVirtualFitter::SetDefaultFitter("Minuit");
  Double_t arglist[100]; //  for(Int_t i=0;i<100;++i){arglist[i]=0.;}

  Double_t offset1=0.;
  for(Int_t i=0;i<NPOINT;++i){
    xx1[i]=rLambda_ms(offset1+STEPSIZE+STEPSIZE*i,2,1,variable)/rLambda_ms(offset1+STEPSIZE+STEPSIZE*i,0,1,variable);
    yy1[i]=rLambda_ms(offset1+STEPSIZE+STEPSIZE*i,3,1,variable)/rLambda_ms(offset1+STEPSIZE+STEPSIZE*i,2,1,variable);
    cout << xx1[i] << "  ,  " << yy1[i] << endl;
  }
// Uncomment if you want to use negative alpha points to fit
//  for(Int_t i=0;i<300;++i){
//    xx1[i]=rLambda_ms(-1-STEPSIZE*i,2,1,variable)/rLambda_ms(-1-STEPSIZE*i,0,1,variable);
//    yy1[i]=rLambda_ms(-1-STEPSIZE*i,3,1,variable)/rLambda_ms(-1-STEPSIZE*i,2,1,variable);
//  }

  for(Int_t i=0;i<NPOINT;++i){
    xx2[i]=rLambda_ms(STEPSIZE+STEPSIZE*i,2,2,variable)/rLambda_ms(STEPSIZE+STEPSIZE*i,0,2,variable);
    yy2[i]=rLambda_ms(STEPSIZE+STEPSIZE*i,3,2,variable)/rLambda_ms(STEPSIZE+STEPSIZE*i,2,2,variable);
  }
  Double_t offset3=0.;
  for(Int_t i=0;i<NPOINT;++i){
    xx3[i]=rLambda_ms(offset3+STEPSIZE+STEPSIZE*i,2,3,variable)/rLambda_ms(offset3+STEPSIZE+STEPSIZE*i,0,3,variable);
    yy3[i]=rLambda_ms(offset3+STEPSIZE+STEPSIZE*i,3,3,variable)/rLambda_ms(offset3+STEPSIZE+STEPSIZE*i,2,3,variable);
  }
  Double_t offset4=0.;
  for(Int_t i=0;i<NPOINT;++i){
    xx4[i]=rLambda_ms(offset4+STEPSIZE+STEPSIZE*i,2,4,variable)/rLambda_ms(offset4+STEPSIZE+STEPSIZE*i,0,4,variable);
    yy4[i]=rLambda_ms(offset4+STEPSIZE+STEPSIZE*i,3,4,variable)/rLambda_ms(offset4+STEPSIZE+STEPSIZE*i,2,4,variable);
  }

  // --------------------------------------------

  TVirtualFitter* minuit1 = TVirtualFitter::Fitter(0,ORDER_P1);
  for(Int_t i=0;i<ORDER_P1;++i){
    minuit1->SetParameter(i,"par",0.0, 0.01,-50., 50.);
  }
  minuit1->SetFCN(minuit_1loop);
  // minimize
  minuit1->ExecuteCommand("MIGRAD",arglist,0);
  minuit1->ExecuteCommand("MINOS",arglist,0);
  // local variables
  Double_t pval1[ORDER_P1];
  // Fit results
  for(Int_t i=0;i<ORDER_P1;++i){
    pval1[i] = minuit1->GetParameter(i);
  }
  cout << "1-loop ";
  for(Int_t i=0;i<ORDER_P1;++i){
    cout << pval1[i] << "  " ;
  }
  cout << endl;

  // --------------------------------------------

  TVirtualFitter* minuit2 = TVirtualFitter::Fitter(0,ORDER_P1);
  for(Int_t i=0;i<ORDER_P1;++i){
    minuit2->SetParameter(i,"par",0.0, 0.01,-50., 50.);
  }
  minuit2->SetFCN(minuit_2loop);
  // minimize
  minuit2->ExecuteCommand("MIGRAD",arglist,0);
  minuit2->ExecuteCommand("MINOS",arglist,0);
  // local variables
  Double_t pval2[ORDER_P1];
  // Fit results
  for(Int_t i=0;i<ORDER_P1;++i){
    pval2[i] = minuit2->GetParameter(i);
  }
  cout << "2-loop ";
  for(Int_t i=0;i<ORDER_P1;++i){
    cout << pval2[i] << "  " ;
  }
  cout << endl;

  // --------------------------------------------

  TVirtualFitter* minuit3 = TVirtualFitter::Fitter(0,ORDER_P1);
  for(Int_t i=0;i<ORDER_P1;++i){
    minuit3->SetParameter(i,"par",0.0, 0.01,-50., 50.);
  }
  minuit3->SetFCN(minuit_3loop);
  // minimize
  minuit3->ExecuteCommand("MIGRAD",arglist,0);
  minuit3->ExecuteCommand("MINOS",arglist,0);
  // local variables
  Double_t pval3[ORDER_P1];
  // Fit results
  for(Int_t i=0;i<ORDER_P1;++i){
    pval3[i] = minuit3->GetParameter(i);
  }
  cout << "3-loop ";
  for(Int_t i=0;i<ORDER_P1;++i){
    cout << pval3[i] << "  " ;
  }
  cout << endl;

  // --------------------------------------------
/*
  TVirtualFitter* minuit4 = TVirtualFitter::Fitter(0,ORDER_P1);
  for(Int_t i=0;i<ORDER_P1;++i){
    minuit4->SetParameter(i,"par",0.0, 0.01,-50., 50.);
  }
  minuit4->SetFCN(minuit_4loop);
  // minimize
  minuit4->ExecuteCommand("MIGRAD",arglist,0);
  minuit4->ExecuteCommand("MINOS",arglist,0);
  // local variables
  Double_t pval4[ORDER_P1];
  // Fit results
  for(Int_t i=0;i<ORDER_P1;++i){
    pval4[i] = minuit4->GetParameter(i);
  }
  cout << "4-loop ";
  for(Int_t i=0;i<ORDER_P1;++i){
    cout << pval4[i] << "  " ;
  }
  cout << endl;
*/
  // --------------------------------------------

  TCanvas* canv = new TCanvas("canv","canv",700,700);
//  canv->SetLogx();

  Double_t minxwindow=0.85;
  Double_t maxxwindow=1.15;

  TH1F* frame=new TH1F("frame","",1000,200,600);
  frame->SetStats(0);
  frame->SetMinimum(0.85);
  frame->SetMaximum(1.15);
  frame->GetXaxis()->SetTitle("#Lambda_{2} / #Lambda_{0}");
  frame->GetXaxis()->SetTickLength(0.02);
  frame->GetXaxis()->SetLabelSize(0.020);
  frame->GetXaxis()->SetLimits(minxwindow,maxxwindow); //minx/maxx limits
  frame->GetYaxis()->SetTitle("#Lambda_{3} / #Lambda_{2}");
  frame->GetYaxis()->SetTickLength(0.02);
  frame->GetYaxis()->SetLabelSize(0.020);
  frame->Draw(" ");

  TF1* func_nl1=new TF1("func_nl1",TF1_Lambda3_over_Lambda2,minxwindow,maxxwindow,2+ORDER_P1);
  func_nl1->SetParameter(0,1);  //nloop
  func_nl1->SetParameter(1,variable);  //var
  for(Int_t i=0;i<(ORDER_P1);++i){
    func_nl1->SetParameter(i+2,pval1[i]);
  }
  func_nl1->SetLineWidth(1.5);
  func_nl1->SetLineColor(2);
  func_nl1->SetLineStyle(2);
  func_nl1->Draw("SAME");

  TGraph* gr1=new TGraph(NPOINT,xx1,yy1);
  gr1->SetMarkerStyle(5);
  gr1->SetMarkerColor(2);
  gr1->Draw("P");

  TF1* func_nl2=new TF1("func_nl2",TF1_Lambda3_over_Lambda2,minxwindow,maxxwindow,2+ORDER_P1);
  func_nl2->SetParameter(0,2);  //nloop
  func_nl2->SetParameter(1,variable);  //var
  for(Int_t i=0;i<(ORDER_P1);++i){
    func_nl2->SetParameter(i+2,pval2[i]);
  }
  func_nl2->SetLineWidth(1.5);
  func_nl2->SetLineColor(3);
  func_nl2->SetLineStyle(3);
  func_nl2->Draw("SAME");

  TGraph* gr2=new TGraph(NPOINT,xx2,yy2);
  gr2->SetMarkerStyle(5);
  gr2->SetMarkerColor(3);
  gr2->Draw("P");


  TF1* func_nl3=new TF1("func_nl3",TF1_Lambda3_over_Lambda2,minxwindow,maxxwindow,2+ORDER_P1);
  func_nl3->SetParameter(0,3);  //nloop
  func_nl3->SetParameter(1,variable);  //var
  for(Int_t i=0;i<(ORDER_P1);++i){
    func_nl3->SetParameter(i+2,pval3[i]);
  }
  func_nl3->SetLineWidth(1.5);
  func_nl3->SetLineColor(4);
  func_nl3->SetLineStyle(4);
  func_nl3->Draw("SAME");

  TGraph* gr3=new TGraph(NPOINT,xx3,yy3);
  gr3->SetMarkerStyle(5);
  gr3->SetMarkerColor(4);
  gr3->Draw("P");

/*
  TF1* func_nl4=new TF1("func_nl4",TF1_Lambda3_over_Lambda2,minxwindow,maxxwindow,2+ORDER_P1);
  func_nl4->SetParameter(0,4);  //nloop
  func_nl4->SetParameter(1,variable);  //var
  for(Int_t i=0;i<(ORDER_P1);++i){
    func_nl4->SetParameter(i+2,pval4[i]);
  }
  func_nl4->SetLineWidth(1.5);
  func_nl4->SetLineColor(1);
  func_nl4->SetLineStyle(1);
  func_nl4->Draw("SAME");
*/

  TLegend* leg = new TLegend(0.15,0.65,0.55,0.85,"Loops");
  leg->SetBorderSize(0);
  leg->SetFillStyle(0); 
  leg->AddEntry(func_nl1,"1","l");
  leg->AddEntry(func_nl2,"2","l");
  leg->AddEntry(func_nl3,"3","l");
//  leg->AddEntry(func_nl4,"4","l");
  leg->Draw();

  Double_t centralx=0.75/0.62;
//  Double_t centralx=0.6803/0.6159;  // = 1.10456
  Double_t xerr=0.04263;
  Double_t meany,yerrp,yerrn,yerr;
  meany=(func_nl1->Eval(centralx)+func_nl2->Eval(centralx)+func_nl3->Eval(centralx))/3.;
  yerrp=meany-(func_nl1->Eval(centralx+xerr)+func_nl2->Eval(centralx+xerr)+func_nl3->Eval(centralx+xerr))/3.;
  yerrn=meany-(func_nl1->Eval(centralx-xerr)+func_nl2->Eval(centralx-xerr)+func_nl3->Eval(centralx-xerr))/3.;
  yerr=TMath::Max(TMath::Abs(yerrp),TMath::Abs(yerrn));
  Double_t sigma = TMath::Sqrt((TMath::Power(meany-func_nl1->Eval(centralx),2.)+TMath::Power(meany-func_nl2->Eval(centralx),2.)+TMath::Power(meany-func_nl3->Eval(centralx),2.))/2.);

  cout << "meany = 1/3( " <<  func_nl1->Eval(centralx) << "+"<< func_nl2->Eval(centralx) << "+"<<   func_nl3->Eval(centralx) << "=" 
  << meany << "+/-" << yerr << "+/-" << sigma << endl;

  canv->SaveAs("rat32_v_rat20_good.pdf");

  cout << endl << endl << t1_qq(0,variable) << "  vs old value  " << t1_qq_new(0,variable) << endl;
  cout << t2_qq(0,variable) << "  vs new value  " << t2_qq_new(0,variable) << endl;


  cout << endl << endl << t1_qq(2,variable) << "  vs old value  " << t1_qq_new(2,variable) << endl;
  cout << t2_qq(2,variable) << "  vs new value  " << t2_qq_new(2,variable) << endl;


  cout << endl << endl << t1_qq(3,variable) << "  vs old value  " << t1_qq_new(3,variable) << endl;
  cout << t2_qq(3,variable) << "  vs new value  " << t2_qq_new(3,variable) << endl << endl << endl;

  BetaFn bf0(4,0,variable); BetaFn bf2(4,2,variable);
  cout << TMath::Exp(t1_qq(2,variable)/(2.*bf2[0]) - t1_qq(0,variable)/(2.*bf0[0]) ) << endl;
  cout << TMath::Exp(t1_qq_new(2,variable)/(2.*bf2[0]) - t1_qq_new(0,variable)/(2.*bf0[0]) ) << endl;

  cout << Lambda2_over_Lambda0(0.5,1,1) << endl;
  cout << Lambda2_over_Lambda0(0.5,2,1) << endl;
  cout << Lambda2_over_Lambda0(0.5,3,1) << endl;




  return 0;

};


