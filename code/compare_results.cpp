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
#include "TLatex.h"

using namespace std;

#define N0 7
#define N2 10
#define N3 6

#define NN02 8
#define NN3 5
#define NN4 1

Int_t main(){

  //  ----------------------------------------------------------------------
  //
  //             COMPARING LAMBDA VALUES
  //
  //  ----------------------------------------------------------------------

  Double_t offsett=0.25;

  TCanvas* canvlambda = new TCanvas("canvlambda","canvlambda",700,700);
  
  TH1F* framelambda=new TH1F("framelambda"," ",1000,200,600);
  framelambda->SetStats(0);
  framelambda->SetMinimum(-1);
  framelambda->SetMaximum(N0+N2+N3+2);
  framelambda->GetXaxis()->SetTitle("r_{0}#Lambda");
  framelambda->GetXaxis()->SetTickLength(0.02);
  framelambda->GetXaxis()->SetLabelSize(0.020);
  framelambda->GetXaxis()->SetLimits(0.4,1.3); //minx/maxx limits
  framelambda->GetYaxis()->SetTitle(" ");
  framelambda->GetYaxis()->SetTickLength(0.);
  framelambda->GetYaxis()->SetLabelSize(0.);
  framelambda->Draw(" ");

  Double_t lambda0vals[N0]; Double_t lambda0errs[N0];
  Double_t lambda2vals[N2]; Double_t lambda2errs[N2];
  Double_t lambda3vals[N3]; Double_t lambda3errs[N3];

  Double_t y0[N0]; Double_t y0e[N0];
  Double_t y2[N2]; Double_t y2e[N2];
  Double_t y3[N3]; Double_t y3e[N3];

  // Lambda 0 -----------------------------

  TLatex latt0_0(1.05,N0-1-offsett,"#it{Inglis-Whalen 2014}"); latt0_0.SetTextSize(0.02); latt0_0.Draw("SAME");
  lambda0vals[0]=0.616; lambda0errs[0]=0.033; y0[0]=N0-1; y0e[0]=0.;
  TLatex latt0_1(1.05,N0-2-offsett,"Sternbeck #it{et al.} 2010"); latt0_1.SetTextSize(0.02); latt0_1.Draw("SAME");
  lambda0vals[1]=0.620; lambda0errs[1]=0.010; y0[1]=N0-2; y0e[1]=0.;
  TLatex latt0_2(1.05,N0-3-offsett,"Brambilla #it{et al.} 2010"); latt0_2.SetTextSize(0.02); latt0_2.Draw("SAME");
  lambda0vals[2]=0.637; lambda0errs[2]=0.032; y0[2]=N0-3; y0e[2]=0.;
  TLatex latt0_3(1.05,N0-4-offsett,"Boucaud #it{et al.} 2009"); latt0_3.SetTextSize(0.02); latt0_3.Draw("SAME");
  lambda0vals[3]=0.590; lambda0errs[3]=0.022; y0[3]=N0-4; y0e[3]=0.;

  TLatex latt0_4(1.05,N0-5-offsett,"Bali #it{et al.} 1992"); latt0_4.SetTextSize(0.02); latt0_4.Draw("SAME");
  lambda0vals[4]=0.661; lambda0errs[4]=0.027; y0[4]=N0-5; y0e[4]=0.;
  TLatex latt0_5(1.05,N0-6-offsett,"UKQCD 1992"); latt0_5.SetTextSize(0.02); latt0_5.Draw("SAME");
  lambda0vals[5]=0.686; lambda0errs[5]=0.054; y0[5]=N0-6; y0e[5]=0.;
  TLatex latt0_6(1.05,N0-7-offsett,"El-Khadra #it{et al.} 1992"); latt0_6.SetTextSize(0.02); latt0_6.Draw("SAME");
  lambda0vals[6]=0.593; lambda0errs[6]=0.025; y0[6]=N0-7; y0e[6]=0.;

  TLatex latnff0(0.390,1.5,"n_{f} = 0"); latnff0.SetTextSize(0.03); latnff0.SetTextAngle(90); latnff0.Draw("SAME");

  TGraphErrors* graphh0=new TGraphErrors(N0,lambda0vals,y0,lambda0errs,y0e);
  graphh0->SetMarkerColor(1);
  graphh0->SetMarkerStyle(33);
  graphh0->SetMarkerSize(1);
  graphh0->Draw("P");

  Double_t linex0[4]={0.60,0.60,0.64,0.64};
  Double_t liney0[4]={-2.,N0,N0,-2.};
  TGraph* linee0=new TGraph(4,linex0,liney0);
  linee0->SetLineColor(1);
  linee0->SetLineStyle(3);
  linee0->Draw("L");

  Double_t linex0x[2]={0.,2.};
  Double_t liney0y[2]={N0,N0};
  TGraph* line0e=new TGraph(2,linex0x,liney0y);
  line0e->SetLineColor(1);
  line0e->SetLineStyle(1);
  line0e->Draw("L");

  // Lambda 2 -----------------------------

  TLatex latt2_0(1.05,N0+N2-offsett,"#it{Inglis-Whalen 2014}"); latt2_0.SetTextSize(0.02); latt2_0.Draw("SAME");
  lambda2vals[0]=0.680; lambda2errs[0]=0.026; y2[0]=N0+N2; y2e[0]=0.;
  TLatex latt2_1(1.05,N0+N2-1-offsett,"Karbstein #it{et al.} 2014"); latt2_1.SetTextSize(0.02); latt2_1.Draw("SAME");
  lambda2vals[1]=0.692; lambda2errs[1]=0.031; y2[1]=N0+N2-1; y2e[1]=0.;
  TLatex latt2_2(1.05,N0+N2-2-offsett,"ALPHA 2012"); latt2_2.SetTextSize(0.02); latt2_2.Draw("SAME");
  lambda2vals[2]=0.789; lambda2errs[2]=0.052; y2[2]=N0+N2-2; y2e[2]=0.;

  TLatex latt2_3(1.05,N0+N2-3-offsett,"ETM 2011"); latt2_3.SetTextSize(0.02); latt2_3.Draw("SAME");
  lambda2vals[3]=0.658; lambda2errs[3]=0.055; y2[3]=N0+N2-3; y2e[3]=0.;
  TLatex latt2_4(1.05,N0+N2-4-offsett,"Sternbeck #it{et al.} 2010"); latt2_4.SetTextSize(0.02); latt2_4.Draw("SAME");
  lambda2vals[4]=0.600; lambda2errs[4]=0.036; y2[4]=N0+N2-4; y2e[4]=0.;
  TLatex latt2_5(1.05,N0+N2-5-offsett,"ETM 2010"); latt2_5.SetTextSize(0.02); latt2_5.Draw("SAME");
  lambda2vals[5]=0.720; lambda2errs[5]=0.050; y2[5]=N0+N2-5; y2e[5]=0.;

  TLatex latt2_6(1.05,N0+N2-6-offsett,"JLQCD/TWQCD 2008"); latt2_6.SetTextSize(0.02); latt2_6.Draw("SAME");
  lambda2vals[6]=0.581; lambda2errs[6]=0.046; y2[6]=N0+N2-6; y2e[6]=0.;
  TLatex latt2_7(1.05,N0+N2-7-offsett,"QCDSF/UKQCD 2005"); latt2_7.SetTextSize(0.02); latt2_7.Draw("SAME");
  lambda2vals[7]=0.617; lambda2errs[7]=0.045; y2[7]=N0+N2-7; y2e[7]=0.;
  TLatex latt2_8(1.05,N0+N2-8-offsett,"ALPHA 2004"); latt2_8.SetTextSize(0.02); latt2_8.Draw("SAME");
  lambda2vals[8]=0.620; lambda2errs[8]=0.03; y2[8]=N0+N2-8; y2e[8]=0.;
  TLatex latt2_9(1.05,N0+N2-9-offsett,"Boucaud #it{et al.} 2001"); latt2_9.SetTextSize(0.02); latt2_9.Draw("SAME");
  lambda2vals[9]=0.669; lambda2errs[9]=0.069; y2[9]=N0+N2-9; y2e[9]=0.;

  TLatex latnff2(0.390,10.5,"n_{f} = 2"); latnff2.SetTextSize(0.03); latnff2.SetTextAngle(90); latnff2.Draw("SAME");

  TGraphErrors* graphh2=new TGraphErrors(N2,lambda2vals,y2,lambda2errs,y2e);
  graphh2->SetMarkerColor(1);
  graphh2->SetMarkerStyle(33);
  graphh2->SetMarkerSize(1);
  graphh2->Draw("P");

  Double_t linex2[4]={0.66,0.66,0.84,0.84};
  Double_t liney2[4]={N0,N0+N2+1,N0+N2+1,N0};
  TGraph* linee2=new TGraph(4,linex2,liney2);
  linee2->SetLineColor(1);
  linee2->SetLineStyle(3);
  linee2->Draw("L");

  Double_t linex2x[2]={0.,2.};
  Double_t liney2y[2]={N0+N2+1,N0+N2+1};
  TGraph* line2e=new TGraph(2,linex2x,liney2y);
  line2e->SetLineColor(1);
  line2e->SetLineStyle(1);
  line2e->Draw("L");

  // Lambda 3 -----------------------------

  TLatex latt3_0(1.05,N0+N2+N3+1-offsett,"#it{Inglis-Whalen 2014}"); latt3_0.SetTextSize(0.02); latt3_0.Draw("SAME");
  lambda3vals[0]=0.727; lambda3errs[0]=0.032; y3[0]=N0+N2+N3+1; y3e[0]=0.;
  TLatex latt3_1(1.05,N0+N2+N3-offsett,"Bazavov #it{et al.} 2012"); latt3_1.SetTextSize(0.02); latt3_1.Draw("SAME");
  lambda3vals[1]=0.720; lambda3errs[1]=0.070; y3[1]=N0+N2+N3; y3e[1]=0.;
  TLatex latt3_2(1.05,N0+N2+N3-1-offsett,"HPQCD 2010"); latt3_2.SetTextSize(0.02); latt3_2.Draw("SAME");
  lambda3vals[2]=0.812; lambda3errs[2]=0.022; y3[2]=N0+N2+N3-1; y3e[2]=0.;
  TLatex latt3_3(1.05,N0+N2+N3-2-offsett,"PACS-CS 2009"); latt3_3.SetTextSize(0.02); latt3_3.Draw("SAME");
  lambda3vals[3]=0.888; lambda3errs[3]=0.074; y3[3]=N0+N2+N3-2; y3e[3]=0.;
  TLatex latt3_4(1.05,N0+N2+N3-3-offsett,"Maltman #it{et al.} 2008"); latt3_4.SetTextSize(0.02); latt3_4.Draw("SAME");
  lambda3vals[4]=0.841; lambda3errs[4]=0.040; y3[4]=N0+N2+N3-3; y3e[4]=0.;
  TLatex latt3_5(1.05,N0+N2+N3-4-offsett,"HPQCD 2008"); latt3_5.SetTextSize(0.02); latt3_5.Draw("SAME");
  lambda3vals[5]=0.763; lambda3errs[5]=0.042; y3[5]=N0+N2+N3-4; y3e[5]=0.;

  TLatex latnff3(0.390,20.0,"n_{f} = 3"); latnff3.SetTextSize(0.03); latnff3.SetTextAngle(90); latnff3.Draw("SAME");

  TGraphErrors* graphh3=new TGraphErrors(N3,lambda3vals,y3,lambda3errs,y3e);
  graphh3->SetMarkerColor(1);
  graphh3->SetMarkerStyle(33);
  graphh3->SetMarkerSize(1);
  graphh3->Draw("P");

  Double_t linex3[4]={0.77,0.77,0.85,0.85};
  Double_t liney3[4]={N0+N2+1,N0+N2+N3+2,N0+N2+N3+2,N0+N2+1};
  TGraph* linee3=new TGraph(4,linex3,liney3);
  linee3->SetLineColor(1);
  linee3->SetLineStyle(3);
  linee3->Draw("L");

  // -------------------------------------

  canvlambda->SaveAs("compare_lambdas.pdf");

  //  ----------------------------------------------------------------------
  //
  //             COMPARING ALPHA VALUES
  //
  //  ----------------------------------------------------------------------

  TCanvas* canvalpha = new TCanvas("canvalpha","canvalpha",700,700);
  
  TH1F* framealpha=new TH1F("framealpha"," ",1000,200,600);
  framealpha->SetStats(0);
  framealpha->SetMinimum(0.25);
  framealpha->SetMaximum(1.);
  framealpha->GetXaxis()->SetTitle("#alpha(m_{Z})");
  framealpha->GetXaxis()->SetTickLength(0.02);
  framealpha->GetXaxis()->SetLabelSize(0.02);
  framealpha->GetXaxis()->SetLimits(0.104,0.135); //minx/maxx limits
  framealpha->GetYaxis()->SetTitle(" ");
  framealpha->GetYaxis()->SetTickLength(0.);
  framealpha->GetYaxis()->SetLabelSize(0.);
  framealpha->Draw(" ");

  Double_t alpha02vals[NN02]; Double_t alpha02errs[NN02];
  Double_t alpha3vals[NN3]; Double_t alpha3errs[NN3];
  Double_t alpha4vals[NN4]; Double_t alpha4errs[NN4];

  Double_t yy02[NN02]; Double_t yy02e[NN02];
  Double_t yy3[NN3]; Double_t yy3e[NN3];
  Double_t yy4[NN4]; Double_t yy4e[NN4];

  Double_t offset=0.01;
  // alpha 02 -----------------------------

  TLatex lat02_0(0.124,0.55-offset,"#it{Inglis-Whalen 2014}"); lat02_0.SetTextSize(0.025); lat02_0.Draw("SAME");
  alpha02vals[0]=0.1146; alpha02errs[0]=0.0013; yy02[0]=0.55; yy02e[0]=0.;
  TLatex lat02_1(0.124,0.50-offset,"QCDSF/UKQCD 2005"); lat02_1.SetTextSize(0.025); lat02_1.Draw("SAME");
  alpha02vals[1]=0.1120; alpha02errs[1]=0.0022; yy02[1]=0.50; yy02e[1]=0.;
  TLatex lat02_2(0.124,0.45-offset,"Boucaud #it{et al.} 2001"); lat02_2.SetTextSize(0.025); lat02_2.Draw("SAME");
  alpha02vals[2]=0.1130; alpha02errs[2]=0.0050; yy02[2]=0.45; yy02e[2]=0.;
  TLatex lat02_3(0.124,0.40-offset,"SESAM 1999"); lat02_3.SetTextSize(0.025); lat02_3.Draw("SAME");
  alpha02vals[3]=0.1118; alpha02errs[3]=0.0017; yy02[3]=0.40; yy02e[3]=0.;
  TLatex lat02_4(0.124,0.35-offset,"Wingate #it{et al.} 1995"); lat02_4.SetTextSize(0.025); lat02_4.Draw("SAME");
  alpha02vals[4]=0.1110; alpha02errs[4]=0.0060; yy02[4]=0.35; yy02e[4]=0.;
  TLatex lat02_5(0.124,0.30-offset,"Davies #it{et al.} 1994"); lat02_5.SetTextSize(0.025); lat02_5.Draw("SAME");
  alpha02vals[5]=0.1150; alpha02errs[5]=0.0017; yy02[5]=0.30; yy02e[5]=0.;

  TLatex latnf02(0.103,0.36,"n_{f} = 0,2"); latnf02.SetTextSize(0.03); latnf02.SetTextAngle(90); latnf02.Draw("SAME");

  TGraphErrors* graph02=new TGraphErrors(NN02,alpha02vals,yy02,alpha02errs,yy02e);
  graph02->SetMarkerColor(1);
  graph02->SetMarkerStyle(33);
  graph02->SetMarkerSize(1);
  graph02->Draw("P");

  // alpha 3 -----------------------------

  TLatex lat3_0(0.124,0.85-offset,"Bazavov #it{et al.} 2014"); lat3_0.SetTextSize(0.025); lat3_0.Draw("SAME");
  alpha3vals[0]=0.1166; alpha3errs[0]=0.0012; yy3[0]=0.85; yy3e[0]=0.;
  TLatex lat3_1(0.124,0.80-offset,"HPQCD 2013"); lat3_1.SetTextSize(0.025); lat3_1.Draw("SAME");
  alpha3vals[1]=0.1183; alpha3errs[1]=0.0006; yy3[1]=0.80; yy3e[1]=0.;
  TLatex lat3_2(0.124,0.75-offset,"PACS-CS 2009"); lat3_2.SetTextSize(0.025); lat3_2.Draw("SAME");
  alpha3vals[2]=0.1205; alpha3errs[2]=0.0019; yy3[2]=0.75; yy3e[2]=0.;
  TLatex lat3_3(0.124,0.70-offset,"Maltman #it{et al.} 2008"); lat3_3.SetTextSize(0.025); lat3_3.Draw("SAME");
  alpha3vals[3]=0.1192; alpha3errs[3]=0.0011; yy3[3]=0.70; yy3e[3]=0.;
  TLatex lat3_4(0.124,0.65-offset,"HPQCD/QCDSF 2005"); lat3_4.SetTextSize(0.025); lat3_4.Draw("SAME");
  alpha3vals[4]=0.1170; alpha3errs[4]=0.0012; yy3[4]=0.65; yy3e[4]=0.;

  TLatex latnf3(0.103,0.71,"n_{f} = 3"); latnf3.SetTextSize(0.03); latnf3.SetTextAngle(90); latnf3.Draw("SAME");

  TGraphErrors* graph3=new TGraphErrors(NN3,alpha3vals,yy3,alpha3errs,yy3e);
  graph3->SetMarkerColor(1);
  graph3->SetMarkerStyle(33);
  graph3->SetMarkerSize(1);
  graph3->Draw("P");

  Double_t line3x[2]={0.0,0.5};
  Double_t line3y[2]={0.6,0.6};
  TGraph* line3=new TGraph(2,line3x,line3y);
  line3->SetLineColor(1);
  line3->SetLineStyle(1);
  line3->Draw("L");  

  // alpha 4 -----------------------------

  TLatex lat4_0(0.124,0.95-offset,"ETM 2014"); lat4_0.SetTextSize(0.025); lat4_0.Draw("SAME");
  alpha4vals[0]=0.1196; alpha4errs[0]=0.0011; yy4[0]=0.95; yy4e[0]=0.;

  TLatex latnf4(0.103,0.91,"n_{f} = 4"); latnf4.SetTextSize(0.03); latnf4.SetTextAngle(90); latnf4.Draw("SAME");

  TGraphErrors* graph4=new TGraphErrors(NN4,alpha4vals,yy4,alpha4errs,yy4e);
  graph4->SetMarkerColor(1);
  graph4->SetMarkerStyle(33);
  graph4->SetMarkerSize(1);
  graph4->Draw("P");

  Double_t line4x[2]={0.0,0.5};
  Double_t line4y[2]={0.9,0.9};
  TGraph* line4=new TGraph(2,line4x,line4y);
  line4->SetLineColor(1);
  line4->SetLineStyle(1);
  line4->Draw("L");  

  // alpha world average -----------------------------

  Double_t linex[4]={0.1179,0.1179,0.1191,0.1191};
  Double_t liney[4]={0.,100.,100.,0.};
  TGraph* line=new TGraph(4,linex,liney);
  line->SetLineColor(1);
  line->SetLineStyle(3);
  line->Draw("L");

  // ---------------------------

  canvalpha->SaveAs("compare_alphas.pdf");

  return 0;

};


