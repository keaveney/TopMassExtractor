/*
This script plots the chi2 results calculated by xfitter for NLO predictions and kfactors applied to NLO predictions (estimated NNLO predictions)
All results are saved to the folder "prefinal".

Note the scale errors are manually extracted by repeating the fits, hence the chi2_up and chi2_down's. All results from xfitter are saved in a folder named "xfitter". 
*/
#include "TH1D.h"
#include "TH1F.h"
#include <TMath.h>

int xfitter_diff(){
  //before k mtt
  double mt[]={168,170,172,173.3,174,176};
  double chi2_b4k[] = {3.58, 1.36,11.38,28.73,46.25,133.78};
  double chi2_b4k_err[] = {0.13,0.46,1.1,1.53,1.78,2.46};
  for (int i =0;i<6;i++){
    chi2_b4k[i] = chi2_b4k[i];
  }
//pT b4k
  double chi2_b4k_pT[] = {4.1,3.3,7.8,14,18,36};
//  double chi2_b4k_err[] = {0.13,0.46,1.1,1.53,1.78,2.46};
  for (int i =0;i<6;i++){
    chi2_b4k_pT[i] = chi2_b4k_pT[i];
  }
  TCanvas *c1 = new TCanvas("c1","chi2_mtt",300,100,700,500);
  gr = new TGraph(6);
  for (int i =0;i<6;i++){
   gr->SetPoint(i, mt[i], chi2_b4k[i]);}
  gr->SetTitle("xfit chi2 before k mtt");
  gr->GetXaxis()->SetTitle("mt [GeV]");
  gr->GetYaxis()->SetTitle("chi2");
  gr->SetMinimum(-1.5);
  gr->SetMarkerStyle(3);
  gr->Draw("AP");
  f_b4k = new TF1("f_b4k","pol2",mt[0],mt[5]);
  f_b4k->SetLineColor(2);
  gr->Fit("f_b4k");
  f_b4k->Draw("same");
  TLegend* legend = new TLegend (0.1, 0.8, 0.25, 0.9);
  legend->AddEntry(gr,"mtt_b4k","p");
  legend->AddEntry(f_b4k,"mtt_b4k_fit","l");
  legend->Draw("same");
  c1->SaveAs("prefinal/xfitter_chi2_b4k_mtt.root");
 
  TCanvas *c2 = new TCanvas("c2","chi2_pT",300,100,700,500);
  gr_pT = new TGraph(6);
 for (int i =0;i<6;i++){
   gr_pT->SetPoint(i,mt[i],chi2_b4k_pT[i]);}
  gr_pT->SetTitle("xfit chi2 before k pT");
  gr_pT->GetXaxis()->SetTitle("mt [GeV]");
  gr_pT->GetYaxis()->SetTitle("chi2");
  gr_pT->SetMinimum(-1.5);
  gr_pT->SetMarkerStyle(34);
  gr_pT->Draw("AP");
  f_b4k_pT = new TF1("f_b4k_pT","pol2",mt[0],mt[5]);
  f_b4k_pT->SetLineColor(4);
  gr_pT->Fit("f_b4k_pT");
  f_b4k_pT->Draw("same");
  TLegend* legend2 = new TLegend (0.1, 0.8, 0.25, 0.9);
  legend2->AddEntry(gr_pT,"pT_b4k","p");
  legend2->AddEntry(f_b4k_pT,"pT_b4k_fit","l");
  legend2->Draw("same");
  c2->SaveAs("prefinal/xfitter_chi2_b4k_pT.root");

//extract best fit mtt
   double xerr_xfit1,xerr_xfit2,minX_xfit,minY;
   minX_xfit = f_b4k->GetMinimumX();
   minY = f_b4k->GetMinimum();
   xerr_xfit1 = f_b4k->GetX(minY+1,0,minX_xfit);
   xerr_xfit2 = f_b4k->GetX(minY+1,minX_xfit,200);
   cout<<"*********************best fit mtt top mass = "<<minX_xfit<<" -/+ "<<minX_xfit-xerr_xfit1<<" & "<<xerr_xfit2-minX_xfit<<"*********************"<<endl;
//extract best fit mtt
   minX_xfit = f_b4k_pT->GetMinimumX();
   minY = f_b4k_pT->GetMinimum();
   xerr_xfit1 = f_b4k_pT->GetX(minY+1,0,minX_xfit);
   xerr_xfit2 = f_b4k_pT->GetX(minY+1,minX_xfit,200);
   cout<<"*********************best fit pT top mass = "<<minX_xfit<<" -/+ "<<minX_xfit-xerr_xfit1<<" & "<<xerr_xfit2-minX_xfit<<"*********************"<<endl;
   
  //c1->SaveAs("prefinal/xfit_chi2_b4k_mtt.root");
  //after k mtt
  double chi2_afterk[] = {10.42, 3.88,4.96,12.74,22.14,74.33};
  double chi2_up_afterk[] = {12.57,5.09,4.01,9.44,16.74,59.81};
  double chi2_down_afterk[] = {6.71,2.68,9.08,22.52,36.85,110.21};
  //double chi2_b4k_err[] = {0.13,0.46,1.1,1.53,1.78,2.46};
  for (int i =0;i<6;i++){
    chi2_afterk[i] = chi2_afterk[i];
    chi2_up_afterk[i] = chi2_up_afterk[i];
    chi2_down_afterk[i] = chi2_down_afterk[i];
  }
//pT after k
 double chi2_afterk_pT[] = {9.49,4.38,3.07,4.29,5.79,13.77};
 double chi2_up_afterk_pT[] = {10.16,4.71,2.92,3.72,5.00,12.13};
 double chi2_down_afterk_pT[] = {7.51,3.36,3.37,5.73,7.84,18.16};
  //double chi2_b4k_err[] = {0.13,0.46,1.1,1.53,1.78,2.46};
  for (int i =0;i<6;i++){
    chi2_afterk_pT[i] = chi2_afterk_pT[i];
    chi2_up_afterk_pT[i] = chi2_up_afterk_pT[i];
    chi2_down_afterk_pT[i] = chi2_down_afterk_pT[i];
  }
  TCanvas *c3 = new TCanvas("c3","chi2_afterk_mtt",300,100,700,500);
  gr2 = new TGraph(6);
  gr2_up = new TGraph(6);
  gr2_down = new TGraph(6);

  for (int i =0;i<6;i++){
   gr2->SetPoint(i, mt[i], chi2_afterk[i]);
   gr2_up->SetPoint(i,mt[i],chi2_up_afterk[i]);
   gr2_down->SetPoint(i,mt[i],chi2_down_afterk[i]);
  }
  gr2->SetTitle("xfit chi2 after k mtt");
  gr2->GetXaxis()->SetTitle("mt [GeV]");
  gr2->GetYaxis()->SetTitle("chi2");
  gr2->SetMinimum(-5);
  gr2->SetMarkerColor(2);
  gr2->SetMarkerStyle(2);
  gr2_up->SetMarkerColor(3);
  gr2_down->SetMarkerColor(4);
  gr2_up->SetMarkerStyle(3);
  gr2_down->SetMarkerStyle(5);
  gr2->Draw("AP");
  gr2_up->Draw("P");
  gr2_down->Draw("P");
  f_afterk = new TF1("f_afterk","pol2",mt[0],mt[5]);
  f_afterk_up= new TF1("f_afterk_up","pol2",mt[0],mt[5]);
  f_afterk_down= new TF1("f_afterk_down","pol2",mt[0],mt[5]);
  f_afterk->SetLineColor(2);
  f_afterk_up->SetLineColor(3);
  f_afterk_down->SetLineColor(4);
  gr2->Fit("f_afterk");
  gr2_up->Fit("f_afterk_up");
  gr2_down->Fit("f_afterk_down");
  f_afterk->Draw("same");
  f_afterk_up->Draw("same");
  f_afterk_down->Draw("same");

   minX_xfit = f_afterk->GetMinimumX();
   minY = f_afterk->GetMinimum();
   xerr_xfit1 = f_afterk->GetX(minY+1,0,minX_xfit);
   xerr_xfit2 = f_afterk->GetX(minY+1,minX_xfit,200);
   double minXup,minXdown,xerr_xfit_up,xerr_xfit_down;
   minXup = f_afterk_up->GetMinimumX();
   minY = f_afterk_up->GetMinimum();
   xerr_xfit_up = f_afterk_up->GetX(minY+1,minXup,200);
   xerr_xfit_up = xerr_xfit_up-minXup;

   minXdown= f_afterk_down->GetMinimumX();
   minY = f_afterk_down->GetMinimum();
   xerr_xfit_down = f_afterk_down->GetX(minY+1,0,minXdown);
   xerr_xfit_down = minXdown-xerr_xfit_down;
 double avg_X=(minX_xfit+minXup+minXdown)/3;
 double X_err_down = avg_X-(minXdown-xerr_xfit_down);
 double X_err_up = (minXup+xerr_xfit_up)-avg_X;
 cout<<"best fit top mass from mtt xfitter = "<<avg_X<<" +/- "<<X_err_up<<"/"<<X_err_down<<endl;
TPaveText *pt1 = new TPaveText(0.7,0.86,0.9,0.9,"NDC");
 pt1->AddText("top mass from mtt xfit = 170.503+/-0.960/0.876 [GeV]");
 pt1->Draw("same");
 TLegend* legend3 = new TLegend (0.1, 0.8, 0.25, 0.9);
 legend3->AddEntry(gr2,"mtt_afterk","p");
 legend3->AddEntry(gr2_up,"mtt_afterk_up","p");
 legend3->AddEntry(gr2_down,"mtt_afterk_down","p");
 legend3->AddEntry(f_afterk,"mtt_afterk_fit","l");
 legend3->AddEntry(f_afterk_up,"mtt_afterk_up_fit","l");
 legend3->AddEntry(f_afterk_down,"mtt_afterk_down_fit","l");
 legend3->Draw("same");

 c3->SaveAs("prefinal/xfitter_chi2_afterk_mtt.root");

  TCanvas *c4 = new TCanvas("c4","chi2_afterk_pT",300,100,700,500);
  gr2_pT = new TGraph(6);
  gr2_pT_up= new TGraph(6);
  gr2_pT_down = new TGraph(6);
  for (int i =0;i<6;i++){
  gr2_pT->SetPoint(i,mt[i],chi2_afterk_pT[i]);
  gr2_pT_up->SetPoint(i,mt[i],chi2_up_afterk_pT[i]);
  gr2_pT_down->SetPoint(i,mt[i],chi2_down_afterk_pT[i]);
  }
  gr2_pT->SetTitle("xfit chi2 after k pT");
  gr2_pT->GetXaxis()->SetTitle("mt [GeV]");
  gr2_pT->GetYaxis()->SetTitle("chi2");
  gr2_pT->SetMarkerColor(2);
  gr2_pT_up->SetMarkerColor(3);
  gr2_pT_down->SetMarkerColor(4);
  gr2_pT->SetMarkerStyle(2);
  gr2_pT_up->SetMarkerStyle(3);
  gr2_pT_down->SetMarkerStyle(5);
  f_afterk_pT = new TF1("f_afterk_pT","pol2",mt[0],mt[5]);
  f_afterk_pT->SetLineColor(2);
  f_afterk_pT_up = new TF1("f_afterk_pT_up","pol2",mt[0],mt[5]);
  f_afterk_pT_up->SetLineColor(3);
  f_afterk_pT_down = new TF1("f_afterk_pT_down","pol2",mt[0],mt[5]);
  f_afterk_pT_down->SetLineColor(4);
  gr2_pT->Fit("f_afterk_pT");
  gr2_pT_up->Fit("f_afterk_pT_up");
  gr2_pT_down->Fit("f_afterk_pT_down");
  gr2_pT->Draw("AP");
  gr2_pT_up->Draw("P");
  gr2_pT_down->Draw("P");
  f_afterk_pT->Draw("same");
  f_afterk_pT_up->Draw("same");
  f_afterk_pT_down->Draw("same");

   minX_xfit = f_afterk_pT->GetMinimumX();
   minY = f_afterk_pT->GetMinimum();
   xerr_xfit1 = f_afterk_pT->GetX(minY+1,0,minX_xfit);
   xerr_xfit2 = f_afterk_pT->GetX(minY+1,minX_xfit,200);

   minXup = f_afterk_pT_up->GetMinimumX();
   minY = f_afterk_pT_up->GetMinimum();
   xerr_xfit_up = f_afterk_pT_up->GetX(minY+1,minXup,200);
   xerr_xfit_up = xerr_xfit_up-minXup;

   minXdown= f_afterk_pT_down->GetMinimumX();
   minY = f_afterk_pT_down->GetMinimum();
   xerr_xfit_down = f_afterk_pT_down->GetX(minY+1,0,minXdown);
   xerr_xfit_down = minXdown-xerr_xfit_down;
 double avg_X=(minX_xfit+minXup+minXdown)/3;
 double X_err_down = avg_X-(minXdown-xerr_xfit_down);
 double X_err_up = (minXup+xerr_xfit_up)-avg_X;
 cout<<"best fit top mass from pT xfitter = "<<avg_X<<" +/- "<<X_err_up<<"/"<<X_err_down<<endl;
TPaveText *pt2 = new TPaveText(0.7,0.86,0.9,0.9,"NDC");
 pt2->AddText("top mass from pT xfit = 171.428+/-1.758/1.786 [GeV]");
 pt2->Draw("same");
 TLegend* legend4 = new TLegend (0.1, 0.8, 0.25, 0.9);
 legend4->AddEntry(gr2_pT,"pT_afterk","p");
 legend4->AddEntry(f_afterk_pT,"pT_afterk_fit","l");
 legend4->AddEntry(gr2_pT_up,"pT_afterk_up","p");
 legend4->AddEntry(f_afterk_pT_up,"pT_afterk_up_fit","l");
 legend4->AddEntry(gr2_pT_down,"pT_afterk_down","p");
 legend4->AddEntry(f_afterk_pT_down,"pT_afterk_down_fit","l");
 legend4->Draw("same");
  c4->SaveAs("prefinal/xfitter_chi2_afterk_pT.root");

  return 0;
}
