/*
This script is used to extract the covariance matrix from the files that contain the statistical correlation coefficients:
/nfs/dust/cms/user/savitsky/ForJames/abs_parton/Unfolding_combined_TtBar_Mass_HypTTBarMass.root 
/nfs/dust/cms/user/savitsky/ForJames/abs_parton/Unfolding_combined_TopAntiQuark_Pt_HypAntiToppT.root

The results are printed to cov_mtt.txt & cov_pT.txt
*/
#include <TFile.h>
#include "TH1D.h"
#include "TH2D.h"
#include "TH2F.h"
#include <TMath.h>
#include <cmath>
#include <fstream>
#include <string>
#include <vector>
#include <TMatrixD.h>

void printCov(TMatrixD *cov,int N){
  for (int i =0;i<N;i++){
     for (int j=0;j<N;j++){
         printf("%.3f,  ",cov(i,j));
     }
     cout<<""<<endl;
  }
}
int extractCov(){
  TFile* f_mtt = TFile::Open("/nfs/dust/cms/user/savitsky/ForJames/abs_parton/Unfolding_combined_TtBar_Mass_HypTTBarMass.root");
  TH2D* h_mtt = (TH2D*) f_mtt->Get("SVD_combined_TtBar_Mass_HypTTBarMass_STATCORR");
  TMatrixD *covmat = new TMatrixD(7,7);
  double cov[7][7];
//extract correlation histogram
  cout<<"***************************mtt******************************"<<endl;
  for(int b = 2; b < h_mtt->GetNbinsX(); b++)
  {
     for (int b2 = b+1;b2<h_mtt->GetNbinsX();b2++){
        double cor = h_mtt->GetBinContent(b, b2)/100;
	cov[b-2][b2-2) = cor;
        cov[b2-2][b-2] = cor;
     }
  }
  for (int i=0;i<7;i++)
	cov[i][i] = 1;
  int N = 7;
    for (int i =0;i<N;i++){
     for (int j=0;j<N;j++){
         printf("%.3f,  ",cov[i][j]);
     }
     cout<<""<<endl;
  }
  //f_mtt->Close();
  TFile* dataF = TFile::Open("data/DiffXS_HypTTBarMass_source.root");
  TH1D* hist_dat = (TH1D*)dataF->Get("mc");
  TGraphAsymmErrors* gr_dat = (TGraphAsymmErrors*) dataF->Get("data");

  double err[7];
  for(int b = 0; b < hist_dat->GetNbinsX(); b++){
     err[b] = gr_dat->GetErrorY(b);
     cout<<err[b]<<endl;}
  
  ofstream myfile;
  myfile.open("cov_mtt.txt");
 for (int i=0;i<N;i++){
     for (int j=0;j<N;j++){
        covmat(i,j) = cov[i][j]*err[i]*err[j];
	myfile<<covmat(i,j)<<"     ";
     }
     cout<<""<<endl;
     myfile<<" "<<endl;
  }
//extract data uncertainties
  printCov(covmat,7);
  TMatrixD *covinv = new TMatrixD(7,7);
  covinv = covmat->Invert();

    for (int i =0;i<N;i++){
     for (int j=0;j<N;j++){
         printf("%.3f,  ",covinv(i,j));
	//myfile<<covmat(i,j)<<"     ";
     }
     cout<<""<<endl;
     //myfile<<" "<<endl;
  }
 f_mtt->Close();
 dataF->Close();
 myfile.close();

//repeat from pT
  cout<<"*********************************pT********************************"<<endl;
  TFile* f_mtt = TFile::Open("/nfs/dust/cms/user/savitsky/ForJames/abs_parton/Unfolding_combined_TopAntiQuark_Pt_HypAntiToppT.root");
  TH2D* h_mtt = (TH2D*) f_mtt->Get("SVD_combined_TopAntiQuark_Pt_HypAntiToppT_STATCORR");
  TMatrixD *covmat_pT = new TMatrixD(6,6);
  double cov_pT[6][6];
//extract correlation histogram
  for(int b = 2; b < h_mtt->GetNbinsX(); b++)
  {
     for (int b2 = b+1;b2<h_mtt->GetNbinsX();b2++){
        double cor = h_mtt->GetBinContent(b, b2)/100;
	cov_pT[b-2][b2-2) = cor;
        cov_pT[b2-2][b-2] = cor;
     }
  }
  for (int i=0;i<6;i++)
	cov_pT[i][i] = 1;
  int N = 6;
    for (int i =0;i<N;i++){
     for (int j=0;j<N;j++){
         printf("%.3f,  ",cov_pT[i][j]);
     }
     cout<<""<<endl;
  }
  //f_mtt->Close();
  TFile* dataF = TFile::Open("data/DiffXS_HypToppT_source.root");
  TH1D* hist_dat = (TH1D*)dataF->Get("mc");
  TGraphAsymmErrors* gr_dat = (TGraphAsymmErrors*) dataF->Get("data");

  double err_pT[6];
  for(int b = 0; b < hist_dat->GetNbinsX(); b++){
     err_pT[b] = gr_dat->GetErrorY(b);
     cout<<err_pT[b]<<endl;}
  
  ofstream myfile;
  myfile.open("cov_pT.txt");
 for (int i=0;i<N;i++){
     for (int j=0;j<N;j++){
        covmat_pT(i,j) = cov_pT[i][j]*err_pT[i]*err_pT[j];
	myfile<<covmat_pT(i,j)<<"     ";
     }
     cout<<""<<endl;
     myfile<<" "<<endl;
  }
//extract data uncertainties
  printCov(covmat_pT,6);

 f_mtt->Close();
 dataF->Close();
 myfile.close();


  return 0;
}
