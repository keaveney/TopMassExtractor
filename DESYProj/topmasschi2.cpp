
#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <cmath>

double calcChi2(vector <double> data, vector <double> pred);

int topmasschi2()
{
  // list of mt variations
  TString names[6] = {"mt168", "mt170", "mt172", "mt1733", "mt174", "mt176"};
  double mt[6] = {168.0, 170.0, 172.0, 173.3, 174.0, 176.0};
  
  // number of iterations times number of calls in each iteration
  const int niter1 = 5e5 * 10;
  const int niter2 = 2e5 * 10;
  
  // list of variables
  TString vars[2] = {"mtt", "ptt"};
  TString vars_dat[2] = {"TTBarMass", "ToppT"};
  TString titles[2] = {"M(ttbar)", "pT(t)"};
  
  //double maxreldiff = 0.0;
  //double averreldiff = 0.0;
  //int nreldiff = 0;
  double x[2][6];
  double y[2][6];
  for(int v = 0; v < 2; v++)
  {
    //extract data
    TFile* dataF = TFile::Open("data/DiffXS_Hyp"+vars_dat[v]+"_source.root");
    TH1D* hist_dat = (TH1D*)dataF->Get("mc");
    TGraphAsymmErrors* gr_dat = (TGraphAsymmErrors*) dataF->Get("data");
        
    double data_bin_x, data_bin_y, data_bin_error;
    double total = 0;
    int binwidth;
    //int size = hist_dat->GetNbinsX();
    double chi2[6];
    //double xsec_data[size];
    vector <double> xsec_data;
    cout<<"    "<<endl;
    cout<<"    "<<endl;
    cout<<"************Data Results************"<<endl;
    for (int b=0;b<hist_dat->GetNbinsX(); b++)
    {
      binwidth = hist_dat->GetBinLowEdge(b + 2) - hist_dat->GetBinLowEdge( b + 1);
      //binwidth = 1;
      gr_dat->GetPoint(b, data_bin_x,data_bin_y);
      xsec_data.push_back(binwidth*data_bin_y);
      total += xsec_data.back();
      cout<<"bin "<<data_bin_x<<"  x-sec = "<<xsec_data.back()<<" pb/GeV"<<endl;
    }
    cout<<"total x-sec data: "<<total<<" pb/GeV"<<endl;

    //extract prediction from histogram
    for(int i =0; i < 6; i++)
    {
      // read M(ttbar) calculated with 500k iterations
      TFile* f500k = TFile::Open("MS_FromSasha_Iter1/MS-" + names[i] + "-500k/grid-TTbar_MS_" + vars[v] + ".root");
      TDirectoryFile* dir500k = (TDirectoryFile*) f500k->Get("grid");
      TH1D* h500k = (TH1D*) dir500k->Get("reference");
      printf("%s\n", std::string(50, '*').c_str());
      printf("*** NLO prediction mt = %.1f GeV***\n", mt[i]);
      double total = 0.0;
      //double xsec_pred[size];
      vector <double> xsec_pred;
      for(int b = 0; b < h500k->GetNbinsX(); b++)
      {
        double binwidth = h500k->GetBinLowEdge(b + 2) - h500k->GetBinLowEdge(b + 1);
        // x-section should be multiplied by the number of iterations
        double xsec = h500k->GetBinContent(b + 1) * niter1;
	xsec_pred.push_back(h500k->GetBinContent(b + 1) * niter1 * binwidth);
        // x-section already divided by bin width: needs to be multiplied to sum up to the total x-section
        total += xsec_pred.back();
        printf("%.1f < %s < %.1f GeV  xsec = %.3f pb/GeV\n", h500k->GetBinLowEdge(b + 1), titles[v].Data(), h500k->GetBinLowEdge(b + 2), xsec_pred.back()/1000);
      }
      // add over/underflow bins to the total x-section: these were not divided by bin width
      total += h500k->GetBinContent(0) * niter1;
      total += h500k->GetBinContent(h500k->GetNbinsX() + 1) * niter1;
      printf("*** Integrated x-section (%.1f < %s < %.1f GeV) = %.3e fb***\n", h500k->GetBinLowEdge(1), titles[v].Data(), h500k->GetBinLowEdge(h500k->GetNbinsX() + 1), h500k->Integral(1, h500k->GetNbinsX(), "width") * niter1);
      printf("*** Total x-section = %.3f pb***\n", total/1000);

      //Calculate chi2 btw data & this particular prediction
      //chi2.push_back(calcChi2(xsec_data,xsec_pred));
      chi2[i] = calcChi2(xsec_data,xsec_pred);
      cout<<"**********************chi2***************************"<<endl;
      cout<<"chi2 btw data & mt = "<<mt[i]<<" is: "<<chi2[i]<<endl;


    /*  // read M(ttbar) calculated with 200k iterations to evaluate uncertainty as difference (500k, 200k)
      TFile* f200k = TFile::Open("MS_FromSasha_Iter1/MS-" + names[i] + "-200k/grid-TTbar_MS_" + vars[v] + ".root");
      TDirectoryFile* dir200k = (TDirectoryFile*) f200k->Get("grid");
      TH1D* h200k = (TH1D*) dir200k->Get("reference");
      for(int b = 0; b < h200k->GetNbinsX(); b++)
      {
        double xsec1 = h500k->GetBinContent(b + 1) * niter1;
        double xsec2 = h200k->GetBinContent(b + 1) * niter2;
        double reldiff = TMath::Abs((xsec2 - xsec1) / xsec1);
        //printf("reldiff = %e bin = %d var = %d\n", reldiff, b, v);
        averreldiff = reldiff * reldiff + averreldiff * averreldiff;
        nreldiff++;
        if(reldiff > maxreldiff)
          maxreldiff = reldiff;
      }*/
      f500k->Close();
    }
      //cout<<"min element in chi2 is: "<<*min_element(chi2.begin(),chi2.end())<<endl;
      //now subtract the minimum chi2
      double min_chi2 = chi2[0];
      for (int i =1;i<6;i++){
        if (chi2[i]<min_chi2) min_chi2 = chi2[i];}
      for (int i = 0;i<6;i++){
        chi2[i] = chi2[i]-min_chi2;
	x[v][i] = mt[i];
	y[v][i] = chi2[i];
	}
      //cout<<"x= "<<x[1][5]<<endl;
      //cout<<"y= "<<y[0][5]<<endl;
      dataF->Close();

  }
  /*printf("*** estimation of intrinsic uncertainties of calculations:  ");
  printf("average relative difference = %e  ", TMath::Sqrt(averreldiff / nreldiff));
  printf("maximum relative difference = %e ***\n", maxreldiff);*/
  
 double xx[6],y1[6],y2[6];
 for (int i=0;i<6;i++){
  xx[i] = x[0][i];
  y1[i] = y[0][i];
  y2[i] = y[1][i];
 }
 //cout<<"y2 = "<<y[2][0]<<endl;
 
 TGraph *gr1 = new TGraph(6,xx,y1);
 TGraph *gr2 = new TGraph(6,xx,y2);
 TCanvas *c = new TCanvas("c","chi2",200,10,700,500);
 gr1->SetMarkerColor(4);
 //gr1->SetMarkerStyle(21);
 gr1->Draw("AC*");
 gr1->SetTitle("chi^2 for ttbar mass");
 gr1->GetXaxis()->SetTitle("t mass (GeV)");
 gr1->GetYaxis()->SetTitle("delta chi^2 [pb]");
 gr2->SetMarkerColor(2);
 gr2->SetMarkerStyle(21);
 gr2->SetLineColor(2);
 gr2->Draw("CP");

 TLegend* legend = new TLegend (0.1, 0.8, 0.25, 0.9);
    
 legend->AddEntry(gr1,"mtt","l");
 legend->AddEntry(gr2,"pT","l");
 legend->Draw("same");
 c->SaveAs("chi2.root");
 return 0;
}

double calcChi2(vector <double> data, vector <double> pred)
{
  int N = data.size();
  double chi2;
  for (int i=0;i<N;i++){
     chi2+= pow((data.at(i)-pred.at(i)/1000),2)/data.at(i);}
  return chi2;
}

	
