// compile:
// g++ -o readhisto `root-config --cflags --libs` readhisto.cpp
// run:
// ./readhisto

#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>

int main()
{
  // list of mt variations
  TString names[6] = {"mt168", "mt170", "mt172", "mt1733", "mt174", "mt176"};
  double mt[6] = {168.0, 170.0, 172.0, 173.3, 174.0, 176.0};
  
  // list of variables
  TString vars[2] = {"mtt", "ptt"};
  TString titles[2] = {"M(ttbar)", "pT(t)"};
  
  double maxreldiff = 0.0;
  double averreldiff = 0.0;
  int nreldiff = 0;
  for(int v = 0; v < 2; v++)
  {
    for(int i =0; i < 6; i++)
    {
      // read M(ttbar) calculated with 500k iterations
      TFile* f500k = TFile::Open("MS-" + names[i] + "-500k/grid-TTbar_MS_" + vars[v] + ".root");
      TDirectoryFile* dir500k = (TDirectoryFile*) f500k->Get("grid");
      TH1D* h500k = (TH1D*) dir500k->Get("reference");
      printf("%s\n", std::string(50, '*').c_str());
      printf("*** NLO prediction mt = %.1f ***\n", mt[i]);
      double total = 0.0;
      for(int b = 0; b < h500k->GetNbinsX(); b++)
      {
        double binwidth = h500k->GetBinLowEdge(b + 2) - h500k->GetBinLowEdge(b + 1);
        // x-section should be multiplied by the number of iterations
        double xsec = h500k->GetBinContent(b + 1) * 5e5;
        // x-section already divided by bin width: needs to be multiplied to sum up to the total x-section
        total += h500k->GetBinContent(b + 1) * 5e5 * binwidth;
        printf("%.1f < %s < %.1f  xsec = %.3e\n", h500k->GetBinLowEdge(b + 1), titles[v].Data(), h500k->GetBinLowEdge(b + 2), xsec);
      }
      // add over/underflow bins to the total x-section: these were not divided by bin width
      total += h500k->GetBinContent(0) * 5e5;
      total += h500k->GetBinContent(h500k->GetNbinsX() + 1) * 5e5;
      printf("*** Integrated x-section (%.1f < %s < %.1f) = %.3e ***\n", h500k->GetBinLowEdge(1), titles[v].Data(), h500k->GetBinLowEdge(h500k->GetNbinsX() + 1), h500k->Integral(1, h500k->GetNbinsX(), "width") * 5e5);
      printf("*** Total x-section = %.3e ***\n", total);

      // read M(ttbar) calculated with 200k iterations to evaluate uncertainty as difference (500k, 200k)
      TFile* f200k = TFile::Open("MS-" + names[i] + "-200k/grid-TTbar_MS_" + vars[v] + ".root");
      TDirectoryFile* dir200k = (TDirectoryFile*) f200k->Get("grid");
      TH1D* h200k = (TH1D*) dir200k->Get("reference");
      for(int b = 0; b < h200k->GetNbinsX(); b++)
      {
        double xsec1 = h500k->GetBinContent(b + 1) * 5e5;
        double xsec2 = h200k->GetBinContent(b + 1) * 2e5;
        double reldiff = TMath::Abs((xsec2 - xsec1) / xsec1);
        //printf("reldiff = %e bin = %d var = %d\n", reldiff, b, v);
        averreldiff = reldiff * reldiff + averreldiff * averreldiff;
        nreldiff++;
        if(reldiff > maxreldiff)
          maxreldiff = reldiff;
      }
      f500k->Close();
      f200k->Close();
    }
  }
  printf("*** estimation of intrinsic uncertainties of calculations:  ");
  printf("average relative difference = %e  ", TMath::Sqrt(averreldiff / nreldiff));
  printf("maximum relative difference = %e ***\n", maxreldiff);
    
  return 0;
}
