/// This is where you write some code to
///   1. Access histogram of *data measurement* (ROOT files in data/ directory) of ttbar differential cross section as a function of a) Pt top, b) Mtt
///   2. Access histograms of predctions for these differntial cross sections for different values of Mt 
///   3. Loop over the predictions and calculate the chi^2 between data and each prediciton
///   4. Make a final plot of chi^2 vs. Mt.
///
/// You should study the code in readhisto.cpp to see how to access histograms from ROOT files.

#include <TFile.h>
#include <TH1F.h>
#include <TMath.h>
#include <iostream>
#include <fstream>

using namespace std;

int topmassextractor()
{	 
  	// list of variables
  	TString vars_dat[2] = {"TTBarMass", "ToppT"};
	TString names[6] = {"mt168", "mt170", "mt172", "mt1733", "mt174", "mt176"};
	TString vars[2] = {"mtt", "ptt"};
	double mt[6] = {168.0, 170.0, 172.0, 173.3, 174.0, 176.0};
	ofstream myfile;
	myfile.open("hist_extractor.txt");
	for (int i = 0; i<2;i++){
		//extract data
		TFile* dataF = TFile::Open("data/DiffXS_Hyp"+vars_dat[i]+"_source.root");
		TH1D* hist_dat = (TH1D*) dataF->Get("mc");
		double total = 0;
		//		cout<<"*****Data Results*****\n";
		myfile<<"*****Data Results*****\n";
		for(int b = 0; b < hist_dat->GetNbinsX(); b++)
      			{
       			 double binwidth = hist_dat->GetBinLowEdge(b + 2) - hist_dat->GetBinLowEdge(b + 1);
        		// x-section should be multiplied by the number of iterations
       			 double xsec = hist_dat->GetBinContent(b + 1);
        		// x-section already divided by bin width: needs to be multiplied to sum up to the total x-section
        		total += hist_dat->GetBinContent(b + 1)* binwidth;
			cout<<setiosflags(ios::fixed)<<hist_dat->GetBinLowEdge(b+1)<<"<"<<vars_dat[i]<<"<"<<hist_dat->GetBinLowEdge(b+2)<<" xsec = "<<setiosflags(ios::scientific)<<xsec<<endl;
			myfile<<setiosflags(ios::fixed)<<hist_dat->GetBinLowEdge(b+1)<<"<"<<vars_dat[i]<<"<"<<hist_dat->GetBinLowEdge(b+2)<<" xsec = "<<setiosflags(ios::scientific)<<xsec<<endl;
			}
			total+=hist_dat->GetBinContent(0);
			total+=hist_dat->GetBinContent(hist_dat->GetNbinsX()+1);
		//extract prediction
		for (int n = 0;n<6;n++)
		{TFile* f500k = TFile::Open("MS_FromSasha_Iter1/MS-" + names[n] + "-500k/grid-TTbar_MS_" + vars[i] + ".root");
      		TDirectoryFile* dir500k = (TDirectoryFile*) f500k->Get("grid");
      		TH1D* h500k = (TH1D*) dir500k->Get("reference");
     		printf("%s\n", std::string(50, '*').c_str());
      		printf("*** NLO prediction mt = %.1f ***\n", mt[n]);
		myfile<<string(50, '*').c_str()<<endl<<"*** NLO prediction mt = "<<mt[n]<<endl;
      		double total_pred = 0.0;
			for(int b = 0; b < h500k->GetNbinsX(); b++){
			double binwidth = h500k->GetBinLowEdge(b + 2) - h500k->GetBinLowEdge(b + 1);
			// x-section should be multiplied by the number of iterations
        		double xsec_pred = h500k->GetBinContent(b + 1) * 5e5;
        		// x-section already divided by bin width: needs to be multiplied to sum up to the total x-section
        		total_pred += h500k->GetBinContent(b + 1) * 5e5 * binwidth;
			printf("%.1f < %s < %.1f  xsec = %.3e\n", h500k->GetBinLowEdge(b + 1), vars_dat[i].Data(), h500k->GetBinLowEdge(b + 2), xsec_pred);
			myfile<<h500k->GetBinLowEdge(b + 1)<<"<"<<vars_dat[i].Data()<<"<"<<h500k->GetBinLowEdge(b + 2)<<", x_sec_pred = "<<xsec_pred<<endl;
			}
      			// add over/underflow bins to the total x-section: these were not divided by bin width
      		total_pred += h500k->GetBinContent(0) * 5e5;
      		total_pred += h500k->GetBinContent(h500k->GetNbinsX() + 1) * 5e5;

		printf("*** Total x-section_pred = %.3e ***\n", total_pred);
	      	printf("*** Integrated x-section_pred (%.1f < %s < %.1f) = %.3e ***\n", h500k->GetBinLowEdge(1), vars_dat[i].Data(), h500k->GetBinLowEdge(h500k->GetNbinsX() + 1), h500k->Integral(1, h500k->GetNbinsX(), "width") * 5e5);
		myfile<<"*** Total x-section_pred = "<<total_pred<<endl<<"Integrated x-sec_pred = "<<h500k->Integral(1, h500k->GetNbinsX(), "width") * 5e5<<endl;
		f500k->Close();
		}
		cout<<"total x sec = "<<setiosflags(ios::scientific)<<total<<endl;
		myfile<<"****Total x-sec from data = " <<setiosflags(ios::scientific)<<total<<endl<<"***************************************************"<<endl;
		dataF->Close();
	}
	myfile.close();
	
	return 0;
}
 
