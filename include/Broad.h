#include "Binning.h"
#include "TStopwatch.h"
#include <iostream>
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TROOT.h"
#include "TH1.h"
#include "TGraphErrors.h"
#include "TMultiGraph.h"
#include "TLegend.h"
#include "TPad.h"
#include "TCanvas.h"
#include "TStyle.h"

char DC[3] = {'D', 'C', '\0' }; char DFe[4] = {'D', 'F', 'e', '\0'}; char DPb[4] = {'D', 'P', 'b', '\0'};

char C[2]  = {'C', '\0'}; char Fe[3]  = {'F', 'e', '\0'}; char Pb[3]  = {'P', 'b', '\0'};

int empty_histo(TH1F* h) {

  int empty = 0;
  for(int i = 1 ; i <= h->GetNbinsX() ; i++) {
    if(h->GetBinContent(i) == 0){ empty++; }
  }
  if(empty == h->GetNbinsX()) { return 1; }
  else { return 0; }
}

void PhiIntegration(TFile* inputFile, TFile* outputFile, char target[]) {

  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
    for(int Q2Counter = 0 ; Q2Counter < N_Q2; Q2Counter++) { // Loops in every Q2 bin
      for(int NuCounter = 0 ; NuCounter < N_Nu; NuCounter++) { // Loops in every Nu bin
        for(int ZhCounter = 0 ; ZhCounter < N_Zh; ZhCounter++) { // Loops in every Zh bin;
          // For every bin of Q2, Nu, Zh and number of pion, generate a histogram  for Pt2
	        TH1F* histPt2 = new TH1F(Form("corr_data_Pt2_%s_%i%i%i_%i_CLEAN_INTERPOLATED", target, Q2Counter, NuCounter, ZhCounter, nPion), "", N_Pt2, Pt2_MIN, Pt2_MAX);
          for(int Pt2Counter = 0 ; Pt2Counter < N_Pt2 ; Pt2Counter++) { // Loops in every Pt2 bin
            // Take the his for this bins
	          TH1F* histPhi = (TH1F*) inputFile->Get(Form("corr_data_%s_%i%i%i%i_%i", target, Q2Counter, NuCounter, ZhCounter, Pt2Counter, nPion));
            // If the histogram is null or empty skip this Pt2 bin
            if(histPhi == NULL){ continue; }
	          if(empty_histo(histPhi) == 1){ continue; }
            // Integrate on PhiPQ
	          // double errors_array[1] = {};
            double error;
	          double integral = histPhi->IntegralAndError(1,N_Phi,error);
	          // Save the value in the Pt2 histogram
	          histPt2->SetBinContent(Pt2Counter + 1, integral);
	          // histPt2->SetBinError(Pt2Counter + 1, errors_array[0]);
            histPt2->SetBinError(Pt2Counter + 1, error);
	        }// End Pt2 loop
          // If the histogram if not empty, save it
	        if(empty_histo(histPt2) == 0) {
	          outputFile->cd();
	          histPt2->Write();
            gROOT->cd();
	        }
          delete histPt2;
        }// End Zh loop
      }// End Nu loop
    }// End Q2 loop
  }// End number pion event loop

}

void NuZhIntegration(TFile* inputFile, TFile* outputFile, char target[]) {

  for(int nPion = 1; nPion <= N_PION ; nPion++) { // Loops in every number of pion
    // Generate a histogram to save Q2 for every number of pion in the final event
    TH1F* histQ2 = new TH1F(Form("corr_data_%s_%i_Q2", target, nPion), "", N_Q2, Q2_BINS);
    for(int Q2Counter = 0 ; Q2Counter < N_Q2 ; Q2Counter++) { // Loops in every Q2 bin
      TH1F* histPt2Integrated = new TH1F(Form("corr_data_Pt2_%s_%i", target, Q2Counter), "", N_Pt2, Pt2_MIN, Pt2_MAX);
      for(int NuCounter = 0 ; NuCounter < N_Nu ; NuCounter++) { // Loops in every Nu bin
        // Starts in the second bin because there is not broadening iz Zh<0.2
	      for(int ZhCounter = 1 ; ZhCounter < N_Zh ; ZhCounter++) { // Loops in every Zh bin
	        TH1F* histPt2 = (TH1F*) inputFile->Get(Form("corr_data_Pt2_%s_%i%i%i_%i_CLEAN_INTERPOLATED", target, Q2Counter, NuCounter, ZhCounter, nPion));
          // Sum the histograms for every bin of Zh and Nu
	        histPt2Integrated->Add(histPt2);
	        delete histPt2;
	      }// End Zh loop
	    }// End Nu loop
      // Take the mean and save in the Q2 histogram

      histQ2->SetBinContent(Q2Counter+1, histPt2Integrated->GetMean());
      histQ2->SetBinError(Q2Counter+1, histPt2Integrated->GetMeanError());
      delete histPt2Integrated;
    }// End Q2 loop

    // Open the direction of the output file and save the histograms
    outputFile->cd();
    histQ2->Write(Form("meanPt2_%s_%i_Q2_CLEAN_INTERPOLATED", target, nPion));
    gROOT->cd();
    delete histQ2;

  }// End number pion event

}

void Q2NuIntegration(TFile* inputFile, TFile* outputFile, char target[]) {

  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
    // Generate a histogram to save Zh for every number of pion in the final event
    TH1F* histZh = new TH1F(Form("corr_data_%s_%i_Zh", target, nPion), "", N_Zh, Zh_BINS);
    for(int ZhCounter = 0 ; ZhCounter < N_Zh ; ZhCounter++) { // Loops in every Zh bin
      // Generate a histogram for every bin of zh
      TH1F* histPt2Integrated = new TH1F(Form("corr_data_%s_Pt2_%i_%i", target, ZhCounter, nPion), "", N_Pt2, Pt2_MIN, Pt2_MAX);
      for(int NuCounter = 0 ; NuCounter < N_Nu ; NuCounter++) { // Loops in every Nu bin
	      for(int Q2Counter = 0 ; Q2Counter < N_Q2 ; Q2Counter++) { // Loops in every Q2 bin
          // Sum the histograms for every bin of Q2 and Nu
	        TH1F* histPt2 = (TH1F*) inputFile->Get(Form("corr_data_Pt2_%s_%i%i%i_%i_CLEAN_INTERPOLATED", target, Q2Counter, NuCounter, ZhCounter, nPion));
          histPt2Integrated->Add(histPt2);
          delete histPt2;
	      }// End Q2 loop
	    }// End Nu loop

      // Take the mean and save in the Zh histogram
      histZh->SetBinContent(ZhCounter+1, histPt2Integrated->GetMean());
      histZh->SetBinError(ZhCounter+1, histPt2Integrated->GetMeanError());

      // Open the direction of the output file and save the data
      delete histPt2Integrated;
    } // End Zh loop

    // Open the direction of the output file and save the data
    outputFile->cd();
    histZh->Write(Form("meanPt2_%s_%i_CLEAN_INTERPOLATED", target, nPion));
    gROOT->cd();
    delete histZh;
  } // End number pion event loop

}

void Q2ZhIntegration(TFile* inputFile, TFile* outputFile, char target[]) {

  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
    // Generate a histogram to save Nu for every number of pion in the final event
    TH1F* histNu = new TH1F(Form("corr_data_%s_Nu", target), "", N_Nu, Nu_BINS);
    for(int NuCounter = 0 ; NuCounter < N_Nu ; NuCounter++) { // Loops in every Nu bin
      TH1F* histPt2Integrated = new TH1F(Form("corr_data_Pt2_%s_%i", target, NuCounter), "", N_Pt2, Pt2_MIN, Pt2_MAX);
      // Starts in the second bin because there is not broadening in Zh<0.2
      for(int ZhCounter = 1 ; ZhCounter < N_Zh ; ZhCounter++) { // Loops in every Zh bin
	      for(int Q2Counter = 0 ; Q2Counter < N_Q2 ; Q2Counter++) { // Loops in every Q2 bin
	        TH1F* histPt2 =  (TH1F*) inputFile->Get(Form("corr_data_Pt2_%s_%i%i%i_%i_CLEAN_INTERPOLATED", target, Q2Counter, NuCounter, ZhCounter, nPion));
	        histPt2Integrated->Add(histPt2);
	        delete histPt2;
	      }// End Q2 loop
	    }// End Zh loop
      // Take the mean and save in the Nu histogram
      histNu->SetBinContent(NuCounter+1, histPt2Integrated->GetMean());
      histNu->SetBinError(NuCounter+1, histPt2Integrated->GetMeanError());
      delete histPt2Integrated;
    }// End Nu loop

    // Open the direction of the output file and save the data
    outputFile->cd();
    histNu->Write(Form("meanPt2_%s_%i_Nu_CLEAN_INTERPOLATED", target, nPion));
    gROOT->cd();
    delete histNu;
  } // End number pion event loop

}

void Q2NuZhIntegration(TFile* inputFile, TFile* outputFile, char target[]) {

  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
    // Generate two histograms for every number of pion in the final event
    // To save the mean of Pt2
    TH1F* hist = new TH1F(Form("corr_data_%s", target), "", 1, 0, 1);
    // To save the sum of the histograms
    TH1F* histPt2Integrated = new TH1F(Form("corr_data_%s_Pt2", target), "", N_Pt2, Pt2_MIN, Pt2_MAX);
    // Starts in the second bin because there is not broadening iz Zh<0.2
    for(int ZhCounter = 1 ; ZhCounter < N_Zh ; ZhCounter++) { // Loops in every Zh bin
      for(int NuCounter = 0 ; NuCounter < N_Nu ; NuCounter++) { // Loops in every Q2 bin
	      for(int Q2Counter = 0 ; Q2Counter < N_Q2 ; Q2Counter++) { // Loops in every Nu bin
          // Sum the histograms for every bin of Q2 and Nu
	        TH1F* histPt2 = (TH1F*) inputFile->Get(Form("corr_data_Pt2_%s_%i%i%i_%i_CLEAN_INTERPOLATED", target, Q2Counter, NuCounter, ZhCounter, nPion));
	        histPt2Integrated->Add(histPt2);
	        delete histPt2;
	      }// End Q2 loop
	    }// End Nu loop
    }// End Zh loop
    // Take the mean and save it
    hist->SetBinContent(1, histPt2Integrated->GetMean());
    hist->SetBinError(1, histPt2Integrated->GetMeanError());

    // Open the direction of the output file and save the data
    outputFile->cd();
    hist->Write(Form("meanPt2_%s_%i_CLEAN_INTERPOLATED",target, nPion));
    gROOT->cd();
    delete hist;
    delete histPt2Integrated;
  } // End number pion event loop


}

void CallPhiIntegration(TString inputDirectory, TString outputDirectory) {

  TFile* inputFile  = new TFile(inputDirectory  + "corr_data_Phi.root", "READ");
  TFile* outputFile = new TFile(outputDirectory + "corr_data_Pt2_processed.root", "RECREATE");
  gROOT->cd();

  PhiIntegration(inputFile, outputFile, DC);
  PhiIntegration(inputFile, outputFile, DFe);
  PhiIntegration(inputFile, outputFile, DPb);

  PhiIntegration(inputFile, outputFile, C);
  PhiIntegration(inputFile, outputFile, Fe);
  PhiIntegration(inputFile, outputFile, Pb);

  inputFile->Close();
  outputFile->Close();

}

void CallNuZhIntegration(TString inputDirectory, TString outputDirectory) {

  TFile* inputFile  = new TFile(inputDirectory  + "corr_data_Pt2_processed.root", "READ");
  TFile* outputFile = new TFile(outputDirectory + Form("meanPt2_Q2_%i.root",N_Q2), "RECREATE");
  gROOT->cd();

  NuZhIntegration(inputFile, outputFile, DC);
  NuZhIntegration(inputFile, outputFile, DFe);
  NuZhIntegration(inputFile, outputFile, DPb);

  NuZhIntegration(inputFile, outputFile, C);
  NuZhIntegration(inputFile, outputFile, Fe);
  NuZhIntegration(inputFile, outputFile, Pb);

  inputFile->Close();
  outputFile->Close();

}

void CallQ2NuIntegration(TString inputDirectory, TString outputDirectory) {

  TFile* inputFile  = new TFile(inputDirectory  + "corr_data_Pt2_processed.root", "READ");
  TFile* outputFile = new TFile(outputDirectory + Form("meanPt2_Zh_%i.root",N_Zh), "RECREATE");
  gROOT->cd();

  Q2NuIntegration(inputFile, outputFile, DC);
  Q2NuIntegration(inputFile, outputFile, DFe);
  Q2NuIntegration(inputFile, outputFile, DPb);

  Q2NuIntegration(inputFile, outputFile, C);
  Q2NuIntegration(inputFile, outputFile, Fe);
  Q2NuIntegration(inputFile, outputFile, Pb);

  inputFile->Close();
  outputFile->Close();

}

void CallQ2ZhIntegration(TString inputDirectory, TString outputDirectory) {

  TFile* inputFile  = new TFile(inputDirectory  + "corr_data_Pt2_processed.root", "READ");
  TFile* outputFile = new TFile(outputDirectory + Form("meanPt2_Nu_%i.root",N_Nu), "RECREATE");
  gROOT->cd();

  Q2ZhIntegration(inputFile, outputFile, DC);
  Q2ZhIntegration(inputFile, outputFile, DFe);
  Q2ZhIntegration(inputFile, outputFile, DPb);

  Q2ZhIntegration(inputFile, outputFile, C);
  Q2ZhIntegration(inputFile, outputFile, Fe);
  Q2ZhIntegration(inputFile, outputFile, Pb);

  inputFile->Close();
  outputFile->Close();

}

void CallQ2NuZhIntegration(TString inputDirectory, TString outputDirectory) {

  TFile* inputFile  = new TFile(inputDirectory  + "corr_data_Pt2_processed.root", "READ");
  TFile* outputFile = new TFile(outputDirectory + "meanPt2.root", "RECREATE");
  gROOT->cd();

  Q2NuZhIntegration(inputFile, outputFile, DC);
  Q2NuZhIntegration(inputFile, outputFile, DFe);
  Q2NuZhIntegration(inputFile, outputFile, DPb);

  Q2NuZhIntegration(inputFile, outputFile, C);
  Q2NuZhIntegration(inputFile, outputFile, Fe);
  Q2NuZhIntegration(inputFile, outputFile, Pb);

  inputFile->Close();
  outputFile->Close();

}

void SetErrorXNull(TGraphErrors* g) {

  double* errors_Y = g->GetEY();
  for(int point = 0 ; point < N_Nu ; point++) {
    g->SetPointError(point,0,errors_Y[point]);
  }

}

void SetXShift(TGraphErrors* g, double shift) {

  double* content_Y = g->GetY();
  double* content_X = g->GetX();

  for(int point = 0; point < N_Nu ; point++) {
    g->SetPoint(point, content_X[point] + shift , content_Y[point]);
  }

}

TLegend* GenTLegendTarget(TGraphErrors* g[N_STARGETS][N_PION]) {

  TLegend* legend = new TLegend(0.22,0.75,0.3,0.92,"","NDC");
  legend->AddEntry(g[0][0],"C","lpf");
  legend->AddEntry(g[1][0],"Fe","lpf");
  legend->AddEntry(g[2][0],"Pb","lpf");
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(62);
  legend->SetTextSize(0.04);
  return legend;

}

TLegend* GenTLegendNPion(TGraphErrors* g[N_STARGETS][N_PION]) {

  TLegend* legend = new TLegend(0.35,0.75,0.3,0.92,"","NDC");
  //legend->AddEntry((TObject*)0, "Final state", "");
  legend->AddEntry(g[2][0],"1 #pi +","lpf");
  legend->AddEntry(g[2][1],"2 #pi +","lpf");
  legend->AddEntry(g[2][2],"3 #pi +","lpf");
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(62);
  legend->SetTextSize(0.04);
  return legend;

}

TGraphErrors* TH1TOTGraph(TH1 *h1){

  if (!h1) std::cout << "TH1TOTGraph: histogram not found !" << std::endl;

  TGraphErrors* g1= new TGraphErrors();
  Double_t x, y, ex, ey;
  for (Int_t i=1 ; i<=h1->GetNbinsX(); i++) {
    y=h1->GetBinContent(i);
    ey=h1->GetBinError(i);
    x=h1->GetBinCenter(i);
    // ex=h1->GetBinWidth(i);
    ex=0.5*h1->GetBinWidth(i);
    g1->SetPoint(i-1,x,y);
    g1->SetPointError(i-1,ex,ey);
 }
 return g1;
}

void PtBroadeningNuIntegrated(TString inputDirectory, TString plotDirectory) {

  TFile* inputFile = new TFile(inputDirectory + Form("meanPt2_Nu_%i.root", N_Nu), "READ");

  TGraphErrors* g[N_STARGETS][N_PION];
  TCanvas* c = new TCanvas("c", "", 800, 600);
  c->Draw();
  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion

    TH1F* histSolid[N_STARGETS];   TH1F* histLiquid[N_STARGETS];  TH1F* histBroadening[N_STARGETS];

   //C
    histSolid[0]  = (TH1F*) inputFile->Get(Form("meanPt2_C_%i_Nu_CLEAN_INTERPOLATED", nPion));
    histLiquid[0] = (TH1F*) inputFile->Get(Form("meanPt2_DC_%i_Nu_CLEAN_INTERPOLATED", nPion));
   //Fe
    histSolid[1]  = (TH1F*) inputFile->Get(Form("meanPt2_Fe_%i_Nu_CLEAN_INTERPOLATED", nPion));
    histLiquid[1] = (TH1F*) inputFile->Get(Form("meanPt2_DFe_%i_Nu_CLEAN_INTERPOLATED", nPion));
   //Pb
    histSolid[2]  = (TH1F*) inputFile->Get(Form("meanPt2_Pb_%i_Nu_CLEAN_INTERPOLATED", nPion));
    histLiquid[2] = (TH1F*) inputFile->Get(Form("meanPt2_DPb_%i_Nu_CLEAN_INTERPOLATED", nPion));

    for(int i = 0 ; i < N_STARGETS; i++) {
      histBroadening[i] = new TH1F(Form("histBroadening_%i", i), "", N_Nu, Nu_BINS);
      histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
    }


    for(int i = 0 ; i < N_STARGETS ; i++){
      g[i][nPion-1] = (TGraphErrors*) TH1TOTGraph(histBroadening[i]);
      SetErrorXNull(g[i][nPion-1]);
    }

    g[0][nPion-1]->SetMarkerColor(kRed);
    g[0][nPion-1]->SetLineColor(kRed);
    g[1][nPion-1]->SetLineColor(kBlue);
    g[1][nPion-1]->SetMarkerColor(kBlue);
    g[2][nPion-1]->SetMarkerColor(kBlack);
    g[2][nPion-1]->SetLineColor(kBlack);

    SetXShift(g[0][nPion-1],-0.035);
    SetXShift(g[2][nPion-1],0.035);

    for(int i = 0; i < 3; i++) {
      delete histSolid[i];
      delete histLiquid[i];
      delete histBroadening[i];
    }
  }// End number pion event loop

  inputFile->Close();

  for(int k = 0; k < N_PION; k++) {
    g[k][0]->SetMarkerStyle(1);
    g[k][1]->SetMarkerStyle(8);
    g[k][1]->SetMarkerSize(.6);
    g[k][2]->SetMarkerStyle(27);
    g[k][2]->SetMarkerSize(1.4);
  }

  TLegend* legendTarget = GenTLegendTarget(g);
  TLegend* legendNPion= GenTLegendNPion(g);

  TPad* p = new TPad("p","p",0,0,1,1);
  p->SetLeftMargin(0.13);
  p->SetTopMargin(0.06);
  p->SetGridx(1);
  p->SetGridy(1);
  p->Draw();
  p->cd();

  TMultiGraph* mg = new TMultiGraph();
  for(int i = 0; i < N_STARGETS; i++) {
    for(int j = 0; j < N_PION; j++) {
        mg->Add(g[i][j]);
    }
  }
  mg->Draw("APE0");

  gStyle->SetTitleFont(62,"XY");
  gStyle->SetTitleSize(0.04,"XY");

  // mg->GetYaxis()->SetRangeUser(0.,0.059);
  mg->GetYaxis()->SetRangeUser(0.,0.11);
  mg->GetXaxis()->SetRangeUser(2.1,4.3);
  mg->GetXaxis()->SetTitle("#nu[GeV]");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  mg->Draw("APE0");
  legendTarget->Draw();
  legendNPion->Draw();

  c->Print(plotDirectory + "Pt_broad_Nu.pdf");
  delete c;

}

void PtBroadeningQ2Integrated(TString inputDirectory, TString plotDirectory) {

  TFile* inputFile = new TFile(inputDirectory + Form("meanPt2_Q2_%i.root", N_Q2), "READ");

  TCanvas* c = new TCanvas("c", "", 800, 600);
  c->Draw();

  TGraphErrors* g[N_STARGETS][N_PION];

  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
    TH1F* histSolid[N_STARGETS];   TH1F* histLiquid[N_STARGETS]; TH1F* histBroadening[N_STARGETS];

    //C
    histSolid[0] = (TH1F*) inputFile->Get(Form("meanPt2_C_%i_Q2_CLEAN_INTERPOLATED",nPion));
    histLiquid[0] = (TH1F*) inputFile->Get(Form("meanPt2_DC_%i_Q2_CLEAN_INTERPOLATED",nPion));
    //Fe
    histSolid[1] = (TH1F*) inputFile->Get(Form("meanPt2_Fe_%i_Q2_CLEAN_INTERPOLATED",nPion));
    histLiquid[1] = (TH1F*) inputFile->Get(Form("meanPt2_DFe_%i_Q2_CLEAN_INTERPOLATED",nPion));
    //Pb
    histSolid[2] = (TH1F*) inputFile->Get(Form("meanPt2_Pb_%i_Q2_CLEAN_INTERPOLATED",nPion));
    histLiquid[2] = (TH1F*) inputFile->Get(Form("meanPt2_DPb_%i_Q2_CLEAN_INTERPOLATED",nPion));


    for(int i = 0 ; i < N_STARGETS ; i++){
      histBroadening[i] = new TH1F(Form("histBroadening_%i",i), "", N_Q2, Q2_BINS);
      histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
    }


    for(int i = 0 ; i < N_STARGETS ; i++){
      g[i][nPion-1] = (TGraphErrors*) TH1TOTGraph(histBroadening[i]);
      SetErrorXNull(g[i][nPion-1]);
    }

    g[0][nPion-1]->SetMarkerColor(kRed);
    g[1][nPion-1]->SetMarkerColor(kBlue);
    g[0][nPion-1]->SetLineColor(kRed);
    g[1][nPion-1]->SetLineColor(kBlue);
    g[2][nPion-1]->SetMarkerColor(kBlack);
    g[2][nPion-1]->SetLineColor(kBlack);

    SetXShift(g[0][nPion-1],-0.035);
    SetXShift(g[2][nPion-1],0.035);

    for(int i = 0; i < N_STARGETS; i++){
      delete histSolid[i];
      delete histLiquid[i];
      delete histBroadening[i];
    }
  } // End number pion event loop

  inputFile->Close();

  TLegend* legendTarget = GenTLegendTarget(g);
  TLegend* legendNPion = GenTLegendNPion(g);

  TPad* p = new TPad("p","p",0,0,1,1);
  p->SetGridx(1);
  p->SetGridy(1);
  p->SetLeftMargin(0.13);
  p->SetTopMargin(0.06);
  p->Draw();
  p->cd();

  for(int k = 0; k < N_STARGETS; k++) {
    g[k][0]->SetMarkerStyle(1);
    g[k][1]->SetMarkerStyle(8);
    g[k][1]->SetMarkerSize(.6);
    g[k][2]->SetMarkerStyle(27);
    g[k][2]->SetMarkerSize(1.2);
  }

  TMultiGraph* mg = new TMultiGraph();
  for(int i = 0; i < N_STARGETS; i++){
    for(int j = 0; j < N_PION; j++){
        mg->Add(g[i][j]);
    }
  }

  mg->Draw("APE0");

  gStyle->SetTitleFont(62,"XY");
  gStyle->SetTitleSize(0.04,"XY");

  //mg->GetYaxis()->SetRangeUser(0.,0.059);
  mg->GetYaxis()->SetRangeUser(0.,0.12);
  mg->GetXaxis()->SetRangeUser(.9,4.1);
  mg->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  mg->Draw("APE0");
  legendTarget->Draw();
  legendNPion->Draw();

  c->Print(plotDirectory + "Pt_broad_Q2.pdf");
  delete c;

}

void PtBroadeningZhIntegrated(TString inputDirectory,  TString plotDirectory) {

  TFile* inputFile = new TFile(inputDirectory + Form("meanPt2_Zh_%i.root", N_Zh), "READ");

  TCanvas* c = new TCanvas("c", "", 800, 600);
  c->Draw();
  TGraphErrors* g[N_STARGETS][N_PION];

  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion

    TH1F* histSolid[N_STARGETS];  TH1F* histLiquid[N_STARGETS];  TH1F* histBroadening[N_STARGETS];

    //C
    histSolid[0]  = (TH1F*) inputFile->Get(Form("meanPt2_C_%i_CLEAN_INTERPOLATED",nPion));
    histLiquid[0] = (TH1F*) inputFile->Get(Form("meanPt2_DC_%i_CLEAN_INTERPOLATED",nPion));
    //Fe
    histSolid[1]  = (TH1F*) inputFile->Get(Form("meanPt2_Fe_%i_CLEAN_INTERPOLATED",nPion));
    histLiquid[1] = (TH1F*) inputFile->Get(Form("meanPt2_DFe_%i_CLEAN_INTERPOLATED",nPion));
    //Pb
    histSolid[2]  = (TH1F*) inputFile->Get(Form("meanPt2_Pb_%i_CLEAN_INTERPOLATED",nPion));
    histLiquid[2] = (TH1F*) inputFile->Get(Form("meanPt2_DPb_%i_CLEAN_INTERPOLATED",nPion));


    for(int i = 0 ; i < N_STARGETS ; i++) {
      histBroadening[i] = new TH1F(Form("histBroadening_%i",i),"", N_Zh, Zh_BINS);
      histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
    }


    for(int i = 0 ; i < N_STARGETS ; i++){
      g[i][nPion-1] = (TGraphErrors*) TH1TOTGraph(histBroadening[i]);
      SetErrorXNull(g[i][nPion-1]);
    }

    g[0][nPion-1]->SetMarkerColor(kRed);
    g[0][nPion-1]->SetLineColor(kRed);
    g[1][nPion-1]->SetMarkerColor(kBlue);
    g[1][nPion-1]->SetLineColor(kBlue);
    g[2][nPion-1]->SetMarkerColor(kBlack);
    g[2][nPion-1]->SetLineColor(kBlack);

    SetXShift(g[0][nPion-1], -0.015);
    SetXShift(g[2][nPion-1], 0.015);

    for(int i = 0; i < N_STARGETS; i++){
      delete histSolid[i];
      delete histLiquid[i];
      delete histBroadening[i];
    }
  } // End number pion event loop

  inputFile->Close();

  for(int k = 0; k < 3; k++) {
    g[k][0]->SetMarkerStyle(1);
    g[k][1]->SetMarkerStyle(8);
    g[k][1]->SetMarkerSize(0.95);
    g[k][2]->SetMarkerStyle(27);
    g[k][2]->SetMarkerSize(1.7);
  }

  TLegend* legendTarget = GenTLegendTarget(g);
  TLegend* legendNPion  = GenTLegendNPion(g);

  TPad* p = new TPad("p", "p", 0, 0, 1, 1);
  p->SetLeftMargin(0.13);
  p->SetTopMargin(0.06);
  p->SetGridx(1);
  p->SetGridy(1);
  p->Draw();
  p->cd();

  TMultiGraph* mg = new TMultiGraph();
  for(int i = 0; i < N_STARGETS; i++) {
    for(int j = 0; j < N_PION; j++) {
        mg->Add(g[i][j]);
    }
  }

  mg->Draw("APE0");

  gStyle->SetTitleFont(62, "XY");
  gStyle->SetTitleSize(0.04, "XY");

  mg->GetYaxis()->SetRangeUser(0.,0.18);
//mg->GetYaxis()->SetRangeUser(0.,0.45);;
  //mg->GetXaxis()->SetRangeUser(0.,0.75);
  mg->GetXaxis()->SetRangeUser(Zh_MIN, Zh_MAX);
  mg->GetXaxis()->SetTitle("Zh_{sum}");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  mg->Draw("APE0");
  legendTarget->Draw();
  legendNPion->Draw();

  c->Print(plotDirectory + "Pt_broad_Zh.pdf");
  delete c;
}

void PtBroadeningFullIntegrated(TString inputDirectory, TString plotDirectory) {

  TFile* inputFile = new TFile(inputDirectory + "meanPt2.root", "READ");

  TGraphErrors* g[N_STARGETS][N_PION];

  for(int nPion = 1; nPion <= N_PION; nPion++) {

    TH1F* histSolid[N_STARGETS]; TH1F* histLiquid[N_STARGETS]; TH1F* histBroadening[N_STARGETS];

    histSolid[0] = (TH1F*) inputFile->Get(Form("meanPt2_C_%i_CLEAN_INTERPOLATED",nPion));
    histSolid[1] = (TH1F*) inputFile->Get(Form("meanPt2_Fe_%i_CLEAN_INTERPOLATED",nPion));
    histSolid[2] = (TH1F*) inputFile->Get(Form("meanPt2_Pb_%i_CLEAN_INTERPOLATED",nPion));

    histLiquid[0] = (TH1F*) inputFile->Get(Form("meanPt2_DC_%i_CLEAN_INTERPOLATED",nPion));
    histLiquid[1] = (TH1F*) inputFile->Get(Form("meanPt2_DFe_%i_CLEAN_INTERPOLATED",nPion));
    histLiquid[2] = (TH1F*) inputFile->Get(Form("meanPt2_DPb_%i_CLEAN_INTERPOLATED",nPion));

    for(int i = 0 ; i < 3 ; i++) {
      //Calculate the Broadening (subtract of the means)
      histBroadening[i] = new TH1F(Form("histBroadening_%i",i),"",1,Zh_MIN,Zh_MAX);
      histBroadening[i]->Add(histSolid[i],histLiquid[i],1,-1);
    }

    // Set the points in TGraphErrors
    g[0][nPion-1] = new TGraphErrors();
    g[0][nPion-1]->SetPoint(1,TMath::Power(12.01,1./3.),histBroadening[0]->GetBinContent(1));
    g[0][nPion-1]->SetPointError(1,0,histBroadening[0]->GetBinError(1));

    g[1][nPion-1] = new TGraphErrors();
    g[1][nPion-1]->SetPoint(1,TMath::Power(55.845,1./3.),histBroadening[1]->GetBinContent(1));
    g[1][nPion-1]->SetPointError(1,0,histBroadening[1]->GetBinError(1));

    g[2][nPion-1] = new TGraphErrors();
    g[2][nPion-1]->SetPoint(1,TMath::Power(207.2,1./3.),histBroadening[2]->GetBinContent(1));
    g[2][nPion-1]->SetPointError(1,0,histBroadening[2]->GetBinError(1));

    //Set a color for each target
    g[0][nPion-1]->SetMarkerColor(kRed);
    g[0][nPion-1]->SetLineColor(kRed);
    g[1][nPion-1]->SetMarkerColor(kBlue);
    g[1][nPion-1]->SetLineColor(kBlue);
    g[2][nPion-1]->SetMarkerColor(kBlack);
    g[2][nPion-1]->SetLineColor(kBlack);

    // SetErrorXNull(g[0][nPion-1]);
    // SetErrorXNull(g[1][nPion-1]);
    // SetErrorXNull(g[2][nPion-1]);

    for(int i = 0; i < N_STARGETS; i++) {
      delete histSolid[i];
      delete histLiquid[i];
      delete histBroadening[i];
    }
  } // End number pion event loop

  inputFile->Close();

  for(int k = 0; k < N_STARGETS; k++) {
    g[k][0]->SetMarkerStyle(1);
    g[k][1]->SetMarkerStyle(8);
    g[k][1]->SetMarkerSize(.6);
    g[k][2]->SetMarkerStyle(27);
    g[k][2]->SetMarkerSize(1.2);
  }

  TCanvas* c = new TCanvas("c", "", 800, 600);
  c->Draw();

  // Create the legends
  TLegend* legendTarget = GenTLegendTarget(g);
  TLegend* legendNPion  = GenTLegendNPion(g);

  TPad* p = new TPad("p", "p", 0, 0, 1, 1);
  p->SetLeftMargin(0.13);
  p->SetTopMargin(0.06);
  p->SetGridx(1);
  p->SetGridy(1);
  p->Draw();
  p->cd();

  TMultiGraph* mg = new TMultiGraph();

  for(int i = 0; i < N_STARGETS; i++) {
    for(int j = 0; j < N_PION; j++) {
        mg->Add(g[i][j]);
    }
  }

  mg->Draw("APE");

  gStyle->SetTitleFont(62,"XY");
  gStyle->SetTitleSize(0.04,"XY");

  mg->GetYaxis()->SetRangeUser(0.,0.1);
  //mg->GetYaxis()->SetRangeUser(0.012,0.08);
  //mg->GetYaxis()->SetRangeUser(0.,0.2);
  mg->GetXaxis()->SetRangeUser(2, 6);
  mg->GetXaxis()->SetTitle("A^{1/3}");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  mg->Draw("APE");
  legendTarget->Draw();
  legendNPion->Draw();

  c->Print(plotDirectory + "Pt_broad_FullIntegrated.pdf");
  delete c;
}
