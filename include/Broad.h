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
#include "TLatex.h"


char DC[3] = {'D', 'C', '\0' }; char DFe[4] = {'D', 'F', 'e', '\0'}; char DPb[4] = {'D', 'P', 'b', '\0'};

char C[2]  = {'C', '\0'}; char Fe[3]  = {'F', 'e', '\0'}; char Pb[3]  = {'P', 'b', '\0'};

// If the histogram if empty return 1 if not return 0
int EmptyHist(TH1F* h) {

  int empty = 0;
  for(int i = 1 ; i <= h->GetNbinsX() ; i++) {
    if(h->GetBinContent(i) == 0){ empty++; }
  }
  if(empty == h->GetNbinsX()) { return 1; }
  else { return 0; }
}

//Integrate the PhiPQ histograms and generate a Pt2 histogram for each Q2, Nu, Zh bin
void PhiIntegration(TFile* inputFile, TFile* outputFile, char target[]) {

  std::cout << Form("Target %s",target) << std::endl;
  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
    for(int Q2Counter = 0 ; Q2Counter < N_Q2; Q2Counter++) { // Loops in every Q2 bin
      for(int NuCounter = 0 ; NuCounter < N_Nu; NuCounter++) { // Loops in every Nu bin
        for(int ZhCounter = 0 ; ZhCounter < N_Zh; ZhCounter++) { // Loops in every Zh bin;
          // For every bin of Q2, Nu, Zh and number of pion, generate a histogram  for Pt2
	        TH1F* histPt2 = new TH1F(Form("corr_data_Pt2_%s_%i%i%i_%i_CLEAN_INTERPOLATED", target, Q2Counter, NuCounter, ZhCounter, nPion), "", N_Pt2, Pt2_MIN, Pt2_MAX);
          for(int Pt2Counter = 0 ; Pt2Counter < N_Pt2 ; Pt2Counter++) { // Loops in every Pt2 bin
            // Take the his for this bins
	          TH1F* histPhi = (TH1F*) inputFile->Get(Form("DataCorr2_%s_%i%i%i%i_%i", target, Q2Counter, NuCounter, ZhCounter, Pt2Counter, nPion));
            // If the histogram is null or empty skip this Pt2 bin
            if(histPhi == NULL){ continue; }
	          if(EmptyHist(histPhi) == 1){ continue; }
            // Integrate on PhiPQ
	          // double errors_array[1] = {};
            double error;
	          double integral = histPhi->IntegralAndError(1, N_Phi, error);
	          // Save the value in the Pt2 histogram
	          histPt2->SetBinContent(Pt2Counter + 1, integral);
	          // histPt2->SetBinError(Pt2Counter + 1, errors_array[0]);
            histPt2->SetBinError(Pt2Counter + 1, error);
            delete histPhi;
	        }// End Pt2 loop
          // If the histogram if not empty, save it
	        if(EmptyHist(histPt2) == 0) {
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

// Integrate the Pt2 histograms for Nu and Zh bins
void NuZhIntegration(TFile* inputFile, TFile* outputFile, char target[]) {

  for(int nPion = 1; nPion <= N_PION ; nPion++) { // Loops in every number of pion
    // Generate a histogram to save Q2 for every number of pion in the final event
    TH1F* histQ2 = new TH1F(Form("corr_data_%s_%i_Q2", target, nPion), "", N_Q2, Q2_BINS);
    for(int Q2Counter = 0 ; Q2Counter < N_Q2 ; Q2Counter++) { // Loops in every Q2 bin
      TH1F* histPt2Integrated = new TH1F(Form("corr_data_Pt2_%s_%i", target, Q2Counter), "", N_Pt2, Pt2_MIN, Pt2_MAX);
      for(int NuCounter = 0 ; NuCounter < N_Nu ; NuCounter++) { // Loops in every Nu bin
        // Starts in the second bin because there is not broadening iz Zh<0.2
	      for(int ZhCounter = ZH_SUM ; ZhCounter < N_Zh ; ZhCounter++) { // Loops in every Zh bin
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

// Integrate the Pt2 histograms for Q2 and Nu bins
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

// Integrate the Pt2 histograms for Q2 and Zh bins
void Q2ZhIntegration(TFile* inputFile, TFile* outputFile, char target[]) {

  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
    // Generate a histogram to save Nu for every number of pion in the final event
    TH1F* histNu = new TH1F(Form("corr_data_%s_Nu", target), "", N_Nu, Nu_BINS);
    for(int NuCounter = 0 ; NuCounter < N_Nu ; NuCounter++) { // Loops in every Nu bin
      TH1F* histPt2Integrated = new TH1F(Form("corr_data_Pt2_%s_%i", target, NuCounter), "", N_Pt2, Pt2_MIN, Pt2_MAX);
      // Starts in the second bin because there is not broadening in Zh<0.2
      for(int ZhCounter = ZH_SUM ; ZhCounter < N_Zh ; ZhCounter++) { // Loops in every Zh bin
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

// Integrate the Pt2 histograms for Q2, Nu and Zh bins
void Q2NuZhIntegration(TFile* inputFile, TFile* outputFile, char target[]) {

  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
    // Generate two histograms for every number of pion in the final event
    // To save the mean of Pt2
    TH1F* hist = new TH1F(Form("corr_data_%s", target), "", 1, 0, 1);
    // To save the sum of the histograms
    TH1F* histPt2Integrated = new TH1F(Form("corr_data_%s_Pt2", target), "", N_Pt2, Pt2_MIN, Pt2_MAX);
    // Starts in the second bin because there is not broadening iz Zh<0.2
    for(int ZhCounter = ZH_SUM ; ZhCounter < N_Zh ; ZhCounter++) { // Loops in every Zh bin
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

void SetErrorXNull(TGraphErrors* g, int Ntarget) {

  double* errors_Y = g->GetEY();
  for(int point = 0 ; point < Ntarget ; point++) {
    g->SetPointError(point,0,errors_Y[point]);
  }

}

void SetXShift(TGraphErrors* g, double shift, int Ntarget) {

  double* content_Y = g->GetY();
  double* content_X = g->GetX();

  for(int point = 0; point < Ntarget ; point++) {
    g->SetPoint(point, content_X[point] + shift , content_Y[point]);
  }

}

TLegend* GenTLegendTarget(TGraphErrors* g[N_STARGETS][N_PION + 1], float pos = 0.17) {

  TLegend* legend = new TLegend(pos, 0.75, pos + 0.13, 0.92, "", "NDC");
  legend->AddEntry(g[0][0],"C","lpf");
  legend->AddEntry(g[1][0],"Fe","lpf");
  legend->AddEntry(g[2][0],"Pb","lpf");
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(62);
  legend->SetTextSize(0.04);
  return legend;

}

TLegend* GenTLegendNPion(TGraphErrors* g[N_STARGETS][N_PION + 1], int m, float pos = 0.28) {

  TLegend* legend = new TLegend(pos, 0.75, pos + 0.07, 0.92, "", "NDC");
  //legend->AddEntry((TObject*)0, "Final state", "");
  legend->AddEntry(g[2][0],"1 #pi +","lpf");
  legend->AddEntry(g[2][1],"2 #pi +","lpf");
  legend->AddEntry(g[2][2],"3 #pi +","lpf");
  if(m == 4){ legend->AddEntry(g[2][3],"All Data","lpf");}
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(62);
  legend->SetTextSize(0.04);
  return legend;

}

// TLegend* GenTLegendTarget(TGraphErrors* g[N_STARGETS][N_PION]) {
//
//   TLegend* legend = new TLegend(0.22,0.75,0.3,0.92,"","NDC");
//   legend->AddEntry(g[0][0],"C","lpf");
//   legend->AddEntry(g[1][0],"Fe","lpf");
//   legend->AddEntry(g[2][0],"Pb","lpf");
//   legend->SetFillStyle(0);
//   legend->SetBorderSize(0);
//   legend->SetTextFont(62);
//   legend->SetTextSize(0.04);
//   return legend;
//
// }
//
TLegend* GenTLegendNPion(TGraphErrors* g[N_STARGETS][N_PION+1]) {

  TLegend* legend = new TLegend(0.35,0.75,0.3,0.92,"","NDC");
  //legend->AddEntry((TObject*)0, "Final state", "");
  legend->AddEntry(g[2][0],"1 #pi +","lpf");
  legend->AddEntry(g[2][1],"2 #pi +","lpf");
  legend->AddEntry(g[2][2],"3 #pi +","lpf");
  legend->AddEntry(g[2][3],"All Data","lpf");
  legend->SetFillStyle(0);
  legend->SetBorderSize(0);
  legend->SetTextFont(62);
  legend->SetTextSize(0.04);
  return legend;

}

TGraphErrors* TH1TOTGraph(TH1 *h1) {

  if (!h1) std::cout << "TH1TOTGraph: histogram not found !" << std::endl;

  TGraphErrors* g1= new TGraphErrors();
  Double_t x, y, ex, ey;
  for (Int_t i = 1; i <= h1->GetNbinsX(); i++) {
    y = h1->GetBinContent(i);
    ey = h1->GetBinError(i);
    x = h1->GetBinCenter(i);
    // ex=h1->GetBinWidth(i);
    ex = 0.5*h1->GetBinWidth(i);
    g1->SetPoint(i-1, x, y);
    g1->SetPointError(i-1, ex, ey);
 }
 return g1;

}

// Calculate the Pt broadening and plot in funcion of the Nu bins
void PtBroadeningNuIntegrated(TString inputDirectory, TString plotDirectory) {

  TFile* inputFile = new TFile(inputDirectory + Form("meanPt2_Nu_%i.root", N_Nu), "READ");
  TFile* inputFileEst = new TFile(inputDirectory + Form("meanPt2_Nu_%i-Esteban.root", N_Nu), "READ");

  TGraphErrors* g[N_STARGETS][N_PION+1];
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
      SetErrorXNull(g[i][nPion-1], N_Nu);
    }

    g[0][nPion-1]->SetMarkerColor(kRed);
    g[0][nPion-1]->SetLineColor(kRed);
    g[1][nPion-1]->SetLineColor(kBlue);
    g[1][nPion-1]->SetMarkerColor(kBlue);
    g[2][nPion-1]->SetMarkerColor(kBlack);
    g[2][nPion-1]->SetLineColor(kBlack);

    // SetXShift(g[0][nPion-1], -0.035, N_Nu);
    // SetXShift(g[2][nPion-1], 0.035,  N_Nu);

    for(int i = 0; i < 3; i++) {
      delete histSolid[i];
      delete histLiquid[i];
      delete histBroadening[i];
    }
  }// End number pion event loop

  inputFile->Close();

  // ***** Add the standar pt2 Broadening pt2 broadening (Esteban results)********/////

  TH1F* histSolid[N_STARGETS]; TH1F* histLiquid[N_STARGETS]; TH1F* histBroadening[N_STARGETS];

  histSolid[0]  = (TH1F*) inputFileEst->Get("meanPt2_C_Nu_CLEAN_INTERPOLATED");
  histSolid[1]  = (TH1F*) inputFileEst->Get("meanPt2_Fe_Nu_CLEAN_INTERPOLATED");
  histSolid[2]  = (TH1F*) inputFileEst->Get("meanPt2_Pb_Nu_CLEAN_INTERPOLATED");

  histLiquid[0] = (TH1F*) inputFileEst->Get("meanPt2_DC_Nu_CLEAN_INTERPOLATED");
  histLiquid[1] = (TH1F*) inputFileEst->Get("meanPt2_DFe_Nu_CLEAN_INTERPOLATED");
  histLiquid[2] = (TH1F*) inputFileEst->Get("meanPt2_DPb_Nu_CLEAN_INTERPOLATED");

  for(int i = 0 ; i < 3 ; i++) {
    //Calculate the Broadening (subtract of the means)
    histBroadening[i] = new TH1F(Form("histBroadening_%i", i), "", N_Nu, Nu_BINS_E);
    histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
  }

  for(int i = 0 ; i < N_STARGETS ; i++){
    g[i][N_PION] = (TGraphErrors*) TH1TOTGraph(histBroadening[i]);
    SetErrorXNull(g[i][N_PION], N_Nu);
  }
  //Set a color for each target
  g[0][N_PION]->SetMarkerColor(kRed);
  g[0][N_PION]->SetLineColor(kRed);
  g[1][N_PION]->SetMarkerColor(kBlue);
  g[1][N_PION]->SetLineColor(kBlue);
  g[2][N_PION]->SetMarkerColor(kBlack);

  SetXShift(g[0][0], -0.075 - 0.030, N_Nu);
  SetXShift(g[1][0], - 0.030, N_Nu);
  SetXShift(g[2][0], 0.075 +  0.040, N_Nu);

  SetXShift(g[0][1], -0.075 , N_Nu);
  SetXShift(g[2][1], 0.075 , N_Nu);

  SetXShift(g[0][2], -0.075  - 0.050, N_Nu);
  SetXShift(g[1][2], -0.030, N_Nu);
  SetXShift(g[2][2], 0.075 ,  N_Nu);


  SetXShift(g[0][N_PION], -0.075, N_Nu);
  SetXShift(g[2][N_PION], 0.055,  N_Nu);

  inputFileEst->Close();
  // ***** Add the standar pt2 Broadening pt2 broadening (Esteban results)********/////

  for(int k = 0; k < N_PION; k++) {
    g[k][0]->SetMarkerStyle(4);
    g[k][0]->SetMarkerSize(.53);
    g[k][1]->SetMarkerStyle(8);
    g[k][1]->SetMarkerSize(.6);
    g[k][2]->SetMarkerStyle(27);
    g[k][2]->SetMarkerSize(1.4);
    g[k][3]->SetMarkerStyle(22);
    g[k][3]->SetMarkerSize(1.);
  }


  TPad* p = new TPad("p","p",0,0,1,1);
  p->SetLeftMargin(0.13);
  p->SetTopMargin(0.055);
  p->SetGridx(1);
  p->SetGridy(1);
  p->Draw();
  p->cd();

  TMultiGraph* mg = new TMultiGraph();
  for(int i = 0; i < N_STARGETS; i++) {
    for(int j = 0; j < N_PION+1; j++) {
        mg->Add(g[i][j]);
    }
  }
  mg->Draw("APE0");

  gStyle->SetTitleFont(62,"XY");
  gStyle->SetTitleSize(0.04,"XY");

  // mg->GetYaxis()->SetRangeUser(0.,0.059);
  mg->GetYaxis()->SetRangeUser(0.,0.049);
  mg->GetXaxis()->SetRangeUser(2.1,4.3);
  mg->GetXaxis()->SetTitle("#nu[GeV]");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  mg->Draw("APE0");

  TLegend* legendTarget = GenTLegendTarget(g);
  TLegend* legendNPion= GenTLegendNPion(g);

  legendTarget->Draw();
  legendNPion->Draw();

  c->Print(plotDirectory + "Pt_broad_Nu.pdf");
  delete c;

}

// Calculate the Pt broadening and plot in funcion of the Q2 bins
void PtBroadeningQ2Integrated(TString inputDirectory, TString plotDirectory) {

  TFile* inputFile = new TFile(inputDirectory + Form("meanPt2_Q2_%i.root", N_Q2), "READ");
  TFile* inputFileEst = new TFile(inputDirectory + Form("meanPt2_Q2_%i-Esteban.root", N_Q2), "READ");

  TCanvas* c = new TCanvas("c", "", 800, 600);
  c->Draw();

  TGraphErrors* g[N_STARGETS][N_PION+1];

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
      SetErrorXNull(g[i][nPion-1], N_Q2);
    }

    g[0][nPion-1]->SetMarkerColor(kRed);
    g[1][nPion-1]->SetMarkerColor(kBlue);
    g[0][nPion-1]->SetLineColor(kRed);
    g[1][nPion-1]->SetLineColor(kBlue);
    g[2][nPion-1]->SetMarkerColor(kBlack);
    g[2][nPion-1]->SetLineColor(kBlack);

    // SetXShift(g[0][nPion-1], -0.03 - 0.020*nPion, N_Q2);
    // SetXShift(g[2][nPion-1], 0.03 + 0.020*nPion,  N_Q2);

    for(int i = 0; i < N_STARGETS; i++){
      delete histSolid[i];
      delete histLiquid[i];
      delete histBroadening[i];
    }
  } // End number pion event loop

  inputFile->Close();

  // ***** Add the standar pt2 Broadening pt2 broadening (Esteban results)********/////

  TH1F* histSolid[N_STARGETS]; TH1F* histLiquid[N_STARGETS]; TH1F* histBroadening[N_STARGETS];

  histSolid[0]  = (TH1F*) inputFileEst->Get("meanPt2_C_Q2_CLEAN_INTERPOLATED");
  histSolid[1]  = (TH1F*) inputFileEst->Get("meanPt2_Fe_Q2_CLEAN_INTERPOLATED");
  histSolid[2]  = (TH1F*) inputFileEst->Get("meanPt2_Pb_Q2_CLEAN_INTERPOLATED");

  histLiquid[0] = (TH1F*) inputFileEst->Get("meanPt2_DC_Q2_CLEAN_INTERPOLATED");
  histLiquid[1] = (TH1F*) inputFileEst->Get("meanPt2_DFe_Q2_CLEAN_INTERPOLATED");
  histLiquid[2] = (TH1F*) inputFileEst->Get("meanPt2_DPb_Q2_CLEAN_INTERPOLATED");

  for(int i = 0 ; i < 3 ; i++) {
    //Calculate the Broadening (subtract of the means)
    histBroadening[i] = new TH1F(Form("histBroadening_%i", i), "", N_Q2, Q2_BINS_E);
    histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
  }

  for(int i = 0 ; i < N_STARGETS ; i++){
    g[i][N_PION] = (TGraphErrors*) TH1TOTGraph(histBroadening[i]);
    SetErrorXNull(g[i][N_PION], N_Q2);
  }

  //Set a color for each target
  g[0][N_PION]->SetMarkerColor(kRed);
  g[0][N_PION]->SetLineColor(kRed);
  g[1][N_PION]->SetMarkerColor(kBlue);
  g[1][N_PION]->SetLineColor(kBlue);
  g[2][N_PION]->SetMarkerColor(kBlack);

  SetXShift(g[0][0], -0.075 - 0.030, N_Q2);
  SetXShift(g[1][0], - 0.030, N_Q2);
  SetXShift(g[2][0], 0.075 +  0.030, N_Q2);

  SetXShift(g[0][1], -0.075 , N_Q2);
  SetXShift(g[2][1], 0.075 , N_Q2);

  SetXShift(g[0][2], -0.075  - 0.050, N_Q2);
  SetXShift(g[1][2], -0.030, N_Q2);
  SetXShift(g[2][2], 0.075 ,  N_Q2);


  SetXShift(g[0][N_PION], -0.055, N_Q2);
  //SetXShift(g[1][N_PION], 0.095,  N_Q2);
  SetXShift(g[2][N_PION], 0.055,  N_Q2);

  inputFileEst->Close();
  // ***** Add the standar pt2 Broadening pt2 broadening (Esteban results)********/////


  TPad* p = new TPad("p","p",0,0,1,1);
  p->SetGridx(1);
  p->SetGridy(1);
  p->SetLeftMargin(0.13);
  p->SetTopMargin(0.06);
  p->Draw();
  p->cd();

  for(int k = 0; k < N_STARGETS; k++) {
    g[k][0]->SetMarkerStyle(4);
    g[k][0]->SetMarkerSize(.53);
    g[k][1]->SetMarkerStyle(8);
    g[k][1]->SetMarkerSize(.6);
    g[k][2]->SetMarkerStyle(27);
    g[k][2]->SetMarkerSize(1.2);
    g[k][3]->SetMarkerStyle(22);
    g[k][3]->SetMarkerSize(1.);
  }

  TMultiGraph* mg = new TMultiGraph();
  for(int i = 0; i < N_STARGETS; i++){
    for(int j = 0; j < N_PION + 1; j++){
        mg->Add(g[i][j]);
    }
  }

  mg->Draw("APE0");

  gStyle->SetTitleFont(62,"XY");
  gStyle->SetTitleSize(0.04,"XY");

  //mg->GetYaxis()->SetRangeUser(0.,0.059);
  mg->GetYaxis()->SetRangeUser(0.,0.048);
  mg->GetXaxis()->SetRangeUser(.9,4.1);
  mg->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
  mg->GetXaxis()->CenterTitle();
  mg->GetXaxis()->SetTitleOffset(1.1);

  mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
  mg->GetYaxis()->CenterTitle();
  mg->GetYaxis()->SetTitleOffset(1.3);

  mg->Draw("APE0");
  TLegend* legendTarget = GenTLegendTarget(g);
  TLegend* legendNPion = GenTLegendNPion(g, 0.42);
  legendTarget->Draw();
  legendNPion->Draw();

  c->Print(plotDirectory + "Pt_broad_Q2.pdf");
  delete c;

}

// Calculate the Pt broadening and plot in funcion of the Zh bins
void PtBroadeningZhIntegrated(TString inputDirectory,  TString plotDirectory) {

  TFile* inputFile = new TFile(inputDirectory + Form("meanPt2_Zh_%i.root", N_Zh), "READ");
  TFile* inputFileEst = new TFile(inputDirectory + Form("meanPt2_Zh_%i-Esteban.root", 8), "READ");

  TCanvas* c = new TCanvas("c", "", 800, 600);
  c->Draw();
  TGraphErrors* g[N_STARGETS][N_PION+1];

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
      SetErrorXNull(g[i][nPion-1], N_Zh);
    }

    g[0][nPion-1]->SetMarkerColor(kRed);
    g[0][nPion-1]->SetLineColor(kRed);
    g[1][nPion-1]->SetMarkerColor(kBlue);
    g[1][nPion-1]->SetLineColor(kBlue);
    g[2][nPion-1]->SetMarkerColor(kBlack);
    g[2][nPion-1]->SetLineColor(kBlack);

    SetXShift(g[0][nPion-1], -0.015, N_Zh);
    SetXShift(g[2][nPion-1], 0.015,  N_Zh);

    for(int i = 0; i < N_STARGETS; i++){
      delete histSolid[i];
      delete histLiquid[i];
      delete histBroadening[i];
    }
  } // End number pion event loop

  inputFile->Close();

  // ***** Add the standar pt2 Broadening pt2 broadening (Esteban results)********/////

  TH1F* histSolid[N_STARGETS]; TH1F* histLiquid[N_STARGETS]; TH1F* histBroadening[N_STARGETS];

  histSolid[0]  = (TH1F*) inputFileEst->Get("meanPt2_C_CLEAN_INTERPOLATED");
  histSolid[1]  = (TH1F*) inputFileEst->Get("meanPt2_Fe_CLEAN_INTERPOLATED");
  histSolid[2]  = (TH1F*) inputFileEst->Get("meanPt2_Pb_CLEAN_INTERPOLATED");

  histLiquid[0] = (TH1F*) inputFileEst->Get("meanPt2_DC_CLEAN_INTERPOLATED");
  histLiquid[1] = (TH1F*) inputFileEst->Get("meanPt2_DFe_CLEAN_INTERPOLATED");
  histLiquid[2] = (TH1F*) inputFileEst->Get("meanPt2_DPb_CLEAN_INTERPOLATED");

  for(int i = 0 ; i < 3 ; i++) {
    //Calculate the Broadening (subtract of the means)
    histBroadening[i] = new TH1F(Form("histBroadening_%i", i), "", 8, Zh_BINS_E);
    histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
    }

  for(int i = 0 ; i < N_STARGETS ; i++){
    g[i][N_PION] = (TGraphErrors*) TH1TOTGraph(histBroadening[i]);
    SetErrorXNull(g[i][N_PION], 8);
  }
        //Set a color for each target
  g[0][N_PION]->SetMarkerColor(kRed);
  g[0][N_PION]->SetLineColor(kRed);
  g[1][N_PION]->SetMarkerColor(kBlue);
  g[1][N_PION]->SetLineColor(kBlue);
  g[2][N_PION]->SetMarkerColor(kBlack);
  g[2][N_PION]->SetLineColor(kBlack);

  SetErrorXNull(g[0][N_PION], N_STARGETS);
  SetErrorXNull(g[1][N_PION], N_STARGETS);
  SetErrorXNull(g[2][N_PION], N_STARGETS);


  inputFileEst->Close();

    // ***** Add the standar pt2 Broadening pt2 broadening (Esteban results)********/////

  for(int k = 0; k < 3; k++) {
    g[k][0]->SetMarkerStyle(4);
    g[k][0]->SetMarkerSize(.53);
    g[k][1]->SetMarkerStyle(8);
    g[k][1]->SetMarkerSize(0.95);
    g[k][2]->SetMarkerStyle(27);
    g[k][2]->SetMarkerSize(1.7);
    g[k][3]->SetMarkerStyle(22);
    g[k][3]->SetMarkerSize(1.2);
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
    for(int j = 0; j < N_PION + 1 ; j++) {
        mg->Add(g[i][j]);
    }
  }

  mg->Draw("APE0");

  gStyle->SetTitleFont(62, "XY");
  gStyle->SetTitleSize(0.04, "XY");

  mg->GetYaxis()->SetRangeUser(0.,0.072);
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

// Calculate the Pt broadening and plot in funcion of the size of the target
void PtBroadeningFullIntegrated(TString inputDirectory, TString plotDirectory) {

  TFile* inputFile = new TFile(inputDirectory + "meanPt2.root", "READ");
  TFile* inputFileEst = new TFile(inputDirectory + "meanPt2-Esteban.root", "READ");

  TGraphErrors* g[N_STARGETS][N_PION + 1];

  for(int nPion = 1; nPion <= N_PION; nPion++) {

    TH1F* histSolid[N_STARGETS]; TH1F* histLiquid[N_STARGETS]; TH1F* histBroadening[N_STARGETS];

    histSolid[0]  = (TH1F*) inputFile->Get(Form("meanPt2_C_%i_CLEAN_INTERPOLATED",   nPion));
    histSolid[1]  = (TH1F*) inputFile->Get(Form("meanPt2_Fe_%i_CLEAN_INTERPOLATED",  nPion));
    histSolid[2]  = (TH1F*) inputFile->Get(Form("meanPt2_Pb_%i_CLEAN_INTERPOLATED",  nPion));

    histLiquid[0] = (TH1F*) inputFile->Get(Form("meanPt2_DC_%i_CLEAN_INTERPOLATED",  nPion));
    histLiquid[1] = (TH1F*) inputFile->Get(Form("meanPt2_DFe_%i_CLEAN_INTERPOLATED", nPion));
    histLiquid[2] = (TH1F*) inputFile->Get(Form("meanPt2_DPb_%i_CLEAN_INTERPOLATED", nPion));

    for(int i = 0 ; i < N_PION ; i++) {
      //Calculate the Broadening (subtract of the means)
      histBroadening[i] = new TH1F(Form("histBroadening_%i", i), "", 1, Zh_MIN, Zh_MAX);
      histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
    }

    // Set the points in TGraphErrors
    g[0][nPion-1] = new TGraphErrors();
    g[0][nPion-1]->SetPoint(1, TMath::Power(12.01,1./3.), histBroadening[0]->GetBinContent(1));
    g[0][nPion-1]->SetPointError(1, 0, histBroadening[0]->GetBinError(1));

    g[1][nPion-1] = new TGraphErrors();
    g[1][nPion-1]->SetPoint(1, TMath::Power(55.845,1./3.), histBroadening[1]->GetBinContent(1));
    g[1][nPion-1]->SetPointError(1, 0,histBroadening[1]->GetBinError(1));

    g[2][nPion-1] = new TGraphErrors();
    g[2][nPion-1]->SetPoint(1, TMath::Power(207.2,1./3.), histBroadening[2]->GetBinContent(1));
    g[2][nPion-1]->SetPointError(1, 0, histBroadening[2]->GetBinError(1));

    //Set a color for each target
    g[0][nPion-1]->SetMarkerColor(kRed);
    g[0][nPion-1]->SetLineColor(kRed);
    g[1][nPion-1]->SetMarkerColor(kBlue);
    g[1][nPion-1]->SetLineColor(kBlue);
    g[2][nPion-1]->SetMarkerColor(kBlack);
    g[2][nPion-1]->SetLineColor(kBlack);

    SetErrorXNull(g[0][nPion-1], N_STARGETS);
    SetErrorXNull(g[1][nPion-1], N_STARGETS);
    SetErrorXNull(g[2][nPion-1], N_STARGETS);

    for(int i = 0; i < N_STARGETS; i++) {
      delete histSolid[i];
      delete histLiquid[i];
      delete histBroadening[i];
    }
  } // End number pion event loop
  inputFile->Close();

  // ***** Add the standar pt broadening (Esteban results)********/////

  TH1F* histSolid[N_STARGETS]; TH1F* histLiquid[N_STARGETS]; TH1F* histBroadening[N_STARGETS];

  histSolid[0]  = (TH1F*) inputFileEst->Get("meanPt2_C_CLEAN_INTERPOLATED");
  histSolid[1]  = (TH1F*) inputFileEst->Get("meanPt2_Fe_CLEAN_INTERPOLATED");
  histSolid[2]  = (TH1F*) inputFileEst->Get("meanPt2_Pb_CLEAN_INTERPOLATED");

  histLiquid[0] = (TH1F*) inputFileEst->Get("meanPt2_DC_CLEAN_INTERPOLATED");
  histLiquid[1] = (TH1F*) inputFileEst->Get("meanPt2_DFe_CLEAN_INTERPOLATED");
  histLiquid[2] = (TH1F*) inputFileEst->Get("meanPt2_DPb_CLEAN_INTERPOLATED");

  for(int i = 0 ; i < 3 ; i++) {
    //Calculate the Broadening (subtract of the means)
    histBroadening[i] = new TH1F(Form("histBroadening_%i",i), "", 1, Zh_MIN ,Zh_MAX);
    histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
  }

  // Set the points in TGraphErrors
  g[0][N_PION] = new TGraphErrors();
  g[0][N_PION]->SetPoint(1,TMath::Power(12.01,1./3.),histBroadening[0]->GetBinContent(1));
  g[0][N_PION]->SetPointError(1,0,histBroadening[0]->GetBinError(1));

  g[1][N_PION] = new TGraphErrors();
  g[1][N_PION]->SetPoint(1,TMath::Power(55.845,1./3.),histBroadening[1]->GetBinContent(1));
  g[1][N_PION]->SetPointError(1,0,histBroadening[1]->GetBinError(1));

  g[2][N_PION] = new TGraphErrors();
  g[2][N_PION]->SetPoint(1,TMath::Power(207.2,1./3.),histBroadening[2]->GetBinContent(1));
  g[2][N_PION]->SetPointError(1,0,histBroadening[2]->GetBinError(1));

  //Set a color for each target
  g[0][N_PION]->SetMarkerColor(kRed);
  g[0][N_PION]->SetLineColor(kRed);
  g[1][N_PION]->SetMarkerColor(kBlue);
  g[1][N_PION]->SetLineColor(kBlue);
  g[2][N_PION]->SetMarkerColor(kBlack);

  SetErrorXNull(g[0][N_PION], N_STARGETS);
  SetErrorXNull(g[1][N_PION], N_STARGETS);
  SetErrorXNull(g[2][N_PION], N_STARGETS);


  inputFileEst->Close();
  // ***** Add the standar pt Broadening (Esteban results)********/////

  for(int k = 0; k < N_STARGETS; k++) {
    g[k][0]->SetMarkerStyle(4);
    g[k][0]->SetMarkerSize(.53);
    g[k][1]->SetMarkerStyle(8);
    g[k][1]->SetMarkerSize(.6);
    g[k][2]->SetMarkerStyle(27);
    g[k][2]->SetMarkerSize(1.2);
    g[k][3]->SetMarkerStyle(22);
    g[k][3]->SetMarkerSize(1.);
    SetXShift(g[k][0], -0.055, N_STARGETS);
    SetXShift(g[k][2], 0.055,  N_STARGETS);
    SetXShift(g[k][N_PION], 0.055, N_STARGETS);
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
    for(int j = 0; j < N_PION + 1; j++) {
        mg->Add(g[i][j]);
    }
  }

  mg->Draw("APE");

  gStyle->SetTitleFont(62,"XY");
  gStyle->SetTitleSize(0.04,"XY");

  mg->GetYaxis()->SetRangeUser(-0.004,0.042);
  mg->GetXaxis()->SetRangeUser(2, 6.3);
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


// void SetAxisMultiGraph(TMultiGraph* mg, char var[]) {
//
//   if (var[0] == 'Q') {
//     mg->GetYaxis()->SetRangeUser(0.,0.048);
//     mg->GetXaxis()->SetRangeUser(.9,4.1);
//     mg->GetXaxis()->SetTitle("Q^{2}[GeV^{2}]");
//     mg->GetXaxis()->CenterTitle();
//     mg->GetXaxis()->SetTitleOffset(1.1);
//     return;
//   }
//   if (var[0] == 'N') {
//     mg->GetYaxis()->SetRangeUser(0.,0.049);
//     mg->GetXaxis()->SetRangeUser(2.1,4.3);
//     mg->GetXaxis()->SetTitle("#nu[GeV]");
//     mg->GetXaxis()->CenterTitle();
//     mg->GetXaxis()->SetTitleOffset(1.1);
//     return;
//   }
//   if (var[0] == 'Z') {
//     mg->GetYaxis()->SetRangeUser(0.,0.072);
//     mg->GetXaxis()->SetRangeUser(Zh_MIN, Zh_MAX);
//     mg->GetXaxis()->SetTitle("Zh_{sum}");
//     mg->GetXaxis()->CenterTitle();
//     mg->GetXaxis()->SetTitleOffset(1.1);
//     return;
//   }
//
// }
//
// // Calculate the Pt broadening and plot in funcion of the Q2 bins
// void PtBroadeningVarIntegrated(char target[], int nTarget, TString inputDirectory, TString plotDirectory, bool AddAllData = false) {
//
//   TFile* inputFile    = new TFile(inputDirectory + Form("meanPt2_%s_%i.root", target, nTarget), "READ");
//
//   TCanvas* c = new TCanvas("c", "", 800, 600);
//   c->Draw();
//
//   int m = N_PION;
//   if(AddAllData) {m += 1;}
//
//   TGraphErrors* g[N_STARGETS][N_PION + 1];
//
//   for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
//     TH1F* histSolid[N_STARGETS];   TH1F* histLiquid[N_STARGETS]; TH1F* histBroadening[N_STARGETS];
//
//     //C
//     histSolid[0] = (TH1F*) inputFile->Get(Form("meanPt2_C_%i_%s_CLEAN_INTERPOLATED",   nPion, target));
//     histLiquid[0] = (TH1F*) inputFile->Get(Form("meanPt2_DC_%i_%s_CLEAN_INTERPOLATED", nPion, target));
//     //Fe
//     histSolid[1] = (TH1F*) inputFile->Get(Form("meanPt2_Fe_%i_%s_CLEAN_INTERPOLATED",   nPion, target));
//     histLiquid[1] = (TH1F*) inputFile->Get(Form("meanPt2_DFe_%i_%s_CLEAN_INTERPOLATED", nPion, target));
//     //Pb
//     histSolid[2] = (TH1F*) inputFile->Get(Form("meanPt2_Pb_%i_%s_CLEAN_INTERPOLATED",   nPion, target));
//     histLiquid[2] = (TH1F*) inputFile->Get(Form("meanPt2_DPb_%i_%s_CLEAN_INTERPOLATED", nPion, target));
//
//
//     for(int i = 0 ; i < N_STARGETS ; i++){
//       histBroadening[i] = new TH1F(Form("histBroadening_%i",i), "", N_Q2, Q2_BINS);
//       histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
//     }
//
//
//     for(int i = 0 ; i < N_STARGETS ; i++){
//       g[i][nPion-1] = (TGraphErrors*) TH1TOTGraph(histBroadening[i]);
//       SetErrorXNull(g[i][nPion-1], N_Q2);
//     }
//
//     g[0][nPion-1]->SetMarkerColor(kRed);
//     g[1][nPion-1]->SetMarkerColor(kBlue);
//     g[0][nPion-1]->SetLineColor(kRed);
//     g[1][nPion-1]->SetLineColor(kBlue);
//     g[2][nPion-1]->SetMarkerColor(kBlack);
//     g[2][nPion-1]->SetLineColor(kBlack);
//
//     // SetXShift(g[0][nPion-1], -0.03 - 0.020*nPion, N_Q2);
//     // SetXShift(g[2][nPion-1], 0.03 + 0.020*nPion,  N_Q2);
//
//     for(int i = 0; i < N_STARGETS; i++){
//       delete histSolid[i];
//       delete histLiquid[i];
//       delete histBroadening[i];
//     }
//   } // End number pion event loop
//
//   inputFile->Close();
//
//   if(AddAllData) { // Add the standar pt2 Broadening pt2 broadening (Esteban results)
//
//   TFile* inputFileEst = new TFile(inputDirectory + Form("meanPt2_%s_%i-Esteban.root",  target, nTarget), "READ");
//
//   TH1F* histSolid[N_STARGETS]; TH1F* histLiquid[N_STARGETS]; TH1F* histBroadening[N_STARGETS];
//
//   histSolid[0]  = (TH1F*) inputFileEst->Get("meanPt2_C_Q2_CLEAN_INTERPOLATED");
//   histSolid[1]  = (TH1F*) inputFileEst->Get("meanPt2_Fe_Q2_CLEAN_INTERPOLATED");
//   histSolid[2]  = (TH1F*) inputFileEst->Get("meanPt2_Pb_Q2_CLEAN_INTERPOLATED");
//
//   histLiquid[0] = (TH1F*) inputFileEst->Get("meanPt2_DC_Q2_CLEAN_INTERPOLATED");
//   histLiquid[1] = (TH1F*) inputFileEst->Get("meanPt2_DFe_Q2_CLEAN_INTERPOLATED");
//   histLiquid[2] = (TH1F*) inputFileEst->Get("meanPt2_DPb_Q2_CLEAN_INTERPOLATED");
//
//   for(int i = 0 ; i < 3 ; i++) {
//     //Calculate the Broadening (subtract of the means)
//     histBroadening[i] = new TH1F(Form("histBroadening_%i", i), "", N_Q2, Q2_BINS_E);
//     histBroadening[i]->Add(histSolid[i], histLiquid[i], 1, -1);
//   }
//
//   for(int i = 0 ; i < N_STARGETS ; i++){
//     g[i][N_PION] = (TGraphErrors*) TH1TOTGraph(histBroadening[i]);
//     SetErrorXNull(g[i][N_PION], N_Q2);
//   }
//
//   //Set a color for each target
//   g[0][N_PION]->SetMarkerColor(kRed);
//   g[0][N_PION]->SetLineColor(kRed);
//   g[1][N_PION]->SetMarkerColor(kBlue);
//   g[1][N_PION]->SetLineColor(kBlue);
//   g[2][N_PION]->SetMarkerColor(kBlack);
//
//   SetXShift(g[0][0], -0.075 - 0.030, N_Q2);
//   SetXShift(g[1][0], - 0.030, N_Q2);
//   SetXShift(g[2][0], 0.075 +  0.030, N_Q2);
//
//   SetXShift(g[0][1], -0.075 , N_Q2);
//   SetXShift(g[2][1], 0.075 , N_Q2);
//
//   SetXShift(g[0][2], -0.075  - 0.050, N_Q2);
//   SetXShift(g[1][2], -0.030, N_Q2);
//   SetXShift(g[2][2], 0.075 ,  N_Q2);
//
//
//   SetXShift(g[0][N_PION], -0.055, N_Q2);
//   //SetXShift(g[1][N_PION], 0.095,  N_Q2);
//   SetXShift(g[2][N_PION], 0.055,  N_Q2);
//
//   inputFileEst->Close();
//   // ***** Add the standar pt2 Broadening pt2 broadening (Esteban results)********/////
//   }
//
//   TPad* p = new TPad("p","p",0,0,1,1);
//   p->SetGridx(1);
//   p->SetGridy(1);
//   p->SetLeftMargin(0.13);
//   p->SetTopMargin(0.06);
//   p->Draw();
//   p->cd();
//
//   for(int k = 0; k < N_STARGETS; k++) {
//     g[k][0]->SetMarkerStyle(4);
//     g[k][0]->SetMarkerSize(.53);
//     g[k][1]->SetMarkerStyle(8);
//     g[k][1]->SetMarkerSize(.6);
//     g[k][2]->SetMarkerStyle(27);
//     g[k][2]->SetMarkerSize(1.2);
//     if(AddAllData) {
//       g[k][3]->SetMarkerStyle(22);
//       g[k][3]->SetMarkerSize(1.);
//     }
//   }
//
//   TMultiGraph* mg = new TMultiGraph();
//   for(int i = 0; i < N_STARGETS; i++){
//     for(int j = 0; j < m; j++){
//         mg->Add(g[i][j]);
//     }
//   }
//
//   mg->Draw("APE0");
//
//   gStyle->SetTitleFont(62, "XY");
//   gStyle->SetTitleSize(0.04, "XY");
//
//   SetAxisMultiGraph(mg, target);
//   mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
//   mg->GetYaxis()->CenterTitle();
//   mg->GetYaxis()->SetTitleOffset(1.3);
//
//   mg->Draw("APE0");
//   TLegend* legendTarget = GenTLegendTarget(g);
//   TLegend* legendNPion = GenTLegendNPion(g, m, 0.42);
//   legendTarget->Draw();
//   legendNPion->Draw();
//
//   c->Print(plotDirectory + Form("Pt_broad_%s.pdf", target));
//   delete c;
//
// }


void MeanPt2(TFile* inputFile, TFile* outputFile, char target[]) {

  TH1F* meanPt2 = new TH1F("meanPt2", "", N_Zh, Zh_BINS);
  for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion
    for(int NuCounter = 0 ; NuCounter < N_Nu ; NuCounter++) { // Loops in every Q2 bin
	    for(int Q2Counter = 0 ; Q2Counter < N_Q2 ; Q2Counter++) { // Loops in every Nu bin

        for(int ZhCounter = 0; ZhCounter < N_Zh ; ZhCounter++) { // Loops in every Zh bin
          TH1F* hist = (TH1F*)inputFile->Get(Form("corr_data_Pt2_%s_%i%i%i_%i_CLEAN_INTERPOLATED", target, Q2Counter, NuCounter, ZhCounter, nPion));
          if(hist == NULL) {
            meanPt2->SetBinContent(ZhCounter + 1, 0);
            meanPt2->SetBinError(ZhCounter + 1, 0);
            continue;
          }
          meanPt2->SetBinContent(ZhCounter + 1, hist->GetMean());
          meanPt2->SetBinError(ZhCounter + 1, hist->GetMeanError());
          delete hist;
        } // End Zh loop
        outputFile->cd();
        meanPt2->Write(Form("meanPt2_%s_%i%i_%i",target, Q2Counter, NuCounter, nPion));
        gROOT->cd();
        meanPt2->Reset();
	    } // End Q2 loop
    } // End Nu loop
  } // End number pion event loop

}

void CallMeanPt2(TString inputDirectory, TString outputDirectory) {

  TFile* inputFile  = new TFile(inputDirectory  + "corr_data_Pt2_processed.root", "READ");
  TFile* outputFile = new TFile(outputDirectory + "multiMeanPt2.root", "RECREATE");
  gROOT->cd();

  MeanPt2(inputFile, outputFile, DC);
  MeanPt2(inputFile, outputFile, DFe);
  MeanPt2(inputFile, outputFile, DPb);

  MeanPt2(inputFile, outputFile, C);
  MeanPt2(inputFile, outputFile, Fe);
  MeanPt2(inputFile, outputFile, Pb);

  inputFile->Close();
  outputFile->Close();

}

void PtbroadeningMultiZh(TString inputDirectory, TString plotDirectory) {

  TFile* fsource  = new TFile(inputDirectory + "multiMeanPt2.root", "READ");
  TCanvas* c = new TCanvas("c","",1500,1500);
  c->Draw();
  c->Divide(3,3,0,0);

  int count = 1;
  for(int NuCounter = 0 ; NuCounter < N_Nu ; NuCounter++) { // Loops in every Nu bin
    for(int Q2Counter = 0 ; Q2Counter < N_Q2 ; Q2Counter++) { // Loops in every Q2 bin
            c->cd(count);count++;
      TGraphErrors* g[N_STARGETS][N_PION];

      for(int nPion = 1; nPion <= N_PION; nPion++) { // Loops in every number of pion

        TH1F* histSolid[3];   TH1F* histLiquid[3];  TH1F* histBroadening[3];

        histSolid[0] = (TH1F*) fsource->Get(Form("meanPt2_C_%i%i_%i", Q2Counter, NuCounter, nPion));
        histSolid[1] = (TH1F*) fsource->Get(Form("meanPt2_Fe_%i%i_%i", Q2Counter, NuCounter, nPion));
        histSolid[2] = (TH1F*) fsource->Get(Form("meanPt2_Pb_%i%i_%i", Q2Counter, NuCounter, nPion));

        histLiquid[0] = (TH1F*) fsource->Get(Form("meanPt2_DC_%i%i_%i", Q2Counter, NuCounter, nPion));
        histLiquid[1] = (TH1F*) fsource->Get(Form("meanPt2_DFe_%i%i_%i", Q2Counter, NuCounter, nPion));
        histLiquid[2] = (TH1F*) fsource->Get(Form("meanPt2_DPb_%i%i_%i", Q2Counter, NuCounter, nPion));

        for(Int_t i = 0 ; i < 3 ; i++){
          histBroadening[i] = new TH1F(Form("histBroadening_%i",i),"",N_Zh,Zh_BINS);
          histBroadening[i]->Add(histSolid[i],histLiquid[i],1,-1);
        }

        for(Int_t i = 0 ; i < 3 ; i++){
          g[i][nPion-1] = (TGraphErrors*) TH1TOTGraph(histBroadening[i]);
          SetErrorXNull(g[i][nPion-1], N_Zh);
        }
        g[0][nPion-1]->SetMarkerColor(kRed);
        g[0][nPion-1]->SetLineColor(kRed);
        g[1][nPion-1]->SetMarkerColor(kBlue);
        g[1][nPion-1]->SetLineColor(kBlue);
        g[2][nPion-1]->SetMarkerColor(kBlack);
        g[2][nPion-1]->SetLineColor(kBlack);

        SetXShift(g[0][nPion-1], -0.03 + (- 0.02 + 0.01*nPion), N_Zh);
        SetXShift(g[1][nPion-1], 0     + (- 0.02 + 0.01*nPion), N_Zh);
        SetXShift(g[2][nPion-1], 0.03  + (- 0.02 + 0.01*nPion), N_Zh);

        for(int i = 0; i < N_STARGETS; i++) {
          delete histSolid[i];
          delete histLiquid[i];
          delete histBroadening[i];
        }
      }
      TLegend* legend = new TLegend(0.15,0.75,0.3,0.92,"","NDC");
      legend->AddEntry(g[0][1],"C","lpf");
      legend->AddEntry(g[1][1],"Fe","lpf");
      legend->AddEntry(g[2][1],"Pb","lpf");
      legend->SetFillStyle(0);
      legend->SetBorderSize(0);
      legend->SetTextFont(62);
      legend->SetTextSize(0.04);

      for(int k = 0; k < 3; k++) {
        g[k][0]->SetMarkerStyle(22);
        g[k][0]->SetMarkerSize(1.2);
        g[k][1]->SetMarkerStyle(8);
        g[k][1]->SetMarkerSize(0.95);
        g[k][2]->SetMarkerStyle(27);
        g[k][2]->SetMarkerSize(1.7);
      }

      TLatex* kinematics = new TLatex();
      kinematics->SetTextFont(62);
      kinematics->SetTextSize(0.06);

      TLatex* kinematics2 = new TLatex();
      kinematics2->SetTextFont(62);
      kinematics2->SetTextSize(0.06);

      TPad* p = new TPad("p", "p", 0, 0, 1, 1);
      if(Q2Counter == 0) { p->SetLeftMargin(0.13); }
      else { p->SetLeftMargin(0.0); }
      if(NuCounter == 0) { p->SetTopMargin(0.13); }
      else { p->SetTopMargin(0.0); }
      if(Q2Counter == 2) { p->SetRightMargin(0.13); }
      else { p->SetRightMargin(0.0); }
      if(NuCounter == 2) { p->SetBottomMargin(0.13); }
      else { p->SetBottomMargin(0.0); }

      p->SetGridx(1);
      p->SetGridy(1);
      p->Draw();
      p->cd();


      TMultiGraph* mg = new TMultiGraph();
      for(int i = 0; i < N_STARGETS; i++){
        for(int j = 0; j < N_PION; j++){
          mg->Add(g[i][j]);
        }
      }


      mg->Draw("APE0");

      gStyle->SetTitleFont(62, "XY");
      gStyle->SetTitleSize(0.04, "XY");

      mg->GetYaxis()->SetRangeUser(0.,0.085);
      if(NuCounter == 2) {
        mg->GetXaxis()->SetRangeUser(Zh_MIN,Zh_MAX);
        mg->GetXaxis()->SetTitle("z_{h}");
        mg->GetXaxis()->CenterTitle();
        mg->GetXaxis()->SetTitleOffset(1.1);
      }
      if (Q2Counter == 0) {
        mg->GetYaxis()->SetTitle("#Delta P^{2}_{T} [GeV^{2}]");
        mg->GetYaxis()->CenterTitle();
        mg->GetYaxis()->SetTitleOffset(1.3);
      }

      //  mg->SetTitle(";z_{h};#Delta P^{2}_{T} [GeV^{2}]");
      mg->Draw("APE0");
      //legend->Draw();
      //kinematics->SetTextSize(1);
      if(NuCounter == 0) {
        kinematics->DrawLatexNDC(p->GetLeftMargin()+0.19,1.024*(1-p->GetTopMargin()),Form("%.2f<Q^{2}[GeV^{2}]<%.2f ", Q2_BINS[Q2Counter], Q2_BINS[Q2Counter+1]));
      }
      kinematics2->SetTextAngle(270);
      if(Q2Counter == 2) {
        kinematics2->DrawLatexNDC(1.024*(0.9-p->GetLeftMargin()),.7,Form("%.2f<#nu[GeV^{2}]<%.2f ", Nu_BINS[NuCounter], Nu_BINS[NuCounter+1]));
      }
      //delete mg;
      // for (int i = 0; i < 3; i++) {
      //   for (int j = 0; j < 3; j++) {
      //     std::cout << "PhiPQ integration4" << i << j<< std::endl;
      //     delete g[i][j];
      //   }
      // }
    } // End Nu loop
  } // End Q2 loop

  c->Print(plotDirectory + "Pt_broad_multiZh.pdf");
  fsource->Close();

}
