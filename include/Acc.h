#include <iostream>
#include <string>
#include "Binning.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TNtuple.h"
#include "TVector2.h"
#include "TStopwatch.h"
#include "TROOT.h"
#include "TH1.h"
#include "TCut.h"
#include "TEventList.h"

// Generated a matrix of the number of events of N generated pion detected as M pion events for each bean
// matrix[generated][detected][Q2Bin][NuBin][ZhBin][Pt2Bin][PhiBin];

void GenMatrixAcc(char targetArr[], int matrix[5][5][N_Q2+1][N_Nu+1][N_Zh+1][N_Pt2+1][N_Phi+1]) {

  std::cout << "Start" << std::endl;

  TString inputName;

  // Create the variables to use inside of the for loops
  int Q2Bin, NuBin, ZhBin, Pt2Bin, PhiBin;
  int genCounter = 0;
  int detCounter = 0;
  int tmpCounter = 0;
  float tmpEvnt;
  float tmpZh[5], tmpPt[5], tmpPhi[5] ;

  for(int folder = 1; folder < 10; folder++) { // Loops in every directory
    for(int sim = 1; sim < 500; sim++) { // Loops in every simulation of the directory
      // Set the name of the file where is the data depends of the target and the folder
      // This are the simulation generated with the code GetSimpleTuple on Hayks simulation
      // https://github.com/utfsm-eg2-data-analysis/GetSimpleTuple
      if(targetArr[0] == 'D' ){
        if(folder < 4) {
          inputName = Form("/eos/user/m/mbarrial/out/GetSimpleTuple_HSim/D2_pb%i/prunedD_%i.root", folder, sim);
        } else {
          inputName = Form("/eos/user/m/mbarrial/out/GetSimpleTuple_HSim/D2_pb%i_yshiftm03/prunedD_%i.root", folder, sim);
        }
      } else {
        if(folder < 4) {
          inputName = Form("/eos/user/m/mbarrial/out/GetSimpleTuple_HSim/%s%i/pruned%s_%i.root", targetArr,folder,targetArr,sim);
        } else {
          inputName = Form("/eos/user/m/mbarrial/out/GetSimpleTuple_HSim/%s%i_yshiftm03/pruned%s_%i.root", targetArr,folder,targetArr,sim);
        }
      }
      //std::cout << "Checking directory " << folder << "  " << sim << std::endl;
      //inputName = Form("/home/matias/proyecto/Piones/Data/Simul/pruned%s_%i.root",targetArr ,sim);
      // Open the file and check if it's exist
      TFile* fSource = new TFile(inputName,"READ");
      if (fSource->IsZombie()) {
        fSource->Close();
        continue;
      }
      // Open the tuple and check if it's exist
      TNtuple* tuple = (TNtuple*)fSource->Get("ntuple_sim");
      if(tuple == NULL) {
        delete tuple;
        fSource->Close();
        continue;
      }
      // Set the variables to read the tuple
      Float_t *vars = new Float_t[8];
      tuple->SetBranchAddress("evnt",&vars[0]);
      tuple->SetBranchAddress("mc_pid",&vars[1]);
      tuple->SetBranchAddress("pid",&vars[2]);
      tuple->SetBranchAddress("mc_Q2",&vars[3]);
      tuple->SetBranchAddress("mc_Nu",&vars[4]);
      tuple->SetBranchAddress("mc_Zh",&vars[5]);
      tuple->SetBranchAddress("mc_Pt2",&vars[6]);
      tuple->SetBranchAddress("mc_PhiPQ",&vars[7]);

      for(int i = 0; i < tuple->GetEntries(); i++) {// Loops in every generated particle
        tuple->GetEntry(i);
        // Check the bin of Q2 for the event
        Q2Bin = N_Q2;
        NuBin = N_Nu;
        for(int j = 0; j < N_Q2; j++ ) { // Loops in every Q2 bin
          if(vars[3] > Q2_BINS[j] && vars[3] < Q2_BINS[j+1]) {
            Q2Bin = j;
            break;
          }
        }
        // Check the bin of Nu for the event
        for(int j = 0; j < N_Nu; j++ ) { // Loops in every Nu bin
          if(vars[4] > Nu_BINS[j] && vars[4] < Nu_BINS[j+1]) {
            NuBin = j;
            break;
          }
        }
        // Check if the generated paricle is a pion+
        if(vars[1] == 211 ){
          // save the angle PhiPQ,Zh and Pt if it's a pion
          tmpZh[0]  = vars[5];
          tmpPt[0]  = TMath::Sqrt(vars[6]);
          tmpPhi[0] = vars[7];
          genCounter++;
        }
        // Check if the detected paricle is a pion+
        if(vars[2] == 211 ){ detCounter++; }
        tmpEvnt = vars[0];
        tuple->GetEntry(i + 1);
        // Check if the next particle cames from the same event
        while(tmpEvnt == vars[0] ) {
          if(vars[1] == 211 ) {
            // save the angle PhiPQ,Zh and Pt of every pion in the event
            tmpZh[genCounter]  = vars[5];
            tmpPt[genCounter]  = TMath::Sqrt(vars[6]);
            tmpPhi[genCounter] = vars[7];
            genCounter++;
          }
          if(vars[2] == 211 ) { detCounter++; }
          tmpCounter++;
          tmpEvnt = vars[0];
          // Go to the next particle
          tuple->GetEntry(i + 1 + tmpCounter);
          if(i + 1 + tmpCounter > tuple->GetEntries() ){
            break;
          }
        }
        vars[5] = 0;
        TVector2* vec = new TVector2(0,0);
        for(int k = 0; k < genCounter; k++) {
          // Calculate de tranvers momentum vector
          TVector2 *tmpVec = new TVector2(tmpPt[k]*TMath::Cos((tmpPhi[k] + 180)*TMath::DegToRad()), tmpPt[k]*TMath::Sin((tmpPhi[k] + 180)*TMath::DegToRad()));
          // Sum the vector and save the sum of Zh
          vars[5] += tmpZh[k];
          *vec += *tmpVec;
          //vecTemp->Print();
          delete tmpVec;
        }
        // Save the Pt2 of the sum vector
        vars[6] = std::pow(vec->Mod(),2);
        // Save the PhiPQ of the sum vector
        vars[7] = vec->Phi()*TMath::RadToDeg()-180;
        delete vec;
        // Check the bin of Zh sum for the event
        ZhBin = N_Zh;
        Pt2Bin = N_Pt2;
        PhiBin = N_Phi;
        for(int j = 0; j < N_Zh; j++ ) { // Loops in every Zhsum bin
          if(vars[5] > Zh_BINS[j] && vars[5] < Zh_BINS[j+1]) {
            ZhBin = j;
            break;
          }
        }
        // Check the bin of Pt2 for the event
        for(int j = 0; j < N_Pt2; j++ ) { // Loops in every Pt2 bin
          if(vars[6] > Pt2_BINS[j] && vars[6] < Pt2_BINS[j+1]) {
            Pt2Bin = j;
            break;
          }
        }
        // Check the bin of Phi for the event
        for(int j = 0; j < N_Phi; j++ ) { // Loops in every PhiPQ bin
          if(vars[7] > Phi_BINS[j] && vars[7] < Phi_BINS[j+1]) {
            PhiBin = j;
            break;
          }
        }
        // Add one to the correspond matrix element
        matrix[genCounter][detCounter][Q2Bin][NuBin][ZhBin][Pt2Bin][PhiBin]++;
        genCounter = 0;
        detCounter = 0;
        // Jump to the next event
        i += tmpCounter;
        tmpCounter = 0;
      }// End Event Loop
      delete tuple;
      fSource->Close();
    }//End sim loop
    std::cout << "Directory " << folder << " checked" << std::endl;
  }//End folder loop
}

// Set all the values of the matrix in 0, this to have defined the inicial value and can use matrix+= later

void SetMatrixValue0(int matrix[5][5][N_Q2+1][N_Nu+1][N_Zh+1][N_Pt2+1][N_Phi+1]) {

  for(int PhiCounter = 0; PhiCounter <= N_Phi; PhiCounter++) {
    for(int Q2Counter = 0; Q2Counter <= N_Q2; Q2Counter++) {
      for(int NuCounter = 0; NuCounter <= N_Nu; NuCounter++) {
        for(int ZhCounter = 0; ZhCounter <= N_Zh; ZhCounter++) {
          for(int Pt2Counter = 0; Pt2Counter <= N_Pt2; Pt2Counter++) {
            for(int i = 0; i < 5; i++) {
              for(int j = 0; j < 5; j++) {
                matrix[i][j][Q2Counter][NuCounter][ZhCounter][Pt2Counter][PhiCounter] = 0;
              }
            }
          }
        }
      }
    }
  }

}

// Generate a PhiPQ histogram storing how many N pions events were generated in the seleced bin
void GenThrown(int matrix[5][5][N_Q2+1][N_Nu+1][N_Zh+1][N_Pt2+1][N_Phi+1], TH1F* histThrown, int gen, int Q2Bin, int NuBin, int ZhBin, int Pt2Bin) {

  int totalGen = 0;

  for(int PhiCounter = 0; PhiCounter < N_Phi; PhiCounter++) { // Loops in every PhiPQ bin
    // Sum all the generated N pion event
    for (int j = 0; j <= N_PION; j++) {
      totalGen += matrix[gen][j][Q2Bin][NuBin][ZhBin][Pt2Bin][PhiCounter];
    }
    histThrown->Fill(-179.5 + Delta_Phi*PhiCounter, totalGen);
    totalGen = 0;
  }

}

// Generate a PhiPQ histogram storing how many N pions events were detected in the seleced bin
void GenTotDetected(int matrix[5][5][N_Q2+1][N_Nu+1][N_Zh+1][N_Pt2+1][N_Phi+1], TH1F* histTotDetected, int dec, int Q2Bin, int NuBin, int ZhBin, int Pt2Bin) {

  int totalDec = 0;

  for(int PhiCounter = 0; PhiCounter < N_Phi; PhiCounter++) { // Loops in every PhiPQ bin
    // Sum all the events detected as an N pion event
    for (int j = 0; j <= N_PION; j++) {
      totalDec += matrix[j][dec][Q2Bin][NuBin][ZhBin][Pt2Bin][PhiCounter];
    }
    histTotDetected->Fill(-179.5 + Delta_Phi*PhiCounter, totalDec);
    totalDec = 0;
  }

}

// Generate a PhiPQ histogram storing how many N pions events were correctly detected as an N pion event in the selected bin
void GenDetected(int matrix[5][5][N_Q2+1][N_Nu+1][N_Zh+1][N_Pt2+1][N_Phi+1], TH1F* histDetected, int nPion, int Q2Bin, int NuBin, int ZhBin, int Pt2Bin) {

  for(int PhiCounter = 0; PhiCounter < N_Phi; PhiCounter++) { // Loops in every PhiPQ bin
    histDetected->Fill(-179.5 + Delta_Phi*PhiCounter, matrix[nPion][nPion][Q2Bin][NuBin][ZhBin][Pt2Bin][PhiCounter]);
  }

}

// If the histogram if empty return 1 if not return 0
int EmptyHist(TH1F* h) {

  int empty = 0;
  for(int i = 1 ; i <= h->GetNbinsX() ; i++) {
    if(h->GetBinContent(i) == 0){ empty++; }
  }
  if(empty == h->GetNbinsX()) { return 1; }
  else { return 0; }

}

// void rec_histo_process(TH1F* h) {
//   /*Here is applied the condition that Naccept>1*/
//   for(Int_t bin = 1; bin <= h->GetNbinsX(); bin++) {
//     if(h->GetBinContent(bin) == 1) {
//       h->SetBinContent(bin, 0);
//       h->SetBinError(bin, 0);
//     }
//   }
// }

// If Acceptance Factor > 1  set it to 0
void AccHist1(TH1F* h) {

  /*Here is applied the condition that acc<1*/
  for(Int_t bin = 1; bin <= h->GetNbinsX(); bin++) {
    if(h->GetBinContent(bin) > 1) {
      h->SetBinContent(bin, 0);
      h->SetBinError(bin, 0);
    }
  }

}
