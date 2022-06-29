// This code generate a tuple with all the events of experimental data
// Saves the electron variables and calculate for the hadrons variables
// calculate vectorial momentum and use it as the hadron momentum for the event
// The code require that you have the number of the event saved in the data tupleName
// if you don't have it you can check by for the paricle has the same Q2 and Nu instead
// It can be compiled with
// g++ -Wall -fPIC -I./include `root-config --cflags` VecSum.cpp -o ./bin/VecSum  `root-config --glibs`
// For the target name use (C,Fe,Pb)

#include <iostream>
#include <string>
#include "Binning.h"
#include "TMath.h"
#include "TString.h"
#include "TFile.h"
#include "TNtuple.h"
#include "TVector2.h"
#include "TStopwatch.h"
#include "TROOT.h"

int main(int argc, char* argv[]) {

  if(argc != 2) {
    std::cout << "Insert (just) the target name as a parameter" << std::endl;
    return 0;
  }

  TStopwatch t;
  std::cout << "Start" << std::endl;

  std::string target = argv[1];
  // Creating a array of chars instead of a string to use Form method
  int n = target.length();
  char targetArr[n + 1];
  strcpy(targetArr, target.c_str());

  TFile* file = new TFile(Form("~/proyecto/Pt2Broadening_multi-pion/Data/PiPlusData_%s.root", targetArr), "READ");
  TNtuple* tuple = (TNtuple*)file->Get("ntuple_data");

  int tmpCounter = 0; // Counts how many partivles there is in the event
  float tmpEvnt, evnt, Q2Evnt, NuEvnt, ZhEvnt, Pt2Evnt, PhiEvnt, YCEvnt, VCEvnt;
  float tmpZh[5], tmpPt[5], tmpPhi[5] ;

  const char* VarList = "Q2:Nu:Zh:Pt2:PhiPQ:YC:VC_TM";
  // Variables to fill the tuple
  float *vars         = new Float_t[7];
  // Read the necesary variables
  tuple->SetBranchAddress("Q2",&Q2Evnt);
  tuple->SetBranchAddress("Nu",&NuEvnt);
  tuple->SetBranchAddress("Zh",&ZhEvnt);
  tuple->SetBranchAddress("Pt2",&Pt2Evnt);
  tuple->SetBranchAddress("PhiPQ",&PhiEvnt);
  tuple->SetBranchAddress("YC",&YCEvnt);
  tuple->SetBranchAddress("VC_TM",&VCEvnt);
  tuple->SetBranchAddress("NEvnt",&evnt);

  gROOT->cd();

  TNtuple* ntuplePion[5];

  for(int i = 0; i < 5; i++) {
    ntuplePion[i] = new TNtuple(Form("ntuple_%i_pion", i + 1),"",VarList);
  }

  for(int i = 0; i < tuple->GetEntries() ; i++) { // Loops in every detected paricle
    tuple->GetEntry(i);
    vars[0] = Q2Evnt;
    vars[1] = NuEvnt;
    vars[2] = ZhEvnt;
    vars[3] = Pt2Evnt;
    vars[4] = PhiEvnt;
    vars[5] = YCEvnt;
    vars[6] = VCEvnt;
    tmpZh[0]  = vars[2];
    tmpPt[0]  = TMath::Sqrt(Pt2Evnt);
    tmpPhi[0] = vars[4];
    tmpEvnt = evnt;
    tuple->GetEntry(i + 1);
    while(tmpEvnt == evnt) { // Check all the paricles in the event
      tmpCounter++;
      tmpZh[tmpCounter]  = ZhEvnt;
      tmpPt[tmpCounter]  = TMath::Sqrt(Pt2Evnt);
      tmpPhi[tmpCounter] = PhiEvnt;
      if(i + 1 + tmpCounter >= tuple->GetEntries() ){ break; }
      tuple->GetEntry(i + 1 + tmpCounter);
    }
    if(tmpCounter == 0) {
      tuple->GetEntry(i);
      ntuplePion[0]->Fill(vars);
    }else {
      tuple->GetEntry(i);
      vars[2] = 0;
      TVector2* vec = new TVector2(0,0);
      for(int k = 0; k <= tmpCounter; k++) {
        // Calculate de tranvers momentum vector
        TVector2 *tmpVec = new TVector2(tmpPt[k]*TMath::Cos((tmpPhi[k] + 180)*TMath::DegToRad()), tmpPt[k]*TMath::Sin((tmpPhi[k] + 180)*TMath::DegToRad()));
        // Sum the vector and save the sum of Zh
        vars[2] += tmpZh[k];
        *vec += *tmpVec;
        //vecTemp->Print();
        delete tmpVec;
      }
      // Save the Pt2 of the sum vector
      vars[3] = std::pow(vec->Mod(),2);
      // Save the PhiPQ of the sum vector
      vars[4] = vec->Phi()*TMath::RadToDeg()-180;
      delete vec;
      ntuplePion[tmpCounter]->Fill(vars);
    }
    // Jump to the next event
    i += tmpCounter;
    tmpCounter = 0;
  } // End paricle loop

  // Save the tuples
  TFile* fOutput = new TFile(Form("~/proyecto/Pt2Broadening_multi-pion/Data/VecSum_%s2.root", targetArr), "RECREATE");
  fOutput->cd();

  for(int i = 0; i < N_PION +1 ; i++) {
    ntuplePion[i]->Write();
  }

  gROOT->cd();
  fOutput->Close();
  std::cout << "Done." << std::endl;
  file->Close();
  t.Print();

}
