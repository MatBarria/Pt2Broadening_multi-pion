#include "Acc.h"

int main(int argc, char* argv[]) {

  if(argc != 2) {
    std::cout << "Insert (just) the target name as a parameter" << std::endl;
    return 0;
  }

  TStopwatch t;

  // For the Target name use (C,Fe,Pb) for the solids targets and (DC,DFe,DPb) for the liquid target
  std::string target = argv[1];
  // Creating a array of chars instead of a string to use Form method
  int n = target.length();
  char targetArr[n + 1];
  strcpy(targetArr, target.c_str());

  int matrix[5][5][N_Q2+1][N_Nu+1][N_Zh+1][N_Pt2+1][N_Phi+1];

  SetMatrixValue0(matrix);
  GenMatrixAcc(targetArr,matrix);

  // Select the solid target

  TFile *fileOutput= new TFile(Form("/eos/user/m/mbarrial/Data/Acc/Gen_Dec_%s.root", targetArr), "RECREATE");
  //TFile *fileOutput= new TFile("Acce.root","RECREATE");
  gROOT->cd();

  // Create all the necessary histograms
  TH1F *histDetected    = new TH1F("histDetected",    "", N_Phi, -180, 180);
  TH1F *histTotDetected = new TH1F("histTotDetected", "", N_Phi, -180, 180);
  TH1F *histThrown      = new TH1F("histThrown",      "", N_Phi, -180, 180);

  // Store the sum of the weights A.K.A the errors
  histThrown->Sumw2();
  histTotDetected->Sumw2();
  histDetected->Sumw2();

  for(int gen = 1; gen <= N_PION ; gen++) { // Loops in every number of generated pions
    for(int Q2Counter = 0; Q2Counter < N_Q2; Q2Counter++) { // Loops in every Q2 bin
      for(int NuCounter = 0; NuCounter < N_Nu; NuCounter++) { // Loops in every Nu bin
        for(int ZhCounter = 0; ZhCounter < N_Zh; ZhCounter++) { // Loops in every Zh bin
          std::cout << "Bin selected: " << gen << Q2Counter << NuCounter << ZhCounter << std::endl;
          for(int Pt2Counter = 0; Pt2Counter < N_Pt2; Pt2Counter++) { // Loops in every Pt2 bin

            // Generate histograms of the all dectected pion, all generated pion, and the pions that was correct dectected
            GenThrown(      matrix, histThrown     , gen, Q2Counter, NuCounter, ZhCounter, Pt2Counter);
            GenDetected(    matrix, histDetected   , gen, Q2Counter, NuCounter, ZhCounter, Pt2Counter);
            GenTotDetected( matrix, histTotDetected, gen, Q2Counter, NuCounter, ZhCounter, Pt2Counter);

            // Save the histograms in the output file
            fileOutput->cd();

            histThrown->Write(Form("histThrown_%s_%i%i%i%i_%i",           targetArr, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));
            histDetected->Write(Form("histDetected_%s_%i%i%i%i_%i",       targetArr, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));
            histTotDetected->Write(Form("histTotDetected_%s_%i%i%i%i_%i", targetArr, Q2Counter, NuCounter, ZhCounter, Pt2Counter, gen));

            gROOT->cd();

            // Set the histograms values to 0
            histThrown->Reset();
            histDetected->Reset();
            histTotDetected->Reset();

          }
        }
      }
    }
  }
  fileOutput->Close();
  t.Print();

}
